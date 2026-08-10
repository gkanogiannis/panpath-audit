use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, Seek, SeekFrom, Write};
use std::os::unix::fs::FileExt;
use std::path::Path;

use crate::sequence::{Orientation, complement, is_ambiguous, is_iupac, walk_steps};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
/// GFA traversal record kind associated with an audited identifier.
pub enum RecordType {
    Path,
    Walk,
}

impl RecordType {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Path => "P",
            Self::Walk => "W",
        }
    }
}

/// Result of comparing all addressed traversal ranges with their source.
pub struct PathComparison {
    pub divergence: Option<usize>,
    pub base_count: usize,
    pub location: Option<GraphLocation>,
    pub record_type: RecordType,
    pub source_range: Option<(usize, usize)>,
    pub range_mismatch: Option<RangeLengthMismatch>,
    pub graph_ambiguous_bases: u64,
    pub not_embedded_bases: u64,
    pub addressed_source_digest: crate::digests::SequenceDigest,
    pub graph_digest: crate::digests::SequenceDigest,
}

/// Declared and reconstructed lengths for an invalid W-line range.
pub struct RangeLengthMismatch {
    pub position: usize,
    pub range_length: usize,
    pub walk_length: usize,
}

/// Materialized graph sequence for one addressed source range.
pub struct TraversalPart {
    pub start: Option<usize>,
    pub end: Option<usize>,
    pub sequence: Vec<u8>,
}

/// Length, ambiguity, and digest measurements for a traversal.
pub struct TraversalMeasure {
    pub bases: u64,
    pub ambiguous_bases: u64,
    pub digest: crate::digests::SequenceDigest,
}

/// Graph coordinate corresponding to a reconstructed traversal position.
pub struct GraphLocation {
    pub segment: String,
    pub orientation: Orientation,
    pub traversal_position: usize,
    pub segment_position: usize,
}

/// Disk-backed graph index and validated traversal metadata.
pub struct Graph {
    pub path_names: BTreeSet<String>,
    pub path_record_count: usize,
    pub walk_record_count: usize,
    pub version: String,
    pub compression: &'static str,
    pub unspecified_overlaps_as_blunt: bool,
    path_names_p: BTreeSet<String>,
    records: BTreeMap<String, Vec<TraversalRef>>,
    segment_index: SegmentIndex,
    traversals: Spool,
    segment_sequences: Spool,
}

/// Stateful reader for graph sequences stored in anonymous spools.
pub struct GraphReader<'a> {
    graph: &'a Graph,
    traversals: File,
    segment_sequences: File,
}

#[derive(Clone, Copy)]
struct RecordRef {
    offset: u64,
    length: u32,
}

#[derive(Clone, Copy)]
struct TraversalRef {
    record: RecordRef,
    start: Option<usize>,
    end: Option<usize>,
}

struct SegmentIndex {
    numeric: Vec<PackedRecordRef>,
    named: HashMap<String, RecordRef>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
struct PackedRecordRef([u8; 12]);

struct Spool {
    file: File,
}

impl Graph {
    /// Validates and indexes a GFA without retaining segment sequences in memory.
    pub fn read(
        path: &Path,
        sources: &BTreeMap<String, usize>,
        temp_dir: &Path,
    ) -> Result<Self, String> {
        let compression = compression(path)?;
        let mut traversals = Spool::create(temp_dir)?;
        let mut segment_sequences = Spool::create(temp_dir)?;
        let mut graph = Self {
            path_names: BTreeSet::new(),
            path_record_count: 0,
            walk_record_count: 0,
            version: "inferred-1.x".into(),
            compression,
            unspecified_overlaps_as_blunt: false,
            path_names_p: BTreeSet::new(),
            records: BTreeMap::new(),
            segment_index: SegmentIndex::new(),
            traversals: Spool::create(temp_dir)?,
            segment_sequences: Spool::create(temp_dir)?,
        };
        let mut pending: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        let mut walk_names = BTreeSet::new();
        let mut walk_ranges: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
        let mut walk_counts: BTreeMap<String, usize> = BTreeMap::new();
        let mut unspecified_walks = BTreeSet::new();
        let mut errors = Vec::new();
        let mut reader = open(path)?;
        let mut line = String::new();
        let mut line_number = 0;

        while read_line(&mut reader, &mut line)? {
            line_number += 1;
            let fields: Vec<_> = line.trim_end_matches(['\n', '\r']).split('\t').collect();
            match fields.as_slice() {
                ["H", tags @ ..] => {
                    if let Some(version) = tags.iter().find_map(|tag| tag.strip_prefix("VN:Z:")) {
                        graph.version = version.into();
                    }
                }
                ["S", name, sequence, ..] => {
                    let duplicate = graph.segment_index.contains(name);
                    if duplicate {
                        errors.push(format!("DUPLICATE_SEGMENT identifier={name}"));
                    } else {
                        graph.segment_index.mark(name);
                    }
                    pending.remove(*name);
                    if *sequence == "*" {
                        errors.push(format!(
                            "SEQUENCELESS_SEGMENT line={line_number} identifier={name}"
                        ));
                    } else if let Some(symbol) = sequence.bytes().find(|base| !is_iupac(*base)) {
                        errors.push(format!(
                            "INVALID_NUCLEOTIDE line={line_number} symbol={}",
                            char::from(symbol)
                        ));
                    } else {
                        let reference = segment_sequences.append(sequence.as_bytes())?;
                        graph.segment_index.insert(name, reference);
                    }
                }
                ["P", name, steps, overlap, ..] => {
                    graph.path_record_count += 1;
                    graph.path_names.insert((*name).into());
                    if !graph.path_names_p.insert((*name).into()) {
                        errors.push(format!("DUPLICATE_PATH identifier={name}"));
                    }
                    let parsed = match parse_path_steps(steps) {
                        Ok(parsed) => parsed,
                        Err(()) => {
                            errors.push(format!("MALFORMED_PATH line={line_number}"));
                            Vec::new()
                        }
                    };
                    let blunt = *overlap == "*"
                        || (parsed.len() > 1
                            && overlap.split(',').count() == parsed.len() - 1
                            && overlap.split(',').all(|value| value == "0M"));
                    if !blunt {
                        errors.push(format!(
                            "UNSUPPORTED_OVERLAP line={line_number} overlap={overlap}"
                        ));
                    }
                    graph.unspecified_overlaps_as_blunt |= *overlap == "*";
                    track_references(name, &parsed, &graph.segment_index, &mut pending);
                    let reference =
                        traversals.append(line.trim_end_matches(['\n', '\r']).as_bytes())?;
                    graph.records.insert(
                        (*name).into(),
                        vec![TraversalRef {
                            record: reference,
                            start: None,
                            end: None,
                        }],
                    );
                }
                ["W", sample, haplotype, sequence, start, end, walk, ..] => {
                    graph.walk_record_count += 1;
                    if haplotype.parse::<u64>().is_err() {
                        errors.push(format!(
                            "INVALID_HAPLOTYPE line={line_number} value={haplotype}"
                        ));
                    }
                    let name = format!("{sample}#{haplotype}#{sequence}");
                    graph.path_names.insert(name.clone());
                    walk_names.insert(name.clone());
                    *walk_counts.entry(name.clone()).or_default() += 1;
                    let parsed = match parse_walk_steps(walk) {
                        Ok(parsed) => parsed,
                        Err(()) => {
                            errors.push(format!("MALFORMED_WALK line={line_number}"));
                            Vec::new()
                        }
                    };
                    track_references(&name, &parsed, &graph.segment_index, &mut pending);
                    let range = match (*start, *end) {
                        ("*", "*") => {
                            unspecified_walks.insert(name.clone());
                            (None, None)
                        }
                        _ => match (start.parse::<usize>(), end.parse::<usize>()) {
                            (Ok(start), Ok(end)) if start < end => {
                                walk_ranges
                                    .entry(name.clone())
                                    .or_default()
                                    .push((start, end));
                                if let Some(source_length) = sources.get(&name)
                                    && end > *source_length
                                {
                                    errors.push(format!("W_RANGE_OUT_OF_BOUNDS identifier={name} end={end} source_length={source_length}"));
                                }
                                (Some(start), Some(end))
                            }
                            _ => {
                                errors.push(format!("INVALID_W_RANGE line={line_number}"));
                                (None, None)
                            }
                        },
                    };
                    let reference =
                        traversals.append(line.trim_end_matches(['\n', '\r']).as_bytes())?;
                    graph.records.entry(name).or_default().push(TraversalRef {
                        record: reference,
                        start: range.0,
                        end: range.1,
                    });
                }
                ["W", ..] => errors.push(format!("MALFORMED_W_RECORD line={line_number}")),
                _ => {}
            }
            line.clear();
        }

        for identifier in graph.path_names_p.intersection(&walk_names) {
            errors.push(format!("DUPLICATE_EMBEDDED_PATH identifier={identifier}"));
        }
        for identifier in unspecified_walks {
            if walk_counts[&identifier] > 1 {
                errors.push(format!("AMBIGUOUS_W_RANGES identifier={identifier}"));
            }
        }
        for (identifier, ranges) in &mut walk_ranges {
            ranges.sort_unstable();
            for pair in ranges.windows(2) {
                if pair[1].0 < pair[0].1 {
                    errors.push(format!("W_RANGE_OVERLAP identifier={identifier}"));
                }
            }
        }
        for records in graph.records.values_mut() {
            records.sort_by_key(|record| record.start);
        }
        for (segment, paths) in pending {
            for path in paths {
                errors.push(format!("UNRESOLVED_SEGMENT path={path} segment={segment}"));
            }
        }
        if !errors.is_empty() {
            return Err(errors.join("\n"));
        }

        std::mem::swap(&mut graph.traversals, &mut traversals);
        std::mem::swap(&mut graph.segment_sequences, &mut segment_sequences);
        Ok(graph)
    }

    pub fn record_type(&self, identifier: &str) -> RecordType {
        if self.path_names_p.contains(identifier) {
            RecordType::Path
        } else {
            RecordType::Walk
        }
    }

    /// Opens independent spool handles for traversal operations.
    pub fn reader(&self) -> Result<GraphReader<'_>, String> {
        Ok(GraphReader {
            graph: self,
            traversals: self
                .traversals
                .file
                .try_clone()
                .map_err(|error| error.to_string())?,
            segment_sequences: self
                .segment_sequences
                .file
                .try_clone()
                .map_err(|error| error.to_string())?,
        })
    }
}

impl GraphReader<'_> {
    /// Streams an indexed traversal and compares its addressed ranges with source bytes.
    pub fn compare_traversal(
        &mut self,
        identifier: &str,
        source: &[u8],
    ) -> Result<PathComparison, String> {
        let records = self.graph.records.get(identifier).ok_or("missing path")?;
        let record_type = self.graph.record_type(identifier);
        let source_range =
            records
                .iter()
                .try_fold((usize::MAX, 0), |(minimum, maximum), record| {
                    Some((minimum.min(record.start?), maximum.max(record.end?)))
                });
        let mut divergence = None;
        let mut location = None;
        let mut range_mismatch = None;
        let mut base_count = 0_usize;
        let mut graph_ambiguous_bases = 0_u64;
        let uses_full_source = record_type == RecordType::Path
            || records
                .iter()
                .any(|record| record.start.is_none() || record.end.is_none());
        let addressed_bases = if uses_full_source {
            source.len() as u64
        } else {
            records
                .iter()
                .map(|record| (record.end.unwrap() - record.start.unwrap()) as u64)
                .sum()
        };
        let mut addressed_source_hasher = crate::digests::SequenceHasher::new();
        if uses_full_source {
            addressed_source_hasher.update(source);
        } else {
            for record in records {
                addressed_source_hasher.update(&source[record.start.unwrap()..record.end.unwrap()]);
            }
        }
        let mut graph_hasher = crate::digests::SequenceHasher::new();

        for traversal in records {
            let record = read_record(&mut self.traversals, traversal.record)?;
            let steps = traversal_steps(&record, record_type)?;
            let start = traversal.start.unwrap_or(0);
            let range_length = traversal
                .end
                .zip(traversal.start)
                .map(|(end, start)| end - start);
            let mut local = 0_usize;
            visit_steps(steps, record_type, |name, orientation| {
                let reference = self
                    .graph
                    .segment_index
                    .get(name)
                    .ok_or("missing segment")?;
                let segment = read_record(&mut self.segment_sequences, reference)?;
                let bases: Box<dyn Iterator<Item = (usize, u8)>> = match orientation {
                    Orientation::Forward => Box::new(segment.bytes().enumerate()),
                    Orientation::Reverse => Box::new(segment.bytes().rev().enumerate()),
                };
                for (index, stored) in bases {
                    let graph_base = match orientation {
                        Orientation::Forward => stored,
                        Orientation::Reverse => complement(stored)?,
                    };
                    graph_ambiguous_bases += u64::from(is_ambiguous(graph_base));
                    graph_hasher.update(&[graph_base]);
                    let source_position = start + local;
                    let outside_range = range_length.is_some_and(|length| local >= length);
                    if divergence.is_none()
                        && (outside_range
                            || source.get(source_position).is_none_or(|source_base| {
                                !source_base.eq_ignore_ascii_case(&graph_base)
                            }))
                    {
                        divergence = Some(source_position);
                        location = Some(GraphLocation {
                            segment: name.into(),
                            orientation,
                            traversal_position: index + 1,
                            segment_position: orientation.segment_position(index, segment.len()),
                        });
                    }
                    local += 1;
                    base_count += 1;
                }
                Ok(())
            })?;
            if let Some(expected) = range_length
                && local != expected
            {
                let mismatch = RangeLengthMismatch {
                    position: start + local.min(expected),
                    range_length: expected,
                    walk_length: local,
                };
                if range_mismatch
                    .as_ref()
                    .is_none_or(|current: &RangeLengthMismatch| {
                        mismatch.position < current.position
                    })
                {
                    range_mismatch = Some(mismatch);
                }
            }
        }
        if record_type == RecordType::Path && divergence.is_none() && source.len() != base_count {
            divergence = Some(source.len().min(base_count));
        }
        if let Some(range) = &range_mismatch {
            divergence = Some(divergence.map_or(range.position, |value| value.min(range.position)));
        }
        Ok(PathComparison {
            divergence,
            base_count,
            location,
            record_type,
            source_range,
            range_mismatch,
            graph_ambiguous_bases,
            not_embedded_bases: source.len() as u64 - addressed_bases,
            addressed_source_digest: addressed_source_hasher.finish(),
            graph_digest: graph_hasher.finish(),
        })
    }

    /// Materializes addressed source/graph pairs for detailed alignment.
    pub fn materialize_parts(
        &mut self,
        identifier: &str,
        expected_bases: usize,
    ) -> Result<Vec<TraversalPart>, String> {
        let records = self.graph.records.get(identifier).ok_or("missing path")?;
        let record_type = self.graph.record_type(identifier);
        let mut parts = Vec::with_capacity(records.len());
        for traversal in records {
            let record = read_record(&mut self.traversals, traversal.record)?;
            let steps = traversal_steps(&record, record_type)?;
            let capacity = traversal
                .end
                .zip(traversal.start)
                .map_or(expected_bases, |(end, start)| end - start);
            let mut sequence = Vec::with_capacity(capacity);
            visit_steps(steps, record_type, |name, orientation| {
                let reference = self
                    .graph
                    .segment_index
                    .get(name)
                    .ok_or("missing segment")?;
                let segment = read_record(&mut self.segment_sequences, reference)?;
                match orientation {
                    Orientation::Forward => sequence.extend_from_slice(segment.as_bytes()),
                    Orientation::Reverse => {
                        for base in segment.bytes().rev() {
                            sequence.push(complement(base)?);
                        }
                    }
                }
                Ok(())
            })?;
            parts.push(TraversalPart {
                start: traversal.start,
                end: traversal.end,
                sequence,
            });
        }
        Ok(parts)
    }

    /// Measures and hashes a traversal without materializing it as one sequence.
    pub fn measure(&mut self, identifier: &str) -> Result<TraversalMeasure, String> {
        let records = self.graph.records.get(identifier).ok_or("missing path")?;
        let record_type = self.graph.record_type(identifier);
        let mut result = TraversalMeasure {
            bases: 0,
            ambiguous_bases: 0,
            digest: crate::digests::SequenceHasher::new().finish(),
        };
        let mut hasher = crate::digests::SequenceHasher::new();
        for traversal in records {
            let record = read_record(&mut self.traversals, traversal.record)?;
            let steps = traversal_steps(&record, record_type)?;
            visit_steps(steps, record_type, |name, orientation| {
                let reference = self
                    .graph
                    .segment_index
                    .get(name)
                    .ok_or("missing segment")?;
                let segment = read_record(&mut self.segment_sequences, reference)?;
                result.bases += segment.len() as u64;
                result.ambiguous_bases +=
                    segment.bytes().filter(|base| is_ambiguous(*base)).count() as u64;
                match orientation {
                    Orientation::Forward => hasher.update(segment.as_bytes()),
                    Orientation::Reverse => {
                        for base in segment.bytes().rev() {
                            hasher.update(&[complement(base)?]);
                        }
                    }
                }
                Ok(())
            })?;
        }
        result.digest = hasher.finish();
        Ok(result)
    }
}

fn read_record(file: &mut File, reference: RecordRef) -> Result<String, String> {
    let mut bytes = vec![0; reference.length as usize];
    file.read_exact_at(&mut bytes, reference.offset)
        .map_err(|error| error.to_string())?;
    String::from_utf8(bytes).map_err(|error| error.to_string())
}

impl Spool {
    fn create(temp_dir: &Path) -> Result<Self, String> {
        tempfile::tempfile_in(temp_dir)
            .map(|file| Self { file })
            .map_err(|error| error.to_string())
    }

    fn append(&mut self, bytes: &[u8]) -> Result<RecordRef, String> {
        let offset = self
            .file
            .seek(SeekFrom::End(0))
            .map_err(|error| error.to_string())?;
        self.file
            .write_all(bytes)
            .map_err(|error| error.to_string())?;
        let length = u32::try_from(bytes.len()).map_err(|_| "temporary record too large")?;
        Ok(RecordRef { offset, length })
    }
}

impl SegmentIndex {
    fn new() -> Self {
        Self {
            numeric: vec![PackedRecordRef::missing()],
            named: HashMap::new(),
        }
    }

    fn insert(&mut self, name: &str, reference: RecordRef) {
        if let Ok(identifier) = name.parse::<usize>() {
            if identifier >= self.numeric.len() {
                self.numeric
                    .resize(identifier + 1, PackedRecordRef::missing());
            }
            self.numeric[identifier] = PackedRecordRef::from_record(reference);
        } else {
            self.named.insert(name.into(), reference);
        }
    }

    fn mark(&mut self, name: &str) {
        self.insert(
            name,
            RecordRef {
                offset: u64::MAX,
                length: 0,
            },
        );
    }

    fn contains(&self, name: &str) -> bool {
        name.parse::<usize>()
            .ok()
            .and_then(|identifier| self.numeric.get(identifier))
            .is_some_and(|reference| !reference.is_missing())
            || self.named.contains_key(name)
    }

    fn get(&self, name: &str) -> Option<RecordRef> {
        name.parse::<usize>()
            .ok()
            .and_then(|identifier| {
                self.numeric
                    .get(identifier)
                    .and_then(PackedRecordRef::record)
            })
            .or_else(|| self.named.get(name).copied())
            .filter(|reference| reference.offset != u64::MAX)
    }
}

impl PackedRecordRef {
    fn missing() -> Self {
        Self([u8::MAX; 12])
    }

    fn from_record(record: RecordRef) -> Self {
        let mut bytes = [0; 12];
        bytes[..8].copy_from_slice(&record.offset.to_le_bytes());
        bytes[8..].copy_from_slice(&record.length.to_le_bytes());
        Self(bytes)
    }

    fn is_missing(&self) -> bool {
        *self == Self::missing()
    }

    fn record(&self) -> Option<RecordRef> {
        if self.is_missing() {
            return None;
        }
        let offset = u64::from_le_bytes(self.0[..8].try_into().unwrap());
        let length = u32::from_le_bytes(self.0[8..].try_into().unwrap());
        Some(RecordRef { offset, length })
    }
}

fn open(path: &Path) -> Result<Box<dyn BufRead>, String> {
    crate::input::open(path).map(|(reader, _)| reader)
}

fn compression(path: &Path) -> Result<&'static str, String> {
    crate::input::open(path).map(|(_, compression)| compression.as_str())
}

fn read_line(reader: &mut Box<dyn BufRead>, line: &mut String) -> Result<bool, String> {
    reader
        .read_line(line)
        .map(|count| count != 0)
        .map_err(|error| error.to_string())
}

fn parse_path_steps(input: &str) -> Result<Vec<(String, Orientation)>, ()> {
    input
        .split(',')
        .map(|step| {
            let split = step.len().checked_sub(1).ok_or(())?;
            let (name, orientation) = step.split_at(split);
            let orientation = orientation
                .chars()
                .next()
                .ok_or(())
                .and_then(|value| Orientation::from_path(value).map_err(|_| ()))?;
            if name.is_empty() {
                return Err(());
            }
            Ok((name.into(), orientation))
        })
        .collect()
}

fn parse_walk_steps(input: &str) -> Result<Vec<(String, Orientation)>, ()> {
    walk_steps(input)
        .map(|steps| {
            steps
                .into_iter()
                .map(|(name, orientation)| (name.into(), orientation))
                .collect()
        })
        .map_err(|_| ())
}

fn traversal_steps(record: &str, record_type: RecordType) -> Result<&str, String> {
    let index = if record_type == RecordType::Path {
        2
    } else {
        6
    };
    record
        .split('\t')
        .nth(index)
        .ok_or_else(|| "malformed traversal record".into())
}

/// Visits every segment step in a traversal record in reconstruction order.
fn visit_steps(
    input: &str,
    record_type: RecordType,
    mut visit: impl FnMut(&str, Orientation) -> Result<(), String>,
) -> Result<(), String> {
    if record_type == RecordType::Path {
        for step in input.split(',') {
            let split = step.len().checked_sub(1).ok_or("empty step")?;
            let (name, orientation) = step.split_at(split);
            let orientation = orientation.chars().next().ok_or("missing orientation")?;
            visit(name, Orientation::from_path(orientation)?)?;
        }
        return Ok(());
    }
    let bytes = input.as_bytes();
    let mut start = 0;
    while start < bytes.len() {
        let orientation = Orientation::from_walk(bytes[start])?;
        let name_start = start + 1;
        let end = bytes[name_start..]
            .iter()
            .position(|byte| matches!(byte, b'>' | b'<'))
            .map_or(bytes.len(), |offset| name_start + offset);
        if end == name_start {
            return Err("empty step".into());
        }
        visit(&input[name_start..end], orientation)?;
        start = end;
    }
    Ok(())
}

/// Marks segment references and reports missing segments during preflight.
fn track_references(
    path: &str,
    steps: &[(String, Orientation)],
    segments: &SegmentIndex,
    pending: &mut BTreeMap<String, BTreeSet<String>>,
) {
    for (segment, _) in steps {
        if !segments.contains(segment) {
            pending
                .entry(segment.clone())
                .or_default()
                .insert(path.into());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::PackedRecordRef;

    #[test]
    fn numeric_reference_is_twelve_bytes() {
        assert_eq!(std::mem::size_of::<PackedRecordRef>(), 12);
    }
}
