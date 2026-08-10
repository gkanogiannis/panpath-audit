#![forbid(unsafe_code)]

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::env;
use std::fs;
use std::io::{BufRead, Read, Seek, SeekFrom, Write};
use std::panic::{self, AssertUnwindSafe};
use std::process::ExitCode;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, mpsc};
use std::thread;

mod digests;
mod input;
mod memory;
mod sequence;
mod statistics;
mod streaming_gfa;
use digests::{SequenceDigest, sequence_digest};
use memory::{MemoryBudget, MemoryPermit, MemoryTelemetry};
use sequence::{Orientation, is_ambiguous, is_iupac};
use statistics::{Alignment, EditCounts};
use streaming_gfa::{Graph, GraphReader, RecordType};

const DEFAULT_ALIGNMENT_MAX_CELLS: u64 = 10_000_000;
const DEFAULT_MEMORY_MIB: u64 = 1024;
const TSV_HEADER: &str = "row_type\tidentifier\tstatus\trecord_type\tcode\tmessage\tsource_start_0based\tsource_end_0based_exclusive\tsource_length\tgraph_length\tdivergence_kind\tsource_position_1based\tgraph_position_1based\tsegment\torientation\ttraversal_position_1based\tsegment_position_1based\tsource_sha256\taddressed_source_sha256\tgraph_sha256\tsource_blake3\taddressed_source_blake3\tgraph_blake3\tmatched\tsubstituted\tsource_only\tgraph_only\tnot_embedded\tmissing_path_bases\tmissing_source_bases\tunaligned_divergent_source\tunaligned_divergent_graph\tedit_distance\talignment_status\ttopology_validated\tunspecified_overlaps_as_blunt";
const HELP: &str = "panpath-audit - audit GFA traversals against source FASTA sequences

Usage:
  panpath-audit [OPTIONS] [FASTA ...] GFA

Inputs:
      --fasta FILE                    Add a FASTA with exact identifiers
      --fasta-pansn SAMPLE HAP FILE   Add FASTA using SAMPLE#HAP#identifier
      --fasta-prefix PREFIX FILE      Add FASTA using PREFIXidentifier

Output:
      --format human|json|tsv         Output format [default: human]
      --stats comprehensive           Add exact base statistics and digests

Resources:
      --threads N                     Worker threads [default: CPUs, max 8]
      --memory-mib N                  Tracked sequence memory [default: 1024]
      --alignment-max-cells N         Alignment limit [default: 10000000]
      --temp-dir DIR                  Scratch directory [default: OS temp]

General:
  -h, --help                          Print help
  -V, --version                       Print version
      --                               Treat remaining arguments as paths

Positional FASTA inputs use exact identifiers. Inputs may be plain or gzip.";

#[derive(Clone, Copy, PartialEq, Eq)]
enum Format {
    Human,
    Json,
    Tsv,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum StatsMode {
    Basic,
    Comprehensive { max_cells: u64 },
}

enum SourceMode {
    Exact,
    Pansn { sample: String, haplotype: String },
    Prefix(String),
}

struct SourceSpec {
    path: std::ffi::OsString,
    mode: SourceMode,
}

struct Config {
    format: Format,
    stats_mode: StatsMode,
    threads: usize,
    memory_mib: u64,
    temp_dir: std::path::PathBuf,
    sources: Vec<SourceSpec>,
    gfa_path: std::ffi::OsString,
}

fn main() -> ExitCode {
    let arguments: Vec<_> = env::args_os().skip(1).collect();
    if arguments
        .iter()
        .any(|argument| argument == "--help" || argument == "-h")
    {
        println!("{HELP}");
        return ExitCode::SUCCESS;
    }
    if arguments
        .iter()
        .any(|argument| argument == "--version" || argument == "-V")
    {
        println!("panpath-audit {}", env!("CARGO_PKG_VERSION"));
        return ExitCode::SUCCESS;
    }
    let json_requested = arguments
        .windows(2)
        .any(|pair| pair[0] == "--format" && pair[1] == "json");
    let tsv_requested = arguments
        .windows(2)
        .any(|pair| pair[0] == "--format" && pair[1] == "tsv");
    let comprehensive_requested = arguments
        .windows(2)
        .any(|pair| pair[0] == "--stats" && pair[1] == "comprehensive");
    match run() {
        Ok(code) => code,
        Err(message) => {
            if json_requested {
                let errors: Vec<_> = message
                    .lines()
                    .map(|raw| {
                        let category = diagnostic_category(raw);
                        let line = clean_error_line(raw);
                        let code = line.split_whitespace().next().unwrap_or("ERROR");
                        serde_json::json!({"category": category, "code": code, "message": line})
                    })
                    .collect();
                let mut document = serde_json::json!({
                    "schema": "panpath-audit-report",
                    "schema_version": env!("CARGO_PKG_VERSION"),
                    "state": "invalid",
                    "summary": {
                        "identical": 0,
                        "divergent": 0,
                        "missing_path": 0,
                        "missing_source": 0
                    },
                    "outcomes": [],
                    "errors": errors,
                    "provenance": {"topology_validated": false},
                    "statistics": if comprehensive_requested { Some(serde_json::Value::Null) } else { None }
                });
                if !comprehensive_requested {
                    document.as_object_mut().unwrap().remove("statistics");
                }
                println!("{document}");
            } else if tsv_requested {
                println!("{TSV_HEADER}");
                for raw in message.lines() {
                    let line = clean_error_line(raw);
                    let code = line.split_whitespace().next().unwrap_or("ERROR");
                    TsvRow::error(code, line).print();
                }
            } else {
                for line in message.lines() {
                    eprintln!("{}", clean_error_line(line));
                }
                eprintln!("Try 'panpath-audit --help' for usage.");
            }
            ExitCode::from(failure_exit(&message))
        }
    }
}

fn clean_error_line(line: &str) -> &str {
    line.strip_prefix("SOURCE_INPUT ")
        .or_else(|| line.strip_prefix("GFA_INPUT "))
        .unwrap_or(line)
}

/// Maps an error set to the documented CLI exit-code precedence.
fn failure_exit(message: &str) -> u8 {
    if message
        .lines()
        .any(|line| line.starts_with("SOURCE_INPUT "))
    {
        return 4;
    }
    if message.lines().any(|line| line.starts_with("GFA_INPUT ")) {
        let ambiguous = [
            "DUPLICATE_PATH",
            "DUPLICATE_EMBEDDED_PATH",
            "AMBIGUOUS_W_RANGES",
            "W_RANGE_OVERLAP",
        ];
        return if message.lines().all(|line| {
            let code = clean_error_line(line)
                .split_whitespace()
                .next()
                .unwrap_or("");
            ambiguous.contains(&code)
        }) {
            2
        } else {
            3
        };
    }
    4
}

/// Classifies a diagnostic for structured machine output.
fn diagnostic_category(line: &str) -> &'static str {
    if line.starts_with("SOURCE_INPUT ") {
        "source"
    } else if line.starts_with("GFA_INPUT ") {
        let code = clean_error_line(line)
            .split_whitespace()
            .next()
            .unwrap_or("");
        if matches!(
            code,
            "DUPLICATE_PATH" | "DUPLICATE_EMBEDDED_PATH" | "AMBIGUOUS_W_RANGES" | "W_RANGE_OVERLAP"
        ) {
            "correspondence"
        } else {
            "gfa"
        }
    } else {
        "operational"
    }
}

/// Runs preflight, auditing, statistics collection, and report rendering.
fn run() -> Result<ExitCode, String> {
    let Config {
        format,
        stats_mode,
        threads,
        memory_mib,
        temp_dir,
        sources,
        gfa_path,
    } = parse_arguments(env::args_os().skip(1).collect())?;

    tempfile::tempfile_in(&temp_dir).map_err(|error| {
        format!(
            "TEMP_STORE path={} error={error}",
            temp_dir.to_string_lossy()
        )
    })?;

    let (source_lengths, fasta_compressions, source_errors) = index_sources(&sources);
    let graph_result = Graph::read(std::path::Path::new(&gfa_path), &source_lengths, &temp_dir);
    if !source_errors.is_empty() || graph_result.is_err() {
        let mut errors: Vec<_> = source_errors
            .into_iter()
            .map(|error| format!("SOURCE_INPUT {error}"))
            .collect();
        if let Err(error) = &graph_result {
            errors.extend(error.lines().map(|line| format!("GFA_INPUT {line}")));
        }
        return Err(errors.join("\n"));
    }
    let graph = graph_result.unwrap();
    let path_record_count = graph.path_record_count;
    let walk_record_count = graph.walk_record_count;
    let mut identical = Vec::new();
    let mut divergent = Vec::new();
    let mut missing_paths = Vec::new();
    let mut missing_sources = Vec::new();
    let mut outcome_statistics = BTreeMap::new();
    let mut outcome_digests = BTreeMap::new();

    let (audited_outcomes, memory_telemetry) = audit_sources(
        &sources,
        &source_lengths,
        &graph,
        threads,
        stats_mode,
        memory_mib,
        &temp_dir,
    )?;
    for audited in audited_outcomes {
        let identifier = match &audited.outcome {
            AuditOutcome::Identical(value) => &value.identifier,
            AuditOutcome::Divergent(value) => &value.identifier,
            AuditOutcome::MissingPath(identifier) => identifier,
        };
        if let Some(statistics) = audited.statistics {
            outcome_statistics.insert(identifier.clone(), statistics);
        }
        outcome_digests.insert(identifier.clone(), audited.digests);
        match audited.outcome {
            AuditOutcome::Identical(value) => identical.push(value),
            AuditOutcome::Divergent(value) => divergent.push(value),
            AuditOutcome::MissingPath(identifier) => missing_paths.push(identifier),
        }
    }
    identical.sort_by(|left, right| left.identifier.cmp(&right.identifier));
    divergent.sort_by(|left, right| left.identifier.cmp(&right.identifier));
    missing_paths.sort();
    let source_names: BTreeSet<_> = source_lengths.keys().cloned().collect();
    missing_sources.extend(graph.path_names.difference(&source_names).cloned());
    if !missing_sources.is_empty() {
        let mut reader = graph.reader()?;
        for identifier in &missing_sources {
            let measure = reader.measure(identifier)?;
            let graph_bases = measure.bases;
            if stats_mode != StatsMode::Basic {
                outcome_statistics.insert(
                    identifier.clone(),
                    OutcomeStatistics::missing_source(graph_bases, measure.ambiguous_bases),
                );
            }
            outcome_digests.insert(
                identifier.clone(),
                OutcomeDigests {
                    source: None,
                    addressed_source: None,
                    graph: Some(measure.digest),
                },
            );
        }
    }

    render_completed(CompletedReport {
        format,
        stats_mode,
        threads,
        sources,
        fasta_compressions,
        gfa_path,
        graph,
        path_record_count,
        walk_record_count,
        identical,
        divergent,
        missing_paths,
        missing_sources,
        outcome_statistics,
        outcome_digests,
        memory_telemetry,
    })
}

/// Owns all data needed to render a completed audit.
struct CompletedReport {
    format: Format,
    stats_mode: StatsMode,
    threads: usize,
    sources: Vec<SourceSpec>,
    fasta_compressions: Vec<&'static str>,
    gfa_path: std::ffi::OsString,
    graph: Graph,
    path_record_count: usize,
    walk_record_count: usize,
    identical: Vec<IdenticalOutcome>,
    divergent: Vec<DivergentOutcome>,
    missing_paths: Vec<String>,
    missing_sources: Vec<String>,
    outcome_statistics: BTreeMap<String, OutcomeStatistics>,
    outcome_digests: BTreeMap<String, OutcomeDigests>,
    memory_telemetry: MemoryTelemetry,
}

impl CompletedReport {
    fn exit_code(&self) -> ExitCode {
        if !self.missing_paths.is_empty() || !self.missing_sources.is_empty() {
            ExitCode::from(2)
        } else if !self.divergent.is_empty() {
            ExitCode::from(1)
        } else {
            ExitCode::SUCCESS
        }
    }
}

/// Renders a completed report in the requested format.
fn render_completed(report: CompletedReport) -> Result<ExitCode, String> {
    let exit_code = report.exit_code();
    let CompletedReport {
        format,
        stats_mode,
        threads,
        sources,
        fasta_compressions,
        gfa_path,
        graph,
        path_record_count,
        walk_record_count,
        identical,
        divergent,
        missing_paths,
        missing_sources,
        outcome_statistics,
        outcome_digests,
        memory_telemetry,
    } = report;
    if format == Format::Tsv {
        render_tsv(
            &identical,
            &divergent,
            &missing_paths,
            &missing_sources,
            &outcome_statistics,
            &outcome_digests,
            &graph,
        );
        return Ok(exit_code);
    }
    if format == Format::Json {
        let mut outcomes = Vec::new();
        for value in &identical {
            let identifier = &value.identifier;
            let mut outcome = serde_json::json!({
                "identifier": identifier,
                "status": "IDENTICAL",
                "record_type": value.record_type.as_str(),
                "source_length": value.source_length,
                "path_length": value.graph_length
            });
            if let Some((start, end)) = value.source_range {
                outcome["source_range"] = serde_json::json!({
                    "start_0based": start,
                    "end_0based_exclusive": end
                });
            }
            if let Some(statistics) = outcome_statistics.get(identifier) {
                attach_statistics(&mut outcome, statistics);
            }
            attach_digests(&mut outcome, &outcome_digests[identifier]);
            outcomes.push(outcome);
        }
        for value in &divergent {
            let identifier = &value.identifier;
            let index = value.index;
            let reconstructed = &value.diagnostic;
            let mut outcome = serde_json::json!({
                "identifier": identifier,
                "status": "DIVERGENT",
                "record_type": reconstructed.record_type.as_str(),
                "source_length": value.source_length,
                "path_length": reconstructed.base_count,
                "source_position_1based": index + 1,
                "path_position_1based": index + 1,
                "divergence": {
                    "source_position_1based": index + 1,
                    "path_position_1based": index + 1
                }
            });
            if let Some((start, end)) = reconstructed.source_range {
                outcome["source_range"] = serde_json::json!({
                    "start_0based": start,
                    "end_0based_exclusive": end
                });
            }
            match reconstructed.divergence(index) {
                Divergence::RangeLength(range) => {
                    outcome["divergence"]["kind"] = "RANGE_LENGTH".into();
                    outcome["divergence"]["range_length"] = range.range_length.into();
                    outcome["divergence"]["walk_length"] = range.walk_length.into();
                }
                Divergence::Located(location) => {
                    outcome["divergence"]["segment"] = location.segment.clone().into();
                    outcome["divergence"]["orientation"] = location.orientation.to_string().into();
                    outcome["divergence"]["traversal_position_1based"] =
                        location.traversal_position.into();
                    outcome["divergence"]["segment_position_1based"] =
                        location.segment_position.into();
                }
                Divergence::GraphEnd => outcome["divergence"]["graph_state"] = "END".into(),
            }
            if let Some(statistics) = outcome_statistics.get(identifier) {
                attach_statistics(&mut outcome, statistics);
            }
            attach_digests(&mut outcome, &outcome_digests[identifier]);
            outcomes.push(outcome);
        }
        for identifier in &missing_paths {
            let mut outcome = serde_json::json!({
                "identifier": identifier,
                "status": "MISSING_PATH"
            });
            if let Some(statistics) = outcome_statistics.get(identifier) {
                attach_statistics(&mut outcome, statistics);
            }
            attach_digests(&mut outcome, &outcome_digests[identifier]);
            outcomes.push(outcome);
        }
        for identifier in &missing_sources {
            let mut outcome = serde_json::json!({
                "identifier": identifier,
                "status": "MISSING_SOURCE",
                "record_type": graph.record_type(identifier).as_str()
            });
            if let Some(statistics) = outcome_statistics.get(identifier) {
                attach_statistics(&mut outcome, statistics);
            }
            attach_digests(&mut outcome, &outcome_digests[identifier]);
            outcomes.push(outcome);
        }
        outcomes.sort_by(|left, right| {
            left["identifier"]
                .as_str()
                .cmp(&right["identifier"].as_str())
        });
        let mut document = serde_json::json!({
            "schema": "panpath-audit-report",
            "schema_version": env!("CARGO_PKG_VERSION"),
            "state": "completed",
            "summary": {
                "identical": identical.len(),
                "divergent": divergent.len(),
                "missing_path": missing_paths.len(),
                "missing_source": missing_sources.len()
            },
            "outcomes": outcomes,
            "errors": [],
            "provenance": {
                "tool_version": env!("CARGO_PKG_VERSION"),
                "source_path": (sources.len() == 1).then(|| sources[0].path.to_string_lossy()),
                "source_paths": sources.iter().map(|source| source.path.to_string_lossy()).collect::<Vec<_>>(),
                "sources": sources.iter().map(SourceSpec::provenance).collect::<Vec<_>>(),
                "graph_path": gfa_path.to_string_lossy(),
                "fasta_compression": (fasta_compressions.len() == 1).then(|| fasta_compressions[0]),
                "fasta_compressions": fasta_compressions,
                "gfa_version": graph.version,
                "gfa_compression": graph.compression,
                "path_records": path_record_count,
                "walk_records": walk_record_count,
                "record_types": match (path_record_count > 0, walk_record_count > 0) {
                    (true, true) => serde_json::json!(["P", "W"]),
                    (true, false) => serde_json::json!(["P"]),
                    (false, true) => serde_json::json!(["W"]),
                    (false, false) => serde_json::json!([]),
                },
                "unspecified_overlaps_as_blunt": graph.unspecified_overlaps_as_blunt,
                "topology_validated": false
                ,"threads": threads,
                "memory_mib": memory_telemetry.limit_bytes / (1024 * 1024),
                "peak_tracked_bytes": memory_telemetry.peak_tracked_bytes,
                "oversized_contigs": memory_telemetry.oversized_contigs
            }
        });
        if stats_mode != StatsMode::Basic {
            document["statistics"] = comprehensive_statistics_json(
                &outcome_statistics,
                &identical,
                &divergent,
                &missing_paths,
                &missing_sources,
                &sources,
                &graph,
            );
            document["provenance"]["statistics"] = serde_json::json!({
                "mode": "comprehensive",
                "alignment": "exact_unit_cost_levenshtein_wavefront",
                "alignment_max_cells": comprehensive(stats_mode).unwrap()
            });
        }
        println!(
            "{}",
            serde_json::to_string(&document).map_err(|error| error.to_string())?
        );
        return Ok(exit_code);
    }

    println!(
        "PROVENANCE topology_validated=false unspecified_overlaps_as_blunt={}",
        graph.unspecified_overlaps_as_blunt
    );
    println!(
        "MEMORY limit_mib={} peak_tracked_bytes={} oversized_contigs={}",
        memory_telemetry.limit_bytes / (1024 * 1024),
        memory_telemetry.peak_tracked_bytes,
        memory_telemetry.oversized_contigs
    );
    if stats_mode != StatsMode::Basic {
        let statistics = comprehensive_statistics_json(
            &outcome_statistics,
            &identical,
            &divergent,
            &missing_paths,
            &missing_sources,
            &sources,
            &graph,
        );
        println!(
            "SOURCE_BASES total={} matched={} substituted={} source_only={} not_embedded={} missing_path={} unaligned_divergent={}",
            statistics["source"]["total_bases"],
            statistics["source"]["bases"]["matched"],
            statistics["source"]["bases"]["substituted"],
            statistics["source"]["bases"]["source_only"],
            statistics["source"]["bases"]["not_embedded"],
            statistics["source"]["bases"]["missing_path"],
            statistics["source"]["bases"]["unaligned_divergent"]
        );
        println!(
            "GRAPH_BASES total={} matched={} substituted={} graph_only={} missing_source={} unaligned_divergent={}",
            statistics["graph"]["total_bases"],
            statistics["graph"]["bases"]["matched"],
            statistics["graph"]["bases"]["substituted"],
            statistics["graph"]["bases"]["graph_only"],
            statistics["graph"]["bases"]["missing_source"],
            statistics["graph"]["bases"]["unaligned_divergent"]
        );
        println!(
            "ALIGNMENT identity_percent={} edit_distance={} unavailable={}",
            statistics["alignment"]["identity_percent"],
            statistics["alignment"]["edit_distance"],
            statistics["alignment"]["unavailable_traversals"]
        );
        if let Some(items) = statistics["by_source"].as_array() {
            for item in items {
                println!(
                    "SOURCE index={} sequences={} bases={} exact_recovery_percent={}",
                    item["source_index"],
                    item["source"]["sequences"],
                    item["source"]["total_bases"],
                    item["source"]["exact_recovery_percent"]
                );
            }
        }
    }

    if !identical.is_empty() {
        println!("IDENTICAL {}", identical.len());
    }
    if !divergent.is_empty() {
        println!("DIVERGENT {}", divergent.len());
    }
    if !missing_paths.is_empty() {
        println!("MISSING_PATH {}", missing_paths.len());
    }
    if !missing_sources.is_empty() {
        println!("MISSING_SOURCE {}", missing_sources.len());
    }
    for value in divergent {
        let identifier = value.identifier;
        let index = value.index;
        let source_length = value.source_length;
        let reconstructed = value.diagnostic;
        let prefix = format!(
            "{identifier} source_position={} path_position={} source_length={source_length} path_length={}",
            index + 1,
            index + 1,
            reconstructed.base_count
        );
        match reconstructed.divergence(index) {
            Divergence::RangeLength(range) => println!(
                "{prefix} range_length={} walk_length={}",
                range.range_length, range.walk_length
            ),
            Divergence::Located(location) => println!(
                "{prefix} segment={}{} traversal_position={} segment_position={}",
                location.segment,
                location.orientation,
                location.traversal_position,
                location.segment_position
            ),
            Divergence::GraphEnd => println!("{prefix} graph=<END>"),
        }
    }
    for identifier in &missing_paths {
        println!("{identifier} MISSING_PATH");
    }
    for identifier in &missing_sources {
        println!("{identifier} MISSING_SOURCE");
    }
    if stats_mode != StatsMode::Basic {
        let identical_names: BTreeSet<_> = identical
            .iter()
            .map(|value| value.identifier.as_str())
            .collect();
        for (identifier, statistics) in &outcome_statistics {
            if !identical_names.contains(identifier.as_str()) {
                println!(
                    "{identifier} BASES source={} graph={} matched={} substituted={} source_only={} graph_only={} alignment={}",
                    statistics.source_bases,
                    statistics.graph_bases,
                    statistics.bases.matched,
                    statistics.bases.substituted,
                    statistics.bases.source_only,
                    statistics.bases.graph_only,
                    statistics.alignment_status.as_str()
                );
            }
        }
    }

    Ok(exit_code)
}

/// Parses CLI arguments and applies resource defaults.
fn parse_arguments(args: Vec<std::ffi::OsString>) -> Result<Config, String> {
    let mut format = Format::Human;
    let mut format_seen = false;
    let mut threads = None;
    let mut stats = None;
    let mut alignment_max_cells = None;
    let mut memory_mib = None;
    let mut temp_dir = None;
    let mut sources = Vec::new();
    let mut positional = Vec::new();
    let mut positional_only = false;
    let mut args = args.into_iter();
    while let Some(argument) = args.next() {
        if positional_only {
            positional.push(argument);
        } else if argument == "--" {
            positional_only = true;
        } else if argument == "--format" {
            if format_seen {
                return Err("format specified more than once".into());
            }
            let value = args.next().ok_or("missing format")?;
            format = match value.to_str() {
                Some("human") => Format::Human,
                Some("json") => Format::Json,
                Some("tsv") => Format::Tsv,
                _ => return Err("format must be human, json, or tsv".into()),
            };
            format_seen = true;
        } else if argument == "--threads" {
            if threads.is_some() {
                return Err("threads specified more than once".into());
            }
            let value = args.next().ok_or("missing thread count")?;
            let value = value
                .to_str()
                .and_then(|value| value.parse::<usize>().ok())
                .filter(|value| *value > 0)
                .ok_or("threads must be a positive integer")?;
            threads = Some(value);
        } else if argument == "--stats" {
            if stats.is_some() {
                return Err("stats specified more than once".into());
            }
            let value = args.next().ok_or("missing stats mode")?;
            if value != "comprehensive" {
                return Err("stats must be comprehensive".into());
            }
            stats = Some(());
        } else if argument == "--alignment-max-cells" {
            if alignment_max_cells.is_some() {
                return Err("alignment max cells specified more than once".into());
            }
            alignment_max_cells = Some(
                args.next()
                    .ok_or("missing alignment max cells")?
                    .to_str()
                    .and_then(|value| value.parse::<u64>().ok())
                    .filter(|value| *value > 0)
                    .ok_or("alignment max cells must be a positive integer")?,
            );
        } else if argument == "--memory-mib" {
            if memory_mib.is_some() {
                return Err("memory specified more than once".into());
            }
            memory_mib = Some(
                args.next()
                    .ok_or("missing memory limit")?
                    .to_str()
                    .and_then(|value| value.parse::<u64>().ok())
                    .filter(|value| *value > 0)
                    .ok_or("memory must be a positive integer in MiB")?,
            );
        } else if argument == "--temp-dir" {
            if temp_dir.is_some() {
                return Err("temporary directory specified more than once".into());
            }
            temp_dir = Some(std::path::PathBuf::from(
                args.next().ok_or("missing temporary directory")?,
            ));
        } else if argument == "--fasta" {
            sources.push(SourceSpec::exact(
                args.next().ok_or("missing --fasta file")?,
            ));
        } else if argument == "--fasta-pansn" {
            let sample = args.next().ok_or("missing --fasta-pansn sample")?;
            let haplotype = args.next().ok_or("missing --fasta-pansn haplotype")?;
            let path = args.next().ok_or("missing --fasta-pansn file")?;
            sources.push(SourceSpec {
                path,
                mode: SourceMode::Pansn {
                    sample: sample.to_string_lossy().into_owned(),
                    haplotype: haplotype.to_string_lossy().into_owned(),
                },
            });
        } else if argument == "--fasta-prefix" {
            let prefix = args.next().ok_or("missing --fasta-prefix prefix")?;
            let path = args.next().ok_or("missing --fasta-prefix file")?;
            sources.push(SourceSpec {
                path,
                mode: SourceMode::Prefix(prefix.to_string_lossy().into_owned()),
            });
        } else if argument.to_string_lossy().starts_with('-') {
            return Err(format!("unknown option {}", argument.to_string_lossy()));
        } else {
            positional.push(argument);
        }
    }
    let gfa = positional.pop().ok_or("missing GFA")?;
    sources.extend(positional.into_iter().map(SourceSpec::exact));
    if sources.is_empty() {
        return Err("missing source FASTA".into());
    }
    if stats.is_none() && alignment_max_cells.is_some() {
        return Err("alignment max cells requires --stats comprehensive".into());
    }
    let stats_mode = if stats.is_some() {
        StatsMode::Comprehensive {
            max_cells: alignment_max_cells.unwrap_or(DEFAULT_ALIGNMENT_MAX_CELLS),
        }
    } else {
        StatsMode::Basic
    };
    let threads = threads.unwrap_or_else(|| {
        thread::available_parallelism()
            .map(|value| value.get())
            .unwrap_or(1)
            .min(8)
    });
    Ok(Config {
        format,
        stats_mode,
        threads,
        memory_mib: memory_mib.unwrap_or(DEFAULT_MEMORY_MIB),
        temp_dir: temp_dir.unwrap_or_else(env::temp_dir),
        sources,
        gfa_path: gfa,
    })
}

impl SourceSpec {
    fn exact(path: std::ffi::OsString) -> Self {
        Self {
            path,
            mode: SourceMode::Exact,
        }
    }

    fn identifier(&self, raw: &str) -> String {
        match &self.mode {
            SourceMode::Exact => raw.into(),
            SourceMode::Pansn { sample, haplotype } => format!("{sample}#{haplotype}#{raw}"),
            SourceMode::Prefix(prefix) => format!("{prefix}{raw}"),
        }
    }

    fn provenance(&self) -> serde_json::Value {
        match &self.mode {
            SourceMode::Exact => {
                serde_json::json!({"path": self.path.to_string_lossy(), "mapping": "exact"})
            }
            SourceMode::Pansn { sample, haplotype } => {
                serde_json::json!({"path": self.path.to_string_lossy(), "mapping": "pansn", "sample": sample, "haplotype": haplotype})
            }
            SourceMode::Prefix(prefix) => {
                serde_json::json!({"path": self.path.to_string_lossy(), "mapping": "prefix", "prefix": prefix})
            }
        }
    }
}

/// Successful reconstruction with identity across all addressed ranges.
struct IdenticalOutcome {
    identifier: String,
    source_length: usize,
    graph_length: usize,
    record_type: RecordType,
    source_range: Option<(usize, usize)>,
}

/// Earliest observed divergence and its traversal context.
struct DivergentOutcome {
    identifier: String,
    index: usize,
    source_length: usize,
    diagnostic: TraversalDiagnostic,
}

/// Streaming comparison metadata used to report a traversal outcome.
struct TraversalDiagnostic {
    record_type: RecordType,
    source_range: Option<(usize, usize)>,
    base_count: usize,
    divergence_location: Option<(usize, Location)>,
    range_mismatch: Option<RangeMismatch>,
    not_embedded_bases: u64,
    graph_ambiguous_bases: u64,
}

impl TraversalDiagnostic {
    fn location(&self, position: usize) -> Option<&Location> {
        self.divergence_location
            .as_ref()
            .filter(|(index, _)| *index == position)
            .map(|(_, location)| location)
    }

    fn divergence(&self, position: usize) -> Divergence<'_> {
        if let Some(range) = self
            .range_mismatch
            .as_ref()
            .filter(|range| range.position == position)
        {
            Divergence::RangeLength(range)
        } else if let Some(location) = self.location(position) {
            Divergence::Located(location)
        } else {
            Divergence::GraphEnd
        }
    }
}

enum Divergence<'a> {
    RangeLength(&'a RangeMismatch),
    Located(&'a Location),
    GraphEnd,
}

#[derive(Clone, Default)]
/// Mutually exclusive base ledgers for one audit outcome.
struct BaseStatistics {
    matched: u64,
    substituted: u64,
    source_only: u64,
    graph_only: u64,
    not_embedded: u64,
    missing_path: u64,
    missing_source: u64,
    unaligned_divergent_source: u64,
    unaligned_divergent_graph: u64,
}

#[derive(Clone)]
/// Detailed measurements for one source/traversal correspondence outcome.
struct OutcomeStatistics {
    source_index: Option<usize>,
    source_bases: u64,
    graph_bases: u64,
    source_ambiguous_bases: u64,
    graph_ambiguous_bases: u64,
    ambiguous_columns: u64,
    length_delta: Option<i64>,
    alignment_status: AlignmentStatus,
    edit_distance: Option<u64>,
    bases: BaseStatistics,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum AlignmentStatus {
    NotApplicable,
    NotRequired,
    Complete,
    LimitExceeded,
    MemoryLimitExceeded,
}

impl AlignmentStatus {
    fn as_str(self) -> &'static str {
        match self {
            Self::NotApplicable => "not_applicable",
            Self::NotRequired => "not_required",
            Self::Complete => "complete",
            Self::LimitExceeded => "limit_exceeded",
            Self::MemoryLimitExceeded => "memory_limit_exceeded",
        }
    }
}

impl OutcomeStatistics {
    fn missing_path(source_index: usize, source: &str) -> Self {
        let source_bases = source.len() as u64;
        Self {
            source_index: Some(source_index),
            source_bases,
            graph_bases: 0,
            source_ambiguous_bases: ambiguous_count(source),
            graph_ambiguous_bases: 0,
            ambiguous_columns: 0,
            length_delta: None,
            alignment_status: AlignmentStatus::NotApplicable,
            edit_distance: None,
            bases: BaseStatistics {
                missing_path: source_bases,
                ..BaseStatistics::default()
            },
        }
    }

    fn missing_source(graph_bases: u64, graph_ambiguous_bases: u64) -> Self {
        Self {
            source_index: None,
            source_bases: 0,
            graph_bases,
            source_ambiguous_bases: 0,
            graph_ambiguous_bases,
            ambiguous_columns: 0,
            length_delta: None,
            alignment_status: AlignmentStatus::NotApplicable,
            edit_distance: None,
            bases: BaseStatistics {
                missing_source: graph_bases,
                ..BaseStatistics::default()
            },
        }
    }

    fn unaligned_divergent(
        source_index: usize,
        source_bases: u64,
        source_ambiguous_bases: u64,
        diagnostic: &TraversalDiagnostic,
    ) -> Result<Self, String> {
        let graph_bases = diagnostic.base_count as u64;
        let addressed_source = source_bases - diagnostic.not_embedded_bases;
        Ok(Self {
            source_index: Some(source_index),
            source_bases,
            graph_bases,
            source_ambiguous_bases,
            graph_ambiguous_bases: diagnostic.graph_ambiguous_bases,
            ambiguous_columns: 0,
            length_delta: Some(signed_delta(graph_bases, addressed_source)?),
            alignment_status: AlignmentStatus::MemoryLimitExceeded,
            edit_distance: None,
            bases: BaseStatistics {
                not_embedded: diagnostic.not_embedded_bases,
                unaligned_divergent_source: addressed_source,
                unaligned_divergent_graph: graph_bases,
                ..BaseStatistics::default()
            },
        })
    }
}

struct AuditedOutcome {
    outcome: AuditOutcome,
    statistics: Option<OutcomeStatistics>,
    digests: OutcomeDigests,
}

#[derive(Clone)]
struct OutcomeDigests {
    source: Option<SequenceDigest>,
    addressed_source: Option<SequenceDigest>,
    graph: Option<SequenceDigest>,
}

/// Source work item admitted under the shared memory budget.
struct SourceJob {
    source_index: usize,
    identifier: String,
    source: String,
    _permit: MemoryPermit,
}

/// Divergent traversal queued for serialized detailed alignment.
struct AlignmentJob {
    audited: AuditedOutcome,
    source_index: usize,
    identifier: String,
    source_length: usize,
    source_ambiguous_bases: u64,
    source: TemporarySequence,
}

/// Anonymous temporary storage for a sequence that cannot remain in memory.
struct TemporarySequence {
    file: fs::File,
}

#[derive(Default)]
/// Aggregates outcome statistics across a selected group.
struct StatisticsAccumulator {
    source_sequences: u64,
    graph_traversals: u64,
    source_bases: u64,
    graph_bases: u64,
    source_ambiguous_bases: u64,
    graph_ambiguous_bases: u64,
    ambiguous_columns: u64,
    bases: BaseStatistics,
    edit_distance: u64,
    aligned_columns: u64,
    completed_alignments: u64,
    unavailable_alignments: u64,
}

impl StatisticsAccumulator {
    fn add(&mut self, value: &OutcomeStatistics) {
        self.source_sequences += u64::from(value.source_index.is_some());
        self.graph_traversals += u64::from(value.graph_bases > 0 || value.source_index.is_none());
        self.source_bases += value.source_bases;
        self.graph_bases += value.graph_bases;
        self.source_ambiguous_bases += value.source_ambiguous_bases;
        self.graph_ambiguous_bases += value.graph_ambiguous_bases;
        self.ambiguous_columns += value.ambiguous_columns;
        self.bases.add(&value.bases);
        if let Some(distance) = value.edit_distance {
            self.edit_distance += distance;
            self.aligned_columns += value.bases.matched
                + value.bases.substituted
                + value.bases.source_only
                + value.bases.graph_only;
        }
        self.completed_alignments += u64::from(value.alignment_status == AlignmentStatus::Complete);
        self.unavailable_alignments += u64::from(matches!(
            value.alignment_status,
            AlignmentStatus::LimitExceeded | AlignmentStatus::MemoryLimitExceeded
        ));
    }
}

impl BaseStatistics {
    fn add(&mut self, other: &Self) {
        self.matched += other.matched;
        self.substituted += other.substituted;
        self.source_only += other.source_only;
        self.graph_only += other.graph_only;
        self.not_embedded += other.not_embedded;
        self.missing_path += other.missing_path;
        self.missing_source += other.missing_source;
        self.unaligned_divergent_source += other.unaligned_divergent_source;
        self.unaligned_divergent_graph += other.unaligned_divergent_graph;
    }
}

impl TemporarySequence {
    fn write(sequence: &[u8], temp_dir: &std::path::Path) -> Result<Self, String> {
        let mut file = tempfile::tempfile_in(temp_dir).map_err(|error| error.to_string())?;
        file.write_all(sequence)
            .map_err(|error| error.to_string())?;
        Ok(Self { file })
    }

    fn read(&mut self) -> Result<String, String> {
        self.file
            .seek(SeekFrom::Start(0))
            .map_err(|error| error.to_string())?;
        let mut sequence = String::new();
        self.file
            .read_to_string(&mut sequence)
            .map_err(|error| error.to_string())?;
        Ok(sequence)
    }
}

enum AuditOutcome {
    Identical(IdenticalOutcome),
    Divergent(DivergentOutcome),
    MissingPath(String),
}

/// Audits all streamed FASTA records through bounded worker queues.
fn audit_sources(
    sources: &[SourceSpec],
    source_lengths: &BTreeMap<String, usize>,
    graph: &Graph,
    threads: usize,
    stats_mode: StatsMode,
    memory_mib: u64,
    temp_dir: &std::path::Path,
) -> Result<(Vec<AuditedOutcome>, MemoryTelemetry), String> {
    let limit_bytes = memory_mib
        .checked_mul(1024 * 1024)
        .ok_or("memory limit overflow")?;
    let budget = MemoryBudget::new(limit_bytes);
    let mut readers = (0..threads)
        .map(|_| graph.reader())
        .collect::<Result<Vec<_>, _>>()?;
    let mut alignment_reader = graph.reader()?;
    let (work_sender, work_receiver) = mpsc::sync_channel::<SourceJob>(1);
    let work_receiver = Mutex::new(work_receiver);
    let (alignment_sender, alignment_receiver) = mpsc::sync_channel::<AlignmentJob>(1);
    let (result_sender, result_receiver) = mpsc::channel();
    let cancelled = AtomicBool::new(false);

    let (producer_result, worker_panicked, aligner_panicked) = thread::scope(|scope| {
        let aligner_result_sender = result_sender.clone();
        let aligner_cancelled = &cancelled;
        let aligner_budget = budget.clone();
        let aligner = scope.spawn(move || {
            for job in alignment_receiver {
                let result = panic::catch_unwind(AssertUnwindSafe(|| {
                    finish_alignment(&mut alignment_reader, &aligner_budget, stats_mode, job)
                }))
                .unwrap_or_else(|_| Err("alignment worker panicked".into()));
                if result.is_err() {
                    aligner_cancelled.store(true, Ordering::Relaxed);
                }
                if aligner_result_sender.send(result).is_err() {
                    break;
                }
            }
        });
        let mut handles = Vec::new();
        for mut reader in readers.drain(..) {
            let result_sender = result_sender.clone();
            let work_receiver = &work_receiver;
            let cancelled = &cancelled;
            let alignment_sender = alignment_sender.clone();
            handles.push(scope.spawn(move || {
                loop {
                    if cancelled.load(Ordering::Relaxed) {
                        break;
                    }
                    let work = match work_receiver.lock() {
                        Ok(receiver) => receiver.recv(),
                        Err(_) => {
                            cancelled.store(true, Ordering::Relaxed);
                            let _ = result_sender.send(Err("audit work queue poisoned".into()));
                            break;
                        }
                    };
                    let Ok(job) = work else {
                        break;
                    };
                    let SourceJob {
                        source_index,
                        identifier,
                        source,
                        _permit,
                    } = job;
                    let result = panic::catch_unwind(AssertUnwindSafe(|| {
                        audit_source(
                            graph,
                            &mut reader,
                            source_index,
                            identifier.clone(),
                            &source,
                            stats_mode,
                        )
                    }))
                    .unwrap_or_else(|_| Err("audit worker panicked".into()));
                    match result {
                        Ok(audited)
                            if comprehensive(stats_mode).is_some()
                                && matches!(audited.outcome, AuditOutcome::Divergent(_)) =>
                        {
                            let source_ambiguous_bases = ambiguous_count(&source);
                            let temporary = TemporarySequence::write(source.as_bytes(), temp_dir);
                            let source_length = source.len();
                            drop(source);
                            drop(_permit);
                            match temporary {
                                Ok(source) => {
                                    if alignment_sender
                                        .send(AlignmentJob {
                                            audited,
                                            source_index,
                                            identifier,
                                            source_length,
                                            source_ambiguous_bases,
                                            source,
                                        })
                                        .is_err()
                                    {
                                        break;
                                    }
                                }
                                Err(error) => {
                                    cancelled.store(true, Ordering::Relaxed);
                                    let _ = result_sender.send(Err(error));
                                    break;
                                }
                            }
                        }
                        result => {
                            if result.is_err() {
                                cancelled.store(true, Ordering::Relaxed);
                            }
                            if result_sender.send(result).is_err() {
                                break;
                            }
                        }
                    }
                }
            }));
        }
        drop(result_sender);
        drop(alignment_sender);
        let producer_result =
            stream_source_jobs(sources, source_lengths, &budget, &cancelled, |job| {
                work_sender
                    .send(job)
                    .map_err(|_| "audit work queue closed".into())
            });
        drop(work_sender);
        let mut worker_panicked = false;
        for handle in handles {
            worker_panicked |= handle.join().is_err();
        }
        let aligner_panicked = aligner.join().is_err();
        (producer_result, worker_panicked, aligner_panicked)
    });

    let mut outcomes = Vec::new();
    let mut worker_error = None;
    for result in result_receiver {
        match result {
            Ok(outcome) => outcomes.push(outcome),
            Err(error) if worker_error.is_none() => worker_error = Some(error),
            Err(_) => {}
        }
    }
    if worker_panicked || aligner_panicked {
        return Err("audit worker panicked".into());
    }
    if let Some(error) = worker_error {
        return Err(error);
    }
    producer_result?;
    Ok((outcomes, budget.telemetry()?))
}

/// Compares one source sequence with its indexed graph traversal.
fn audit_source(
    graph: &Graph,
    reader: &mut GraphReader<'_>,
    source_index: usize,
    identifier: String,
    source: &str,
    stats_mode: StatsMode,
) -> Result<AuditedOutcome, String> {
    let full_source_digest = sequence_digest(source.as_bytes());
    if !graph.path_names.contains(&identifier) {
        let statistics = comprehensive(stats_mode)
            .map(|_| OutcomeStatistics::missing_path(source_index, source));
        return Ok(AuditedOutcome {
            outcome: AuditOutcome::MissingPath(identifier),
            statistics,
            digests: OutcomeDigests {
                source: Some(full_source_digest),
                addressed_source: None,
                graph: None,
            },
        });
    }
    let comparison = reader.compare_traversal(&identifier, source.as_bytes())?;
    let digests = OutcomeDigests {
        source: Some(full_source_digest),
        addressed_source: Some(comparison.addressed_source_digest.clone()),
        graph: Some(comparison.graph_digest.clone()),
    };
    let divergence = comparison.divergence;
    let diagnostic = TraversalDiagnostic {
        record_type: comparison.record_type,
        source_range: comparison.source_range,
        base_count: comparison.base_count,
        divergence_location: comparison.location.as_ref().map(|location| {
            (
                divergence.unwrap_or(0),
                Location {
                    segment: location.segment.clone(),
                    orientation: location.orientation,
                    traversal_position: location.traversal_position,
                    segment_position: location.segment_position,
                },
            )
        }),
        range_mismatch: comparison
            .range_mismatch
            .as_ref()
            .map(|range| RangeMismatch {
                position: range.position,
                range_length: range.range_length,
                walk_length: range.walk_length,
            }),
        not_embedded_bases: comparison.not_embedded_bases,
        graph_ambiguous_bases: comparison.graph_ambiguous_bases,
    };
    if let Some(index) = divergence {
        Ok(AuditedOutcome {
            outcome: AuditOutcome::Divergent(DivergentOutcome {
                identifier,
                index,
                source_length: source.len(),
                diagnostic,
            }),
            statistics: None,
            digests,
        })
    } else {
        let statistics = comprehensive(stats_mode)
            .map(|_| OutcomeStatistics::identical(source_index, source, &comparison))
            .transpose()?;
        Ok(AuditedOutcome {
            outcome: AuditOutcome::Identical(IdenticalOutcome {
                identifier,
                source_length: source.len(),
                graph_length: comparison.base_count,
                record_type: comparison.record_type,
                source_range: comparison.source_range,
            }),
            statistics,
            digests,
        })
    }
}

/// Materializes and aligns a divergent traversal under memory and cell limits.
fn finish_alignment(
    reader: &mut GraphReader<'_>,
    budget: &MemoryBudget,
    stats_mode: StatsMode,
    mut job: AlignmentJob,
) -> Result<AuditedOutcome, String> {
    let diagnostic = match &job.audited.outcome {
        AuditOutcome::Divergent(value) => &value.diagnostic,
        _ => return Err("alignment job is not divergent".into()),
    };
    let tracked_bytes = (job.source_length as u64)
        .checked_add(diagnostic.base_count as u64)
        .ok_or("alignment memory overflow")?;
    let Some(_permit) = budget.acquire_within_limit(tracked_bytes)? else {
        job.audited.statistics = Some(OutcomeStatistics::unaligned_divergent(
            job.source_index,
            job.source_length as u64,
            job.source_ambiguous_bases,
            diagnostic,
        )?);
        return Ok(job.audited);
    };
    let source = job.source.read()?;
    job.audited.statistics = Some(OutcomeStatistics::aligned_divergent(
        reader,
        &job.identifier,
        job.source_index,
        &source,
        diagnostic.base_count,
        comprehensive(stats_mode).ok_or("missing comprehensive statistics mode")?,
    )?);
    Ok(job.audited)
}

fn comprehensive(mode: StatsMode) -> Option<u64> {
    match mode {
        StatsMode::Basic => None,
        StatsMode::Comprehensive { max_cells } => Some(max_cells),
    }
}

impl OutcomeStatistics {
    /// Computes exact statistics for a divergent traversal when budgets permit.
    fn aligned_divergent(
        reader: &mut GraphReader<'_>,
        identifier: &str,
        source_index: usize,
        source: &str,
        expected_graph_bases: usize,
        max_cells: u64,
    ) -> Result<Self, String> {
        let parts = reader.materialize_parts(identifier, expected_graph_bases)?;
        let mut pairs = Vec::new();
        let mut addressed = 0_u64;
        for part in parts {
            let slice = match (part.start, part.end) {
                (Some(start), Some(end)) => {
                    addressed += (end - start) as u64;
                    &source.as_bytes()[start..end]
                }
                _ => {
                    addressed = source.len() as u64;
                    source.as_bytes()
                }
            };
            pairs.push((slice, part.sequence));
        }
        let not_embedded = source.len() as u64 - addressed;

        let source_bases = source.len() as u64;
        let graph_bases = pairs.iter().map(|(_, graph)| graph.len() as u64).sum();
        let addressed_source = source_bases - not_embedded;
        let graph_ambiguous_bases = pairs
            .iter()
            .flat_map(|(_, graph)| graph.iter())
            .filter(|base| is_ambiguous(**base))
            .count() as u64;
        let mut bases = BaseStatistics {
            not_embedded,
            ..BaseStatistics::default()
        };
        let mut counts = EditCounts::default();
        let mut cells = 0_u64;
        let mut available = true;
        for (source_part, graph_part) in &pairs {
            let remaining = max_cells.saturating_sub(cells);
            match statistics::align(source_part, graph_part, remaining) {
                Alignment::Complete {
                    counts: part,
                    cells: used,
                } => {
                    cells += used;
                    counts.add(part);
                }
                Alignment::LimitExceeded => {
                    available = false;
                    break;
                }
            }
        }
        if available {
            bases.matched = counts.matched;
            bases.substituted = counts.substituted;
            bases.source_only = counts.source_only;
            bases.graph_only = counts.graph_only;
        } else {
            bases.unaligned_divergent_source = addressed_source;
            bases.unaligned_divergent_graph = graph_bases;
            counts = EditCounts::default();
        }
        Ok(Self {
            source_index: Some(source_index),
            source_bases,
            graph_bases,
            source_ambiguous_bases: ambiguous_count(source),
            graph_ambiguous_bases,
            ambiguous_columns: counts.ambiguous_columns,
            length_delta: Some(signed_delta(graph_bases, addressed_source)?),
            alignment_status: if available {
                AlignmentStatus::Complete
            } else {
                AlignmentStatus::LimitExceeded
            },
            edit_distance: available.then(|| counts.edit_distance()),
            bases,
        })
    }

    /// Constructs exact statistics for an identical addressed traversal.
    fn identical(
        source_index: usize,
        source: &str,
        comparison: &streaming_gfa::PathComparison,
    ) -> Result<Self, String> {
        let source_bases = source.len() as u64;
        let graph_bases = comparison.base_count as u64;
        let not_embedded = source_bases - graph_bases;
        let graph_ambiguous_bases = comparison.graph_ambiguous_bases;
        Ok(Self {
            source_index: Some(source_index),
            source_bases,
            graph_bases,
            source_ambiguous_bases: ambiguous_count(source),
            graph_ambiguous_bases,
            ambiguous_columns: graph_ambiguous_bases,
            length_delta: Some(0),
            alignment_status: AlignmentStatus::NotRequired,
            edit_distance: Some(0),
            bases: BaseStatistics {
                matched: graph_bases,
                not_embedded,
                ..BaseStatistics::default()
            },
        })
    }
}

fn ambiguous_count(sequence: &str) -> u64 {
    sequence.bytes().filter(|base| is_ambiguous(*base)).count() as u64
}

fn signed_delta(graph: u64, source: u64) -> Result<i64, String> {
    i64::try_from(graph as i128 - source as i128).map_err(|_| "length delta overflow".into())
}

fn percentage(numerator: u64, denominator: u64) -> serde_json::Value {
    if denominator == 0 {
        serde_json::Value::Null
    } else {
        serde_json::json!(
            ((numerator as f64 * 100_000_000.0 / denominator as f64).round()) / 1_000_000.0
        )
    }
}

fn base_statistics_json(value: &OutcomeStatistics) -> serde_json::Value {
    serde_json::json!({
        "source": {
            "matched": value.bases.matched,
            "substituted": value.bases.substituted,
            "source_only": value.bases.source_only,
            "not_embedded": value.bases.not_embedded,
            "missing_path": value.bases.missing_path,
            "unaligned_divergent": value.bases.unaligned_divergent_source
        },
        "graph": {
            "matched": value.bases.matched,
            "substituted": value.bases.substituted,
            "graph_only": value.bases.graph_only,
            "missing_source": value.bases.missing_source,
            "unaligned_divergent": value.bases.unaligned_divergent_graph
        }
    })
}

fn attach_statistics(outcome: &mut serde_json::Value, statistics: &OutcomeStatistics) {
    outcome["base_statistics"] = base_statistics_json(statistics);
    outcome["alignment_status"] = statistics.alignment_status.as_str().into();
    outcome["edit_distance"] = statistics.edit_distance.into();
    outcome["length_delta"] = statistics.length_delta.into();
    outcome["ambiguity"] = serde_json::json!({
        "source_bases": statistics.source_ambiguous_bases,
        "graph_bases": statistics.graph_ambiguous_bases,
        "aligned_columns": statistics.ambiguous_columns
    });
}

fn attach_digests(outcome: &mut serde_json::Value, digests: &OutcomeDigests) {
    fn value(digest: &Option<SequenceDigest>, algorithm: &str) -> serde_json::Value {
        digest
            .as_ref()
            .map_or(serde_json::Value::Null, |digest| match algorithm {
                "sha256" => digest.sha256.clone().into(),
                _ => digest.blake3.clone().into(),
            })
    }
    outcome["digests"] = serde_json::json!({
        "sha256": {
            "source": value(&digests.source, "sha256"),
            "addressed_source": value(&digests.addressed_source, "sha256"),
            "graph": value(&digests.graph, "sha256")
        },
        "blake3": {
            "source": value(&digests.source, "blake3"),
            "addressed_source": value(&digests.addressed_source, "blake3"),
            "graph": value(&digests.graph, "blake3")
        }
    });
}

#[derive(Default)]
struct TsvDivergence {
    kind: String,
    source_position: String,
    graph_position: String,
    segment: String,
    orientation: String,
    traversal_position: String,
    segment_position: String,
}

struct TsvDigests {
    source_sha256: String,
    addressed_source_sha256: String,
    graph_sha256: String,
    source_blake3: String,
    addressed_source_blake3: String,
    graph_blake3: String,
}

struct TsvStatistics {
    fields: [String; 11],
}

struct TsvRow {
    row_type: String,
    identifier: String,
    status: String,
    record_type: String,
    code: String,
    message: String,
    source_start: String,
    source_end: String,
    source_length: String,
    graph_length: String,
    divergence: TsvDivergence,
    digests: TsvDigests,
    statistics: TsvStatistics,
    topology_validated: String,
    unspecified_overlaps_as_blunt: String,
}

impl TsvDigests {
    fn from_outcome(value: &OutcomeDigests) -> Self {
        let digest = |value: &Option<SequenceDigest>, algorithm: &str| {
            value.as_ref().map_or_else(String::new, |value| {
                if algorithm == "sha256" {
                    value.sha256.clone()
                } else {
                    value.blake3.clone()
                }
            })
        };
        Self {
            source_sha256: digest(&value.source, "sha256"),
            addressed_source_sha256: digest(&value.addressed_source, "sha256"),
            graph_sha256: digest(&value.graph, "sha256"),
            source_blake3: digest(&value.source, "blake3"),
            addressed_source_blake3: digest(&value.addressed_source, "blake3"),
            graph_blake3: digest(&value.graph, "blake3"),
        }
    }

    fn empty() -> Self {
        Self::from_outcome(&OutcomeDigests {
            source: None,
            addressed_source: None,
            graph: None,
        })
    }
}

impl TsvStatistics {
    fn from_outcome(value: Option<&OutcomeStatistics>) -> Self {
        let fields = value.map_or_else(
            || std::array::from_fn(|_| String::new()),
            |value| {
                [
                    value.bases.matched.to_string(),
                    value.bases.substituted.to_string(),
                    value.bases.source_only.to_string(),
                    value.bases.graph_only.to_string(),
                    value.bases.not_embedded.to_string(),
                    value.bases.missing_path.to_string(),
                    value.bases.missing_source.to_string(),
                    value.bases.unaligned_divergent_source.to_string(),
                    value.bases.unaligned_divergent_graph.to_string(),
                    value
                        .edit_distance
                        .map_or_else(String::new, |value| value.to_string()),
                    value.alignment_status.as_str().into(),
                ]
            },
        );
        Self { fields }
    }
}

impl TsvRow {
    fn error(code: &str, message: &str) -> Self {
        Self {
            row_type: "error".into(),
            identifier: String::new(),
            status: String::new(),
            record_type: String::new(),
            code: code.into(),
            message: message.replace('\t', " "),
            source_start: String::new(),
            source_end: String::new(),
            source_length: String::new(),
            graph_length: String::new(),
            divergence: TsvDivergence::default(),
            digests: TsvDigests::empty(),
            statistics: TsvStatistics::from_outcome(None),
            topology_validated: String::new(),
            unspecified_overlaps_as_blunt: String::new(),
        }
    }

    fn print(self) {
        let mut fields = vec![
            self.row_type,
            self.identifier,
            self.status,
            self.record_type,
            self.code,
            self.message,
            self.source_start,
            self.source_end,
            self.source_length,
            self.graph_length,
            self.divergence.kind,
            self.divergence.source_position,
            self.divergence.graph_position,
            self.divergence.segment,
            self.divergence.orientation,
            self.divergence.traversal_position,
            self.divergence.segment_position,
            self.digests.source_sha256,
            self.digests.addressed_source_sha256,
            self.digests.graph_sha256,
            self.digests.source_blake3,
            self.digests.addressed_source_blake3,
            self.digests.graph_blake3,
        ];
        fields.extend(self.statistics.fields);
        fields.push(self.topology_validated);
        fields.push(self.unspecified_overlaps_as_blunt);
        debug_assert_eq!(fields.len(), TSV_HEADER.split('\t').count());
        println!("{}", fields.join("\t"));
    }
}

/// Emits one TSV row per correspondence outcome.
fn render_tsv(
    identical: &[IdenticalOutcome],
    divergent: &[DivergentOutcome],
    missing_paths: &[String],
    missing_sources: &[String],
    statistics: &BTreeMap<String, OutcomeStatistics>,
    digests: &BTreeMap<String, OutcomeDigests>,
    graph: &Graph,
) {
    println!("{TSV_HEADER}");
    let emit = |identifier: &str,
                status: &str,
                record_type: Option<RecordType>,
                source_range: Option<(usize, usize)>,
                source_length: Option<usize>,
                graph_length: Option<usize>,
                divergence: Option<(usize, &TraversalDiagnostic)>| {
        let mut location = TsvDivergence::default();
        if let Some((index, diagnostic)) = divergence {
            location.source_position = (index + 1).to_string();
            location.graph_position = (index + 1).to_string();
            match diagnostic.divergence(index) {
                Divergence::RangeLength(range) => {
                    location.kind = "RANGE_LENGTH".into();
                    location.segment =
                        format!("range={} walk={}", range.range_length, range.walk_length);
                }
                Divergence::Located(value) => {
                    location.kind = "BASE".into();
                    location.segment = value.segment.clone();
                    location.orientation = value.orientation.to_string();
                    location.traversal_position = value.traversal_position.to_string();
                    location.segment_position = value.segment_position.to_string();
                }
                Divergence::GraphEnd => location.kind = "END".into(),
            }
        }
        TsvRow {
            row_type: "outcome".into(),
            identifier: identifier.into(),
            status: status.into(),
            record_type: record_type.map_or_else(String::new, |value| value.as_str().into()),
            code: String::new(),
            message: String::new(),
            source_start: source_range.map_or_else(String::new, |value| value.0.to_string()),
            source_end: source_range.map_or_else(String::new, |value| value.1.to_string()),
            source_length: source_length.map_or_else(String::new, |value| value.to_string()),
            graph_length: graph_length.map_or_else(String::new, |value| value.to_string()),
            divergence: location,
            digests: TsvDigests::from_outcome(&digests[identifier]),
            statistics: TsvStatistics::from_outcome(statistics.get(identifier)),
            topology_validated: "false".into(),
            unspecified_overlaps_as_blunt: graph.unspecified_overlaps_as_blunt.to_string(),
        }
        .print();
    };
    for value in identical {
        emit(
            &value.identifier,
            "IDENTICAL",
            Some(value.record_type),
            value.source_range,
            Some(value.source_length),
            Some(value.graph_length),
            None,
        );
    }
    for value in divergent {
        emit(
            &value.identifier,
            "DIVERGENT",
            Some(value.diagnostic.record_type),
            value.diagnostic.source_range,
            Some(value.source_length),
            Some(value.diagnostic.base_count),
            Some((value.index, &value.diagnostic)),
        );
    }
    for identifier in missing_paths {
        emit(identifier, "MISSING_PATH", None, None, None, None, None);
    }
    for identifier in missing_sources {
        emit(
            identifier,
            "MISSING_SOURCE",
            Some(graph.record_type(identifier)),
            None,
            None,
            None,
            None,
        );
    }
}

fn accumulator_json(value: &StatisticsAccumulator) -> serde_json::Value {
    let source_categories = [
        ("matched", value.bases.matched),
        ("substituted", value.bases.substituted),
        ("source_only", value.bases.source_only),
        ("not_embedded", value.bases.not_embedded),
        ("missing_path", value.bases.missing_path),
        (
            "unaligned_divergent",
            value.bases.unaligned_divergent_source,
        ),
    ];
    let graph_categories = [
        ("matched", value.bases.matched),
        ("substituted", value.bases.substituted),
        ("graph_only", value.bases.graph_only),
        ("missing_source", value.bases.missing_source),
        ("unaligned_divergent", value.bases.unaligned_divergent_graph),
    ];
    let source_counts: serde_json::Map<_, _> = source_categories
        .iter()
        .map(|(name, count)| ((*name).into(), (*count).into()))
        .collect();
    let source_percentages: serde_json::Map<_, _> = source_categories
        .iter()
        .map(|(name, count)| ((*name).into(), percentage(*count, value.source_bases)))
        .collect();
    let graph_counts: serde_json::Map<_, _> = graph_categories
        .iter()
        .map(|(name, count)| ((*name).into(), (*count).into()))
        .collect();
    let graph_percentages: serde_json::Map<_, _> = graph_categories
        .iter()
        .map(|(name, count)| ((*name).into(), percentage(*count, value.graph_bases)))
        .collect();
    let embedded = value.source_bases - value.bases.not_embedded - value.bases.missing_path;
    let sourced_graph = value.graph_bases - value.bases.missing_source;
    serde_json::json!({
        "source": {
            "sequences": value.source_sequences,
            "total_bases": value.source_bases,
            "bases": source_counts,
            "percentages": source_percentages,
            "embedded_coverage_percent": percentage(embedded, value.source_bases),
            "exact_recovery_percent": percentage(value.bases.matched, value.source_bases)
        },
        "graph": {
            "traversals": value.graph_traversals,
            "total_bases": value.graph_bases,
            "bases": graph_counts,
            "percentages": graph_percentages,
            "source_coverage_percent": percentage(sourced_graph, value.graph_bases)
        },
        "alignment": {
            "identity_percent": percentage(value.bases.matched, value.aligned_columns),
            "edit_distance": value.edit_distance,
            "completed_traversals": value.completed_alignments,
            "unavailable_traversals": value.unavailable_alignments,
            "aligned_columns": value.aligned_columns
        },
        "ambiguity": {
            "source_bases": value.source_ambiguous_bases,
            "graph_bases": value.graph_ambiguous_bases,
            "aligned_columns": value.ambiguous_columns
        }
    })
}

/// Builds aggregate comprehensive statistics and source-level breakdowns.
fn comprehensive_statistics_json(
    statistics: &BTreeMap<String, OutcomeStatistics>,
    identical: &[IdenticalOutcome],
    divergent: &[DivergentOutcome],
    missing_paths: &[String],
    missing_sources: &[String],
    sources: &[SourceSpec],
    graph: &Graph,
) -> serde_json::Value {
    let mut statuses = HashMap::new();
    for value in identical {
        statuses.insert(value.identifier.as_str(), "identical");
    }
    for value in divergent {
        statuses.insert(value.identifier.as_str(), "divergent");
    }
    for identifier in missing_paths {
        statuses.insert(identifier.as_str(), "missing_path");
    }
    for identifier in missing_sources {
        statuses.insert(identifier.as_str(), "missing_source");
    }
    let mut global = StatisticsAccumulator::default();
    let mut by_source: Vec<_> = (0..sources.len())
        .map(|_| StatisticsAccumulator::default())
        .collect();
    let mut by_record_type: BTreeMap<RecordType, StatisticsAccumulator> = BTreeMap::new();
    let mut global_statuses: BTreeMap<&str, u64> = BTreeMap::new();
    let mut source_statuses: Vec<BTreeMap<&str, u64>> =
        (0..sources.len()).map(|_| BTreeMap::new()).collect();
    let mut record_statuses: BTreeMap<RecordType, BTreeMap<&str, u64>> = BTreeMap::new();
    for (identifier, value) in statistics {
        global.add(value);
        let status = statuses
            .get(identifier.as_str())
            .copied()
            .unwrap_or("unknown");
        *global_statuses.entry(status).or_default() += 1;
        if let Some(index) = value.source_index {
            by_source[index].add(value);
            *source_statuses[index].entry(status).or_default() += 1;
        }
        if graph.path_names.contains(identifier) {
            let record_type = graph.record_type(identifier);
            by_record_type.entry(record_type).or_default().add(value);
            *record_statuses
                .entry(record_type)
                .or_default()
                .entry(status)
                .or_default() += 1;
        }
    }
    let source_total = global.source_sequences;
    let graph_total =
        identical.len() as u64 + divergent.len() as u64 + missing_sources.len() as u64;
    let mut document = accumulator_json(&global);
    document["source"]["statuses"] = status_json(
        &global_statuses,
        source_total,
        &["identical", "divergent", "missing_path"],
    );
    document["graph"]["statuses"] = status_json(
        &global_statuses,
        graph_total,
        &["identical", "divergent", "missing_source"],
    );
    document["by_source"] = serde_json::Value::Array(
        by_source
            .iter()
            .enumerate()
            .map(|(index, value)| {
                let mut item = accumulator_json(value);
                item["source_index"] = index.into();
                item["path"] = sources[index].path.to_string_lossy().into_owned().into();
                item["source"]["statuses"] = status_json(
                    &source_statuses[index],
                    value.source_sequences,
                    &["identical", "divergent", "missing_path"],
                );
                item
            })
            .collect(),
    );
    document["by_record_type"] = serde_json::Value::Object(
        by_record_type
            .iter()
            .map(|(record_type, value)| {
                let mut item = accumulator_json(value);
                item["graph"]["statuses"] = status_json(
                    &record_statuses[record_type],
                    value.graph_traversals,
                    &["identical", "divergent", "missing_source"],
                );
                (record_type.as_str().into(), item)
            })
            .collect(),
    );
    document
}

fn status_json(
    counts: &BTreeMap<&str, u64>,
    denominator: u64,
    names: &[&str],
) -> serde_json::Value {
    serde_json::Value::Object(
        names
            .iter()
            .map(|name| {
                let count = counts.get(name).copied().unwrap_or(0);
                (
                    (*name).into(),
                    serde_json::json!({"count": count, "percent": percentage(count, denominator)}),
                )
            })
            .collect(),
    )
}

/// Indexes identifiers and lengths while collecting all detectable FASTA errors.
fn index_sources(
    sources: &[SourceSpec],
) -> (BTreeMap<String, usize>, Vec<&'static str>, Vec<String>) {
    let mut index = BTreeMap::new();
    let mut compressions = Vec::new();
    let mut errors = Vec::new();
    for source in sources {
        let path = std::path::Path::new(&source.path);
        match fasta_compression(path) {
            Ok(compression) => compressions.push(compression),
            Err(error) => {
                errors.push(format!("FASTA_READ path={} error={error}", path.display()));
                continue;
            }
        }
        let result = stream_fasta(path, |identifier, sequence| {
            let identifier = source.identifier(identifier);
            if let Some((position, symbol)) = sequence
                .bytes()
                .enumerate()
                .find(|(_, symbol)| !is_iupac(*symbol))
            {
                errors.push(format!(
                    "INVALID_SOURCE_NUCLEOTIDE path={} identifier={identifier} position={} symbol={}",
                    path.display(),
                    position + 1,
                    char::from(symbol)
                ));
            }
            match index.entry(identifier) {
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(sequence.len());
                }
                std::collections::btree_map::Entry::Occupied(entry) => {
                    errors.push(format!(
                        "DUPLICATE_SEQUENCE_IDENTIFIER identifier={}",
                        entry.key()
                    ));
                }
            }
            Ok(())
        });
        if let Err(error) = result {
            errors.push(format!(
                "INVALID_FASTA path={} error={error}",
                path.display()
            ));
        }
    }
    (index, compressions, errors)
}

/// Streams source records into bounded worker jobs without retaining the FASTA set.
fn stream_source_jobs(
    sources: &[SourceSpec],
    source_lengths: &BTreeMap<String, usize>,
    budget: &MemoryBudget,
    cancelled: &AtomicBool,
    mut visit: impl FnMut(SourceJob) -> Result<(), String>,
) -> Result<(), String> {
    for (source_index, source) in sources.iter().enumerate() {
        let mut reader = open_fasta(std::path::Path::new(&source.path))?;
        let mut line = String::new();
        let mut current: Option<(String, String, MemoryPermit)> = None;
        while reader
            .read_line(&mut line)
            .map_err(|error| error.to_string())?
            != 0
        {
            let value = line.trim_end_matches(['\n', '\r']);
            if let Some(header) = value.strip_prefix('>') {
                if let Some((identifier, sequence, permit)) = current.take() {
                    visit(SourceJob {
                        source_index,
                        identifier,
                        source: sequence,
                        _permit: permit,
                    })?;
                }
                if cancelled.load(Ordering::Relaxed) {
                    return Err("audit worker stopped".into());
                }
                let raw = header
                    .split_whitespace()
                    .next()
                    .ok_or("invalid FASTA header")?;
                let identifier = source.identifier(raw);
                let length = *source_lengths
                    .get(&identifier)
                    .ok_or("source index changed between passes")?;
                let permit = budget.acquire(length as u64)?;
                current = Some((identifier, String::with_capacity(length), permit));
            } else if let Some((_, sequence, _)) = &mut current {
                sequence.push_str(value);
            }
            line.clear();
        }
        let (identifier, sequence, permit) = current.ok_or("empty FASTA")?;
        visit(SourceJob {
            source_index,
            identifier,
            source: sequence,
            _permit: permit,
        })?;
    }
    Ok(())
}

/// Parses FASTA records sequentially and invokes a callback for each sequence.
fn stream_fasta(
    path: &std::path::Path,
    mut visit: impl FnMut(&str, &str) -> Result<(), String>,
) -> Result<(), String> {
    let mut reader = open_fasta(path)?;
    let mut line = String::new();
    let mut name = None;
    let mut sequence = String::new();
    while reader
        .read_line(&mut line)
        .map_err(|error| error.to_string())?
        != 0
    {
        let value = line.trim_end_matches(['\n', '\r']);
        if let Some(header) = value.strip_prefix('>') {
            if let Some(previous) = name.replace(
                header
                    .split_whitespace()
                    .next()
                    .ok_or("invalid FASTA header")?
                    .to_owned(),
            ) {
                visit(&previous, &sequence)?;
                sequence.clear();
            }
        } else {
            sequence.push_str(value);
        }
        line.clear();
    }
    let name = name.ok_or("empty FASTA")?;
    visit(&name, &sequence)
}

fn open_fasta(path: &std::path::Path) -> Result<Box<dyn BufRead>, String> {
    input::open(path).map(|(reader, _)| reader)
}

fn fasta_compression(path: &std::path::Path) -> Result<&'static str, String> {
    input::open(path).map(|(_, compression)| compression.as_str())
}

struct Location {
    segment: String,
    orientation: Orientation,
    traversal_position: usize,
    segment_position: usize,
}

struct RangeMismatch {
    position: usize,
    range_length: usize,
    walk_length: usize,
}
