# Usage guide

`panpath-audit` compares source FASTA sequences with the named `P` paths and `W`
walks embedded in a GFA file. It accepts files, not standard input, because FASTA
inputs are read twice and the GFA is indexed in temporary storage.

## Command line

```text
panpath-audit [OPTIONS] [FASTA ...] GFA
```

The final positional path is the GFA. Any earlier positional paths are FASTA
sources using exact identifier matching:

```bash
panpath-audit assembly-1.fa assembly-2.fa graph.gfa
```

Explicit source options may be mixed with positional FASTAs:

```text
Input mapping:
  --fasta FILE                    Exact identifiers
  --fasta-pansn SAMPLE HAP FILE   SAMPLE#HAP#identifier
  --fasta-prefix PREFIX FILE      PREFIXidentifier

Output:
  --format human|json|tsv         Report format (default: human)
  --stats comprehensive           Add exact edit-aligned base statistics

Resources:
  --threads N                     Worker threads (default: available CPUs, max 8)
  --memory-mib N                  Tracked sequence memory (default: 1024)
  --alignment-max-cells N         Alignment work limit (default: 10000000)
  --temp-dir DIR                  Scratch directory (default: OS temporary directory)

General:
  -h, --help                      Print the option summary
  -V, --version                   Print the package version
  --                              Treat all remaining arguments as paths
```

Every numeric resource option must be a positive integer. An explicit
`--alignment-max-cells` requires `--stats comprehensive`. Options that take one
value may not be repeated, except for the three FASTA source options, which add
sources to the same audit.

Use `--` before positional filenames that begin with `-`:

```bash
panpath-audit -- -source.fa -graph.gfa
```

## FASTA inputs and identifier mapping

FASTA and GFA compression is detected from gzip content, not the filename, so
`.gz` is conventional but not required. Sequences may contain the DNA IUPAC
symbols `A C G T R Y S W K M B D H V N` in either letter case. Other symbols are
input errors.

The raw FASTA identifier is the first whitespace-delimited token after `>`.
Descriptions after that token are ignored. Identifier matching is case-sensitive,
whereas sequence comparison normalizes nucleotide symbols to uppercase.

Choose a mapping independently for each FASTA file:

| Input | FASTA header | Resulting identifier |
|---|---|---|
| `--fasta ref.fa` or positional `ref.fa` | `>chr1 description` | `chr1` |
| `--fasta-pansn CHM13 1 ref.fa` | `>chr1 description` | `CHM13#1#chr1` |
| `--fasta-prefix GRCh38. ref.fa` | `>chr1 description` | `GRCh38.chr1` |

For example, two FASTAs that both contain `>chr1` can be disambiguated as
different source assemblies:

```bash
panpath-audit \
  --fasta-pansn CHM13 1 chm13.fa \
  --fasta-prefix GRCh38. grch38.fa \
  graph.gfa
```

All mapped identifiers must be unique across the complete source set. A duplicate
is rejected even if the sequences are identical.

## Supported GFA contract

The audit reads these GFA1-family records:

- `S` records provide segment names and inline nucleotide sequences. A sequence
  value of `*`, an invalid IUPAC symbol, or a duplicate segment name is rejected.
- `P` records define named traversals. Steps may use forward (`+`) or reverse
  (`-`) orientation. The overlap field must be `*`, or contain one `0M` value for
  every segment boundary; other CIGAR overlaps are unsupported.
- `W` records define a traversal name as `sample#haplotype#sequence`. The
  haplotype must be numeric. A range may be `* *` or a valid half-open
  `[start, end)` interval within the paired source sequence.

Multiple `W` records for the same identifier are combined in coordinate order
when their explicit ranges do not overlap. Gaps between ranges are allowed: those
source bases are classified as `not_embedded` and are not treated as divergence.
The reconstructed length of each walk must equal its declared range length. An
unspecified range cannot be combined with another walk for the same identifier.

The entire input is preflighted before sequence auditing. Unresolved segments,
duplicate paths, a `P` and `W` claiming the same identifier, overlapping or
ambiguous walk ranges, malformed traversal syntax, and unsupported overlaps are
reported as errors rather than partially audited.

Records unrelated to reconstruction are ignored. In particular, panpath-audit
does not validate links or graph topology. Completed reports record
`topology_validated=false`; use a topology validator separately when that property
matters.

## Correspondence and outcomes

After identifier mapping, the source and graph name sets are compared:

| Outcome | Meaning |
|---|---|
| `IDENTICAL` | Source and reconstructed traversal are equal across every addressed range. |
| `DIVERGENT` | Symbols differ, or one addressed sequence ends before the other. |
| `MISSING_PATH` | A mapped source identifier has no graph traversal. |
| `MISSING_SOURCE` | A graph traversal has no mapped source sequence. |

Sequence equality is case-insensitive but symbol-exact. For example, `N` matches
`n`, but `N` does not match `A`. Reverse steps use the complete DNA IUPAC
complement mapping.

Positions in diagnostics are one-based. Walk ranges in JSON and TSV retain the
GFA convention: `source_start_0based` is inclusive and
`source_end_0based_exclusive` is exclusive. A base divergence may include the
source position, reconstructed traversal position, segment name and orientation,
position within the traversal, and position within the segment. Exhaustion at the
end of the graph has no segment coordinate.

## Report formats

### Human

The default report is intended for terminals and CI logs. It prints provenance
and memory telemetry, outcome counts, and one line for every non-identical
correspondence. With comprehensive statistics it also prints aggregate source,
graph, and alignment ledgers.

```bash
panpath-audit source.fa graph.gfa
```

Human output is a concise summary, not a stable interchange format. Use JSON or
TSV for automation.

### JSON

JSON output is one document on standard output, including when preflight fails:

```bash
panpath-audit --format json source.fa graph.gfa > report.json
```

The document identifies itself with `schema: "panpath-audit-report"` and uses
the panpath-audit package version as `schema_version`. A successful preflight has
`state: "completed"`, a summary, sorted per-identifier outcomes, an empty error
list, and provenance. A failed preflight has `state: "invalid"`, no outcomes, and
structured errors with `category`, `code`, and `message`.

Each completed outcome contains the fields meaningful for its status. Digests use
`null` when the corresponding sequence does not exist. Comprehensive mode adds
per-outcome base statistics and an aggregate `statistics` object; on an invalid
preflight that object is `null`.

### TSV

TSV output always begins with the same header. Completed audits emit one
`row_type=outcome` row per identifier. Invalid preflight emits `row_type=error`
rows using the same columns, leaving non-applicable values empty.

```bash
panpath-audit --format tsv source.fa graph.gfa > report.tsv
```

The columns cover status and record type, diagnostics, source ranges, sequence
lengths, divergence coordinates, digests, comprehensive base counts, alignment
status, and the topology/overlap provenance flags. Read columns by header name
rather than fixed numeric position.

## Digests and comprehensive statistics

JSON and TSV always report SHA-256 and BLAKE3 digests where the relevant sequence
exists:

- `source` covers the complete mapped FASTA sequence.
- `addressed_source` covers only the source ranges represented by the traversal.
- `graph` covers the reconstructed graph traversal.

Digest input is normalized to uppercase. Consequently, files that differ only in
nucleotide letter case produce the same sequence digests.

`--stats comprehensive` adds an exact unit-cost Levenshtein alignment for
divergent traversals when resource limits permit:

```bash
panpath-audit --format json --stats comprehensive source.fa graph.gfa
```

Source bases are partitioned into `matched`, `substituted`, `source_only`,
`not_embedded`, `missing_path`, and `unaligned_divergent`. Graph bases are
partitioned separately into `matched`, `substituted`, `graph_only`,
`missing_source`, and `unaligned_divergent`. Percentages use the corresponding
source or graph denominator.

Alignment work is bounded by `--alignment-max-cells` (default 10,000,000). If the
limit is exhausted, the audit remains `DIVERGENT`, `alignment_status` is
`limit_exceeded`, and addressed bases are assigned to the unaligned-divergent
buckets. If the source and reconstructed traversal cannot simultaneously fit the
memory budget, the status is `memory_limit_exceeded`. In both cases edit distance
is unavailable, but first-divergence diagnostics and digests remain valid.

## Resource behavior

`--threads` defaults to the available CPU count capped at eight. Set
`--threads 1` for sequential comparison or increase it explicitly when the
storage and memory budget can support more concurrent work.

`--memory-mib` defaults to 1024 and limits tracked sequence bytes admitted to
concurrent comparison and alignment. One source contig larger than the limit is
allowed to run alone so the audit can make progress; reports count such contigs as
`oversized_contigs`. Parser metadata, temporary-file metadata, allocator overhead,
compression buffers, and operating-system caches are not included in the tracked
peak.

FASTA inputs are streamed twice. Graph segments and traversals are stored in
anonymous scratch files rather than retained as one in-memory graph. Scratch uses
the operating-system temporary directory unless `--temp-dir DIR` selects another
writable filesystem. Anonymous stores are reclaimed when the process exits,
including after forced termination.

## Exit codes

| Code | Meaning |
|---|---|
| `0` | All paired sequences are identical and correspondence is complete. |
| `1` | At least one paired sequence is divergent. |
| `2` | Correspondence is missing or ambiguous. |
| `3` | The GFA is invalid or unsupported. |
| `4` | Source input, usage, scratch storage, or another operational failure. |

When several conditions coexist, the highest code takes precedence. Inspect the
exit status even when consuming JSON or TSV.

## Troubleshooting

- `DUPLICATE_SEQUENCE_IDENTIFIER`: choose PanSN or prefix mappings so every
  source identifier is unique.
- `MISSING_PATH` or `MISSING_SOURCE`: compare mapped FASTA identifiers with `P`
  names and the derived `sample#haplotype#sequence` names of `W` records.
- `UNSUPPORTED_OVERLAP`: export blunt paths with `*` or explicit `0M` overlaps;
  panpath-audit does not trim overlap CIGARs.
- `SEQUENCELESS_SEGMENT`: embed segment sequences on `S` records before auditing.
- `W_RANGE_OUT_OF_BOUNDS`, `W_RANGE_OVERLAP`, or `AMBIGUOUS_W_RANGES`: correct
  coordinate ranges or split the intended sources into unambiguous identifiers.
- `alignment_status=limit_exceeded`: increase `--alignment-max-cells` if exact
  edit statistics are required.
- `alignment_status=memory_limit_exceeded`: increase `--memory-mib`; the primary
  audit result itself remains valid.
- `TEMP_STORE`: select a writable scratch filesystem with sufficient space using
  `--temp-dir`.
