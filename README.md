<p align="center">
  <img src="https://raw.githubusercontent.com/gkanogiannis/panpath-audit/main/assets/panpath-audit-icon.png" alt="panpath-audit" width="140">
</p>

<h1 align="center">panpath-audit</h1>

<p align="center">
  <em>Verify that a pangenome graph preserved your sequences.</em>
</p>

<p align="center">
  <a href="https://crates.io/crates/panpath-audit"><img src="https://img.shields.io/crates/v/panpath-audit.svg" alt="crates.io version"></a>
  <a href="https://github.com/gkanogiannis/panpath-audit/actions/workflows/ci.yml"><img src="https://github.com/gkanogiannis/panpath-audit/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://github.com/gkanogiannis/panpath-audit/blob/main/LICENSE"><img src="https://img.shields.io/badge/License-Apache_2.0-blue.svg" alt="License: Apache-2.0"></a>
  <img src="https://img.shields.io/badge/rust-1.89%2B-orange.svg" alt="Rust 1.89+">
</p>

`panpath-audit` checks that every genome or haplotype traversal embedded in a
pangenome graph reconstructs **exactly** the source FASTA sequence it is supposed
to represent, and tells you precisely where it does not.

## Why

Pangenome variation graphs (GFA) represent input genomes as paths through
sequence-labelled segments. Graph validators such as `odgi validate` check that a
graph is internally consistent, but a structurally valid graph can still encode
the *wrong* sequence. Reconstruction steps, format conversions, or accidentally
redistributing a different FASTA can silently damage a path while leaving the
topology valid.

panpath-audit answers one binary question:

> Did this graph preserve my genomes?

It reconstructs each embedded traversal, compares it base-for-base against the
source FASTA, and reports the outcome with localized diagnostics and CI-grade exit
codes. It does **not** validate graph topology; that is a separate concern. See
[ADR-0003](https://github.com/gkanogiannis/panpath-audit/blob/main/docs/adr/0003-do-not-validate-graph-topology.md).

## Features

- **GFA1 / GFA1.1**: inline `S` segments, `P` paths, and coordinate-aware `W` walks
- **Plain or gzip** FASTA and GFA input, detected from content
- **Exact IUPAC reverse-complement** for reverse-oriented segments (R↔Y, K↔M, B↔V, D↔H, …)
- **Flexible identifier mapping**: exact, PanSN (`sample#haplotype#contig`), or a user-supplied prefix
- **Blunt-path contract**: `*` and explicit `0M` overlaps are treated as oriented concatenation
- **Streaming with bounded memory** (`--memory-mib`) and **parallel** comparison (`--threads`)
- **Human, JSON, and TSV reports** with stable outcome names and CI-friendly exit codes
- **SHA-256 and BLAKE3 digests** over source, addressed source ranges, and reconstructed traversals
- **Comprehensive statistics** with separate source/graph base ledgers and bounded exact edit-distance alignment

## Install

The current release is available from
[crates.io](https://crates.io/crates/panpath-audit) and requires Rust 1.89 or
newer:

```bash
cargo install panpath-audit
```

To build a repository checkout instead:

```bash
git clone https://github.com/gkanogiannis/panpath-audit.git
cd panpath-audit
cargo build --release              # binary: target/release/panpath-audit
# or: cargo install --path .
```

Confirm the installed version with `panpath-audit --version`.

## Quick start

Use positional arguments when FASTA identifiers already match graph traversal
names. The final positional path is always the GFA:

```bash
panpath-audit source.fa graph.gfa
```

For a FASTA record named `chr01` and a graph traversal named
`SAMPLE#1#chr01`, apply a PanSN mapping:

```bash
panpath-audit --fasta-pansn SAMPLE 1 source.fa graph.gfa
```

For automation, request JSON and optionally comprehensive base statistics:

```bash
panpath-audit --format json --stats comprehensive source.fa.gz graph.gfa.gz \
  > audit.json
```

FASTA identifiers are the first whitespace-delimited token after `>`. Matching
is case-sensitive; nucleotide symbols are compared case-insensitively. Multiple
FASTA inputs form one source set, and every resulting identifier must be unique.
See the [usage guide](https://github.com/gkanogiannis/panpath-audit/blob/main/docs/usage.md)
for mapping examples and the complete option reference.

## Example

```text
$ panpath-audit assemblies.fa.gz pangenome.gfa
PROVENANCE topology_validated=false unspecified_overlaps_as_blunt=true
MEMORY limit_mib=1024 peak_tracked_bytes=558010206 oversized_contigs=0
IDENTICAL 22
DIVERGENT 1
SampleA#1#chr02 source_position=1 path_position=1 source_length=55739602 path_length=55739602 segment=5+ traversal_position=1 segment_position=1
```

A divergence is localized to both the source coordinate and the graph segment
(including orientation and within-segment position), rather than merely reporting
that two sequences differ.

Use the process status as a CI gate:

```bash
panpath-audit assemblies.fa.gz pangenome.gfa
echo "exit: $?"   # 0 means complete correspondence and sequence identity
```

## Audit contract

| Outcome | Meaning |
|---|---|
| **Identical** | The paired source and traversal have sequence identity across all addressed ranges after case normalization. Ambiguity codes remain distinct and must match exactly. |
| **Divergent** | The paired source and traversal differ, or one is exhausted, at the first reported position. |
| **Missing path** | A source sequence has no paired embedded traversal. |
| **Missing source** | An embedded traversal has no paired source sequence. |

`P` records must use `*` or one `0M` overlap per segment boundary. `W` records
are paired by the PanSN name `sample#haplotype#sequence`; multiple non-overlapping
ranges for that name are audited in coordinate order. Bases in gaps between those
ranges are reported as not embedded, not divergent. Segment sequences must be
present inline on `S` records.

The tool intentionally does not validate graph topology, external segment
sequence references, or non-blunt path overlaps. Invalid or ambiguous input is
rejected before auditing starts.

## Exit codes

| Code | Meaning |
|---|---|
| `0` | Complete correspondence; all paired sequences are identical |
| `1` | At least one paired sequence is divergent |
| `2` | Missing or ambiguous source/traversal correspondence |
| `3` | Invalid or unsupported GFA |
| `4` | Invalid source, command usage, scratch storage, or another operational failure |

When conditions coexist, the highest code takes precedence. JSON and TSV remain
valid machine-readable reports for detected preflight failures.

## How it works

<p align="center">
  <img src="https://raw.githubusercontent.com/gkanogiannis/panpath-audit/main/docs/assets/pangenome-gfa-terms.png" alt="Pangenome / GFA terms" width="620">
</p>

FASTA identifiers are indexed, then each graph traversal is streamed, its oriented
segment sequences reconstructed, and compared with the paired source sequence.
Source FASTAs are read twice so sequence memory is bounded by concurrent contigs;
graph segments and traversals stream through anonymous temporary stores.

For operational details, read the
[usage guide](https://github.com/gkanogiannis/panpath-audit/blob/main/docs/usage.md).
The repository also records its
[design decisions](https://github.com/gkanogiannis/panpath-audit/tree/main/docs/adr)
and [terminology](https://github.com/gkanogiannis/panpath-audit/blob/main/CONTEXT.md).

## Development

```bash
cargo test
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo build --release
```

## Citation

If you use panpath-audit in your work, please cite it using
[`CITATION.cff`](https://github.com/gkanogiannis/panpath-audit/blob/main/CITATION.cff).
GitHub also exposes this metadata through its “Cite this repository” button.

## License

Apache-2.0. See
[`LICENSE`](https://github.com/gkanogiannis/panpath-audit/blob/main/LICENSE) and
[`NOTICE`](https://github.com/gkanogiannis/panpath-audit/blob/main/NOTICE).
