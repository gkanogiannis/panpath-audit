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

## Test data

The large integration-test data used while developing `panpath-audit` are not
stored in this repository. They are public, but the complete human set is about
100 GB compressed and the tomato preparation needs additional working space.
Everything below is written relative to the repository root and leaves the data
under the git-ignored `data/` directory.

The commands require `curl`, `awk`, `gzip`, and `sha256sum`. Tomato preparation
also requires [samtools](https://www.htslib.org/), `bgzip`, and
[GNU parallel](https://www.gnu.org/software/parallel/). Use `wget` instead of
`curl` if preferred.

### Human (HPRC)

The human test uses the HPRC Year 1 haplotype assemblies and the GRCh38-based
Minigraph-Cactus graphs. The official assembly index contains both anonymous S3
URLs and SHA-256 checksums:

```bash
mkdir -p data/human-HPRC/{assemblies,pangenome}

index=data/human-HPRC/assemblies/Year1_assemblies_v2_genbank.index
curl -fL --retry 5 \
  -o "$index" \
  https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index

# Columns 2 and 3 are the paternal and maternal S3 URLs; columns 6 and 7
# contain their SHA-256 digests. Blank maternal entries are skipped.
tail -n +2 "$index" |
  awk -F '\t' '{if ($2 != "") print $2 "\t" $6; if ($3 != "") print $3 "\t" $7}' |
  while IFS=$'\t' read -r s3_url digest; do
    https_url=${s3_url/s3:\/\/human-pangenomics/https:\/\/s3-us-west-2.amazonaws.com\/human-pangenomics}
    output=data/human-HPRC/assemblies/${https_url##*/}
    curl -fL -C - --retry 5 --speed-limit 1024 --speed-time 60 \
      -o "$output" "$https_url"
    printf '%s  %s\n' "$digest" "$output" | sha256sum -c -
  done
```

Download the two graphs used for the clipped and full-graph checks:

```bash
base=https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus

curl -fL -C - --retry 5 --speed-limit 1024 --speed-time 60 \
  -o data/human-HPRC/pangenome/hprc-v1.0-mc-grch38.gfa.gz \
  "$base/hprc-v1.0-mc-grch38/hprc-v1.0-mc-grch38.gfa.gz"
curl -fL -C - --retry 5 --speed-limit 1024 --speed-time 60 \
  -o data/human-HPRC/pangenome/hprc-v1.1-mc-grch38.full.gfa.gz \
  "$base/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.full.gfa.gz"

gzip -t data/human-HPRC/pangenome/*.gfa.gz
```

For a smaller endurance check, one assembly and the v1.0 graph are enough:

```bash
panpath-audit \
  --fasta-pansn HG00438 1 \
  data/human-HPRC/assemblies/HG00438.paternal.f1_assembly_v2_genbank.fa.gz \
  data/human-HPRC/pangenome/hprc-v1.0-mc-grch38.gfa.gz
```

The clipped graph intentionally contains gaps in some `W` walks, so the report
can include `not_embedded` bases. See the
[HPRC data-use policy](https://humanpangenome.org/data-use/) before redistributing
derived human data.

### Tomato (TGG)

The tomato test follows the public 23-assembly preparation from the
[PGGB paper workflow](https://github.com/pangenome/pggb-paper/blob/main/workflows/0.Preparation.md).
The assemblies and their MD5 files are hosted by the
[Sol Genomics Network](https://solgenomics.net/ftp/genomes/TGG/genome/):

```bash
mkdir -p data/tomato-tgg/{assemblies,sources-chr2,pangenome}
curl -fL --retry 5 \
  -o data/tomato-tgg/tomato23.urls.tsv \
  https://raw.githubusercontent.com/pangenome/pggb-paper/main/data/tomato23.urls.tsv

cut -f 2 data/tomato-tgg/tomato23.urls.tsv |
  parallel -j 4 --halt soon,fail=1 \
    'curl -fL -C - --retry 5 -O --output-dir data/tomato-tgg/assemblies {}'

# Fetch and check the publisher-provided MD5 file for every assembly.
cut -f 2 data/tomato-tgg/tomato23.urls.tsv |
  parallel -j 4 --halt soon,fail=1 \
    'curl -fL --retry 5 -O --output-dir data/tomato-tgg/assemblies {}.md5.txt'
(cd data/tomato-tgg/assemblies && for check in *.md5.txt; do md5sum -c "$check"; done)
```

The audit-ready chromosome-2 FASTA and GFA are derived files, not publisher
downloads. Prepare PanSN names and extract chromosome 2 as described upstream:

```bash
cd data/tomato-tgg/assemblies
for input in *.fasta.gz; do
  sample=${input%%.*}
  gzip -dc "$input" |
    sed "s/^>/>${sample}#1#/" |
    bgzip -@ 4 > "${sample}.pansn.fa.gz"
  samtools faidx "${sample}.pansn.fa.gz"
done

: > ../sources-chr2/tomato23_chr2.sources.fa
for input in *.pansn.fa.gz; do
  sequence=$(cut -f 1 "${input}.fai" | grep -E '#(chr|ch)0?2$' | head -n 1)
  test -n "$sequence"
  samtools faidx "$input" "$sequence" >> ../sources-chr2/tomato23_chr2.sources.fa
done
bgzip -@ 4 ../sources-chr2/tomato23_chr2.sources.fa
samtools faidx ../sources-chr2/tomato23_chr2.sources.fa.gz
cd ../../..
```

Build the chromosome graph with PGGB using the parameters from its published
tomato workflow, then copy the final GFA to
`data/tomato-tgg/pangenome/tomato23_chr2.cleaned.gfa.gz`. The graph must contain
all 23 `P` paths with the same PanSN identifiers as the prepared FASTA, inline
segment sequences, and only `*` or `0M` path overlaps. Validate those conditions
before using a graph produced by a different PGGB release.

Run the prepared test with:

```bash
panpath-audit \
  --fasta data/tomato-tgg/sources-chr2/tomato23_chr2.sources.fa.gz \
  data/tomato-tgg/pangenome/tomato23_chr2.cleaned.gfa.gz
# Expected: IDENTICAL 23 and exit status 0.
```

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
