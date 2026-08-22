# panpath-audit manuscript

This directory contains the reproducible manuscript package for the
panpath-audit bioRxiv preprint. The large FASTA and GFA inputs are not copied
here; the workflow reads them through the repository-level `data` symlink.

## Reproduce the analyses

From the repository root:

```bash
paper/scripts/verify_inputs.sh
paper/scripts/run_tomato.sh
paper/scripts/run_human.sh
python3 paper/scripts/summarize.py
make -C paper
```

`verify_inputs.sh` is a gate, not just a checksum pass. Besides the FASTA and GFA
inputs it requires two download records next to the human assemblies:
`manifest.tsv` with one row per file and `verify.tsv` with one `PASS` row per
file, 96 rows each covering the 94 haplotypes, CHM13, and a GRCh38 analysis set.
The acquisition loop in the repository README verifies each digest as it
downloads but does not persist those two files, so write them from that loop or
reconstruct them with `sha256sum` before running the gate.

The GRCh38 the audit compares against is
`GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz`, the analysis
set that includes the `hs38d1` decoys. The repository README explains why the
no-decoy set gives a different and wrong answer.

The experiment scripts are resumable: a case with a JSON report and exit-code
file is skipped. Set `FORCE=1` to rerun it. Resource settings can be overridden
with `THREADS` (default 8) and `MEMORY_MIB` (default 4096).

Export `SCRATCH_DIR` once, before `verify_inputs.sh`, and every later step will
agree with it. It needs to be a fast local filesystem with room for the
decompressed graphs, which reach 48 GB and 58 GB; budget at least 64 GB free.
A default tmpfs is RAM backed and is the wrong choice at that size, and a
Windows-mounted scratch directory makes random segment access much slower.
`record_environment` reads `SCRATCH_DIR` when it writes `results/environment.tsv`,
and the manuscript reports the filesystem it names, so setting it late means the
paper describes a directory the audits never used.

Raw reports, stderr, exit codes, and GNU time telemetry are stored in
`results/raw/`. Compact derived tables and LaTeX fragments are stored in
`results/derived/`. Raw reports, scratch files, and build products are ignored
by Git; manuscript sources, scripts, input checksums, environment metadata, and
compact derived results are intended to be versioned.

`results/derived/summary.tsv` is the compact run-level record, while
`results/derived/hprc_divergences.tsv` retains the sequence and graph coordinates
for every HPRC divergence discussed in the manuscript.
`results/derived/grch38_iupac_divergences.tsv` is the supporting evidence for the
discussion of those divergences: for each affected GRCh38 chromosome it records
the first divergent position, the source symbol, the symbol the graph carries
there, the graph segment and its length, and the IUPAC ambiguity-code census of
that chromosome.

`results/reproduction/` holds the corresponding artifacts from an earlier,
independent execution of the same pipeline on different hardware, copied
verbatim rather than regenerated. It is what the manuscript's reproduction
sentence rests on: the two executions differed in processor count, operating
system, kernel, filesystem, and Rust minor version, and comparing the files
there against `results/input_checksums.tsv` and `results/derived/summary.tsv`
shows byte-identical inputs and identical outcome classifications, counts, and
exit codes. Per-sequence digests and full base ledgers were compared against the
raw JSON reports as well and were unchanged for every identifier;
`results/reproduction/README.txt` records that. The raw reports themselves are
not versioned because of their size.

The two HPRC graphs name their GRCh38 traversals differently, so
`human_source_arguments` takes the mapping as an argument: `prefix` for the v1.0
clipped graph (`P` paths named `GRCh38.chr1`) and `pansn` for the v1.1 full graph
(`W` walks named `GRCh38#0#chr1`). Using the wrong one does not fail loudly; it
reports every GRCh38 contig as missing.

## Regenerating the environment paragraph

The manuscript's execution-environment sentence and the results-table caption are
generated, not hand written. `record_environment` in `scripts/common.sh` records
the toolchain, thread count, memory budget, logical CPUs, host memory, platform,
and the filesystem types behind the data and scratch directories;
`summarize.py` turns those into `\Env...` macros in
`results/derived/macros.tex`. Rerunning the pipeline on different hardware
therefore updates the paper by itself.

`summarize.py` refuses to run if `results/environment.tsv` is missing or predates
those fields, rather than emitting an unknown value into a methods section.
If it reports a missing key, rerun `verify_inputs.sh`.
