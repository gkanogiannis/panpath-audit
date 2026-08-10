# Usage

```text
panpath-audit [--format human|json|tsv] [--threads N] [--memory-mib N] \
  [--temp-dir DIR] \
  [--stats comprehensive] \
  [--alignment-max-cells N] [--fasta FILE]... \
  [--fasta-pansn SAMPLE HAP FILE]... \
  [--fasta-prefix PREFIX FILE]... GFA
```

Use `panpath-audit --help` for the option summary and `--version` for the tool
version. `--` permits positional paths beginning with `-`.

FASTA and GFA inputs may be plain or gzip-compressed. Source and segment
sequences must use DNA IUPAC symbols. Sequence identifiers must be unique across
all FASTA inputs. FASTAs are streamed twice; source memory is bounded by the
largest contig.

Positional `FASTA... GFA` remains shorthand for exact identifiers.

Scratch defaults to the OS temporary directory. `--temp-dir` selects another
filesystem. Stores are anonymous and reclaimed after forced termination.

`--threads` defaults to available CPUs capped at 8. Use `--threads 1` for
sequential comparison.

`--memory-mib` limits concurrent tracked sequence bytes and defaults to 1024.
One oversized source contig runs alone. Reports include the limit, tracked peak,
and oversized-contig count. Parser metadata and OS buffers are not counted.

`--stats comprehensive` adds exact edit-aligned base statistics. JSON reports
the panpath-audit package version as its schema version. Detailed alignments are
serialized.
The wavefront limit defaults to 10,000,000 cells. Limit exhaustion keeps the
audit result and reports the traversal as unaligned.
Alignment also reports `memory_limit_exceeded` when both sequences cannot fit
the memory budget.

Source bases are partitioned into matched, substituted, source-only,
not-embedded, missing-path, and unaligned-divergent bases. Graph traversal bases
use matched, substituted, graph-only, missing-source, and unaligned-divergent
categories. Percentages use separate source and graph denominators.

JSON and TSV include SHA-256 and BLAKE3 for the full source, addressed source
ranges, and reconstructed graph traversal. Symbols are normalized to uppercase.
Human output remains summary-only.

Exit codes are 0 clean, 1 divergent, 2 missing or ambiguous correspondence,
3 invalid GFA, and 4 source, usage, scratch, or operational failure.
When conditions coexist, the highest code takes precedence.
