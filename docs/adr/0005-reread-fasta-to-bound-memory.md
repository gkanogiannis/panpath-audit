# Reread FASTA to bound memory

Panpath Audit reads every plain or gzip FASTA twice: once to index identifiers and lengths, then to audit one contig at a time. Repeated sequential decompression is preferred to retaining assemblies in memory. Source memory is bounded by the largest contig.

GFA segments and traversal records are streamed into anonymous temporary stores
in the configured scratch directory. The OS reclaims their open files after
normal or forced termination.
Workers admit source contigs through a tracked byte budget and read one graph
segment at a time. Divergent sources are temporarily spooled for one serialized
alignment.

P traversals are parsed and compared step-by-step from the traversal store. A complete P record, reconstructed chromosome, and per-base location map are never retained. W records remain grouped per source contig to preserve range stitching.

Source comparison uses a bounded worker pool. Each worker owns independent graph-store readers; the queue holds at most one pending contig per worker. Results are sorted before reporting.
