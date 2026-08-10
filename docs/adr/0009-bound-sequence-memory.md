# Bound sequence memory

Comparison uses a configurable tracked-sequence budget, capacity-one queues, and
streaming `P`/`W` traversal readers. Only first-divergence coordinates persist.
The numeric segment index uses safe 12-byte entries. An oversized source contig
runs alone; an oversized alignment is reported without materialization.
