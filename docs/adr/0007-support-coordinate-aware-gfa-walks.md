# Support coordinate-aware GFA walks

GFA 1.1 `W` records are grouped as `sample#haplotype#sequence` and ordered by their half-open source ranges. Each traversal is compared only with its declared `[seqStart, seqEnd)` range. Prefix, internal, and suffix gaps are unembedded source bases, not divergences. A mismatch or traversal-length disagreement within a declared range is a divergence. Overlapping, ambiguous, or invalid ranges are invalid input, and a derived `W` identifier may not collide with a `P` name.
