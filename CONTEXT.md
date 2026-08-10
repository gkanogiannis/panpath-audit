# Panpath Audit

Panpath Audit establishes whether genome sequences embedded as paths in a pangenome graph preserve their source sequences.

## Language

**Source sequence**:
A named biological sequence supplied as the reference truth against which an embedded path is audited.
_Avoid_: Original sequence, expected sequence

**Sequence identifier**:
The case-sensitive name that uniquely pairs one source sequence with one embedded path; repeated identifiers make the pairing ambiguous and invalid.
_Avoid_: Header, path label, sequence name

**Sequence identity**:
Equality of two nucleotide-symbol sequences after normalizing letter case; ambiguity symbols remain distinct symbols and must match exactly.
_Avoid_: Byte identity, case-sensitive identity

**Divergence**:
The earliest position within an addressed range at which a source sequence and its embedded path have unequal symbols or either addressed sequence is exhausted, without implying a particular biological or graph edit.
_Avoid_: Mutation, variant, edit

**Position**:
A one-based location in a source sequence, embedded path, or segment; position 1 denotes its first nucleotide symbol.
_Avoid_: Index, zero-based offset

**Nucleotide sequence**:
A sequence over the DNA IUPAC symbols A, C, G, T, R, Y, S, W, K, M, B, D, H, V, and N, independent of letter case.
_Avoid_: DNA string, base string

**Embedded path**:
A named traversal of graph segments intended to represent a source sequence.
_Avoid_: Genome path, reconstructed genome

**Walk range**:
A half-open interval of a source sequence represented by one coordinate-bearing graph walk; non-overlapping ranges for one source are combined in coordinate order.
_Avoid_: Walk fragment, W-line

**Segment**:
A uniquely named nucleotide sequence that an embedded path may traverse in forward or reverse orientation.
_Avoid_: Node, vertex, sequence record

**Complete correspondence**:
A one-to-one pairing in which every source sequence has exactly one embedded path and every embedded path has exactly one source sequence.
_Avoid_: Coverage, matched subset

## Audit outcomes

**Identical**:
A source sequence and its paired embedded path have sequence identity across all
addressed ranges; unembedded W-range gaps are reported separately.
_Avoid_: OK, pass

**Divergent**:
A source sequence and its paired embedded path have a divergence.
_Avoid_: Mismatch, mutation

**Missing path**:
A source sequence has no paired embedded path.
_Avoid_: Absent sequence

**Missing source**:
An embedded path has no paired source sequence.
_Avoid_: Unexpected path, extra path

**Blunt path**:
An embedded path whose consecutive oriented segment sequences are joined without overlapping or trimming bases.
_Avoid_: Zero-overlap path, concatenated path
