# Do not validate graph topology

Panpath Audit validates the records and references needed to reconstruct embedded sequences but does not require consecutive path steps to be supported by GFA links. Topology validation is a separate responsibility already served by graph validators; keeping it outside the audit preserves the narrow claim that the recorded path traversal reconstructs its source sequence, and reports must disclose that topology was not validated.
