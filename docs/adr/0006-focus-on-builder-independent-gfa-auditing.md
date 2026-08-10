# Focus on builder-independent GFA auditing

_Superseded in part by ADR-0007._

Panpath Audit proceeds as a standalone auditor for source-FASTA equivalence of blunt GFA1 embedded paths, with deterministic CI results and traversal-localized diagnostics. Pangraph already verifies reconstruction for its native JSON graph format, while ODGI validates graph topology rather than equivalence to source FASTA; the product therefore makes no claim that reconstruction verification is unprecedented and differentiates itself through builder-independent auditing of distributed GFA artifacts.
