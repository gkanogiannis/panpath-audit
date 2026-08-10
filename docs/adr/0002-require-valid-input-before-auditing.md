# Require valid input before auditing

Panpath Audit performs an all-or-nothing preflight: it reports all detectable violations of the accepted input contract and produces no sequence-audit results if any are found. Once preflight succeeds, it audits every source/path correspondence and reports all audit failures rather than stopping early; this prevents partial comparisons from being mistaken for proof that the complete inputs were preserved.
