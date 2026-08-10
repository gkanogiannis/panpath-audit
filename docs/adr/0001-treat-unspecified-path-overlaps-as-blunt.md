# Treat unspecified path overlaps as blunt

Panpath Audit v0.1 treats a GFA `P` record's `*` overlap field as an instruction to concatenate its oriented segment sequences without trimming, and discloses that interpretation in its report. This departs from the field's technical meaning of "unspecified" because the blunt pangenome graphs targeted by the tool commonly use `*`; rejecting it would exclude the primary workflows, while accepting sequence-consuming overlaps would make reconstruction ambiguous and is therefore unsupported.
