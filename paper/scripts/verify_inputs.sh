#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "$0")/common.sh"

require_file "$BINARY"
require_file "$DATA_DIR/human-HPRC/assemblies/manifest.tsv"
require_file "$DATA_DIR/human-HPRC/assemblies/verify.tsv"
require_file "$DATA_DIR/human-HPRC/assemblies/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
require_file "$DATA_DIR/human-HPRC/pangenome/hprc-v1.0-mc-grch38.gfa.gz"
require_file "$DATA_DIR/human-HPRC/pangenome/hprc-v1.1-mc-grch38.full.gfa.gz"
require_file "$DATA_DIR/tomato-tgg/sources-chr2/tomato23_chr2.sources.fa.gz"
require_file "$DATA_DIR/tomato-tgg/pangenome/tomato23_chr2.cleaned.gfa.gz"
require_file "$DATA_DIR/tomato-tgg/fixtures/divergent_MM.fa.gz"
require_file "$DATA_DIR/tomato-tgg/fixtures/missing_SL5.fa.gz"

[[ $(wc -l < "$DATA_DIR/human-HPRC/assemblies/manifest.tsv") -eq 96 ]]
[[ $(awk '$1=="PASS"{count++} END{print count+0}' \
    "$DATA_DIR/human-HPRC/assemblies/verify.tsv") -eq 96 ]]
[[ $(wc -l < "$DATA_DIR/tomato-tgg/assemblies/manifest.tsv") -eq 23 ]]

gzip -t "$DATA_DIR/human-HPRC/pangenome/hprc-v1.0-mc-grch38.gfa.gz"
gzip -t "$DATA_DIR/human-HPRC/pangenome/hprc-v1.1-mc-grch38.full.gfa.gz"
gzip -t "$DATA_DIR/tomato-tgg/sources-chr2/tomato23_chr2.sources.fa.gz"
gzip -t "$DATA_DIR/tomato-tgg/pangenome/tomato23_chr2.cleaned.gfa.gz"
gzip -t "$DATA_DIR/tomato-tgg/fixtures/divergent_MM.fa.gz"
gzip -t "$DATA_DIR/tomato-tgg/fixtures/missing_SL5.fa.gz"

checksums="$PAPER_DIR/results/input_checksums.tsv"
{
    printf 'sha256\tbytes\tpath\n'
    for path in \
        "$DATA_DIR/human-HPRC/assemblies/manifest.tsv" \
        "$DATA_DIR/human-HPRC/assemblies/verify.tsv" \
        "$DATA_DIR/human-HPRC/assemblies/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz" \
        "$DATA_DIR/human-HPRC/pangenome/hprc-v1.0-mc-grch38.gfa.gz" \
        "$DATA_DIR/human-HPRC/pangenome/hprc-v1.1-mc-grch38.full.gfa.gz" \
        "$DATA_DIR/tomato-tgg/sources-chr2/tomato23_chr2.sources.fa.gz" \
        "$DATA_DIR/tomato-tgg/pangenome/tomato23_chr2.cleaned.gfa.gz" \
        "$DATA_DIR/tomato-tgg/fixtures/divergent_MM.fa.gz" \
        "$DATA_DIR/tomato-tgg/fixtures/missing_SL5.fa.gz"; do
        digest=$(sha256sum "$path" | awk '{print $1}')
        size=$(stat -c %s "$path")
        printf '%s\t%s\t%s\n' "$digest" "$size" "${path#"$REPO_DIR/"}"
    done
} > "$checksums"

record_environment
printf 'verified input manifests and compressed artifacts\n'
