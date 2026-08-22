#!/usr/bin/env bash
set -euo pipefail
export SCRATCH_DIR=${SCRATCH_DIR:-/tmp/panpath-audit-paper-tomato}
source "$(dirname -- "$0")/common.sh"

require_file "$BINARY"
graph="$DATA_DIR/tomato-tgg/pangenome/tomato23_chr2.cleaned.gfa.gz"

run_audit tomato_baseline 0 \
    "$DATA_DIR/tomato-tgg/sources-chr2/tomato23_chr2.sources.fa.gz" "$graph"
run_audit tomato_divergent_mm 1 \
    "$DATA_DIR/tomato-tgg/fixtures/divergent_MM.fa.gz" "$graph"
run_audit tomato_missing_sl5 2 \
    "$DATA_DIR/tomato-tgg/fixtures/missing_SL5.fa.gz" "$graph"
