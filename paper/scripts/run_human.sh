#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "$0")/common.sh"

require_file "$BINARY"
mapfile -d '' -t clipped_sources < <(human_source_arguments prefix)
[[ ${#clipped_sources[@]} -eq 195 ]]
mapfile -d '' -t full_sources < <(human_source_arguments pansn)
[[ ${#full_sources[@]} -eq 196 ]]

run_audit hprc_clipped any "${clipped_sources[@]}" \
    "$DATA_DIR/human-HPRC/pangenome/hprc-v1.0-mc-grch38.gfa.gz"
run_audit hprc_full any "${full_sources[@]}" \
    "$DATA_DIR/human-HPRC/pangenome/hprc-v1.1-mc-grch38.full.gfa.gz"
