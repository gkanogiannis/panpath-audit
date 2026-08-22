#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PAPER_DIR=$(cd -- "$SCRIPT_DIR/.." && pwd)
REPO_DIR=$(cd -- "$PAPER_DIR/.." && pwd)
DATA_DIR="$REPO_DIR/data"
RESULT_DIR="$PAPER_DIR/results/raw"
DERIVED_DIR="$PAPER_DIR/results/derived"
SCRATCH_DIR=${SCRATCH_DIR:-"$PAPER_DIR/scratch"}
BINARY="$REPO_DIR/target/release/panpath-audit"
THREADS=${THREADS:-8}
MEMORY_MIB=${MEMORY_MIB:-4096}

mkdir -p "$RESULT_DIR" "$DERIVED_DIR" "$SCRATCH_DIR"

require_file() {
    if [[ ! -f "$1" ]]; then
        printf 'required file not found: %s\n' "$1" >&2
        exit 1
    fi
}

# Filesystem type backing a path. The manuscript reports these, because input
# and scratch storage dominate wall-clock time far more than CPU does.
fstype_of() {
    findmnt -no FSTYPE --target "$1" 2>/dev/null \
        || stat -f -c %T "$1" 2>/dev/null \
        || printf 'unknown'
}

detect_platform() {
    case "$(uname -r)" in
        *microsoft-standard-WSL2*) printf 'WSL2' ;;
        *[Mm]icrosoft*) printf 'WSL' ;;
        *) printf 'Linux' ;;
    esac
}

# Every value the manuscript's environment paragraph needs is recorded here and
# turned into a LaTeX macro by summarize.py, so a rerun on different hardware
# updates the paper instead of leaving a stale hand-written description.
#
# Export SCRATCH_DIR before calling this if the audits will use a scratch
# location other than the default, otherwise scratch_fstype describes a
# directory the runs never touched.
record_environment() {
    local output="$PAPER_DIR/results/environment.tsv"
    {
        printf 'key\tvalue\n'
        printf 'timestamp_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'git_commit\t%s\n' "$(git -C "$REPO_DIR" rev-parse HEAD)"
        printf 'panpath_audit\t%s\n' "$($BINARY --version)"
        printf 'rustc\t%s\n' "$(rustc --version)"
        printf 'cargo\t%s\n' "$(cargo --version)"
        printf 'kernel\t%s\n' "$(uname -srmo)"
        printf 'kernel_release\t%s\n' "$(uname -r)"
        printf 'platform\t%s\n' "$(detect_platform)"
        printf 'logical_cpus\t%s\n' "$(nproc)"
        printf 'threads\t%s\n' "$THREADS"
        printf 'memory_mib\t%s\n' "$MEMORY_MIB"
        printf 'host_memory_kib\t%s\n' "$(awk '/^MemTotal:/{print $2}' /proc/meminfo)"
        printf 'data_fstype\t%s\n' "$(fstype_of "$DATA_DIR")"
        printf 'scratch_dir\t%s\n' "$SCRATCH_DIR"
        printf 'scratch_fstype\t%s\n' "$(fstype_of "$SCRATCH_DIR")"
    } > "$output"
}

run_audit() {
    local case_name=$1
    local expected_exit=$2
    shift 2
    local json="$RESULT_DIR/$case_name.json"
    local stderr="$RESULT_DIR/$case_name.stderr"
    local timing="$RESULT_DIR/$case_name.time.txt"
    local exit_file="$RESULT_DIR/$case_name.exit_code"

    if [[ ${FORCE:-0} != 1 && -s "$json" && -s "$exit_file" ]]; then
        printf 'skip %s (existing result)\n' "$case_name"
        return 0
    fi

    printf 'run %s\n' "$case_name"
    set +e
    /usr/bin/time -v -o "$timing" \
        "$BINARY" --format json --stats comprehensive \
        --threads "$THREADS" --memory-mib "$MEMORY_MIB" \
        --temp-dir "$SCRATCH_DIR" "$@" > "$json" 2> "$stderr"
    local status=$?
    set -e
    printf '%s\n' "$status" > "$exit_file"

    python3 -m json.tool "$json" >/dev/null
    if [[ "$expected_exit" != any && "$status" -ne "$expected_exit" ]]; then
        printf '%s: expected exit %s, observed %s\n' \
            "$case_name" "$expected_exit" "$status" >&2
        return 1
    fi
    printf 'done %s (exit %s)\n' "$case_name" "$status"
}

# GRCh38 is named differently in the two HPRC graphs: the v1.0 clipped graph
# embeds it as P paths named GRCh38.chrN, while the v1.1 full graph embeds it as
# PanSN walks named GRCh38#0#chrN. Callers pass the mapping that matches the
# graph being audited; using the wrong one silently reclassifies every GRCh38
# contig as missing rather than comparing it.
human_source_arguments() {
    local grch38_mapping=${1:?GRCh38 mapping required: prefix|pansn}
    local assembly_dir="$DATA_DIR/human-HPRC/assemblies"
    local path
    for path in "$assembly_dir"/*.paternal.f1_assembly_v2_genbank.fa.gz; do
        printf '%s\0%s\0' --fasta "$path"
    done
    for path in "$assembly_dir"/*.maternal.f1_assembly_v2_genbank.fa.gz; do
        printf '%s\0%s\0' --fasta "$path"
    done
    printf '%s\0%s\0%s\0%s\0' --fasta-pansn CHM13 0 \
        "$assembly_dir/chm13.draft_v1.1.fasta.gz"
    local grch38="$assembly_dir/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
    case "$grch38_mapping" in
        prefix)
            printf '%s\0%s\0%s\0' --fasta-prefix GRCh38. "$grch38"
            ;;
        pansn)
            printf '%s\0%s\0%s\0%s\0' --fasta-pansn GRCh38 0 "$grch38"
            ;;
        *)
            printf 'unknown GRCh38 mapping: %s\n' "$grch38_mapping" >&2
            return 1
            ;;
    esac
}
