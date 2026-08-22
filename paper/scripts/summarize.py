#!/usr/bin/env python3
"""Create compact, manuscript-ready summaries from panpath-audit JSON reports."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path


PAPER = Path(__file__).resolve().parent.parent
RAW = PAPER / "results" / "raw"
DERIVED = PAPER / "results" / "derived"
ENVIRONMENT = PAPER / "results" / "environment.tsv"


@dataclass
class Case:
    key: str
    label: str
    graph_label: str


CASES = [
    Case("tomato_baseline", "Tomato baseline", "PGGB chr. 2"),
    Case("tomato_divergent_mm", "Tomato altered MM", "PGGB chr. 2"),
    Case("tomato_missing_sl5", "Tomato omitted SL5", "PGGB chr. 2"),
    Case("hprc_clipped", "HPRC clipped", "Minigraph--Cactus v1.0"),
    Case("hprc_full", "HPRC full", "Minigraph--Cactus v1.1"),
]


# The manuscript's environment paragraph is generated from these, so rerunning
# the pipeline on different hardware updates the paper instead of leaving a stale
# hand-written description behind. Every key is required: a silent "unknown" in a
# methods section is worse than a failed build.
REQUIRED_ENVIRONMENT = (
    "panpath_audit",
    "rustc",
    "threads",
    "memory_mib",
    "logical_cpus",
    "host_memory_kib",
    "platform",
    "data_fstype",
    "scratch_fstype",
)


def load_environment() -> dict:
    if not ENVIRONMENT.exists():
        raise SystemExit(f"missing {ENVIRONMENT}; run scripts/verify_inputs.sh first")
    values: dict[str, str] = {}
    with ENVIRONMENT.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) >= 2:
                values[row[0]] = row[1]
    missing = [key for key in REQUIRED_ENVIRONMENT if not values.get(key)]
    if missing:
        raise SystemExit(
            f"{ENVIRONMENT} lacks {', '.join(missing)}; "
            "rerun scripts/verify_inputs.sh with the current common.sh"
        )
    return values


def environment_macros(values: dict) -> list[str]:
    """LaTeX macros for the manuscript's execution-environment paragraph."""
    host_gib = int(values["host_memory_kib"]) / (1024 * 1024)
    return [
        "% Execution environment, from results/environment.tsv.",
        macro("EnvToolVersion", values["panpath_audit"].split()[-1]),
        macro("EnvRustVersion", values["rustc"].split()[1]),
        macro("EnvThreads", values["threads"]),
        macro("EnvMemoryMiB", values["memory_mib"]),
        macro("EnvLogicalCpus", values["logical_cpus"]),
        macro("EnvHostMemory", f"{host_gib:.1f}"),
        macro("EnvPlatform", values["platform"]),
        macro("EnvKernelRelease", values.get("kernel_release", "")),
        macro("EnvDataFilesystem", values["data_fstype"]),
        macro("EnvScratchFilesystem", values["scratch_fstype"]),
    ]


def load_case(case: Case) -> dict:
    path = RAW / f"{case.key}.json"
    if not path.exists():
        raise SystemExit(f"missing result: {path}")
    with path.open(encoding="utf-8") as handle:
        document = json.load(handle)
    if document.get("state") != "completed":
        raise SystemExit(f"result is not completed: {path}")
    return document


def parse_timing(case: Case) -> tuple[float, int]:
    text = (RAW / f"{case.key}.time.txt").read_text(encoding="utf-8")
    elapsed_match = re.search(
        r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)", text
    )
    rss_match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    if not elapsed_match or not rss_match:
        raise SystemExit(f"cannot parse timing for {case.key}")
    parts = [float(value) for value in elapsed_match.group(1).split(":")]
    if len(parts) == 2:
        seconds = parts[0] * 60 + parts[1]
    elif len(parts) == 3:
        seconds = parts[0] * 3600 + parts[1] * 60 + parts[2]
    else:
        raise SystemExit(f"unexpected elapsed time for {case.key}")
    return seconds, int(rss_match.group(1))


def format_integer(value: int) -> str:
    return f"{value:,}"


def result_text(summary: dict) -> str:
    parts = []
    for key, label in [
        ("identical", "identical"),
        ("divergent", "divergent"),
        ("missing_path", "missing path"),
        ("missing_source", "missing source"),
    ]:
        count = int(summary.get(key, 0))
        if count:
            parts.append(f"{format_integer(count)} {label}")
    return "; ".join(parts) if parts else "none"


def format_percent(value: object) -> str:
    return f"{float(value):.2f}"


def base_share(document: dict, key: str) -> str:
    """Percentage of total source bases sitting in one ledger bucket."""
    source = document["statistics"]["source"]
    return format_percent(100 * source["bases"][key] / source["total_bases"])


def format_duration(seconds: float) -> str:
    if seconds < 600:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.2f} h"


def latex_escape(value: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
    }
    return "".join(replacements.get(character, character) for character in value)


def macro(name: str, value: object) -> str:
    return rf"\newcommand{{\{name}}}{{{latex_escape(str(value))}}}"


def main() -> None:
    DERIVED.mkdir(parents=True, exist_ok=True)
    rows = []
    documents: dict[str, dict] = {}
    timings: dict[str, tuple[float, int]] = {}
    for case in CASES:
        document = load_case(case)
        elapsed, rss_kib = parse_timing(case)
        documents[case.key] = document
        timings[case.key] = (elapsed, rss_kib)
        stats = document["statistics"]
        record_types = "+".join(document["provenance"]["record_types"])
        row = {
            "case": case.key,
            "label": case.label,
            "graph": case.graph_label,
            "record_types": record_types,
            "source_sequences": stats["source"]["sequences"],
            "graph_traversals": stats["graph"]["traversals"],
            "source_bases": stats["source"]["total_bases"],
            "graph_bases": stats["graph"]["total_bases"],
            "not_embedded_bases": stats["source"]["bases"]["not_embedded"],
            "result": result_text(document["summary"]),
            "exit_code": int((RAW / f"{case.key}.exit_code").read_text().strip()),
            "elapsed_seconds": elapsed,
            "max_rss_kib": rss_kib,
            "tracked_peak_bytes": document["provenance"]["peak_tracked_bytes"],
        }
        rows.append(row)

    with (DERIVED / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    divergence_fields = [
        "case",
        "identifier",
        "record_type",
        "source_length",
        "path_length",
        "edit_distance",
        "substituted_bases",
        "first_source_position_1based",
        "first_path_position_1based",
        "first_segment",
        "first_orientation",
        "first_segment_position_1based",
    ]
    with (DERIVED / "hprc_divergences.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=divergence_fields, delimiter="\t")
        writer.writeheader()
        for case_key in ("hprc_clipped", "hprc_full"):
            for outcome in documents[case_key]["outcomes"]:
                if outcome["status"] != "DIVERGENT":
                    continue
                first = outcome["divergence"]
                writer.writerow(
                    {
                        "case": case_key,
                        "identifier": outcome["identifier"],
                        "record_type": outcome["record_type"],
                        "source_length": outcome["source_length"],
                        "path_length": outcome["path_length"],
                        "edit_distance": outcome["edit_distance"],
                        "substituted_bases": outcome["base_statistics"]["source"]["substituted"],
                        "first_source_position_1based": first["source_position_1based"],
                        "first_path_position_1based": first["path_position_1based"],
                        "first_segment": first.get("segment", ""),
                        "first_orientation": first.get("orientation", ""),
                        "first_segment_position_1based": first.get("segment_position_1based", ""),
                    }
                )

    table_lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Sequence audits on the two public pangenome datasets and prepared tomato controls. Runtime and resident memory describe this single \EnvPlatform{} environment and are not comparative benchmarks.}",
        r"\label{tab:results}",
        r"\small",
        r"\begin{tabularx}{\textwidth}{@{}l l c r r X r r@{}}",
        r"\toprule",
        r"Case & Graph & GFA & Sources & Traversals & Outcome & Time & RSS \\",
        r"\midrule",
    ]
    for row in rows:
        table_lines.append(
            " & ".join(
                [
                    latex_escape(row["label"]),
                    latex_escape(row["graph"]),
                    latex_escape(row["record_types"]),
                    format_integer(int(row["source_sequences"])),
                    format_integer(int(row["graph_traversals"])),
                    latex_escape(str(row["result"])),
                    format_duration(float(row["elapsed_seconds"])),
                    f'{int(row["max_rss_kib"]) / 1024:.0f} MiB',
                ]
            )
            + r" \\"
        )
    table_lines.extend([r"\bottomrule", r"\end{tabularx}", r"\end{table*}"])
    (DERIVED / "results_table.tex").write_text(
        "\n".join(table_lines) + "\n", encoding="utf-8"
    )

    tomato = documents["tomato_baseline"]
    divergent = documents["tomato_divergent_mm"]
    missing = documents["tomato_missing_sl5"]
    clipped = documents["hprc_clipped"]
    full = documents["hprc_full"]
    divergence = next(
        outcome for outcome in divergent["outcomes"] if outcome["status"] == "DIVERGENT"
    )
    missing_outcome = next(
        outcome for outcome in missing["outcomes"] if outcome["status"] == "MISSING_SOURCE"
    )
    clipped_divergences = [
        outcome for outcome in clipped["outcomes"] if outcome["status"] == "DIVERGENT"
    ]
    full_divergences = [
        outcome for outcome in full["outcomes"] if outcome["status"] == "DIVERGENT"
    ]
    macro_lines = [
        "% Generated by scripts/summarize.py; do not edit.",
        macro("TomatoIdentical", tomato["summary"]["identical"]),
        macro("TomatoSourceBases", format_integer(tomato["statistics"]["source"]["total_bases"])),
        macro("TomatoGraphBases", format_integer(tomato["statistics"]["graph"]["total_bases"])),
        macro("TomatoDivergentIdentifier", divergence["identifier"]),
        macro("TomatoDivergenceSourcePosition", format_integer(divergence["divergence"]["source_position_1based"])),
        macro("TomatoDivergencePathPosition", format_integer(divergence["divergence"]["path_position_1based"])),
        macro("TomatoDivergenceSegment", divergence["divergence"].get("segment", "end")),
        macro("TomatoDivergenceSegmentPosition", format_integer(divergence["divergence"]["segment_position_1based"])),
        macro("TomatoDivergentSourceLength", format_integer(divergence["source_length"])),
        macro("TomatoMissingIdentifier", missing_outcome["identifier"]),
        macro("HprcClippedIdentical", format_integer(clipped["summary"]["identical"])),
        macro("HprcClippedDivergent", format_integer(clipped["summary"]["divergent"])),
        macro("HprcClippedMissingPath", format_integer(clipped["summary"]["missing_path"])),
        macro("HprcClippedMissingSource", format_integer(clipped["summary"]["missing_source"])),
        macro("HprcClippedSourceSequences", format_integer(clipped["statistics"]["source"]["sequences"])),
        macro("HprcClippedGraphTraversals", format_integer(clipped["statistics"]["graph"]["traversals"])),
        macro("HprcClippedSourceBases", format_integer(clipped["statistics"]["source"]["total_bases"])),
        macro("HprcClippedGraphBases", format_integer(clipped["statistics"]["graph"]["total_bases"])),
        macro("HprcNotEmbedded", format_integer(clipped["statistics"]["source"]["bases"]["not_embedded"])),
        macro("HprcEmbeddedPercent", format_percent(clipped["statistics"]["source"]["embedded_coverage_percent"])),
        macro("HprcClippedSubstituted", format_integer(clipped["statistics"]["source"]["bases"]["substituted"])),
        macro("HprcClippedMissingPathPercent", base_share(clipped, "missing_path")),
        macro("HprcClippedNotEmbeddedPercent", base_share(clipped, "not_embedded")),
        macro("HprcClippedMissingPathBases", format_integer(clipped["statistics"]["source"]["bases"]["missing_path"])),
        macro("HprcClippedMissingSourceBases", format_integer(clipped["statistics"]["graph"]["bases"]["missing_source"])),
        macro("HprcClippedDivergentIdentifiers", ", ".join(outcome["identifier"] for outcome in clipped_divergences)),
        macro("HprcFullIdentical", format_integer(full["summary"]["identical"])),
        macro("HprcFullDivergent", format_integer(full["summary"]["divergent"])),
        macro("HprcFullMissingPath", format_integer(full["summary"]["missing_path"])),
        macro("HprcFullMissingSource", format_integer(full["summary"]["missing_source"])),
        macro("HprcFullSourceSequences", format_integer(full["statistics"]["source"]["sequences"])),
        macro("HprcFullGraphTraversals", format_integer(full["statistics"]["graph"]["traversals"])),
        macro("HprcFullSourceBases", format_integer(full["statistics"]["source"]["total_bases"])),
        macro("HprcFullGraphBases", format_integer(full["statistics"]["graph"]["total_bases"])),
        macro("HprcFullNotEmbedded", format_integer(full["statistics"]["source"]["bases"]["not_embedded"])),
        macro("HprcFullEmbeddedPercent", format_percent(full["statistics"]["source"]["embedded_coverage_percent"])),
        macro("HprcFullSubstituted", format_integer(full["statistics"]["source"]["bases"]["substituted"])),
        macro("HprcFullMissingPathPercent", base_share(full, "missing_path")),
        macro("HprcFullMissingPathBases", format_integer(full["statistics"]["source"]["bases"]["missing_path"])),
        macro("HprcFullMissingSourceBases", format_integer(full["statistics"]["graph"]["bases"]["missing_source"])),
        macro("HprcFullDivergentIdentifiers", ", ".join(outcome["identifier"] for outcome in full_divergences)),
    ]
    macro_lines.extend(environment_macros(load_environment()))
    (DERIVED / "macros.tex").write_text(
        "\n".join(macro_lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
