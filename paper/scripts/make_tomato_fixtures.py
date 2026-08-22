#!/usr/bin/env python3
"""Regenerate the two prepared tomato controls from the audit-ready source FASTA.

Both controls are derived files rather than publisher downloads, so they have to
be reproducible on any machine that holds the source set:

    divergent_MM.fa.gz   one T->A substitution at the midpoint of MM#1#chr02
    missing_SL5.fa.gz    the SL5#1#ch02 record removed entirely

Record layout is copied verbatim, so each fixture differs from the source only in
the intended way. Gzip bytes depend on the compressor, so a regenerated fixture
can have a different SHA-256 from an earlier build while carrying identical
content; verify_inputs.sh records whatever is present rather than comparing
against a stored digest.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
REPO = PAPER.parent
SOURCE = REPO / "data/tomato-tgg/sources-chr2/tomato23_chr2.sources.fa.gz"
FIXTURES = REPO / "data/tomato-tgg/fixtures"

# The substitution sits at the midpoint of the 55,739,602 bp chromosome so that
# the audit has to traverse roughly 1.12M segments before reporting it, and so
# that the reported graph coordinate carries a non-trivial within-segment offset.
# A first-base control localizes trivially and demonstrates neither.
DIVERGENT_RECORD = b"MM#1#chr02"
DIVERGENT_POSITION = 27_869_801
DIVERGENT_FROM = b"T"
DIVERGENT_TO = b"A"
MISSING_RECORD = b"SL5#1#ch02"


def write_divergent(source: Path, target: Path, level: int) -> None:
    done = False
    with gzip.open(source, "rb") as reader, gzip.open(target, "wb", level) as writer:
        in_record = False
        consumed = 0
        for line in reader:
            if line[:1] == b">":
                in_record = line[1:].strip() == DIVERGENT_RECORD
                consumed = 0
                writer.write(line)
                continue
            if in_record and not done:
                body = line.rstrip(b"\r\n")
                tail = line[len(body):]
                if consumed < DIVERGENT_POSITION <= consumed + len(body):
                    offset = DIVERGENT_POSITION - consumed - 1
                    found = body[offset:offset + 1]
                    if found != DIVERGENT_FROM:
                        raise SystemExit(
                            f"expected {DIVERGENT_FROM.decode()} at "
                            f"{DIVERGENT_RECORD.decode()}:{DIVERGENT_POSITION}, "
                            f"found {found.decode() or '<end of record>'}"
                        )
                    line = body[:offset] + DIVERGENT_TO + body[offset + 1:] + tail
                    done = True
                consumed += len(body)
            writer.write(line)
    if not done:
        raise SystemExit(f"{DIVERGENT_RECORD.decode()} shorter than {DIVERGENT_POSITION}")


def write_missing(source: Path, target: Path, level: int) -> None:
    dropped = False
    with gzip.open(source, "rb") as reader, gzip.open(target, "wb", level) as writer:
        skipping = False
        for line in reader:
            if line[:1] == b">":
                skipping = line[1:].strip() == MISSING_RECORD
                dropped = dropped or skipping
            if not skipping:
                writer.write(line)
    if not dropped:
        raise SystemExit(f"{MISSING_RECORD.decode()} not present in {source}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--out-dir", type=Path, default=FIXTURES)
    parser.add_argument("--level", type=int, default=6, help="gzip level (default 6)")
    parser.add_argument("--force", action="store_true", help="overwrite existing fixtures")
    arguments = parser.parse_args()

    if not arguments.source.exists():
        raise SystemExit(f"missing source set: {arguments.source}")
    arguments.out_dir.mkdir(parents=True, exist_ok=True)

    for name, builder in (
        ("divergent_MM.fa.gz", write_divergent),
        ("missing_SL5.fa.gz", write_missing),
    ):
        target = arguments.out_dir / name
        if target.exists() and not arguments.force:
            print(f"skip {target} (exists; pass --force to rebuild)")
            continue
        builder(arguments.source, target, arguments.level)
        print(f"wrote {target} ({target.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
