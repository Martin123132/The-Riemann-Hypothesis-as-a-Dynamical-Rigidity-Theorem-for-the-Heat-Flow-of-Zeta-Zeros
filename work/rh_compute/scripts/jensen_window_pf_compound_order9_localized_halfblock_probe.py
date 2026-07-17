#!/usr/bin/env python3
"""Probe one rigorous half-unit order-nine localized continuation block."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
import flint  # noqa: E402
from jensen_window_pf_compound_order9_localized_final_gap_interval_core import (  # noqa: E402
    PRECISION_BITS,
    localized_sixth_continuation_row,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_h_rows(path: Path) -> list[dict]:
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    rows = []
    for expected, record in enumerate(records):
        if record.get("index") != expected or record.get("passed") is not True:
            raise RuntimeError(f"invalid exact-t source row {expected}")
        if set(record.get("h_derivatives", {})) != {
            str(order) for order in range(2, 23)
        }:
            raise RuntimeError(f"source row {expected} does not contain H2-H22")
        rows.append(
            {
                "target_t_left": Fraction(record["target_t_left"]),
                "target_t_right": Fraction(record["target_t_right"]),
                "H": {
                    order: compact.interval_from_text(
                        record["h_derivatives"][str(order)]
                    )
                    for order in range(2, 23)
                },
            }
        )
    return rows


def load_point_h_source(path: Path) -> dict[Fraction, tuple[list, dict]]:
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    source = {}
    for record in records:
        if record.get("passed") is not True:
            raise RuntimeError("point H source contains a failed row")
        target = Fraction(record["target_t"])
        derivatives = record.get("h_derivatives", {})
        if set(derivatives) != {str(order) for order in range(9)}:
            raise RuntimeError(f"point H source at {target} is not H0-H8")
        source[target] = (
            [
                compact.interval_from_text(derivatives[str(order)])
                / math.factorial(order)
                for order in range(9)
            ],
            {
                "target_t": str(target),
                "mode_bracket": [record["mode_left"], record["mode_right"]],
                "maximum_panel_error_upper": record[
                    "maximum_panel_error_upper"
                ],
                "maximum_tail_moment_upper": record[
                    "maximum_tail_moment_upper"
                ],
                "minimum_tail_slope_lower": record[
                    "minimum_tail_slope_lower"
                ],
            },
        )
    return source


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor", type=Fraction, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--point-cache", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    flint.ctx.prec = PRECISION_BITS
    source = args.cache.resolve()
    rows = load_h_rows(source)
    point_source = (
        load_point_h_source(args.point_cache.resolve())
        if args.point_cache is not None
        else None
    )
    right = args.anchor + Fraction(1, 2)
    result = localized_sixth_continuation_row(
        args.anchor,
        right,
        rows,
        point_order=6,
        remainder_order=7,
        point_h_source=point_source,
    )
    artifact = {
        "kind": "jensen_window_pf_compound_order9_localized_halfblock_probe",
        "date": "2026-07-14",
        "status": "rigorous single-block probe; not a global curvature theorem",
        "proof_boundary": (
            "This artifact certifies only its displayed half-unit t block. "
            "It does not certify the untested intervals between probe anchors."
        ),
        "source": {
            "cache": source.relative_to(REPO_ROOT).as_posix(),
            "sha256": sha256(source),
            "rows": len(rows),
            "h_derivative_orders": [2, 22],
            "point_cache": (
                args.point_cache.resolve().relative_to(REPO_ROOT).as_posix()
                if args.point_cache is not None
                else None
            ),
            "point_cache_sha256": (
                sha256(args.point_cache.resolve())
                if args.point_cache is not None
                else None
            ),
        },
        "result": result,
    }
    payload = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
