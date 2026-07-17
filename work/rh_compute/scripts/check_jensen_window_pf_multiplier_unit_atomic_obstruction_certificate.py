#!/usr/bin/env python3
"""Validate the interval-certified unit-atomic multiplier obstruction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_multiplier_unit_atomic_obstruction_certificate as cert  # noqa: E402


EXPECTED = (
    "validated Jensen-window PF multiplier unit-atomic obstruction certificate: "
    "8 rows, 0 issues, 1 atom cutoff, 1 ratio cap, "
    "1 unit-atomic route ruled out, 0 open requirements"
)


def validate(stored: dict, note_path: Path) -> list[str]:
    issues: list[str] = []
    rebuilt = cert.build_payload()
    if stored != rebuilt:
        issues.append("stored payload differs from independent Arb reconstruction")
    if stored.get("kind") != "jensen_window_pf_multiplier_unit_atomic_obstruction_certificate":
        issues.append(f"bad kind: {stored.get('kind')!r}")
    expected_summary = {
        "certificate_rows": 8,
        "source_difference_orders": 56,
        "cutoff_order": 6,
        "atom_cutoff": "4.863538496",
        "ratio_order": 6,
        "unit_atomic_routes_ruled_out": 1,
        "arbitrary_positive_weight_routes_ruled_out": 0,
        "open_requirements": 0,
        "target_closing": True,
    }
    if stored.get("summary") != expected_summary:
        issues.append("bad summary")
    rows = {row.get("id"): row for row in stored.get("rows", [])}
    required_ids = {f"muao_{index:02d}_{suffix}" for index, suffix in (
        (1, "atom_kernel"),
        (2, "unit_sum"),
        (3, "cutoff_exclusion"),
        (4, "ratio_monotonicity"),
        (5, "weighted_ratio_cap"),
        (6, "interval_ratio_contradiction"),
        (7, "unit_product_rejected"),
        (8, "route_boundary"),
    )}
    if set(rows) != required_ids:
        issues.append("bad row ids")
    cutoff_gap = rows.get("muao_03_cutoff_exclusion", {}).get("diagnostics", {}).get("gap_ball", "")
    ratio_gap = rows.get("muao_06_interval_ratio_contradiction", {}).get("diagnostics", {}).get("gap_ball", "")
    if "3.634019802524695" not in cutoff_gap:
        issues.append("cutoff gap drifted")
    if "0.000768250097512900" not in ratio_gap:
        issues.append("ratio contradiction gap drifted")
    for key in (
        "source_leading_bound",
        "source_hausdorff_bridge",
        "source_target",
        "generator",
        "checker",
    ):
        ref = stored.get(key)
        if not isinstance(ref, str) or not (cert.REPO_ROOT / ref).exists():
            issues.append(f"missing ref {key}: {ref!r}")
    if not note_path.exists():
        issues.append(f"missing note: {note_path}")
    else:
        text = note_path.read_text(encoding="utf-8")
        required = (
            EXPECTED,
            "Status: interval-certified unit-atomic obstruction; not a proof",
            "R_m'(alpha)=Cov_alpha(1-r,log r)<0.",
            "a_7/a_6<R_6(4.863538496).",
            "gap=[0.000768250097512900",
            "do not admit the proposed convergent",
            "does not rule",
        )
        issues.extend(
            f"missing note text: {item}" for item in required if item not in text
        )
    boundary = str(stored.get("proof_boundary", "")).lower()
    for marker in (
        "rules out the specific",
        "does not rule out every multiplier sequence",
        "other all-order",
        "does not prove pf-infinity",
        "rh",
        "lambda <= 0",
    ):
        if marker not in boundary:
            issues.append(f"weak proof boundary: {marker}")
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=cert.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=cert.DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    stored = json.loads(args.result.read_text(encoding="utf-8"))
    issues = validate(stored, args.note)
    for item in issues:
        print(f"ISSUE {item}")
    print(EXPECTED.replace("0 issues", f"{len(issues)} issues"))
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
