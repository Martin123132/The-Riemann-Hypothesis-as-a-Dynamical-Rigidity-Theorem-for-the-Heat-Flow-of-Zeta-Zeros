#!/usr/bin/env python3
"""Validate the conditional multiplier leading-atom interval bounds."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_multiplier_leading_atom_bound_certificate as cert  # noqa: E402


EXPECTED = (
    "validated Jensen-window PF multiplier leading-atom bound certificate: "
    "8 rows, 0 issues, 56 difference orders, "
    "beta6 in (4.863538496,4.863538497), alpha_min>4.863538496, "
    "N(11/2)<=1, 1 open existence handoff"
)


def validate_symbolic() -> list[str]:
    alpha, r = sp.symbols("alpha r", positive=True)
    m = sp.symbols("m", integer=True, nonnegative=True)
    integrand = r ** (alpha - 1) * (1 - r) ** (m + 2) / (-sp.log(r))
    derivative = sp.simplify(sp.diff(integrand, alpha))
    expected = -r ** (alpha - 1) * (1 - r) ** (m + 2)
    if sp.simplify(derivative - expected) != 0:
        return ["bad atom-kernel alpha derivative"]
    return []


def validate(stored: dict, note_path: Path) -> list[str]:
    issues: list[str] = []
    rebuilt = cert.build_payload()
    if stored != rebuilt:
        issues.append("stored payload differs from independent Arb reconstruction")
    if stored.get("kind") != "jensen_window_pf_multiplier_leading_atom_bound_certificate":
        issues.append(f"bad kind: {stored.get('kind')!r}")
    expected_summary = {
        "certificate_rows": 8,
        "difference_orders": 56,
        "root_order": 6,
        "root_bracket": ["4.863538496", "4.863538497"],
        "conditional_alpha_min_lower_bound": "4.863538496",
        "other_orders_certified_weaker": 55,
        "count_cutoff": "11/2",
        "count_order": 3,
        "max_atoms_at_or_below_cutoff": 1,
        "open_existence_handoffs": 1,
        "target_closing": False,
    }
    if stored.get("summary") != expected_summary:
        issues.append("bad summary")
    rows = {row.get("id"): row for row in stored.get("rows", [])}
    required_ids = {f"mlab_{index:02d}_{suffix}" for index, suffix in (
        (1, "atom_difference_kernel"),
        (2, "kernel_monotonicity"),
        (3, "conditional_sum"),
        (4, "root_lower_endpoint"),
        (5, "root_upper_endpoint"),
        (6, "strongest_finite_order"),
        (7, "low_atom_count"),
        (8, "open_existence_handoff"),
    )}
    if set(rows) != required_ids:
        issues.append("bad row ids")
    lower_gap = rows.get("mlab_04_root_lower_endpoint", {}).get("diagnostics", {}).get("gap_ball", "")
    upper_gap = rows.get("mlab_05_root_upper_endpoint", {}).get("diagnostics", {}).get("gap_ball", "")
    if "e-14" not in lower_gap or "e-13" not in upper_gap:
        issues.append("root endpoint margins drifted")
    strongest = rows.get("mlab_06_strongest_finite_order", {}).get("diagnostics", {})
    if strongest.get("checked_other_orders") != 55 or strongest.get("weakest_other_order") != 55:
        issues.append("strongest-order audit drifted")
    ratio = rows.get("mlab_07_low_atom_count", {}).get("diagnostics", {}).get("ratio_ball", "")
    if "1.631788551345843" not in ratio:
        issues.append("count ratio drifted")
    issues.extend(validate_symbolic())
    for key in (
        "source_hausdorff_bridge",
        "source_frontier",
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
            "Status: conditional interval leading-atom bounds; not a proof",
            "partial_alpha g_m(alpha)=-Beta(alpha,m+3)<0.",
            "4.863538496<beta_6<4.863538497",
            "N(11/2)<=1.",
            "does not prove that an atom lies below `11/2`",
        )
        issues.extend(
            f"missing note text: {item}" for item in required if item not in text
        )
    boundary = str(stored.get("proof_boundary", "")).lower()
    for marker in (
        "assuming the still-open",
        "alpha_min>4.863538496",
        "at most one atom",
        "does not prove that any atom",
        "pf-infinity",
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
