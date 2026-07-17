#!/usr/bin/env python3
"""Validate the exact multiplier Hausdorff-uniqueness bridge."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_multiplier_hausdorff_uniqueness_bridge as bridge  # noqa: E402


EXPECTED = (
    "validated Jensen-window PF multiplier Hausdorff-uniqueness bridge: "
    "10 rows, 0 issues, 1 Hausdorff measure theorem, "
    "1 unit-atomic characterization, 1 interpolation guard, "
    "1 open recovery handoff"
)
REQUIRED_IDS = {f"mhb_{index:02d}_{suffix}" for index, suffix in (
    (1, "atom_laplace_kernel"),
    (2, "counting_laplace_sum"),
    (3, "hausdorff_pushforward"),
    (4, "hausdorff_uniqueness"),
    (5, "complete_monotonicity"),
    (6, "unit_atomic_characterization"),
    (7, "atom_uniqueness"),
    (8, "finite_frontier_input"),
    (9, "periodic_interpolation_guard"),
    (10, "open_recovery_handoff"),
)}


def validate_symbolic_identities() -> list[str]:
    issues: list[str] = []
    z, t, r, s = sp.symbols("z t r s", positive=True)
    exponential_kernel = sp.exp(-z * t) * 2 * (sp.cosh(t) - 1)
    frullani_split = (
        sp.exp(-(z - 1) * t) + sp.exp(-(z + 1) * t) - 2 * sp.exp(-z * t)
    )
    if sp.simplify(exponential_kernel - frullani_split) != 0:
        issues.append("bad Frullani exponential split")
    rational_product = sp.factor((z / (z - 1)) * (z / (z + 1)))
    if sp.simplify(rational_product - 1 / (1 - z ** -2)) != 0:
        issues.append("bad Frullani logarithm product")
    for order in range(17):
        alternating_difference = sum(
            (-1) ** j * sp.binomial(order, j) * r**j
            for j in range(order + 1)
        )
        if sp.expand(alternating_difference - (1 - r) ** order) != 0:
            issues.append(f"bad alternating finite-difference identity at order {order}")
    periodic = sp.sin(2 * sp.pi * s)
    if sp.simplify(periodic.subs(s, s + 1) - 2 * periodic + periodic.subs(s, s - 1)) != 0:
        issues.append("bad periodic interpolation guard")
    return issues


def validate(stored: dict, note_path: Path) -> list[str]:
    issues: list[str] = []
    rebuilt = bridge.build_payload()
    if stored != rebuilt:
        issues.append("stored payload differs from exact reconstruction")
    if stored.get("kind") != "jensen_window_pf_multiplier_hausdorff_uniqueness_bridge":
        issues.append(f"bad kind: {stored.get('kind')!r}")
    ids = {row.get("id") for row in stored.get("rows", [])}
    if ids != REQUIRED_IDS:
        issues.append(f"bad row ids: missing={sorted(REQUIRED_IDS-ids)}, extra={sorted(ids-REQUIRED_IDS)}")
    expected_summary = {
        "bridge_rows": 10,
        "exact_rows": 7,
        "hausdorff_measure_theorems": 1,
        "unit_atomic_characterizations": 1,
        "finite_frontier_rows": 1,
        "interpolation_guards": 1,
        "open_recovery_handoffs": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    if stored.get("summary") != expected_summary:
        issues.append("bad summary")
    issues.extend(validate_symbolic_identities())
    for key in (
        "source_frontier",
        "source_target",
        "source_mellin_guard",
        "generator",
        "checker",
    ):
        ref = stored.get(key)
        if not isinstance(ref, str) or not (bridge.REPO_ROOT / ref).exists():
            issues.append(f"missing ref {key}: {ref!r}")
    if not note_path.exists():
        issues.append(f"missing note: {note_path}")
    else:
        text = note_path.read_text(encoding="utf-8")
        required = (
            EXPECTED,
            "Status: exact Hausdorff uniqueness and unit-atomic characterization.",
            "dnu(r)=q(-log r)*S(-log r)dr.",
            "r^(k-1)*(1-r)^m dnu(r)>=0.",
            "unit integer multiplicities",
            "f(s)=sin(2*pi*s)",
            "unique Hausdorff measure nu",
        )
        issues.extend(
            f"missing note text: {item}" for item in required if item not in text
        )
    boundary = str(stored.get("proof_boundary", "")).lower()
    for marker in (
        "proves uniqueness",
        "does not construct",
        "natural mellin interpolation",
        "pf-infinity",
        "rh",
        "lambda <= 0",
    ):
        if marker not in boundary:
            issues.append(f"weak proof boundary: {marker}")
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=bridge.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=bridge.DEFAULT_NOTE)
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
