#!/usr/bin/env python3
"""Check the exact alternate quartic length-13 survivor gate."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import sympy as sp


sys.set_int_max_str_digits(0)

REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.py"
)


def load_builder():
    spec = importlib.util.spec_from_file_location("quartic_alt13_builder", BUILDER)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load alternate survivor builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def direct_terminal_signs(exact: dict) -> dict[str, int]:
    """Independently reconstruct the stored rational tail and terminal signs."""
    fixed = {
        int(name.split("_")[1]): sp.Rational(value)
        for name, value in exact["fixed_outer_contact"].items()
    }
    parameters = {
        int(name.split("_")[1]): sp.Rational(value)
        for name, value in exact["corridor_parameters"].items()
    }
    x = dict(fixed)

    def d(index: int) -> sp.Expr:
        return sp.cancel(1 - x[index])

    def gap(index: int) -> sp.Expr:
        return sp.cancel(
            d(index + 2) ** 2
            - x[index + 2] ** 2 * d(index + 1) * d(index + 3)
        )

    def cap(index: int) -> sp.Expr:
        return sp.cancel(
            gap(index - 4) ** 2
            / (x[index - 2] ** 3 * gap(index - 5))
        )

    for index in range(6, 14):
        target = sp.cancel(parameters[index] * cap(index))
        d_new = sp.cancel(
            (d(index - 1) ** 2 - target)
            / (x[index - 1] ** 2 * d(index - 2))
        )
        x[index] = sp.cancel(1 - d_new)
        stored = exact["contractions"][f"x_{index}"]
        if hashlib.sha256(str(x[index]).encode("ascii")).hexdigest() != stored[
            "sha256"
        ]:
            raise RuntimeError(f"stored x_{index} hash did not reproduce")

    def delta(index: int) -> sp.Expr:
        repeated = sp.cancel(
            d(index - 1) ** 2
            - x[index - 1] ** 2 * d(index - 2) * d(index - 1)
        )
        return sp.cancel(cap(index) - repeated)

    return {
        "positive_order3_gaps": sum(
            1 for index in range(1, 11) if gap(index) > 0
        ),
        "positive_order4_margins": sum(
            1
            for index in range(6, 14)
            if cap(index) - gap(index - 3) > 0
        ),
        "delta_13_sign": int(sp.sign(delta(13))),
        "delta_14_sign": int(sp.sign(delta(14))),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    issues: list[str] = []
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    expected_kind = (
        "jensen_window_pf_quartic_outer_branch_"
        "alternate_length13_survivor_gate"
    )
    if payload.get("kind") != expected_kind:
        issues.append("artifact kind changed")
    if payload.get("status") != (
        "exact alternate length-thirteen survivor and fixed-tail obstruction"
    ):
        issues.append("artifact status changed")
    if len(payload.get("rows", [])) != 14:
        issues.append("expected 14 gate rows")

    builder = load_builder()
    parent = json.loads(builder.PARENT_RESULT.read_text(encoding="utf-8"))
    rebuilt = builder.build_exact(parent)
    if payload.get("exact") != rebuilt:
        issues.append("stored exact survivor data did not reproduce")
    parent_hash = hashlib.sha256(builder.PARENT_RESULT.read_bytes()).hexdigest()
    if payload.get("parent", {}).get("sha256") != parent_hash:
        issues.append("parent artifact hash changed")

    direct = direct_terminal_signs(rebuilt)
    if direct != {
        "positive_order3_gaps": 10,
        "positive_order4_margins": 8,
        "delta_13_sign": 1,
        "delta_14_sign": -1,
    }:
        issues.append(f"independent terminal signs changed: {direct}")

    expected_summary = {
        "rows": 14,
        "fixed_outer_contacts": 1,
        "corridor_parameters": 8,
        "contractions": 12,
        "strict_adjacent_scalar_steps": 11,
        "positive_order3_gaps": 10,
        "positive_order4_margins": 8,
        "signed_order2_minors": 364,
        "signed_order3_minors": 715,
        "signed_order4_minors": 792,
        "positive_length13_compatibilities": 1,
        "negative_fixed_tail_length14_compatibilities": 1,
        "uniform_length14_theorems": 0,
        "repaired_scope_handoffs": 1,
    }
    if payload.get("summary") != expected_summary:
        issues.append("summary counts changed")

    scalar = rebuilt["scalar_corridors"]
    if (
        scalar["strict_increase_count"] != 11
        or scalar["strict_upper_count"] != 12
        or scalar["strict_point_wall_count"] != 12
        or scalar["strict_scaled_defect_increase_count"] != 11
        or scalar["strict_cubic_frontier_count"] != 11
        or len(scalar["reciprocal_certificates"]) != 11
    ):
        issues.append("scalar corridor counts changed")
    branch_counts = {
        branch: sum(
            row["branch"] == branch
            for row in scalar["reciprocal_certificates"]
        )
        for branch in {
            row["branch"] for row in scalar["reciprocal_certificates"]
        }
    }
    if branch_counts != {
        "automatic_nonpositive_comparison": 1,
        "positive_comparison_one_squaring": 10,
    }:
        issues.append(f"reciprocal proof branches changed: {branch_counts}")

    minors = rebuilt["signed_layers"]["arbitrary_column_minors"]
    for order, count in {"2": 364, "3": 715, "4": 792}.items():
        if minors[order]["count"] != count:
            issues.append(f"order-{order} minor count changed")
        if minors[order]["strict_positive"] != count:
            issues.append(f"order-{order} minor sign count changed")
        if len(minors[order]["sha256"]) != 64:
            issues.append(f"order-{order} determinant hash changed shape")

    scout = rebuilt["bounded_delta14_scout"]
    if scout.get("status") != "numerical_evidence_only":
        issues.append("bounded scout was promoted beyond numerical evidence")
    if scout.get("uniform_samples") != 800000:
        issues.append("bounded scout sample count changed")

    note = args.note.read_text(encoding="utf-8")
    required_markers = [
        "same outer quartic contact",
        "different tail",
        "order 4: 792 strict signed minors",
        "Delta_13",
        "Delta_14",
        "prefix-specific",
        "numerical evidence only",
        "not promoted to a",
        "uniform theorem",
        "PF-infinity",
        "`Lambda <= 0`",
        "RH",
    ]
    for marker in required_markers:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated alternate quartic length-13 survivor gate: "
        "14 rows, 0 issues, 8 exact corridor parameters, "
        "12 contractions, 11 strict adjacent scalar steps, "
        "10 positive order-three gaps, 8 positive order-four margins, "
        "364 signed order-two minors, 715 signed order-three minors, "
        "792 signed order-four minors, 1 positive length-13 compatibility, "
        "1 negative fixed-tail length-14 compatibility, "
        "0 uniform length-14 theorems, 1 repaired scope handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
