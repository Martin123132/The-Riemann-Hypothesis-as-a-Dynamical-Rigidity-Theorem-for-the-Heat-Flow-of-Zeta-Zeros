#!/usr/bin/env python3
"""Validate the multiplier counting-measure theorem target."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_multiplier_counting_measure_target as target


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_multiplier_counting_measure_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_multiplier_counting_measure_target.md"
EXPECTED = (
    "validated Jensen-window PF multiplier counting-measure target: "
    "10 rows, 0 issues, 4 exact rows, 1 finite evidence row, 3 countermodel rows, "
    "1 live route, 0 ready-to-apply rows"
)


def exact_issues() -> list[str]:
    issues: list[str] = []
    alpha, k, z = sp.symbols("alpha k z", positive=True)
    q = alpha / (alpha + 1)
    atom = q ** (k - 1) * (alpha + k) / (alpha + 1)
    atom_next = atom.subs(k, k + 1)
    atom_prev = atom.subs(k, k - 1)
    contraction = sp.factor(atom_next * atom_prev / atom**2)
    if sp.factor(contraction - (1 - 1 / (k + alpha) ** 2)) != 0:
        issues.append("elementary contraction identity failed")
    if sp.simplify(atom.subs(k, 0) - 1) != 0 or sp.simplify(atom.subs(k, 1) - 1) != 0:
        issues.append("elementary normalization failed")
    s = sp.symbols("s", positive=True)
    kernel_log = 2 * sp.log(s) - sp.log(s - 1) - sp.log(s + 1)
    kernel_derivative = sp.factor(sp.diff(kernel_log, s))
    laplace_derivative = sp.factor(-(1 / (s - 1) + 1 / (s + 1) - 2 / s))
    if sp.factor(kernel_derivative - laplace_derivative) != 0:
        issues.append("Laplace kernel derivative identity failed")
    radical_upper = -432 * sp.Rational(707, 500) - 324 * sp.Rational(433, 250) + 378 + 324 * sp.Rational(49, 20)
    if radical_upper != sp.Rational(-27, 125):
        issues.append("fractional-weight guard bound drifted")
    return issues


def validate_payload(stored: dict) -> list[str]:
    issues: list[str] = []
    if stored != target.build_payload():
        return ["stored counting-measure target differs from reconstruction"]
    if stored.get("status") != "open theorem target":
        issues.append("target status drifted")
    summary = stored.get("summary", {})
    expected_summary = {
        "target_rows": 10,
        "exact_rows": 4,
        "finite_evidence_rows": 1,
        "countermodel_rows": 3,
        "conditional_theorem_rows": 1,
        "open_statement_rows": 1,
        "ready_to_apply_rows": 0,
        "live_routes": 1,
        "target_closing": False,
    }
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(f"summary {key} drifted")
    for key in ("source_boundary_family", "source_defect_scout", "source_bridge_target", "source_mellin_obstruction"):
        ref = stored.get(key)
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(f"missing source reference {key}")
    rows = stored.get("target_rows", [])
    if sum(row.get("readiness") == "ready_to_apply" for row in rows) != 0:
        issues.append("target has a ready-to-apply row")
    issues.extend(exact_issues())
    return issues


def validate_note(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    required = (
        EXPECTED,
        "Status: open theorem target.",
        "B_k=A_k*A_0^(k-1)/A_1^k, B_0=B_1=1.",
        "EGF(M^(alpha))=(1+z/(alpha+1))*exp(alpha*z/(alpha+1)).",
        "unit integer multiplicity",
        "-log x_k=sum_j -log(1-1/(k+alpha_j)^2)",
        "positive Laplace representation",
        "discriminant `<-27/125`",
        "discriminant `-937/3456`",
        "No such construction is currently known.",
        "sufficient subclass theorem, not a claimed characterization",
        "strictly negative",
        "rule out equality asserted only for integer `k`",
        "interpolation-uniqueness theorem",
    )
    return [f"missing note text: {value}" for value in required if value not in text]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    stored = json.loads(args.result.read_text(encoding="utf-8"))
    issues = validate_payload(stored)
    issues.extend(validate_note(args.note))
    for value in issues:
        print(f"ISSUE {value}")
    print(EXPECTED.replace("0 issues", f"{len(issues)} issues"))
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
