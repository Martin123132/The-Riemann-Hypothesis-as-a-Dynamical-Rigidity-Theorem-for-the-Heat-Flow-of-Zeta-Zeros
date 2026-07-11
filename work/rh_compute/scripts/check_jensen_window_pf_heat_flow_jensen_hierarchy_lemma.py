#!/usr/bin/env python3
"""Validate the exact Jensen-window heat-flow hierarchy lemma."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_heat_flow_jensen_hierarchy_lemma as lemma


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.md"
EXPECTED = (
    "validated Jensen-window PF heat-flow Jensen hierarchy lemma: "
    "9 rows, 0 issues, 5 exact hierarchy identities, 2 cubic countermodels, "
    "1 open higher-minor handoff, 0 ready-to-apply rows"
)


def symbolic_issues() -> list[str]:
    issues: list[str] = []
    d, n, j = sp.symbols("d n j", integer=True, nonnegative=True)
    ode_lhs = 2 * (2 * (n + j) + 1)
    ode_rhs = (4 * n + 2) + 4 * j
    if sp.expand(ode_lhs - ode_rhs) != 0:
        issues.append("lambda hierarchy coefficient identity failed")
    if sp.simplify(sp.binomial(d + 1, j) - sp.binomial(d, j) - sp.binomial(d, j - 1)) != 0:
        issues.append("Pascal shift-degree identity failed")
    if sp.simplify(j * sp.binomial(d, j) - d * sp.binomial(d - 1, j - 1)) != 0:
        issues.append("z-derivative identity failed")
    x, y, z = sp.symbols("x y z")
    polynomial = 1 + 3 * z + 3 * x * z**2 + x**2 * y * z**3
    discriminant = sp.factor(sp.discriminant(polynomial, z))
    expected = -27 * x**2 * (x**2 * y**2 - 6 * x * y + 4 * x + 4 * y - 3)
    if sp.expand(discriminant - expected) != 0:
        issues.append("cubic discriminant identity failed")
    if sp.factor(discriminant.subs({x: sp.Rational(1, 2), y: 1})) != sp.Rational(-27, 16):
        issues.append("cubic countermodel value failed")
    n = sp.Integer(0)
    xb = sp.Rational(5, 9)
    yb = sp.Rational(21, 25)
    zb = sp.Rational(589, 625)
    frontier = x**2 * y**2 - 6 * x * y + 4 * x + 4 * y - 3
    frontier_x = sp.diff(frontier, x)
    frontier_y = sp.diff(frontier, y)
    x_dot_over_r = 2 * x * ((2 * n + 5) * x * y - 2 * (2 * n + 3) * x + (2 * n + 1))
    y_dot_over_r = 2 * x * y * ((2 * n + 7) * y * z - 2 * (2 * n + 5) * y + (2 * n + 3))
    frontier_dot_over_r = sp.factor(frontier_x * x_dot_over_r + frontier_y * y_dot_over_r)
    substitutions = {x: xb, y: yb, z: zb}
    if sp.factor(frontier.subs(substitutions)) != 0:
        issues.append("dynamic cubic witness is not on the boundary")
    if sp.factor(frontier_dot_over_r.subs(substitutions)) != sp.Rational(329728, 2109375):
        issues.append("dynamic cubic boundary derivative drifted")
    d1, d2, d3 = 1 - xb, 1 - yb, 1 - zb
    if sp.factor(d1 * d3 - d2**2) != 0 or sp.factor(d2 / d1) != sp.Rational(9, 25):
        issues.append("one-atom Hausdorff defect witness drifted")
    return issues


def validate_payload(stored: dict) -> list[str]:
    issues: list[str] = []
    rebuilt = lemma.build_payload()
    if stored != rebuilt:
        return ["stored hierarchy payload differs from reconstruction"]
    summary = stored.get("summary", {})
    expected = {
        "rows": 9,
        "exact_identity_rows": 5,
        "exact_definition_rows": 1,
        "exact_cubic_countermodels": 2,
        "open_hierarchy_handoffs": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(f"summary {key} drifted")
    for key in ("source_flow_certificate", "source_bridge_target", "source_defect_scout"):
        ref = stored.get(key)
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(f"missing source reference {key}")
    issues.extend(symbolic_issues())
    return issues


def validate_note(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    required = (
        EXPECTED,
        "Status: exact heat-flow hierarchy lemma and open higher-minor handoff.",
        "J_(d+1,n)=J_(d,n)+z*J_(d,n+1)",
        "partial_z J_(d,n)=d*J_(d-1,n+1)",
        "partial_lambda J_(d,n)=(4*n+2)*J_(d,n+1)+4*z*partial_z J_(d,n+1)",
        "discriminant is `-27/16`",
        "(partial_lambda F)/r_0=329728/2109375>0",
        "one-atom Hausdorff moment",
        "heat derivative is",
        "negative at this boundary point",
        "needs an additional shift-coupled invariant",
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
