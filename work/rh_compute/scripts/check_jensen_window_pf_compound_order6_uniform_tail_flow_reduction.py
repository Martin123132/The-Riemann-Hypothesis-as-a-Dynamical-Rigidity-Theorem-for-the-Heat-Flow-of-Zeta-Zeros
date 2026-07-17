#!/usr/bin/env python3
"""Validate the uniform signed order-six tail and flow reduction."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order6_uniform_tail_flow_reduction import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    HEAT_TILT_SOURCE,
    HIGHER_THETA_SOURCE,
    ORDER5_SOURCE,
    XI_RATIO_SOURCE,
    build_artifact,
)


@dataclass(frozen=True)
class ReductionIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> ReductionIssue:
    return ReductionIssue(section, code, str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_artifact(path: Path) -> tuple[dict, list[ReductionIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = load_json(path)
    issues: list[ReductionIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order6_uniform_tail_flow_reduction"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "one open lambda=-100 entry" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

    expected_summary = {
        "rows": 13,
        "ready_to_apply_rows": 11,
        "conditional_rows": 1,
        "open_rows": 1,
        "suitability_coefficients_checked": 17,
        "lambert_derivative_orders_checked": 16,
        "newton_coefficients_checked": 10,
        "determinant_permutations_checked": 720,
        "weighted_monomials_checked": 684,
        "uniform_eventual_tail_theorems": 1,
        "signed_condensation_recursions": 1,
        "cooperative_flow_recursions": 1,
        "conditional_forward_theorems": 1,
        "lower_layer_countermodels": 1,
        "open_entry_targets": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    rows = artifact.get("rows", [])
    if len(rows) != 13:
        issues.append(issue("rows", "bad-count", len(rows)))
    ids = [row.get("id") for row in rows]
    if len(ids) != len(set(ids)):
        issues.append(issue("rows", "duplicate-id", ids))
    if sum(row.get("readiness") == "ready_to_apply" for row in rows) != 11:
        issues.append(issue("rows", "bad-ready-count", rows))
    if sum(row.get("readiness") == "conditional_on_open_input" for row in rows) != 1:
        issues.append(issue("rows", "bad-conditional-count", rows))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or "lambda=-100" not in open_rows[0].get("claim", ""):
        issues.append(issue("rows", "bad-open-entry", open_rows))

    suitability = artifact.get("uniform_suitability_extension", {})
    coefficients = suitability.get("coefficient_checks", [])
    if [row.get("degree") for row in coefficients] != list(range(17)):
        issues.append(issue("suitability", "bad-degrees", coefficients))
    for row in coefficients:
        degree = row.get("degree")
        if row.get("T_degree", -1) > degree or row.get("v_degree", -1) > degree:
            issues.append(issue("suitability", "degree-bound", row))

    lambert_rows = artifact.get("lambert_extension", {}).get("rows", [])
    if [row.get("order") for row in lambert_rows] != list(range(1, 17)):
        issues.append(issue("lambert", "bad-orders", lambert_rows))
    for row in lambert_rows:
        if sp.Integer(str(row.get("limit_w_to_infinity", "0"))) == 0:
            issues.append(issue("lambert", "zero-limit", row))

    newton_rows = artifact.get("newton_transfer", {}).get("coefficient_rows", [])
    if [row.get("degree") for row in newton_rows] != list(range(1, 11)):
        issues.append(issue("newton", "bad-degrees", newton_rows))
    for row in newton_rows:
        degree = row["degree"]
        if sp.Rational(row.get("diagonal", "0")) != sp.Rational(
            1, sp.factorial(degree)
        ):
            issues.append(issue("newton", "bad-diagonal", row))
        if not row.get("difference_orders") or min(row["difference_orders"]) != degree:
            issues.append(issue("newton", "bad-support", row))

    cancellation = artifact.get("determinant_cancellation", {})
    if cancellation.get("coefficients_h0_through_h15") != (
        ["0"] * 15 + ["-1132462080*G2**15"]
    ):
        issues.append(issue("determinant", "bad-coefficients", cancellation))
    if cancellation.get("first_nonzero_term") != "-1132462080*G2^15*h^15":
        issues.append(issue("determinant", "bad-main-term", cancellation))
    if cancellation.get("permutation_monomial_checks") != 492480:
        issues.append(issue("determinant", "bad-check-count", cancellation))

    flow = artifact.get("exact_flow", {})
    if flow.get("symbolic_residuals") != {
        "affine_heat": "0",
        "condensation": "0",
        "plucker": "0",
    }:
        issues.append(issue("flow", "nonzero-residual", flow))
    if "a_n=(4*n+42)*H_(5,n)/H_(5,n+1)>0" not in flow.get(
        "order6_cooperative_flow", ""
    ):
        issues.append(issue("flow", "missing-positive-off-diagonal", flow))
    if "Q_(m,n)*Q_(m-2,n+2)" not in flow.get(
        "signed_condensation_recursion", ""
    ):
        issues.append(issue("flow", "missing-general-recursion", flow))

    countermodel = artifact.get("countermodel", {})
    if countermodel.get("H6_n0") != (
        "3247/6004150801779916800000000"
    ):
        issues.append(issue("countermodel", "bad-H6", countermodel))
    if countermodel.get("Q6_n0") != (
        "-3247/6004150801779916800000000"
    ):
        issues.append(issue("countermodel", "bad-Q6", countermodel))
    lower = countermodel.get("strict_signed_lower_layers", {})
    if sorted(lower) != ["1", "2", "3", "4", "5"]:
        issues.append(issue("countermodel", "bad-lower-layers", lower))
    for order, values in lower.items():
        if not all(sp.Rational(value) > 0 for value in values):
            issues.append(issue("countermodel", f"nonpositive-order-{order}", values))

    conclusions = artifact.get("conclusions", {})
    if "there exists N_6" not in conclusions.get("uniform_eventual_tail", ""):
        issues.append(issue("conclusions", "missing-uniform-tail", conclusions))
    if "Q_(6,n)(-100)>0" not in conclusions.get("conditional_forward", ""):
        issues.append(issue("conclusions", "missing-entry-condition", conclusions))

    for source in (ORDER5_SOURCE, HEAT_TILT_SOURCE, HIGHER_THETA_SOURCE, XI_RATIO_SOURCE):
        if not source.exists():
            issues.append(issue("sources", "missing-source", source))

    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover - diagnostic path
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[ReductionIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: exact uniform eventual signed order-six tail",
        "This is not a proof",
        "-1132462080*G_2^15",
        "1132462080*G_2(lambda,n)^15",
        "Q_(m,n)*Q_(m-2,n+2)",
        "a_n=(4*n+42)*H_(5,n)/H_(5,n+1)>0",
        "Q_(6,0)=-3247/6004150801779916800000000<0",
        "all-shift signed order-six entry",
        "analytic all-shift tail remains open",
        "outputs/formal_core.md",
    )
    return [
        issue("note", "missing-marker", marker)
        for marker in required
        if marker not in text
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    for finding in issues:
        print(f"ORDER6-TAIL-FLOW {finding.section} [{finding.code}] {finding.detail}")
    if issues:
        print(f"order-six uniform tail and flow reduction: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-six uniform tail and flow reduction: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['suitability_coefficients_checked']} suitability coefficients, "
        f"{summary['lambert_derivative_orders_checked']} Lambert orders, "
        f"{summary['newton_coefficients_checked']} Newton coefficients, "
        f"{summary['determinant_permutations_checked']} determinant permutations, "
        f"{summary['weighted_monomials_checked']} weighted monomials, "
        f"{summary['uniform_eventual_tail_theorems']} uniform signed tail theorem, "
        f"{summary['signed_condensation_recursions']} condensation recursion, "
        f"{summary['cooperative_flow_recursions']} cooperative recursion, "
        f"{summary['conditional_forward_theorems']} conditional forward theorem, "
        f"{summary['open_entry_targets']} open lambda=-100 entry"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
