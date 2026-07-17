#!/usr/bin/env python3
"""Validate the uniform order-five tail and cooperative-flow reduction."""

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

from jensen_window_pf_compound_order5_uniform_tail_flow_reduction import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    HEAT_TILT_SOURCE,
    HIGHER_THETA_SOURCE,
    ORDER4_SOURCE,
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


def without_runtime_fields(value: object) -> object:
    if isinstance(value, dict):
        return {
            key: without_runtime_fields(item)
            for key, item in value.items()
            if key != "elapsed_seconds"
        }
    if isinstance(value, list):
        return [without_runtime_fields(item) for item in value]
    return value


def validate_artifact(path: Path) -> tuple[dict, list[ReductionIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = load_json(path)
    issues: list[ReductionIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order5_uniform_tail_flow_reduction"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    status = str(artifact.get("status", ""))
    if "one open lambda=-100 entry" not in status:
        issues.append(issue("artifact", "bad-status", status))

    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 12,
        "ready_to_apply_rows": 10,
        "conditional_rows": 1,
        "open_rows": 1,
        "suitability_coefficients_checked": 12,
        "lambert_derivative_orders_checked": 11,
        "newton_coefficients_checked": 8,
        "determinant_permutations_checked": 120,
        "uniform_eventual_tail_theorems": 1,
        "cooperative_flow_theorems": 1,
        "conditional_forward_theorems": 1,
        "lower_layer_countermodels": 1,
        "open_entry_targets": 1,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    rows = artifact.get("rows", [])
    if len(rows) != 12:
        issues.append(issue("rows", "bad-count", len(rows)))
    ids = [row.get("id") for row in rows]
    if len(ids) != len(set(ids)):
        issues.append(issue("rows", "duplicate-id", ids))
    if sum(row.get("readiness") == "ready_to_apply" for row in rows) != 10:
        issues.append(issue("rows", "bad-ready-count", rows))
    if sum(row.get("readiness") == "conditional_on_open_input" for row in rows) != 1:
        issues.append(issue("rows", "bad-conditional-count", rows))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or "lambda=-100" not in open_rows[0].get("claim", ""):
        issues.append(issue("rows", "bad-open-entry", open_rows))

    suitability = artifact.get("uniform_suitability_extension", {})
    coefficients = suitability.get("coefficient_checks", [])
    if [row.get("degree") for row in coefficients] != list(range(12)):
        issues.append(issue("suitability", "bad-degrees", coefficients))
    for row in coefficients:
        degree = row.get("degree")
        if row.get("T_degree", -1) > degree or row.get("v_degree", -1) > degree:
            issues.append(issue("suitability", "degree-bound", row))

    lambert = artifact.get("lambert_extension", {})
    lambert_rows = lambert.get("rows", [])
    if [row.get("order") for row in lambert_rows] != list(range(1, 12)):
        issues.append(issue("lambert", "bad-orders", lambert_rows))
    for row in lambert_rows:
        limit = sp.Integer(str(row.get("limit_w_to_infinity", "0")))
        if limit == 0:
            issues.append(issue("lambert", "zero-limit", row))

    newton = artifact.get("newton_transfer", {})
    newton_rows = newton.get("coefficient_rows", [])
    if [row.get("degree") for row in newton_rows] != list(range(1, 9)):
        issues.append(issue("newton", "bad-degrees", newton_rows))
    for row in newton_rows:
        degree = row["degree"]
        if sp.Rational(row.get("diagonal", "0")) != sp.Rational(
            1, sp.factorial(degree)
        ):
            issues.append(issue("newton", "bad-diagonal", row))
        support = row.get("difference_orders", [])
        if not support or min(support) != degree:
            issues.append(issue("newton", "bad-support", row))

    cancellation = artifact.get("determinant_cancellation", {})
    expected_coefficients = ["0"] * 10 + ["294912*G2**10"]
    if cancellation.get("coefficients_h0_through_h10") != expected_coefficients:
        issues.append(
            issue(
                "determinant",
                "bad-coefficients",
                cancellation.get("coefficients_h0_through_h10"),
            )
        )
    if cancellation.get("first_nonzero_term") != "294912*G2^10*h^10":
        issues.append(
            issue(
                "determinant",
                "bad-main-term",
                cancellation.get("first_nonzero_term"),
            )
        )
    if cancellation.get("permutations_checked") != 120:
        issues.append(
            issue(
                "determinant",
                "bad-permutation-count",
                cancellation.get("permutations_checked"),
            )
        )

    flow = artifact.get("exact_flow", {})
    residuals = flow.get("symbolic_residuals", {})
    if residuals != {"affine_heat": "0", "condensation": "0", "plucker": "0"}:
        issues.append(issue("flow", "nonzero-residual", residuals))
    if "a_n=(4*n+34)*H_(4,n)/H_(4,n+1)>0" not in flow.get(
        "cooperative_flow", ""
    ):
        issues.append(issue("flow", "missing-positive-off-diagonal", flow))

    countermodel = artifact.get("countermodel", {})
    if countermodel.get("H5_n0") != "-1/3657830400000":
        issues.append(issue("countermodel", "bad-H5", countermodel.get("H5_n0")))
    if countermodel.get("H4_log_concavity_margin") != (
        "-1/3792438558720000000"
    ):
        issues.append(
            issue(
                "countermodel",
                "bad-H4-margin",
                countermodel.get("H4_log_concavity_margin"),
            )
        )
    lower = countermodel.get("strict_signed_lower_layers", {})
    if sorted(lower) != ["1", "2", "3", "4"]:
        issues.append(issue("countermodel", "bad-lower-layers", lower))
    for order, values in lower.items():
        if not all(sp.Rational(value) > 0 for value in values):
            issues.append(issue("countermodel", f"nonpositive-order-{order}", values))

    conclusions = artifact.get("conclusions", {})
    if "there exists N_5" not in conclusions.get("uniform_eventual_tail", ""):
        issues.append(issue("conclusions", "missing-uniform-tail", conclusions))
    if "H_(5,n)(-100)>0" not in conclusions.get("conditional_forward", ""):
        issues.append(issue("conclusions", "missing-entry-condition", conclusions))

    for source in (ORDER4_SOURCE, HEAT_TILT_SOURCE, HIGHER_THETA_SOURCE, XI_RATIO_SOURCE):
        if not source.exists():
            issues.append(issue("sources", "missing-source", source))

    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover - diagnostic path
        issues.append(issue("rebuild", "exception", exc))
    else:
        canonical = without_runtime_fields(canonical)
        comparable = without_runtime_fields(artifact)
        if comparable != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[ReductionIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: exact uniform eventual order-five tail and conditional forward",
        "This is not a proof of",
        "[h^0,...,h^10] det K=[0,0,0,0,0,0,0,0,0,0,294912*G_2^10]",
        "2^10*(1!*2!*3!*4!)=294912",
        "there exists N_5 such that H_(5,n)(lambda)>0",
        "H_(5,n)*H_(3,n+2)=H_(4,n)*H_(4,n+2)-H_(4,n+1)^2",
        "H_(4,n+1)(-100)^2>H_(4,n)(-100)*H_(4,n+2)(-100)",
        "a_n=(4*n+34)*H_(4,n)/H_(4,n+1)>0",
        "The flow needs no order-six sign",
        "H_(5,0)=-1/3657830400000<0",
        "prove H_(5,n)(-100)>0 for every integer n>=0",
        "sole order-five propagation input",
        "outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
    )
    issues: list[ReductionIssue] = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "all-shift order five is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER5-TAIL-FLOW {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-five uniform tail and flow reduction: {len(issues)} issues")
        return 1

    summary = artifact["summary"]
    print(
        "validated order-five uniform tail and flow reduction: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['suitability_coefficients_checked']} suitability coefficients, "
        f"{summary['lambert_derivative_orders_checked']} Lambert orders, "
        f"{summary['newton_coefficients_checked']} Newton coefficients, "
        f"{summary['determinant_permutations_checked']} determinant permutations, "
        f"{summary['uniform_eventual_tail_theorems']} uniform tail theorem, "
        f"{summary['cooperative_flow_theorems']} cooperative flow theorem, "
        f"{summary['conditional_forward_theorems']} conditional forward theorem, "
        f"{summary['open_entry_targets']} open lambda=-100 entry"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
