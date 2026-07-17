#!/usr/bin/env python3
"""Validate the first-summand curvature bridge for the order-five tail."""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import dataclass
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order5_first_summand_curvature_bridge import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    DOMINANCE_SOURCE,
    FIRST_CONTINUOUS_CONSTANT,
    FIRST_CONTINUOUS_T,
    FIRST_DISCRETE_CONSTANT,
    FULL_TRANSFER_CONSTANT,
    ORDER4_BRIDGE_SOURCE,
    ORDER4_ENTRY_SOURCE,
    SCOUT_T,
    TAIL_FIRST_K,
    TAIL_REDUCTION_SOURCE,
    build_artifact,
)


@dataclass(frozen=True)
class BridgeIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> BridgeIssue:
    return BridgeIssue(section, code, str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def without_scout(artifact: dict) -> dict:
    stripped = deepcopy(artifact)
    stripped["scout"] = None
    stripped["summary"]["scout_rows"] = 0
    for row in stripped.get("rows", []):
        if row.get("id") == "co5fscb_08_open_continuous_target":
            row["diagnostics"] = {
                "scout_scaled_range": None,
                "formal_limit_scout": None,
            }
    return stripped


def validate_positive_polynomial(
    diagnostics: dict, expected_start: int, expected_degree: int
) -> list[BridgeIssue]:
    issues: list[BridgeIssue] = []
    if diagnostics.get("start") != expected_start:
        issues.append(issue("exact", "bad-polynomial-start", diagnostics))
    if diagnostics.get("degree") != expected_degree:
        issues.append(issue("exact", "bad-polynomial-degree", diagnostics))
    coefficients = diagnostics.get("shifted_coefficients", [])
    if len(coefficients) != expected_degree + 1:
        issues.append(issue("exact", "bad-polynomial-count", len(coefficients)))
        return issues
    try:
        parsed = [sp.Integer(value) for value in coefficients]
    except Exception as exc:
        issues.append(issue("exact", "bad-polynomial-coefficient", exc))
        return issues
    if any(value <= 0 for value in parsed):
        issues.append(issue("exact", "nonpositive-polynomial-coefficient", diagnostics))
    if diagnostics.get("leading_coefficient") != str(parsed[0]):
        issues.append(issue("exact", "bad-leading-coefficient", diagnostics))
    if diagnostics.get("constant_coefficient") != str(parsed[-1]):
        issues.append(issue("exact", "bad-constant-coefficient", diagnostics))
    return issues


def validate_artifact(path: Path) -> tuple[dict, list[BridgeIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = load_json(path)
    issues: list[BridgeIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order5_first_summand_curvature_bridge"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "one open continuous first-summand ceiling" not in str(
        artifact.get("status", "")
    ):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

    expected_parameters = {
        "tail_first_k": TAIL_FIRST_K,
        "continuous_first_t": FIRST_CONTINUOUS_T,
        "full_transfer_constant": FULL_TRANSFER_CONSTANT,
        "first_discrete_constant": FIRST_DISCRETE_CONSTANT,
        "first_continuous_constant": FIRST_CONTINUOUS_CONSTANT,
    }
    if artifact.get("parameters") != expected_parameters:
        issues.append(issue("parameters", "drift", artifact.get("parameters")))

    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 8,
        "ready_rows": 5,
        "conditional_rows": 2,
        "open_rows": 1,
        "exact_factorizations": 1,
        "positive_floor_theorems": 2,
        "full_kernel_transfer_theorems": 1,
        "scout_rows": len(SCOUT_T),
        "open_continuous_targets": 1,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    exact = artifact.get("exact", {})
    if exact.get("stable_derivative_residuals") != ["0", "0"]:
        issues.append(issue("exact", "bad-derivative-residuals", exact))
    if exact.get("continuous_target") != (
        "q_1''(t)<=60/t^2 for every real t>=320"
    ):
        issues.append(issue("exact", "bad-continuous-target", exact))
    if "<=37/k^2" not in exact.get("full_transfer", ""):
        issues.append(issue("exact", "bad-full-transfer", exact))
    if "<63/k^2" not in exact.get("tent_transfer", ""):
        issues.append(issue("exact", "bad-tent-transfer", exact))
    if "C_n<=100/k^2" not in exact.get("conditional_full_ceiling", ""):
        issues.append(issue("exact", "bad-composition", exact))

    issues.extend(
        validate_positive_polynomial(
            exact.get("gap_floor_polynomial", {}), 319, 29
        )
    )
    issues.extend(
        validate_positive_polynomial(
            exact.get("first_R_floor_polynomial", {}), 320, 2
        )
    )
    issues.extend(
        validate_positive_polynomial(
            exact.get("full_R_floor_polynomial", {}), 320, 2
        )
    )
    issues.extend(
        validate_positive_polynomial(
            exact.get("transfer_polynomial", {}), TAIL_FIRST_K, 52
        )
    )
    if exact.get("first_R_floor_polynomial", {}).get(
        "shifted_coefficients"
    ) != ["1", "597", "88622"]:
        issues.append(issue("exact", "bad-first-R-floor", exact))
    if exact.get("full_R_floor_polynomial", {}).get(
        "shifted_coefficients"
    ) != ["53", "31570", "4674200"]:
        issues.append(issue("exact", "bad-full-R-floor", exact))

    scout = artifact.get("scout", {})
    if scout.get("sample_t") != list(SCOUT_T):
        issues.append(issue("scout", "bad-sample-grid", scout))
    rows = scout.get("rows", [])
    if len(rows) != len(SCOUT_T):
        issues.append(issue("scout", "bad-row-count", len(rows)))
    for row in rows:
        if not row.get("below_target") or not row.get("positive_stable_coordinates"):
            issues.append(issue("scout", "failed-diagnostic-row", row))
        try:
            scaled = sp.Float(row.get("t2_q_second"), 40)
        except Exception as exc:
            issues.append(issue("scout", "bad-scaled-value", exc))
        else:
            if not (0 < scaled < FIRST_CONTINUOUS_CONSTANT):
                issues.append(issue("scout", "scaled-value-out-of-range", row))

    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    if len(rows) != 8 or len(ids) != len(set(ids)):
        issues.append(issue("rows", "count-or-uniqueness", ids))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or "q_1''(t)<=60/t^2" not in open_rows[0].get(
        "formula", ""
    ):
        issues.append(issue("rows", "bad-open-row", open_rows))

    for source in (
        TAIL_REDUCTION_SOURCE,
        ORDER4_ENTRY_SOURCE,
        ORDER4_BRIDGE_SOURCE,
        DOMINANCE_SOURCE,
    ):
        if not source.exists():
            issues.append(issue("sources", "missing-source", source))

    try:
        canonical = build_artifact(include_scout=False)
    except Exception as exc:  # pragma: no cover - diagnostic path
        issues.append(issue("rebuild", "exception", exc))
    else:
        if without_scout(artifact) != canonical:
            issues.append(issue("rebuild", "exact-artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[BridgeIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: exact order-five first/full curvature transfer with one open",
        "This is not a proof of order-five",
        "f(t)=g(t)^2-x(t)^3*g(t-1)*g(t+1)",
        "q(t)=log(f(t)/d(t))",
        "q''=2*h''-ell''+phi(R)*R''-chi(R)*(R')^2",
        "min(J_j,J_j^(1))>=1/(8*j)",
        "min(R_j,R_j^(1))>=7/(5*j)",
        "|C_n-C_n^(1)|<=E_(k-1)+2*E_k+E_(k+1)<=37/k^2",
        "q_1''(t)<=60/t^2 for every real t>=320",
        "<63/k^2, k>=321",
        "37+63=100",
        "The displayed continuous theorem is open.",
        "factor-ten reserve",
        "outputs/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.md",
        "outputs/formal_core.md",
    )
    issues: list[BridgeIssue] = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "the displayed continuous theorem is proved",
        "order-five entry is complete",
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
                f"ORDER5-FIRST-CURVATURE {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-five first-summand curvature bridge: {len(issues)} issues")
        return 1

    summary = artifact["summary"]
    print(
        "validated order-five first-summand curvature bridge: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['positive_floor_theorems']} positive floors, "
        f"{summary['full_kernel_transfer_theorems']} full-kernel transfer, "
        f"{summary['scout_rows']} scout rows, "
        f"{summary['open_continuous_targets']} open continuous target, "
        "budgets 37+63=100"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
