#!/usr/bin/env python3
"""Validate the signed order-eight tail and cooperative-flow reduction."""

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

from jensen_window_pf_compound_order8_uniform_tail_flow_reduction import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
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
    if artifact.get("kind") != "jensen_window_pf_compound_order8_uniform_tail_flow_reduction":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    status = str(artifact.get("status", ""))
    if "one open lambda=-100 entry" not in status:
        issues.append(issue("artifact", "bad-status", status))

    expected_summary = {
        "rows": 9,
        "ready_to_apply_rows": 6,
        "conditional_rows": 2,
        "open_rows": 1,
        "universal_tail_specializations": 1,
        "signed_condensation_coordinates": 1,
        "cooperative_flow_recursions": 1,
        "conditional_forward_theorems": 1,
        "conditional_arbitrary_column_theorems": 1,
        "lower_layer_countermodels": 1,
        "open_entry_targets": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    source = artifact.get("source_contract", {})
    if source.get("order8_leading_degree") != 28:
        issues.append(issue("source", "bad-leading-degree", source))
    if source.get("order8_raw_leading_constant") != 33664847019245568000:
        issues.append(issue("source", "bad-leading-constant", source))
    if source.get("order8_positive_signed_constant") != 33664847019245568000:
        issues.append(issue("source", "bad-signed-constant", source))

    exact = artifact.get("exact", {})
    required_exact = {
        "order8_orientation": "epsilon_8=1, Q_(8,n)=H_(8,n)",
        "leading_term": (
            "det K=33664847019245568000*G_2^28*h^28+o(h^28)"
        ),
        "signed_condensation": (
            "Q_(8,n)*Q_(6,n+2)=Q_(7,n+1)^2-"
            "Q_(7,n)*Q_(7,n+2)"
        ),
        "affine_derivative": "Q_(8,n)'=(4*n+58)*delta(Q_(8,n))",
    }
    for key, expected in required_exact.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))
    if "a_n=(4*n+58)*Q_(7,n)/Q_(7,n+1)>0" not in exact.get(
        "cooperative_flow", ""
    ):
        issues.append(issue("exact", "bad-cooperative-flow", exact))

    countermodel = artifact.get("countermodel", {})
    lower = countermodel.get("strict_signed_lower_layers", {})
    if sorted(lower) != [str(order) for order in range(1, 8)]:
        issues.append(issue("countermodel", "bad-lower-orders", lower))
    for order, values in lower.items():
        if not all(sp.Rational(value) > 0 for value in values):
            issues.append(issue("countermodel", f"nonpositive-order-{order}", values))
    if sp.Rational(countermodel.get("Q7_log_concavity_margin", "0")) >= 0:
        issues.append(issue("countermodel", "bad-q7-margin", countermodel))
    if sp.Rational(countermodel.get("Q8_n0", "0")) >= 0:
        issues.append(issue("countermodel", "bad-q8", countermodel))
    if countermodel.get("Q8_n0") != (
        "-463/69210459277322870541707704182767616000000000000000"
    ):
        issues.append(issue("countermodel", "changed-q8", countermodel.get("Q8_n0")))
    if countermodel.get("condensation_residual") != "0":
        issues.append(issue("countermodel", "bad-condensation", countermodel))

    rows = artifact.get("rows", [])
    if len(rows) != 9 or len({row.get("id") for row in rows}) != 9:
        issues.append(issue("rows", "bad-rows", rows))
    if sum(row.get("readiness") == "ready_to_apply" for row in rows) != 6:
        issues.append(issue("rows", "bad-ready-count", rows))
    if sum(row.get("readiness") == "conditional_on_open_input" for row in rows) != 2:
        issues.append(issue("rows", "bad-conditional-count", rows))
    if sum(row.get("readiness") == "not_ready_to_apply" for row in rows) != 1:
        issues.append(issue("rows", "bad-open-count", rows))

    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover
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
        "Status: exact uniform eventual signed order-eight tail",
        "This is not a",
        "D=binom(8,2)=28",
        "33664847019245568000*G_2^28*h^28",
        "Q_(8,n)*Q_(6,n+2)=Q_(7,n+1)^2-Q_(7,n)*Q_(7,n+2)",
        "Q_(8,n)'=(4*n+58)*delta(Q_(8,n))",
        "a_n=(4*n+58)*Q_(7,n)/Q_(7,n+1)>0",
        "Q_(7,1)^2-Q_(7,0)Q_(7,2)",
        "-463/69210459277322870541707704182767616000000000000000<0",
        "uniform-in-order theorem",
        "outputs/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "all-shift order eight is proved",
        "pf-infinity follows",
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
                f"ORDER8-TAIL-FLOW {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-eight uniform tail and flow reduction: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-eight uniform tail and flow reduction: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['universal_tail_specializations']} universal-tail specialization, "
        f"{summary['signed_condensation_coordinates']} condensation coordinate, "
        f"{summary['cooperative_flow_recursions']} cooperative recursion, "
        f"{summary['conditional_forward_theorems']} conditional forward theorem, "
        f"{summary['lower_layer_countermodels']} lower-cone countermodel, "
        f"{summary['open_entry_targets']} open lambda=-100 entry"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
