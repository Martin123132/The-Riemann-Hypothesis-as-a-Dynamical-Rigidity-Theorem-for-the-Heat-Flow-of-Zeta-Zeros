#!/usr/bin/env python3
"""Validate the signed order-nine tail and cooperative-flow reduction."""

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

from jensen_window_pf_compound_order9_uniform_tail_flow_reduction import (  # noqa: E402
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


def validate_artifact(path: Path) -> tuple[dict, list[ReductionIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[ReductionIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order9_uniform_tail_flow_reduction":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "open lambda=-100 entry" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

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

    exact = artifact.get("exact", {})
    expected_exact = {
        "order9_orientation": "epsilon_9=1, Q_(9,n)=H_(9,n)",
        "signed_condensation": (
            "Q_(9,n)*Q_(7,n+2)=Q_(8,n+1)^2-Q_(8,n)*Q_(8,n+2)"
        ),
        "affine_derivative": "Q_(9,n)'=(4*n+66)*delta(Q_(9,n))",
        "open_entry": (
            "Q_(9,n)(-100)>0 for every n>=0, equivalently strict "
            "log-concavity of Q_(8,n)(-100)"
        ),
    }
    for key, expected in expected_exact.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))

    countermodel = artifact.get("countermodel", {})
    try:
        values = tuple(sp.Rational(value) for value in countermodel["sequence"])
        q9 = sp.Rational(countermodel["Q9_n0"])
        margin = sp.Rational(countermodel["Q8_log_concavity_margin"])
        residual = sp.Rational(countermodel["condensation_residual"])
    except Exception as exc:
        issues.append(issue("countermodel", "parse-failure", exc))
    else:
        if len(values) != 17 or values[-1] != sp.Rational(1, 20926400000000):
            issues.append(issue("countermodel", "bad-sequence", values))
        if q9 >= 0 or margin >= 0 or residual != 0:
            issues.append(issue("countermodel", "bad-obstruction", countermodel))
        for rows in countermodel.get("strict_signed_lower_layers", {}).values():
            if not rows or any(sp.Rational(value) <= 0 for value in rows):
                issues.append(issue("countermodel", "bad-lower-layer", rows))
                break

    rows = artifact.get("rows", [])
    if len(rows) != 9 or len({row.get("id") for row in rows}) != 9:
        issues.append(issue("rows", "bad-rows", rows))
    try:
        canonical = build_artifact()
    except Exception as exc:
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
        "Status: exact uniform eventual signed order-nine tail and conditional",
        "This is not a",
        "D=36",
        "epsilon_9=1",
        "Q_(9,n)*Q_(7,n+2)=Q_(8,n+1)^2-Q_(8,n)*Q_(8,n+2)",
        "a_n=(4*n+66)*Q_(8,n)/Q_(8,n+1)>0",
        "1/20926400000000",
        "Q_(9,0)=",
        "Q_(9,n)(-100)>0 for every n>=0",
        "outputs/formal_core.md",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "all-shift order nine is proved",
        "pf-infinity follows",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER9-TAIL-FLOW {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-nine uniform tail and flow reduction: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-nine uniform tail and flow reduction: "
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
