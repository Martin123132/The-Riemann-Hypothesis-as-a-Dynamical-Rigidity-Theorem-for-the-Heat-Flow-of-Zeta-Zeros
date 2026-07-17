#!/usr/bin/env python3
"""Validate inverse-thirteenth-power first-summand dominance."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_first_summand_power13_rebalanced_dominance_extension import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_artifact,
)


@dataclass(frozen=True)
class DominanceIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> DominanceIssue:
    return DominanceIssue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[DominanceIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[DominanceIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_first_summand_power13_rebalanced_dominance_extension":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "inverse-thirteenth-power" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))
    expected_summary = {
        "rows": 6,
        "ready_to_apply_rows": 6,
        "positive_analytic_gates": 14,
        "full_tail_power": 13,
        "tail_start_k": 340,
        "dominance_theorems": 1,
        "rebalanced_splits": 1,
    }
    for key, expected in expected_summary.items():
        if artifact.get("summary", {}).get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", artifact.get("summary", {}).get(key)))
    diagnostics = artifact.get("diagnostics", {})
    required = {
        "high_region_bound": "epsilon(a(k))<=17*exp(-3*pi*k^(2/5))<k^(-13), k>=340",
        "low_region_bound": "epsilon(0)*P_k(u<a(k))<k^(-13), k>=340",
        "full_tail_relative_bound": "0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^13 for every integer k>=340",
        "moment_log_error": "0<=log(1+delta_k)<2/k^13, k>=340",
        "adjacent_B_error": "|B_j-B_j^(1)|<=a_j=2*((j-1)^(-13)+2*j^(-13)+(j+1)^(-13)), j>=341",
    }
    for key, expected in required.items():
        if diagnostics.get(key) != expected:
            issues.append(issue("diagnostics", f"bad-{key}", diagnostics.get(key)))
    gates = diagnostics.get("positive_gates", {})
    if len(gates) != 14 or not all(gates.values()):
        issues.append(issue("diagnostics", "failed-gates", gates))
    rows = artifact.get("rows", [])
    if len(rows) != 6 or len({row.get("id") for row in rows}) != 6:
        issues.append(issue("rows", "bad-rows", rows))
    if any(row.get("readiness") != "ready_to_apply" for row in rows):
        issues.append(issue("rows", "bad-readiness", rows))
    try:
        rebuilt = build_artifact()
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != rebuilt:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[DominanceIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous inverse-thirteenth-power first-summand dominance",
        "This is not a proof of order eleven",
        "a(k)=log(k)/10",
        "all fourteen Arb endpoint and derivative gates are strict",
        "<k^(-13), k>=340",
        "<2/k^13 for every integer k>=340",
    )
    return [issue("note", "missing-marker", marker) for marker in required if marker not in text]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(f"POWER13-DOMINANCE {finding.section} [{finding.code}] {finding.detail}")
        print(f"power-thirteen first-summand dominance: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(f"validated power-thirteen rebalanced first-summand dominance: {summary['rows']} rows, 0 issues, {summary['positive_analytic_gates']} positive analytic gates, 1 dominance theorem")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
