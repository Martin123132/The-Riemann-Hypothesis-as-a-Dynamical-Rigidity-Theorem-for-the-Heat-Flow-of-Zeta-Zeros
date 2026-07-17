#!/usr/bin/env python3
"""Validate the order-nine lambda=-100 tail-curvature reduction."""

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

from jensen_window_pf_compound_order9_m100_tail_curvature_reduction import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    TAIL_FIRST_K,
    build_artifact,
)


@dataclass(frozen=True)
class TailIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> TailIssue:
    return TailIssue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[TailIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[TailIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order9_m100_tail_curvature_reduction":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "sixth-nested curvature ceiling" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))
    expected_summary = {
        "rows": 7,
        "ready_to_apply_rows": 6,
        "exact_factorizations": 1,
        "exact_curvature_reductions": 1,
        "coefficient_positive_comparisons": 1,
        "conditional_tail_theorems": 1,
        "open_curvature_targets": 1,
        "required_top_derivative_order": 16,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    exact = artifact.get("exact", {})
    expected_exact = {
        "sixth_gap": "V(t)=7*B(t)-s(t-1)+2*s(t)-s(t+1)",
        "order9_coordinate": "w(t)=2*s(t)-r(t)+log(1-exp(-V(t)))",
        "canonical_factorization": "Q_(8,n)=A_(n+7)^8*exp(w(n+7))",
        "tail_index": "k=n+8, so n>=1241 iff k>=1249",
        "sufficient_ceiling": "Y_k<=4900/k^2 for every real/integer k>=1249",
    }
    for key, expected in expected_exact.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))
    try:
        coefficients = [sp.Integer(value) for value in exact["shifted_coefficients"]]
    except Exception as exc:
        issues.append(issue("exact", "bad-shifted-coefficients", exc))
    else:
        if not coefficients or any(value <= 0 for value in coefficients):
            issues.append(issue("exact", "nonpositive-shifted-coefficient", coefficients))
    if TAIL_FIRST_K != 1249:
        issues.append(issue("constants", "bad-tail-first-k", TAIL_FIRST_K))
    try:
        canonical = build_artifact()
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[TailIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: exact endpoint-tail reduction with one open sixth-nested",
        "This is not a proof",
        "V(t)=7*B(t)-s(t-1)+2*s(t)-s(t+1)",
        "Q_(8,n)=A_(n+7)^8*exp(w(n+7))",
        "k=n+8, so n>=1241 iff k>=1249",
        "4900/k^2<1004/(125*(2*k+1)), k>=1249",
        "Prove Y_k<=4900/k^2 for every k>=1249",
        "through order sixteen",
        "outputs/formal_core.md",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
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
                f"ORDER9-TAIL-CURVATURE {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-nine tail curvature reduction: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-nine tail curvature reduction: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['exact_factorizations']} exact factorization, "
        f"{summary['exact_curvature_reductions']} exact curvature reduction, "
        f"{summary['coefficient_positive_comparisons']} positive comparison, "
        f"{summary['conditional_tail_theorems']} conditional tail theorem, "
        f"{summary['open_curvature_targets']} open curvature target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
