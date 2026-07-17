#!/usr/bin/env python3
"""Validate the order-seven endpoint-tail curvature reduction."""

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

from jensen_window_pf_compound_order7_m100_tail_curvature_reduction import (  # noqa: E402
    CURVATURE_CONSTANT,
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


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_artifact(path: Path) -> tuple[dict, list[TailIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = load_json(path)
    issues: list[TailIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order7_m100_tail_curvature_reduction":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    status = str(artifact.get("status", ""))
    if "one open fourth-nested curvature ceiling" not in status:
        issues.append(issue("artifact", "bad-status", status))

    expected_summary = {
        "rows": 7,
        "ready_to_apply_rows": 6,
        "exact_factorizations": 1,
        "exact_curvature_reductions": 1,
        "coefficient_positive_comparisons": 1,
        "conditional_tail_theorems": 1,
        "open_curvature_targets": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    exact = artifact.get("exact", {})
    required_exact = {
        "fourth_gap": "T(t)=5*B(t)-p(t-1)+2*p(t)-p(t+1)",
        "order7_coordinate": (
            "r(t)=2*p(t)-q(t)+log(1-exp(-T(t)))"
        ),
        "canonical_factorization": (
            "Q_(6,n)=A_(n+5)^6*exp(r(n+5))"
        ),
        "canonical_factorization_residual": "0",
        "tail_index": "k=n+6, so n>=315 iff k>=321",
        "sufficient_ceiling": "R_k<=900/k^2 for every real/integer k>=321",
        "log_buffer": (
            "-6*log(x_k)>=6*d_k>=753/(125*(2*k+1)), k>=320"
        ),
        "rational_comparison": (
            "900/k^2<753/(125*(2*k+1)), k>=321"
        ),
    }
    for key, expected in required_exact.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))
    if exact.get("shifted_coefficients") != ["753", "258426", "5252373"]:
        issues.append(issue("exact", "bad-shifted-coefficients", exact.get("shifted_coefficients")))

    m = sp.symbols("m", nonnegative=True, integer=True)
    shifted = sp.sympify(
        exact.get("shifted_polynomial_k_321_plus_m", "0"), locals={"m": m}
    )
    expected_shifted = sp.expand(
        sp.Integer(753) * (TAIL_FIRST_K + m) ** 2
        - sp.Integer(125) * CURVATURE_CONSTANT * (2 * (TAIL_FIRST_K + m) + 1)
    )
    if sp.expand(shifted - expected_shifted) != 0:
        issues.append(issue("exact", "bad-shifted-polynomial", shifted))

    rows = artifact.get("rows", [])
    if len(rows) != 7 or len({row.get("id") for row in rows}) != 7:
        issues.append(issue("rows", "bad-rows", rows))
    if sum(row.get("readiness") == "ready_to_apply" for row in rows) != 6:
        issues.append(issue("rows", "bad-ready-count", rows))
    if sum(row.get("readiness") == "not_ready_to_apply" for row in rows) != 1:
        issues.append(issue("rows", "bad-open-count", rows))

    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover - diagnostic path
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
        "Status: exact endpoint-tail reduction with one open fourth-nested",
        "This is not a proof of order-seven entry",
        "T(t)=5*B(t)-p(t-1)+2*p(t)-p(t+1)",
        "r(t)=2*p(t)-q(t)+log(1-exp(-T(t)))",
        "Q_(6,n)=A_(n+5)^6*exp(r(n+5))",
        "E_n=log(Q_(6,n)*Q_(6,n+2)/Q_(6,n+1)^2)",
        "R_k<=900/k^2 for every real/integer k>=321",
        "900/k^2<753/(125*(2*k+1))",
        "753*m**2 + 258426*m + 5252373>0",
        "Prove R_k<=900/k^2 for every k>=321",
        "common `t+-5` cover",
        "derivatives through order",
        "outputs/jensen_window_pf_compound_order7_m100_prefix_certificate.md",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "order-seven entry is proved",
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
                f"ORDER7-TAIL-CURVATURE {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-seven tail curvature reduction: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-seven tail curvature reduction: "
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
