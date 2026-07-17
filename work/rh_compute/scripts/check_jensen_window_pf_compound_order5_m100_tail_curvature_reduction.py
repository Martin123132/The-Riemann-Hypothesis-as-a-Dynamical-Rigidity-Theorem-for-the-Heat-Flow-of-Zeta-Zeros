#!/usr/bin/env python3
"""Validate the scalar curvature reduction for the order-five endpoint tail."""

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

from jensen_window_pf_compound_order5_m100_tail_curvature_reduction import (  # noqa: E402
    CURVATURE_CONSTANT,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    ORDER4_ENTRY_SOURCE,
    ORDER5_REDUCTION,
    PREFIX_SOURCE,
    TAIL_FIRST_K,
    TAIL_FIRST_N,
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
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order5_m100_tail_curvature_reduction"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "one open stable log-curvature ceiling" not in str(
        artifact.get("status", "")
    ):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 8,
        "ready_to_apply_rows": 5,
        "conditional_rows": 2,
        "open_rows": 1,
        "exact_identity_rows": 2,
        "rational_comparison_rows": 1,
        "conditional_tail_theorems": 1,
        "open_curvature_targets": 1,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    exact = artifact.get("exact", {})
    if exact.get("stable_ratio_residual") != "0":
        issues.append(issue("exact", "bad-stable-residual", exact))
    if exact.get("shifted_coefficients") != ["502", "297284", "43689082"]:
        issues.append(
            issue("exact", "bad-shifted-coefficients", exact.get("shifted_coefficients"))
        )
    if exact.get("sufficient_ceiling") != (
        "C_n<=100/k^2 for every k=n+4>=321"
    ):
        issues.append(issue("exact", "bad-ceiling", exact.get("sufficient_ceiling")))
    if exact.get("sign_equivalence") != "J_n>0 iff C_n<-4*log(x_k)":
        issues.append(issue("exact", "bad-sign-equivalence", exact.get("sign_equivalence")))

    k, m = sp.symbols("k m", integer=True, nonnegative=True)
    polynomial = sp.expand(
        502 * k**2 - 125 * CURVATURE_CONSTANT * (2 * k + 1)
    )
    shifted = sp.expand(polynomial.subs(k, TAIL_FIRST_K + m))
    if str(polynomial) != exact.get("cleared_polynomial"):
        issues.append(issue("exact", "cleared-polynomial-drift", polynomial))
    if str(shifted) != exact.get("shifted_polynomial_k_321_plus_n"):
        issues.append(issue("exact", "shifted-polynomial-drift", shifted))
    if any(value <= 0 for value in sp.Poly(shifted, m).all_coeffs()):
        issues.append(issue("exact", "nonpositive-shifted-coefficient", shifted))

    source_contract = artifact.get("source_contract", {})
    if source_contract.get("prefix_theorem") != (
        "H_(5,n)(-100)>0 for every 0<=n<=316"
    ):
        issues.append(issue("sources", "bad-prefix-contract", source_contract))
    if source_contract.get("defect_anchor") != (
        "d_k>=251/(250*(2*k+1)), k>=320"
    ):
        issues.append(issue("sources", "bad-defect-anchor", source_contract))
    for source in (PREFIX_SOURCE, ORDER5_REDUCTION, ORDER4_ENTRY_SOURCE):
        if not source.exists():
            issues.append(issue("sources", "missing-source", source))

    conclusions = artifact.get("conclusions", {})
    if "n>=317" not in conclusions.get("conditional_tail", ""):
        issues.append(issue("conclusions", "bad-tail", conclusions))
    if conclusions.get("open_target") != (
        "prove C_n<=100/(n+4)^2 for every n>=317 at lambda=-100"
    ):
        issues.append(issue("conclusions", "bad-open-target", conclusions))

    rows = artifact.get("rows", [])
    if len(rows) != 8:
        issues.append(issue("rows", "bad-count", len(rows)))
    ids = [row.get("id") for row in rows]
    if len(ids) != len(set(ids)):
        issues.append(issue("rows", "duplicate-id", ids))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or "C_n<=100" not in open_rows[0].get("formula", ""):
        issues.append(issue("rows", "bad-open-row", open_rows))

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
        "Status: exact order-five `lambda=-100` tail reduction to one open",
        "This is not a proof of the analytic tail",
        "C_n=log(F_n*F_(n+2)/F_(n+1)^2)",
        "J_n>0 iff C_n<-4*log(x_k)",
        "C_n=Delta^2 log(F_n)-Delta^2 log(d_(n+3))",
        "C_n<=100/k^2 for every k=n+4>=321",
        "100/k^2<502/(125*(2*k+1)), k>=321",
        "502*m**2 + 297284*m + 43689082>0",
        "(n+4)^2*C_n=3.5869277550969014082",
        "factor-above-27 reserve",
        "prove C_n<=100/(n+4)^2 for every n>=317 at lambda=-100",
        "outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
        "outputs/formal_core.md",
    )
    issues: list[TailIssue] = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "the curvature ceiling is proved",
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
                f"ORDER5-TAIL-CURVATURE {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-five lambda=-100 tail curvature reduction: {len(issues)} issues")
        return 1

    summary = artifact["summary"]
    print(
        "validated order-five lambda=-100 tail curvature reduction: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['exact_identity_rows']} exact identities, "
        f"{summary['rational_comparison_rows']} rational comparison, "
        f"{summary['conditional_tail_theorems']} conditional tail theorem, "
        f"{summary['open_curvature_targets']} open curvature target, "
        f"cap {CURVATURE_CONSTANT}/k^2 from k={TAIL_FIRST_K}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
