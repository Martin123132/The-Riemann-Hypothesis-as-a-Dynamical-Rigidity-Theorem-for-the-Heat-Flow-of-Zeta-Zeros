#!/usr/bin/env python3
"""Validate uniform compact-heat contiguous order-four forward invariance."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
import math
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "co4uhfi_01_lambda0_ratio_input",
    "co4uhfi_02_heat_tilt",
    "co4uhfi_03_kernel_remainder",
    "co4uhfi_04_uniform_ratio",
    "co4uhfi_05_uniform_eventual_tail",
    "co4uhfi_06_entry_and_flow",
    "co4uhfi_07_finite_confinement",
    "co4uhfi_08_forward_invariance",
    "co4uhfi_09_lambda_zero",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Compound Order-Four Uniform-Heat Forward Invariance",
    "all-shift contiguous order-four forward invariance",
    "same G_2 limit 1",
    "exact Newton interpolation",
    "one finite, non-effective `N`",
    "variation of constants",
    "H_(4,n)(lambda)>0 for every integer n>=0",
    "H_(4,n)(0)>0 for every integer n>=0",
    "not RH",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class InvarianceIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> InvarianceIssue:
    return InvarianceIssue(section=section, issue=name, detail=str(detail))


def validate_ref(ref: object) -> list[InvarianceIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def validate_artifact(artifact: dict) -> list[InvarianceIssue]:
    findings: list[InvarianceIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate"
    ):
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("date") != "2026-07-13":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for ref in artifact.get("sources", []):
        findings.extend(validate_ref(ref))
    for key in ("generator", "checker"):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "proves h_(4,n)(lambda)>0",
        "does not prove",
        "noncontiguous",
        "arbitrary-column",
        "order five",
        "all-degree",
        "rh",
        "lambda<=0",
    ):
        if required not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", required))

    expected_summary = {
        "rows": 9,
        "ready_to_apply_rows": 9,
        "open_analytic_rows": 0,
        "uniform_ratio_theorems": 1,
        "uniform_eventual_tail_theorems": 1,
        "finite_confinement_theorems": 1,
        "forward_invariance_theorems": 1,
        "lambda_zero_all_shift_theorems": 1,
        "global_coefficient_ceilings_used": 0,
        "newton_coefficients_checked": 6,
    }
    summary = artifact.get("summary", {})
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))

    source_diagnostics = artifact.get("source_diagnostics", {})
    checks = source_diagnostics.get("checks", {})
    if len(checks) != 7 or not all(checks.values()):
        findings.append(issue("sources", "bad-contract-checks", checks))
    if source_diagnostics.get("symbolic_main_term") != "768*G2^6*h^6":
        findings.append(issue("sources", "bad-main-term", source_diagnostics))

    exact = artifact.get("exact", {})
    expected_exact = {
        "all_interval_theorem": (
            "H_(4,n)(lambda)>0 for every integer n>=0 and every lambda in [-100,0]"
        ),
        "lambda_zero_theorem": "H_(4,n)(0)>0 for every integer n>=0",
        "cooperative_flow": "Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n>0",
    }
    for key, value in expected_exact.items():
        if exact.get(key) != value:
            findings.append(issue("exact", f"bad-{key}", exact.get(key)))
    exact_text = " ".join(str(value) for value in exact.values())
    for required in (
        "same G_2 limit 1",
        "768*G_2(lambda,n)^6",
        "there exists N",
        "integral_(-100)^lambda",
        "closed indirectly",
        "c_r(T,M)=O(log(M)/M^r)",
    ):
        if required not in exact_text:
            findings.append(issue("exact", "missing-formula", required))

    newton = artifact.get("newton_transfer", {})
    coefficient_rows = newton.get("coefficient_rows", [])
    if len(coefficient_rows) != 6:
        findings.append(issue("newton", "bad-row-count", len(coefficient_rows)))
    for degree, row in enumerate(coefficient_rows, start=1):
        if row.get("degree") != degree:
            findings.append(issue("newton", "bad-degree", row))
        if row.get("difference_orders", [None])[0] != degree:
            findings.append(issue("newton", "nontriangular-support", row))
        if row.get("diagonal") != str(Fraction(1, math.factorial(degree))):
            findings.append(issue("newton", "bad-diagonal", row))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(value) for value in ids)))
    if any(row.get("readiness") != "ready_to_apply" for row in rows):
        findings.append(issue("rows", "non-ready-row", rows))
    if any(not row.get("claim") or not row.get("proof_boundary") for row in rows):
        findings.append(issue("rows", "incomplete-row", rows))

    try:
        recomputed = build_artifact()
        for key in ("source_hashes", "source_diagnostics", "exact", "newton_transfer", "rows", "summary"):
            if artifact.get(key) != recomputed.get(key):
                findings.append(issue("recompute", f"mismatch-{key}", "stored differs"))
    except Exception as exc:
        findings.append(issue("recompute", "exception", exc))
    return findings


def validate_note(path: Path) -> list[InvarianceIssue]:
    if not path.exists():
        return [issue("note", "missing", path)]
    text = path.read_text(encoding="utf-8")
    return [
        issue("note", "missing-text", required)
        for required in REQUIRED_NOTE_STRINGS
        if required not in text
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    findings: list[InvarianceIssue] = []
    if not args.input.exists():
        findings.append(issue("artifact", "missing", args.input))
    else:
        try:
            artifact = json.loads(args.input.read_text(encoding="utf-8"))
            findings.extend(validate_artifact(artifact))
        except Exception as exc:
            findings.append(issue("artifact", "exception", exc))
    findings.extend(validate_note(args.note))

    if args.json:
        print(json.dumps([asdict(value) for value in findings], indent=2))
    elif findings:
        print(f"order-four uniform-heat forward invariance: {len(findings)} issues")
        for finding in findings:
            print(f"- [{finding.section}] {finding.issue}: {finding.detail}")
    else:
        print(
            "validated order-four uniform-heat forward invariance: "
            "9 rows, 0 issues, 9 ready rows, 1 uniform tail theorem, "
            "1 forward theorem, 1 lambda-zero all-shift theorem, 0 open rows"
        )
    raise SystemExit(1 if findings else 0)


if __name__ == "__main__":
    main()
