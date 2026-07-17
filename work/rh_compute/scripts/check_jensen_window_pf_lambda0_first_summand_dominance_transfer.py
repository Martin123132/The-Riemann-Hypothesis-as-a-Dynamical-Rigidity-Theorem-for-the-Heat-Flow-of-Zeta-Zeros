#!/usr/bin/env python3
"""Validate the lambda-zero first-summand dominance transfer."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_lambda0_first_summand_dominance_transfer import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "l0fsdt_01_tail_expectation",
    "l0fsdt_02_covariance_derivative",
    "l0fsdt_03_opposite_monotonicity",
    "l0fsdt_04_uniform_heat_transfer",
    "l0fsdt_05_lambda_zero_dominance",
    "l0fsdt_06_order4_penalty_transfer",
    "l0fsdt_07_curvature_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Lambda-Zero First-Summand Dominance Transfer",
    "Status: exact heat-parameter dominance transfer",
    "delta_k'(T)=-Cov_(mu_(k,T))(epsilon(U),U^2)",
    "0<=delta_k(T)<=delta_k(100)<=2/k^6",
    "k=n+3=504",
    "|P_n-P_n^(1)|<=10112*(k+1)^2/(k-3)^6<=2/(5*k^2)",
    "K_1(t)<=(7/2)/t^2 for every real t>=503",
    "does not supply that curvature theorem",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class TransferIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> TransferIssue:
    return TransferIssue(section=section, issue=name, detail=str(detail))


def validate_ref(ref: object) -> list[TransferIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def validate_artifact(artifact: dict) -> list[TransferIssue]:
    findings: list[TransferIssue] = []
    if artifact.get("kind") != "jensen_window_pf_lambda0_first_summand_dominance_transfer":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("date") != "2026-07-13":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in ("source", "generator", "checker"):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "first-summand dominance",
        "does not prove",
        "curvature",
        "all-shift order-four",
        "rh",
        "lambda<=0",
    ):
        if required not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", required))

    expected_summary = {
        "rows": 7,
        "exact_rows": 6,
        "ready_to_apply_rows": 6,
        "open_analytic_rows": 1,
        "uniform_heat_intervals": 1,
        "lambda_zero_dominance_theorems": 1,
        "order4_penalty_transfers": 1,
        "tail_start_k": 300,
        "lambda0_splice_k": 504,
    }
    summary = artifact.get("summary", {})
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))

    exact = artifact.get("exact", {})
    expected_exact = {
        "heat_monotonicity": "delta_k'(T)>=0 on 0<=T<=100",
        "lambda_zero_dominance": "0<=delta_k(0)<=2/k^6, every integer k>=300",
        "splice": "lambda=0 prefix ends at n=500; tail starts at n=501, k=n+3=504",
    }
    for key, value in expected_exact.items():
        if exact.get(key) != value:
            findings.append(issue("exact", f"bad-{key}", exact.get(key)))
    coefficients = exact.get("perturbation_coefficients_descending", [])
    if len(coefficients) != 7 or any(not isinstance(value, int) or value <= 0 for value in coefficients):
        findings.append(issue("exact", "bad-perturbation-coefficients", coefficients))
    for required in (
        "-Cov",
        "Cov(epsilon(U),U^2)<=0",
        "delta_k(100)<=2/k^6",
        "<=2/(5*k^2)",
    ):
        if required not in " ".join(str(value) for value in exact.values()):
            findings.append(issue("exact", "missing-formula", required))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(value) for value in ids)))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or open_rows[0].get("id") != "l0fsdt_07_curvature_handoff":
        findings.append(issue("rows", "bad-open-row", open_rows))
    if any(not row.get("claim") or not row.get("proof_boundary") for row in rows):
        findings.append(issue("rows", "incomplete-row", rows))

    try:
        recomputed = build_artifact()
        for key in ("source_sha256", "exact", "rows", "summary"):
            if artifact.get(key) != recomputed.get(key):
                findings.append(issue("recompute", f"mismatch-{key}", "stored differs"))
    except Exception as exc:
        findings.append(issue("recompute", "exception", exc))
    return findings


def validate_note(path: Path) -> list[TransferIssue]:
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
    findings: list[TransferIssue] = []
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
        print(f"lambda-zero first-summand dominance transfer: {len(findings)} issues")
        for finding in findings:
            print(f"- [{finding.section}] {finding.issue}: {finding.detail}")
    else:
        print(
            "validated lambda-zero first-summand dominance transfer: "
            "7 rows, 0 issues, 6 exact rows, 1 lambda-zero dominance theorem, "
            "1 order-four penalty transfer, 1 open curvature tail"
        )
    raise SystemExit(1 if findings else 0)


if __name__ == "__main__":
    main()
