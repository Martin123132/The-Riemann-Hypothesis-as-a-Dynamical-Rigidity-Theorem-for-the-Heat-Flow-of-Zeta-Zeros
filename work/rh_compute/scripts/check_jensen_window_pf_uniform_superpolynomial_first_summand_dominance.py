#!/usr/bin/env python3
"""Validate uniform superpolynomial first-summand dominance."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_uniform_superpolynomial_first_summand_dominance import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "usfsd_01_low_region_asymptotic",
    "usfsd_02_high_region_asymptotic",
    "usfsd_03_m100_superpolynomial",
    "usfsd_04_uniform_heat_transfer",
    "usfsd_05_local_log_differences",
    "usfsd_06_uniform_ratio_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Uniform Superpolynomial First-Summand Dominance",
    "uniform compact-heat superpolynomial first-summand dominance",
    "lim_(k->infinity) (B_low(k)+p*log(k))*log(k)/k=-8/5",
    "17*exp(-3*pi*sqrt(k))",
    "delta_k(T)<=delta_k(100)",
    "|nabla^m e_k(T)|<=2^m",
    "remaining asymptotic problem is",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class DominanceIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> DominanceIssue:
    return DominanceIssue(section=section, issue=name, detail=str(detail))


def validate_ref(ref: object) -> list[DominanceIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def validate_artifact(artifact: dict) -> list[DominanceIssue]:
    findings: list[DominanceIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_uniform_superpolynomial_first_summand_dominance"
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
        "superpolynomial",
        "does not prove",
        "first-summand heat-tilt",
        "uniform order-four",
        "rh",
        "lambda<=0",
    ):
        if required not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", required))

    expected_summary = {
        "rows": 6,
        "exact_rows": 6,
        "ready_to_apply_rows": 6,
        "open_analytic_rows": 0,
        "superpolynomial_tail_theorems": 1,
        "uniform_heat_transfers": 1,
        "local_log_difference_theorems": 1,
        "higher_theta_handoffs_closed": 1,
    }
    summary = artifact.get("summary", {})
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))

    exact = artifact.get("exact", {})
    if exact.get("symbolic_low_scaled_limit") != "-8/5":
        findings.append(issue("exact", "bad-low-limit", exact.get("symbolic_low_scaled_limit")))
    if exact.get("symbolic_high_scaled_limit") != "-3*pi":
        findings.append(issue("exact", "bad-high-limit", exact.get("symbolic_high_scaled_limit")))
    exact_text = " ".join(str(value) for value in exact.values())
    for required in (
        "for every p>0",
        "sup_(0<=T<=100)",
        "2^m*max",
        "O_p,m(k^(-p))",
        "o(k^-N)",
    ):
        if required not in exact_text:
            findings.append(issue("exact", "missing-formula", required))

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
        for key in ("source_hashes", "source_diagnostics", "exact", "rows", "summary"):
            if artifact.get(key) != recomputed.get(key):
                findings.append(issue("recompute", f"mismatch-{key}", "stored differs"))
    except Exception as exc:
        findings.append(issue("recompute", "exception", exc))
    return findings


def validate_note(path: Path) -> list[DominanceIssue]:
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
    findings: list[DominanceIssue] = []
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
        print(f"uniform superpolynomial first-summand dominance: {len(findings)} issues")
        for finding in findings:
            print(f"- [{finding.section}] {finding.issue}: {finding.detail}")
    else:
        print(
            "validated uniform superpolynomial first-summand dominance: "
            "6 rows, 0 issues, 6 exact rows, 1 superpolynomial theorem, "
            "1 local-difference theorem, 0 open rows"
        )
    raise SystemExit(1 if findings else 0)


if __name__ == "__main__":
    main()
