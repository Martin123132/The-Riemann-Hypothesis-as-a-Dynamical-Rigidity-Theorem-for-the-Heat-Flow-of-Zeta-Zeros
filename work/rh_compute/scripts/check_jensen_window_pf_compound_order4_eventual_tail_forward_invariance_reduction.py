#!/usr/bin/env python3
"""Validate the order-four eventual-tail forward-invariance reduction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "co4etfi_01_entry",
    "co4etfi_02_cooperative_flow",
    "co4etfi_03_uniform_tail_target",
    "co4etfi_04_finite_confinement",
    "co4etfi_05_variation_of_constants",
    "co4etfi_06_backward_induction",
    "co4etfi_07_forward_invariance",
    "co4etfi_08_lambda_zero_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Compound Order-Four Eventual-Tail Forward-Invariance Reduction",
    "exact finite-confinement forward-invariance reduction",
    "Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n>0",
    "exists N such that Q_n(lambda)>0",
    "Variation of constants gives",
    "backward induction",
    "needs no supremum bound",
    "endpoint statement does not yet control the interior heat interval",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class ReductionIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> ReductionIssue:
    return ReductionIssue(section=section, issue=name, detail=str(detail))


def validate_ref(ref: object) -> list[ReductionIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def validate_artifact(artifact: dict) -> list[ReductionIssue]:
    findings: list[ReductionIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction"
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
        "uniform eventual",
        "does not prove",
        "unconditional",
        "arbitrary-column",
        "rh",
        "lambda<=0",
    ):
        if required not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", required))

    expected_summary = {
        "rows": 8,
        "exact_or_input_rows": 7,
        "ready_to_apply_rows": 4,
        "conditional_rows": 3,
        "open_analytic_rows": 1,
        "finite_confinement_reductions": 1,
        "variation_of_constants_identities": 1,
        "conditional_forward_theorems": 1,
        "global_coefficient_ceilings_required": 0,
    }
    summary = artifact.get("summary", {})
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))

    exact = artifact.get("exact", {})
    exact_text = " ".join(str(value) for value in exact.values()).lower()
    for required in (
        "Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n>0",
        "exists N",
        "integral_(-100)^lambda",
        "successively",
        "every n>=0",
        "eventual positivity at lambda=0 alone",
    ):
        if required.lower() not in exact_text:
            findings.append(issue("exact", "missing-formula", required))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(value) for value in ids)))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or open_rows[0].get("id") != "co4etfi_03_uniform_tail_target":
        findings.append(issue("rows", "bad-open-row", open_rows))
    conditional = [row for row in rows if row.get("readiness") == "conditional"]
    if len(conditional) != 3:
        findings.append(issue("rows", "bad-conditional-count", len(conditional)))
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


def validate_note(path: Path) -> list[ReductionIssue]:
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
    findings: list[ReductionIssue] = []
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
        print(f"order-four eventual-tail forward invariance: {len(findings)} issues")
        for finding in findings:
            print(f"- [{finding.section}] {finding.issue}: {finding.detail}")
    else:
        print(
            "validated order-four eventual-tail forward invariance: "
            "8 rows, 0 issues, 7 exact/input rows, 1 finite-confinement reduction, "
            "1 conditional forward theorem, 1 open uniform tail"
        )
    raise SystemExit(1 if findings else 0)


if __name__ == "__main__":
    main()
