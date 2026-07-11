#!/usr/bin/env python3
"""Validate the Jensen-window PF heat-flow monotone-closure scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

from jensen_window_pf_heat_flow_monotone_closure_scout import (
    DEFAULT_ENCLOSURE_JSONL,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    finite_diagnostics,
)


EXPECTED_FORMULAS = {
    "hfmc_01_coefficient_flow_identity": "dA_k/dlambda = 2*(2*k+1)*A_{k+1}",
    "hfmc_02_contraction_log_flow": "dlog(x_k)/dlambda = 2*r_k*((2*k+3)*x_{k+1} + (2*k-1)/x_k - 2*(2*k+1))",
    "hfmc_03_monotone_gap_flow": "dlog(x_{k+1}/x_k)/dlambda = 2*r_k*((2*k+5)*x_{k+1}*x_{k+2} - 3*(2*k+3)*x_{k+1} + 3*(2*k+1) - (2*k-1)/x_k)",
    "hfmc_04_boundary_threshold_factorization": "if x_k=x_{k+1}=q and x_{k+2}>=q, bracket >= ((q-1)^2*((2*k+5)*q-(2*k-1)))/q",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Heat-Flow Monotone-Closure Scout",
    "Status: finite heat-flow closure scout",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_heat_flow_monotone_closure_scout`",
    "work/rh_compute/results/jensen_window_pf_heat_flow_monotone_closure_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_heat_flow_monotone_closure_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_monotone_closure_scout.py",
    "validated Jensen-window PF heat-flow monotone closure scout: 4 exact rows, 315 threshold rows, 305 flow-bracket rows, 0 issues",
    "dA_k/dlambda = 2*(2*k+1)*A_{k+1}",
    "dlog(x_{k+1}/x_k)/dlambda",
    "((q-1)^2*((2*k+5)*q-(2*k-1)))/q",
    "q >= (2*k-1)/(2*k+5)",
    "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
    "validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues",
    "threshold rows: 315",
    "flow-bracket rows: 305",
    "3.822433850353900366E-2",
    "2.665474275226736438E-4",
    "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
    "outputs/jensen_window_pf_monotone_contraction_stress.md",
)


@dataclass(frozen=True)
class HeatFlowClosureIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> HeatFlowClosureIssue:
    return HeatFlowClosureIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[HeatFlowClosureIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(scout: dict) -> list[HeatFlowClosureIssue]:
    issues: list[HeatFlowClosureIssue] = []
    if scout.get("kind") != "jensen_window_pf_heat_flow_monotone_closure_scout":
        issues.append(issue("<scout>", "bad-kind", repr(scout.get("kind"))))
    if scout.get("status") != "finite_heat_flow_closure_scout":
        issues.append(issue("<scout>", "bad-status", repr(scout.get("status"))))
    for key in ("source_theorem_target", "source_boundary_threshold", "source_stress", "source_stress_summary"):
        issues.extend(validate_ref("<scout>", scout.get(key)))
    for ref in scout.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<scout>", ref))
    boundary = str(scout.get("proof_boundary", "")).lower()
    for required in ("finite arb", "does not prove", "closed differential", "cauchy-binet", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<scout>", "weak-proof-boundary", required))
    return issues


def validate_exact_rows(scout: dict) -> list[HeatFlowClosureIssue]:
    rows = scout.get("exact_rows", [])
    issues: list[HeatFlowClosureIssue] = []
    if not isinstance(rows, list):
        return [issue("exact_rows", "bad-rows", repr(type(rows)))]
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for row_id, formula in EXPECTED_FORMULAS.items():
        row = by_id.get(row_id)
        if row is None:
            issues.append(issue(row_id, "missing-row", row_id))
            continue
        if row.get("formula") != formula:
            issues.append(issue(row_id, "bad-formula", repr(row.get("formula"))))
        text = " ".join(str(row.get(key, "")) for key in ("role", "derivation", "proof_boundary")).lower()
        required_terms = ("boundary", "not", "proof") if row_id == "hfmc_04_boundary_threshold_factorization" else ("exact", "not", "proof")
        for required in required_terms:
            if required not in text:
                issues.append(issue(row_id, "weak-row-boundary", required))
    threshold = by_id.get("hfmc_04_boundary_threshold_factorization", {})
    if threshold.get("sufficient_threshold") != "q >= (2*k-1)/(2*k+5)":
        issues.append(issue("hfmc_04_boundary_threshold_factorization", "bad-threshold", repr(threshold.get("sufficient_threshold"))))
    return issues


def validate_finite_diagnostics(scout: dict) -> tuple[list[HeatFlowClosureIssue], int, int]:
    issues: list[HeatFlowClosureIssue] = []
    recorded = scout.get("finite_diagnostics", {})
    max_index = int(recorded.get("max_coefficient_index", 0))
    if max_index != 64:
        issues.append(issue("finite_diagnostics", "bad-max-coefficient-index", repr(max_index)))
        max_index = 64
    recomputed = asdict(finite_diagnostics(list(DEFAULT_ENCLOSURE_JSONL), max_index))
    if recorded != recomputed:
        issues.append(issue("finite_diagnostics", "recorded-diagnostics-stale", f"{recorded!r} != {recomputed!r}"))
    expected_counts = {
        "contraction_threshold_rows": 315,
        "contraction_threshold_positive_rows": 315,
        "flow_bracket_rows": 305,
        "flow_bracket_positive_rows": 305,
    }
    for key, value in expected_counts.items():
        if recorded.get(key) != value:
            issues.append(issue("finite_diagnostics", f"bad-{key}", f"{recorded.get(key)!r} != {value!r}"))
    if recorded.get("min_threshold_margin", {}).get("sample") != "3.822433850353900366E-2":
        issues.append(issue("finite_diagnostics", "bad-min-threshold", repr(recorded.get("min_threshold_margin"))))
    if recorded.get("min_flow_bracket", {}).get("sample") != "2.665474275226736438E-4":
        issues.append(issue("finite_diagnostics", "bad-min-flow", repr(recorded.get("min_flow_bracket"))))
    return issues, int(recorded.get("contraction_threshold_rows", 0)), int(recorded.get("flow_bracket_rows", 0))


def validate_summary(scout: dict) -> list[HeatFlowClosureIssue]:
    summary = scout.get("summary", {})
    issues: list[HeatFlowClosureIssue] = []
    expected = {
        "exact_rows": 4,
        "threshold_rows": 315,
        "threshold_positive_rows": 315,
        "flow_bracket_rows": 305,
        "flow_bracket_positive_rows": 305,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("heat-flow", "boundary threshold", "finite zeta", "not a closed analytic"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in scout.get("invariants", [])).lower()
    for required in ("finite", "threshold", "forbidden", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[HeatFlowClosureIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[HeatFlowClosureIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "closed differential inequality is proved",
        "monotone contractions are proved for all zeta windows",
        "analytic monotone-contraction theorem is proved",
        "cauchy-binet identity is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(scout_path: Path, note_path: Path) -> tuple[list[HeatFlowClosureIssue], int, int, int]:
    scout = load_json(scout_path)
    issues: list[HeatFlowClosureIssue] = []
    issues.extend(validate_top_level(scout))
    issues.extend(validate_exact_rows(scout))
    finite_issues, threshold_rows, flow_rows = validate_finite_diagnostics(scout)
    issues.extend(finite_issues)
    issues.extend(validate_summary(scout))
    issues.extend(validate_note(note_path))
    exact_rows = len(scout.get("exact_rows", []))
    return issues, exact_rows, threshold_rows, flow_rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scout", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, exact_rows, threshold_rows, flow_rows = validate(args.scout, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "exact_rows": exact_rows,
                    "threshold_rows": threshold_rows,
                    "flow_bracket_rows": flow_rows,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-HEAT-FLOW-CLOSURE {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF heat-flow monotone closure scout: "
            f"{exact_rows} exact rows, {threshold_rows} threshold rows, "
            f"{flow_rows} flow-bracket rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
