#!/usr/bin/env python3
"""Validate the Jensen-window PF heat-flow boundary-threshold lemma."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from dataclasses import dataclass
import json
from pathlib import Path

from jensen_window_pf_heat_flow_boundary_threshold_lemma import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    finite_diagnostics,
)
from jensen_window_pf_heat_flow_monotone_closure_scout import DEFAULT_ENCLOSURE_JSONL  # noqa: E402


EXPECTED_FORMULAS = {
    "hfbtl_01_phi_positive_half_line": "Phi(u)=sum_{n>=1} pi*n^2*exp(5*u)*(2*pi*n^2*exp(4*u)-3)*exp(-pi*n^2*exp(4*u)) > 0 for u>=0",
    "hfbtl_02_raw_moment_logconvexity": "mu_{2k}(lambda)^2 <= mu_{2k-2}(lambda)*mu_{2k+2}(lambda)",
    "hfbtl_03_normalized_strong_threshold": "x_k = ((2*k-1)/(2*k+1))*(mu_{2k+2}*mu_{2k-2}/mu_{2k}^2) >= (2*k-1)/(2*k+1)",
    "hfbtl_04_heat_flow_threshold": "x_k >= (2*k-1)/(2*k+1) > (2*k-1)/(2*k+5)",
    "hfbtl_05_boundary_bracket_nonnegative": "if x_k=x_{k+1}=q, x_{k+2}>=q, and q>=((2*k-1)/(2*k+5)), then bracket >= ((q-1)^2*((2*k+5)*q-(2*k-1)))/q >= 0",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Heat-Flow Boundary-Threshold Lemma",
    "Status: exact boundary-threshold lemma",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_heat_flow_boundary_threshold_lemma`",
    "work/rh_compute/results/jensen_window_pf_heat_flow_boundary_threshold_lemma.json",
    "python work/rh_compute/scripts/jensen_window_pf_heat_flow_boundary_threshold_lemma.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_boundary_threshold_lemma.py",
    "validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues",
    "Phi(u)=sum_{n>=1}",
    "2*pi-3 > 0",
    "mu_{2k}(lambda)^2 <= mu_{2k-2}(lambda)*mu_{2k+2}(lambda)",
    "x_k >= (2*k-1)/(2*k+1) > (2*k-1)/(2*k+5)",
    "boundary pointing condition",
    "strong-threshold rows: 315",
    "heat-threshold rows: 315",
    "8.170843282044743874E-3",
    "3.822433850353900366E-2",
    "outputs/jensen_window_pf_heat_flow_monotone_closure_scout.md",
    "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
)


@dataclass(frozen=True)
class BoundaryThresholdIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> BoundaryThresholdIssue:
    return BoundaryThresholdIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[BoundaryThresholdIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(lemma: dict) -> list[BoundaryThresholdIssue]:
    issues: list[BoundaryThresholdIssue] = []
    if lemma.get("kind") != "jensen_window_pf_heat_flow_boundary_threshold_lemma":
        issues.append(issue("<lemma>", "bad-kind", repr(lemma.get("kind"))))
    if lemma.get("status") != "available_exact_boundary_threshold_lemma":
        issues.append(issue("<lemma>", "bad-status", repr(lemma.get("status"))))
    for key in ("source_heat_flow_closure", "source_theorem_target"):
        issues.extend(validate_ref("<lemma>", lemma.get(key)))
    for ref in lemma.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<lemma>", ref))
    boundary = str(lemma.get("proof_boundary", "")).lower()
    for required in ("exact boundary-threshold", "does not prove", "adjacent log-concavity", "closed differential", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<lemma>", "weak-proof-boundary", required))
    return issues


def validate_exact_rows(lemma: dict) -> tuple[list[BoundaryThresholdIssue], int]:
    rows = lemma.get("exact_rows", [])
    issues: list[BoundaryThresholdIssue] = []
    if not isinstance(rows, list):
        return [issue("exact_rows", "bad-rows", repr(type(rows)))], 0
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for row_id, formula in EXPECTED_FORMULAS.items():
        row = by_id.get(row_id)
        if row is None:
            issues.append(issue(row_id, "missing-row", row_id))
            continue
        if row.get("formula") != formula:
            issues.append(issue(row_id, "bad-formula", repr(row.get("formula"))))
        text = " ".join(str(row.get(key, "")) for key in ("role", "argument", "proof_boundary")).lower()
        for required in ("exact", "proof_boundary"):
            if required == "proof_boundary":
                if "not" not in text and "only" not in text:
                    issues.append(issue(row_id, "weak-row-boundary", "not/only"))
            elif required not in text:
                issues.append(issue(row_id, "weak-row-boundary", required))
    return issues, len(rows)


def validate_finite_diagnostics(lemma: dict) -> tuple[list[BoundaryThresholdIssue], int, int]:
    issues: list[BoundaryThresholdIssue] = []
    recorded = lemma.get("finite_diagnostics", {})
    max_index = int(recorded.get("max_coefficient_index", 0))
    if max_index != 64:
        issues.append(issue("finite_diagnostics", "bad-max-coefficient-index", repr(max_index)))
        max_index = 64
    recomputed = asdict(finite_diagnostics(list(DEFAULT_ENCLOSURE_JSONL), max_index))
    if recorded != recomputed:
        issues.append(issue("finite_diagnostics", "recorded-diagnostics-stale", f"{recorded!r} != {recomputed!r}"))
    expected_counts = {
        "strong_threshold_rows": 315,
        "strong_threshold_positive_rows": 315,
        "heat_flow_threshold_rows": 315,
        "heat_flow_threshold_positive_rows": 315,
    }
    for key, value in expected_counts.items():
        if recorded.get(key) != value:
            issues.append(issue("finite_diagnostics", f"bad-{key}", f"{recorded.get(key)!r} != {value!r}"))
    if recorded.get("min_strong_threshold_margin", {}).get("sample") != "8.170843282044743874E-3":
        issues.append(issue("finite_diagnostics", "bad-min-strong", repr(recorded.get("min_strong_threshold_margin"))))
    if recorded.get("min_heat_flow_threshold_margin", {}).get("sample") != "3.822433850353900366E-2":
        issues.append(issue("finite_diagnostics", "bad-min-heat", repr(recorded.get("min_heat_flow_threshold_margin"))))
    return issues, int(recorded.get("strong_threshold_rows", 0)), int(recorded.get("heat_flow_threshold_rows", 0))


def validate_summary(lemma: dict) -> list[BoundaryThresholdIssue]:
    summary = lemma.get("summary", {})
    issues: list[BoundaryThresholdIssue] = []
    expected = {
        "exact_rows": 5,
        "strong_threshold_rows": 315,
        "strong_threshold_positive_rows": 315,
        "heat_flow_threshold_rows": 315,
        "heat_flow_threshold_positive_rows": 315,
        "boundary_threshold_closed": True,
        "monotone_contraction_target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("not an open subtarget", "cauchy-schwarz", "remaining heat-flow route", "initial or terminal"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in lemma.get("invariants", [])).lower()
    for required in ("threshold is closed", "monotone-contraction theorem remains open", "finite", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[BoundaryThresholdIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[BoundaryThresholdIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "monotone contractions are proved for all zeta windows",
        "analytic monotone-contraction theorem is proved",
        "closed differential inequality is proved",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(lemma_path: Path, note_path: Path) -> tuple[list[BoundaryThresholdIssue], int, int, int]:
    lemma = load_json(lemma_path)
    issues: list[BoundaryThresholdIssue] = []
    issues.extend(validate_top_level(lemma))
    exact_issues, exact_rows = validate_exact_rows(lemma)
    issues.extend(exact_issues)
    finite_issues, strong_rows, heat_rows = validate_finite_diagnostics(lemma)
    issues.extend(finite_issues)
    issues.extend(validate_summary(lemma))
    issues.extend(validate_note(note_path))
    return issues, exact_rows, strong_rows, heat_rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lemma", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, exact_rows, strong_rows, heat_rows = validate(args.lemma, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "exact_rows": exact_rows,
                    "strong_threshold_rows": strong_rows,
                    "heat_threshold_rows": heat_rows,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-HEAT-FLOW-THRESHOLD {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF heat-flow boundary threshold lemma: "
            f"{exact_rows} exact rows, {strong_rows} strong-threshold rows, "
            f"{heat_rows} heat-threshold rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
