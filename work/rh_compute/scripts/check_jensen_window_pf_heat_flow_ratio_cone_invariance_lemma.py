#!/usr/bin/env python3
"""Validate the Jensen-window PF heat-flow ratio-cone invariance lemma."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

from jensen_window_pf_heat_flow_monotone_closure_scout import DEFAULT_ENCLOSURE_JSONL  # noqa: E402
from jensen_window_pf_heat_flow_ratio_cone_invariance_lemma import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    finite_diagnostics,
)


EXPECTED_FORMULAS = {
    "hfrci_01_ratio_flow_identity": "dlog(x_k)/dlambda = 2*r_k*((2*k+3)*x_{k+1} + (2*k-1)/x_k - 2*(2*k+1))",
    "hfrci_02_lower_wall_inward": "at x_k=(2*k-1)/(2*k+1), F_k=(2*k+3)*x_{k+1}-(2*k+1)>=0 if x_{k+1}>=(2*k+1)/(2*k+3)",
    "hfrci_03_upper_wall_inward": "at x_k=1, F_k=(2*k+3)*(x_{k+1}-1)<=0 if x_{k+1}<=1",
    "hfrci_04_monotone_wall_inward": "at x_{k+1}=x_k=q, dlog(x_{k+1}/x_k)/dlambda >= 0 if x_{k+2}>=q and q>=(2*k-1)/(2*k+5)",
    "hfrci_05_finite_collar_requirement": "checking x_1..x_K requires collar variables x_{K+1} for lower/upper walls and x_{K+2} for monotone walls",
    "hfrci_06_conditional_forward_invariance": "if the full ratio cone holds at lambda0 and the positive heat-flow ratio ODE is well posed, then the cone is forward-invariant while the hypotheses persist",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Heat-Flow Ratio-Cone Invariance Lemma",
    "Status: exact conditional ratio-cone invariance lemma",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_heat_flow_ratio_cone_invariance_lemma`",
    "work/rh_compute/results/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.json",
    "python work/rh_compute/scripts/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.py",
    "validated Jensen-window PF heat-flow ratio cone invariance lemma: 6 exact rows, 315 lower rows, 315 upper rows, 310 monotone rows, 0 issues",
    "a_k=(2*k-1)/(2*k+1)",
    "a_k <= x_k <= 1",
    "x_{k+1} >= x_k",
    "F_k=(2*k+3)*x_{k+1} + (2*k-1)/x_k - 2*(2*k+1)",
    "A finite prefix",
    "needs collar variables",
    "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
    "validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues",
    "8.170843282044743874E-3",
    "7.573304754109889330E-3",
    "9.895636183563448458E-5",
    "1.037697096819682472E+0",
    "2.951966135520726476E-1",
    "2.665474275226736438E-4",
)


@dataclass(frozen=True)
class RatioConeIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> RatioConeIssue:
    return RatioConeIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[RatioConeIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(lemma: dict) -> list[RatioConeIssue]:
    issues: list[RatioConeIssue] = []
    if lemma.get("kind") != "jensen_window_pf_heat_flow_ratio_cone_invariance_lemma":
        issues.append(issue("<lemma>", "bad-kind", repr(lemma.get("kind"))))
    if lemma.get("status") != "available_exact_conditional_invariance_lemma":
        issues.append(issue("<lemma>", "bad-status", repr(lemma.get("status"))))
    for key in ("source_boundary_threshold", "source_heat_flow_closure", "source_theorem_target"):
        issues.extend(validate_ref("<lemma>", lemma.get(key)))
    for ref in lemma.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<lemma>", ref))
    boundary = str(lemma.get("proof_boundary", "")).lower()
    for required in ("exact conditional", "inward-pointing", "collared finite", "does not prove", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<lemma>", "weak-proof-boundary", required))
    return issues


def validate_exact_rows(lemma: dict) -> tuple[list[RatioConeIssue], int]:
    rows = lemma.get("exact_rows", [])
    issues: list[RatioConeIssue] = []
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
        if "only" not in text and "conditional" not in text:
            issues.append(issue(row_id, "weak-row-boundary", text))
    return issues, len(rows)


def validate_finite_diagnostics(lemma: dict) -> tuple[list[RatioConeIssue], int, int, int]:
    issues: list[RatioConeIssue] = []
    recorded = lemma.get("finite_diagnostics", {})
    max_index = int(recorded.get("max_coefficient_index", 0))
    if max_index != 64:
        issues.append(issue("finite_diagnostics", "bad-max-coefficient-index", repr(max_index)))
        max_index = 64
    recomputed = asdict(finite_diagnostics(list(DEFAULT_ENCLOSURE_JSONL), max_index))
    if recorded != recomputed:
        issues.append(issue("finite_diagnostics", "recorded-diagnostics-stale", f"{recorded!r} != {recomputed!r}"))
    expected_counts = {
        "coordinate_lower_rows": 315,
        "coordinate_lower_positive_rows": 315,
        "coordinate_upper_rows": 315,
        "coordinate_upper_positive_rows": 315,
        "coordinate_monotone_rows": 310,
        "coordinate_monotone_positive_rows": 310,
        "lower_wall_rows": 310,
        "lower_wall_positive_rows": 310,
        "upper_wall_rows": 310,
        "upper_wall_positive_rows": 310,
        "monotone_wall_rows": 305,
        "monotone_wall_positive_rows": 305,
    }
    for key, value in expected_counts.items():
        if recorded.get(key) != value:
            issues.append(issue("finite_diagnostics", f"bad-{key}", f"{recorded.get(key)!r} != {value!r}"))
    expected_mins = {
        "min_coordinate_lower_margin": "8.170843282044743874E-3",
        "min_coordinate_upper_margin": "7.573304754109889330E-3",
        "min_coordinate_monotone_margin": "9.895636183563448458E-5",
        "min_lower_wall_margin": "1.037697096819682472E+0",
        "min_upper_wall_margin": "2.951966135520726476E-1",
        "min_monotone_wall_margin": "2.665474275226736438E-4",
    }
    for key, value in expected_mins.items():
        if recorded.get(key, {}).get("sample") != value:
            issues.append(issue("finite_diagnostics", f"bad-{key}", repr(recorded.get(key))))
    return (
        issues,
        int(recorded.get("coordinate_lower_rows", 0)),
        int(recorded.get("coordinate_upper_rows", 0)),
        int(recorded.get("coordinate_monotone_rows", 0)),
    )


def validate_summary(lemma: dict) -> list[RatioConeIssue]:
    summary = lemma.get("summary", {})
    issues: list[RatioConeIssue] = []
    expected = {
        "exact_rows": 6,
        "coordinate_lower_rows": 315,
        "coordinate_upper_rows": 315,
        "coordinate_monotone_rows": 310,
        "lower_wall_rows": 310,
        "upper_wall_rows": 310,
        "monotone_wall_rows": 305,
        "conditional_invariance_available": True,
        "zeta_initial_condition_available": False,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("conditional", "forward-invariance", "remaining theorem gap", "enters the full infinite cone"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in lemma.get("invariants", [])).lower()
    for required in ("conditional", "finite zeta", "no initial", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[RatioConeIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[RatioConeIssue] = []
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


def validate(lemma_path: Path, note_path: Path) -> tuple[list[RatioConeIssue], int, int, int, int]:
    lemma = load_json(lemma_path)
    issues: list[RatioConeIssue] = []
    issues.extend(validate_top_level(lemma))
    exact_issues, exact_rows = validate_exact_rows(lemma)
    issues.extend(exact_issues)
    finite_issues, lower_rows, upper_rows, monotone_rows = validate_finite_diagnostics(lemma)
    issues.extend(finite_issues)
    issues.extend(validate_summary(lemma))
    issues.extend(validate_note(note_path))
    return issues, exact_rows, lower_rows, upper_rows, monotone_rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lemma", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, exact_rows, lower_rows, upper_rows, monotone_rows = validate(args.lemma, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "exact_rows": exact_rows,
                    "lower_rows": lower_rows,
                    "upper_rows": upper_rows,
                    "monotone_rows": monotone_rows,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-HEAT-FLOW-RATIO-CONE {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF heat-flow ratio cone invariance lemma: "
            f"{exact_rows} exact rows, {lower_rows} lower rows, "
            f"{upper_rows} upper rows, {monotone_rows} monotone rows, "
            f"{len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
