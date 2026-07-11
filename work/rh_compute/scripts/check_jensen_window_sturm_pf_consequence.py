#!/usr/bin/env python3
"""Validate the finite Sturm-to-PF Jensen-window consequence note.

This gate records a finite consequence of the promoted Arb/Sturm diagnostics:
for the checked windows only, the certified positive-root count for
Q(y)=P(-y) gives real nonpositive zeros for P.  By the finite
Polya-frequency characterization, those checked binomial windows are finite
PF-infinity sequences.  This is not an all-degree/all-shift theorem.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_sturm_pf_consequence.md"


@dataclass(frozen=True)
class ManifestSpec:
    summary: Path
    rows: Path
    degrees: tuple[int, ...]
    expected_rows: int
    needed_max_k: int


SPECS = (
    ManifestSpec(
        summary=REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520_summary.json",
        rows=REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520.jsonl",
        degrees=(3, 4),
        expected_rows=210,
        needed_max_k=24,
    ),
    ManifestSpec(
        summary=REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520_summary.json",
        rows=REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520.jsonl",
        degrees=(5,),
        expected_rows=105,
        needed_max_k=25,
    ),
    ManifestSpec(
        summary=REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d6_d12_dps520_summary.json",
        rows=REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d6_d12_dps520.jsonl",
        degrees=(6, 7, 8, 9, 10, 11, 12),
        expected_rows=735,
        needed_max_k=32,
    ),
)


REQUIRED_NOTE_STRINGS = (
    "Status: finite consequence artifact",
    "This is not a proof of all-degree or all-shift Jensen hyperbolicity",
    "Q_{d,n,lambda}(y)=P_{d,n,lambda}(-y)",
    "finite Polya-frequency characterization",
    "1050/1050 checked Jensen windows",
    "degree d = 3..12",
    "lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}",
    "shifts n = 0..20",
    "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520_summary.json",
    "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520_summary.json",
    "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d6_d12_dps520_summary.json",
    "python work/rh_compute/scripts/check_jensen_window_sturm_pf_consequence.py",
    "not all-minor Jensen-window PF-infinity",
)


@dataclass(frozen=True)
class ConsequenceIssue:
    section: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_note(path: Path) -> list[ConsequenceIssue]:
    if not path.exists():
        return [ConsequenceIssue("<note>", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ConsequenceIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(ConsequenceIssue("<note>", "missing-text", required))
    forbidden = (
        "This proves the Jensen-window PF bridge",
        "therefore RH",
        "therefore `Lambda <= 0`",
        "we have proved Lambda <= 0",
    )
    lowered = text.lower()
    for phrase in forbidden:
        if phrase.lower() in lowered:
            issues.append(ConsequenceIssue("<note>", "forbidden-text", phrase))
    return issues


def validate_spec(spec: ManifestSpec) -> tuple[list[ConsequenceIssue], int]:
    issues: list[ConsequenceIssue] = []
    summary = load_json(spec.summary)
    if summary.get("kind") != "arb_jensen_window_sturm_hyperbolicity_summary":
        issues.append(ConsequenceIssue(str(spec.summary), "bad-kind", repr(summary.get("kind"))))
    if tuple(summary.get("degrees", ())) != spec.degrees:
        issues.append(ConsequenceIssue(str(spec.summary), "bad-degrees", repr(summary.get("degrees"))))
    if summary.get("rows") != spec.expected_rows:
        issues.append(ConsequenceIssue(str(spec.summary), "bad-rows", repr(summary.get("rows"))))
    if summary.get("ok") != spec.expected_rows or summary.get("failed_or_inconclusive") != 0:
        issues.append(
            ConsequenceIssue(
                str(spec.summary),
                "not-all-ok",
                f"ok={summary.get('ok')} failed={summary.get('failed_or_inconclusive')}",
            )
        )
    if summary.get("needed_max_k") != spec.needed_max_k:
        issues.append(ConsequenceIssue(str(spec.summary), "bad-needed-max-k", repr(summary.get("needed_max_k"))))
    if "not all-degree or all-shift Jensen hyperbolicity" not in str(summary.get("proof_boundary", "")):
        issues.append(ConsequenceIssue(str(spec.summary), "weak-boundary", str(summary.get("proof_boundary", ""))))

    row_count = 0
    with spec.rows.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            row_count += 1
            row = json.loads(raw)
            row_id = f"lambda={row.get('lam')} n={row.get('shift_n')} d={row.get('degree_d')}"
            degree = row.get("degree_d")
            if degree not in spec.degrees:
                issues.append(ConsequenceIssue(row_id, "unexpected-degree", repr(degree)))
            if row.get("ok") is not True or row.get("positive_root_count") != degree:
                issues.append(ConsequenceIssue(row_id, "bad-positive-root-count", repr(row)))
            if row.get("signs_at_zero", [None])[0] != 1:
                issues.append(ConsequenceIssue(row_id, "nonpositive-constant-term", repr(row.get("signs_at_zero"))))
            if any(sign is None for sign in row.get("signs_at_zero", [])):
                issues.append(ConsequenceIssue(row_id, "inconclusive-zero-sign", repr(row.get("signs_at_zero"))))
            if any(sign is None for sign in row.get("signs_at_infinity", [])):
                issues.append(ConsequenceIssue(row_id, "inconclusive-infinity-sign", repr(row.get("signs_at_infinity"))))
    if row_count != spec.expected_rows:
        issues.append(ConsequenceIssue(str(spec.rows), "bad-row-count", repr(row_count)))
    return issues, row_count


def validate(note: Path) -> list[ConsequenceIssue]:
    issues = validate_note(note)
    total_rows = 0
    for spec in SPECS:
        spec_issues, rows = validate_spec(spec)
        issues.extend(spec_issues)
        total_rows += rows
    if total_rows != 1050:
        issues.append(ConsequenceIssue("<total>", "bad-total-window-count", repr(total_rows)))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues = validate(args.note)
    if args.json:
        print(json.dumps({"ok": not issues, "issues": [asdict(issue) for issue in issues]}, indent=2, sort_keys=True))
    else:
        for issue in issues:
            print(f"JENSEN-STURM-PF {issue.section} [{issue.issue}] {issue.detail}")
        print(f"validated 1050 finite Sturm-to-PF Jensen-window consequences with {len(issues)} issues")
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
