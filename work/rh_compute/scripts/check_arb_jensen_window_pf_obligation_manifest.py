#!/usr/bin/env python3
"""Validate the Arb Jensen-window PF obligation finite diagnostic manifest."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/arb_jensen_window_pf_obligations_lamgrid_n0_n20_d3_m1_m8_d4_m1_m6_dps520_summary.json"
)
DEFAULT_ROWS = (
    REPO_ROOT
    / "work/rh_compute/results/arb_jensen_window_pf_obligations_lamgrid_n0_n20_d3_m1_m8_d4_m1_m6_dps520.jsonl"
)


@dataclass(frozen=True)
class ManifestIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate(summary: dict, rows_path: Path) -> list[ManifestIssue]:
    issues: list[ManifestIssue] = []
    if summary.get("kind") != "arb_jensen_window_pf_obligation_summary":
        issues.append(ManifestIssue("<summary>", "bad-kind", repr(summary.get("kind"))))
    boundary = str(summary.get("proof_boundary", ""))
    if "not all-minor Jensen-window PF-infinity" not in boundary:
        issues.append(ManifestIssue("<summary>", "weak-proof-boundary", boundary))
    if summary.get("lambdas") != ["0", "1e-6", "1e-4", "1e-2", "1e-1"]:
        issues.append(ManifestIssue("<summary>", "bad-lambdas", repr(summary.get("lambdas"))))
    if summary.get("shifts") != list(range(21)):
        issues.append(ManifestIssue("<summary>", "bad-shifts", repr(summary.get("shifts"))))
    if summary.get("degree_sizes") != {"3": list(range(1, 9)), "4": list(range(1, 7))}:
        issues.append(ManifestIssue("<summary>", "bad-degree-sizes", repr(summary.get("degree_sizes"))))
    if summary.get("dps") != 520:
        issues.append(ManifestIssue("<summary>", "bad-dps", repr(summary.get("dps"))))
    if summary.get("needed_max_k") != 24:
        issues.append(ManifestIssue("<summary>", "bad-needed-max-k", repr(summary.get("needed_max_k"))))
    if summary.get("summary_rows") != 210:
        issues.append(ManifestIssue("<summary>", "bad-summary-row-count", repr(summary.get("summary_rows"))))
    if summary.get("rows") != 1470:
        issues.append(ManifestIssue("<summary>", "bad-row-count", repr(summary.get("rows"))))
    if summary.get("ok") != 1470 or summary.get("failed_or_inconclusive") != 0:
        issues.append(
            ManifestIssue(
                "<summary>",
                "not-all-ok",
                f"ok={summary.get('ok')} failed={summary.get('failed_or_inconclusive')}",
            )
        )
    if summary.get("all_ok") is not True:
        issues.append(ManifestIssue("<summary>", "all-ok-false", repr(summary.get("all_ok"))))

    summaries = summary.get("summaries", [])
    if len(summaries) != 210:
        issues.append(ManifestIssue("<summary>", "bad-subsummary-count", repr(len(summaries))))
    for sub in summaries:
        row_id = f"lambda={sub.get('lam')} n={sub.get('shift_n')} d={sub.get('degree_d')}"
        if sub.get("all_ok") is not True:
            issues.append(ManifestIssue(row_id, "subsummary-not-ok", repr(sub)))
        degree = sub.get("degree_d")
        expected_tests = 8 if degree == 3 else 6 if degree == 4 else None
        if sub.get("tests") != expected_tests or sub.get("positive") != expected_tests:
            issues.append(ManifestIssue(row_id, "bad-subsummary-counts", repr(sub)))

    row_count = 0
    bad_rows = 0
    with rows_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row_count += 1
            row = json.loads(line)
            if row.get("kind") != "arb_jensen_window_pf_obligation_row":
                bad_rows += 1
                continue
            if row.get("classification") != "positive" or row.get("contains_zero") or row.get("ok") is not True:
                bad_rows += 1
    if row_count != 1470:
        issues.append(ManifestIssue("<rows>", "bad-jsonl-row-count", repr(row_count)))
    if bad_rows:
        issues.append(ManifestIssue("<rows>", "bad-jsonl-rows", repr(bad_rows)))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues = validate(load_json(args.summary), args.rows)
    if args.json:
        print(json.dumps({"ok": not issues, "issues": [asdict(issue) for issue in issues]}, indent=2, sort_keys=True))
    else:
        for issue in issues:
            print(f"ARB-JENSEN-WINDOW-PF {issue.row_id} [{issue.issue}] {issue.detail}")
        print(f"validated 1470 Arb Jensen-window PF obligation finite diagnostics with {len(issues)} issues")
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
