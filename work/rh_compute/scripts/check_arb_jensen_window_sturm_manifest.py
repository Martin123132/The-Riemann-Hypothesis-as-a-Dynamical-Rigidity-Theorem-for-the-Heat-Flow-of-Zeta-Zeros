#!/usr/bin/env python3
"""Validate the Arb Jensen-window Sturm finite diagnostic manifest."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SUMMARY = REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520_summary.json"
DEFAULT_ROWS = REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520.jsonl"
DEFAULT_LAMBDAS = ["0", "1e-6", "1e-4", "1e-2", "1e-1"]
DEFAULT_SHIFTS = list(range(21))


@dataclass(frozen=True)
class ManifestIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_int_range(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if ".." in part:
            left, right = part.split("..", 1)
            start = int(left.strip())
            stop = int(right.strip())
            if stop < start:
                raise ValueError(f"descending range is not supported: {part}")
            out.extend(range(start, stop + 1))
        else:
            out.append(int(part))
    return sorted(set(out))


def validate(summary: dict, rows_path: Path, expected_degrees: list[int]) -> list[ManifestIssue]:
    issues: list[ManifestIssue] = []
    expected_rows = len(DEFAULT_LAMBDAS) * len(DEFAULT_SHIFTS) * len(expected_degrees)
    expected_rows_by_degree = {str(degree): len(DEFAULT_LAMBDAS) * len(DEFAULT_SHIFTS) for degree in expected_degrees}
    expected_max_k = max(DEFAULT_SHIFTS) + max(expected_degrees)

    if summary.get("kind") != "arb_jensen_window_sturm_hyperbolicity_summary":
        issues.append(ManifestIssue("<summary>", "bad-kind", repr(summary.get("kind"))))
    boundary = str(summary.get("proof_boundary", ""))
    if "not all-degree or all-shift Jensen hyperbolicity" not in boundary:
        issues.append(ManifestIssue("<summary>", "weak-proof-boundary", boundary))
    if summary.get("lambdas") != DEFAULT_LAMBDAS:
        issues.append(ManifestIssue("<summary>", "bad-lambdas", repr(summary.get("lambdas"))))
    if summary.get("shifts") != DEFAULT_SHIFTS:
        issues.append(ManifestIssue("<summary>", "bad-shifts", repr(summary.get("shifts"))))
    if summary.get("degrees") != expected_degrees:
        issues.append(ManifestIssue("<summary>", "bad-degrees", repr(summary.get("degrees"))))
    if summary.get("dps") != 520:
        issues.append(ManifestIssue("<summary>", "bad-dps", repr(summary.get("dps"))))
    if summary.get("needed_max_k") != expected_max_k:
        issues.append(ManifestIssue("<summary>", "bad-needed-max-k", repr(summary.get("needed_max_k"))))
    if summary.get("rows") != expected_rows:
        issues.append(ManifestIssue("<summary>", "bad-row-count", repr(summary.get("rows"))))
    if summary.get("ok") != expected_rows or summary.get("failed_or_inconclusive") != 0:
        issues.append(
            ManifestIssue(
                "<summary>",
                "not-all-ok",
                f"ok={summary.get('ok')} failed={summary.get('failed_or_inconclusive')}",
            )
        )
    if summary.get("rows_by_degree") != expected_rows_by_degree:
        issues.append(ManifestIssue("<summary>", "bad-rows-by-degree", repr(summary.get("rows_by_degree"))))
    if summary.get("ok_by_degree") != expected_rows_by_degree:
        issues.append(ManifestIssue("<summary>", "bad-ok-by-degree", repr(summary.get("ok_by_degree"))))
    if summary.get("all_ok") is not True:
        issues.append(ManifestIssue("<summary>", "all-ok-false", repr(summary.get("all_ok"))))

    row_count = 0
    bad_rows = 0
    with rows_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row_count += 1
            row = json.loads(line)
            if row.get("kind") != "arb_jensen_window_sturm_hyperbolicity_row":
                bad_rows += 1
                continue
            degree = row.get("degree_d")
            row_id = f"lambda={row.get('lam')} n={row.get('shift_n')} d={degree}"
            if degree not in expected_degrees:
                bad_rows += 1
                issues.append(ManifestIssue(row_id, "unexpected-degree", repr(degree)))
            if row.get("positive_root_count") != degree or row.get("ok") is not True:
                bad_rows += 1
                issues.append(ManifestIssue(row_id, "bad-root-count", repr(row)))
            if any(sign is None for sign in row.get("signs_at_zero", [])):
                bad_rows += 1
                issues.append(ManifestIssue(row_id, "inconclusive-zero-sign", repr(row.get("signs_at_zero"))))
            if any(sign is None for sign in row.get("signs_at_infinity", [])):
                bad_rows += 1
                issues.append(ManifestIssue(row_id, "inconclusive-infinity-sign", repr(row.get("signs_at_infinity"))))
    if row_count != expected_rows:
        issues.append(ManifestIssue("<rows>", "bad-jsonl-row-count", repr(row_count)))
    if bad_rows:
        issues.append(ManifestIssue("<rows>", "bad-jsonl-rows", repr(bad_rows)))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--expected-degrees", default="3,4")
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    expected_degrees = parse_int_range(args.expected_degrees)
    expected_rows = len(DEFAULT_LAMBDAS) * len(DEFAULT_SHIFTS) * len(expected_degrees)
    issues = validate(load_json(args.summary), args.rows, expected_degrees)
    if args.json:
        print(json.dumps({"ok": not issues, "issues": [asdict(issue) for issue in issues]}, indent=2, sort_keys=True))
    else:
        for issue in issues:
            print(f"ARB-JENSEN-WINDOW-STURM {issue.row_id} [{issue.issue}] {issue.detail}")
        print(f"validated {expected_rows} Arb Jensen-window Sturm hyperbolicity finite diagnostics with {len(issues)} issues")
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
