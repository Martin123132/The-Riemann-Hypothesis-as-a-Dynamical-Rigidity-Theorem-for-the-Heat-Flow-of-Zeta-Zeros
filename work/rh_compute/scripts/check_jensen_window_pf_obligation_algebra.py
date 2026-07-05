#!/usr/bin/env python3
"""Validate the exact Jensen-window PF obligation algebra gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_obligation_algebra.json"


@dataclass(frozen=True)
class ObligationIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def require(condition: bool, issues: list[ObligationIssue], row_id: str, issue: str, detail: str) -> None:
    if not condition:
        issues.append(ObligationIssue(row_id, issue, detail))


def validate(payload: dict) -> list[ObligationIssue]:
    issues: list[ObligationIssue] = []
    require(
        payload.get("kind") == "jensen_window_pf_obligation_algebra",
        issues,
        "<result>",
        "bad-kind",
        repr(payload.get("kind")),
    )
    require(
        "not a proof of Jensen-window PF-infinity" in str(payload.get("proof_boundary", "")),
        issues,
        "<result>",
        "weak-proof-boundary",
        str(payload.get("proof_boundary")),
    )

    degree2 = payload.get("degree2", {})
    require(
        degree2.get("matches_signed_hankel_m1") is True,
        issues,
        "degree2",
        "missing-signed-hankel-contact",
        repr(degree2),
    )
    require(
        degree2.get("pf_threshold") == "a1**2 - a0*a2 >= 0",
        issues,
        "degree2",
        "bad-pf-threshold",
        repr(degree2.get("pf_threshold")),
    )

    degree3 = payload.get("degree3", {})
    d3_formulas = {
        row.get("determinant")
        for row in degree3.get("selected_toeplitz_minors", [])
        if isinstance(row, dict)
    }
    for expected in (
        "-3*(a0*a2 - 3*a1**2)",
        "-a0*a3 + 9*a1*a2",
        "-3*(a1*a3 - 3*a2**2)",
        "a0**2*a3 - 18*a0*a1*a2 + 27*a1**3",
        "a0*a3**2 - 18*a1*a2*a3 + 27*a2**3",
    ):
        require(expected in d3_formulas, issues, "degree3", "missing-formula", expected)
    require(
        "-27*(" in str(degree3.get("cubic_discriminant", ""))
        and "a0**2*a3**2" in str(degree3.get("cubic_discriminant", "")),
        issues,
        "degree3",
        "missing-cubic-discriminant",
        repr(degree3.get("cubic_discriminant")),
    )

    degree4 = payload.get("degree4", {})
    d4_formulas = {
        row.get("determinant")
        for row in degree4.get("selected_toeplitz_minors", [])
        if isinstance(row, dict)
    }
    for expected in (
        "-2*(3*a0*a2 - 8*a1**2)",
        "-4*(a0*a3 - 6*a1*a2)",
        "-a0*a4 + 16*a1*a3",
        "4*(a0**2*a3 - 12*a0*a1*a2 + 16*a1**3)",
        "-2*(3*a0*a2*a4 - 8*a0*a3**2 - 8*a1**2*a4 + 96*a1*a2*a3 - 108*a2**3)",
        "4*(a1*a4**2 - 12*a2*a3*a4 + 16*a3**3)",
    ):
        require(expected in d4_formulas, issues, "degree4", "missing-formula", expected)
    require(
        "a0**3*a4**3" in str(degree4.get("quartic_discriminant", "")),
        issues,
        "degree4",
        "missing-quartic-discriminant",
        repr(degree4.get("quartic_discriminant")),
    )

    counter = payload.get("finite_countermodel", {})
    require(
        counter.get("sequence_A0_to_A4") == ["1", "33/40", "429/640", "4719/12800", "4719/35840"],
        issues,
        "countermodel",
        "bad-sequence",
        repr(counter.get("sequence_A0_to_A4")),
    )
    require(
        counter.get("d3_selected_toeplitz_minors_positive") is True,
        issues,
        "countermodel",
        "d3-selected-not-positive",
        repr(counter.get("d3_selected_toeplitz_minors")),
    )
    require(
        counter.get("d3_cubic_discriminant") == "-2476526481/3276800000",
        issues,
        "countermodel",
        "bad-d3-discriminant",
        repr(counter.get("d3_cubic_discriminant")),
    )
    require(
        counter.get("d3_cubic_hyperbolicity_breaks") is True,
        issues,
        "countermodel",
        "d3-does-not-break",
        repr(counter.get("d3_cubic_hyperbolicity_breaks")),
    )
    require(
        counter.get("d3_first_negative_contiguous_toeplitz_minor", {}).get("size") == 8,
        issues,
        "countermodel",
        "bad-d3-first-negative-size",
        repr(counter.get("d3_first_negative_contiguous_toeplitz_minor")),
    )
    require(
        counter.get("d3_first_negative_contiguous_toeplitz_minor", {}).get("determinant")
        == "-435846079534239/104857600000000",
        issues,
        "countermodel",
        "bad-d3-first-negative-det",
        repr(counter.get("d3_first_negative_contiguous_toeplitz_minor")),
    )
    require(
        counter.get("d4_selected_toeplitz_minors_positive") is True,
        issues,
        "countermodel",
        "d4-selected-not-positive",
        repr(counter.get("d4_selected_toeplitz_minors")),
    )
    require(
        counter.get("d4_quartic_discriminant") == "-668519580649275927/359661568000000000",
        issues,
        "countermodel",
        "bad-d4-discriminant",
        repr(counter.get("d4_quartic_discriminant")),
    )
    require(
        counter.get("d4_quartic_hyperbolicity_breaks") is True,
        issues,
        "countermodel",
        "d4-does-not-break",
        repr(counter.get("d4_quartic_hyperbolicity_breaks")),
    )
    require(
        counter.get("d4_first_negative_contiguous_toeplitz_minor", {}).get("size") == 6,
        issues,
        "countermodel",
        "bad-d4-first-negative-size",
        repr(counter.get("d4_first_negative_contiguous_toeplitz_minor")),
    )
    require(
        counter.get("d4_first_negative_contiguous_toeplitz_minor", {}).get("determinant")
        == "-229760849637/28672000000",
        issues,
        "countermodel",
        "bad-d4-first-negative-det",
        repr(counter.get("d4_first_negative_contiguous_toeplitz_minor")),
    )
    require(
        "low-order Jensen-window Toeplitz minors" in str(counter.get("blocked_promotion", "")),
        issues,
        "countermodel",
        "missing-blocked-promotion",
        repr(counter.get("blocked_promotion")),
    )
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues = validate(load_json(args.result))
    if args.json:
        print(json.dumps({"ok": not issues, "issues": [asdict(issue) for issue in issues]}, indent=2, sort_keys=True))
    else:
        for issue in issues:
            print(f"JENSEN-WINDOW-PF-OBLIGATION {issue.row_id} [{issue.issue}] {issue.detail}")
        print(f"validated Jensen-window PF obligation algebra with {len(issues)} issues")
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
