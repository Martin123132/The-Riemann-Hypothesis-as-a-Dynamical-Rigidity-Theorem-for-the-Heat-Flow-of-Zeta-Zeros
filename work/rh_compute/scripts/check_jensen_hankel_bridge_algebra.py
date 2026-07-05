#!/usr/bin/env python3
"""Validate the exact Jensen/Hankel bridge algebra countermodel gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results/jensen_hankel_bridge_algebra.json"


@dataclass(frozen=True)
class AlgebraIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate(payload: dict) -> list[AlgebraIssue]:
    issues: list[AlgebraIssue] = []
    if payload.get("kind") != "jensen_hankel_bridge_algebra":
        issues.append(AlgebraIssue("<result>", "bad-kind", repr(payload.get("kind"))))
    if "not an all-order sign-consistency theorem" not in str(payload.get("proof_boundary", "")):
        issues.append(AlgebraIssue("<result>", "weak-proof-boundary", str(payload.get("proof_boundary"))))

    degree2 = payload.get("degree2_identity", {})
    if degree2.get("exact_identity_verified") is not True:
        issues.append(AlgebraIssue("degree2", "identity-not-verified", repr(degree2)))
    if "-4*det" not in str(degree2.get("hankel_relation", "")):
        issues.append(AlgebraIssue("degree2", "missing-hankel-relation", repr(degree2.get("hankel_relation"))))

    counter = payload.get("finite_countermodel", {})
    expected_sequence = ["1", "33/40", "429/640", "4719/12800", "4719/35840"]
    if counter.get("sequence_A0_to_A4") != expected_sequence:
        issues.append(
            AlgebraIssue(
                "countermodel",
                "bad-sequence",
                f"{counter.get('sequence_A0_to_A4')!r} != {expected_sequence!r}",
            )
        )
    if counter.get("all_coefficients_positive") is not True:
        issues.append(AlgebraIssue("countermodel", "nonpositive-coefficient", repr(counter.get("all_coefficients_positive"))))
    if counter.get("finite_reshaped_signs_pass") is not True:
        issues.append(AlgebraIssue("countermodel", "reshaped-signs-do-not-pass", repr(counter.get("finite_reshaped_signs_pass"))))
    if counter.get("degree3_jensen_discriminant") != "-2476526481/3276800000":
        issues.append(AlgebraIssue("countermodel", "bad-cubic-discriminant", repr(counter.get("degree3_jensen_discriminant"))))
    if counter.get("degree3_jensen_hyperbolicity_breaks") is not True:
        issues.append(AlgebraIssue("countermodel", "cubic-does-not-break", repr(counter.get("degree3_jensen_hyperbolicity_breaks"))))

    k2_rows = counter.get("reshaped_k2_minors_N3", [])
    if len(k2_rows) != 3:
        issues.append(AlgebraIssue("countermodel", "bad-k2-row-count", repr(len(k2_rows))))
    for row in k2_rows:
        row_id = f"k2:{row.get('columns')}"
        if row.get("expected_sign") != -1 or row.get("signed_positive") is not True:
            issues.append(AlgebraIssue(row_id, "bad-k2-sign", repr(row)))

    k3_rows = counter.get("reshaped_k3_minors_N3", [])
    if len(k3_rows) != 1:
        issues.append(AlgebraIssue("countermodel", "bad-k3-row-count", repr(len(k3_rows))))
    for row in k3_rows:
        row_id = f"k3:{row.get('columns')}"
        if row.get("determinant") != "-281710143/9175040000":
            issues.append(AlgebraIssue(row_id, "bad-k3-det", repr(row.get("determinant"))))
        if row.get("expected_sign") != -1 or row.get("signed_positive") is not True:
            issues.append(AlgebraIssue(row_id, "bad-k3-sign", repr(row)))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = load_json(args.result)
    issues = validate(payload)
    if args.json:
        print(json.dumps({"ok": not issues, "issues": [asdict(issue) for issue in issues]}, indent=2, sort_keys=True))
    else:
        for issue in issues:
            print(f"JENSEN-HANKEL-ALGEBRA {issue.row_id} [{issue.issue}] {issue.detail}")
        print(
            "validated Jensen/Hankel bridge algebra gate: "
            "degree2 identity and degree3 finite countermodel "
            f"with {len(issues)} issues"
        )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
