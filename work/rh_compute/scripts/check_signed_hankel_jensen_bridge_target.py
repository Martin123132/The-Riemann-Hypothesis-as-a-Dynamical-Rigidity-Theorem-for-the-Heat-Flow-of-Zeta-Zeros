#!/usr/bin/env python3
"""Validate the signed-Hankel/Jensen bridge target specification."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_NOTE = REPO_ROOT / "outputs/signed_hankel_jensen_bridge_target.md"


@dataclass(frozen=True)
class TargetIssue:
    section: str
    issue: str
    detail: str


REQUIRED_STRINGS = (
    "Status: theorem target",
    "This is not a proof of RH or `Lambda <= 0`",
    "for every d >= 1 and n >= 0",
    "R_{k,n}(j_1,...,j_k)",
    "(-1)^(k(k-1)/2)",
    "n = 0",
    "k = 2..7",
    "689,795/689,795 finite minors positive",
    "n = 0..20",
    "k = 2..5",
    "1,322,685/1,322,685 finite minors positive",
    "k = 6",
    "840,840/840,840 finite minors positive",
    "k = 7",
    "675,675/675,675 finite minors positive",
    "All certificates remain finite",
    "Candidate Theorem B-Star",
    "Degree 2 is exact",
    "Degree 3 is the first nontrivial bridge obstruction",
    "low-order finite sign checks cannot be promoted into Jensen hyperbolicity",
    "Proof Obligations",
    "Kill Gates",
    "This theorem target is open",
)


FORBIDDEN_STRINGS = (
    "therefore RH",
    "therefore `Lambda <= 0`",
    "we have proved Lambda <= 0",
    "the bridge is proved",
)


def validate_note(path: Path) -> list[TargetIssue]:
    issues: list[TargetIssue] = []
    if not path.exists():
        return [TargetIssue("<file>", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    for required in REQUIRED_STRINGS:
        if required not in text:
            issues.append(TargetIssue("<required>", "missing-text", required))
    lowered = text.lower()
    for forbidden in FORBIDDEN_STRINGS:
        if forbidden.lower() in lowered:
            issues.append(TargetIssue("<forbidden>", "forbidden-text", forbidden))

    proof_obligation_count = sum(
        1 for line in text.splitlines() if line.startswith(("1. ", "2. ", "3. ", "4. ", "5. "))
    )
    if proof_obligation_count < 5:
        issues.append(TargetIssue("Proof Obligations", "too-few-obligations", str(proof_obligation_count)))

    if "outputs/jensen_hankel_bridge_algebra.md" not in text:
        issues.append(TargetIssue("Degree 3", "missing-algebra-note-ref", "outputs/jensen_hankel_bridge_algebra.md"))
    if "work/rh_compute/results/jensen_hankel_bridge_algebra.json" not in text:
        issues.append(TargetIssue("Degree 3", "missing-algebra-json-ref", "work/rh_compute/results/jensen_hankel_bridge_algebra.json"))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues = validate_note(args.note)
    if args.json:
        print(json.dumps({"ok": not issues, "issues": [asdict(issue) for issue in issues]}, indent=2, sort_keys=True))
    else:
        for issue in issues:
            print(f"SIGNED-HANKEL-JENSEN-TARGET {issue.section} [{issue.issue}] {issue.detail}")
        print(f"validated signed-Hankel/Jensen bridge target specification with {len(issues)} issues")
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
