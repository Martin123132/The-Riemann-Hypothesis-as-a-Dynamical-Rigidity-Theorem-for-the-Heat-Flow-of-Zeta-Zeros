#!/usr/bin/env python3
"""Validate the Jensen-window PF bridge target specification."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_bridge_target.md"


@dataclass(frozen=True)
class TargetIssue:
    section: str
    issue: str
    detail: str


REQUIRED_STRINGS = (
    "Status: theorem target artifact",
    "This is not a proof of PF-infinity",
    "B^{d,n,lambda}_j = binom(d,j) A_{n+j}(lambda)",
    "for every d >= 1 and n >= 0",
    "the finite sequence B^{d,n,0}_0,...,B^{d,n,0}_d is PF-infinity",
    "Toeplitz total-positivity language",
    "infinite banded Toeplitz matrix",
    "T_{d,n}(r,c) = B^{d,n,0}_{c-r}",
    "Exact Relation To Jensen",
    "Contact With Signed Hankel",
    "d = 2",
    "A_{n+1}^2 - A_n A_{n+2} >= 0",
    "m = 1",
    "For `d >= 3`",
    "outputs/jensen_window_pf_obligation_algebra.md",
    "work/rh_compute/results/jensen_window_pf_obligation_algebra.json",
    "python work/rh_compute/scripts/jensen_window_pf_obligation_algebra.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_obligation_algebra.py",
    "selected degree-3 and degree-4",
    "outputs/jensen_hankel_bridge_algebra.md",
    "work/rh_compute/results/jensen_hankel_bridge_algebra.json",
    "python work/rh_compute/scripts/check_jensen_hankel_bridge_algebra.py",
    "Required Bridge Theorem",
    "all-order shifted reshaped-Hankel sign-consistency for A_k(0)",
    "every binomially weighted Jensen window B^{d,n,0} is PF-infinity",
    "python work/rh_compute/scripts/check_arb_shifted_hankel_staircase_manifest.py",
    "3,154,515/3,154,515 finite shifted reshaped-Hankel minors",
    "c_k(lambda) = A_k(lambda)/k!",
    "Kill Gates",
    "The Jensen-window PF bridge is open",
)


FORBIDDEN_STRINGS = (
    "therefore RH",
    "therefore `Lambda <= 0`",
    "we have proved Lambda <= 0",
    "the bridge is proved",
    "This proves the Jensen-window PF bridge",
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

    kill_gate_count = 0
    in_kill_gates = False
    for line in text.splitlines():
        if line.startswith("## Kill Gates"):
            in_kill_gates = True
            continue
        if in_kill_gates and line.startswith("## "):
            break
        if in_kill_gates and line.strip().endswith(";"):
            kill_gate_count += 1
    if kill_gate_count < 6:
        issues.append(TargetIssue("Kill Gates", "too-few-kill-gates", str(kill_gate_count)))

    required_refs = (
        "outputs/signed_hankel_jensen_bridge_target.md",
        "outputs/sign_regularity_theorem_fit_matrix.md",
        "work/rh_compute/scripts/check_signed_hankel_jensen_bridge_target.py",
        "work/rh_compute/scripts/check_sign_regularity_theorem_fit_matrix.py",
    )
    for ref in required_refs:
        if ref not in text:
            issues.append(TargetIssue("Integration Points", "missing-ref", ref))
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
            print(f"JENSEN-WINDOW-PF-TARGET {issue.section} [{issue.issue}] {issue.detail}")
        print(f"validated Jensen-window PF bridge target with {len(issues)} issues")
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
