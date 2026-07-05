#!/usr/bin/env python3
"""Validate the sign-regularity theorem fit/misfit matrix.

This is a theorem-search hygiene gate.  It checks that the fit matrix keeps
known total-positivity, PF, Jensen, Hankel, and signed-Hankel routes separated
and that finite evidence is not presented as an all-order bridge theorem.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_NOTE = REPO_ROOT / "outputs/sign_regularity_theorem_fit_matrix.md"


@dataclass(frozen=True)
class FitMatrixIssue:
    section: str
    issue: str
    detail: str


REQUIRED_STRINGS = (
    "Status: theorem-search artifact",
    "This is not a proof of PF-infinity, Laguerre-Polya membership, RH, or `Lambda <= 0`.",
    "ASW/Edrei identifies the exact all-order condition needed",
    "The finite Toeplitz certificate ledger proves ASW/Edrei hypotheses.",
    "Edrei-Thoma / Schur-Positive Specializations",
    "direct kernel PF function",
    "!= coefficient PF sequence c_k",
    "!= signed Hankel sequence A_k",
    "The criterion needs all degrees and shifts.",
    "Sign-regular Hankel data alone does not imply Toeplitz PF-infinity",
    "not a Level-4 interval",
    "certificate and not an all-order sign-consistency proof",
    "work/rh_compute/results/arb_shifted_hankel_sign_consistency_lamgrid_n0_n20_k2_k5_N18_dps520_summary.json",
    "work/rh_compute/results/arb_shifted_hankel_sign_consistency_lamgrid_n0_n20_k2_k5_N18_dps520.jsonl.gz",
    "work/rh_compute/results/arb_shifted_hankel_sign_consistency_lamgrid_n0_n20_k6_N16_dps520_summary.json",
    "work/rh_compute/results/arb_shifted_hankel_sign_consistency_lamgrid_n0_n20_k6_N16_dps520.jsonl.gz",
    "Using rigorous `A_ball` coefficient enclosures through `A_41`",
    "n = 0..20",
    "1,322,685/1,322,685 finite shifted minors positive and separated from zero",
    "The next-order shifted slice checks:",
    "needed max coefficient index = 40",
    "840,840/840,840 finite shifted minors positive and separated from zero",
    "degree 2:",
    "degree 3:",
    "low-order reshaped-Hankel signs do not supply the missing all-order bridge",
    "Current Best Theorem Targets",
    "Target 1: Positive Schur Specialization",
    "Target 2: Edrei Log-Power Representation",
    "Target 3: Determinantal Integral Formula",
    "Target 4: Signed-Hankel/Jensen Bridge",
    "Kill Gates",
    "uses only finitely many coefficients",
    "treats a finite shifted-principal Hankel grid as all-order sign consistency",
    "No checked theorem currently closes the signed-Hankel/Jensen route.",
)


FORBIDDEN_OUTSIDE_FENCES = (
    "we have proved RH",
    "therefore RH",
    "we have proved Lambda <= 0",
    "the bridge is proved",
    "finite evidence proves",
    "finite certificates prove",
    "the finite Toeplitz certificate ledger proves ASW/Edrei hypotheses",
    "the shifted Arb certificate proves all-shift",
)


def read_text(path: Path) -> str:
    with path.open("r", encoding="utf-8") as handle:
        return handle.read()


def unfenced_lines(text: str) -> list[tuple[int, str]]:
    lines: list[tuple[int, str]] = []
    in_fence = False
    for lineno, line in enumerate(text.splitlines(), start=1):
        if line.strip().startswith("```"):
            in_fence = not in_fence
            continue
        if not in_fence:
            lines.append((lineno, line))
    return lines


def validate(path: Path) -> list[FitMatrixIssue]:
    if not path.exists():
        return [FitMatrixIssue("<file>", "missing-note", str(path))]

    text = read_text(path)
    issues: list[FitMatrixIssue] = []

    for required in REQUIRED_STRINGS:
        if required not in text:
            issues.append(FitMatrixIssue("<required>", "missing-text", required))

    lowered_unfenced = [(lineno, line.lower()) for lineno, line in unfenced_lines(text)]
    for forbidden in FORBIDDEN_OUTSIDE_FENCES:
        needle = forbidden.lower()
        for lineno, line in lowered_unfenced:
            if needle in line:
                issues.append(FitMatrixIssue(f"line {lineno}", "forbidden-outside-fence", forbidden))

    target_count = sum(1 for _, line in unfenced_lines(text) if line.startswith("### Target "))
    if target_count < 4:
        issues.append(FitMatrixIssue("Current Best Theorem Targets", "too-few-targets", str(target_count)))

    kill_gate_count = sum(1 for line in text.splitlines() if line.strip().startswith(tuple(f"{i}. " for i in range(1, 8))))
    if kill_gate_count < 7:
        issues.append(FitMatrixIssue("Kill Gates", "too-few-kill-gates", str(kill_gate_count)))

    for ref in (
        "work/rh_compute/scripts/countermodel_gate_examples.py",
        "work/rh_compute/scripts/check_jensen_hankel_bridge_algebra.py",
        "work/rh_compute/scripts/check_signed_hankel_jensen_bridge_target.py",
    ):
        if not (REPO_ROOT / ref).exists():
            issues.append(FitMatrixIssue("references", "missing-ref", ref))

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
            print(f"FIT-MATRIX {issue.section} [{issue.issue}] {issue.detail}")
        print(f"validated sign-regularity theorem fit matrix with {len(issues)} issues")
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
