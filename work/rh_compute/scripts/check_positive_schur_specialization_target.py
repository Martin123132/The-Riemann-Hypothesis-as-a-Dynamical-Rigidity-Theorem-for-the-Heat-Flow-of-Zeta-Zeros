#!/usr/bin/env python3
"""Validate the positive Schur-specialization theorem target note."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_NOTE = REPO_ROOT / "outputs/positive_schur_specialization_target.md"


@dataclass(frozen=True)
class SchurTargetIssue:
    section: str
    issue: str
    detail: str


REQUIRED_STRINGS = (
    "Status: theorem target note",
    "This is not a proof of PF-infinity, Laguerre-Polya membership, RH, or `Lambda <= 0`.",
    "construct a positive Schur/Edrei-Thoma specialization h_k -> d_k(0)",
    "c_k(lambda) = A_k(lambda) / k! = mu_{2k}(lambda) / (2k)!",
    "d_k(lambda) = c_k(lambda) / c_0(lambda)",
    "c_0^m s_{lambda/mu}(d)",
    "h_k -> d_k",
    "python work/rh_compute/scripts/check_toeplitz_jacobi_trudi_map.py",
    "validated Toeplitz/Jacobi-Trudi reindexing: N=10, orders<=5, 124129 minors, 39094 nonzero maps, 85035 structural zeros",
    "This validates the translation. It does not prove positivity",
    "S1: Positive Specialization",
    "h_k(alpha) = d_k(0) for every k >= 0",
    "S2: Planar Network",
    "P_{r,q} = d_{q-r}(0)",
    "S3: Production Matrix Or Continued Fraction",
    "S4: Positive Determinant Integral",
    "Evidence That May Guide The Search",
    "These artifacts are evidence and falsification pressure. They are not a proof",
    "Kill Gates",
    "uses only finitely many Toeplitz, Schur, Hankel, or Edrei-log checks",
    "invokes ASW/Edrei after assuming PF-infinity or real-negative zeros",
    "proves only Stieltjes or Hankel moment positivity",
    "assumes Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0",
    "leaves a signed determinant or signed integral with no independent positivity theorem",
    "python work/rh_compute/scripts/countermodel_gate_examples.py",
    "No positive specialization, planar network, production matrix, continued",
    "Target S remains an open theorem target.",
)


FORBIDDEN_STRINGS = (
    "therefore RH",
    "therefore Lambda <= 0",
    "we have proved PF-infinity",
    "we have proved RH",
    "Target S is proved",
    "finite certificates prove Target S",
    "finite evidence proves Target S",
)


REQUIRED_REFS = (
    "outputs/toeplitz_jacobi_trudi_bridge_note.md",
    "outputs/finite_toeplitz_certificate_status.md",
    "outputs/edrei_log_sign_diagnostic.md",
    "outputs/missing_coefficient_pf_theorem.md",
    "outputs/sign_regularity_theorem_fit_matrix.md",
    "outputs/countermodel_library.md",
    "work/rh_compute/scripts/check_toeplitz_certificate_manifest.py",
    "work/rh_compute/scripts/check_toeplitz_jacobi_trudi_map.py",
    "work/rh_compute/scripts/check_edrei_log_sign_manifest.py",
    "work/rh_compute/scripts/check_edrei_power_hankel_manifest.py",
    "work/rh_compute/scripts/countermodel_gate_examples.py",
)


def validate(path: Path) -> list[SchurTargetIssue]:
    if not path.exists():
        return [SchurTargetIssue("<file>", "missing-note", str(path))]

    text = path.read_text(encoding="utf-8")
    issues: list[SchurTargetIssue] = []

    for required in REQUIRED_STRINGS:
        if required not in text:
            issues.append(SchurTargetIssue("<required>", "missing-text", required))

    lowered = text.lower()
    for forbidden in FORBIDDEN_STRINGS:
        if forbidden.lower() in lowered:
            issues.append(SchurTargetIssue("<forbidden>", "forbidden-text", forbidden))

    for ref in REQUIRED_REFS:
        if ref not in text:
            issues.append(SchurTargetIssue("references", "missing-note-ref", ref))
        elif not (REPO_ROOT / ref).exists():
            issues.append(SchurTargetIssue("references", "missing-path", ref))

    sufficient_count = sum(1 for line in text.splitlines() if line.startswith("### S"))
    if sufficient_count < 4:
        issues.append(SchurTargetIssue("Target S", "too-few-sufficient-routes", str(sufficient_count)))

    kill_count = sum(1 for line in text.splitlines() if line.strip().startswith(tuple(f"{i}. " for i in range(1, 8))))
    if kill_count < 7:
        issues.append(SchurTargetIssue("Kill Gates", "too-few-kill-gates", str(kill_count)))

    acceptance_count = sum(1 for line in text.splitlines() if line.strip().startswith(tuple(f"{i}. " for i in range(1, 6))))
    if acceptance_count < 5:
        issues.append(SchurTargetIssue("Acceptance Condition", "too-few-acceptance-items", str(acceptance_count)))

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
            print(f"SCHUR-TARGET {issue.section} [{issue.issue}] {issue.detail}")
        print(f"validated positive Schur-specialization target note with {len(issues)} issues")
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
