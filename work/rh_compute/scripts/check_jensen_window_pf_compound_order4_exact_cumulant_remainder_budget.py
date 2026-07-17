#!/usr/bin/env python3
"""Validate the order-four exact-cumulant remainder budget."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from jensen_window_pf_compound_order4_exact_cumulant_remainder_budget import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
    sha256,
)


EXPECTED_SUMMARY = (
    "validated order-four exact cumulant remainder budget: 9 rows, 0 issues, "
    "8 exact rows, finite epsilon-ten budget 9/1000, "
    "ray epsilon-ten budget 1/(100u), "
    "1 open central-tail theorem"
)


@dataclass(frozen=True)
class Finding:
    section: str
    issue: str
    detail: str


def finding(section: str, issue: str, detail: object) -> Finding:
    return Finding(section=section, issue=issue, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate(artifact_path: Path, note_path: Path) -> list[Finding]:
    findings: list[Finding] = []
    artifact = load_json(artifact_path)
    try:
        rebuilt = build_artifact()
    except Exception as exc:
        findings.append(finding("rebuild", "failed", repr(exc)))
    else:
        if rebuilt != artifact:
            findings.append(finding("rebuild", "artifact-mismatch", "rebuilt JSON differs"))

    expected_summary = {
        "rows": 9,
        "exact_rows": 8,
        "open_analytic_rows": 1,
        "finite_exact_remainder_budget": "1/100",
        "ray_exact_remainder_budget": "1/(50u)",
        "global_epsilon_eight_layer_closed": True,
        "global_epsilon_ten_layer_closed": True,
        "finite_exact_minus_epsilon_ten_budget": "9/1000",
        "ray_exact_minus_epsilon_ten_budget": "1/(100u)",
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    finite = artifact.get("finite_budget", {})
    if finite.get("reserved_final_corridor_margin") != "7/1000":
        findings.append(finding("finite", "bad-final-margin", finite))
    ray = artifact.get("asymptotic_budget", {})
    if ray.get("reserved_final_corridor_margin") != "1/(50u)":
        findings.append(finding("ray", "bad-final-margin", ray))
    epsilon_ten = artifact.get("epsilon_ten_budget", {})
    if epsilon_ten.get("finite", {}).get("reserved_final_corridor_margin") != (
        "79999/10000000"
    ):
        findings.append(finding("epsilon-ten-finite", "bad-final-margin", epsilon_ten))
    if epsilon_ten.get("ray", {}).get("reserved_final_corridor_margin") != (
        "29/(1000u)"
    ):
        findings.append(finding("epsilon-ten-ray", "bad-final-margin", epsilon_ten))

    for ref, expected_hash in artifact.get("source_hashes", {}).items():
        path = REPO_ROOT / ref
        if not path.exists():
            findings.append(finding("source", "missing", ref))
        elif sha256(path) != expected_hash:
            findings.append(finding("source", "sha256-mismatch", ref))
    for ref in artifact.get("sources", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(finding("artifact", "missing-source", ref))
    generator = artifact.get("generator")
    if not isinstance(generator, str) or not (REPO_ROOT / generator).exists():
        findings.append(finding("artifact", "missing-generator", generator))

    text = note_path.read_text(encoding="utf-8")
    required = (
        "Status: exact formal composition and sharpened exact-density remainder target.",
        "scaled |kappa_r^[8]-kappa_r^[6]|<1/1000",
        "scaled |kappa_r-kappa_r^[8]|<1/100, 2<=u<=20",
        "scaled |kappa_r-kappa_r^[8]|<1/(50u), u>=20",
        "scaled |kappa_r^[10]-kappa_r^[8]|<1/10000000",
        "scaled |kappa_r^[10]-kappa_r^[8]|<1/(1000u)",
        "scaled |kappa_r-kappa_r^[10]|<9/1000",
        "scaled |kappa_r-kappa_r^[10]|<1/(100u)",
        "sufficient budgets, not exact-density estimates",
        "dependency widths destroy the q-scaled cancellations",
        "not a proof of the exact cumulant corridors",
    )
    for marker in required:
        if marker not in text:
            findings.append(finding("note", "missing-marker", marker))
    return findings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings = validate(args.artifact, args.note)
    if findings:
        for row in findings:
            print(f"{row.section}: {row.issue}: {row.detail}")
        print(f"order-four exact cumulant remainder budget: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
