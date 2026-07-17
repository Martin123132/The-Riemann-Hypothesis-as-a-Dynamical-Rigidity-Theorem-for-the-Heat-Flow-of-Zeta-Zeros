#!/usr/bin/env python3
"""Validate the formal-cumulant asymptotic-ray certificate."""

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

from jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
)


EXPECTED_SUMMARY = (
    "validated order-four formal cumulant asymptotic ray: 8 rows, 0 issues, "
    "7 exact rows, 7 buffered corridors, 14 jet-remainder sign gates, "
    "1 open exact-density remainder"
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
        "rows": 8,
        "exact_rows": 7,
        "open_analytic_rows": 1,
        "leading_corridor_gates": 7,
        "jet_remainder_sign_gates": 14,
        "formal_ray_closed": True,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    row_ids = {row.get("id") for row in artifact.get("rows", [])}
    if len(row_ids) != 8 or "co4fcarc_07_formal_ray_theorem" not in row_ids:
        findings.append(finding("rows", "bad-ids", sorted(str(value) for value in row_ids)))
    if artifact.get("scalar_geometry", {}).get("formal_transfer") != (
        "|R_r^[6]-F_r|<=22000000*u^6/q<1/(20u), u>=20"
    ):
        findings.append(finding("geometry", "bad-transfer", artifact.get("scalar_geometry")))
    for ref in artifact.get("sources", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(finding("artifact", "missing-source", ref))
    generator = artifact.get("generator")
    if not isinstance(generator, str) or not (REPO_ROOT / generator).exists():
        findings.append(finding("artifact", "missing-generator", generator))

    text = note_path.read_text(encoding="utf-8")
    required = (
        "Status: exact analytic theorem for the epsilon-six formal cumulant ray.",
        "corridor floor+1/(10u) <= F_r <= corridor ceiling-1/(10u)",
        "|V^(r)-q*P_r(u)|<=1000*u^8",
        "|R_r^[6]-F_r|<=22000000*u^6/q<1/(20u)",
        "formal cumulant polynomial satisfies every candidate",
        "model is now closed for every `u>=2`",
        "not a proof of the exact cumulant ray",
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
        print(f"order-four formal cumulant asymptotic ray: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
