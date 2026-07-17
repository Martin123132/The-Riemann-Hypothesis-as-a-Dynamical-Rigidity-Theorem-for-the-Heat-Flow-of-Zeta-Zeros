#!/usr/bin/env python3
"""Validate the next-parity asymptotic-ray coefficient certificate."""

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

from jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
    sha256,
)


EXPECTED_SUMMARY = (
    "validated order-four next-parity asymptotic ray: 8 rows, 0 issues, "
    "7 exact rows, 7 global coefficient bounds, 14 leading buffer gates, "
    "4 new jet-remainder sign gates, 1 open exact-density remainder"
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
        "coefficient_bounds": 7,
        "leading_buffer_gates": 14,
        "odd_nonvanishing_gates": 3,
        "normalized_jet_cap_gates": 8,
        "reused_jet_remainder_sign_gates": 14,
        "new_jet_remainder_sign_gates": 4,
        "global_coefficient_layer_closed": True,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    if artifact.get("scalar_geometry", {}).get("coefficient_transfer") != (
        "|C_r-C_r^(infinity)|<=64000000*eta(u)<1/(200u)"
    ):
        findings.append(finding("geometry", "bad-transfer", artifact.get("scalar_geometry")))
    if artifact.get("jet_remainders", {}).get("new_sign_gates") != 4:
        findings.append(finding("jet-remainders", "bad-new-gate-count", artifact.get("jet_remainders")))

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
        "Status: exact analytic theorem for the global next-parity coefficient layer.",
        "coefficient floor+1/(100u)",
        "0<L_r^(infinity)<8/5",
        "|V^(r)-q*P_r(u)|<=10000*u^10",
        "|C_r-C_r^(infinity)|<=64000000*eta(u)<1/(200u)",
        "hold for every `u>=2`",
        "No exact cumulant",
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
        print(f"order-four next-parity asymptotic ray: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
