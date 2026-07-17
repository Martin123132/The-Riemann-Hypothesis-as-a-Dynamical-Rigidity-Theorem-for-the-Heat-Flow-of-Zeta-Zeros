#!/usr/bin/env python3
"""Validate the exact-corridor localized-curvature ray certificate."""

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

from jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
    sha256,
)


EXPECTED_SUMMARY = (
    "validated order-four exact-corridor curvature ray: 8 rows, 0 issues, "
    "8 exact rows, 7 normalized H boxes, 6 defect bounds, "
    "5 localized gates, global corridor-to-curvature closed"
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
        "exact_rows": 8,
        "open_analytic_rows": 0,
        "normalized_h_boxes": 7,
        "defect_derivative_bounds": 6,
        "localized_scaled_gates": 5,
        "positive_rational_comparisons": 1,
        "ray_corridor_to_curvature_closed": True,
        "global_corridor_to_curvature_closed": True,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))

    scalar = artifact.get("scalar_comparison", {})
    if scalar.get("scaled_U_upper") != "3011223637/866377000":
        findings.append(finding("scalar", "unexpected-upper", scalar))
    if scalar.get("margin_lower") != "21095863/866377000":
        findings.append(finding("scalar", "unexpected-margin", scalar))
    gates = artifact.get("localized_scaled_gates", {}).get("gates", {})
    expected_gates = {
        "t2_ell_second_upper": "23/20",
        "A_minus_lower": "193/100",
        "A_plus_upper": "201/100",
        "C_plus_upper": "401/100",
        "P_minus_lower": "191/100",
        "j_plus_cap": "1/1000",
    }
    for key, value in expected_gates.items():
        if gates.get(key) != value:
            findings.append(finding("gates", f"unexpected-{key}", gates.get(key)))

    boxes = artifact.get("normalized_h_boxes", {}).get("boxes", {})
    if len(boxes) != 7:
        findings.append(finding("boxes", "wrong-count", len(boxes)))
    if boxes.get("2", {}).get("required_floor") != "97/100":
        findings.append(finding("boxes", "x2-floor", boxes.get("2")))
    if boxes.get("3", {}).get("required_floor") != "24/25":
        findings.append(finding("boxes", "x3-floor", boxes.get("3")))

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
        "Status: rigorous exact-corridor localized-curvature theorem",
        "`u>=20`",
        "`u>=2`",
        "q >= 1000000000000000000000000000000000",
        "x_2 >= 97/100",
        "x_3 >= 24/25",
        "t^2*ell'' <= 23/20",
        "negative square term is retained",
        "3011223637/866377000",
        "21095863/866377000",
        "now complete for every `u>=2`",
        "not a proof of order-four entry",
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
        print(f"order-four exact-corridor curvature ray: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
