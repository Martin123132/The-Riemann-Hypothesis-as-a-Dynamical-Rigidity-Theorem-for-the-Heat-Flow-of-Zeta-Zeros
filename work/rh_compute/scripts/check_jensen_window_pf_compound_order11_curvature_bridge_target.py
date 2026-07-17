#!/usr/bin/env python3
"""Validate the conditional order-eleven curvature bridge target."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_curvature_bridge_target as target  # noqa: E402


@dataclass(frozen=True)
class Finding:
    section: str
    issue: str
    detail: str


def finding(section: str, issue: str, detail: object) -> Finding:
    return Finding(section, issue, str(detail))


def validate(artifact_path: Path, note_path: Path) -> list[Finding]:
    findings: list[Finding] = []
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    try:
        rebuilt = target.build_artifact()
    except Exception as exc:
        findings.append(finding("rebuild", "failed", repr(exc)))
    else:
        if rebuilt != artifact:
            findings.append(finding("rebuild", "artifact-mismatch", "rebuilt JSON differs"))
    if artifact.get("kind") != "jensen_window_pf_compound_order11_curvature_bridge_target":
        findings.append(finding("artifact", "bad-kind", artifact.get("kind")))
    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 9,
        "ready_rows": 8,
        "open_rows": 1,
        "exact_factorizations": 1,
        "power_envelope_rows": 18,
        "conditional_transfer_theorems": 1,
        "conditional_endpoint_tail_theorems": 1,
        "open_continuum_targets": 1,
        "finite_prefix_theorems": 1,
        "conditional_heat_handoffs": 1,
    }
    if summary != expected_summary:
        findings.append(finding("summary", "mismatch", summary))
    exact = artifact.get("exact", {})
    required = {
        "canonical_factorization": "Q_(10,n)=A_(n+9)^10*exp(y(n+9))",
        "first_X_floor": "X_j^(1)>=9/(2*j+1)-4201/j^2>1/j, j>=1252",
        "full_X_floor": "X_j>=2259/(250*(2*j+1))-4211/j^2>1/j, j>=1252",
        "full_transfer": target.Y_FULL_TRANSFER,
        "continuous_target": target.Y_CONTINUOUS_TARGET,
        "conditional_full_ceiling": target.Y_FULL_CEILING,
        "rational_comparison": "6100/k^2<2510/(250*(2*k+1)), k>=1253",
        "finite_prefix_theorem": target.ORDER11_FINITE_PREFIX,
        "delayed_heat_theorem": target.ORDER11_DELAYED_HANDOFF,
    }
    for key, expected in required.items():
        if exact.get(key) != expected:
            findings.append(finding("exact", f"bad-{key}", exact.get(key)))
    envelope = exact.get("power_envelope", {})
    rows = envelope.get("rows", [])
    if len(rows) != 18:
        findings.append(finding("envelope", "bad-row-count", len(rows)))
    if [row.get("power") for row in rows] != [13, 12, 12, 11, 11, 10, 10, 9, 9, 8, 8, 7, 7, 6, 6, 5, 5, 4]:
        findings.append(finding("envelope", "bad-power-ladder", [row.get("power") for row in rows]))
    if [row.get("start") for row in rows] != [341, 341, 342, 342, 343, 343, 344, 344, 345, 345, 346, 1249, 1250, 1250, 1251, 1251, 1252, 1252]:
        findings.append(finding("envelope", "bad-start-ladder", [row.get("start") for row in rows]))
    scaled = Fraction(envelope.get("transfer_scaled_exact", "999"))
    if not scaled < 37:
        findings.append(finding("envelope", "transfer-exhausted", scaled))
    if not Fraction(36) < scaled:
        findings.append(finding("envelope", "unexpected-transfer-scale", scaled))
    for key in ("first_floor_polynomial", "full_floor_polynomial"):
        polynomial = exact.get(key, {})
        coefficients = [int(value) for value in polynomial.get("shifted_coefficients", [])]
        if polynomial.get("start") != 1252 or not coefficients or min(coefficients) <= 0:
            findings.append(finding("floor", f"bad-{key}", polynomial))
    if [int(value) for value in exact.get("shifted_coefficients", [])] != [2510, 3240060, 117547590]:
        findings.append(finding("comparison", "bad-shifted-coefficients", exact.get("shifted_coefficients")))
    rows_artifact = artifact.get("rows", [])
    if len(rows_artifact) != 9 or sum(row.get("readiness") == "not_ready_to_apply" for row in rows_artifact) != 1:
        findings.append(finding("rows", "readiness-mismatch", rows_artifact))
    note = note_path.read_text(encoding="utf-8")
    for marker in (
        "exact scaled transfer=",
        "<37",
        "eighteen exact rational rows",
        "y_1''(t)<=6000/t^2",
        target.ORDER11_FINITE_PREFIX,
        target.ORDER11_DELAYED_HANDOFF,
        "Only the continuum target",
        "not all-shift order eleven",
        "PF-infinity, RH, or `Lambda<=0`",
    ):
        if marker not in note:
            findings.append(finding("note", "missing-marker", marker))
    return findings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=target.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=target.DEFAULT_NOTE)
    args = parser.parse_args()
    findings = validate(args.artifact, args.note)
    if findings:
        for row in findings:
            print(f"{row.section}: {row.issue}: {row.detail}")
        print(f"order-eleven curvature bridge target: {len(findings)} issues")
        return 1
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    scaled = artifact["exact"]["power_envelope"]["transfer_scaled_decimal"]
    print(f"validated order-eleven curvature bridge target: 18 envelope rows, scaled transfer {scaled}<37, 0 issues")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
