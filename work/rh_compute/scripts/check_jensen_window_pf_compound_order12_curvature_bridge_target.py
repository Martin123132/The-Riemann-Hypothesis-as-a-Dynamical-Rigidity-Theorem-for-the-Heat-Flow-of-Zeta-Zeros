#!/usr/bin/env python3
"""Validate the conditional order-twelve curvature bridge target."""

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

import jensen_window_pf_compound_order12_curvature_bridge_target as target  # noqa: E402


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
    if artifact.get("kind") != "jensen_window_pf_compound_order12_curvature_bridge_target":
        findings.append(finding("artifact", "bad-kind", artifact.get("kind")))
    expected_summary = {
        "rows": 11,
        "ready_rows": 10,
        "open_rows": 1,
        "exact_factorizations": 1,
        "power_envelope_rows": 20,
        "conditional_transfer_theorems": 1,
        "conditional_endpoint_tail_theorems": 1,
        "open_continuum_targets": 1,
        "open_finite_endpoint_targets": 0,
        "lambda0_prefix_theorems": 1,
        "open_lambda0_prefix_targets": 0,
        "conditional_heat_handoffs": 1,
        "orders_above_12": 0,
        "rh_claims": 0,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    exact = artifact.get("exact", {})
    required = {
        "canonical_factorization": "Q_(11,n)=A_(n+10)^11*exp(v(n+10))",
        "first_U_floor": "U_j^(1)>=10/(2*j+1)-6001/j^2>1/j, j>=1503",
        "full_U_floor": "U_j>=2510/(250*(2*j+1))-6038/j^2>1/j, j>=1503",
        "continuous_target": target.V_CONTINUOUS_TARGET,
        "full_transfer": target.V_FULL_TRANSFER,
        "conditional_full_ceiling": target.V_FULL_CEILING,
        "certified_endpoint_sign_chart": "Q_(12,n)(-100)<0 for n=0,1,2,3 and Q_(12,n)(-100)>0 for every 4<=n<=1240",
        "finite_prefix_theorem": "Q_(12,n)(-100)>0 for every 1241<=n<=1492",
        "lambda0_prefix_target": target.ORDER12_LAMBDA0_PREFIX_TARGET,
        "delayed_heat_theorem": target.ORDER12_DELAYED_HANDOFF,
    }
    for key, expected in required.items():
        if exact.get(key) != expected:
            findings.append(finding("exact", f"bad-{key}", exact.get(key)))
    envelope = exact.get("power_envelope", {})
    rows = envelope.get("rows", [])
    if len(rows) != 20:
        findings.append(finding("envelope", "bad-row-count", len(rows)))
    scaled = Fraction(envelope.get("transfer_scaled_exact", "999"))
    if not Fraction(99) < scaled < 100:
        findings.append(finding("envelope", "unexpected-transfer-scale", scaled))
    for key in ("first_floor_polynomial", "full_floor_polynomial"):
        polynomial = exact.get(key, {})
        coefficients = [int(value) for value in polynomial.get("shifted_coefficients", [])]
        if polynomial.get("start") != 1503 or not coefficients or min(coefficients) <= 0:
            findings.append(finding("floor", f"bad-{key}", polynomial))
    coefficients = [int(value) for value in exact.get("shifted_coefficients", [])]
    if not coefficients or min(coefficients) <= 0:
        findings.append(finding("comparison", "bad-shifted-coefficients", coefficients))
    rows_artifact = artifact.get("rows", [])
    if len(rows_artifact) != 11:
        findings.append(finding("rows", "bad-row-count", len(rows_artifact)))
    if sum(row.get("readiness") == "not_ready_to_apply" for row in rows_artifact) != 1:
        findings.append(finding("rows", "readiness-mismatch", rows_artifact))
    note = note_path.read_text(encoding="utf-8")
    for marker in (
        "exact conditional reduction with one open continuum target",
        "This is not a proof of order twelve",
        "U(t)=10*B(t)",
        "twenty exact rational rows",
        "exact scaled transfer=",
        "<100",
        target.V_CONTINUOUS_TARGET,
        target.ORDER12_FINITE_PREFIX_THEOREM,
        target.ORDER12_LAMBDA0_PREFIX_TARGET,
        "No order above twelve",
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
        print(f"order-twelve curvature bridge target: {len(findings)} issues")
        return 1
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    scaled = artifact["exact"]["power_envelope"]["transfer_scaled_decimal"]
    print(
        "validated order-twelve curvature bridge target: "
        f"20 envelope rows, scaled transfer {scaled}<100, 1 open target, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
