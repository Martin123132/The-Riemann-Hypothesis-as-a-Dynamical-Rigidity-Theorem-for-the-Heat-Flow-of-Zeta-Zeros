#!/usr/bin/env python3
"""Validate the positive-boundary diagonal compact-exhaustion gate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate as gate


EXPECTED_IDS = [
    "npbdeg_01_boundary_attainment",
    "npbdeg_02_abstract_exhaustion",
    "npbdeg_03_linear_exhaustion",
    "npbdeg_04_compact_shell",
    "npbdeg_05_compact_minimum",
    "npbdeg_06_independent_rates",
    "npbdeg_07_arbitrary_center_field",
    "npbdeg_08_shifted_boundary_flow",
    "npbdeg_09_local_field_guard",
    "npbdeg_10_proof_handoff",
]


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != gate.STEM:
        issues.append("artifact kind mismatch")
    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if artifact.get("exact") != gate.build_exact():
        issues.append("exact payload drifted")
    if artifact.get("source_audit") != gate.source_audit():
        issues.append("source audit drifted")
    if artifact.get("linear_examples") != gate.linear_examples():
        issues.append("linear examples drifted")
    if artifact.get("field_examples") != gate.field_examples():
        issues.append("field examples drifted")

    expected_roles = {
        "published_input": 1,
        "exact_equivalence": 2,
        "exact_composition": 1,
        "exact_reduction": 1,
        "nonpromotion_gate": 2,
        "exact_countermodel": 2,
        "open_theorem_target": 1,
    }
    role_counts = {
        role: sum(row.get("role") == role for row in rows)
        for role in expected_roles
    }
    if role_counts != expected_roles:
        issues.append(f"role counts drifted: {role_counts}")

    guard_rows = [
        row
        for row in rows
        if row.get("role") in {"nonpromotion_gate", "exact_countermodel"}
    ]
    if len(guard_rows) != 4 or any(
        row.get("readiness") != "guard_validated" for row in guard_rows
    ):
        issues.append("diagonal/local-field guards are not active")
    open_rows = [
        row for row in rows if row.get("role") == "open_theorem_target"
    ]
    if len(open_rows) != 1 or open_rows[0].get("readiness") != (
        "not_ready_to_apply"
    ):
        issues.append("independent-rate target was promoted")

    exact_text = json.dumps(artifact.get("exact", {}))
    for marker in (
        "delta_j->0",
        "R_j->infinity",
        "R_j=38+j",
        "38<|x|<=38+j",
        "R_j=38+log(1+j)",
        "not required",
        "a_c^2=c^2*(3-b*c)/(1-b*c)",
        "B_c=1/c+2*c/(c^2-a_c^2)=b",
        "K_c=1/(2*c^2)",
        "pi^2/64",
        "Newman-style boundary lambda",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    field_rows = artifact.get("field_examples", [])
    if len(field_rows) != 4:
        issues.append("field example count mismatch")
    else:
        expected_field = -math.pi / 8
        for row in field_rows:
            if abs(float(row.get("field", 0.0)) - expected_field) > 1e-15:
                issues.append(f"field specialization drifted at c={row.get('c')}")
        if not (
            float(field_rows[-1]["stiffness"])
            < float(field_rows[0]["stiffness"])
            and abs(
                float(field_rows[-1]["stiffness"]) - math.pi**2 / 64
            )
            < 0.002
        ):
            issues.append("arbitrary-height stiffness limit diagnostic failed")

    status = artifact.get("status", "")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "independent diagonal",
        "remains open",
        "not a proof",
    ):
        if marker not in status:
            issues.append(f"status marker missing: {marker}")
    for marker in (
        "unnecessary exponential coupling",
        "does not certify",
        "positive-time simplicity",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = gate.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Diagonal Exhaustion",
        "Quantifier Correction",
        "Arbitrary-Height Field Guard",
        "Live Handoff",
        "38+log(1+j)",
        "not a proof of",
        "`Lambda <= 0` or RH",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, default=gate.DEFAULT_OUT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman positive-boundary diagonal-exhaustion gate: "
        f"{len(artifact['rows'])} rows, 0 issues, "
        "5 exact exhaustion reductions, 1 compact-shell composition, "
        "1 arbitrary-height classical-field boundary countermodel, "
        "1 independent-rate open target"
    )


if __name__ == "__main__":
    main()
