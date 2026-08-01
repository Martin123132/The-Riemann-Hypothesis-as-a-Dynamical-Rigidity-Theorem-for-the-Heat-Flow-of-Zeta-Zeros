#!/usr/bin/env python3
"""Validate the first-order Wronskian crossing reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_wronskian_crossing_reduction as builder


EXPECTED_IDS = [
    "nfowc_01_coordinates",
    "nfowc_02_contact_wronskian",
    "nfowc_03_contact_disjunction",
    "nfowc_04_crossing_margin",
    "nfowc_05_upward_sign",
    "nfowc_06_exceptional_zero",
    "nfowc_07_count_decomposition",
    "nfowc_08_successor_link",
    "nfowc_09_phase_sign_guard",
    "nfowc_10_frequency_target",
    "nfowc_11_parabolic_target",
    "nfowc_12_shoulder_target",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    x_part, y_part, u_part, v_part = sp.symbols("X Y U V", real=True)
    wronskian = sp.expand_complex(
        sp.im(
            (u_part + sp.I * v_part)
            * sp.conjugate(x_part + sp.I * y_part)
        )
    )
    if sp.simplify(wronskian - (v_part * x_part - u_part * y_part)) != 0:
        issues.append("independent Wronskian coordinate audit failed")
    if sp.simplify(wronskian.subs(x_part, 0) + u_part * y_part) != 0:
        issues.append("independent crossing Wronskian audit failed")
    if sp.sqrt(builder.HALF_SQUARED_BUDGET) != 50_000 * sp.sqrt(5):
        issues.append("independent crossing threshold audit failed")


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != builder.STEM:
        issues.append("artifact kind mismatch")
    if artifact.get("source_sha256") != builder.source_hashes():
        issues.append("source hash chain drifted")
    if artifact.get("source_audit") != builder.source_audit():
        issues.append("source audit drifted")
    if artifact.get("exact") != builder.build_exact():
        issues.append("exact payload drifted")
    if artifact.get("constants") != {
        "eta_0": 100_000,
        "eta_1": 200_000,
        "half_squared_budget": 12_500_000_000,
        "crossing_threshold": "50000*sqrt(5)",
    }:
        issues.append("constant payload drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    open_rows = [
        row for row in rows if row.get("role") == "open_theorem_target"
    ]
    if len(open_rows) != 3:
        issues.append(f"open target count drifted: {len(open_rows)}")
    if any(row.get("readiness") != "not_ready_to_apply" for row in open_rows):
        issues.append("an open target was promoted")

    independent_symbolic_audit(issues)
    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "W_[1]=V*X-U*Y",
        "2*|W_[1]|<exp(-5L/4)",
        "50000*sqrt(5)",
        "W_[1]*Y<0",
        "N_(E_[1]=0,U>0)",
        "kappa_j=N_up(t_j;R_j)",
        "oriented count is <1",
        "inertia (1,1,m-2)",
        "q=2tL^2>=1",
        "q<1",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    proof_boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove",
        "q>=1",
        "q<1",
        "finite shoulders",
        "one-sided successor winding",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in proof_boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "First-Order Wronskian",
        "Boundary Crossings",
        "Successor Composition",
        "Route Guards",
        "Remaining Xi Input",
        "Live Handoff",
        "explicit `E_[1]=0` class",
        "not a proof of RH",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=builder.DEFAULT_OUT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Newman first-order Wronskian crossing reduction: "
        "12 rows, eta_0=100000, eta_1=200000, "
        "2 exact crossing classes, 3 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
