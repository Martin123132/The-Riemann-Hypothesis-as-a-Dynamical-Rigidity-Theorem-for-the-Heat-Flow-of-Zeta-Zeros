#!/usr/bin/env python3
"""Validate the Newman first-jet winding and boundary-flux gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import jensen_window_pf_newman_first_jet_winding_gate as gate
import sympy as sp


EXPECTED_IDS = [
    "njwg_01_universal_heat_polynomial",
    "njwg_02_discriminant",
    "njwg_03_root_count_jump",
    "njwg_04_local_index",
    "njwg_05_global_winding",
    "njwg_06_phase_flux",
    "njwg_07_xi_half_rectangles",
    "njwg_08_known_edge_composition",
    "njwg_09_winding_guard",
    "njwg_10_proof_handoff",
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
    if artifact.get("hermite_audits") != gate.hermite_audits():
        issues.append("Hermite audits drifted")

    f_x, f_xx, f_xt = sp.symbols("f_x f_xx f_xt", real=True)
    f_t = -f_xx
    independent_x_t = sp.Matrix(
        [[f_x, f_t], [f_xx, f_xt]]
    ).det().subs(f_x, 0)
    independent_t_x = sp.Matrix(
        [[f_t, f_x], [f_xt, f_xx]]
    ).det().subs(f_x, 0)
    orientation = artifact.get("exact", {}).get("orientation_audit", {})
    if independent_x_t != f_xx**2:
        issues.append("independent standard-orientation determinant failed")
    if independent_t_x != -f_xx**2:
        issues.append("independent reversed-orientation determinant failed")
    if orientation.get("determinant_x_t") != str(independent_x_t):
        issues.append("stored standard-orientation determinant drifted")
    if orientation.get("determinant_t_x") != str(independent_t_x):
        issues.append("stored reversed-orientation determinant drifted")
    if orientation.get("positive_bottom_edge") != "t=delta, x:0->R":
        issues.append("positive boundary convention drifted")

    expected_roles = {
        "exact_local_model": 1,
        "exact_all_multiplicity_identity": 1,
        "exact_all_multiplicity_reduction": 1,
        "exact_topological_lemma": 2,
        "exact_differential_identity": 1,
        "exact_equivalence": 1,
        "exact_composition": 1,
        "nonpromotion_gate": 1,
        "open_theorem_target": 1,
    }
    role_counts = {
        role: sum(row.get("role") == role for row in rows)
        for role in expected_roles
    }
    if role_counts != expected_roles:
        issues.append(f"role counts drifted: {role_counts}")

    open_rows = [
        row for row in rows if row.get("role") == "open_theorem_target"
    ]
    if len(open_rows) != 1 or open_rows[0].get("readiness") != (
        "not_ready_to_apply"
    ):
        issues.append("Xi edge target was promoted")
    guards = [
        row for row in rows if row.get("role") == "nonpromotion_gate"
    ]
    if len(guards) != 1 or guards[0].get("readiness") != "guard_validated":
        issues.append("winding nonpromotion guard is not active")

    audits = artifact.get("hermite_audits", [])
    if len(audits) != 9:
        issues.append("Hermite audit count mismatch")
    for row in audits:
        m = int(row.get("multiplicity", 0))
        expected_pairs = m // 2
        if row.get("pair_births") != expected_pairs:
            issues.append(f"pair-birth count drifted at m={m}")
        if row.get("local_index") != expected_pairs:
            issues.append(f"local index drifted at m={m}")
        if row.get("discriminant_tau_power") != m * (m - 1) // 2:
            issues.append(f"discriminant exponent drifted at m={m}")

    exact_text = json.dumps(artifact.get("exact", {}))
    for marker in (
        "Disc_u P_m",
        "+floor(m/2)",
        "det D_(x,t)(F,F_x)",
        "F_xx^2>0",
        "reversed (t,x)",
        "wind((F+iF_x)(partial Omega),0)",
        "partial_x arg(H+iH_x)=-L/Q",
        "partial_t arg(H+iH_x)=L_x/Q",
        "Q_j=[1/(5j),1/4]x[0,38+j]",
        "bottom edge",
        "right edge",
        "2pi*wind",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    status = artifact.get("status", "")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "signed winding theorem",
        "remain open",
        "not a proof",
    ):
        if marker not in status:
            issues.append(f"status marker missing: {marker}")
    for marker in (
        "does not prove",
        "bottom/right Xi edge",
        "zero winding",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = gate.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Universal Contact Charge",
        "Laguerre Phase Flux",
        "Xi Half-Rectangle Criterion",
        "Scope Guard",
        "Live Handoff",
        "Hidden contacts cannot",
        "not a proof of `Lambda <= 0` or RH",
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
        "validated Newman first-jet winding gate: "
        f"{len(artifact['rows'])} rows, 0 issues, "
        f"{len(artifact['hermite_audits'])} Hermite audits, "
        "5 exact local/global index reductions, 1 signed winding theorem, "
        "1 half-rectangle flux criterion, 1 open Xi edge target"
    )


if __name__ == "__main__":
    main()
