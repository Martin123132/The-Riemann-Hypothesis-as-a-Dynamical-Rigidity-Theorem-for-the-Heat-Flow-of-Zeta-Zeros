#!/usr/bin/env python3
"""Validate the residual-placement and cumulative-current handoff gate."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_residual_placement_c1_cumulative_handoff_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "residual_placement_c1_cumulative_handoff_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    "rpc_01_whole_jet",
    "rpc_02_c1_box",
    "rpc_03_main_components",
    "rpc_04_allocation_gauge",
    "rpc_05_four_defects",
    "rpc_06_allocation_witness",
    "rpc_07_derivative_order",
    "rpc_08_c2_nonrequirement",
    "rpc_09_chart_covariance",
    "rpc_10_cutoff_roles",
    "rpc_11_retained_edge",
    "rpc_12_polarization",
    "rpc_13_diagonal_reserve",
    "rpc_14_sum_guard",
    "rpc_15_invariant_handoff",
    "rpc_16_boundary",
]
EXPECTED_COUNTS = {
    "rows": 16,
    "global_remainder_coordinates": 2,
    "allocation_gauge_parameters": 1,
    "allocation_dependent_edge_current_guards": 1,
    "required_c2_for_boundary_transfer": 0,
    "required_c2_for_allocated_component_current": 1,
    "whole_jet_c1_homotopies": 1,
    "retained_model_signed_edge_families": 1,
    "retained_model_signed_interior_families": 1,
    "exact_cumulative_polarizations": 1,
    "open_cross_current_or_abel_targets": 1,
    "xi_level_edge_signs": 0,
    "contact_exclusions": 0,
}
SUMMARY = (
    "validated residual-placement C1 cumulative handoff gate: "
    "16 rows, 2 C1 residual coordinates, 1 allocation-sign guard, "
    "0 C2 boundary requirements, 1 whole-jet homotopy, "
    "1 cumulative polarization, 1 open cross-current/Abel target"
)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing result: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid result: {exc}")
        return {}


def independent_allocation_audit(issues: list[str]) -> None:
    c, d, u, w = sp.symbols("c d u w", real=True)
    c_x, d_x, u_x, w_x = sp.symbols(
        "c_x d_x u_x w_x", real=True
    )
    lam = sp.symbols("lambda", real=True)
    base = c * d_x - d * c_x
    allocated = (c + lam * u) * (d_x + lam * w_x) - (
        d + lam * w
    ) * (c_x + lam * u_x)
    expected = (
        base
        + lam * (c * w_x + u * d_x - d * u_x - w * c_x)
        + lam**2 * (u * w_x - w * u_x)
    )
    if sp.simplify(allocated - expected) != 0:
        issues.append("independent allocation identity failed")

    witness = allocated.subs(
        {
            c: 1,
            d: 0,
            c_x: 0,
            d_x: -1,
            u: 0,
            w: 1,
            u_x: -2,
            w_x: 0,
        }
    )
    if sp.simplify(witness - (-1 + 2 * lam**2)) != 0:
        issues.append("independent allocation sign witness failed")


def independent_polarization_audit(issues: list[str]) -> None:
    c = sp.symbols("c0:4", real=True)
    d = sp.symbols("d0:4", real=True)
    c_x = sp.symbols("cx0:4", real=True)
    d_x = sp.symbols("dx0:4", real=True)
    total = sp.expand(sum(c) * sum(d_x) - sum(d) * sum(c_x))
    diagonal = sum(c[i] * d_x[i] - d[i] * c_x[i] for i in range(4))
    cross = sum(
        c[i] * d_x[j]
        + c[j] * d_x[i]
        - d[i] * c_x[j]
        - d[j] * c_x[i]
        for i in range(4)
        for j in range(i + 1, 4)
    )
    if sp.simplify(total - diagonal - cross) != 0:
        issues.append("independent four-atom polarization failed")

    theta = sp.symbols("theta", real=True)
    v1 = sp.Matrix([sp.cos(theta), -sp.sin(theta)])
    v2 = sp.Matrix([sp.cos(3 * theta), -sp.sin(3 * theta)])
    for vector, target in ((v1, -1), (v2, -3)):
        derivative = vector.diff(theta)
        current = vector[0] * derivative[1] - vector[1] * derivative[0]
        if sp.simplify(current - target) != 0:
            issues.append("independent clockwise component failed")
    if sp.simplify((v1 + v2).subs(theta, sp.pi / 2)) != sp.zeros(2, 1):
        issues.append("independent clockwise sum guard failed")


def validate(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    if payload.get("date") != "2026-07-31":
        issues.append("artifact date drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    if [row.get("id") for row in payload.get("rows", [])] != EXPECTED_IDS:
        issues.append("row ids or order drifted")

    for token in (
        "residual-placement",
        "C1 whole-jet",
        "cumulative-current",
    ):
        if token not in payload.get("status", ""):
            issues.append(f"status missing token: {token}")
    for token in (
        "noninvariant",
        "no cross-current bound",
        "Abel gap",
        "contact exclusion",
        "Lambda<=0",
        "RH",
    ):
        if token not in payload.get("proof_boundary", ""):
            issues.append(f"proof boundary missing token: {token}")

    hashes = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    try:
        if payload.get("symbolic_certificate") != gate.symbolic_certificate():
            issues.append("stored symbolic certificate drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"symbolic certificate recomputation failed: {exc}")

    exact = payload.get("exact", {})
    for token in (
        "V_Z=(Z,partial_x Z/L)=V_J+V_r",
        "C_edge,D_edge",
        "Delta V_r=-Delta V_J",
        "C_cross<R_diag",
        "no derivative of V_r is used",
    ):
        if token not in json.dumps(exact):
            issues.append(f"exact payload missing token: {token}")

    independent_allocation_audit(issues)
    independent_polarization_audit(issues)


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Allocation Gauge",
        "Derivative Order",
        "Cutoff Covariance",
        "Cumulative Current",
        "The four-defect algebra is correct",
        "not invariant",
        "No derivative of `V_r`",
        "C_cross<R_diag",
        "No cumulative cross-current bound",
    ):
        if token not in text:
            issues.append(f"note marker missing: {token}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
        validate(payload, issues)
    validate_note(issues)
    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        return 1
    print(SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
