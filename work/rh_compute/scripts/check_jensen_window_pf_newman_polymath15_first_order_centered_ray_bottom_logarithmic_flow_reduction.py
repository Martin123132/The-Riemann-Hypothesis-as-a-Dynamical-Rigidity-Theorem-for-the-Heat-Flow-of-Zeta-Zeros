#!/usr/bin/env python3
"""Validate the centered ray-bottom logarithmic-flow reduction."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_ray_bottom_logarithmic_flow_reduction"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "rectangular_boundary": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_rectangular_boundary_degree_reduction.json"
    ),
    "normalized_phase_flux": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_normalized_prefix_phase_flux_reduction.json"
    ),
    "abel_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_carrier_kernel_abel_prefix_reduction.json"
    ),
    "abel_shear": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_abel_scalar_shear_flux_reduction.json"
    ),
    "target_reconciliation": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_target_reconciliation_gate.json"
    ),
    "connector_cap": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_dominant_ray_connector_phase_cap.json"
    ),
}
EXPECTED_IDS = [
    "rlf_01_ray_schedule",
    "rlf_02_hard_interval",
    "rlf_03_cutoff_cells",
    "rlf_04_equality_ownership",
    "rlf_05_log_chain_rule",
    "rlf_06_terminal_distance",
    "rlf_07_carrier_rates",
    "rlf_08_inward_current",
    "rlf_09_ordered_rotation",
    "rlf_10_prefix_flow",
    "rlf_11_endpoint_value_flow",
    "rlf_12_endpoint_slope_flow",
    "rlf_13_projector_flow",
    "rlf_14_abel_scalar_flow",
    "rlf_15_ray_proxy",
    "rlf_16_exact_phase_flux",
    "rlf_17_orientation_defect",
    "rlf_18_scaled_orientation",
    "rlf_19_ray_intersections",
    "rlf_20_cutoff_join",
    "rlf_21_local_flux_target",
    "rlf_22_recrossing_guard",
    "rlf_23_joined_ledger",
    "rlf_24_handoff",
]


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, label: str, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing {label}: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid {label}: {exc}")
        return {}


def validate_ray_geometry(issues: list[str]) -> None:
    ell, t, j = sp.symbols("L t j", positive=True)
    n = sp.symbols("N", positive=True, integer=True)
    pi = sp.pi
    x = 4 * pi * sp.exp(ell)
    a_sq = sp.exp(ell) + t / 16
    t_0 = x / 2 + pi * t / 8
    u_n = sp.log(a_sq) / 2 - sp.log(n)
    expected_u = sp.exp(ell) / (2 * a_sq)
    if sp.simplify(sp.diff(u_n, ell) - expected_u) != 0:
        issues.append("terminal-distance L derivative failed")
    if sp.simplify(expected_u - x / (4 * t_0)) != 0:
        issues.append("terminal-distance T_0 identity failed")

    cutoff = sp.log(n**2 - t / 16)
    if sp.simplify(a_sq.subs(ell, cutoff) - n**2) != 0:
        issues.append("cutoff equality ownership failed")

    t_j = sp.Rational(25, 1) / (100 + j)
    q = 2 * t_j * ell**2
    q_one = sp.sqrt((100 + j) / 50)
    c_edge = sp.symbols("c_edge", positive=True)
    low_c = c_edge * (100 + j) / 25
    if sp.simplify(q.subs(ell, q_one) - 1) != 0:
        issues.append("ray q=1 endpoint failed")
    if sp.simplify((t_j * ell).subs(ell, low_c) - c_edge) != 0:
        issues.append("ray low-c endpoint failed")


def validate_flux_algebra(issues: list[str]) -> None:
    ell = sp.symbols("L", positive=True, real=True)
    x_value = sp.symbols("x", positive=True, real=True)
    x_fun = sp.Function("X")(ell)
    c_fun = sp.Function("C")(ell)
    denominator = ell**2 * x_fun**2 + c_fun**2
    numerator = (
        ell * x_fun * sp.diff(c_fun, ell)
        - x_fun * c_fun
        - ell * c_fun * sp.diff(x_fun, ell)
    )
    real = x_fun
    imag = c_fun / ell
    direct = sp.simplify(
        (real * sp.diff(imag, ell) - imag * sp.diff(real, ell))
        / (real**2 + imag**2)
    )
    if sp.simplify(direct - numerator / denominator) != 0:
        issues.append("ray phase-flux identity failed")

    defect = sp.symbols("E", real=True)
    decomposed = (
        -ell * x_value * c_fun**2
        + x_fun * (ell * sp.diff(c_fun, ell) - c_fun)
        - ell * c_fun * defect
    )
    substituted = numerator.subs(
        sp.diff(x_fun, ell), x_value * c_fun + defect
    )
    if sp.simplify(substituted - decomposed) != 0:
        issues.append("orientation-defect decomposition failed")

    crossing_c, crossing_x_l = sp.symbols(
        "C_cross X_L_cross", real=True, nonzero=True
    )
    crossing = sp.simplify(
        (-ell * crossing_c * crossing_x_l) / crossing_c**2
    )
    expected_crossing = -ell * crossing_x_l / crossing_c
    if sp.simplify(crossing - expected_crossing) != 0:
        issues.append("crossing-speed identity failed")


def validate_projector_and_shear(issues: list[str]) -> None:
    x_0, y_0, p_0, q_0, omega = sp.symbols(
        "X Y P_0 Q_0 Omega", real=True
    )
    rotated_value = x_0 + sp.I * y_0
    rotated_fixed_frame_rate = p_0 + sp.I * q_0
    rotated_total_rate = (
        rotated_fixed_frame_rate + sp.I * omega * rotated_value
    )
    if sp.simplify(sp.re(rotated_total_rate) - (p_0 - omega * y_0)) != 0:
        issues.append("moving-projector real derivative failed")
    if sp.simplify(sp.im(rotated_total_rate) - (q_0 + omega * x_0)) != 0:
        issues.append("moving-projector imaginary derivative failed")

    c, u, c_l, u_l = sp.symbols("c u c_L u_L", real=True)
    alpha = c * u
    alpha_l = c_l * u + c * u_l
    if sp.expand(alpha_l - (c_l * u + c * u_l)) != 0:
        issues.append("terminal-shear product rule failed")

    s, ell = sp.symbols("s L", real=True, positive=True)
    shear = sp.Matrix([[1, 0], [s * alpha / ell, 1]])
    if sp.simplify(shear.det() - 1) != 0:
        issues.append("terminal-shear determinant failed")

    b, b_l, p_b, q_b, p_b_l = sp.symbols(
        "b b_L P_B Q_B P_B_L", real=True
    )
    x_l = p_0 - omega * y_0
    y_l = q_0 + omega * x_0
    projected_b_l = p_b_l - omega * q_b
    mathsf_a = p_b - b * u * y_0 + c * u * x_0
    mathsf_a_l = (
        projected_b_l
        - (b_l * u + b * u_l) * y_0
        - b * u * y_l
        + alpha_l * x_0
        + alpha * x_l
    )
    abel_from_shear = sp.expand(mathsf_a_l - alpha_l * x_0 - alpha * x_l)
    abel_direct = sp.expand(
        projected_b_l
        - (b_l * u + b * u_l) * y_0
        - b * u * y_l
    )
    if sp.simplify(abel_from_shear - abel_direct) != 0:
        issues.append("direct and shear Abel derivatives disagree")


def validate_sources(
    stored: dict, sources: dict[str, dict], issues: list[str]
) -> None:
    expected_hashes = {
        key: file_hash(path)
        for key, path in SOURCES.items()
        if path.is_file()
    }
    if stored.get("source_audit", {}).get("source_sha256") != expected_hashes:
        issues.append("source hashes drifted")

    rectangular = sources["rectangular_boundary"].get("exact", {})
    if rectangular.get("ray_bottom", {}).get("hard_interval") != (
        "max(B_epsilon,sqrt((100+j)/50))<=L"
        "<=min(L_j,(c_*+epsilon)*(100+j)/25)"
    ):
        issues.append("source hard interval drifted")

    normalized = sources["normalized_phase_flux"].get("exact", {})
    if "gamma_n=q_(n,x)/q_n" not in normalized.get(
        "coefficient_current", ""
    ):
        issues.append("source carrier current drifted")
    if "Z_(A,x)=" not in normalized.get("normalized_flux_derivative", ""):
        issues.append("source normalized derivative drifted")

    shear = sources["abel_shear"].get("exact", {})
    if "2.03e-14" not in shear.get("orientation_ratio", ""):
        issues.append("source orientation ratio drifted")
    if "partial_s arg(Psi_ell)" not in shear.get("general_flux", ""):
        issues.append("source proxy flux drifted")

    connector = sources["connector_cap"].get("exact", {})
    if "H_j<3*pi/2" not in connector.get(
        "reduced_horizontal_target", ""
    ):
        issues.append("source horizontal target drifted")


def validate() -> list[str]:
    issues: list[str] = []
    stored = load_json(RESULT, "reduction result", issues)
    sources = {
        key: load_json(path, f"{key} source", issues)
        for key, path in SOURCES.items()
    }
    if issues:
        return issues

    if stored.get("kind") != STEM:
        issues.append("result kind drifted")
    if stored.get("date") != "2026-07-29":
        issues.append("result date drifted")
    if [row.get("id") for row in stored.get("rows", [])] != EXPECTED_IDS:
        issues.append("row ids or order drifted")

    expected_summary = {
        "rows": 24,
        "exact_cutoff_partitions": 1,
        "exact_logarithmic_flow_identities": 8,
        "exact_phase_flux_identities": 2,
        "orientation_transfers": 1,
        "cutoff_join_contracts": 1,
        "nonpromotion_guards": 2,
        "open_Xi_targets": 2,
        "proved_abel_gaps": 0,
        "proved_horizontal_phase_bounds": 0,
        "pointwise_contact_exclusions": 0,
    }
    if stored.get("summary") != expected_summary:
        issues.append("summary drifted")

    exact = stored.get("exact", {})
    if exact.get("logarithmic_derivative", {}).get("operator") != (
        "D_j=d/dL|_(t=t_j,N fixed)=x*partial_x, "
        "x=4*pi*exp(L)>0"
    ):
        issues.append("stored logarithmic operator drifted")
    if exact.get("ray_proxy", {}).get("orientation_defect") != (
        "E_j=D_j mathsf_X-x*mathcal_C_N"
    ):
        issues.append("stored orientation defect drifted")
    if "mathsf_X*(L*D_j mathcal_C_N-mathcal_C_N)" not in (
        exact.get("ray_proxy", {}).get("defect_flux", "")
    ):
        issues.append("stored endpoint-complete defect drifted")
    if "H_j<3*pi/2" not in (
        exact.get("conditional_targets", {}).get("joined_horizontal", "")
    ):
        issues.append("stored joined target drifted")

    boundary = stored.get("proof_boundary", "")
    for required in (
        "does not prove the Xi Abel gap",
        "aggregate defect",
        "H_j<3*pi/2",
        "contact exclusion",
        "RH",
    ):
        if required not in boundary:
            issues.append(f"proof boundary missing: {required}")

    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
    else:
        note = NOTE.read_text(encoding="utf-8")
        for required in (
            "# First-Order Ray-Bottom Logarithmic-Flow Reduction",
            "## Ray And Cutoff Cells",
            "## Endpoint-Complete Derivatives",
            "## Signed Ray Flux",
            "0 Abel gaps",
            "0 horizontal phase bounds",
        ):
            if required not in note:
                issues.append(f"note marker missing: {required}")

    validate_ray_geometry(issues)
    validate_flux_algebra(issues)
    validate_projector_and_shear(issues)
    validate_sources(stored, sources, issues)
    return issues


def main() -> int:
    issues = validate()
    for issue in issues:
        print(f"ERROR: {issue}")
    if issues:
        return 1
    print(
        "validated Newman centered ray-bottom logarithmic-flow reduction: "
        "24 rows, 1 exact cutoff partition, 8 logarithmic-flow identities, "
        "2 phase-flux identities, 1 orientation transfer, 1 cutoff-join "
        "contract, 2 nonpromotion guards, 2 open Xi targets, "
        "0 Abel gaps, 0 horizontal phase bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
