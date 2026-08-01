#!/usr/bin/env python3
"""Independently validate the first-order target reconciliation gate."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_target_reconciliation_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
SOURCES = {
    "zeroth_order_remainder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_C1_global_remainder_certificate.json"
    ),
    "first_order_remainder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_first_order_global_remainder_certificate.json"
    ),
    "first_order_signed_contact": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_first_order_signed_contact_reduction.json"
    ),
    "wronskian_phase": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_wronskian_phase_reduction.json"
    ),
    "component_wronskian": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_component_wronskian_gate.json"
    ),
    "absolute_phase_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_absolute_phase_anchor_reduction.json"
    ),
    "real_residual": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_real_residual_reduction.json"
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
    "oscillatory_splice": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_oscillatory_spliced_outer_collar_reduction.json"
    ),
}
EXPECTED_IDS = [
    "fotr_01_chronology",
    "fotr_02_current_main",
    "fotr_03_current_remainder",
    "fotr_04_local_lift",
    "fotr_05_cutoff_quantifiers",
    "fotr_06_cartesian_energy",
    "fotr_07_phase_amplitude",
    "fotr_08_phase_energy",
    "fotr_09_cross_terms",
    "fotr_10_rank_cancellation_guard",
    "fotr_11_contact_box",
    "fotr_12_radial_target",
    "fotr_13_threshold_improvement",
    "fotr_14_box_gauge",
    "fotr_15_signed_band",
    "fotr_16_radial_overstrength",
    "fotr_17_absolute_anchor",
    "fotr_18_abel_shear",
    "fotr_19_abel_implies_band",
    "fotr_20_low_c_splice",
    "fotr_21_inner_guard",
    "fotr_22_route_decision",
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


def validate_phase_identity(issues: list[str]) -> None:
    count = 3
    amplitudes = sp.symbols("a0:3", real=True)
    phases = sp.symbols("theta0:3", real=True)
    radial = sp.symbols("u0:3", real=True)
    angular = sp.symbols("v0:3", real=True)
    ell = sp.symbols("L", positive=True)

    c = [
        amplitudes[j] * sp.cos(phases[j])
        for j in range(count)
    ]
    d = [
        amplitudes[j]
        * (
            radial[j] * sp.cos(phases[j])
            - angular[j] * sp.sin(phases[j])
        )
        for j in range(count)
    ]
    compact = sp.expand(sum(c) ** 2 + sum(d) ** 2 / ell**2)
    expanded = sum(
        c[j] ** 2 + d[j] ** 2 / ell**2
        for j in range(count)
    ) + 2 * sum(
        c[j] * c[k] + d[j] * d[k] / ell**2
        for j in range(count)
        for k in range(j + 1, count)
    )
    if sp.simplify(compact - expanded) != 0:
        issues.append("phase-energy cross-term expansion failed")

    a, theta, u, v = sp.symbols("a theta u v", real=True)
    z = a * sp.exp(sp.I * theta)
    modeled_derivative = (u + sp.I * v) * z
    expected_real = a * (
        u * sp.cos(theta) - v * sp.sin(theta)
    )
    if sp.simplify(sp.re(modeled_derivative).expand(complex=True) - expected_real) != 0:
        issues.append("component derivative phase identity failed")


def validate_thresholds(issues: list[str]) -> None:
    delta = 50_000
    gamma = 100_000
    radial = 4 * (delta**2 + gamma**2)
    if radial != 50_000_000_000:
        issues.append("first-order radial constant failed")
    ratio = Fraction(radial, 32_000_000)
    if ratio != Fraction(3125, 2):
        issues.append("old/new threshold ratio failed")
    if Fraction(3125, 2) >= 2**50:
        issues.append("L=50 exponential-improvement guard failed")

    x = Fraction(0)
    derivative = Fraction(21, 10)
    if max(abs(x), abs(derivative) / 2) <= 1:
        issues.append("box-gauge overstrength witness failed")
    if x**2 + derivative**2 >= 5:
        issues.append("radial overstrength witness failed")

    positive_slack = Fraction(1) - Fraction(1, 10**6)
    if positive_slack <= 0:
        issues.append("Abel-to-band residual slack failed")


def validate_cancellation_guard(issues: list[str]) -> None:
    amplitude = sp.symbols("A", nonzero=True, real=True)
    w_1 = sp.Matrix([amplitude, 0])
    w_2 = sp.Matrix([-amplitude, 0])
    diagonal = (w_1.dot(w_1) + w_2.dot(w_2)).simplify()
    aggregate = ((w_1 + w_2).dot(w_1 + w_2)).simplify()
    if diagonal != 2 * amplitude**2 or aggregate != 0:
        issues.append("rank-two cancellation guard failed")


def validate_sources(
    stored: dict, source_payloads: dict[str, dict], issues: list[str]
) -> None:
    stored_hashes = (
        stored.get("source_audit", {}).get("source_sha256", {})
    )
    expected_hashes = {
        key: file_hash(path)
        for key, path in SOURCES.items()
        if path.is_file()
    }
    if stored_hashes != expected_hashes:
        issues.append("source hashes drifted")

    old = source_payloads["zeroth_order_remainder"]
    if old.get("exact", {}).get("remaining_target") != (
        "Prove T_L[J]>32000000*exp(-3L/2) "
        "for the corrected finite main"
    ):
        issues.append("historical target drifted")

    first = source_payloads["first_order_remainder"]
    if first.get("exact", {}).get("global_remainder") != (
        "For every critical disk, including cutoff crossings, "
        "|r_[1](x)|<100000*exp(-5L/4) and "
        "|partial_x r_[1](x)|<200000*L*exp(-5L/4)"
    ):
        issues.append("first-order remainder drifted")

    signed = source_payloads["first_order_signed_contact"]
    if signed.get("exact", {}).get("corrected_complex_main") != (
        "E_[1]=g_0+sum_(n=1)^N f_n, J_[1]=2*Re(E_[1])"
    ):
        issues.append("first-order main drifted")
    if signed.get("exact", {}).get("certified_contact_box") != (
        "|X_[1]|<50000*exp(-5L/4), "
        "|U_[1]|<100000*L*exp(-5L/4)"
    ):
        issues.append("first-order contact box drifted")

    phase = source_payloads["wronskian_phase"]
    if phase.get("exact", {}).get("cartesian") != (
        "For E=X+iY and E'=U+iV, "
        "T_L[J]=4*(X^2+(U/L)^2)"
    ):
        issues.append("Wronskian Cartesian identity drifted")

    component = source_payloads["component_wronskian"]
    if "genuine double zero" not in (
        component.get("exact", {}).get("ordered_speed_countermodel", "")
    ):
        issues.append("ordered-speed guard drifted")

    anchor = source_payloads["absolute_phase_anchor"]
    if "X/|f_1|=Re(eta*Z_0)" not in (
        anchor.get("exact", {}).get("real_projection", "")
    ):
        issues.append("absolute phase anchor drifted")

    residual = source_payloads["real_residual"]
    if residual.get("exact", {}).get("core_scalar", {}).get(
        "band_approximation"
    ) != (
        "For |X|<=50000*exp(-5L/4), "
        "|U-A_a|<1e-6*exp(-5L/4)."
    ):
        issues.append("real-residual band approximation drifted")

    prefix = source_payloads["abel_prefix"]
    if "mathsf_A=Re(W_A)=mathcal_C_N+c*u_N*mathsf_X" not in (
        prefix.get("exact", {}).get("contact_scalar", "")
    ):
        issues.append("Abel terminal shear drifted")
    if "delta_L=50000*exp(-5L/4)/|f_1|" not in (
        prefix.get("exact", {}).get("live_target", "")
    ):
        issues.append("Abel contact band drifted")

    shear = source_payloads["abel_shear"]
    if "|mathcal_C_N|>A_L+epsilon_term" not in (
        shear.get("exact", {}).get("q_ge_1_pointwise_target", "")
    ):
        issues.append("Abel pointwise target drifted")

    splice = source_payloads["oscillatory_splice"]
    outer = splice.get("exact", {}).get("open_outer_target", "")
    if "0<tL<=c_*+epsilon" not in outer or "q=2tL^2>=1" not in outer:
        issues.append("oscillatory-spliced outer target drifted")


def validate() -> list[str]:
    issues: list[str] = []
    stored = load_json(RESULT, "reconciliation result", issues)
    source_payloads = {
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
        "rows": 22,
        "superseded_handoffs": 1,
        "cutoff_uniform_first_order_remainders": 1,
        "exact_energy_identities": 4,
        "contact_box_coordinates": 2,
        "direct_radial_targets": 1,
        "box_optimal_targets": 1,
        "normalized_abel_targets": 1,
        "nonpromotion_guards": 2,
        "open_outer_abel_gaps": 1,
        "pointwise_contact_exclusions": 0,
        "cofinal_descendant_theorem": False,
    }
    if stored.get("summary") != expected_summary:
        issues.append("summary drifted")

    exact = stored.get("exact", {})
    if exact.get("division_free_components", {}).get("energy") != (
        "T_L[J_[1]]=4[X^2+(U/L)^2]"
    ):
        issues.append("stored Cartesian energy drifted")
    if exact.get("direct_radial_target", {}).get("constant_identity") != (
        "4*(50000^2+100000^2)=50000000000"
    ):
        issues.append("stored radial constant drifted")
    if exact.get("box_optimal_target", {}).get("exact_target") != (
        "G_L(X,U)>1, equivalently "
        "|X|<=delta_0(L) => |U|>gamma_0(L)"
    ):
        issues.append("stored box-optimal target drifted")
    if "implies the box-optimal signed band" not in (
        exact.get("abel_target", {}).get("conclusion", "")
    ):
        issues.append("stored Abel implication drifted")
    if exact.get("outer_domain_reduction", {}).get("live_outer_wedge") != (
        "L>=B_epsilon, q=2tL^2>=1, "
        "0<tL<=c_*+epsilon"
    ):
        issues.append("stored low-c wedge drifted")
    if "U=u_N*S+sum_k log((k+1)/k)F_k" not in (
        exact.get("open_handoff", "")
    ):
        issues.append("stored Abel-prefix handoff drifted")

    boundary = stored.get("proof_boundary", "")
    for required in (
        "does not prove the Xi Abel gap",
        "q<1",
        "Lambda<=0",
        "RH",
    ):
        if required not in boundary:
            issues.append(f"proof boundary missing: {required}")

    validate_phase_identity(issues)
    validate_thresholds(issues)
    validate_cancellation_guard(issues)
    validate_sources(stored, source_payloads, issues)
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Newman first-order target reconciliation gate: "
        "22 rows, 1 superseded handoff, "
        "1 cutoff-uniform first-order remainder, "
        "4 exact energy identities, 2 contact-box coordinates, "
        "1 radial target, 1 box-optimal target, "
        "1 normalized Abel target, 2 nonpromotion guards, "
        "0 pointwise contact exclusions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
