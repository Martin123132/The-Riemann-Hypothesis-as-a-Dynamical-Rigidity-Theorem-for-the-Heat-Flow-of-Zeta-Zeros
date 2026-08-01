#!/usr/bin/env python3
"""Build the corrected-component Wronskian reduction and resonance guard."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

import jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout as base
import jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction as wronskian_gate


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_polymath15_critical_component_wronskian_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_newman_polymath15_critical_component_wronskian_gate.md"
)
COARSE_DPS = 55
FINE_DPS = 75
ROOT_SPECS = (
    (("517.1", "517.3"), 6, "ordinary_negative"),
    (("519.6", "519.9"), 6, "ordinary_positive"),
    (("14010.12", "14010.13"), 33, "lehmer_left"),
    (("14010.19", "14010.21"), 33, "lehmer_right"),
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def build_exact() -> dict:
    a, b, u_1, u_2, v_1, v_2, delta = sp.symbols(
        "a b u_1 u_2 v_1 v_2 delta", real=True
    )
    e_1 = a * sp.exp(sp.I * delta)
    e_2 = b
    e_1_prime = (u_1 + sp.I * v_1) * e_1
    e_2_prime = (u_2 + sp.I * v_2) * e_2
    direct_pair = sp.simplify(
        sp.im((e_1_prime + e_2_prime) * sp.conjugate(e_1 + e_2))
    )
    expected_pair = (
        a**2 * v_1
        + b**2 * v_2
        + a
        * b
        * ((u_1 - u_2) * sp.sin(delta) + (v_1 + v_2) * sp.cos(delta))
    )
    if sp.trigsimp(direct_pair - expected_pair) != 0:
        raise RuntimeError("two-component Wronskian identity failed")

    x, y, r_x, r_y, u_star, v_star = sp.symbols(
        "x y r_x r_y u_star v_star", real=True
    )
    value = x + sp.I * y
    residual = r_x + sp.I * r_y
    first = (u_star + sp.I * v_star) * value + residual
    reference_wronskian = sp.simplify(sp.im(first * sp.conjugate(value)))
    expected_reference = v_star * (x**2 + y**2) + r_y * x - r_x * y
    if sp.simplify(reference_wronskian - expected_reference) != 0:
        raise RuntimeError("reference-rate Wronskian identity failed")
    crossing_first = sp.simplify(sp.re(first).subs(x, 0))
    if sp.simplify(crossing_first - (r_x - v_star * y)) != 0:
        raise RuntimeError("reference-rate crossing identity failed")

    ell, beta_first, s_real, s_imag = sp.symbols(
        "ell beta_first s_real s_imag", real=True
    )
    component_rate = sp.I * beta_first - (s_real + sp.I * s_imag) * ell
    if sp.re(component_rate) != -s_real * ell:
        raise RuntimeError("Dirichlet radial-rate identity failed")
    if sp.im(component_rate) != beta_first - s_imag * ell:
        raise RuntimeError("Dirichlet phase-rate identity failed")

    u_3, v_3 = sp.symbols("u_3 v_3", real=True)
    rates = [
        u_1 + sp.I * v_1,
        u_2 + sp.I * v_2,
        u_3 + sp.I * v_3,
    ]
    rate_matrix = sp.Matrix(
        3,
        3,
        lambda row, column: (
            rates[column] - sp.conjugate(rates[row])
        )
        / (2 * sp.I),
    )
    if any(
        sp.simplify(
            rate_matrix[row, column]
            - sp.conjugate(rate_matrix[column, row])
        )
        != 0
        for row in range(3)
        for column in range(3)
    ):
        raise RuntimeError("component-rate matrix is not Hermitian")
    if sp.simplify(rate_matrix.det()) != 0:
        raise RuntimeError("three-component rate matrix is not rank at most two")
    rate_sum = sum(rates, sp.Integer(0))
    rate_square_sum = sum(
        (
            sp.re(rate) ** 2 + sp.im(rate) ** 2
            for rate in rates
        ),
        sp.Integer(0),
    )
    second_symmetric = sp.simplify(
        (
            sp.trace(rate_matrix) ** 2
            - sp.trace(rate_matrix * rate_matrix)
        )
        / 2
    )
    expected_second_symmetric = sp.simplify(
        -(
            3 * rate_square_sum
            - (
                sp.re(rate_sum) ** 2
                + sp.im(rate_sum) ** 2
            )
        )
        / 4
    )
    if (
        sp.simplify(
            second_symmetric - expected_second_symmetric
        )
        != 0
    ):
        raise RuntimeError("component-rate matrix inertia identity failed")

    variable = sp.symbols("variable", real=True)
    toy = (
        sp.exp(-3 * sp.I * variable)
        - sp.exp(-4 * sp.I * variable)
        + sp.I * sp.exp(-sp.I * variable)
        - sp.I * sp.exp(-2 * sp.I * variable) / 2
    )
    toy_value = sp.simplify(toy.subs(variable, 0))
    toy_first = sp.simplify(sp.diff(toy, variable).subs(variable, 0))
    toy_second_real = sp.simplify(
        sp.re(sp.diff(toy, variable, 2).subs(variable, 0))
    )
    toy_wronskian = sp.simplify(
        sp.im(toy_first * sp.conjugate(toy_value))
    )
    if (
        toy_value != sp.I / 2
        or toy_first != sp.I
        or toy_second_real != 7
        or toy_wronskian != 0
    ):
        raise RuntimeError("ordered-speed double-crossing countermodel failed")

    return {
        "component_split": (
            "Write E=sum_(j=0)^N e_j, where e_n=exp(i*beta_t)"
            "*exp((t/4)log(n)^2-s_*(x)log(n)) for 1<=n<=N and "
            "e_0=-q_(N,t) is the sharp endpoint pseudo-component"
        ),
        "component_rates": (
            "For e_j=a_j exp(i*theta_j), write e_j'=(u_j+i*v_j)e_j, "
            "u_j=(log a_j)' and v_j=theta_j'"
        ),
        "pairwise_wronskian": (
            "W_E=Im(E'*conj(E))=sum_j a_j^2 v_j"
            "+sum_(j<k)a_j*a_k*((u_j-u_k)sin(theta_j-theta_k)"
            "+(v_j+v_k)cos(theta_j-theta_k))"
        ),
        "wronskian_matrix": (
            "For z=(e_0,...,e_N)^T and r_j=u_j+i*v_j, "
            "W_E=z^*Kz with K_(jk)=(r_k-conj(r_j))/(2i); "
            "K is Hermitian and rank(K)<=2"
        ),
        "matrix_inertia": (
            "If m=N+1, the product of the two possible nonzero eigenvalues is "
            "lambda_+*lambda_-=-(m*sum_j|r_j|^2-|sum_j r_j|^2)/4"
            "=-(m/4)*sum_j|r_j-r_avg|^2; unequal component rates give "
            "inertia (1,1,m-2)"
        ),
        "dirichlet_rates": (
            "For 1<=n<=N, u_n=-Re(s_*')log(n) and "
            "v_n=beta_t'-Im(s_*')log(n); at t=0, "
            "u_n=0 and v_n=beta_0'+log(n)/2"
        ),
        "reference_rate": (
            "For any real u_* and v_*, set "
            "R_*=E'-(u_*+i*v_*)E. Then "
            "W_E=v_*|E|^2+Im(R_*conj(E))"
        ),
        "reference_crossing": (
            "At Re(E)=0 with E=iY, Re(E')=Re(R_*)-v_*Y; "
            "a nonzero double crossing is equivalent to Re(R_*)=v_*Y"
        ),
        "corrected_crossing_system": (
            "Re(E)=sum_j a_j cos(theta_j)=0 and "
            "Re(E')=sum_j a_j*(u_j cos(theta_j)-v_j sin(theta_j))=0"
        ),
        "ordered_speed_countermodel": (
            "E_toy(x)=exp(-3ix)-exp(-4ix)+i*exp(-ix)"
            "-(i/2)exp(-2ix) has positive component amplitudes and distinct "
            "negative phase speeds {-4,-3,-2,-1}, but "
            "E_toy(0)=i/2, E_toy'(0)=i, Re(E_toy) has a genuine double zero "
            "with second derivative 7, and W_Etoy(0)=0"
        ),
        "small_ball_target": (
            "Prove X^2+(U/L)^2>8000000*exp(-3L/2), "
            "X=sum_j a_j cos(theta_j), "
            "U=sum_j a_j*(u_j cos(theta_j)-v_j sin(theta_j)), "
            "on L>=50 and 0<tL<=c_*+o(1)"
        ),
        "proof_handoff": (
            "The missing input is a joint arithmetic small-ball or "
            "phase-critical-value exclusion for the logarithmic Riemann-Siegel "
            "phases; componentwise phase ordering, absolute upper bounds, and "
            "a positive-semidefinite argument based only on the instantaneous "
            "component-rate matrix cannot supply it"
        ),
    }


def corrected_components(
    x: mp.mpf, time: mp.mpf, cutoff: int
) -> list[mp.mpc]:
    s = (1 - mp.j * x) / 2
    alpha = base.alpha(s)
    shifted_s = s + time * alpha / 2
    normalizer = base.m_time(s, time)
    phase = normalizer / abs(normalizer)
    components = [
        phase
        * mp.exp(
            time * mp.log(n) ** 2 / 4 - shifted_s * mp.log(n)
        )
        for n in range(1, cutoff + 1)
    ]

    saddle_time = x / 2 + mp.pi * time / 8
    saddle = mp.sqrt(saddle_time / (2 * mp.pi))
    p = 1 - 2 * (saddle - cutoff)
    u_rs = mp.exp(
        -mp.j
        * (
            saddle_time / 2 * mp.log(saddle_time / (2 * mp.pi))
            - saddle_time / 2
            - mp.pi / 8
        )
    )
    endpoint_half = (
        (-1) ** cutoff
        * mp.exp(time * mp.pi**2 / 64)
        * base.m_zero(mp.j * saddle_time)
        * base.c_zero(p)
        * u_rs
        * mp.exp(mp.j * mp.pi / 8)
        / abs(normalizer)
    )
    components.append(-endpoint_half)
    return components


def component_audit(
    root: mp.mpf, time: mp.mpf, cutoff: int, label: str
) -> dict:
    components = corrected_components(root, time, cutoff)
    derivatives = [
        mp.diff(
            lambda value, index=index: corrected_components(
                value, time, cutoff
            )[index],
            root,
        )
        for index in range(len(components))
    ]
    value = mp.fsum(components)
    first = mp.fsum(derivatives)
    direct_wronskian = mp.im(first * mp.conj(value))
    rates = [derivative / component for component, derivative in zip(
        components, derivatives, strict=True
    )]
    radial_rates = [mp.re(rate) for rate in rates]
    phase_speeds = [mp.im(rate) for rate in rates]

    diagonal_terms = [
        abs(component) ** 2 * phase_speed
        for component, phase_speed in zip(
            components, phase_speeds, strict=True
        )
    ]
    radial_pair_terms: list[mp.mpf] = []
    phase_pair_terms: list[mp.mpf] = []
    for left in range(len(components)):
        for right in range(left + 1, len(components)):
            amplitude_product = abs(components[left]) * abs(components[right])
            phase_difference = mp.arg(components[left]) - mp.arg(
                components[right]
            )
            radial_pair_terms.append(
                amplitude_product
                * (radial_rates[left] - radial_rates[right])
                * mp.sin(phase_difference)
            )
            phase_pair_terms.append(
                amplitude_product
                * (phase_speeds[left] + phase_speeds[right])
                * mp.cos(phase_difference)
            )
    diagonal_sum = mp.fsum(diagonal_terms)
    radial_pair_sum = mp.fsum(radial_pair_terms)
    phase_pair_sum = mp.fsum(phase_pair_terms)
    reconstructed = diagonal_sum + radial_pair_sum + phase_pair_sum
    scale = max(mp.mpf(1), abs(direct_wronskian))
    identity_error = abs(reconstructed - direct_wronskian) / scale
    all_terms = diagonal_terms + radial_pair_terms + phase_pair_terms
    cancellation_factor = mp.fsum(abs(term) for term in all_terms) / abs(
        direct_wronskian
    )
    dirichlet_speeds = phase_speeds[:-1]
    return {
        "label": label,
        "time": mp.nstr(time, 20, strip_zeros=False),
        "cutoff": cutoff,
        "component_count": len(components),
        "root": mp.nstr(root, 40, strip_zeros=False),
        "abs_corrected_crossing_residual": mp.nstr(
            abs(2 * mp.re(value)), 25, strip_zeros=False
        ),
        "corrected_first_derivative": mp.nstr(
            2 * mp.re(first), 35, strip_zeros=False
        ),
        "direct_wronskian": mp.nstr(
            direct_wronskian, 35, strip_zeros=False
        ),
        "aggregate_phase_velocity": mp.nstr(
            direct_wronskian / abs(value) ** 2, 35, strip_zeros=False
        ),
        "diagonal_sum": mp.nstr(diagonal_sum, 35, strip_zeros=False),
        "radial_pair_sum": mp.nstr(
            radial_pair_sum, 35, strip_zeros=False
        ),
        "phase_pair_sum": mp.nstr(
            phase_pair_sum, 35, strip_zeros=False
        ),
        "pairwise_identity_relative_error": mp.nstr(
            identity_error, 15, strip_zeros=False
        ),
        "absolute_term_cancellation_factor": mp.nstr(
            cancellation_factor, 25, strip_zeros=False
        ),
        "all_component_phase_speeds_negative": all(
            speed < 0 for speed in phase_speeds
        ),
        "dirichlet_phase_speeds_strictly_increasing": all(
            left < right
            for left, right in zip(
                dirichlet_speeds, dirichlet_speeds[1:]
            )
        ),
        "endpoint_is_maximum_phase_speed": (
            phase_speeds[-1] == max(phase_speeds)
        ),
        "minimum_component_phase_speed": mp.nstr(
            min(phase_speeds), 25, strip_zeros=False
        ),
        "maximum_component_phase_speed": mp.nstr(
            max(phase_speeds), 25, strip_zeros=False
        ),
        "endpoint_radial_rate": mp.nstr(
            radial_rates[-1], 25, strip_zeros=False
        ),
        "endpoint_phase_speed": mp.nstr(
            phase_speeds[-1], 25, strip_zeros=False
        ),
    }


def diagnostics(dps: int) -> dict:
    mp.mp.dps = dps
    time = mp.mpf(0)
    rows = []
    for bracket, cutoff, label in ROOT_SPECS:
        crossing = lambda value: 2 * mp.re(
            wronskian_gate.corrected_complex_main(value, time, cutoff)
        )
        root = mp.findroot(
            crossing, tuple(mp.mpf(endpoint) for endpoint in bracket)
        )
        rows.append(component_audit(root, time, cutoff, label))
    if not all(row["all_component_phase_speeds_negative"] for row in rows):
        raise RuntimeError("component phase-speed sign diagnostic drifted")
    if not all(
        row["dirichlet_phase_speeds_strictly_increasing"] for row in rows
    ):
        raise RuntimeError("Dirichlet phase-speed ordering drifted")
    if not all(row["endpoint_is_maximum_phase_speed"] for row in rows):
        raise RuntimeError("endpoint phase-speed ordering drifted")
    wronskians = [mp.mpf(row["direct_wronskian"]) for row in rows]
    if not (min(wronskians) < 0 < max(wronskians)):
        raise RuntimeError("aggregate Wronskian sign diagnostic drifted")
    maximum_cancellation = max(
        mp.mpf(row["absolute_term_cancellation_factor"]) for row in rows
    )
    if maximum_cancellation <= 500:
        raise RuntimeError("Lehmer cancellation stress unexpectedly weak")
    return {
        "role": "finite_route_shaping_diagnostics_only",
        "proof_boundary": (
            "The four t=0 corrected crossings audit the exact pairwise identity "
            "and expose severe cancellation. They do not prove any L>=50 "
            "small-ball bound or a positive-time theorem."
        ),
        "dps": dps,
        "rows": rows,
        "maximum_absolute_term_cancellation_factor": mp.nstr(
            maximum_cancellation, 25, strip_zeros=False
        ),
    }


def compare_diagnostics(coarse: dict, fine: dict) -> dict:
    maximum_root_delta = mp.mpf(0)
    maximum_wronskian_delta = mp.mpf(0)
    maximum_cancellation_delta = mp.mpf(0)
    for left, right in zip(coarse["rows"], fine["rows"], strict=True):
        if left["label"] != right["label"]:
            raise RuntimeError("component diagnostic ordering drifted")
        maximum_root_delta = max(
            maximum_root_delta,
            abs(mp.mpf(left["root"]) - mp.mpf(right["root"])),
        )
        maximum_wronskian_delta = max(
            maximum_wronskian_delta,
            abs(
                mp.mpf(left["direct_wronskian"])
                - mp.mpf(right["direct_wronskian"])
            ),
        )
        maximum_cancellation_delta = max(
            maximum_cancellation_delta,
            abs(
                mp.mpf(left["absolute_term_cancellation_factor"])
                - mp.mpf(right["absolute_term_cancellation_factor"])
            ),
        )
    if (
        maximum_root_delta >= mp.mpf("1e-35")
        or maximum_wronskian_delta >= mp.mpf("1e-30")
        or maximum_cancellation_delta >= mp.mpf("1e-20")
    ):
        raise RuntimeError("component Wronskian diagnostics are unstable")
    return {
        "maximum_abs_root_delta": mp.nstr(
            maximum_root_delta, 25, strip_zeros=False
        ),
        "maximum_abs_wronskian_delta": mp.nstr(
            maximum_wronskian_delta, 25, strip_zeros=False
        ),
        "maximum_abs_cancellation_factor_delta": mp.nstr(
            maximum_cancellation_delta, 25, strip_zeros=False
        ),
    }


def build_artifact() -> dict:
    exact = build_exact()
    coarse = diagnostics(COARSE_DPS)
    fine = diagnostics(FINE_DPS)
    convergence = compare_diagnostics(coarse, fine)
    rows = [
        GateRow(
            id="np15cwg_01_component_split",
            role="exact_definition",
            readiness="ready_to_apply",
            claim="The corrected complex main is one finite component sum after treating the endpoint lift as a pseudo-component.",
            formula=f"{exact['component_split']}; {exact['component_rates']}",
            proof_boundary="One fixed analytic cutoff cell on the real axis.",
        ),
        GateRow(
            id="np15cwg_02_pairwise_wronskian",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The corrected phase Wronskian has an exact diagonal-plus-pair decomposition.",
            formula=(
                f"{exact['pairwise_wronskian']}; "
                f"{exact['wronskian_matrix']}; "
                f"{exact['matrix_inertia']}"
            ),
            proof_boundary=(
                "The instantaneous component-rate matrix is generically "
                "rank-two indefinite, so no sign is asserted."
            ),
        ),
        GateRow(
            id="np15cwg_03_dirichlet_rates",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Every Dirichlet component has an explicit radial rate and phase speed.",
            formula=exact["dirichlet_rates"],
            proof_boundary="The endpoint pseudo-component has its own explicit rates.",
        ),
        GateRow(
            id="np15cwg_04_reference_rate",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Subtracting any common complex rate isolates the component defects without changing the exact Wronskian bookkeeping.",
            formula=f"{exact['reference_rate']}; {exact['reference_crossing']}",
            proof_boundary="An identity, not a bound on the residual.",
        ),
        GateRow(
            id="np15cwg_05_crossing_system",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="A corrected double crossing is a simultaneous two-observable trigonometric cancellation.",
            formula=exact["corrected_crossing_system"],
            proof_boundary="The certified remainder still requires the robust collar margin.",
        ),
        GateRow(
            id="np15cwg_06_ordered_speed_guard",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Positive amplitudes and strictly ordered same-sign component phase speeds do not exclude a nonzero double crossing.",
            formula=exact["ordered_speed_countermodel"],
            proof_boundary="Finite exponential-polynomial countermodel; not an Xi counterexample.",
        ),
        GateRow(
            id="np15cwg_07_corrected_diagnostics",
            role="finite_diagnostics",
            readiness="diagnostic_only",
            claim="Corrected ordinary and Lehmer-stress crossings verify the pair identity and exhibit both aggregate Wronskian signs despite ordered negative component speeds.",
            formula="four t=0 corrected crossings at N=6 and N=33",
            proof_boundary=fine["proof_boundary"],
            diagnostics={
                "convergence": convergence,
                "rows": fine["rows"],
                "maximum_absolute_term_cancellation_factor": fine[
                    "maximum_absolute_term_cancellation_factor"
                ],
            },
        ),
        GateRow(
            id="np15cwg_08_cancellation_guard",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Absolute component bounds and phase-speed ordering lose the pair cancellation that controls corrected transversality.",
            formula=exact["proof_handoff"],
            proof_boundary="The finite cancellation factors do not establish an asymptotic rate.",
        ),
        GateRow(
            id="np15cwg_09_small_ball_target",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim="The high-frequency problem is a joint arithmetic small-ball exclusion for the corrected value and first derivative.",
            formula=exact["small_ball_target"],
            proof_boundary="Open on the residual positive scaled-time layer.",
        ),
        GateRow(
            id="np15cwg_10_proof_boundary",
            role="proof_guard",
            readiness="guard_validated",
            claim="The component reduction and diagnostics do not prove corrected transversality.",
            formula="exact pairwise reduction plus finite diagnostics != Lambda<=0",
            proof_boundary="No positive-time simplicity, Lambda<=0, RH, or Clay-prize result.",
        ),
    ]
    return {
        "kind": (
            "jensen_window_pf_newman_polymath15_"
            "critical_component_wronskian_gate"
        ),
        "date": "2026-07-24",
        "status": (
            "exact corrected-component Wronskian reduction with an ordered-speed "
            "countermodel and finite cancellation diagnostics; the arithmetic "
            "small-ball theorem remains open"
        ),
        "proof_boundary": (
            "This artifact proves the component, pairwise Wronskian, rate, "
            "reference-frame, and crossing identities and validates an exact "
            "same-sign ordered-speed countermodel. Its numerical rows are "
            "route-shaping diagnostics only. It does not prove the residual "
            "small-ball target, corrected transversality, positive-time "
            "simplicity, Lambda<=0, RH, or a Clay-prize result."
        ),
        "exact": exact,
        "diagnostics": fine,
        "convergence": convergence,
        "rows": [asdict(row) for row in rows],
        "sources": [
            "outputs/jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_transversality_target.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.md",
            "https://arxiv.org/abs/1904.12438",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    diagnostics_payload = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF Newman Polymath-15 Critical Component Wronskian Gate",
        "",
        "Date: 2026-07-24",
        "",
        "Status: exact corrected-component Wronskian reduction, exact route",
        "guard, and finite cancellation diagnostics. The arithmetic small-ball",
        "theorem remains open; this is not a proof of `Lambda <= 0` or RH.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_component_wronskian_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_component_wronskian_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_component_wronskian_gate.py",
        "```",
        "",
        "## Pairwise Identity",
        "",
        "```text",
        exact["component_split"],
        exact["component_rates"],
        exact["pairwise_wronskian"],
        "```",
        "",
        "Equivalently,",
        "",
        "```text",
        exact["wronskian_matrix"],
        exact["matrix_inertia"],
        "```",
        "",
        "Thus the unrestricted instantaneous rate form has one positive and",
        "one negative direction whenever the component rates are not all",
        "equal. A positive-semidefinite Gram argument cannot arise from this",
        "matrix alone; arithmetic phase placement or additional Xi structure",
        "must remove its negative direction on the relevant crossing set.",
        "",
        "For the actual Dirichlet components,",
        "",
        "```text",
        exact["dirichlet_rates"],
        "```",
        "",
        "At `t=0` the Dirichlet amplitudes have zero radial rate and their phase",
        "speeds increase toward the saddle endpoint. The endpoint correction is",
        "one additional explicit component, so no cross term is discarded.",
        "",
        "## Reference Rate",
        "",
        "```text",
        exact["reference_rate"],
        exact["reference_crossing"],
        "```",
        "",
        "This makes endpoint-relative formulations exact, but it does not make",
        "the defect term one-signed.",
        "",
        "## Double-Crossing System",
        "",
        "```text",
        exact["corrected_crossing_system"],
        "```",
        "",
        "The two equations require joint cancellation of different trigonometric",
        "observables. Upper bounds for either sum separately do not exclude their",
        "simultaneous smallness.",
        "",
        "## Ordered-Speed Countermodel",
        "",
        "```text",
        exact["ordered_speed_countermodel"],
        "```",
        "",
        "Thus positive component amplitudes, distinct ordered phase speeds, and",
        "a common speed sign are not a transversality theorem.",
        "",
        "## Corrected Diagnostics",
        "",
        "| label | N | crossing | W_E | phase speed | diagonal | pair radial | pair phase | cancellation |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in diagnostics_payload["rows"]:
        lines.append(
            "| {label} | {cutoff} | {root} | {direct_wronskian} | "
            "{aggregate_phase_velocity} | {diagonal_sum} | "
            "{radial_pair_sum} | {phase_pair_sum} | "
            "{absolute_term_cancellation_factor} |".format(**row)
        )
    lines.extend(
        [
            "",
            "All component phase speeds in these four rows are negative, the",
            "Dirichlet speeds are strictly increasing, and the endpoint speed is",
            "the largest. Nevertheless the aggregate Wronskian has both signs.",
            "The largest absolute-term cancellation factor is",
            f"`{diagnostics_payload['maximum_absolute_term_cancellation_factor']}`.",
            "These are finite `t=0` diagnostics, not an asymptotic theorem.",
            "",
            "## Live Target",
            "",
            "```text",
            exact["small_ball_target"],
            "```",
            "",
            exact["proof_handoff"] + ".",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 critical component Wronskian gate: "
        f"{len(artifact['rows'])} rows, 5 exact identities/reductions, "
        "1 exact ordered-speed countermodel, 4 corrected diagnostics, "
        "1 open arithmetic small-ball target"
    )


if __name__ == "__main__":
    main()
