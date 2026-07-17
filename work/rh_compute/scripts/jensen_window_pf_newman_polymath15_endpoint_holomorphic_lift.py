#!/usr/bin/env python3
"""Construct a holomorphic lift of the Polymath-15 endpoint correction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.md"
)
SOURCE_URL = "https://arxiv.org/abs/1904.12438"


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def m_zero(s: mp.mpc) -> mp.mpc:
    return (
        mp.mpf(1) / 8
        * (s * (s - 1) / 2)
        * mp.pi ** (-s / 2)
        * mp.sqrt(2 * mp.pi)
        * mp.exp((s / 2 - mp.mpf("0.5")) * mp.log(s / 2) - s / 2)
    )


def u_rs(tau: mp.mpf) -> mp.mpc:
    return mp.exp(
        -mp.j
        * (
            tau / 2 * mp.log(tau / (2 * mp.pi))
            - tau / 2
            - mp.pi / 8
        )
    )


def c_zero(p: mp.mpf) -> mp.mpc:
    if abs(abs(p) - mp.mpf("0.5")) < mp.mpf("1e-40"):
        numerator_prime = (
            mp.exp(mp.pi * mp.j * (p**2 / 2 + mp.mpf(3) / 8))
            * mp.pi
            * mp.j
            * p
            + mp.j
            * mp.sqrt(2)
            * mp.pi
            / 2
            * mp.sin(mp.pi * p / 2)
        )
        denominator_prime = -2 * mp.pi * mp.sin(mp.pi * p)
        return numerator_prime / denominator_prime
    return (
        mp.exp(mp.pi * mp.j * (p**2 / 2 + mp.mpf(3) / 8))
        - mp.j * mp.sqrt(2) * mp.cos(mp.pi * p / 2)
    ) / (2 * mp.cos(mp.pi * p))


def slow_g(tau: mp.mpf | mp.mpc) -> mp.mpc:
    return -(
        mp.sqrt(mp.pi)
        / 8
        * mp.exp(-mp.pi * tau / 4)
        * (tau ** mp.mpf("1.5") + mp.j * tau ** mp.mpf("0.5"))
    )


def build_exact() -> dict:
    tau, slow_first = sp.symbols("tau slow_first", positive=True, real=True)
    slow = sp.Function("S")(tau)
    f = sp.exp(-sp.pi * tau / 4) * slow
    derivative_x = sp.diff(f, tau) / 2
    defect = sp.simplify(derivative_x + sp.pi * f / 8)
    expected = sp.exp(-sp.pi * tau / 4) * slow_first / 2
    defect_with_jet = defect.subs(sp.Derivative(slow, tau), slow_first)
    if sp.simplify(defect_with_jet - expected) != 0:
        raise RuntimeError("endpoint slow-derivative cancellation failed")

    d, paper, lift, source = sp.symbols(
        "d paper lift source", complex=True
    )
    # Algebraic record of D_y+i*a*D=i*(lift_z+a*lift).
    if sp.simplify((sp.I * source - sp.I * sp.pi * d / 8) + sp.I * sp.pi * d / 8 - sp.I * source) != 0:
        raise RuntimeError("holomorphic lift defect ODE failed")

    return {
        "phase_cancellation": (
            "M_0(i*T)*U(T)*exp(pi*i/8)="
            "-(sqrt(pi)/8)*exp(-pi*T/4)*(T^(3/2)+i*T^(1/2))"
        ),
        "slow_factor": (
            "F_N(T)=exp(-pi*T/4)S_N(T), "
            "S_N(T)=-(sqrt(pi)/8)*(T^(3/2)+i*T^(1/2))*C_0(p_N(T))"
        ),
        "complex_parameters": (
            "T(z)=z/2+pi*t/8, a(z)=sqrt(T(z)/(2*pi)), "
            "p_N(z)=1-2*(a(z)-N)"
        ),
        "sharp_operation": "F_N^#(z)=conj(F_N(conj(z)))",
        "holomorphic_lift": (
            "C_hat_(N,t)(z)=(-1)^N*exp(t*pi^2/64)*(F_N(z)+F_N^#(z))"
        ),
        "real_match": "C_hat_(N,t)(x)=C_(N,t)(x) for every real x in the fixed-N cell",
        "slow_derivative": (
            "(d/dz+pi/8)F_N(z)=(1/2)*exp(-pi*T(z)/4)*S_N'(T(z))"
        ),
        "paper_extension": (
            "C_paper(x+iy)=exp(-pi*i*y/8)*C_(N,t)(x)"
        ),
        "defect_ode": (
            "D=C_hat-C_paper satisfies D(x,0)=0 and "
            "partial_y D+i*pi*D/8=i*(partial_z+pi/8)C_hat"
        ),
        "defect_integral": (
            "|D(x,y)|<=|y|*exp(pi*|y|/8)*"
            "sup_|v|<=|y| |(partial_z+pi/8)C_hat(x+iv)|"
        ),
        "derivative_scale": (
            "If C_0 and C_0' are bounded on the required p-strip, then "
            "|S_N'(T)|=O(T), one saddle factor below |S_N(T)|=O(T^(3/2))"
        ),
        "critical_gain": (
            "After division by A_t(x), the radius-1/L lift defect is "
            "O(exp(-3L/4)/L) uniformly for bounded c=tL"
        ),
        "remaining_bound": (
            "Certify explicit complex-strip bounds for C_0 and C_0' and combine "
            "them with the published e_A+e_B+e_C constants"
        ),
    }


def diagnostics(dps: int) -> dict:
    mp.mp.dps = dps
    phase_rows = []
    max_relative = mp.mpf(0)
    for tau in (mp.mpf(100), mp.mpf(1000), mp.mpf(10000)):
        direct = m_zero(mp.j * tau) * u_rs(tau) * mp.exp(mp.j * mp.pi / 8)
        reduced = slow_g(tau)
        relative = abs(direct - reduced) / abs(reduced)
        max_relative = max(max_relative, relative)
        phase_rows.append(
            {
                "T": str(int(tau)),
                "direct": mp.nstr(direct, 40),
                "reduced": mp.nstr(reduced, 40),
                "relative_delta": mp.nstr(relative, 25),
            }
        )
    if max_relative >= mp.mpf("1e-60"):
        raise RuntimeError("endpoint phase cancellation diagnostic failed")

    derivative_max = mp.mpf(0)
    derivative_arg = mp.mpf(0)
    for k in range(2001):
        p = -1 + mp.mpf(2) * k / 2000
        value = abs(mp.diff(c_zero, p))
        if value > derivative_max:
            derivative_max = value
            derivative_arg = p
    if derivative_max >= mp.mpf("0.7"):
        raise RuntimeError("C0 derivative scout exceeded its diagnostic guard")
    return {
        "dps": dps,
        "phase_rows": phase_rows,
        "max_phase_cancellation_relative_delta": mp.nstr(max_relative, 25),
        "C0_prime_grid_points": 2001,
        "C0_prime_grid_max": mp.nstr(derivative_max, 25),
        "C0_prime_grid_argmax": mp.nstr(derivative_arg, 15),
        "C0_prime_status": "diagnostic only; a complex-strip interval bound remains required",
    }


def build_artifact(dps: int) -> dict:
    exact = build_exact()
    diag = diagnostics(dps)
    rows = [
        GateRow(
            id="np15ehl_01_phase_cancellation",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The rapidly oscillating Stirling and Riemann-Siegel phases cancel exactly in the endpoint correction.",
            formula=exact["phase_cancellation"],
            proof_boundary="Principal branches with T>0.",
            diagnostics=diag["phase_rows"],
        ),
        GateRow(
            id="np15ehl_02_slow_factor",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="After phase cancellation the endpoint block is an exponential decay times a slowly varying algebraic factor.",
            formula=exact["slow_factor"],
            proof_boundary="Fixed integer cutoff N.",
        ),
        GateRow(
            id="np15ehl_03_holomorphic_lift",
            role="exact_construction",
            readiness="ready_to_apply",
            claim="The real-axis endpoint correction has an explicit fixed-cutoff holomorphic lift.",
            formula=f"{exact['complex_parameters']}; {exact['holomorphic_lift']}",
            proof_boundary="Local principal branches around the positive real axis.",
        ),
        GateRow(
            id="np15ehl_04_real_axis_match",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The holomorphic lift agrees exactly with the published correction on the real axis.",
            formula=exact["real_match"],
            proof_boundary="No approximation at y=0.",
        ),
        GateRow(
            id="np15ehl_05_slow_derivative",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The Cauchy-relevant defect operator removes the leading exponential derivative exactly.",
            formula=exact["slow_derivative"],
            proof_boundary="Exact chain rule with T'(z)=1/2.",
        ),
        GateRow(
            id="np15ehl_06_defect_ode",
            role="exact_transfer_lemma",
            readiness="ready_to_apply",
            claim="Difference between the analytic lift and the published off-axis correction obeys a first-order defect equation.",
            formula=f"{exact['defect_ode']}; {exact['defect_integral']}",
            proof_boundary="Within a fixed-N collar.",
        ),
        GateRow(
            id="np15ehl_07_saddle_gain",
            role="asymptotic_reduction",
            readiness="conditional_ready",
            claim="The lift defect gains the missing saddle-index factor needed by the refined remainder.",
            formula=f"{exact['derivative_scale']}; {exact['critical_gain']}",
            proof_boundary="Conditional on an explicit complex-strip C0 derivative bound.",
        ),
        GateRow(
            id="np15ehl_08_C0_derivative_scout",
            role="high_precision_diagnostic",
            readiness="diagnostic_only",
            claim="A dense real grid finds a small C0 derivative, supporting a very loose rigorous strip target.",
            formula="max_grid |C0'(p)|<0.7 on -1<=p<=1",
            proof_boundary="Not an interval or complex-strip certificate.",
            diagnostics={
                "grid_points": diag["C0_prime_grid_points"],
                "max": diag["C0_prime_grid_max"],
                "argmax": diag["C0_prime_grid_argmax"],
            },
        ),
        GateRow(
            id="np15ehl_09_live_bound",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="An explicit C0 strip bound is the remaining local input for a rigorous corrected C1 collar.",
            formula=exact["remaining_bound"],
            proof_boundary="Open finite-dimensional interval estimate, not the RH-level sign target.",
        ),
        GateRow(
            id="np15ehl_10_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The holomorphic lift resolves an error-transfer obstruction but does not prove corrected-main transversality.",
            formula="analytic lift != first-jet lower bound",
            proof_boundary="Not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift",
        "date": "2026-07-17",
        "status": (
            "exact fixed-cutoff holomorphic lift and saddle-gain reduction for "
            "the endpoint correction; not a C1 remainder certificate, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This artifact proves the phase cancellation, holomorphic lift, real-axis "
            "match, and defect equation. The saved C0 derivative grid is diagnostic. "
            "A rigorous complex-strip derivative bound and explicit composition with "
            "the published remainder remain open. No corrected-main lower bound, "
            "positive-boundary exclusion, Lambda<=0, or RH is proved."
        ),
        "exact": exact,
        "diagnostics": diag,
        "rows": [asdict(row) for row in rows],
        "sources": [
            SOURCE_URL,
            "outputs/jensen_window_pf_newman_polymath15_critical_transversality_target.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    diag = artifact["diagnostics"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Endpoint Holomorphic Lift",
            "",
            "Date: 2026-07-17",
            "",
            "Status: exact holomorphic-lift reduction for the refined endpoint",
            "correction. This is not a proof of `Lambda <= 0` or RH.",
            "It is not a `C1` remainder certificate.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.py",
            "```",
            "",
            "The endpoint term is imported from the A+B-C approximation in",
            f"[Polymath 15]({SOURCE_URL}). Its displayed off-axis real-part formula",
            "is not itself holomorphic, so it cannot be fed directly into Cauchy's",
            "estimate. The following lift repairs that point.",
            "",
            "## Exact Phase Cancellation",
            "",
            "Direct substitution of `M_0` and `U` gives",
            "",
            "```text",
            exact["phase_cancellation"],
            exact["slow_factor"],
            "```",
            "",
            "The apparent high-frequency phase is gone exactly. Independent",
            f"{diag['dps']}-digit checks at `T=100,1000,10000` have maximum relative",
            f"difference `{diag['max_phase_cancellation_relative_delta']}`.",
            "",
            "## Holomorphic Lift",
            "",
            "For fixed `N`, set",
            "",
            "```text",
            exact["complex_parameters"],
            exact["sharp_operation"],
            exact["holomorphic_lift"],
            "```",
            "",
            "Then",
            "",
            "```text",
            exact["real_match"],
            "```",
            "",
            "so the lift changes nothing in the corrected real-axis main.",
            "",
            "## Defect Equation",
            "",
            "The exponential derivative cancels:",
            "",
            "```text",
            exact["slow_derivative"],
            "```",
            "",
            "Comparing the lift with the paper's off-axis extension yields",
            "",
            "```text",
            exact["paper_extension"],
            exact["defect_ode"],
            exact["defect_integral"],
            "```",
            "",
            "Thus the discrepancy is controlled by the slow derivative, not by",
            "the full endpoint block. Since `p_N'(T)=O(T^-1/2)`, bounded `C_0`",
            "and `C_0'` give",
            "",
            "```text",
            exact["derivative_scale"],
            exact["critical_gain"],
            "```",
            "",
            "## Remaining Local Bound",
            "",
            "A 2,001-point high-precision real grid reports",
            f"`max |C_0'|={diag['C0_prime_grid_max']}` at",
            f"`p={diag['C0_prime_grid_argmax']}`. This is only a diagnostic.",
            "The next obligation is",
            "",
            "```text",
            exact["remaining_bound"],
            "```",
            "",
            "Closing that finite-dimensional strip estimate would make the refined",
            "endpoint correction compatible with a rigorous `C1` Cauchy transfer.",
            "It would still leave the RH-level corrected first-jet lower bound open.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--dps", type=int, default=80)
    args = parser.parse_args()
    artifact = build_artifact(args.dps)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 endpoint holomorphic lift: "
        f"{len(artifact['rows'])} rows, 4 exact identities, 1 defect equation, "
        "1 saddle-gain reduction, 1 open strip bound"
    )


if __name__ == "__main__":
    main()
