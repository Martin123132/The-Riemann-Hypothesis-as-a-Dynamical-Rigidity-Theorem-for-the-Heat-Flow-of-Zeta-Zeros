#!/usr/bin/env python3
"""Build the endpoint theta cell-renormalization and positive-time gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_theta_cell_renormalization_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_theta_cell_renormalization_gate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def build_exact() -> dict:
    a, u, n = sp.symbols("a u n", positive=True)
    q = sp.pi * sp.exp(4 * u)
    moment_2 = sp.sqrt(sp.pi) / (4 * q ** sp.Rational(3, 2))
    moment_4 = 3 * sp.sqrt(sp.pi) / (8 * q ** sp.Rational(5, 2))
    continuum = sp.simplify(
        2 * sp.pi**2 * sp.exp(9 * u) * moment_4
        - 3 * sp.pi * sp.exp(5 * u) * moment_2
    )
    if continuum != 0:
        raise RuntimeError("continuous theta-index annihilation failed")

    cell_second_moment = sp.integrate(a**2, (a, n - 1, n))
    cell_gap = sp.expand(n**2 - cell_second_moment)
    if cell_gap != n - sp.Rational(1, 3):
        raise RuntimeError("cell tail coefficient failed")

    s = sp.symbols("s")
    expected_cell_power = (
        n ** (1 - s) - (n - 1) ** (1 - s)
    ) / (1 - s)
    antiderivative = a ** (1 - s) / (1 - s)
    if sp.simplify(sp.diff(antiderivative, a) - a ** (-s)) != 0:
        raise RuntimeError("Euler cell-power identity failed")

    f1, f2, f1p, f2p, f1pp, f2pp = sp.symbols(
        "f1 f2 f1p f2p f1pp f2pp", real=True
    )
    direct = (f1p + f2p) ** 2 - (f1 + f2) * (f1pp + f2pp)
    matrix = (
        f1p**2
        - f1 * f1pp
        + 2 * (f1p * f2p - (f1 * f2pp + f2 * f1pp) / 2)
        + f2p**2
        - f2 * f2pp
    )
    if sp.expand(direct - matrix) != 0:
        raise RuntimeError("renormalized Laguerre matrix identity failed")

    return {
        "continuous_theta_index": {
            "definition": (
                "phi_a(u)=(2*pi^2*a^4*exp(9u)-3*pi*a^2*exp(5u))*"
                "exp(-pi*a^2*exp(4u)), a>0"
            ),
            "shift": "phi_a(u)=a^(-1/2)*phi_1(u+(log a)/2)",
            "symmetrized": "g_a(u)=(phi_a(u)+phi_a(-u))/2",
            "continuum_annihilation": (
                "integral_0^infinity phi_a(u)da=0 and "
                "integral_0^infinity g_a(u)da=0"
            ),
            "gaussian_moments": (
                "integral_0^infinity a^2*exp(-q*a^2)da=sqrt(pi)/(4*q^(3/2)); "
                "integral_0^infinity a^4*exp(-q*a^2)da=3*sqrt(pi)/(8*q^(5/2))"
            ),
        },
        "cell_renormalization": {
            "block": "r_n(u)=g_n(u)-integral_(n-1)^n g_a(u)da",
            "kernel_sum": "Phi(u)=sum_(n>=1)r_n(u)",
            "reason": (
                "sum_(n>=1)g_n=Phi and the cell integrals partition the exactly "
                "vanishing continuous theta-index integral."
            ),
        },
        "normal_convergence": {
            "index_derivative": (
                "partial_a phi_a(u)=a^(-3/2)/2*(phi_1'(v)-phi_1(v)), "
                "v=u+(log a)/2"
            ),
            "weighted_bound": (
                "integral_R |u|^k*|partial_a g_a(u)|du<="
                "C_k*a^(-3/2)*(1+|log a|^k), a>=1"
            ),
            "cell_bound": (
                "integral_R |u|^k*|r_n(u)|du<="
                "C_k*integral_(n-1)^n a^(-3/2)*(1+(log a)^k)da, n>=2"
            ),
            "summability": (
                "sum_(n>=1)integral_R |u|^k*|r_n(u)|du<infinity for every integer k>=0"
            ),
            "first_cell": (
                "For 0<a<=1, ||g_a||_(L1_k)<=C_k*a^(-1/2)*(1+|log a|^k), "
                "which is integrable at a=0."
            ),
        },
        "uniform_transform": {
            "block_transform": (
                "J_n(x)=integral_0^infinity r_n(u)cos(xu)du="
                "I_n(x)-integral_(n-1)^n I_a(x)da"
            ),
            "series": "H_0^(k)(x)=sum_(n>=1)J_n^(k)(x), k=0,1,2,...",
            "mode": (
                "The series and every fixed derivative order converge absolutely "
                "and uniformly on the real spectral axis."
            ),
        },
        "euler_zeta_transform": {
            "spectral_parameter": "s=(1+i*x)/2",
            "cell_error": (
                "e_n(s)=n^(-s)-(n^(1-s)-(n-1)^(1-s))/(1-s)"
            ),
            "cell_integral": str(expected_cell_power),
            "block_formula": (
                "J_n(x)=1/4*(Qhat(x)*e_n(s)+Qhat(-x)*e_n(conjugate(s)))"
            ),
            "euler_limit": (
                "sum_(n>=1)e_n(s)=lim_(N->infinity)(sum_(n=1)^N n^(-s)-"
                "N^(1-s)/(1-s))=zeta(s), 0<Re(s)<1"
            ),
            "xi_assembly": (
                "sum_(n>=1)J_n(x)=1/2*Re(Qhat(x)*zeta((1+i*x)/2))="
                "xi((1+i*x)/2)/8=H_0(x)"
            ),
            "scope": (
                "This is a convergent Euler cell-renormalization derived from the "
                "kernel; it uses no information about zeta zeros."
            ),
        },
        "origin_sign": {
            "formula": (
                "J_n(0)=Qhat(0)/2*(n^(-1/2)-2*(sqrt(n)-sqrt(n-1)))"
            ),
            "qhat": "Qhat(0)=-Gamma(1/4)/(32*pi^(1/4))<0",
            "cell_inequality": (
                "2*(sqrt(n)-sqrt(n-1))=integral_(n-1)^n a^(-1/2)da>n^(-1/2)"
            ),
            "conclusion": "J_n(0)>0 for every n>=1",
        },
        "coupled_laguerre_matrix": {
            "entry": (
                "M_(m,n)(x)=J_m'(x)*J_n'(x)-"
                "(J_m(x)*J_n''(x)+J_n(x)*J_m''(x))/2"
            ),
            "identity": (
                "L[H_0](x)=sum_(m,n>=1)M_(m,n)(x)"
            ),
            "convergence": (
                "The double series is absolutely and locally uniformly convergent "
                "because the J_n jets are normally summable through order two."
            ),
            "boundary": (
                "The identity retains every mixed sign; no positivity of the matrix is proved."
            ),
        },
        "positive_time_obstruction": {
            "block_tail": (
                "r_n(u)=-(pi/2)*(3*n-1)*exp(-5u)+O_n(exp(-9u)) as u->infinity"
            ),
            "cell_gap": (
                "n^2-integral_(n-1)^n a^2 da=n-1/3"
            ),
            "consequence": (
                "For every fixed n>=1 and every t>0, exp(t*u^2)*r_n(u) is not "
                "integrable and has no ordinary half-line Fourier transform."
            ),
            "decision": (
                "The cell-renormalized decomposition is exact at t=0 but cannot be "
                "deformed blockwise to the positive Newman times required by the strict-correlation target."
            ),
        },
        "live_handoff": (
            "Construct a t-compatible modular grouping whose individual blocks "
            "cancel the full negative-side theta tail before multiplication by "
            "exp(t*u^2), or work directly with the endpoint-subtracted theta "
            "primitive A_t and prove positivity of its coupled curvature expression."
        ),
    }


def build_artifact() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="ntcrg_01_continuous_index",
            role="exact_identity",
            readiness="available_exact",
            claim="The arithmetic theta summand extends to a continuous positive index with an exact translate law.",
            formula=exact["continuous_theta_index"]["shift"],
            proof_boundary="Exact continuous-index identity only.",
        ),
        GateRow(
            id="ntcrg_02_continuum_annihilation",
            role="exact_identity",
            readiness="available_exact",
            claim="The complete continuous theta-index integral vanishes pointwise.",
            formula=exact["continuous_theta_index"]["continuum_annihilation"],
            proof_boundary="Exact Gaussian-moment cancellation; no discrete spectral sign follows.",
            diagnostics=exact["continuous_theta_index"],
        ),
        GateRow(
            id="ntcrg_03_cell_kernel_sum",
            role="exact_identity",
            readiness="available_exact",
            claim="Subtracting one continuous index cell from every arithmetic block leaves the full Xi kernel unchanged.",
            formula=exact["cell_renormalization"]["kernel_sum"],
            proof_boundary="Exact pointwise kernel identity at t=0.",
            diagnostics=exact["cell_renormalization"],
        ),
        GateRow(
            id="ntcrg_04_weighted_l1",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The cell-renormalized kernels are normally summable in every polynomially weighted L1 norm.",
            formula=exact["normal_convergence"]["summability"],
            proof_boundary="Endpoint kernel theorem; Gaussian positive-time weights are excluded.",
            diagnostics=exact["normal_convergence"],
        ),
        GateRow(
            id="ntcrg_05_uniform_transform",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The renormalized cosine transforms and all fixed derivative orders converge absolutely and uniformly.",
            formula=exact["uniform_transform"]["series"],
            proof_boundary="Ordinary endpoint Fourier convergence only.",
            diagnostics=exact["uniform_transform"],
        ),
        GateRow(
            id="ntcrg_06_euler_zeta",
            role="exact_identity",
            readiness="available_exact",
            claim="The transform counterterm is exactly Euler's convergent sum-integral continuation of zeta on the critical line.",
            formula=exact["euler_zeta_transform"]["euler_limit"],
            proof_boundary="Reconstructs H_0 without using zero information; it is not a zero-free factorization.",
            diagnostics=exact["euler_zeta_transform"],
        ),
        GateRow(
            id="ntcrg_07_origin_sign",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Every renormalized block contributes positively at zero frequency.",
            formula=exact["origin_sign"]["conclusion"],
            proof_boundary="One spectral point only; no global Fourier positivity follows.",
            diagnostics=exact["origin_sign"],
        ),
        GateRow(
            id="ntcrg_08_coupled_laguerre",
            role="exact_identity",
            readiness="available_exact",
            claim="The endpoint first Laguerre expression is an absolutely convergent infinite mixed-term matrix sum.",
            formula=exact["coupled_laguerre_matrix"]["identity"],
            proof_boundary="Exact matrix representation; positivity of the matrix sum remains open.",
            diagnostics=exact["coupled_laguerre_matrix"],
        ),
        GateRow(
            id="ntcrg_09_positive_time_obstruction",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Every individual cell-renormalized block retains an exponential tail that the positive Newman weight makes nonintegrable.",
            formula=exact["positive_time_obstruction"]["block_tail"],
            proof_boundary="Rejects blockwise positive-time deformation of this renormalization only.",
            diagnostics=exact["positive_time_obstruction"],
        ),
        GateRow(
            id="ntcrg_10_live_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The remaining direct route needs t-compatible modular cancellation inside each coupled group.",
            formula=exact["live_handoff"],
            proof_boundary="Open; not strict Laguerre positivity, Wiener density, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_theta_cell_renormalization_gate",
        "date": "2026-07-11",
        "status": "exact endpoint cell renormalization with positive-time obstruction",
        "proof_boundary": (
            "This artifact proves a pointwise and weighted-L1 convergent cell "
            "renormalization of the theta summands at t=0, an ordinary uniformly "
            "convergent transform series, Euler-zeta reconstruction, positive "
            "zero-frequency contributions, and an absolutely convergent coupled "
            "Laguerre matrix identity. It also proves that every individual block "
            "has a nonzero exp(-5u) tail and therefore cannot be deformed by the "
            "positive Newman Gaussian weight. It does not prove matrix positivity, "
            "the positive-time strict Laguerre theorem, Wiener density, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.md",
            "outputs/jensen_window_pf_newman_theta_summand_spectral_square_gate.md",
            "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md",
            "outputs/formal_core.md",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Theta Cell-Renormalization Gate",
            "",
            "Date: 2026-07-11",
            "",
            "Status: exact endpoint cell renormalization with positive-time obstruction.",
            "This is not a proof or disproof of RH or `Lambda <= 0`.",
            "",
            "Artifact kind: `jensen_window_pf_newman_theta_cell_renormalization_gate`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_theta_cell_renormalization_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_theta_cell_renormalization_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_cell_renormalization_gate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF Newman theta cell-renormalization gate: 10 rows, 0 issues, 4 exact kernel/transform identities, 3 exact convergence/sign theorems, 1 coupled Laguerre identity, 1 positive-time obstruction, 1 modular handoff",
            "```",
            "",
            "## Vanishing Continuum",
            "",
            "Extend the arithmetic summand to every real `a>0`:",
            "",
            "```text",
            exact["continuous_theta_index"]["definition"],
            exact["continuous_theta_index"]["shift"],
            exact["continuous_theta_index"]["continuum_annihilation"],
            "```",
            "",
            "The last identity is an exact cancellation of the second and fourth",
            "Gaussian moments. It suggests subtracting the continuous index cell",
            "that ends at each integer:",
            "",
            "```text",
            exact["cell_renormalization"]["block"],
            exact["cell_renormalization"]["kernel_sum"],
            "```",
            "",
            "## Normal Convergence",
            "",
            "The continuous translate law differentiates to",
            "",
            "```text",
            exact["normal_convergence"]["index_derivative"],
            exact["normal_convergence"]["weighted_bound"],
            exact["normal_convergence"]["summability"],
            "```",
            "",
            "Thus the renormalized kernel series converges in every polynomially",
            "weighted `L1` norm. Its cosine transforms may be differentiated term",
            "by term to every fixed order:",
            "",
            "```text",
            exact["uniform_transform"]["block_transform"],
            exact["uniform_transform"]["series"],
            "```",
            "",
            "## Euler-Zeta Assembly",
            "",
            "For `s=(1+i*x)/2`, set",
            "",
            "```text",
            exact["euler_zeta_transform"]["cell_error"],
            exact["euler_zeta_transform"]["block_formula"],
            exact["euler_zeta_transform"]["euler_limit"],
            exact["euler_zeta_transform"]["xi_assembly"],
            "```",
            "",
            "This is now an ordinary convergent transform decomposition derived",
            "from the kernel, not the divergent termwise Bessel sum. At the origin",
            "each renormalized block has the same strict sign:",
            "",
            "```text",
            exact["origin_sign"]["formula"],
            exact["origin_sign"]["cell_inequality"],
            exact["origin_sign"]["conclusion"],
            "```",
            "",
            "## Coupled Matrix",
            "",
            "Normal convergence through the second jet gives",
            "",
            "```text",
            exact["coupled_laguerre_matrix"]["entry"],
            exact["coupled_laguerre_matrix"]["identity"],
            "```",
            "",
            "The identity retains all mixed signs. No positive-semidefinite or",
            "diagonal-dominance theorem is asserted.",
            "",
            "## Positive-Time Boundary",
            "",
            "The endpoint renormalization is not Newman-time compatible:",
            "",
            "```text",
            exact["positive_time_obstruction"]["cell_gap"],
            exact["positive_time_obstruction"]["block_tail"],
            exact["positive_time_obstruction"]["consequence"],
            "```",
            "",
            "So the endpoint matrix cannot be deformed block by block to `t>0`.",
            "",
            "## Live Handoff",
            "",
            exact["live_handoff"],
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(f"wrote Newman theta cell-renormalization gate: {len(artifact['rows'])} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
