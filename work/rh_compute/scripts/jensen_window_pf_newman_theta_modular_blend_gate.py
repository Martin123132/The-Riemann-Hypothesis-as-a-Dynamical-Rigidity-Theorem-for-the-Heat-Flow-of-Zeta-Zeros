#!/usr/bin/env python3
"""Build a positive-time-compatible modular theta-block decomposition."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from math import comb, factorial, pi
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import erf
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_theta_modular_blend_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_newman_theta_modular_blend_gate.md"
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


def quadrature(node_count: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = leggauss(node_count)
    panels = (
        (0.0, 0.1),
        (0.1, 0.2),
        (0.2, 0.35),
        (0.35, 0.55),
        (0.55, 0.8),
        (0.8, 1.2),
        (1.2, 2.0),
        (2.0, 4.0),
    )
    u_rows: list[np.ndarray] = []
    weight_rows: list[np.ndarray] = []
    for left, right in panels:
        u_rows.append((right - left) * nodes / 2 + (right + left) / 2)
        weight_rows.append((right - left) * weights / 2)
    return np.concatenate(u_rows), np.concatenate(weight_rows)


def theta_summand(u: np.ndarray, n: int) -> np.ndarray:
    e4u = np.exp(4 * u)
    return (
        2 * pi**2 * n**4 * np.exp(9 * u)
        - 3 * pi * n**2 * np.exp(5 * u)
    ) * np.exp(-pi * n**2 * e4u)


def blend_weight(u: np.ndarray) -> np.ndarray:
    return (1 + erf(3 * np.sinh(4 * u))) / 2


def blended_block(u: np.ndarray, n: int) -> np.ndarray:
    omega = blend_weight(u)
    return omega * theta_summand(u, n) + (1 - omega) * theta_summand(-u, n)


def transform_witnesses(node_count: int) -> list[dict]:
    u, weights = quadrature(node_count)
    rows: list[dict] = []
    for t in (0.0, 0.5):
        heat = np.exp(t * u * u)
        for n, x in ((1, 35.0), (2, 25.0), (3, 15.0)):
            kernel = blended_block(u, n)
            value = float(np.dot(weights, heat * kernel * np.cos(x * u)))
            rows.append(
                {
                    "t": format(t, ".17g"),
                    "n": n,
                    "x": format(x, ".17g"),
                    "block_transform": format(value, ".17e"),
                }
            )
    return rows


def jensen_witnesses(node_count: int) -> list[dict]:
    u, weights = quadrature(node_count)
    kernel = blended_block(u, 1)
    degree = 9
    rows: list[dict] = []
    for t in (0.0, 0.5):
        weighted_kernel = weights * np.exp(t * u * u) * kernel
        moments = np.array(
            [np.dot(weighted_kernel, u ** (2 * k)) for k in range(degree + 1)]
        )
        gamma = np.array(
            [factorial(k) * moments[k] / factorial(2 * k) for k in range(degree + 1)]
        )
        coefficients = np.array(
            [comb(degree, j) * gamma[j] for j in range(degree + 1)]
        )
        variable_scale = (coefficients[0] / coefficients[-1]) ** (1 / degree)
        scaled = coefficients * np.array(
            [variable_scale**j for j in range(degree + 1)]
        )
        scaled /= np.max(np.abs(scaled))
        roots = np.roots(scaled[::-1])
        nonreal = [root for root in roots if root.imag > 1e-8]
        if len(nonreal) != 1:
            raise RuntimeError(f"expected one upper-half-plane Jensen root, got {roots}")
        root = nonreal[0]
        rows.append(
            {
                "t": format(t, ".17g"),
                "block_n": 1,
                "jensen_degree": degree,
                "jensen_shift": 0,
                "positive_variable_scale": format(variable_scale, ".17e"),
                "scaled_nonreal_root_real": format(root.real, ".17e"),
                "scaled_nonreal_root_imag": format(root.imag, ".17e"),
            }
        )
    return rows


def build_numerics() -> dict:
    coarse = transform_witnesses(220)
    fine = transform_witnesses(420)
    max_relative_delta = 0.0
    for left, right in zip(coarse, fine, strict=True):
        left_value = float(left["block_transform"])
        right_value = float(right["block_transform"])
        if right_value >= -1e-12:
            raise RuntimeError(f"negative block-transform witness failed: {right}")
        max_relative_delta = max(
            max_relative_delta,
            abs(left_value - right_value) / max(abs(right_value), 1e-300),
        )
    if max_relative_delta >= 1e-8:
        raise RuntimeError("modular-blend quadrature convergence failed")

    coarse_jensen = jensen_witnesses(220)
    fine_jensen = jensen_witnesses(420)
    max_jensen_root_delta = 0.0
    for left, right in zip(coarse_jensen, fine_jensen, strict=True):
        left_root = complex(
            float(left["scaled_nonreal_root_real"]),
            float(left["scaled_nonreal_root_imag"]),
        )
        right_root = complex(
            float(right["scaled_nonreal_root_real"]),
            float(right["scaled_nonreal_root_imag"]),
        )
        if right_root.imag <= 1e-3:
            raise RuntimeError(f"degree-nine Jensen witness failed: {right}")
        max_jensen_root_delta = max(
            max_jensen_root_delta, abs(left_root - right_root) / abs(right_root)
        )
    if max_jensen_root_delta >= 1e-8:
        raise RuntimeError("modular-blend Jensen-root convergence failed")
    return {
        "method": (
            "Composite Gauss-Legendre quadrature on [0,4] with analytic "
            "Newman weight and the exact modular blend"
        ),
        "coarse_nodes_per_panel": 220,
        "fine_nodes_per_panel": 420,
        "max_relative_coarse_fine_delta": format(max_relative_delta, ".17e"),
        "transform_rows": fine,
        "jensen_rows": fine_jensen,
        "max_relative_jensen_root_delta": format(
            max_jensen_root_delta, ".17e"
        ),
        "jensen_definition": (
            "For G_t(z)=B_(1,t)(i*sqrt(z))=sum_(k>=0)gamma_k*z^k/k!, "
            "gamma_k=k!*mu_(2k)/(2k)!. The degree-nine shift-zero Jensen "
            "polynomial is sum_(j=0)^9 binom(9,j)*gamma_j*X^j."
        ),
        "scope": (
            "Numerical non-promotion witnesses only: positivity of the blended "
            "kernel blocks does not make their cosine transforms nonnegative or "
            "put the individual transforms in the Laguerre-Polya class."
        ),
    }


def build_exact() -> dict:
    u = sp.symbols("u", nonnegative=True)
    slack = (
        10 * u
        - sp.Rational(3, 2) * (sp.exp(8 * u) - 1)
        + 9 * sp.sinh(4 * u) ** 2
    )
    expected_slack = (
        10 * u
        + sp.Rational(3, 4) * sp.exp(8 * u)
        + sp.Rational(9, 4) * sp.exp(-8 * u)
        - 3
    )
    if sp.simplify(sp.expand(slack.rewrite(sp.exp)) - expected_slack) != 0:
        raise RuntimeError("blend positivity slack identity failed")

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
        raise RuntimeError("positive-time coupled Laguerre identity failed")

    f, fp, fpp, r, rp, rpp = sp.symbols(
        "f fp fpp r rp rpp", real=True
    )
    full_laguerre = (fp + rp) ** 2 - (f + r) * (fpp + rpp)
    partial_laguerre = fp**2 - f * fpp
    tail_difference = 2 * fp * rp + rp**2 - f * rpp - r * fpp - r * rpp
    if sp.expand(full_laguerre - partial_laguerre - tail_difference) != 0:
        raise RuntimeError("finite-to-infinite Laguerre error identity failed")

    return {
        "theta_summand": {
            "definition": (
                "phi_n(u)=(2*pi^2*n^4*exp(9u)-3*pi*n^2*exp(5u))*"
                "exp(-pi*n^2*exp(4u))"
            ),
            "kernel_series": "Phi(u)=sum_(n>=1)phi_n(u)=Phi(-u)",
            "scope": "The second equality is the theta modular identity.",
        },
        "modular_switch": {
            "omega": "omega(u)=(1+erf(3*sinh(4u)))/2",
            "reflection": "omega(-u)=1-omega(u)",
            "entire": "omega is entire because erf and sinh are entire",
            "tail_bound": (
                "(1-omega(u))/omega(u)<=erfc(3*sinh(4u))<="
                "exp(-9*sinh(4u)^2), u>=0"
            ),
            "erfc_proof": (
                "exp(-x^2)-erfc(x)>=0 for x>=0: it starts and ends at zero, "
                "and its derivative changes sign once at x=1/sqrt(pi)."
            ),
        },
        "blended_blocks": {
            "definition": (
                "b_n(u)=omega(u)*phi_n(u)+(1-omega(u))*phi_n(-u)"
            ),
            "evenness": "b_n(-u)=b_n(u)",
            "smoothness": "b_n is an even entire function",
        },
        "strict_positivity": {
            "positive_side": "phi_n(u)>0 for n>=1 and u>=0",
            "negative_case": (
                "If phi_n(-u)<0 and q=pi*n^2, then q*exp(-4u)<3/2."
            ),
            "ratio": (
                "phi_n(u)/(-phi_n(-u))=exp(10u)*"
                "(2q*exp(4u)-3)/(3-2q*exp(-4u))*"
                "exp(-q*(exp(4u)-exp(-4u)))"
            ),
            "ratio_lower_bound": (
                "log(phi_n(u)/(-phi_n(-u)))>"
                "10u-(3/2)*(exp(8u)-1)>-9*sinh(4u)^2"
            ),
            "slack_identity": (
                "10u-(3/2)*(exp(8u)-1)+9*sinh(4u)^2="
                "10u+(3/4)*exp(8u)+(9/4)*exp(-8u)-3>0 in the negative case"
            ),
            "negative_case_threshold": (
                "The negative case forces exp(4u)>2*pi/3>2 and hence "
                "exp(8u)>4; the displayed slack is then strictly positive."
            ),
            "conclusion": "b_n(u)>0 for every n>=1 and every real u",
        },
        "exact_partition": {
            "identity": "sum_(n>=1)b_n(u)=Phi(u)",
            "proof": (
                "The sum equals omega(u)*Phi(u)+(1-omega(u))*Phi(-u), "
                "then theta modular evenness applies."
            ),
        },
        "positive_time_compatibility": {
            "block_tail": (
                "For u->+infinity, omega*phi_n(u) is O_n(exp(9u-pi*n^2*exp(4u))) "
                "and (1-omega)*phi_n(-u) is O_n(exp(-5u-9*sinh(4u)^2))."
            ),
            "schwartz": (
                "For every finite T, exp(t*u^2)*b_n(u) is Schwartz uniformly "
                "for 0<=t<=T."
            ),
            "normal_moments": (
                "sum_n integral_R |u|^k*exp(t*u^2)*b_n(u)du="
                "integral_R |u|^k*exp(t*u^2)*Phi(u)du<infinity"
            ),
        },
        "transform_series": {
            "block_transform": (
                "B_(n,t)(x)=integral_0^infinity exp(t*u^2)*b_n(u)*cos(xu)du"
            ),
            "series": "H_t^(k)(x)=sum_(n>=1)B_(n,t)^(k)(x), k=0,1,2,...",
            "convergence": (
                "For every fixed derivative order and finite T, the series is "
                "absolutely and uniformly convergent on R_x times [0,T]."
            ),
        },
        "coupled_laguerre": {
            "entry": (
                "M_(m,n,t)=B_m'*B_n'-(B_m*B_n''+B_n*B_m'')/2"
            ),
            "identity": "L_t(x)=sum_(m,n>=1)M_(m,n,t)(x)",
            "convergence": (
                "The double series is absolutely and uniformly convergent for "
                "x in R and t in every compact positive-time interval."
            ),
        },
        "tail_enclosure": {
            "partial_sum": (
                "S_(N,t)(x)=sum_(1<=n<=N)B_(n,t)(x), "
                "R_(N,t)(x)=H_t(x)-S_(N,t)(x)"
            ),
            "tail_moments": (
                "mu_(N,j)(T)=integral_0^infinity u^j*exp(T*u^2)*"
                "sum_(n>N)b_n(u)du, j=0,1,2"
            ),
            "jet_bound": (
                "|R_(N,t)^(j)(x)|<=mu_(N,j)(T) for every real x and "
                "0<=t<=T"
            ),
            "laguerre_error": (
                "|L[H_t]-L[S_(N,t)]|<=2*|S'|*mu_1+|S|*mu_2+"
                "|S''|*mu_0+mu_1^2+mu_0*mu_2"
            ),
            "criterion": (
                "If L[S_(N,t)] exceeds the displayed error, then L_t(x)>0."
            ),
            "scope": (
                "The bound is uniform in frequency and uses only positivity of "
                "the omitted blended blocks; numerical evaluation of the moments "
                "must still be enclosed rigorously before promotion."
            ),
        },
        "moment_floor_guard": {
            "positive_floor": (
                "For every fixed finite N, mu_(N,j)(T)>0 for j=0,1,2 because "
                "the omitted blended tail is strictly positive."
            ),
            "riemann_lebesgue": (
                "S_(N,t)(x), S_(N,t)'(x), and S_(N,t)''(x) tend to zero as "
                "|x| tends to infinity."
            ),
            "consequence": (
                "The moment-only error majorant tends to the positive constant "
                "mu_1^2+mu_0*mu_2 while L[S_(N,t)] tends to zero."
            ),
            "decision": (
                "For every fixed N, the three-moment criterion is necessarily a "
                "compact-frequency certificate and cannot close the global tail."
            ),
        },
        "decaying_tail_enclosure": {
            "tail_kernel": "r_N(u)=sum_(n>N)b_n(u)=Phi(u)-sum_(n<=N)b_n(u)",
            "derivative_budget": (
                "d_(N,j,m)(T)=(1/2)*sup_(0<=t<=T)||partial_u^m[(i*u)^j*"
                "exp(t*u^2)*r_N(u)]||_L1(R)"
            ),
            "fourier_bound": (
                "|R_(N,t)^(j)(x)|<=d_(N,j,m)(T)/|x|^m for x!=0, "
                "j=0,1,2"
            ),
            "combined_budget": (
                "epsilon_(N,j,m)(x,T)=min(mu_(N,j)(T),"
                "d_(N,j,m)(T)/|x|^m)"
            ),
            "laguerre_use": (
                "Replace each mu_j in the finite-to-infinite Laguerre error by "
                "epsilon_j; the same algebraic error inequality remains valid."
            ),
            "finiteness": (
                "Every d_(N,j,m)(T) is finite because the blended tail and all "
                "of its derivatives remain Schwartz after bounded Newman time."
            ),
            "scope": (
                "This repairs the positive moment floor. Rigorous useful constants "
                "and a lower bound for the finite partial Laguerre expression remain open."
            ),
        },
        "boundary": {
            "decision": (
                "The decomposition removes both the cusp obstruction of reflected "
                "finite theta truncations and the exponential-tail obstruction of "
                "Euler cell blocks. It does not make the block transforms or their "
                "mixed Laguerre entries nonnegative."
            ),
            "live_handoff": (
                "Enclose derivative L1 budgets for the modular tail and compare "
                "their frequency-decaying Laguerre error with a finite blended "
                "partial sum. Termwise Fourier positivity and individual block "
                "Laguerre-Polya membership are rejected."
            ),
        },
    }


def build_artifact() -> dict:
    exact = build_exact()
    numerics = build_numerics()
    rows = [
        GateRow(
            id="ntmbg_01_switch_reflection",
            role="exact_identity",
            readiness="available_exact",
            claim="The entire super-exponential switch exchanges with its complement under reflection.",
            formula=exact["modular_switch"]["reflection"],
            proof_boundary="Exact scalar identity.",
        ),
        GateRow(
            id="ntmbg_02_even_smooth_blocks",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Every modularly blended theta block is even and entire.",
            formula=exact["blended_blocks"]["evenness"],
            proof_boundary="Block regularity only; no spectral sign follows.",
            diagnostics=exact["blended_blocks"],
        ),
        GateRow(
            id="ntmbg_03_strict_block_positivity",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Every modularly blended theta block is strictly positive on the real axis.",
            formula=exact["strict_positivity"]["conclusion"],
            proof_boundary="Kernel positivity only; not Fourier positivity.",
            diagnostics=exact["strict_positivity"],
        ),
        GateRow(
            id="ntmbg_04_exact_partition",
            role="exact_identity",
            readiness="available_exact",
            claim="The positive blended blocks sum exactly to the Xi Fourier kernel.",
            formula=exact["exact_partition"]["identity"],
            proof_boundary="Uses theta modular evenness, not zeta-zero information.",
            diagnostics=exact["exact_partition"],
        ),
        GateRow(
            id="ntmbg_05_positive_time_schwartz",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Every block remains rapidly decreasing after every bounded positive Newman deformation.",
            formula=exact["positive_time_compatibility"]["schwartz"],
            proof_boundary="Kernel-space compatibility only.",
            diagnostics=exact["positive_time_compatibility"],
        ),
        GateRow(
            id="ntmbg_06_normal_transform_series",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The positive-time block transforms and every fixed spectral jet form a normal series.",
            formula=exact["transform_series"]["series"],
            proof_boundary="Absolute uniform convergence; no sign of the sum is asserted.",
            diagnostics=exact["transform_series"],
        ),
        GateRow(
            id="ntmbg_07_coupled_laguerre",
            role="exact_identity",
            readiness="available_exact",
            claim="The strict-Laguerre target has an absolutely convergent positive-time block matrix expansion.",
            formula=exact["coupled_laguerre"]["identity"],
            proof_boundary="The mixed matrix entries retain both signs.",
            diagnostics=exact["coupled_laguerre"],
        ),
        GateRow(
            id="ntmbg_08_tail_enclosure",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Three positive tail moments give a uniform finite-to-infinite Laguerre enclosure.",
            formula=exact["tail_enclosure"]["laguerre_error"],
            proof_boundary="Exact conditional enclosure; no finite partial sum is certified globally here.",
            diagnostics=exact["tail_enclosure"],
        ),
        GateRow(
            id="ntmbg_09_moment_floor_guard",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The positive three-moment enclosure is necessarily compact in frequency for every finite block count.",
            formula=exact["moment_floor_guard"]["consequence"],
            proof_boundary="Rejects global promotion of the moment-only enclosure, not the blended decomposition.",
            diagnostics=exact["moment_floor_guard"],
        ),
        GateRow(
            id="ntmbg_10_decaying_tail_enclosure",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Derivative L1 budgets give arbitrarily high algebraic decay for every remainder jet.",
            formula=exact["decaying_tail_enclosure"]["fourier_bound"],
            proof_boundary="Exact tail theorem; no useful global constants or Laguerre lower bound are certified here.",
            diagnostics=exact["decaying_tail_enclosure"],
        ),
        GateRow(
            id="ntmbg_11_termwise_spectral_guard",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim=(
                "Strict positivity of the blended kernel blocks implies neither "
                "nonnegative transforms nor individual Laguerre-Polya membership."
            ),
            formula=(
                "Six block transforms are negative, and the first block's "
                "degree-nine Jensen polynomial has a nonreal pair at two times."
            ),
            proof_boundary=(
                "Numerical witnesses only; rejects termwise transform positivity "
                "and the naive individual Laguerre-Polya promotion."
            ),
            diagnostics=numerics,
        ),
        GateRow(
            id="ntmbg_12_live_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The surviving use of the blend is a coupled spectral coercivity or tail-dominance estimate.",
            formula=exact["boundary"]["live_handoff"],
            proof_boundary="Open; not strict Laguerre positivity, Lambda<=0, or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_theta_modular_blend_gate",
        "date": "2026-07-16",
        "status": (
            "exact positive entire t-compatible modular theta blend with termwise "
            "spectral non-promotion guard"
        ),
        "proof_boundary": (
            "This artifact constructs an exact, strictly positive, entire, even "
            "theta-block decomposition of Phi whose blocks remain Schwartz under "
            "every bounded positive Newman deformation. It proves normal transform "
            "and coupled Laguerre matrix expansions and records explicit negative "
            "block-transform and nonreal Jensen witnesses. It does not prove the coupled matrix is "
            "positive, strict Laguerre positivity, Lambda<=0, or RH."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_theta_summand_spectral_square_gate.md",
            "outputs/jensen_window_pf_newman_theta_cell_renormalization_gate.md",
            "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md",
            "outputs/formal_core.md",
        ],
        "exact": exact,
        "numerics": numerics,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    numerics = artifact["numerics"]
    witness_lines = [
        (
            f"t={row['t']}, n={row['n']}, x={row['x']}, "
            f"B={row['block_transform']}"
        )
        for row in numerics["transform_rows"]
    ]
    jensen_lines = [
        (
            f"t={row['t']}, degree={row['jensen_degree']}, "
            f"scaled root={row['scaled_nonreal_root_real']}"
            f" + {row['scaled_nonreal_root_imag']}*i"
        )
        for row in numerics["jensen_rows"]
    ]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Theta Modular-Blend Gate",
            "",
            "Date: 2026-07-16",
            "",
            "Status: exact positive-time-compatible modular theta blend with a",
            "termwise spectral non-promotion guard. This is not a proof of RH or",
            "`Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_theta_modular_blend_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_theta_modular_blend_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_blend_gate.py",
            "```",
            "",
            "## Modular Blend",
            "",
            "For the standard theta summands, set",
            "",
            "```text",
            exact["theta_summand"]["definition"],
            exact["modular_switch"]["omega"],
            exact["blended_blocks"]["definition"],
            "```",
            "",
            "The switch is entire and satisfies `omega(-u)=1-omega(u)`, so every",
            "`b_n` is even and entire. Theta modular evenness gives the exact partition",
            "",
            "```text",
            exact["exact_partition"]["identity"],
            "```",
            "",
            "## Strict Positivity",
            "",
            "For `u>=0`, `phi_n(u)>0`. If the reflected summand is negative,",
            "write `q=pi*n^2`; then `q*exp(-4u)<3/2` and",
            "",
            "```text",
            exact["strict_positivity"]["ratio"],
            exact["modular_switch"]["tail_bound"],
            exact["strict_positivity"]["ratio_lower_bound"],
            exact["strict_positivity"]["slack_identity"],
            "```",
            "",
            "Thus the positive forward summand dominates the reflected negative",
            "piece by more than the exact blend ratio, proving",
            "",
            "```text",
            exact["strict_positivity"]["conclusion"],
            "```",
            "",
            "## Positive Newman Time",
            "",
            "The unwanted reflected tail is multiplied by a double-exponential",
            "switch:",
            "",
            "```text",
            exact["positive_time_compatibility"]["block_tail"],
            exact["positive_time_compatibility"]["schwartz"],
            "```",
            "",
            "Positivity and Tonelli then give the exact normal-moment identity",
            "",
            "```text",
            exact["positive_time_compatibility"]["normal_moments"],
            "```",
            "",
            "Consequently, for every fixed spectral derivative and bounded time",
            "interval,",
            "",
            "```text",
            exact["transform_series"]["series"],
            exact["coupled_laguerre"]["identity"],
            "```",
            "",
            "Both series are absolutely and uniformly convergent. This removes the",
            "cusp obstruction of hard-reflected finite summands and the exponential",
            "tail obstruction of the smooth Euler-cell blocks.",
            "",
            "## Tail Enclosure",
            "",
            "For a finite blended sum, define",
            "",
            "```text",
            exact["tail_enclosure"]["partial_sum"],
            exact["tail_enclosure"]["tail_moments"],
            "```",
            "",
            "Positivity of every omitted block gives a frequency-uniform jet bound",
            "and hence the exact a posteriori criterion",
            "",
            "```text",
            exact["tail_enclosure"]["jet_bound"],
            exact["tail_enclosure"]["laguerre_error"],
            exact["tail_enclosure"]["criterion"],
            "```",
            "",
            "This turns the decomposition into a usable finite-to-infinite bridge,",
            "but the partial-sum inequality and tail moments still require rigorous",
            "uniform enclosures on any claimed region.",
            "",
            "## Frequency Tail",
            "",
            "The three-moment bound cannot itself be global. For every fixed `N`,",
            "the omitted positive tail makes all three moments strictly positive,",
            "whereas Riemann-Lebesgue gives",
            "",
            "```text",
            exact["moment_floor_guard"]["riemann_lebesgue"],
            exact["moment_floor_guard"]["consequence"],
            exact["moment_floor_guard"]["decision"],
            "```",
            "",
            "The exact repair is to retain oscillatory decay. Since the entire",
            "blended remainder is Schwartz after bounded Newman time, define",
            "",
            "```text",
            exact["decaying_tail_enclosure"]["derivative_budget"],
            exact["decaying_tail_enclosure"]["fourier_bound"],
            exact["decaying_tail_enclosure"]["combined_budget"],
            "```",
            "",
            "Replacing each `mu_j` by `epsilon_j` in the Laguerre error is valid",
            "and removes the positive floor. The remaining hard task is quantitative:",
            "obtain useful rigorous derivative budgets and compare them with a",
            "global lower profile for the finite blended Laguerre expression.",
            "",
            "## Spectral Guard",
            "",
            "Positive smooth blocks still do not have nonnegative cosine transforms.",
            "Independent coarse/fine quadrature gives",
            "",
            "```text",
            *witness_lines,
            "```",
            "",
            "The stronger individual Laguerre-Polya shortcut also fails numerically.",
            "For the first block, form `G_t(z)=B_(1,t)(i*sqrt(z))`; its degree-nine",
            "shift-zero Jensen polynomial has the stable nonreal roots",
            "",
            "```text",
            *jensen_lines,
            "```",
            "",
            "So termwise Fourier positivity and individual block Laguerre-Polya",
            "promotion are closed. The surviving use of this",
            "decomposition is a coupled matrix estimate, a relative coercivity",
            "bound, or a transform-tail dominance theorem that retains the infinite",
            "arithmetic cancellation. None of those signs is proved here.",
            "",
        ]
    )


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
        "built Jensen-window PF Newman theta modular-blend gate: "
        f"{len(artifact['rows'])} rows, "
        f"{len(artifact['numerics']['transform_rows'])} transform witnesses, "
        f"{len(artifact['numerics']['jensen_rows'])} Jensen witnesses"
    )


if __name__ == "__main__":
    main()
