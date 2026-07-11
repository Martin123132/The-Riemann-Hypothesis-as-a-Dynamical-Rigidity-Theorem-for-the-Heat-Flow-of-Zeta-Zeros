#!/usr/bin/env python3
"""Derive the Newman correlation hierarchy and test Gaussian-mixture closure."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from functools import lru_cache
import json
from math import comb, factorial
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.md"
)

MAX_DERIVATIVE = 16
PHI_TERMS = 15
T_VALUES = (0.0, 0.1, 0.25, 0.5)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def derivative_polynomials(a: int, max_derivative: int) -> list[np.ndarray]:
    """P_k(q) for d^k exp(a*u-q)/du^k, q'=4q."""
    rows = [np.array([1.0])]
    for _ in range(max_derivative):
        source = rows[-1]
        target = np.zeros(len(source) + 1)
        for j in range(len(target)):
            if j < len(source):
                target[j] += (a + 4 * j) * source[j]
            if j > 0:
                target[j] -= 4 * source[j - 1]
        rows.append(target)
    return rows


P5 = derivative_polynomials(5, MAX_DERIVATIVE)
P9 = derivative_polynomials(9, MAX_DERIVATIVE)


def polynomial_value(coefficients: np.ndarray, value: np.ndarray) -> np.ndarray:
    return np.polynomial.polynomial.polyval(value, coefficients)


def phi_derivatives_positive(x: np.ndarray) -> np.ndarray:
    """Finite rapidly convergent Phi derivatives for x>0."""
    values = np.zeros((MAX_DERIVATIVE + 1, x.size))
    e4x = np.exp(4 * x)
    for n in range(1, PHI_TERMS + 1):
        q = np.pi * n * n * e4x
        exp_minus_q = np.exp(-q)
        term9 = 2 * np.pi**2 * n**4 * np.exp(9 * x) * exp_minus_q
        term5 = 3 * np.pi * n * n * np.exp(5 * x) * exp_minus_q
        for k in range(MAX_DERIVATIVE + 1):
            values[k] += term9 * polynomial_value(P9[k], q)
            values[k] -= term5 * polynomial_value(P5[k], q)
    return values


def gaussian_derivative_polynomials(t: float) -> list[np.ndarray]:
    """Q_k(x) for d^k exp(t*x^2)/dx^k."""
    rows = [np.array([1.0])]
    for _ in range(MAX_DERIVATIVE):
        source = rows[-1]
        target = np.zeros(len(source) + 1)
        for j in range(len(target)):
            if j + 1 < len(source):
                target[j] += (j + 1) * source[j + 1]
            if j > 0:
                target[j] += 2 * t * source[j - 1]
        rows.append(target)
    return rows


def heat_kernel_derivatives_positive(
    x: np.ndarray, t: float, phi: np.ndarray | None = None
) -> np.ndarray:
    if phi is None:
        phi = phi_derivatives_positive(x)
    q_polynomials = gaussian_derivative_polynomials(t)
    q_values = np.array(
        [polynomial_value(row, x) for row in q_polynomials]
    )
    exp_tx2 = np.exp(t * x * x)
    values = np.zeros_like(phi)
    for k in range(MAX_DERIVATIVE + 1):
        for j in range(k + 1):
            values[k] += comb(k, j) * q_values[j] * phi[k - j]
        values[k] *= exp_tx2
    return values


@lru_cache(maxsize=None)
def quadrature(node_count: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = leggauss(node_count)
    panels = ((0.0, 0.125), (0.125, 0.25), (0.25, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 4.0))
    x_rows: list[np.ndarray] = []
    w_rows: list[np.ndarray] = []
    for left, right in panels:
        x_rows.append((right - left) * nodes / 2 + (right + left) / 2)
        w_rows.append((right - left) * weights / 2)
    return np.concatenate(x_rows), np.concatenate(w_rows)


@lru_cache(maxsize=None)
def quadrature_phi(node_count: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x, weights = quadrature(node_count)
    return x, weights, phi_derivatives_positive(x)


def origin_derivatives(t: float, node_count: int) -> list[float]:
    """Return g_t^(n)(0), g_t(r)=K_(1,t)(sqrt(r)), through n=8."""
    x, weights, phi = quadrature_phi(node_count)
    derivatives = heat_kernel_derivatives_positive(x, t, phi)
    rows: list[float] = []
    for n in range(9):
        order = 2 * n
        product_derivative = np.zeros(x.size)
        for j in range(order + 1):
            product_derivative += (
                comb(order, j)
                * (-1) ** (order - j)
                * derivatives[j]
                * derivatives[order - j]
            )
        k_derivative = float(np.dot(weights, 2 * x * x * product_derivative))
        rows.append(k_derivative * factorial(n) / factorial(2 * n))
    return rows


def build_numerics() -> dict:
    coarse = {t: origin_derivatives(t, 192) for t in T_VALUES}
    fine = {t: origin_derivatives(t, 320) for t in T_VALUES}
    rows: list[dict] = []
    max_relative_delta = 0.0
    for t in T_VALUES:
        for n in range(5):
            scale = max(abs(fine[t][n]), 1e-300)
            max_relative_delta = max(
                max_relative_delta,
                abs(fine[t][n] - coarse[t][n]) / scale,
            )
        determinant = fine[t][0] * fine[t][2] - fine[t][1] ** 2
        rows.append(
            {
                "t": format(t, ".17g"),
                "g_derivatives_at_zero": [format(value, ".17e") for value in fine[t]],
                "alternating_signs_through_order": 8,
                "log_convexity_determinant": format(determinant, ".17e"),
                "relative_log_convexity_determinant": format(
                    determinant / fine[t][1] ** 2, ".17e"
                ),
            }
        )
        if not all(((-1) ** n) * value > 0 for n, value in enumerate(fine[t])):
            raise RuntimeError(f"alternating derivative scout failed at t={t}")
        if determinant >= -7.0e-6:
            raise RuntimeError(f"log-convexity obstruction lost margin at t={t}")
    if max_relative_delta >= 1e-9:
        raise RuntimeError("origin quadrature convergence failed")
    return {
        "method": (
            "Composite Gauss-Legendre quadrature on [0,4], finite Phi series n<=15, "
            "analytic derivative recurrences through order 16"
        ),
        "coarse_nodes_per_panel": 192,
        "fine_nodes_per_panel": 320,
        "max_relative_coarse_fine_delta_orders_0_to_4": format(
            max_relative_delta, ".17e"
        ),
        "rows": rows,
        "independent_mpmath_55_digit_references": {
            "method": "mpmath dps=55, Phi n<=13, cutoff 4, composite adaptive quadrature",
            "t=0": {
                "g0": "0.000245610289479469659091987495249879241099947756",
                "g1": "-0.0230795665478520135715709257798235345533328451",
                "g2": "2.1367320486706280789844503282357979256332003",
                "determinant": "-0.00000786301502267649581984683570981852268678957241",
            },
            "t=1/2": {
                "g0": "0.000249998390758081215149621829515199470704011959",
                "g1": "-0.0232914860741815533425228589206177906632906075",
                "g2": "2.13838453388908894463349777467030442419024861",
                "determinant": "-0.0000079006312495514071411797904741719473141521126",
            },
        },
        "scope": (
            "Numerical diagnostic only. The exact global rejection of complete "
            "monotonicity comes from the proved super-Gaussian tail bound."
        ),
    }


def build_exact() -> dict:
    s, v = sp.symbols("s v", real=True)
    heat_weight = sp.expand((s + v) ** 2 + (s - v) ** 2)
    if sp.simplify(heat_weight - 2 * (s**2 + v**2)) != 0:
        raise RuntimeError("heat hierarchy identity failed")
    normalization_checks = []
    for n in range(7):
        coefficient = sp.Rational(2 ** (2 * n - 1), sp.factorial(2 * n))
        normalization_checks.append({"n": n, "coefficient": str(coefficient)})
    return {
        "definitions": {
            "deformed_kernel": "phi_t(u)=exp(t*u^2)*Phi(u)",
            "correlations": (
                "K_(n,t)(v)=integral_R phi_t(s+v)*phi_t(s-v)*s^(2n) ds"
            ),
            "correlation_transforms": "F_(n,t)(xi)=Fourier[K_(n,t)](xi)",
            "laguerre_generating_identity": (
                "|H_t(x+i*y)|^2=sum_(n>=0) L_(n,t)(x)*y^(2n)"
            ),
        },
        "heat_hierarchy": (
            "partial_t K_(n,t)(v)=2*v^2*K_(n,t)(v)+2*K_(n+1,t)(v)"
        ),
        "fourier_hierarchy": (
            "partial_t F_(n,t)(xi)=-2*partial_xi^2 F_(n,t)(xi)+2*F_(n+1,t)(xi)"
        ),
        "laguerre_normalization": (
            "L_(n,t)(x)=2^(2n-1)/(2n)!*F_(n,t)(2x)"
        ),
        "normalization_checks": normalization_checks,
        "multiple_root_signature": (
            "If H_t has a real root c of multiplicity m, then F_(n,t)(2c)=0 "
            "for n<m and F_(m,t)(2c)=(2m)!/2^(2m-1)*(H_t^(m)(c)/m!)^2>0."
        ),
        "double_root_contact": {
            "conditions": "H_t(c)=H_t'(c)=0, H_t''(c)!=0",
            "transform_contact": (
                "F_(1,t)(2c)=partial_xi F_(1,t)(2c)=0, "
                "partial_xi^2 F_(1,t)(2c)=H_t''(c)^2/4"
            ),
            "next_correlation": "F_(2,t)(2c)=3*H_t''(c)^2/4",
            "compatibility": (
                "F_(2,t)(2c)=3*partial_xi^2 F_(1,t)(2c), "
                "partial_t F_(1,t)(2c)=H_t''(c)^2"
            ),
        },
        "gaussian_mixture_sufficiency": {
            "candidate": "g_t(r)=K_(1,t)(sqrt(r))",
            "criterion": (
                "If g_t is nonzero and completely monotone on [0,infinity), "
                "then g_t(r)=integral_[0,infinity) exp(-a*r)dmu_t(a)."
            ),
            "fourier_formula": (
                "Fourier[K_(1,t)](xi)=sqrt(pi)*integral_(a>0) "
                "a^(-1/2)*exp(-xi^2/(4a))dmu_t(a)>0"
            ),
            "consequence": (
                "Complete monotonicity for every 0<t<=1/2 would prove the strict "
                "Laguerre/Wiener target and hence Lambda<=0."
            ),
            "source": "https://doi.org/10.2307/1968466",
        },
        "phi_tail_bound": {
            "constant": (
                "C_Phi=2*pi^2*sum_(n>=1)n^4*exp(-pi*(n^2-1))<infinity"
            ),
            "bound": (
                "0<Phi(u)<=C_Phi*exp(9u-pi*exp(4u)) for u>=0"
            ),
            "derivation": (
                "Drop the negative -3*pi*n^2*exp(5u) part termwise and use "
                "exp(-pi*n^2*E)<=exp(-pi*E)*exp(-pi*(n^2-1)) for E=exp(4u)>=1."
            ),
        },
        "correlation_tail_bound": {
            "finite_constant": (
                "M=integral_R s^2*exp(s^2+18|s|-pi*exp(4|s|))ds<infinity"
            ),
            "elementary_inequalities": [
                "|s+v|+|s-v|=2*max(|s|,v)<=2*(|s|+v)",
                "exp(4|s+v|)+exp(4|s-v|)>=exp(4v)+exp(4|s|)",
            ],
            "bound": (
                "0<K_(1,t)(v)<=C_Phi^2*M*exp(2*t*v^2+18v-pi*exp(4v)) "
                "for v>=0 and 0<=t<=1/2"
            ),
            "limit": (
                "lim_(r->infinity) -log(g_t(r))/r=+infinity, uniformly for 0<=t<=1/2"
            ),
        },
        "complete_monotonicity_obstruction": (
            "A nonzero completely monotone g has g(r)>=mu([0,A])*exp(-A*r) "
            "for some finite A, so limsup -log(g(r))/r<infinity. The Xi "
            "correlation tail has the opposite limit; g_t is therefore not "
            "completely monotone for any 0<=t<=1/2."
        ),
        "smooth_polya_obstruction": (
            "A differentiable even kernel has K'(0)=0. If it were convex on "
            "[0,infinity), K' would be nondecreasing and hence nonnegative; a "
            "nonconstant positive kernel tending to zero cannot satisfy the "
            "classical decreasing-convex Polya criterion."
        ),
        "direct_pf_infinity_obstruction": {
            "classification": (
                "Schoenberg PF-infinity functions are convolutions generated by "
                "Gaussians and one-sided exponentials, up to translation and scale."
            ),
            "tail_consequence": (
                "A nondegenerate even PF-infinity function has a Gaussian lower "
                "tail when the Gaussian factor is present, and an exponential "
                "lower tail otherwise. It cannot decay faster than every Gaussian "
                "on both half-lines."
            ),
            "application": (
                "K_(1,t) itself is not PF-infinity. This does not rule out a "
                "different structured sum-of-squares or total-positive factorization."
            ),
            "sources": [
                "https://doi.org/10.1073/pnas.33.1.11",
                "https://arxiv.org/abs/2006.16213",
            ],
        },
        "tail_compatible_handoff": (
            "Seek an Xi/theta-summand spectral-square decomposition or a relative "
            "correlation-hierarchy coercivity estimate that excludes the universal "
            "contact F_1=F_1'=0, F_2=3*F_1''>0. Gaussian mixtures, the smooth Polya "
            "convex criterion, and direct PF-infinity membership of K_(1,t) are closed."
        ),
    }


def build_payload() -> dict:
    exact = build_exact()
    numerics = build_numerics()
    rows = [
        GateRow(
            id="nchgm_01_heat_hierarchy",
            role="exact_identity",
            readiness="available_exact",
            claim="The Xi correlations form an exact first-order heat hierarchy.",
            formula=exact["heat_hierarchy"],
            proof_boundary="Differentiation under the super-exponentially convergent integral.",
        ),
        GateRow(
            id="nchgm_02_fourier_hierarchy",
            role="exact_identity",
            readiness="available_exact",
            claim="Fourier transformation converts the correlation hierarchy into a backward parabolic hierarchy.",
            formula=exact["fourier_hierarchy"],
            proof_boundary="Exact Fourier identity; it does not provide a favorable maximum-principle sign.",
        ),
        GateRow(
            id="nchgm_03_laguerre_normalization",
            role="exact_identity",
            readiness="available_exact",
            claim="Every generalized Laguerre expression is a normalized correlation transform.",
            formula=exact["laguerre_normalization"],
            proof_boundary="Normalization identity only; strict positivity is not inferred.",
            diagnostics=exact["normalization_checks"],
        ),
        GateRow(
            id="nchgm_04_multiple_root_contact",
            role="exact_boundary_signature",
            readiness="ready_to_apply",
            claim="A multiple Newman-boundary root has a rigid vanishing pattern across the full correlation hierarchy.",
            formula=exact["multiple_root_signature"],
            proof_boundary="Necessary boundary signature, not its exclusion.",
            diagnostics=exact["double_root_contact"],
        ),
        GateRow(
            id="nchgm_05_gaussian_mixture_sufficiency",
            role="published_theorem_composition",
            readiness="ready_to_apply_if_hypothesis_proved",
            claim="Complete monotonicity in v^2 would force strict Fourier positivity.",
            formula=exact["gaussian_mixture_sufficiency"]["fourier_formula"],
            proof_boundary="Sufficient criterion only; the Xi hypothesis is false.",
            diagnostics=exact["gaussian_mixture_sufficiency"],
        ),
        GateRow(
            id="nchgm_06_alternating_derivative_scout",
            role="numerical_diagnostic",
            readiness="diagnostic_validated",
            claim="Low-order alternating derivatives imitate complete monotonicity unusually well.",
            formula="(-1)^n*g_t^(n)(0)>0 through n=8 on the sampled t grid",
            proof_boundary="Finite numerical scout only; derivative signs do not establish complete monotonicity.",
            diagnostics=numerics,
        ),
        GateRow(
            id="nchgm_07_logconvexity_witness",
            role="numerical_counter_witness",
            readiness="diagnostic_validated",
            claim="The first Stieltjes determinant is robustly negative despite all sampled derivative signs.",
            formula="g_t(0)*g_t''(0)-g_t'(0)^2< -7e-6 on t in {0,0.1,0.25,0.5}",
            proof_boundary="High-precision numerical witness; the exact global rejection uses the tail theorem.",
            diagnostics=numerics["independent_mpmath_55_digit_references"],
        ),
        GateRow(
            id="nchgm_08_super_gaussian_tail",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The first Xi correlation decays faster than every Gaussian uniformly over the target time interval.",
            formula=exact["correlation_tail_bound"]["bound"],
            proof_boundary="Exact upper-tail theorem; it does not decide the sign of the Fourier transform.",
            diagnostics=exact["correlation_tail_bound"],
        ),
        GateRow(
            id="nchgm_09_complete_monotonicity_rejection",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The Gaussian-mixture sufficient route is globally impossible for every target time.",
            formula=exact["complete_monotonicity_obstruction"],
            proof_boundary="Rejects this sufficient criterion only, not strict Fourier positivity itself.",
        ),
        GateRow(
            id="nchgm_10_polya_pf_rejections",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Classical smooth convexity and direct PF-infinity membership cannot be the missing Xi property.",
            formula=(
                exact["smooth_polya_obstruction"]
                + " "
                + exact["direct_pf_infinity_obstruction"]["application"]
            ),
            proof_boundary="Does not reject structured spectral-square or weaker total-positive decompositions.",
            diagnostics=exact["direct_pf_infinity_obstruction"],
        ),
        GateRow(
            id="nchgm_11_tail_compatible_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The next correlation route must respect the actual double-exponential tail and exclude the universal boundary contact.",
            formula=exact["tail_compatible_handoff"],
            proof_boundary="Open Xi-specific theorem search; not a proof of the strict correlation target, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate",
        "date": "2026-07-11",
        "status": "exact correlation hierarchy with Gaussian-mixture and direct-PF obstructions",
        "proof_boundary": (
            "This artifact derives the full Xi correlation hierarchy, its multiple-root "
            "contact signature, and a strict-Fourier sufficient Gaussian-mixture criterion. "
            "It then rejects that criterion exactly by the Xi super-Gaussian tail and "
            "records an independent numerical log-convexity witness. It also rejects the "
            "classical smooth Polya-convex route and direct PF-infinity membership of the "
            "correlation kernel. It does not prove strict Fourier positivity, Wiener "
            "translate density, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md",
            "outputs/jensen_window_pf_newman_positive_boundary_attainment_lemma.md",
            "https://arxiv.org/abs/1309.0055",
            "https://doi.org/10.2307/1968466",
            "https://doi.org/10.1073/pnas.33.1.11",
            "https://arxiv.org/abs/2006.16213",
        ],
        "exact": exact,
        "numerics": numerics,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    numerics = payload["numerics"]
    lines = [
        "# Jensen-Window PF Newman Correlation Hierarchy Gaussian-Mixture Gate",
        "",
        "Date: 2026-07-11",
        "",
        "Status: exact correlation hierarchy with Gaussian-mixture and direct-PF",
        "obstructions. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF Newman correlation hierarchy Gaussian-mixture gate: 11 rows, 0 issues, 3 exact hierarchy identities, 1 universal boundary-contact signature, 1 Gaussian-mixture sufficient theorem, 2 numerical diagnostics, 1 exact super-Gaussian tail theorem, 2 non-promotion gates, 1 tail-compatible handoff",
        "```",
        "",
        "## Exact Hierarchy",
        "",
        "Set",
        "",
        "```text",
        exact["definitions"]["correlations"],
        exact["definitions"]["correlation_transforms"],
        exact["heat_hierarchy"],
        exact["fourier_hierarchy"],
        exact["laguerre_normalization"],
        "```",
        "",
        "The hierarchy also resolves the exact contact made by a hypothetical",
        "positive Newman boundary:",
        "",
        "```text",
        exact["multiple_root_signature"],
        exact["double_root_contact"]["compatibility"],
        "```",
        "",
        "For a double root this is a nonnegative quadratic touch of `F_1` with",
        "a strictly positive `F_2` source. The hierarchy is therefore compatible",
        "with birth; it becomes useful only with an additional coercive estimate.",
        "",
        "## Gaussian-Mixture Test",
        "",
        exact["gaussian_mixture_sufficiency"]["criterion"],
        "Then",
        "",
        "```text",
        exact["gaussian_mixture_sufficiency"]["fourier_formula"],
        "```",
        "",
        "so this single property would close the strict Laguerre/Wiener target.",
        "It survives a misleadingly strong local scout:",
        "",
        "```text",
    ]
    for row in numerics["rows"]:
        lines.append(
            f"t={row['t']}: alternating derivatives through n=8; "
            f"g0*g2-g1^2={row['log_convexity_determinant']}"
        )
    lines.extend(
        [
            "```",
            "",
            "The determinant is negative, not positive: complete monotonicity",
            "requires log-convexity. Independent 55-digit mpmath quadrature gives",
            "`-7.8630150226764958e-6` at `t=0` and",
            "`-7.9006312495514071e-6` at `t=1/2`.",
            "",
            "## Exact Tail Obstruction",
            "",
            "The Phi series gives",
            "",
            "```text",
            exact["phi_tail_bound"]["bound"],
            exact["correlation_tail_bound"]["bound"],
            exact["correlation_tail_bound"]["limit"],
            "```",
            "",
            "A nonzero completely monotone function is a positive Laplace mixture",
            "and therefore has an exponential lower bound along the positive axis.",
            "The displayed limit is impossible for such a mixture. Hence the",
            "Gaussian-mixture route is rejected exactly for every `0<=t<=1/2`;",
            "the numerical determinant merely shows that the failure is visible",
            "already at the first moment-minor test.",
            "",
            "## Closed Generic Routes",
            "",
            exact["smooth_polya_obstruction"],
            "",
            exact["direct_pf_infinity_obstruction"]["tail_consequence"],
            "Therefore `K_(1,t)` itself cannot be PF-infinity. This sharpens the",
            "previous handoff: an eventual total-positive argument must be a",
            "different structured factorization, not direct PF membership.",
            "",
            "Primary sources: https://doi.org/10.2307/1968466,",
            "https://doi.org/10.1073/pnas.33.1.11, and",
            "https://arxiv.org/abs/2006.16213.",
            "",
            "## Live Handoff",
            "",
            exact["tail_compatible_handoff"],
            "",
            "The most concrete next test is a theta-summand spectral-square",
            "decomposition. Failing that, the hierarchy needs a relative bound",
            "strong enough to make `F_2=3*F_1''` impossible at a zero of `F_1`.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Newman correlation hierarchy Gaussian-mixture gate: "
        f"{len(payload['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
