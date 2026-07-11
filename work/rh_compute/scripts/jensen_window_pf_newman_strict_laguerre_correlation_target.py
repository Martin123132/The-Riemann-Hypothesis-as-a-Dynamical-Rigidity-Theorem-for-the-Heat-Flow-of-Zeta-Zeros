#!/usr/bin/env python3
"""Recast the Newman upper direction as a strict Laguerre correlation target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_strict_laguerre_correlation_target.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md"


@dataclass(frozen=True)
class TargetRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def build_exact() -> dict:
    a, sigma, xi = sp.symbols("a sigma xi", positive=True)
    triangular_transform = 4 * sp.sin(a * xi / 2) ** 2 / xi**2
    gaussian_transform = sp.sqrt(2 * sp.pi) * sigma * sp.exp(
        -sigma**2 * xi**2 / 2
    )
    guard_transform = sp.factor(triangular_transform * gaussian_transform)

    xi0 = 2 * sp.pi
    guard_a1_sigma2 = sp.simplify(guard_transform.subs({a: 1, sigma: 2}))
    if sp.simplify(guard_a1_sigma2.subs(xi, xi0)) != 0:
        raise RuntimeError("guard Fourier zero failed")
    if sp.simplify(sp.diff(guard_a1_sigma2, xi).subs(xi, xi0)) != 0:
        raise RuntimeError("guard double-zero slope failed")
    second_at_zero = sp.simplify(
        sp.limit(sp.diff(guard_a1_sigma2, xi, 2), xi, xi0)
    )
    if second_at_zero == 0:
        raise RuntimeError("guard double-zero curvature failed")

    curvature_bound = sp.factor((a**2 - sigma**2) / sigma**4)
    guard_curvature_bound = curvature_bound.subs({a: 1, sigma: 2})
    if guard_curvature_bound != -sp.Rational(3, 16):
        raise RuntimeError("guard log-curvature bound failed")

    return {
        "definitions": {
            "deformed_kernel": "phi_t(u)=exp(t*u^2)*Phi(u), u in R",
            "heat_function": "H_t(x)=integral_0^infinity phi_t(u)*cos(x*u)du",
            "first_laguerre": "L_t(x)=H_t'(x)^2-H_t(x)*H_t''(x)",
            "correlation": (
                "K_(1,t)(v)=integral_R phi_t(s+v)*phi_t(s-v)*s^2 ds"
            ),
        },
        "strict_laguerre_equivalence": (
            "Lambda<=0 if and only if L_t(x)>0 for every real x and every "
            "0<t<=1/2. The forward implication uses the simple real-zero "
            "factorization for t>Lambda; the reverse implication uses the finite "
            "multiple zero of H_Lambda forced when Lambda>0."
        ),
        "factorization_identity": (
            "When H_t has simple real roots r_j, L_t(x)/H_t(x)^2="
            "sum_j 1/(x-r_j)^2 away from roots, and L_t(r_j)=H_t'(r_j)^2>0."
        ),
        "origin_sign": (
            "H_t(0)>0, H_t'(0)=0, H_t''(0)=-integral_0^infinity "
            "u^2*phi_t(u)du<0, hence L_t(0)>0."
        ),
        "correlation_identity": (
            "L_t(x)=integral_R K_(1,t)(v)*cos(2*x*v)dv; equivalently, "
            "Fourier[K_(1,t)](xi)=L_t(xi/2)."
        ),
        "correlation_normalization": {
            "full_transform": "F_t(x)=integral_R phi_t(u)*exp(i*x*u)du=2*H_t(x)",
            "full_laguerre": "L_1[F_t](x)=4*L_t(x)",
            "jacobian_identity": (
                "L_1[F_t](x)=4*integral_R K_(1,t)(v)*cos(2*x*v)dv"
            ),
        },
        "wiener_density": {
            "source": "https://arxiv.org/abs/1606.05011",
            "statement": (
                "For K in L1(R), the linear span of its translates is dense in "
                "L1(R) if and only if its Fourier transform has no real zero."
            ),
            "application": (
                "Translations of K_(1,t) are dense iff L_t has no real zero. "
                "Since L_t(0)>0, zero-freeness is equivalent to L_t(x)>0 for all x."
            ),
        },
        "density_rh_equivalence": (
            "Lambda<=0 if and only if, for every 0<t<=1/2, the translations "
            "of K_(1,t) are dense in L1(R)."
        ),
        "positive_definite_scope": (
            "At a hypothetical boundary H_Lambda is Laguerre-Polya, so L_Lambda>=0 "
            "and K_(1,Lambda) is positive definite; a multiple zero c still gives "
            "Fourier[K_(1,Lambda)](2c)=L_Lambda(c)=0. Positive definiteness alone "
            "does not provide the required zero-freeness."
        ),
        "published_correlation_source": {
            "source": "https://arxiv.org/abs/1309.0055",
            "theorem": "Theorem 3.7 and equations (3.7), (3.14)",
            "scope": (
                "Correlation kernels represent generalized Laguerre expressions; "
                "positive definiteness gives nonnegativity, not strict positivity."
            ),
        },
        "strict_logconcavity_guard": {
            "parameters": "a=1, sigma=2",
            "triangle": "T_a(y)=(a-|y|)_+ = 1_[-a/2,a/2]*1_[-a/2,a/2]",
            "gaussian": "G_sigma(x)=exp(-x^2/(2*sigma^2))",
            "kernel": (
                "K_(a,sigma)(x)=integral_(-a)^a (a-|y|)*"
                "exp(-(x-y)^2/(2*sigma^2))dy"
            ),
            "conditional_curvature": (
                "With p_x(y) proportional to T_a(y)G_sigma(x-y), "
                "(log K)''=Var_x(Y)/sigma^4-1/sigma^2"
            ),
            "strict_bound": (
                "Since Y lies in [-a,a], Var_x(Y)<=a^2; for sigma>a, "
                "(log K)''<=(a^2-sigma^2)/sigma^4<0."
            ),
            "fourier_transform": str(guard_transform),
            "a1_sigma2_transform": str(guard_a1_sigma2),
            "double_zeros": "xi=2*pi*n/a for every nonzero integer n",
            "positive_definite": (
                "The transform is nonnegative, so K_(a,sigma) is positive definite."
            ),
            "translation_failure": (
                "The double Fourier zeros make its translates non-dense by Wiener."
            ),
            "curvature_bound_a1_sigma2": str(guard_curvature_bound),
            "second_derivative_at_2pi": str(second_at_zero),
            "scope_warning": (
                "This Gaussian-tail model blocks promotion from strict log-concavity "
                "and positive definiteness alone; it does not reproduce the Xi tail."
            ),
        },
        "nonpromotion_decision": (
            "Pointwise positivity, evenness, smoothness, strict log-concavity, and "
            "positive definiteness of a kernel do not force a zero-free Fourier "
            "transform. The Xi route needs a stronger correlation property such as "
            "Wiener translate density, strict Fourier positivity, or an Xi-specific "
            "total-positivity factorization."
        ),
        "open_handoff": (
            "Prove uniformly for every 0<t<=1/2 that Fourier[K_(1,t)] has no real "
            "zero, equivalently that translations of K_(1,t) are dense in L1(R). "
            "A viable proof must use Xi-specific structure beyond generic strict "
            "log-concavity and positive definiteness."
        ),
        "checks": {
            "triangular_transform": str(triangular_transform),
            "gaussian_transform": str(gaussian_transform),
            "guard_transform": str(guard_transform),
            "guard_zero": str(sp.simplify(guard_a1_sigma2.subs(xi, xi0))),
            "guard_zero_derivative": str(
                sp.simplify(sp.diff(guard_a1_sigma2, xi).subs(xi, xi0))
            ),
            "guard_zero_second_derivative": str(second_at_zero),
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        TargetRow(
            id="nslc_01_boundary_laguerre_zero",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="A positive Newman boundary forces a zero of the first Laguerre expression.",
            formula="H_Lambda(c)=H_Lambda'(c)=0 => L_Lambda(c)=0",
            proof_boundary="Uses the finite positive-boundary attainment lemma.",
        ),
        TargetRow(
            id="nslc_02_strict_laguerre_equivalence",
            role="exact_equivalence",
            readiness="ready_to_apply",
            claim="Strict first-Laguerre positivity throughout positive time is equivalent to the Newman upper direction.",
            formula=exact["strict_laguerre_equivalence"],
            proof_boundary="All real x and the full continuum 0<t<=1/2 are required.",
            diagnostics={
                "factorization": exact["factorization_identity"],
                "origin": exact["origin_sign"],
            },
        ),
        TargetRow(
            id="nslc_03_correlation_identity",
            role="exact_identity",
            readiness="available_exact",
            claim="The first Laguerre expression is exactly the Fourier transform of one positive Xi correlation kernel.",
            formula=exact["correlation_identity"],
            proof_boundary="Fubini and the midpoint/difference change of variables.",
            diagnostics=exact["correlation_normalization"],
        ),
        TargetRow(
            id="nslc_04_wiener_density",
            role="published_theorem",
            readiness="available_published",
            claim="Wiener translate density is equivalent to zero-freeness of the correlation transform.",
            formula=exact["wiener_density"]["application"],
            proof_boundary="Classical Wiener L1 Tauberian theorem; K_(1,t) is integrable.",
            diagnostics=exact["wiener_density"],
        ),
        TargetRow(
            id="nslc_05_density_rh_equivalence",
            role="exact_equivalence",
            readiness="ready_to_apply",
            claim="The missing Newman direction is equivalent to an all-positive-time correlation-density theorem.",
            formula=exact["density_rh_equivalence"],
            proof_boundary="Equivalent target only; density is not proved.",
        ),
        TargetRow(
            id="nslc_06_positive_definite_scope",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Positive definiteness yields nonnegative Laguerre expressions but does not exclude a boundary zero.",
            formula=exact["positive_definite_scope"],
            proof_boundary="Separates nonnegativity from strict zero-freeness.",
            diagnostics=exact["published_correlation_source"],
        ),
        TargetRow(
            id="nslc_07_strict_logconcavity_guard",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="A smooth positive even strictly log-concave positive-definite kernel can have double Fourier zeros.",
            formula=exact["strict_logconcavity_guard"]["a1_sigma2_transform"],
            proof_boundary="Exact Gaussian-tail convolution model, not the Xi kernel.",
            diagnostics=exact["strict_logconcavity_guard"],
        ),
        TargetRow(
            id="nslc_08_generic_kernel_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Generic kernel shape properties cannot be promoted to Wiener density.",
            formula=exact["nonpromotion_decision"],
            proof_boundary="Does not rule out stronger Xi tail or total-positivity structure.",
        ),
        TargetRow(
            id="nslc_09_first_correlation_target",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim="It is sufficient to prove zero-freeness of the first Xi correlation transform for all positive times up to 1/2.",
            formula="Fourier[K_(1,t)](xi)=L_t(xi/2)>0 for all xi in R, 0<t<=1/2",
            proof_boundary="Open Xi-specific strict correlation theorem.",
        ),
        TargetRow(
            id="nslc_10_xi_density_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The next route should attack translate density using Xi-specific correlation structure.",
            formula=exact["open_handoff"],
            proof_boundary="Not a proof of density, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_strict_laguerre_correlation_target",
        "date": "2026-07-11",
        "status": "exact strict-Laguerre/Wiener equivalence with generic-kernel guard",
        "proof_boundary": (
            "This artifact composes positive-boundary attainment with the first "
            "Laguerre correlation identity and Wiener's theorem to give an RH-equivalent "
            "translation-density target for every 0<t<=1/2. It supplies an exact smooth "
            "strict-log-concavity and positive-definiteness countermodel with double "
            "Fourier zeros. It does not prove the Xi correlation transform zero-free, "
            "translation density, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_positive_boundary_attainment_lemma.md",
            "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
            "outputs/formal_core.md",
            "https://arxiv.org/abs/1309.0055",
            "https://arxiv.org/abs/1606.05011",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    guard = exact["strict_logconcavity_guard"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Strict-Laguerre Correlation Target",
            "",
            "Date: 2026-07-11",
            "",
            "Status: exact strict-Laguerre/Wiener equivalence with a generic-kernel",
            "guard. This is not a proof of RH or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_strict_laguerre_correlation_target.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_strict_laguerre_correlation_target.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_strict_laguerre_correlation_target.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF Newman strict-Laguerre correlation target: 10 rows, 0 issues, 1 strict-Laguerre equivalence, 1 exact correlation identity, 1 Wiener-density equivalence, 1 RH-equivalent density target, 1 exact strict-log-concavity/positive-definiteness countermodel, 2 non-promotion gates, 1 open Xi handoff",
            "```",
            "",
            "## Strict Laguerre Target",
            "",
            "Set",
            "",
            "```text",
            exact["definitions"]["deformed_kernel"],
            exact["definitions"]["heat_function"],
            exact["definitions"]["first_laguerre"],
            exact["definitions"]["correlation"],
            "```",
            "",
            "The positive-boundary attainment theorem turns the first Laguerre",
            "inequality into an exact endgame:",
            "",
            "```text",
            exact["strict_laguerre_equivalence"],
            "```",
            "",
            "This is weaker in form than proving every `H_t` is Laguerre-Polya:",
            "only one strict real-axis inequality is requested, but it must hold",
            "for the complete positive-time continuum.",
            "",
            "## Correlation Identity",
            "",
            "The midpoint/difference change of variables gives exactly",
            "",
            "```text",
            exact["correlation_identity"],
            "```",
            "",
            "There is no missing normalization factor. Wiener's theorem now yields",
            "",
            "```text",
            exact["wiener_density"]["application"],
            exact["density_rh_equivalence"],
            "```",
            "",
            "Primary sources: https://arxiv.org/abs/1309.0055 and",
            "https://arxiv.org/abs/1606.05011.",
            "",
            "At a hypothetical positive boundary the correlation kernel is already",
            "positive definite, but its transform touches zero at the multiple root.",
            "The missing property is therefore zero-freeness or translate density,",
            "not ordinary nonnegativity.",
            "",
            "## Exact Shape Guard",
            "",
            "Let",
            "",
            "```text",
            guard["triangle"],
            guard["gaussian"],
            guard["kernel"],
            "```",
            "",
            "For `a=1`, `sigma=2`, conditional differentiation gives",
            "",
            "```text",
            guard["conditional_curvature"],
            guard["strict_bound"],
            f"(log K)''<={guard['curvature_bound_a1_sigma2']}",
            "```",
            "",
            "so `K` is smooth, positive, even, and strictly log-concave. Its",
            "Fourier transform is nonnegative, hence `K` is positive definite, but",
            "",
            "```text",
            guard["a1_sigma2_transform"],
            guard["double_zeros"],
            "```",
            "",
            "The zeros are double, so the translates of `K` are not dense. This",
            "blocks promotion from strict log-concavity plus positive definiteness.",
            "It does not block stronger Xi-specific tail or correlation structure.",
            "",
            "## Live Handoff",
            "",
            exact["open_handoff"],
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Jensen-window PF Newman strict-Laguerre correlation target: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
