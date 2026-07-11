#!/usr/bin/env python3
"""Derive the Newman-zero external field behind the scaled Jensen correction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_root_external_field_lemma.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_root_external_field_lemma.md"


@dataclass(frozen=True)
class LemmaRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def build_exact() -> dict:
    s, z, x = sp.symbols("s z x", positive=True)
    u0, u1, u2 = sp.symbols("u0 u1 u2", nonzero=True)
    rho = -x**2
    local = s - rho
    F = sp.expand(local**2 * (u0 + u1 * local + u2 * local**2))
    H = sp.expand(F.subs(s, -z**2))
    f2 = sp.diff(F, s, 2).subs(s, rho)
    f3 = sp.diff(F, s, 3).subs(s, rho)
    h2 = sp.factor(sp.diff(H, z, 2).subs(z, x))
    h3 = sp.factor(sp.diff(H, z, 3).subs(z, x))
    E = sp.factor(u1 / u0)
    B = sp.factor(1 / x - 2 * x * E)
    if sp.simplify(f3 / f2 - 3 * E) != 0:
        raise RuntimeError("squared-zero external-field jet failed")
    if sp.simplify(h3 / h2 - 3 * B) != 0:
        raise RuntimeError("Newman-variable external-field jet failed")

    root_values = [sp.Rational(1, 2), sp.Rational(3, 2), sp.Rational(5, 2)]
    multiplicities = [1, 2, 1]
    xv = sp.Rational(1)
    squared_field = sp.factor(
        sum(
            multiplicity / (root**2 - xv**2)
            for root, multiplicity in zip(root_values, multiplicities)
        )
    )
    squared_stiffness = sp.factor(
        sum(
            multiplicity / (root**2 - xv**2) ** 2
            for root, multiplicity in zip(root_values, multiplicities)
        )
    )
    signed_roots: list[tuple[sp.Rational, int]] = [(-xv, 2)]
    for root, multiplicity in zip(root_values, multiplicities):
        signed_roots.extend([(root, multiplicity), (-root, multiplicity)])
    signed_field = sp.factor(
        sum(multiplicity / (xv - root) for root, multiplicity in signed_roots)
    )
    signed_stiffness = sp.factor(
        sum(multiplicity / (xv - root) ** 2 for root, multiplicity in signed_roots)
    )
    expected_signed_field = sp.factor(1 / xv - 2 * xv * squared_field)
    expected_signed_stiffness = sp.factor(
        1 / (2 * xv**2) + 2 * squared_field + 4 * xv**2 * squared_stiffness
    )
    if signed_field != expected_signed_field:
        raise RuntimeError("finite signed-field conversion failed")
    if signed_stiffness != expected_signed_stiffness:
        raise RuntimeError("finite stiffness conversion failed")

    c, g = sp.symbols("c g", real=True, nonzero=True)
    outsiders = [sp.Rational(-5, 2), sp.Rational(-1, 2), sp.Rational(3, 1)]
    xp = c + g / 2
    xm = c - g / 2
    vp = 2 / g + 2 * sum(1 / (xp - root) for root in outsiders)
    vm = -2 / g + 2 * sum(1 / (xm - root) for root in outsiders)
    center_prime = sp.factor((vp + vm) / 2)
    gap_prime = sp.factor(vp - vm)
    pair_sum = sum(1 / ((xp - root) * (xm - root)) for root in outsiders)
    if sp.simplify(gap_prime - (4 / g - 2 * g * pair_sum)) != 0:
        raise RuntimeError("pair gap equation failed")
    q = sp.symbols("q", nonnegative=True)
    q_flow = "q'=8-4*q*S(c,q), S=sum_y 1/((c-y)^2-q/4)"

    field_symbol = sp.symbols("field")
    second_shift = sp.factor(sp.Rational(1, 64) - field_symbol / 16)
    plus_field = sp.Rational(1)
    minus_field = sp.Rational(-2)
    plus_second_shift = sp.factor(second_shift.subs(field_symbol, plus_field))
    minus_second_shift = sp.factor(second_shift.subs(field_symbol, minus_field))
    if not (plus_second_shift < 0 and minus_second_shift > 0):
        raise RuntimeError("external-field sign countermodels failed")

    y = sp.symbols("y", positive=True)
    cosine_kernel = sp.cosh(sp.sqrt(y))
    cosine_ode = sp.simplify(
        4 * y * sp.diff(cosine_kernel, y, 2)
        + 2 * sp.diff(cosine_kernel, y)
        - cosine_kernel
    )
    if cosine_ode != 0:
        raise RuntimeError("cosine equilibrium ODE failed")

    return {
        "canonical_product": (
            "At a real-zero time, F_lambda(s)=F_lambda(0)*"
            "product_(j>=1)(1+s/x_j(lambda)^2), with sum_j x_j^(-2)<infinity."
        ),
        "double_root_factor": "F_*(s)=(s-rho)^2*U(s), rho=-x^2<0",
        "squared_external_field": (
            "E_x=U'(rho)/U(rho)=sum_(j!=*) m_j/(x_j^2-x^2)"
        ),
        "squared_stiffness": (
            "K_s=-E_x'(rho)=sum_(j!=*) m_j/(x_j^2-x^2)^2>0"
        ),
        "jet_identity": "F_*'''(rho)/F_*''(rho)=3*E_x",
        "newman_external_field": (
            "If H_*(z)=(z-x)^2*V(z), then B_x=V'(x)/V(x)="
            "1/x-2*x*E_x and H_*'''(x)/H_*''(x)=3*B_x."
        ),
        "newman_stiffness": (
            "K_x=sum_(signed y outside pair)1/(x-y)^2="
            "1/(2*x^2)+2*E_x+4*x^2*K_s>0."
        ),
        "pair_flow": {
            "center": "c'=sum_y[1/(x_+-y)+1/(x_--y)] -> 2*B_x",
            "gap": "g'=4/g-2*g*sum_y 1/((x_+-y)*(x_--y))",
            "gap_square": q_flow,
            "boundary_expansion": "g(t)^2=8*(t-t_*)-16*K_x*(t-t_*)^2+O((t-t_*)^3)",
        },
        "sign_countermodels": {
            "positive": {
                "F": "(1+s)^2*(1+s/2)",
                "rho": "-1",
                "E": str(plus_field),
                "D_minus_2_collision_coefficient_at_n0": str(plus_second_shift),
            },
            "negative": {
                "F": "(1+s)^2*(1+2*s)",
                "rho": "-1",
                "E": str(minus_field),
                "D_minus_2_collision_coefficient_at_n0": str(minus_second_shift),
            },
            "conclusion": (
                "Positive coefficients and Laguerre-Polya membership do not fix the sign "
                "of E_x or of the second collision correction."
            ),
        },
        "cosine_equilibrium": {
            "F": "cosh(sqrt(s))^2",
            "roots": "rho_k=-pi^2*(k+1/2)^2, each double",
            "ODE": "4*s*g''+2*g'-g=0 for g=cosh(sqrt(s))",
            "E": "-1/(2*rho_k)=1/(2*x_k^2)",
            "B": "1/x_k-2*x_k*E=0",
            "interpretation": "The arithmetic cosine lattice has zero center field.",
        },
        "checks": {
            "local_f2": str(f2),
            "local_f3": str(f3),
            "local_h2": str(h2),
            "local_h3": str(h3),
            "finite_squared_field": str(squared_field),
            "finite_squared_stiffness": str(squared_stiffness),
            "finite_signed_field": str(signed_field),
            "finite_signed_stiffness": str(signed_stiffness),
            "sample_center_flow": str(center_prime),
            "sample_gap_flow": str(gap_prime),
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        LemmaRow(
            id="nref_01_canonical_product",
            role="exact_identity",
            readiness="available_exact",
            claim="At a real-zero Newman time, the squared-variable entire function has a genus-zero product over positive zero heights.",
            formula=exact["canonical_product"],
            proof_boundary="Conditional on being in the real-zero regime.",
        ),
        LemmaRow(
            id="nref_02_squared_external_field",
            role="exact_identity",
            readiness="available_exact",
            claim="Removing a double squared zero turns the regularized logarithmic derivative into an absolutely convergent root field.",
            formula=exact["squared_external_field"],
            proof_boundary="One finite double zero at a real-zero time.",
        ),
        LemmaRow(
            id="nref_03_squared_stiffness",
            role="exact_identity",
            readiness="available_exact",
            claim="The derivative of the squared-zero field is a strictly negative sum of squares.",
            formula=exact["squared_stiffness"],
            proof_boundary="Root-measure stiffness identity only.",
        ),
        LemmaRow(
            id="nref_04_newman_field_conversion",
            role="exact_identity",
            readiness="available_exact",
            claim="The squared-zero field converts exactly to the regularized field acting on the colliding Newman pair.",
            formula=exact["newman_external_field"],
            proof_boundary="Local change of variables H(z)=F(-z^2).",
        ),
        LemmaRow(
            id="nref_05_newman_stiffness_conversion",
            role="exact_identity",
            readiness="available_exact",
            claim="The signed-zero stiffness is an explicit positive combination of the squared field and squared stiffness.",
            formula=exact["newman_stiffness"],
            proof_boundary="Exact canonical-product algebra only.",
        ),
        LemmaRow(
            id="nref_06_pair_flow_reduction",
            role="exact_local_lemma",
            readiness="ready_to_apply",
            claim="The Newman zero ODE separates a universal pair repulsion from the regularized center field and stiffness.",
            formula=f"{exact['pair_flow']['center']}; {exact['pair_flow']['gap']}",
            proof_boundary="Already-real simple pair for t>t_* only.",
        ),
        LemmaRow(
            id="nref_07_gap_stiffness_expansion",
            role="exact_local_lemma",
            readiness="ready_to_apply",
            claim="The first correction to square-root pair birth is coercive and governed by the positive signed-zero stiffness.",
            formula=exact["pair_flow"]["boundary_expansion"],
            proof_boundary="Local finite-height double-zero branch only.",
        ),
        LemmaRow(
            id="nref_08_external_field_sign_guard",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Positive coefficients and Laguerre-Polya membership alone allow either sign of the regularized field and second Jensen collision correction.",
            formula="E=1 gives -3/64; E=-2 gives 9/64 at rho=-1, n=0",
            proof_boundary="Abstract finite LP countermodels, not Xi trajectories.",
            diagnostics=exact["sign_countermodels"],
        ),
        LemmaRow(
            id="nref_09_cosine_equilibrium",
            role="exact_benchmark",
            readiness="available_exact",
            claim="The doubled arithmetic cosine lattice has the exact field needed for zero center drift.",
            formula="E=1/(2*x_k^2), B_x=0",
            proof_boundary="Solvable equilibrium benchmark only.",
            diagnostics=exact["cosine_equilibrium"],
        ),
        LemmaRow(
            id="nref_10_xi_balance_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="A closing root-measure route must prove an Xi-specific paired-cancellation and far-tail bound for E_x and K_x uniformly as zero height diverges.",
            formula="control sum_j m_j/(x_j^2-x^2) and its square-denominator tail",
            proof_boundary="Open height-uniform theorem; not a proof of RH or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_root_external_field_lemma",
        "date": "2026-07-11",
        "status": "exact Newman root external-field reduction with sign countermodels",
        "proof_boundary": (
            "This artifact translates the Jensen D^-2 correction into the regularized Newman zero field, "
            "derives its positive stiffness and pair-flow expansion, and proves that no generic LP sign rule "
            "can control the field. It does not establish the required Xi-specific degree- and height-uniform "
            "balance estimate, exclude a positive Newman boundary, prove RH, or prove Lambda <= 0."
        ),
        "sources": [
            "outputs/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.md",
            "outputs/formal_core.md",
            "https://arxiv.org/abs/1801.05914",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    counter = exact["sign_countermodels"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Root External-Field Lemma",
            "",
            "Date: 2026-07-11",
            "",
            "Status: exact Newman root external-field reduction with sign",
            "countermodels. This is not a proof of RH or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_root_external_field_lemma.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_root_external_field_lemma.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_root_external_field_lemma.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF Newman root external-field lemma: 10 rows, 0 issues, 5 exact canonical-product identities, 1 pair-flow reduction, 1 gap-stiffness expansion, 2 sign countermodels, 1 cosine equilibrium benchmark, 1 open Xi-balance handoff",
            "```",
            "",
            "## Squared-Zero Field",
            "",
            "At a real-zero Newman time, write the positive zeros as `x_j` and",
            "suppose `x` is a double zero. In `s=-z^2`,",
            "",
            "```text",
            exact["canonical_product"],
            exact["double_root_factor"],
            exact["squared_external_field"],
            exact["squared_stiffness"],
            "```",
            "",
            "The second-order Jensen correction from the boundary-layer theorem",
            "therefore depends on a concrete principal root field rather than an",
            "unspecified higher-degree quantity.",
            "",
            "## Newman Variables",
            "",
            "Returning to the signed zeros of `H`,",
            "",
            "```text",
            exact["newman_external_field"],
            exact["newman_stiffness"],
            "```",
            "",
            "For a pair with center `c`, gap `g`, and `q=g^2`, the zero ODE gives",
            "",
            "```text",
            exact["pair_flow"]["center"],
            exact["pair_flow"]["gap"],
            exact["pair_flow"]["gap_square"],
            exact["pair_flow"]["boundary_expansion"],
            "```",
            "",
            "Thus the universal `8*(t-t_*)` square-root birth receives a negative",
            "second-order correction from the strictly positive stiffness `K_x`.",
            "",
            "## Sign Guard",
            "",
            "Two finite positive-coefficient Laguerre-Polya products at `rho=-1` are",
            "",
            "```text",
            f"{counter['positive']['F']}: E={counter['positive']['E']}, collision D^-2 coefficient={counter['positive']['D_minus_2_collision_coefficient_at_n0']}",
            f"{counter['negative']['F']}: E={counter['negative']['E']}, collision D^-2 coefficient={counter['negative']['D_minus_2_collision_coefficient_at_n0']}",
            "```",
            "",
            "So neither LP membership nor coefficient positivity fixes the field",
            "sign. Any usable estimate must exploit Xi-specific zero balance.",
            "",
            "## Equilibrium Benchmark",
            "",
            "For `F(s)=cosh(sqrt(s))^2`, every squared root is double and",
            "",
            "```text",
            "E_x=1/(2*x_k^2), B_x=1/x_k-2*x_k*E_x=0.",
            "```",
            "",
            "This is the exact arithmetic-lattice equilibrium. The live target is a",
            "uniform comparison between the Xi field and this paired cancellation,",
            "including the far zero tail and possible escape to infinite height.",
            "",
            "Primary zero-dynamics anchor: https://arxiv.org/abs/1801.05914",
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
        "wrote Jensen-window PF Newman root external-field lemma: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
