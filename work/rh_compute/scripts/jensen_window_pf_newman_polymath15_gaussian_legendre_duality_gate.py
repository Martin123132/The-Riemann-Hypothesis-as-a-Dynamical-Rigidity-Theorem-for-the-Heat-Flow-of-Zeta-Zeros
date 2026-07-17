#!/usr/bin/env python3
"""Build the Gaussian/Legendre duality gate for the Newman scaled layer."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

from jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem import (
    exponent_pair,
)
from jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate import (
    C_STAR,
    R_STAR,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate.md"
)

POLYMATH_SOURCE_URL = "https://arxiv.org/abs/1904.12438"
EXPONENT_PAIR_SOURCE_URL = "https://arxiv.org/abs/2306.05599"
CURRENT_PHASE = Fraction(220_633, 310_306)
C_TWO = Fraction(2, 1)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def fraction_record(value: Fraction) -> dict:
    return {
        "exact": f"{value.numerator}/{value.denominator}",
        "decimal": f"{float(value):.15f}",
    }


def dual_point(scaled_time: Fraction) -> dict:
    radius = R_STAR
    sigma = Fraction(1, 2) + scaled_time / 4
    shift = scaled_time * radius / 4
    line = sigma - shift
    partial_sum_exponent = CURRENT_PHASE - line * radius
    gaussian_cost = 2 * shift**2 / scaled_time
    return {
        "scaled_time": fraction_record(scaled_time),
        "radius": fraction_record(radius),
        "line": fraction_record(line),
        "gaussian_shift": fraction_record(shift),
        "partial_sum_exponent": fraction_record(partial_sum_exponent),
        "gaussian_cost": fraction_record(gaussian_cost),
        "net_exponent": fraction_record(partial_sum_exponent - gaussian_cost),
    }


def build_exact() -> dict:
    c, r, u = sp.symbols("c r u", positive=True, real=True)
    dual_quadratic = u * r - 2 * u**2 / c
    saddle = c * r / 4
    saddle_value = sp.simplify(dual_quadratic.subs(u, saddle))
    if saddle_value != c * r**2 / 8:
        raise RuntimeError("Gaussian Legendre saddle identity failed")
    if sp.diff(dual_quadratic, u, 2) != -4 / c:
        raise RuntimeError("Gaussian saddle concavity failed")

    star = dual_point(C_STAR)
    c_two = dual_point(C_TWO)
    if Fraction(star["net_exponent"]["exact"]) != 0:
        raise RuntimeError("c_star dual equality failed")
    if Fraction(c_two["net_exponent"]["exact"]) != Fraction(
        3_133_668_399, 48_144_906_818
    ):
        raise RuntimeError("c=2 dual deficit failed")

    slope_two = exponent_pair(2)[1] - exponent_pair(2)[0]
    slope_one = exponent_pair(1)[1] - exponent_pair(1)[0]
    q_star = Fraction(star["line"]["exact"])
    q_two = Fraction(c_two["line"]["exact"])
    if not (slope_two >= q_star >= slope_one):
        raise RuntimeError("q_star does not expose the pair-2/pair-1 corner")
    if not (slope_two >= q_two >= slope_one):
        raise RuntimeError("q_2 does not expose the pair-2/pair-1 corner")

    return {
        "gaussian_identity": (
            "exp(t*y^2/4)=(pi*t)^(-1/2)*integral_R "
            "exp(-u^2/t+u*y)du; hence D_(N,t)(s_*)=(pi*t)^(-1/2)*"
            "integral_R exp(-u^2/t)S_N(s_*-u)du, "
            "S_N(z)=sum_(n<=N)n^(-z)"
        ),
        "moment_identity": (
            "The same finite interchange gives D_k=(pi*t)^(-1/2)*"
            "integral_R exp(-u^2/t)S_(N,k)(s_*-u)du for k=0,1,2"
        ),
        "scaled_gaussian": (
            "For t=c/L and L=2log(N)+o(1), exp(-u^2/t)="
            "N^(-2u^2/c+o(1))"
        ),
        "partial_sum_profile": (
            "If Phi(r) bounds the phase sum on M=N^r, define "
            "E(q)=sup_(0<=r<=1)(Phi(r)-q*r)"
        ),
        "dual_identity": (
            "With sigma=1/2+c/4, sup_u{E(sigma-u)-2u^2/c}="
            "sup_r{Phi(r)-r/2-c*r*(2-r)/8}"
        ),
        "saddle_identity": (
            "For fixed r, sup_u{u*r-2u^2/c}=c*r^2/8 at u=c*r/4"
        ),
        "threshold_equivalence": (
            "A Gaussian partial-sum proof using the same pointwise envelope "
            "has strict decay exactly when the weighted dyadic frontier has "
            "strict decay; neither formulation is stronger"
        ),
        "critical_star": star,
        "critical_c2": c_two,
        "corner_slopes": {
            "pair_2": fraction_record(slope_two),
            "pair_1": fraction_record(slope_one),
        },
        "current_nonpromotion": (
            "At c=c_* the pair-2/pair-1 corner gives equality, not decay. "
            "Rewriting the current eleven-pair hull as a Gaussian heat average "
            "therefore cannot lower c_*"
        ),
        "required_upgrade": (
            "Lowering c_* requires a genuinely smaller partial-sum profile "
            "E(q) near q=0.6207094699015768 and across the exposed q-range, "
            "or an argument exploiting cancellation beyond pointwise partial-sum "
            "majorants; reaching c=2 still leaves the fixed c<2 phase wall"
        ),
    }


def build_artifact() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="np15gld_01_gaussian_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The finite heat-weighted Dirichlet sum has an exact Gaussian shift representation.",
            formula=exact["gaussian_identity"],
            proof_boundary="Finite-sum Gaussian integration; no asymptotic estimate is used.",
        ),
        GateRow(
            id="np15gld_02_moment_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The Gaussian representation transfers termwise through the two logarithmic moments needed for the curvature jet.",
            formula=exact["moment_identity"],
            proof_boundary="Exact for the finite cutoff only.",
        ),
        GateRow(
            id="np15gld_03_scaled_cost",
            role="exact_asymptotic_geometry",
            readiness="ready_to_apply",
            claim="In scaled Newman coordinates the Gaussian shift has a quadratic N-exponent cost.",
            formula=exact["scaled_gaussian"],
            proof_boundary="Leading exponent with bounded c and L tending to infinity.",
        ),
        GateRow(
            id="np15gld_04_partial_sum_profile",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="The full partial-sum bound is the Legendre profile of the dyadic phase envelope.",
            formula=exact["partial_sum_profile"],
            proof_boundary="Pointwise exponent bookkeeping; logarithmic losses are suppressed by eta.",
        ),
        GateRow(
            id="np15gld_05_saddle_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Each dyadic radius selects one exact Gaussian shift saddle.",
            formula=exact["saddle_identity"],
            proof_boundary="Exact concave quadratic maximization.",
        ),
        GateRow(
            id="np15gld_06_dual_equivalence",
            role="exact_lemma",
            readiness="ready_to_apply",
            claim="The Gaussian and weighted-dyadic decay frontiers are exactly equal.",
            formula=exact["dual_identity"],
            proof_boundary="Equality of suprema for the same pointwise phase envelope.",
        ),
        GateRow(
            id="np15gld_07_cstar_contact",
            role="exact_frontier_contact",
            readiness="ready_to_apply",
            claim="The known c_* obstruction is an exact Gaussian/partial-sum contact at one exposed critical-strip line.",
            formula=json.dumps(exact["critical_star"], sort_keys=True),
            proof_boundary="Exact rational contact for the current published hull.",
            diagnostics=exact["critical_star"],
        ),
        GateRow(
            id="np15gld_08_c2_deficit",
            role="exact_frontier_target",
            readiness="ready_to_apply",
            claim="At c=2 the same exposed corner retains the exact previously identified positive exponent deficit.",
            formula=json.dumps(exact["critical_c2"], sort_keys=True),
            proof_boundary="Necessary local target only, not an all-line cancellation theorem.",
            diagnostics=exact["critical_c2"],
        ),
        GateRow(
            id="np15gld_09_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="A Gaussian heat-semigroup rewrite using the current partial-sum bounds cannot improve c_*.",
            formula=exact["current_nonpromotion"],
            proof_boundary="Route-equivalence gate only; it does not rule out stronger cancellation or Xi-specific structure.",
        ),
        GateRow(
            id="np15gld_10_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The dual picture identifies the precise partial-sum profile that must improve, while preserving the separate inner phase wall.",
            formula=exact["required_upgrade"],
            proof_boundary="The required analytic improvement remains open; not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate",
        "date": "2026-07-17",
        "status": (
            "exact Gaussian/Legendre equivalence for the Newman cancellation "
            "frontier; one route nonpromotion gate and one open analytic "
            "handoff, not a proof of Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "This artifact proves an exact finite Gaussian identity and the "
            "Legendre equivalence of two uses of the same phase-sum envelope. "
            "It does not improve a partial-sum estimate, establish strict decay "
            "at c_*, close c=2, cross the fixed c<2 phase wall, prove "
            "Wronskian separation, prove Lambda<=0, or prove RH."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "sources": [
            POLYMATH_SOURCE_URL,
            EXPONENT_PAIR_SOURCE_URL,
            "outputs/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.md",
            "outputs/jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    star = exact["critical_star"]
    c_two = exact["critical_c2"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Gaussian/Legendre Duality Gate",
            "",
            "Date: 2026-07-17",
            "",
            "Status: exact equivalence and route nonpromotion gate. The required",
            "cancellation improvement remains open. This is not a proof of `Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate.py",
            "```",
            "",
            "## Exact Gaussian Rewrite",
            "",
            "For every finite cutoff, completing the Gaussian square gives",
            "",
            "```text",
            exact["gaussian_identity"],
            exact["moment_identity"],
            "```",
            "",
            "This is tempting as an alternative to dyadic exponent pairs, but its",
            "scaled saddle contains exactly the same optimization.",
            "",
            "## Legendre Equivalence",
            "",
            "```text",
            exact["scaled_gaussian"],
            exact["partial_sum_profile"],
            exact["saddle_identity"],
            exact["dual_identity"],
            "```",
            "",
            "The equality follows by interchanging the two suprema and maximizing",
            "the displayed concave quadratic in `u`. Therefore",
            "",
            "```text",
            exact["threshold_equivalence"],
            "```",
            "",
            "## Exact Contact",
            "",
            "At the current threshold, the exposed pair-2/pair-1 corner maps to",
            "",
            "```text",
            f"r_*={star['radius']['exact']}",
            f"q_*={star['line']['exact']}={star['line']['decimal']}...",
            f"u_*={star['gaussian_shift']['exact']}",
            f"E(q_*)={star['partial_sum_exponent']['exact']}",
            f"2u_*^2/c_*={star['gaussian_cost']['exact']}",
            f"net exponent={star['net_exponent']['exact']}",
            "```",
            "",
            "Thus `c_*` is a true equality point for this proof architecture, not an artifact",
            "of whether the estimates are written blockwise or as a",
            "Gaussian heat average.",
            "",
            "At `c=2`, the same corner gives",
            "",
            "```text",
            f"q_2={c_two['line']['exact']}={c_two['line']['decimal']}...",
            f"E(q_2)={c_two['partial_sum_exponent']['exact']}",
            f"Gaussian cost={c_two['gaussian_cost']['exact']}",
            f"deficit={c_two['net_exponent']['exact']}={c_two['net_exponent']['decimal']}...",
            "```",
            "",
            "This reproduces the cancellation/zero-free wall gate independently.",
            "",
            "## Nonpromotion And Handoff",
            "",
            "```text",
            exact["current_nonpromotion"],
            exact["required_upgrade"],
            "```",
            "",
            "A useful next theorem must therefore improve the actual partial-sum",
            "profile, exploit cancellation between estimates that the pointwise",
            "profile discards, or attack the corrected Xi phase directly. Merely",
            "rewriting the same hull as a heat-semigroup estimate does not move the",
            "boundary.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 Gaussian/Legendre duality gate: "
        f"{len(artifact['rows'])} rows, 1 exact Gaussian identity, "
        "1 exact Legendre equivalence, 1 c_* equality point, "
        "1 c=2 deficit, 1 nonpromotion gate"
    )


if __name__ == "__main__":
    main()
