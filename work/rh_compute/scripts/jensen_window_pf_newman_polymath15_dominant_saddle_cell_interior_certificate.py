#!/usr/bin/env python3
"""Certify exact Newman Laguerre positivity on dominant-saddle cell interiors."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

from jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate import (  # noqa: E402
    PRECISION_BITS,
    SOURCE_URL,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate.md"
)
RHO = Fraction(1, 4)
CELL_DISTANCE = Fraction(1, 2)
L_MIN = Fraction(50)
TL_MIN = Fraction(25)
DELTA = Fraction(1, 8000)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits)


def build_exact() -> dict:
    r0, r1, r2, ell, potential = sp.symbols(
        "r0 r1 r2 ell potential", real=True
    )
    e0 = r0
    e1 = r1 - ell * r0
    e2 = r2 - 2 * ell * r1 + (ell**2 + potential) * r0
    if sp.simplify(e1 - (r1 - ell * e0)) != 0:
        raise RuntimeError("normalized first remainder jet failed")
    if sp.simplify(
        e2 - (r2 - 2 * ell * r1 + (ell**2 + potential) * e0)
    ) != 0:
        raise RuntimeError("normalized second remainder jet failed")

    p0, p1, p2 = sp.symbols("p0 p1 p2", real=True)
    curvature = lambda f0, f1, f2: (
        f1**2 - f0 * f2 + potential * f0**2
    )
    delta = sp.expand(
        curvature(p0 + e0, p1 + e1, p2 + e2)
        - curvature(p0, p1, p2)
    )
    expected = (
        2 * p1 * e1
        + e1**2
        - p0 * e2
        - e0 * p2
        - e0 * e2
        + 2 * potential * p0 * e0
        + potential * e0**2
    )
    if sp.simplify(delta - expected) != 0:
        raise RuntimeError("normalized remainder curvature identity failed")

    return {
        "region": (
            "0<t<=1/2, L=log(x/(4*pi)), t*L>=25, "
            "N=floor(sqrt(x/(4*pi)+t/16)), and x is at distance at least "
            "1/2 from both boundaries of its fixed-N cell"
        ),
        "cell": (
            "4*pi*(N^2-t/16)<=x<4*pi*((N+1)^2-t/16), "
            "dist(x,cell boundary)>=1/2"
        ),
        "raw_main": (
            "S_(N,t)(z)=B_t(z)*f_(N,t)(z), R_(N,t)(z)=H_t(z)-S_(N,t)(z)"
        ),
        "real_axis": (
            "S_(N,t)(x)=A_t(x)*P_(N,t)(x), A_t=|B_t|, "
            "Z_t-P_(N,t)=R_(N,t)/A_t"
        ),
        "reflection": (
            "S_(N,t) is holomorphic and real on the real axis, so "
            "R_(N,t)(conj(z))=conj(R_(N,t)(z))"
        ),
        "scalar_collar": (
            "sup_(|z-x|=1/4)|R_(N,t)(z)|/A_t(x)<1/8000"
        ),
        "cauchy_raw": (
            "|R|/A<=delta, |R'|/A<=4*delta, |R''|/A<=32*delta"
        ),
        "normalized_jets": (
            "eps0<=delta, eps1<=0.35*L*delta, "
            "eps2<=0.13*L^2*delta, delta=1/8000"
        ),
        "main_jets": (
            "|P|<=2.256, |P'|<=0.59L, |P''|<=0.154L^2"
        ),
        "remainder_curvature": (
            "|C[Z_t]-C[P_(N,t)]|<10^-3*L^2"
        ),
        "theorem": (
            "C_t[Z_t](x)>3/40*L^2 and L_t(x)=A_t(x)^2*C_t[Z_t](x)>0"
        ),
    }


def arb_collar_data() -> dict:
    flint.ctx.prec = PRECISION_BITS
    l0 = flint.arb(50)
    x0 = 4 * flint.arb.pi() * l0.exp()
    if not x0 > flint.arb(10) ** 22:
        raise RuntimeError("L=50 x-scale lower bound failed")

    quarter_power = (flint.arb(3).log() / 4).exp()
    l_guard = flint.arb(49)
    x_guard = 4 * flint.arb.pi() * l_guard.exp()
    n_lower = (l_guard / 2).exp() - flint.arb(9) / 8
    positive_correction = (
        (flint.arb(31) / 25)
        * (quarter_power + 1 / quarter_power)
        / n_lower
        + (
            3 * (l_guard**2 + flint.arb.pi() ** 2 / 4).sqrt()
            + flint.arb(261) / 25
        )
        / (x_guard - 12)
    )
    if not positive_correction < flint.arb(1) / 100:
        raise RuntimeError("published eC positive correction bound failed")

    u_bound = 2 * (((l0 + 1) ** 2) / 64 + flint.arb(313) / 500) / x0
    if not u_bound < 1:
        raise RuntimeError("published eAB exponential argument bound failed")
    raw_ab_bound = 8 * u_bound * (l0 / 16 + flint.arb(1) / 4).exp()
    raw_c_bound = (-3 * l0 / 16 + flint.arb(3) / 10).exp()
    if not raw_ab_bound < flint.arb(1) / 200000:
        raise RuntimeError("raw eAB collar bound failed")
    if not raw_c_bound < 3 * flint.arb(1) / 25000:
        raise RuntimeError("raw eC collar bound failed")
    if not raw_ab_bound + raw_c_bound < flint.arb(1) / 8000:
        raise RuntimeError("combined raw collar bound failed")
    return {
        "precision_bits": PRECISION_BITS,
        "x_at_L50_ball": arb_text(x0),
        "x_at_L50_gt": "10^22",
        "positive_eC_correction_ball": arb_text(positive_correction),
        "positive_eC_correction_lt": "1/100",
        "eAB_argument_ball": arb_text(u_bound),
        "raw_eAB_over_center_amplitude_ball": arb_text(raw_ab_bound),
        "raw_eAB_over_center_amplitude_lt": "1/200000",
        "raw_eC_over_center_amplitude_ball": arb_text(raw_c_bound),
        "raw_eC_over_center_amplitude_lt": "3/25000",
        "combined_raw_collar_ball": arb_text(raw_ab_bound + raw_c_bound),
        "combined_raw_collar_lt": "1/8000",
    }


def rational_budget() -> dict:
    delta = DELTA
    main0 = Fraction(282, 125)
    main1 = Fraction(59, 100)
    main2 = Fraction(77, 500)
    eps1 = Fraction(35, 100) * delta
    eps2 = Fraction(13, 100) * delta
    potential = Fraction(1, 10**21)
    l2_min = L_MIN**2
    error = (
        2 * main1 * eps1
        + eps1**2
        + main0 * eps2
        + delta * main2
        + delta * eps2
        + potential * (2 * main0 * delta + delta**2) / l2_min
    )
    error_cap = Fraction(1, 1000)
    if not error < error_cap:
        raise RuntimeError("normalized remainder curvature budget failed")
    finite_main_margin = Fraction(19, 250)
    final_margin = finite_main_margin - error_cap
    if final_margin != Fraction(3, 40):
        raise RuntimeError("unexpected full cell-interior margin")
    return {
        "collar": {
            "rho": fraction_text(RHO),
            "cell_boundary_distance": fraction_text(CELL_DISTANCE),
            "delta": fraction_text(delta),
        },
        "normalized_jet_caps": {
            "eps0": fraction_text(delta),
            "eps1_over_L": fraction_text(eps1),
            "eps2_over_L2": fraction_text(eps2),
        },
        "main_jet_caps": {
            "P": fraction_text(main0),
            "P1_over_L": fraction_text(main1),
            "P2_over_L2": fraction_text(main2),
        },
        "curvature": {
            "remainder_error_over_L2_raw": fraction_text(error),
            "remainder_error_over_L2": fraction_text(error_cap),
            "finite_main_margin_over_L2": fraction_text(finite_main_margin),
            "full_margin_over_L2": fraction_text(final_margin),
        },
    }


def build_artifact() -> dict:
    exact = build_exact()
    arb_data = arb_collar_data()
    budget = rational_budget()
    rows = [
        GateRow(
            id="np15dscic_01_fixed_n_remainder",
            role="published_theorem_transfer",
            readiness="ready_to_apply",
            claim="Inside one cutoff cell the published approximation defines a fixed holomorphic raw remainder.",
            formula=exact["raw_main"],
            proof_boundary="Imports the scalar error bounds from Polymath 15, Theorem 1.3.",
        ),
        GateRow(
            id="np15dscic_02_schwarz_collar",
            role="exact_lemma",
            readiness="ready_to_apply",
            claim="Reality on the real axis reflects the upper-half published bound across the full Cauchy circle.",
            formula=exact["reflection"],
            proof_boundary="Requires one fixed N throughout the disk.",
        ),
        GateRow(
            id="np15dscic_03_scalar_collar",
            role="interval_certificate",
            readiness="certified",
            claim="The raw theorem remainder is uniformly tiny relative to the center amplitude on every eligible collar.",
            formula=exact["scalar_collar"],
            proof_boundary="Uses the explicit published eA+eB+eC0 bounds and L>=50.",
            diagnostics=arb_data,
        ),
        GateRow(
            id="np15dscic_04_cauchy_raw_jets",
            role="exact_lemma",
            readiness="ready_to_apply",
            claim="Cauchy's estimate transfers the raw scalar collar through two derivatives.",
            formula=exact["cauchy_raw"],
            proof_boundary="Exact for disk radius 1/4.",
        ),
        GateRow(
            id="np15dscic_05_normalized_jets",
            role="rational_certificate",
            readiness="certified",
            claim="Dividing by the real normalizer preserves an explicit small C2 error.",
            formula=exact["normalized_jets"],
            proof_boundary="Uses |(log A)'|<=0.27L and |V|<=10^-21.",
            diagnostics=budget["normalized_jet_caps"],
        ),
        GateRow(
            id="np15dscic_06_main_jets",
            role="rational_certificate",
            readiness="certified",
            claim="The dominant-saddle arithmetic certificate also supplies uniform main-sum jet caps.",
            formula=exact["main_jets"],
            proof_boundary="Direct consequence of the prior infinite-zeta tail budgets.",
            diagnostics=budget["main_jet_caps"],
        ),
        GateRow(
            id="np15dscic_07_remainder_curvature",
            role="rational_certificate",
            readiness="certified",
            claim="The transferred theorem remainder consumes less than one thousandth of L-squared curvature.",
            formula=exact["remainder_curvature"],
            proof_boundary="Exact rational perturbation budget.",
            diagnostics=budget["curvature"],
        ),
        GateRow(
            id="np15dscic_08_exact_h_cell_theorem",
            role="proved_theorem",
            readiness="ready_to_apply",
            claim="The exact Newman heat-flow function has strict first Laguerre positivity on every dominant-saddle cutoff-cell interior.",
            formula=exact["theorem"],
            proof_boundary="Applies only under the stated tL and cell-distance hypotheses.",
        ),
        GateRow(
            id="np15dscic_09_cutoff_collar_gap",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Width-one neighborhoods of every discrete Riemann-Siegel cutoff transition remain outside this theorem.",
            formula="dist(x,cell boundary)<1/2",
            proof_boundary="Needs adjacent-N overlap or a direct transition-correction estimate.",
        ),
        GateRow(
            id="np15dscic_10_boundary_layer_gap",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The coupled near-zero-time layer remains the principal global obstruction.",
            formula="0<t*log(x/(4*pi))<25",
            proof_boundary="Open; not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate",
        "date": "2026-07-16",
        "status": (
            "strict Laguerre positivity for the exact Newman H_t on dominant-saddle "
            "fixed-cutoff cell interiors; not a proof of Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "This artifact proves L_t(x)>0 when 0<t<=1/2, "
            "t*log(x/(4*pi))>=25, and x is at least 1/2 from both boundaries "
            "of its Riemann-Siegel cutoff cell. It does not cover the transition "
            "collars, the layer t*log(x/(4*pi))<25, all real x and t, Lambda<=0, "
            "or RH."
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "rho": fraction_text(RHO),
            "cell_distance": fraction_text(CELL_DISTANCE),
            "L_min": fraction_text(L_MIN),
            "tL_min": fraction_text(TL_MIN),
            "delta": fraction_text(DELTA),
        },
        "exact": exact,
        "arb": arb_data,
        "rational_budget": budget,
        "rows": [asdict(row) for row in rows],
        "sources": [
            SOURCE_URL,
            "outputs/jensen_window_pf_newman_polymath15_normalized_laguerre_bridge.md",
            "outputs/jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    arb_data = artifact["arb"]
    budget = artifact["rational_budget"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Dominant-Saddle Cell-Interior Certificate",
            "",
            "Date: 2026-07-16",
            "",
            "Status: strict first Laguerre positivity for the exact Newman heat-flow",
            "function on dominant-saddle fixed-cutoff cell interiors. This is not a",
            "proof of `Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate.py",
            "```",
            "",
            "The scalar approximation is imported from D. H. J. Polymath,",
            f"[Effective approximation of heat flow evolution of the Riemann xi function]({SOURCE_URL}),",
            "Theorem 1.3. All derivative and curvature transfers below are explicit.",
            "",
            "## Region",
            "",
            "```text",
            exact["region"],
            exact["cell"],
            "```",
            "",
            "Since `t<=1/2`, the first condition forces `L>=50` and",
            f"`x>{arb_data['x_at_L50_ball']}`, in particular `x>10^22`.",
            "This deliberately large threshold buys a clean theorem; it does not",
            "indicate failure below it.",
            "",
            "## Raw Collar",
            "",
            "With `N` fixed on the radius-`1/4` disk, define",
            "",
            "```text",
            exact["raw_main"],
            exact["real_axis"],
            "```",
            "",
            "The main sum is holomorphic and real on the real axis. Schwarz",
            "reflection therefore transfers the upper-half Polymath-15 bound to",
            "the lower semicircle. The explicit published errors give",
            "",
            "```text",
            f"positive eC exponent correction < {arb_data['positive_eC_correction_lt']}",
            f"raw eA+eB / A(center) < {arb_data['raw_eAB_over_center_amplitude_lt']}",
            f"raw eC0 / A(center) < {arb_data['raw_eC_over_center_amplitude_lt']}",
            exact["scalar_collar"],
            "```",
            "",
            "The saved outward-rounded combined bound is",
            f"`{arb_data['combined_raw_collar_ball']}`.",
            "",
            "## Derivatives",
            "",
            "Cauchy's estimate first gives",
            "",
            "```text",
            exact["cauchy_raw"],
            "```",
            "",
            "For `ell=(log A)'` and `V=-(log A)''`, exact differentiation of",
            "`e=R/A=Z-P` gives",
            "",
            "```text",
            "e'=R'/A-ell*e",
            "e''=R''/A-2*ell*R'/A+(ell^2+V)*e",
            exact["normalized_jets"],
            "```",
            "",
            "The arithmetic ray theorem supplies the companion main-jet caps",
            "",
            "```text",
            exact["main_jets"],
            "```",
            "",
            "## Exact Heat Flow",
            "",
            "Substitution into the exact curvature perturbation identity gives",
            "the rational budget",
            "",
            "```text",
            f"finite-main margin / L^2 = {budget['curvature']['finite_main_margin_over_L2']}",
            f"remainder error / L^2  < {budget['curvature']['remainder_error_over_L2']}",
            f"exact-H margin / L^2   = {budget['curvature']['full_margin_over_L2']}",
            "```",
            "",
            "Therefore",
            "",
            "```text",
            exact["theorem"],
            "```",
            "",
            "## Remaining Gap",
            "",
            "The exact high-frequency theorem now fails only in two explicitly",
            "identified places: width-one neighborhoods around the discrete cutoff",
            "transitions, and the near-zero-time layer `0<tL<25`. The first is a",
            "local adjacent-`N` problem; the second retains the full RH-level",
            "arithmetic coupling.",
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
        "built Newman Polymath-15 dominant-saddle cell-interior certificate: "
        f"{len(artifact['rows'])} rows, exact-H margin "
        f"{artifact['rational_budget']['curvature']['full_margin_over_L2']}*L^2"
    )


if __name__ == "__main__":
    main()
