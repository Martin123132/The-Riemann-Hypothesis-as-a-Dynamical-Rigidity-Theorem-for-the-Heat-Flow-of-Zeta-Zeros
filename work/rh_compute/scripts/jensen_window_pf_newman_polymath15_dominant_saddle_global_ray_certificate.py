#!/usr/bin/env python3
"""Close cutoff transitions on the exact Newman dominant-saddle ray."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

from jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate import (  # noqa: E402
    DELTA,
    L_MIN,
    PRECISION_BITS,
    RHO,
    SOURCE_URL,
    arb_collar_data,
    rational_budget as cell_rational_budget,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate.md"
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


def fraction_text(value) -> str:
    return f"{value.numerator}/{value.denominator}"


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits)


def build_exact() -> dict:
    return {
        "region": (
            "0<t<=1/2, L=log(x/(4*pi)), t*L>=25, "
            "N_x=floor(sqrt(x/(4*pi)+t/16))"
        ),
        "boundary_spacing": (
            "x_(N+1)-x_N=4*pi*(2N+1); for L>=50 a radius-1/4 disk "
            "crosses at most one cutoff boundary"
        ),
        "fixed_center": (
            "Use S_(N_x,t)(z)=B_t(z)f_(N_x,t)(z) throughout the disk, even "
            "where the theorem-prescribed cutoff N_z differs from N_x"
        ),
        "adjacent_difference": (
            "|f_(N_z,t)(z)-f_(N_x,t)(z)|<=3*n^(-7/2), "
            "n>=exp(L/2)-1"
        ),
        "raw_jump": (
            "|S_(N_z,t)(z)-S_(N_x,t)(z)|/A_t(x)<10^-30"
        ),
        "global_collar": (
            "sup_(|z-x|=1/4)|H_t(z)-S_(N_x,t)(z)|/A_t(x)<1/8000"
        ),
        "theorem": (
            "For every real x on the ray, C_t[H_t/A_t](x)>3/40*L^2 and "
            "L_t(x)>0"
        ),
        "remaining_layer": "0<t*log(x/(4*pi))<25",
    }


def arb_transition_data() -> dict:
    flint.ctx.prec = PRECISION_BITS
    collar = arb_collar_data()
    l0 = flint.arb(50)
    n_lower = (l0 / 2).exp() - 1
    if not n_lower > flint.arb(10) ** 10:
        raise RuntimeError("adjacent cutoff index lower bound failed")
    raw_jump = (
        3
        * n_lower ** (-flint.arb(7) / 2)
        * (l0 / 16 + flint.arb(1) / 4).exp()
    )
    if not raw_jump < flint.arb(1) / flint.arb(10) ** 30:
        raise RuntimeError("adjacent cutoff raw jump bound failed")
    prior_raw = flint.arb(collar["combined_raw_collar_ball"])
    combined = prior_raw + raw_jump
    if not combined < flint.arb(1) / 8000:
        raise RuntimeError("global transition collar bound failed")
    spacing = 4 * flint.arb.pi() * (2 * n_lower + 1)
    if not spacing > flint.arb(1) / 2:
        raise RuntimeError("cutoff boundary spacing bound failed")
    return {
        "precision_bits": PRECISION_BITS,
        "adjacent_index_lower_ball": arb_text(n_lower),
        "adjacent_index_gt": "10^10",
        "cutoff_spacing_ball": arb_text(spacing),
        "cutoff_spacing_gt": "1/2",
        "raw_adjacent_jump_ball": arb_text(raw_jump),
        "raw_adjacent_jump_lt": "10^-30",
        "prior_fixed_cell_collar_ball": collar["combined_raw_collar_ball"],
        "global_collar_ball": arb_text(combined),
        "global_collar_lt": "1/8000",
    }


def build_artifact() -> dict:
    exact = build_exact()
    transition = arb_transition_data()
    budget = cell_rational_budget()
    rows = [
        GateRow(
            id="np15dsgrc_01_boundary_spacing",
            role="exact_geometry",
            readiness="ready_to_apply",
            claim="A Cauchy disk on the dominant ray can meet no more than one Riemann-Siegel cutoff transition.",
            formula=exact["boundary_spacing"],
            proof_boundary="Exact cutoff spacing plus the Arb lower bound on N.",
            diagnostics={
                "spacing_ball": transition["cutoff_spacing_ball"],
                "spacing_gt": transition["cutoff_spacing_gt"],
            },
        ),
        GateRow(
            id="np15dsgrc_02_fixed_center_cutoff",
            role="exact_handoff",
            readiness="ready_to_apply",
            claim="One may keep the center cutoff fixed as an analytic main sum throughout a transition-crossing disk.",
            formula=exact["fixed_center"],
            proof_boundary="The fixed finite sum remains holomorphic; only its comparison with the prescribed sum changes.",
        ),
        GateRow(
            id="np15dsgrc_03_adjacent_block",
            role="analytic_theorem",
            readiness="ready_to_apply",
            claim="The discrepancy between adjacent prescribed cutoffs is one strongly suppressed two-saddle block.",
            formula=exact["adjacent_difference"],
            proof_boundary="Uses tL>=25, the published Re(s_*) bound, and the gamma/kappa bounds.",
        ),
        GateRow(
            id="np15dsgrc_04_arb_jump",
            role="interval_certificate",
            readiness="certified",
            claim="Arb bounds the raw adjacent-cutoff jump far below the available collar slack.",
            formula=exact["raw_jump"],
            proof_boundary="Outward-rounded evaluation at the monotone endpoint L=50.",
            diagnostics=transition,
        ),
        GateRow(
            id="np15dsgrc_05_global_collar",
            role="interval_certificate",
            readiness="certified",
            claim="The same raw Cauchy collar bound now holds whether or not the disk crosses a cutoff.",
            formula=exact["global_collar"],
            proof_boundary="Combines the fixed-cell published remainder with at most one adjacent block.",
        ),
        GateRow(
            id="np15dsgrc_06_global_exact_h_ray",
            role="proved_theorem",
            readiness="ready_to_apply",
            claim="Strict first Laguerre positivity holds for the exact Newman heat flow on the entire dominant-saddle ray.",
            formula=exact["theorem"],
            proof_boundary="No cutoff-distance exception remains; the tL threshold remains essential.",
            diagnostics=budget["curvature"],
        ),
        GateRow(
            id="np15dsgrc_07_cutoff_gap_closed",
            role="closure_record",
            readiness="closed",
            claim="The discrete cutoff-transition collars are no longer an open high-frequency gap.",
            formula="adjacent-N correction absorbed into delta=1/8000",
            proof_boundary="Closure only on tL>=25.",
        ),
        GateRow(
            id="np15dsgrc_08_live_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The sole remaining global Newman-tail regime is the shrinking near-zero-time arithmetic layer.",
            formula=exact["remaining_layer"],
            proof_boundary="Open; not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate",
        "date": "2026-07-17",
        "status": (
            "strict first Laguerre positivity for exact H_t on the complete ray "
            "t*log(x/(4*pi))>=25, including all cutoff transitions; not a proof "
            "of Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "This artifact proves L_t(x)>0 for every real x satisfying "
            "0<t<=1/2 and t*log(x/(4*pi))>=25, with no cutoff-distance "
            "restriction. It does not cover 0<t*log(x/(4*pi))<25, all real x "
            "and t, Lambda<=0, or RH."
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "rho": fraction_text(RHO),
            "L_min": fraction_text(L_MIN),
            "delta": fraction_text(DELTA),
            "tL_min": "25/1",
        },
        "exact": exact,
        "arb": transition,
        "rational_budget": budget,
        "rows": [asdict(row) for row in rows],
        "sources": [
            SOURCE_URL,
            "outputs/jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate.md",
            "outputs/jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    arb_data = artifact["arb"]
    curvature = artifact["rational_budget"]["curvature"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Dominant-Saddle Global Ray Certificate",
            "",
            "Date: 2026-07-17",
            "",
            "Status: strict first Laguerre positivity for the exact Newman heat-flow",
            "function on the complete dominant-saddle ray, including every cutoff",
            "transition. This is not a proof of `Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate.py",
            "```",
            "",
            "The effective approximation is imported from D. H. J. Polymath,",
            f"[Effective approximation of heat flow evolution of the Riemann xi function]({SOURCE_URL}),",
            "Theorem 1.3.",
            "",
            "## Region",
            "",
            "```text",
            exact["region"],
            "```",
            "",
            "The prior certificate proved the exact theorem whenever a radius-`1/4`",
            "disk stayed inside one cutoff cell. This certificate removes that",
            "restriction.",
            "",
            "## One-Block Repair",
            "",
            "Consecutive boundaries have spacing",
            "",
            "```text",
            exact["boundary_spacing"],
            "```",
            "",
            "At `L=50` the relevant adjacent index is already enclosed below by",
            f"`{arb_data['adjacent_index_lower_ball']}`. Thus a radius-`1/4` disk",
            "crosses at most one boundary. Keep the center cutoff fixed throughout",
            "the disk. Where the theorem prescribes its neighbor, the two main sums",
            "differ by one block and",
            "",
            "```text",
            exact["adjacent_difference"],
            "```",
            "",
            "After restoring the raw normalizer, Arb gives",
            "",
            "```text",
            f"adjacent raw jump / A(center) = {arb_data['raw_adjacent_jump_ball']}",
            exact["raw_jump"],
            "```",
            "",
            "This is negligible beside the existing fixed-cell collar. Their",
            "outward-rounded sum is",
            "",
            "```text",
            f"{arb_data['global_collar_ball']} < 1/8000",
            "```",
            "",
            "so all previous Cauchy, normalized-jet, and curvature budgets apply",
            "unchanged at the transition itself.",
            "",
            "## Global Ray",
            "",
            "The exact rational margins remain",
            "",
            "```text",
            f"finite-main margin / L^2 = {curvature['finite_main_margin_over_L2']}",
            f"remainder error / L^2  < {curvature['remainder_error_over_L2']}",
            f"exact-H margin / L^2   = {curvature['full_margin_over_L2']}",
            "```",
            "",
            "Therefore, with no cutoff exception,",
            "",
            "```text",
            exact["theorem"],
            "```",
            "",
            "## Remaining Gap",
            "",
            "The high-frequency dominant-saddle branch is now closed. The remaining",
            "global target is exactly the nonuniform layer",
            "",
            "```text",
            exact["remaining_layer"],
            "```",
            "",
            "where the finite Dirichlet terms no longer form a small perturbation",
            "of the first saddle and the RH-level arithmetic cancellation must be",
            "retained.",
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
        "built Newman Polymath-15 dominant-saddle global ray certificate: "
        f"{len(artifact['rows'])} rows, cutoff gap closed, exact-H margin 3/40*L^2"
    )


if __name__ == "__main__":
    main()
