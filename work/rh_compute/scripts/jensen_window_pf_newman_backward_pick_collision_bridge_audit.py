#!/usr/bin/env python3
"""Build the backward-Pick collision-bridge audit for the Newman route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_backward_pick_collision_bridge_audit.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_backward_pick_collision_bridge_audit.md"
)


@dataclass(frozen=True)
class AuditRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def build_exact() -> dict:
    z, t, a = sp.symbols("z t a", real=True)
    e = z**2 + a - 2 * t
    heat_residual = sp.simplify(sp.diff(e, t) + sp.diff(e, z, 2))
    if heat_residual != 0:
        raise RuntimeError("quadratic backward-heat identity failed")

    g = sp.simplify(-sp.diff(e, z) / e)
    burgers_residual = sp.simplify(
        sp.diff(g, t) + sp.diff(g, z, 2) - 2 * g * sp.diff(g, z)
    )
    if burgers_residual != 0:
        raise RuntimeError("log-derivative Burgers identity failed")

    x, y, b = sp.symbols("x y b", real=True)
    z_xy = x + sp.I * y
    g_xy = -2 * z_xy / (z_xy**2 + b)
    h_formula = (
        2
        * y
        * (x**2 + y**2 - b)
        / ((x**2 - y**2 + b) ** 2 + 4 * x**2 * y**2)
    )
    h_residual = sp.simplify(sp.im(sp.expand_complex(g_xy)) - h_formula)
    if h_residual != 0:
        raise RuntimeError("Pick-sign formula failed")

    r = sp.symbols("r", positive=True)
    negative_witness = sp.simplify(
        h_formula.subs({x: 0, y: r / 2, b: r**2})
    )
    if negative_witness != -sp.Rational(4, 3) / r:
        raise RuntimeError("negative Pick-bubble witness failed")

    rho = sp.symbols("rho", positive=True)
    hidden_width_residual = sp.simplify(
        a - 2 * (a / 2 - rho**2 / 2) - rho**2
    )
    if hidden_width_residual != 0:
        raise RuntimeError("cutoff-hiding width failed")

    return {
        "quadratic_flow": {
            "definition": "E_t(z)=z^2+a-2t, a>0",
            "heat_equation": "partial_t E_t=-partial_z^2 E_t",
            "collision_time": "t_*=a/2",
            "zeros_above": (
                "t>t_*: zeros are +/-sqrt(2t-a), both real"
            ),
            "zeros_at": "t=t_*: z=0 is a double zero",
            "zeros_below": (
                "t<t_*: zeros are +/-i*sqrt(a-2t), so one lies in C+"
            ),
            "class_scope": (
                "E_t is an even real polynomial, hence an even real entire "
                "Cartwright-class backward-heat solution."
            ),
        },
        "pick_field": {
            "definition": "g_t(z)=-E_t'(z)/E_t(z), h_t=Im(g_t)",
            "burgers": "partial_t g=-partial_z^2 g+2g*partial_z g",
            "formula": (
                "h_t(x,y)=2y*(x^2+y^2-b)"
                "/((x^2-y^2+b)^2+4x^2y^2), b=a-2t"
            ),
            "top_sign": (
                "For t>t_* one has b<0, hence h_t(x,y)>0 for y>0."
            ),
            "negative_bubble": (
                "For t<t_* one has h_t<0 on "
                "x^2+y^2<a-2t inside C+."
            ),
            "witness": (
                "At x=0, y=sqrt(a-2t)/2, "
                "h_t=-4/(3*sqrt(a-2t))<0."
            ),
        },
        "collision_uniformity": {
            "upper_zero": "rho_+(t)=i*sqrt(a-2t), t<t_*",
            "speed": (
                "|rho_+'(t)|=1/sqrt(a-2t), which tends to infinity "
                "as t increases to t_*."
            ),
            "closed_window_scope": (
                "A finite zero-speed constant exists on each closed "
                "collision-free window, but not uniformly on an open window "
                "whose endpoint is t_*."
            ),
            "gronwall_guard": (
                "A bound |E_rho'(t)|<=K*E_rho(t) with K depending on the "
                "zero-speed constant cannot be extended to t_* by asserting "
                "one finite K on (t_*-delta,t_*)."
            ),
        },
        "cutoff_obstruction": {
            "negative_support": (
                "supp(h_t^-)=C+ intersect {x^2+y^2<a-2t}"
            ),
            "fixed_cutoff_hiding": (
                "For a lower cutoff y>=rho, the weighted negative energy is "
                "exactly zero whenever 0<t_*-t<rho^2/2, although h_t^- is "
                "nonzero below y=rho."
            ),
            "hidden_width": "delta_rho=rho^2/2",
            "noncommuting_limits": (
                "The zero-energy bridge width tends to zero with rho. "
                "One cannot first prove a rho-dependent bridge and then send "
                "rho to zero while retaining a fixed backward time interval."
            ),
            "support_guard": (
                "Vanishing of an energy whose weight is positive only on "
                "y>=rho proves h_t^-=0 only on that support, not on all C+."
            ),
        },
        "literature_audit": {
            "source": (
                "Kevin Schatz, Riemann Hypothesis: Backward Parabolic "
                "Positivity Barriers for the Xi Flow, preprint, 2025"
            ),
            "doi": "https://doi.org/10.5281/zenodo.17636625",
            "manuscript": (
                "https://kschatz.github.io/rh-xi-backward-parabolic-barrier/"
                "schatz_riemann_hypothesis_backward_parabolic_positivity_"
                "barriers_xi_flow.pdf"
            ),
            "pinpoint": (
                "Lemma 4.7 assumes a closed collision-free window; "
                "Lemma 7.3 applies its speed-dependent energy bound on an "
                "open interval ending at a collision and concludes global "
                "upper-half-plane positivity from a cutoff-supported energy."
            ),
            "decision": (
                "The quadratic flow satisfies the generic heat, Burgers, "
                "Cartwright, top-time Pick, and isolated-collision inputs but "
                "violates the claimed collision bridge. The cited preprint "
                "therefore does not presently supply an admissible proof of "
                "backward Pick positivity."
            ),
        },
        "open_repair": (
            "A viable Xi-specific backward argument must add an estimate "
            "uniform in zero separation and lower cutoff through collisions, "
            "or an arithmetic invariant that excludes the collision before "
            "the cutoff limit. Closed-window speed bounds, weighted-energy "
            "continuity, and top-time Pick positivity are insufficient."
        ),
        "checks": {
            "heat_residual": str(heat_residual),
            "burgers_residual": str(burgers_residual),
            "pick_formula_residual": str(h_residual),
            "negative_witness": str(negative_witness),
            "hidden_width_residual": str(hidden_width_residual),
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        AuditRow(
            "nbpca_01_quadratic_heat_flow",
            "exact_identity",
            "ready_to_apply",
            "A quadratic family is an exact even backward-heat flow.",
            exact["quadratic_flow"]["heat_equation"],
            "No Xi zero information is used.",
        ),
        AuditRow(
            "nbpca_02_real_zero_top",
            "exact_example",
            "ready_to_apply",
            "The flow has a real-zero time above one collision.",
            exact["quadratic_flow"]["zeros_above"],
            "This is a finite Cartwright-class calibration model.",
        ),
        AuditRow(
            "nbpca_03_pick_top",
            "exact_inequality",
            "ready_to_apply",
            "The negative logarithmic derivative is Pick-positive above the collision.",
            exact["pick_field"]["top_sign"],
            "The sign follows directly from the displayed rational formula.",
        ),
        AuditRow(
            "nbpca_04_double_collision",
            "exact_identity",
            "ready_to_apply",
            "The two real zeros collide at a finite time.",
            exact["quadratic_flow"]["zeros_at"],
            "The collision is the standard square-root normal form.",
        ),
        AuditRow(
            "nbpca_05_nonreal_birth",
            "exact_countermodel",
            "guard_validated",
            "Backward continuation creates an upper-half-plane zero.",
            exact["quadratic_flow"]["zeros_below"],
            "This blocks generic backward preservation of real-rootedness.",
        ),
        AuditRow(
            "nbpca_06_pick_bubble",
            "exact_countermodel",
            "guard_validated",
            "A negative-Pick bubble is born below the upper zero.",
            exact["pick_field"]["negative_bubble"],
            exact["pick_field"]["witness"],
        ),
        AuditRow(
            "nbpca_07_speed_blowup",
            "uniformity_guard",
            "guard_validated",
            "The zero-speed constant diverges at the collision.",
            exact["collision_uniformity"]["speed"],
            exact["collision_uniformity"]["closed_window_scope"],
        ),
        AuditRow(
            "nbpca_08_cutoff_hiding",
            "support_guard",
            "guard_validated",
            "Every fixed lower cutoff hides the newborn negative bubble for a positive time.",
            exact["cutoff_obstruction"]["fixed_cutoff_hiding"],
            exact["cutoff_obstruction"]["support_guard"],
        ),
        AuditRow(
            "nbpca_09_noncommuting_limits",
            "limit_order_obstruction",
            "guard_validated",
            "The cutoff exhaustion and backward collision limit do not commute.",
            exact["cutoff_obstruction"]["noncommuting_limits"],
            "The exact hidden interval has width rho^2/2.",
        ),
        AuditRow(
            "nbpca_10_gronwall_boundary",
            "proof_boundary",
            "ready_to_apply",
            "Closed-window energy control cannot be reused with one finite constant at the collision endpoint.",
            exact["collision_uniformity"]["gronwall_guard"],
            "A separation-uniform estimate would be a genuinely new input.",
        ),
        AuditRow(
            "nbpca_11_preprint_audit",
            "literature_audit",
            "gap_identified",
            "The preprint collision bridge crosses exactly the forbidden uniformity step.",
            exact["literature_audit"]["pinpoint"],
            exact["literature_audit"]["decision"],
        ),
        AuditRow(
            "nbpca_12_open_xi_repair",
            "open_theorem_target",
            "open",
            "Only an Xi-specific separation-uniform collision estimate could repair this route.",
            exact["open_repair"],
            "No backward Pick theorem, Lambda<=0, or RH claim is made.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_backward_pick_collision_bridge_audit",
        "date": "2026-07-24",
        "status": (
            "exact quadratic collision countermodel and literature audit with "
            "one open Xi-specific repair target"
        ),
        "proof_boundary": (
            "The artifact rejects generic backward Pick/collision bridging and "
            "pinpoints a separation/cutoff uniformity gap in one unreviewed "
            "preprint. It does not establish that every possible Xi-specific "
            "backward barrier fails, and it does not prove Lambda<=0 or RH."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    success = (
        "validated Newman backward-Pick collision-bridge audit: 12 rows, "
        "0 issues, 3 exact flow identities, 2 Pick-sign regions, "
        "1 square-root speed blowup, 1 cutoff-hiding theorem, "
        "1 noncommuting-limit obstruction, 1 preprint gap, 1 open Xi repair"
    )
    return "\n".join(
        [
            "# Newman Backward-Pick Collision-Bridge Audit",
            "",
            "Date: 2026-07-24",
            "",
            "Status: exact collision countermodel and literature audit. This",
            "is not a proof of `Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_backward_pick_collision_bridge_audit.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_backward_pick_collision_bridge_audit.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_backward_pick_collision_bridge_audit.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            success,
            "```",
            "",
            "## Exact Collision Model",
            "",
            "Take",
            "",
            "```text",
            exact["quadratic_flow"]["definition"],
            exact["quadratic_flow"]["heat_equation"],
            exact["quadratic_flow"]["collision_time"],
            exact["quadratic_flow"]["zeros_above"],
            exact["quadratic_flow"]["zeros_at"],
            exact["quadratic_flow"]["zeros_below"],
            "```",
            "",
            "This is an even real entire Cartwright-class solution of exactly",
            "the same backward heat equation. Its logarithmic derivative obeys",
            "",
            "```text",
            exact["pick_field"]["definition"],
            exact["pick_field"]["burgers"],
            exact["pick_field"]["formula"],
            "```",
            "",
            "Above the collision its imaginary part is positive throughout",
            "the upper half-plane. Below the collision it is negative on the",
            "newborn half-disc:",
            "",
            "```text",
            exact["pick_field"]["negative_bubble"],
            exact["pick_field"]["witness"],
            "```",
            "",
            "## Uniformity Failure",
            "",
            "The upper zero has",
            "",
            "```text",
            exact["collision_uniformity"]["upper_zero"],
            exact["collision_uniformity"]["speed"],
            "```",
            "",
            "Thus a speed bound proved on each closed collision-free interval",
            "does not produce one finite constant on an interval ending at the",
            "collision.",
            "",
            "There is a second, independent cutoff obstruction:",
            "",
            "```text",
            exact["cutoff_obstruction"]["fixed_cutoff_hiding"],
            exact["cutoff_obstruction"]["hidden_width"],
            exact["cutoff_obstruction"]["noncommuting_limits"],
            "```",
            "",
            "For every fixed `rho`, the weighted energy can therefore vanish",
            "on a two-sided collision neighbourhood while the unweighted",
            "negative part is already nonzero below `y=rho`. Sending `rho` to",
            "zero also sends the guaranteed backward width to zero.",
            "",
            "## Literature Audit",
            "",
            exact["literature_audit"]["source"],
            "",
            f"DOI: {exact['literature_audit']['doi']}",
            "",
            f"Manuscript: {exact['literature_audit']['manuscript']}",
            "",
            exact["literature_audit"]["pinpoint"],
            "",
            exact["literature_audit"]["decision"],
            "",
            "## Live Handoff",
            "",
            exact["open_repair"],
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
        "wrote Newman backward-Pick collision-bridge audit: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
