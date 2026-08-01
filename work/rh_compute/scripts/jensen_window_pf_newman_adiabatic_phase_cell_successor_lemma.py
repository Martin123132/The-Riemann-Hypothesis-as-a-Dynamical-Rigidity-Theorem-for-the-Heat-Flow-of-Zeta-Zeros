#!/usr/bin/env python3
"""Build the exact adiabatic phase-cell successor lemma."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{STEM}.md"
DIAGNOSTIC_PATH = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_cofinal_phase_cell_scaling_diagnostics.json"
)
Q208_PATH = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_q208_closed_boundary_winding_certificate.json"
)
BOUNDARY_TARGET_PATH = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.json"
)
DATE = "2026-07-25"


@dataclass(frozen=True)
class LemmaRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def exact_derivations() -> dict:
    j = sp.symbols("j", integer=True, positive=True)
    t_j = sp.Rational(1, 5) / j
    t_next = sp.Rational(1, 5) / (j + 1)
    delta = sp.factor(t_j - t_next)
    expected_delta = 1 / (5 * j * (j + 1))
    if sp.simplify(delta - expected_delta) != 0:
        raise RuntimeError("successor time step identity failed")
    right = j + 38
    right_next = j + 39
    if sp.simplify(right_next - right) != 1:
        raise RuntimeError("successor right step identity failed")
    if sp.simplify(t_j * right - (sp.Rational(1, 5) + 38 / (5 * j))) != 0:
        raise RuntimeError("bottom hyperbola identity failed")

    h_xx, h_xxx, ell = sp.symbols(
        "h_xx h_xxx ell",
        real=True,
        positive=True,
    )
    heat_jet = sp.Matrix([-h_xx, -h_xxx / ell])
    heat_jet_norm_squared = sp.expand(heat_jet.dot(heat_jet))
    if heat_jet_norm_squared != h_xx**2 + h_xxx**2 / ell**2:
        raise RuntimeError("heat-jet transport identity failed")

    margin, transport, s = sp.symbols(
        "margin transport s",
        positive=True,
    )
    reverse_triangle_floor = sp.factor(margin - s * transport)
    if sp.simplify(
        reverse_triangle_floor.subs(s, 1) - (margin - transport)
    ) != 0:
        raise RuntimeError("reverse-triangle homotopy identity failed")

    return {
        "shell": {
            "low_rectangle": "P_j=[t_j,1/5]x[0,R_j]",
            "parameters": "t_j=1/(5*j), R_j=j+38",
            "bottom_collar": (
                "C_j=[t_(j+1),t_j]x[0,R_j]"
            ),
            "right_strip": (
                "S_j=[t_(j+1),1/5]x[R_j,R_(j+1)]"
            ),
            "decomposition": "P_(j+1)=P_j union C_j union S_j",
            "time_step": str(delta),
            "right_step": str(sp.factor(right_next - right)),
            "time_times_right": str(sp.factor(t_j * right)),
        },
        "heat_jet": {
            "heat_equation": "partial_t H=-partial_x^2 H",
            "scaled_jet": "V_ell=(H,H_x/ell(x)), ell(x)>0",
            "time_derivative": (
                "partial_t V_ell=(-H_xx,-H_xxx/ell)"
            ),
            "time_derivative_norm_squared": str(
                heat_jet_norm_squared
            ),
            "scale_guard": (
                "ell is independent of t on each collar; a t-dependent "
                "scale requires its extra derivative term."
            ),
        },
        "transport": {
            "pointwise_bound": (
                "||V_ell(t,x)-V_ell(t_j,x)||_2 "
                "<=(t_j-t)*M_(j,k)"
            ),
            "panel_bound": (
                "M_(j,k)>=sup_(C_j over I_(j,k)) "
                "sqrt(H_xx^2+(H_xxx/ell)^2)"
            ),
            "strict_gate": (
                "delta_j*M_(j,k)<d_(j,k)"
            ),
            "homotopy_floor": (
                "||V_ell(t,x)||_2>=d_(j,k)-"
                "delta_j*M_(j,k)>0"
            ),
            "symbolic_floor": str(reverse_triangle_floor),
        },
        "phase_cells": {
            "reference_cell": (
                "V_ell(t_j,I_(j,k)) subset C_(j,k), "
                "where C_(j,k) is compact and convex"
            ),
            "distance": "d_(j,k)=dist(0,C_(j,k))>0",
            "transported_cell": (
                "V_ell(C_j over I_(j,k)) subset "
                "C_(j,k)+closed_ball(0,delta_j*M_(j,k))"
            ),
            "origin_avoidance": (
                "dist(0,C_(j,k)+closed_ball(0,rho))"
                ">=d_(j,k)-rho>0"
            ),
            "consequence": (
                "The actual heat-flow path supplies an origin-free "
                "homotopy between the two bottom paths."
            ),
        },
        "right_half_plane": {
            "gate": (
                "There is a fixed nonzero q_j in R^2 and eta_j>0 "
                "such that q_j dot V_ell(t,x)>=eta_j on S_j."
            ),
            "consequence": (
                "V_ell is nonzero on S_j and the complete strip image "
                "lies in one open half-plane."
            ),
            "q208_calibration": (
                "For the Q208 transformed right strip one may take "
                "q=(0,1), because F_x>0 in all 20 stored cells."
            ),
        },
    }


def build_payload() -> dict:
    exact = exact_derivations()
    diagnostics = json.loads(DIAGNOSTIC_PATH.read_text(encoding="utf-8"))
    q208 = json.loads(Q208_PATH.read_text(encoding="utf-8"))
    boundary = json.loads(
        BOUNDARY_TARGET_PATH.read_text(encoding="utf-8")
    )
    comparison = diagnostics["comparisons"]["q207_q208"]
    turning = {
        row["scale"]: row
        for row in diagnostics["witness_turning"]["bottom"]
    }
    rows = [
        LemmaRow(
            id="napcs_01_shell_decomposition",
            role="exact_geometry",
            readiness="proved_exact",
            claim=(
                "Each successor low rectangle is the old rectangle "
                "plus one global bottom collar and one local right strip."
            ),
            formula=exact["shell"]["decomposition"],
            proof_boundary=(
                "This identity does not show either added region is "
                "contact-free."
            ),
            diagnostics=exact["shell"],
        ),
        LemmaRow(
            id="napcs_02_heat_jet_transport",
            role="exact_identity",
            readiness="proved_exact",
            claim=(
                "For a time-independent positive derivative scale, the "
                "time derivative of the scaled first jet is controlled "
                "exactly by the second and third spatial derivatives."
            ),
            formula=exact["heat_jet"]["time_derivative"],
            proof_boundary=exact["heat_jet"]["scale_guard"],
            diagnostics=exact["heat_jet"],
        ),
        LemmaRow(
            id="napcs_03_reverse_triangle_collar",
            role="exact_lemma",
            readiness="proved_exact",
            claim=(
                "A reference first-jet margin larger than the integrated "
                "time displacement excludes the origin throughout a "
                "bottom collar."
            ),
            formula=exact["transport"]["homotopy_floor"],
            proof_boundary=(
                "Both the reference margin and the collar derivative "
                "upper bound must be proved for Xi."
            ),
            diagnostics=exact["transport"],
        ),
        LemmaRow(
            id="napcs_04_phase_cell_thickening",
            role="exact_lemma",
            readiness="proved_exact",
            claim=(
                "A compact convex origin-free phase cell remains "
                "origin-free after any Euclidean thickening smaller "
                "than its distance from the origin."
            ),
            formula=exact["phase_cells"]["origin_avoidance"],
            proof_boundary=(
                "The Q208 cells supply finite reference examples, not "
                "all-j distances or transport bounds."
            ),
            diagnostics=exact["phase_cells"],
        ),
        LemmaRow(
            id="napcs_05_right_half_plane",
            role="exact_lemma",
            readiness="proved_exact",
            claim=(
                "A fixed strict half-plane functional on the new right "
                "strip excludes contacts there and contributes zero "
                "local strip degree."
            ),
            formula=exact["right_half_plane"]["gate"],
            proof_boundary=(
                "Q208 has derivative positivity; an all-j half-plane "
                "theorem remains open."
            ),
            diagnostics=exact["right_half_plane"],
        ),
        LemmaRow(
            id="napcs_06_successor_theorem",
            role="exact_composition",
            readiness="proved_conditional",
            claim=(
                "If P_j is contact-free, every bottom phase cell passes "
                "the strict transport gate, and the new right strip "
                "passes the half-plane gate, then P_(j+1) is "
                "contact-free and has the same zero first-jet degree."
            ),
            formula=(
                "P_j contact-free and "
                "forall k: delta_j*M_(j,k)<d_(j,k) and "
                "q_j dot V_ell>=eta_j>0 on S_j "
                "imply P_(j+1) contact-free"
            ),
            proof_boundary=(
                "This is an induction interface. It does not prove the "
                "Xi inequalities in its antecedent."
            ),
        ),
        LemmaRow(
            id="napcs_07_cofinal_induction",
            role="exact_composition",
            readiness="proved_conditional",
            claim=(
                "A finite base and the successor hypotheses for every "
                "later j certify every low rectangle; the published "
                "positive-time cap then yields Lambda<=0 and RH."
            ),
            formula=(
                "P_208 plus successor(j) for all j>=208 "
                "imply union_j P_j=(0,1/5]x[0,infinity)"
            ),
            proof_boundary=(
                "The all-j successor hypotheses are the open RH-level "
                "arithmetic theorem."
            ),
        ),
        LemmaRow(
            id="napcs_08_q208_calibration",
            role="finite_calibration",
            readiness="proved_finite_only",
            claim=(
                "Q208 supplies a rigorous finite base, 20 derivative-"
                "positive right-strip cells, and stable half-unit phase "
                "geometry for calibrating the successor interface."
            ),
            formula=(
                "Q1..Q208 certified; Q208 boundary has 1005 cells and "
                "winding zero"
            ),
            proof_boundary=(
                "Q208 was proved directly. The adiabatic Q207-to-Q208 "
                "transport inequality has not been certified."
            ),
            diagnostics={
                "q207_q208_branch_agreement": (
                    f"{comparison['branch_agreements']}/"
                    f"{comparison['comparable_panels']}"
                ),
                "bottom_unit_max_panel_turn": turning["unit"][
                    "maximum_absolute_panel_turn"
                ],
                "bottom_half_log_max_panel_turn": turning["half_log"][
                    "maximum_absolute_panel_turn"
                ],
                "q208_conclusion": q208["conclusion"],
            },
        ),
        LemmaRow(
            id="napcs_09_asymptotic_budget",
            role="open_handoff",
            readiness="open",
            claim=(
                "The live analytic task is a panelwise lower bound for "
                "the scaled first jet and an upper bound for its "
                "second/third-derivative transport whose ratio beats "
                "delta_j=1/(5*j*(j+1))."
            ),
            formula=(
                "sup_k delta_j*M_(j,k)/d_(j,k)<1, together with a "
                "uniform right-strip half-plane cone"
            ),
            proof_boundary=(
                "The dominant-saddle t*L>=25 theorem cannot supply this "
                "bottom estimate because t_j*L_j tends to zero."
            ),
            diagnostics={
                "selected_route": diagnostics["route_decision"][
                    "selected"
                ],
                "next_exact_obligation": diagnostics["route_decision"][
                    "next_exact_obligation"
                ],
            },
        ),
        LemmaRow(
            id="napcs_10_nonpromotion",
            role="nonpromotion_gate",
            readiness="active",
            claim=(
                "Finite branch agreement, small panel turning, large "
                "tail domination, and Q208 do not establish the "
                "adiabatic or right-strip hypotheses for all j."
            ),
            formula=(
                "finite diagnostics != all-j Xi transport theorem"
            ),
            proof_boundary=(
                "Q209, every later successor, the cofinal theorem, "
                "Lambda<=0, RH, and the Clay prize remain open."
            ),
        ),
    ]
    return {
        "kind": STEM,
        "date": DATE,
        "status": (
            "exact conditional adiabatic phase-cell successor lemma "
            "with finite Q208 calibration; Xi all-j antecedent open"
        ),
        "source_sha256": {
            "scaling_diagnostics": file_hash(DIAGNOSTIC_PATH),
            "q208_closed_boundary": file_hash(Q208_PATH),
            "cofinal_boundary_target": file_hash(
                BOUNDARY_TARGET_PATH
            ),
        },
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "conditional_endgame": (
            boundary["exact"]["cofinal_contract"][
                "conditional_endgame"
            ]
        ),
        "next_exact_obligation": (
            "Construct a rigorous collar certifier for "
            "sqrt(H_xx^2+(H_xxx/ell)^2), first on the Q207-to-Q208 "
            "collar as a calibration, and compare "
            "delta_j*M_(j,k) with exact phase-cell distances. Then seek "
            "an analytic all-j bound before using Q209 as anything "
            "beyond a falsification case."
        ),
        "proof_boundary": (
            "The shell decomposition, heat-jet transport identity, "
            "cell-thickening criterion, half-plane criterion, and "
            "successor induction are exact. No Q207-to-Q208 transport "
            "budget, no all-j bottom-collar estimate, no all-j "
            "right-strip cone, no cofinal Xi theorem, no Lambda<=0, "
            "and no RH proof is supplied."
        ),
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    return "\n".join(
        [
            "# Newman Adiabatic Phase-Cell Successor Lemma",
            "",
            f"Date: {DATE}",
            "",
            "Status: exact conditional successor theorem with finite Q208",
            "calibration. The Xi all-j antecedent remains open; this is",
            "not a proof of `Lambda<=0` or RH.",
            "",
            "## Shell Decomposition",
            "",
            "```text",
            exact["shell"]["low_rectangle"],
            exact["shell"]["bottom_collar"],
            exact["shell"]["right_strip"],
            exact["shell"]["decomposition"],
            f"delta_j={exact['shell']['time_step']}",
            "```",
            "",
            "A successor is not merely a new corner: lowering the bottom",
            "time adds a collar across the complete old x-range. The",
            "one-unit right strip is the only spatially local addition.",
            "",
            "## Heat-Jet Transport",
            "",
            "For a positive scale `ell(x)` independent of time,",
            "",
            "```text",
            exact["heat_jet"]["scaled_jet"],
            exact["heat_jet"]["time_derivative"],
            (
                "||partial_t V_ell||_2^2="
                f"{exact['heat_jet']['time_derivative_norm_squared']}"
            ),
            "```",
            "",
            "On a phase cell `I_(j,k)`, prove",
            "",
            "```text",
            "V_ell(t_j,I_(j,k)) subset C_(j,k)",
            "d_(j,k)=dist(0,C_(j,k))>0",
            "M_(j,k)>=sup sqrt(H_xx^2+(H_xxx/ell)^2) on the collar",
            "delta_j*M_(j,k)<d_(j,k)",
            "```",
            "",
            "Then reverse triangle and convex-cell thickening give",
            "",
            "```text",
            "||V_ell(t,x)||_2>=d_(j,k)-delta_j*M_(j,k)>0.",
            "```",
            "",
            "Thus the actual heat flow is an origin-free homotopy between",
            "the old and new bottom paths.",
            "",
            "## Right Strip",
            "",
            "It is enough to find a fixed vector `q_j` and `eta_j>0`",
            "with",
            "",
            "```text",
            "q_j dot V_ell(t,x)>=eta_j on S_j.",
            "```",
            "",
            "The strip then lies in an open half-plane and contains no",
            "contact. Q208 calibrates this with `q=(0,1)` and `F_x>0`",
            "on all 20 stored right-strip cells.",
            "",
            "## Successor Theorem",
            "",
            "```text",
            "P_j contact-free",
            "and every bottom cell satisfies delta_j*M_(j,k)<d_(j,k)",
            "and the new right strip satisfies one strict half-plane gate",
            "imply P_(j+1) contact-free with zero first-jet degree.",
            "```",
            "",
            "A finite base plus these hypotheses for every later stage",
            "would certify the cofinal exhaustion and hence `Lambda<=0`.",
            "",
            "## Live Obligation",
            "",
            payload["next_exact_obligation"],
            "",
            "## Proof Boundary",
            "",
            payload["proof_boundary"],
            "",
        ]
    )


def main() -> None:
    payload = build_payload()
    RESULT_PATH.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    NOTE_PATH.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Newman adiabatic phase-cell successor lemma: "
        f"{len(payload['rows'])} rows, Xi all-j antecedent open"
    )


if __name__ == "__main__":
    main()
