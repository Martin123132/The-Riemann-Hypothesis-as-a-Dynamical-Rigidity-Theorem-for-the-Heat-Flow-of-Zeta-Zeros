#!/usr/bin/env python3
"""Build the ray-aligned energy cell-localization and ownership gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_ray_aligned_energy_cell_localization_gate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "energy_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_parabolic_frequency_energy_current_gate.json"
    ),
    "endpoint_margin": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_ray_aligned_energy_endpoint_margin_gate.json"
    ),
    "dominant_global": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate.json"
    ),
    "oscillatory_handoff": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.json"
    ),
    "critical_c1": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.json"
    ),
    "successor_lemma": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma.json"
    ),
    "finite_transport": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_q207_q208_forward_adiabatic_successor_certificate.json"
    ),
    "q208_bottom": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_q208_bottom_phase_cell_certificate.json"
    ),
    "q208_top": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_q208_top_phase_cell_certificate.json"
    ),
}


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }


def audit_sources(payloads: dict[str, dict]) -> dict:
    energy = payloads["energy_current"]
    energy_exact = energy.get("exact", {})
    if energy_exact.get("fixed_slice_balance") != (
        "d/dt integral_a^b E_pf dx+2[J_pf(t,b)-J_pf(t,a)]"
        "=2 integral_a^b B_pf dx"
    ):
        raise RuntimeError("energy fixed-slice balance drifted")
    bridge = energy_exact.get("spatial_gradient_bridge", {})
    if "B_pf=|partial_x V_pf|^2" not in bridge.get("domination", ""):
        raise RuntimeError("energy gradient domination drifted")
    if energy.get("summary", {}).get("pointwise_contact_exclusion"):
        raise RuntimeError("energy gate was unexpectedly promoted")

    endpoint = payloads["endpoint_margin"]
    endpoint_summary = endpoint.get("summary", {})
    if endpoint_summary.get("left_raw_margins") != 1:
        raise RuntimeError("left endpoint margin ownership drifted")
    if endpoint_summary.get("right_normalized_margins") != 1:
        raise RuntimeError("right endpoint margin ownership drifted")
    if endpoint_summary.get("pointwise_contact_exclusions") != 0:
        raise RuntimeError("endpoint gate was unexpectedly promoted")
    endpoint_exact = endpoint.get("exact", {})
    if endpoint_exact.get("left_source", {}).get("full_margin") != (
        "|J_t'(38)|>3*38^3/2800=20577/350"
    ):
        raise RuntimeError("left endpoint raw margin drifted")
    if endpoint_exact.get("right_energy_margin", {}).get("conclusion") != (
        "sqrt(E_pf(t,R_j))>(99/1000)A_t(R_j)"
    ):
        raise RuntimeError("right endpoint energy margin drifted")

    dominant = payloads["dominant_global"]
    if dominant.get("exact", {}).get("region") != (
        "0<t<=1/2, L=log(x/(4*pi)), t*L>=25, "
        "N_x=floor(sqrt(x/(4*pi)+t/16))"
    ):
        raise RuntimeError("dominant-ray region drifted")
    if "C_t[H_t/A_t](x)>3/40*L^2" not in (
        dominant.get("exact", {}).get("theorem", "")
    ):
        raise RuntimeError("dominant-ray theorem drifted")

    oscillatory = payloads["oscillatory_handoff"]
    oscillatory_theorem = oscillatory.get("exact", {}).get("theorem", "")
    if "there exists L_epsilon" not in oscillatory_theorem:
        raise RuntimeError("oscillatory existential threshold drifted")
    if "does not provide a practical L_epsilon" not in (
        oscillatory.get("proof_boundary", "")
    ):
        raise RuntimeError("oscillatory effectivity boundary drifted")

    critical = payloads["critical_c1"]
    critical_exact = critical.get("exact", {})
    if critical_exact.get("global_c1") != (
        "|r|<2500*exp(-3L/4), |r'|<5000*L*exp(-3L/4), "
        "r^2+(r'/L)^2<32000000*exp(-3L/2)"
    ):
        raise RuntimeError("critical C1 remainder bound drifted")
    if critical_exact.get("remaining_target") != (
        "Prove T_L[J]>32000000*exp(-3L/2) "
        "for the corrected finite main"
    ):
        raise RuntimeError("critical corrected-main target drifted")

    successor = payloads["successor_lemma"]
    transport = successor.get("exact", {}).get("transport", {})
    if transport.get("strict_gate") != "delta_j*M_(j,k)<d_(j,k)":
        raise RuntimeError("phase-cell successor gate drifted")
    if "all-j antecedent open" not in successor.get("status", ""):
        raise RuntimeError("phase-cell all-j boundary drifted")

    finite = payloads["finite_transport"]
    finite_summary = finite.get("summary", {})
    if finite_summary.get("finite_successor_instances") != 1:
        raise RuntimeError("finite successor count drifted")
    if finite_summary.get("certified_transport_panels") != 525:
        raise RuntimeError("finite transport panel count drifted")
    if finite_summary.get("transport_cover") != ["38", "245"]:
        raise RuntimeError("finite transport cover drifted")

    bottom = payloads["q208_bottom"]
    if bottom.get("contract", {}).get("time") != "1/1040":
        raise RuntimeError("Q208 bottom time drifted")
    if bottom.get("contract", {}).get("x_domain") != ["0", "246"]:
        raise RuntimeError("Q208 bottom domain drifted")
    if not bottom.get("summary", {}).get("complete_bottom_edge_certificate"):
        raise RuntimeError("Q208 bottom certificate is incomplete")
    if bottom.get("summary", {}).get("unresolved_boxes") != 0:
        raise RuntimeError("Q208 bottom unresolved-box count drifted")

    top = payloads["q208_top"]
    if top.get("contract", {}).get("time") != "1/5":
        raise RuntimeError("Q208 top time drifted")
    if top.get("contract", {}).get("x_domain") != ["0", "246"]:
        raise RuntimeError("Q208 top domain drifted")
    if not top.get("summary", {}).get("complete_top_edge_certificate"):
        raise RuntimeError("Q208 top certificate is incomplete")
    if top.get("summary", {}).get("unresolved_boxes") != 0:
        raise RuntimeError("Q208 top unresolved-box count drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCES.items()
        },
        "outer_stage_uniform_anchors": [
            {
                "node": "x_0=38",
                "margin": "sqrt(E_pf(t,38))>=m_(L,j)>0",
                "time_domain": "t_(j+1)<=t<=t_j, j>=25",
            },
            {
                "node": "x_N=R_j",
                "margin": (
                    "sqrt(E_pf(t,R_j))>"
                    "(99/1000)A_t(R_j)>0"
                ),
                "time_domain": "t_(j+1)<=t<=t_j, j>=25",
            },
        ],
        "finite_phase_sources": {
            "q208_bottom": {
                "time": bottom["contract"]["time"],
                "x_domain": bottom["contract"]["x_domain"],
                "cells": bottom["summary"]["certified_leaf_cells"],
                "unresolved_boxes": bottom["summary"]["unresolved_boxes"],
            },
            "q208_top": {
                "time": top["contract"]["time"],
                "x_domain": top["contract"]["x_domain"],
                "cells": top["summary"]["certified_leaf_cells"],
                "unresolved_boxes": top["summary"]["unresolved_boxes"],
            },
            "q207_q208_transport": {
                "instances": finite_summary["finite_successor_instances"],
                "panels": finite_summary["certified_transport_panels"],
                "x_cover": finite_summary["transport_cover"],
            },
        },
        "cofinal_internal_anchor_registry": False,
        "inherited_boundary": (
            "The sources own two cofinal outer endpoints and several "
            "finite fixed-time paths, but no stage-uniform current-time "
            "anchor set meeting every bounded interior cell."
        ),
    }


def build_exact() -> dict:
    return {
        "partition": {
            "nodes": "38=x_0<x_1<...<x_N=R_j",
            "cells": "I_k=[x_(k-1),x_k], 1<=k<=N",
            "widths": "h_k=x_k-x_(k-1)>0",
            "endpoint_norms": (
                "e_i(t)=sqrt(E_pf(t,x_i))"
                "=||(H_t,s_pf H_(t,x))(x_i)||_2"
            ),
            "certified_margins": "0<=m_i(t)<=e_i(t)",
        },
        "local_balance": {
            "definition": (
                "D_k(t)=d/dt integral_(x_(k-1))^(x_k) E_pf dx"
                "+2[J_pf(t,x_k)-J_pf(t,x_(k-1))]"
            ),
            "identity": (
                "D_k(t)=2 integral_(I_k) B_pf(t,x) dx>=0"
            ),
            "reason": (
                "B_pf>=|partial_x V_pf|^2>=0, "
                "V_pf=(H_t,s_pf H_(t,x))"
            ),
        },
        "cell_contact_floor": {
            "exact": (
                "V_pf(t,c)=0 in I_k => "
                "D_k(t)>=2[e_(k-1)(t)+e_k(t)]^2/h_k"
            ),
            "certified": (
                "V_pf(t,c)=0 in I_k => "
                "D_k(t)>=T_k(t), "
                "T_k=2[m_(k-1)(t)+m_k(t)]^2/h_k"
            ),
            "one_sided": (
                "m_(k-1)>0, m_k=0 => T_k=2m_(k-1)^2/h_k; "
                "the right-sided analogue is identical"
            ),
            "minimizer": (
                "min_(0<u<h) [e_L^2/u+e_R^2/(h-u)]"
                "=(e_L+e_R)^2/h"
            ),
        },
        "local_criterion": {
            "statement": (
                "If D_k(t)<T_k(t) for every k, then V_pf(t,x)"
                " has no zero on [38,R_j]."
            ),
            "input": (
                "Each local numerator needs either its two boundary fluxes "
                "and energy derivative, or a direct source upper bound for "
                "the local bulk."
            ),
        },
        "telescoping": {
            "identity": (
                "sum_(k=1)^N D_k(t)"
                "=d/dt integral_38^R_j E_pf dx"
                "+2[J_pf(t,R_j)-J_pf(t,38)]"
                "=2 integral_38^R_j B_pf dx=:D_total(t)"
            ),
            "cancellation": (
                "Every internal current J_pf(t,x_k), 1<=k<N, "
                "appears once with sign + and once with sign -."
            ),
            "nonnegative_cells": "D_k(t)>=0 for every k",
        },
        "global_minimum_criterion": {
            "statement": (
                "If T_k(t)>0 for every k and "
                "D_total(t)<min_(1<=k<=N) T_k(t), "
                "then V_pf(t,x) has no zero on [38,R_j]."
            ),
            "proof": (
                "A contact in I_r gives D_r>=T_r, while "
                "D_total=sum_k D_k>=D_r because every D_k>=0."
            ),
            "anchor_condition": (
                "T_k>0 exactly when m_(k-1)+m_k>0; hence the "
                "certified node set must meet every partition edge."
            ),
        },
        "sum_threshold_guard": {
            "false_shortcut": (
                "D_total<sum_k T_k does not exclude a contact."
            ),
            "exact_budget_counterexample": (
                "T_1=1, T_2=100, D_1=1, D_2=0 gives "
                "D_total=1<T_1+T_2=101 while D_1>=T_1."
            ),
            "boundary": (
                "This is an exact implication counterexample for the "
                "budget logic, not a claimed realization by the Xi field."
            ),
        },
        "anchor_graph": {
            "graph": (
                "The partition is the path graph on nodes x_0,...,x_N; "
                "owned positive margins must form a vertex cover."
            ),
            "current_uniform_set": "{x_0=38,x_N=R_j}",
            "consequence": (
                "The two owned outer nodes cover both cells only when N<=2. "
                "For N>=3 at least one interior edge has neither endpoint "
                "owned; bounded-width localization therefore needs a "
                "cofinal interior anchor skeleton."
            ),
        },
        "flux_localization_tradeoff": {
            "local_side": (
                "Keeping the thresholds cellwise preserves their 1/h_k "
                "strength but also keeps the internal fluxes in D_k."
            ),
            "summed_side": (
                "Summing cancels every internal flux but compresses all "
                "contact information to the weakest threshold min_k T_k."
            ),
            "conclusion": (
                "Telescoping is exact bookkeeping, not a substitute for "
                "current-time interior endpoint ownership."
            ),
        },
        "ownership_audit": {
            "compact": (
                "At x=38 and 0<t<=1/5, the compact certificate gives "
                "|J_t'(38)|>20577/350 and the inherited margin "
                "sqrt(E_pf(t,38))>=m_(L,j)>0."
            ),
            "dominant": (
                "The dominant theorem proves strict normalized Laguerre "
                "curvature on tL>=25. The endpoint gate exposes the "
                "quantitative margin at x=R_j. Uniformly over a descendant "
                "collar, the dominant region shrinks to x=R_j at the bottom "
                "time, so it does not populate a bounded-cell skeleton "
                "through the old interior."
            ),
            "oscillatory": (
                "For each epsilon>0 the oscillatory theorem gives an "
                "existential L_epsilon above c_*+epsilon, "
                "c_*=4911678521/1933561194, but no practical L_epsilon or "
                "effective cofinal endpoint-margin registry."
            ),
            "critical_c1": (
                "For L>=50 and 0<tL<25, the corrected global C1 remainder "
                "satisfies r^2+(r'/L)^2<32000000 exp(-3L/2), but the "
                "finite-main lower bound "
                "T_L[J]>32000000 exp(-3L/2) remains open."
            ),
            "finite_phase": (
                "The Q208 top and bottom certificates own fixed paths "
                "t=1/5 and t=1/1040 on 0<=x<=246. The Q207-to-Q208 "
                "transport owns one finite successor on 38<=x<=245."
            ),
            "conditional_successor": (
                "The all-stage phase-cell lemma would transfer margins if "
                "delta_j M_(j,k)<d_(j,k) on every cell, but that Xi "
                "antecedent is precisely open."
            ),
        },
        "nonpromotion": {
            "snapshot": (
                "A finite fixed-time phase chain cannot be reused as a "
                "current-time cofinal anchor registry."
            ),
            "positivity": (
                "Qualitative strict positivity with an ineffective height "
                "does not provide the effective positive m_i needed by the "
                "minimum-threshold criterion."
            ),
            "partition": (
                "Choosing more cells improves 1/h_k only after supplying "
                "more certified anchors; partition refinement alone proves "
                "nothing pointwise."
            ),
        },
        "route_decision": (
            "Retain the local and globally telescoped criteria as exact "
            "conditional lemmas, but park the energy-localization route. "
            "The audited sources do not supply a cofinal current-time "
            "interior anchor vertex cover, and attempting to transport one "
            "reintroduces the open all-j phase-cell estimate. Return to the "
            "critical contact-normal C1 arithmetic target."
        ),
        "open_handoff": (
            "Prove the corrected finite-main inequality "
            "T_L[J]>32000000 exp(-3L/2) for L>=50 and 0<tL<25, with "
            "cutoff transitions retained. First derive a cutoff-stable "
            "phase/amplitude decomposition of T_L[J] and identify the "
            "smallest signed cancellation statement that would force the "
            "lower bound; do not divide by the unknown first-jet norm."
        ),
    }


def build_payload() -> dict:
    payloads = load_sources()
    source_audit = audit_sources(payloads)
    exact = build_exact()
    rows = [
        GateRow(
            id="raecl_01_partition",
            role="exact_geometry",
            readiness="proved_exact",
            claim="Fix an arbitrary finite partition of the old descendant collar.",
            formula=(
                f"{exact['partition']['nodes']}; "
                f"{exact['partition']['cells']}; "
                f"{exact['partition']['widths']}"
            ),
            proof_boundary="No regular spacing or numerical mesh is assumed.",
        ),
        GateRow(
            id="raecl_02_local_balance",
            role="exact_identity",
            readiness="proved_exact",
            claim="The energy/current identity localizes exactly to every cell.",
            formula=(
                f"{exact['local_balance']['definition']}; "
                f"{exact['local_balance']['identity']}"
            ),
            proof_boundary=exact["local_balance"]["reason"],
        ),
        GateRow(
            id="raecl_03_cell_nonnegativity",
            role="exact_inequality",
            readiness="proved_exact",
            claim="Every local energy/current numerator is nonnegative.",
            formula=exact["local_balance"]["identity"],
            proof_boundary="This uses the proved outer-domain bulk coercivity.",
        ),
        GateRow(
            id="raecl_04_two_sided_contact_floor",
            role="exact_inequality",
            readiness="proved_exact",
            claim="A contact in one cell forces a two-sided local bulk floor.",
            formula=exact["cell_contact_floor"]["exact"],
            proof_boundary=exact["cell_contact_floor"]["minimizer"],
        ),
        GateRow(
            id="raecl_05_certified_threshold",
            role="exact_composition",
            readiness="available_exact",
            claim="Certified endpoint margins give a usable cell threshold.",
            formula=exact["cell_contact_floor"]["certified"],
            proof_boundary="The margins must be derived independently of contact exclusion.",
        ),
        GateRow(
            id="raecl_06_one_sided_threshold",
            role="exact_composition",
            readiness="available_exact",
            claim="One positive endpoint anchor is enough to make a cell threshold positive.",
            formula=exact["cell_contact_floor"]["one_sided"],
            proof_boundary="This is weaker than the two-sided floor but still exact.",
        ),
        GateRow(
            id="raecl_07_local_no_contact",
            role="conditional_criterion",
            readiness="not_ready_to_apply",
            claim="Strictly beating every local threshold excludes all contacts.",
            formula=exact["local_criterion"]["statement"],
            proof_boundary=exact["local_criterion"]["input"],
        ),
        GateRow(
            id="raecl_08_flux_telescoping",
            role="exact_identity",
            readiness="proved_exact",
            claim="All internal currents cancel when the cell balances are summed.",
            formula=exact["telescoping"]["identity"],
            proof_boundary=exact["telescoping"]["cancellation"],
        ),
        GateRow(
            id="raecl_09_global_minimum_criterion",
            role="conditional_criterion",
            readiness="not_ready_to_apply",
            claim="The summed numerator can exclude contact only below the weakest cell threshold.",
            formula=exact["global_minimum_criterion"]["statement"],
            proof_boundary=exact["global_minimum_criterion"]["proof"],
        ),
        GateRow(
            id="raecl_10_sum_threshold_guard",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The sum of cell thresholds is not a valid global contact threshold.",
            formula=exact["sum_threshold_guard"]["exact_budget_counterexample"],
            proof_boundary=exact["sum_threshold_guard"]["boundary"],
        ),
        GateRow(
            id="raecl_11_anchor_vertex_cover",
            role="exact_graph_reduction",
            readiness="proved_exact",
            claim="Positive owned endpoint margins must meet every partition cell.",
            formula=exact["global_minimum_criterion"]["anchor_condition"],
            proof_boundary=exact["anchor_graph"]["graph"],
        ),
        GateRow(
            id="raecl_12_outer_anchor_obstruction",
            role="source_audit",
            readiness="not_ready_to_apply",
            claim="The two cofinal outer anchors do not cover a bounded-width interior partition.",
            formula=exact["anchor_graph"]["consequence"],
            proof_boundary="This concerns currently audited stage-uniform ownership only.",
        ),
        GateRow(
            id="raecl_13_flux_localization_tradeoff",
            role="exact_route_guard",
            readiness="guard_validated",
            claim="Flux cancellation and strong local thresholds cannot be obtained for free at the same time.",
            formula=exact["flux_localization_tradeoff"]["conclusion"],
            proof_boundary=(
                f"{exact['flux_localization_tradeoff']['local_side']} "
                f"{exact['flux_localization_tradeoff']['summed_side']}"
            ),
        ),
        GateRow(
            id="raecl_14_compact_ownership",
            role="source_audit",
            readiness="certified",
            claim="The compact source owns the left outer node at every active time.",
            formula=exact["ownership_audit"]["compact"],
            proof_boundary="Inherited from the checked endpoint-margin gate.",
        ),
        GateRow(
            id="raecl_15_dominant_ownership",
            role="source_audit",
            readiness="certified_boundary_only",
            claim="The dominant source owns the moving right end but not a uniform old-interior skeleton.",
            formula=exact["ownership_audit"]["dominant"],
            proof_boundary="At the collar bottom the dominant tail begins at R_j.",
        ),
        GateRow(
            id="raecl_16_oscillatory_effectivity",
            role="source_audit",
            readiness="not_effective_for_grid",
            claim="The oscillatory theorem is asymptotic and does not expose effective cell anchors.",
            formula=exact["ownership_audit"]["oscillatory"],
            proof_boundary="No practical L_epsilon is present in the source theorem.",
        ),
        GateRow(
            id="raecl_17_critical_c1_target",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The critical remainder is controlled, but the corrected finite-main first-jet floor is open.",
            formula=exact["ownership_audit"]["critical_c1"],
            proof_boundary="This is the live arithmetic transversality obligation.",
        ),
        GateRow(
            id="raecl_18_finite_phase_ownership",
            role="source_audit",
            readiness="finite_only",
            claim="The rigorous phase-cell margins currently cover fixed finite paths and one successor.",
            formula=exact["ownership_audit"]["finite_phase"],
            proof_boundary="They do not quantify every later descendant collar.",
        ),
        GateRow(
            id="raecl_19_successor_and_snapshot_guard",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Finite snapshots and a conditional successor lemma do not provide a cofinal anchor registry.",
            formula=(
                f"{exact['ownership_audit']['conditional_successor']} "
                f"{exact['nonpromotion']['snapshot']}"
            ),
            proof_boundary=(
                f"{exact['nonpromotion']['positivity']} "
                f"{exact['nonpromotion']['partition']}"
            ),
        ),
        GateRow(
            id="raecl_20_route_decision",
            role="route_decision",
            readiness="available_exact",
            claim="Park localization and return to the corrected finite-main C1 target.",
            formula=exact["open_handoff"],
            proof_boundary=exact["route_decision"],
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-29",
        "status": (
            "exact cell-local energy criteria, flux telescoping, and "
            "current-time endpoint-ownership obstruction"
        ),
        "proof_boundary": (
            "This artifact proves the cell-local balance, local contact "
            "floors, one-sided threshold, exact internal-flux telescoping, "
            "global minimum-threshold criterion, vertex-cover ownership "
            "condition, and source ownership audit. It does not prove a "
            "cofinal interior anchor registry, any strict Xi local or "
            "global bulk upper budget, the corrected finite-main C1 lower "
            "bound, Q209, the cofinal descendant theorem, Lambda<=0, "
            "PF-infinity, RH, or a Clay-prize conclusion."
        ),
        "source_audit": source_audit,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_local_balances": 2,
            "exact_contact_floors": 2,
            "conditional_no_contact_criteria": 2,
            "outer_stage_uniform_anchors": 2,
            "cofinal_internal_anchor_registries": 0,
            "finite_successor_instances": 1,
            "nonpromotion_guards": 2,
            "open_corrected_main_targets": 1,
            "pointwise_contact_exclusions": 0,
            "cofinal_descendant_theorem": False,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    audit = payload["source_audit"]
    finite = audit["finite_phase_sources"]
    return "\n".join(
        [
            "# Ray-Aligned Energy Cell-Localization Gate",
            "",
            "Date: 2026-07-29",
            "",
            "Status: exact localization and endpoint-ownership audit. This",
            "is not a proof of contact exclusion, Lambda<=0, or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated ray-aligned energy cell-localization gate: 20 rows, "
            "2 exact local balances, 2 exact contact floors, "
            "2 conditional no-contact criteria, 2 outer anchors, "
            "0 cofinal internal anchor registries, "
            "1 finite successor instance, 2 nonpromotion guards, "
            "0 pointwise contact exclusions",
            "```",
            "",
            "## Cell Identity",
            "",
            "For",
            "",
            "```text",
            exact["partition"]["nodes"],
            exact["partition"]["cells"],
            exact["partition"]["widths"],
            "```",
            "",
            "the fixed-slice balance gives",
            "",
            "```text",
            exact["local_balance"]["definition"],
            exact["local_balance"]["identity"],
            "```",
            "",
            "A contact in `I_k` therefore forces",
            "",
            "```text",
            exact["cell_contact_floor"]["exact"],
            exact["cell_contact_floor"]["certified"],
            exact["cell_contact_floor"]["one_sided"],
            "```",
            "",
            "The `1/h_k` gain is real. It is useful only when at least one",
            "endpoint of that cell has an independently certified positive",
            "current-time margin.",
            "",
            "## Telescoping",
            "",
            "Summing the local identities gives",
            "",
            "```text",
            exact["telescoping"]["identity"],
            "```",
            "",
            exact["telescoping"]["cancellation"],
            "",
            "The exact summed criterion is",
            "",
            "```text",
            exact["global_minimum_criterion"]["statement"],
            "```",
            "",
            exact["global_minimum_criterion"]["proof"],
            "",
            "The tempting replacement of the minimum by a sum is false:",
            "",
            "```text",
            exact["sum_threshold_guard"]["exact_budget_counterexample"],
            "```",
            "",
            exact["sum_threshold_guard"]["boundary"],
            "",
            "## Endpoint Ownership",
            "",
            "The positive-margin nodes form a vertex-cover problem on the",
            "partition path:",
            "",
            "```text",
            exact["anchor_graph"]["current_uniform_set"],
            exact["anchor_graph"]["consequence"],
            "```",
            "",
            "Source audit:",
            "",
            f"- Compact: {exact['ownership_audit']['compact']}",
            f"- Dominant: {exact['ownership_audit']['dominant']}",
            f"- Oscillatory: {exact['ownership_audit']['oscillatory']}",
            f"- Critical C1: {exact['ownership_audit']['critical_c1']}",
            f"- Finite phase cells: {exact['ownership_audit']['finite_phase']}",
            f"- Conditional successor: {exact['ownership_audit']['conditional_successor']}",
            "",
            "The finite inventory is exact:",
            "",
            "```text",
            f"Q208 bottom: t={finite['q208_bottom']['time']}, "
            f"{finite['q208_bottom']['cells']} cells, "
            f"{finite['q208_bottom']['unresolved_boxes']} unresolved",
            f"Q208 top: t={finite['q208_top']['time']}, "
            f"{finite['q208_top']['cells']} cells, "
            f"{finite['q208_top']['unresolved_boxes']} unresolved",
            f"Q207-Q208 transport: {finite['q207_q208_transport']['instances']} "
            f"successor, {finite['q207_q208_transport']['panels']} panels",
            "```",
            "",
            "These are finite certificates, not a cofinal current-time",
            "anchor registry.",
            "",
            "## Exact Tradeoff",
            "",
            exact["flux_localization_tradeoff"]["local_side"],
            "",
            exact["flux_localization_tradeoff"]["summed_side"],
            "",
            exact["flux_localization_tradeoff"]["conclusion"],
            "",
            "## Route Decision",
            "",
            exact["route_decision"],
            "",
            "Open handoff:",
            "",
            exact["open_handoff"],
            "",
            "Pi provenance is unchanged: every occurrence comes from the",
            "established coordinate `L=log(x/(4*pi))`, where `pi` is the",
            "ordinary circle constant.",
            "",
            "## Boundary",
            "",
            payload["proof_boundary"],
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
    args.out.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "built ray-aligned energy cell-localization gate: 20 rows, "
        "2 exact local balances, 2 exact contact floors, "
        "2 conditional no-contact criteria, 2 outer anchors, "
        "0 cofinal internal anchor registries, "
        "1 finite successor instance, 2 nonpromotion guards, "
        "0 pointwise contact exclusions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
