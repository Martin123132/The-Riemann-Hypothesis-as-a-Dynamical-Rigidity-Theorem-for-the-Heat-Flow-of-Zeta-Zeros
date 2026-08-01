#!/usr/bin/env python3
"""Build the exact Q208-base/Q209 successor-shell coverage gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_q208_base_q209_shell_coverage_gate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "cofinal_target": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.json"
    ),
    "q208_closed": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_q208_closed_boundary_winding_certificate.json"
    ),
    "q208_bottom": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_q208_bottom_phase_cell_certificate.json"
    ),
    "q208_right": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_q208_selected_boundary_pilot.json"
    ),
    "q207_q208_forward": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_q207_q208_forward_adiabatic_successor_certificate.json"
    ),
    "compact_core": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_theta_compact_transversality_interval_certificate.json"
    ),
    "successor_lemma": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma.json"
    ),
    "scaled_successor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_time_dependent_scaled_successor_lemma.json"
    ),
    "ray_aligned": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_ray_aligned_parabolic_frequency_reduction.json"
    ),
}


@dataclass(frozen=True)
class CoverageRow:
    id: str
    region: str
    status: str
    source: str
    obligation: str
    proof_boundary: str


@dataclass(frozen=True)
class BoundaryRow:
    id: str
    edge: str
    parameterization: str
    arc: str
    status: str
    reason: str


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


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


def source_audit(payloads: dict[str, dict]) -> dict:
    cofinal = payloads["cofinal_target"]["exact"]["cofinal_contract"]
    if "Q_207=" not in cofinal["current_base"]:
        raise RuntimeError("historical cofinal source no longer records Q207")
    if "j>=208" not in cofinal["open_boundary_theorem"]:
        raise RuntimeError("historical cofinal open index drifted")

    q208 = payloads["q208_closed"]
    if q208.get("exact_winding") != 0:
        raise RuntimeError("Q208 winding source is not zero")
    if q208.get("domain", {}).get("q208") != "[1/1040,1/4]x[0,246]":
        raise RuntimeError("Q208 domain source drifted")
    if "Q_208 is contact-free" not in q208.get("conclusion", ""):
        raise RuntimeError("Q208 finite theorem source missing")

    bottom = payloads["q208_bottom"]
    contract = bottom.get("contract", {})
    summary = bottom.get("summary", {})
    if contract.get("time") != "1/1040":
        raise RuntimeError("Q208 bottom time drifted")
    if contract.get("x_domain") != ["0", "246"]:
        raise RuntimeError("Q208 bottom domain drifted")
    if contract.get("initial_x_step") != "1/2":
        raise RuntimeError("Q208 bottom phase-cell step drifted")
    if summary.get("tasks_total") != 492:
        raise RuntimeError("Q208 bottom cell count drifted")
    if not summary.get("complete_bottom_edge_certificate"):
        raise RuntimeError("Q208 bottom certificate is incomplete")

    right = payloads["q208_right"].get("summary", {}).get("right_strip", {})
    if right.get("domain") != "[1/1040,1/5]x[245,246]":
        raise RuntimeError("Q208 right-strip source domain drifted")
    if not right.get("full_derivative_positive"):
        raise RuntimeError("Q208 right-strip source lost positivity")

    forward = payloads["q207_q208_forward"]
    forward_contract = forward.get("contract", {})
    if forward_contract.get("old_time") != "1/1035":
        raise RuntimeError("Q207 forward old time drifted")
    if forward_contract.get("new_time") != "1/1040":
        raise RuntimeError("Q208 forward new time drifted")
    if forward_contract.get("new_radius") != "246":
        raise RuntimeError("Q208 forward radius drifted")
    if not forward.get("summary", {}).get("complete_forward_collar"):
        raise RuntimeError("Q207-to-Q208 forward collar is incomplete")

    compact = payloads["compact_core"].get("exact", {})
    if compact.get("proved_rectangle") != (
        "0<=t<=1/5 and 1/4<=x<=38"
    ):
        raise RuntimeError("compact transversality rectangle drifted")
    if "|x|<=1/4" not in compact.get("origin_overlap", ""):
        raise RuntimeError("compact origin overlap missing")

    successor = payloads["successor_lemma"]
    shell = successor.get("exact", {}).get("shell", {})
    if shell.get("parameters") != "t_j=1/(5*j), R_j=j+38":
        raise RuntimeError("linear successor schedule drifted")
    if shell.get("decomposition") != "P_(j+1)=P_j union C_j union S_j":
        raise RuntimeError("linear successor decomposition drifted")
    if "delta_j*M_(j,k)<d_(j,k)" not in successor.get(
        "exact", {}
    ).get("transport", {}).get("strict_gate", ""):
        raise RuntimeError("successor transport gate missing")
    if "fixed nonzero q_j" not in successor.get(
        "exact", {}
    ).get("right_half_plane", {}).get("gate", ""):
        raise RuntimeError("successor right half-plane gate missing")
    if "before using Q209 as anything beyond a falsification case" not in (
        successor.get("next_exact_obligation", "")
    ):
        raise RuntimeError("finite-Q209 nonpromotion source missing")

    scaled = payloads["scaled_successor"].get("exact_summary", {})
    if not scaled.get("relative_transport_available"):
        raise RuntimeError("relative scaled transport source missing")
    if scaled.get("all_j_xi_theorem"):
        raise RuntimeError("scaled source overpromoted to an all-j theorem")

    ray = payloads["ray_aligned"]
    ray_summary = ray.get("summary", {})
    if ray_summary.get("new_strip_open_antecedents") != 0:
        raise RuntimeError("ray-aligned new strips are no longer closed")
    if ray_summary.get("old_collar_open_antecedents") != 1:
        raise RuntimeError("ray-aligned old-collar count drifted")
    if ray.get("exact", {}).get("scale") != (
        "s_pf(t,x)=(L(x)^2+(2t)^(-1))^(-1/2)"
    ):
        raise RuntimeError("parabolic-frequency scale drifted")

    return {
        "source_kinds": {
            key: payload.get("kind", "")
            for key, payload in payloads.items()
        },
        "source_sha256": {
            key: file_hash(path)
            for key, path in SOURCES.items()
        },
        "historical_handoff_correction": (
            "The cofinal target's Q207/j>=208 handoff is retained as a "
            "historical artifact. The later rigorous Q208 closed-boundary "
            "certificate supersedes only that finite base: Q208 is proved "
            "and the first unresolved linear stage is Q209."
        ),
    }


def build_exact() -> dict:
    old_index = 208
    new_index = 209
    t_old = Fraction(1, 5 * old_index)
    t_new = Fraction(1, 5 * new_index)
    low_top = Fraction(1, 5)
    full_top = Fraction(1, 4)
    radius_old = old_index + 38
    radius_new = new_index + 38
    compact_radius = 38
    delta = t_old - t_new
    half_step = Fraction(1, 2)
    outer_width = radius_old - compact_radius
    outer_cells = int(Fraction(outer_width, 1) / half_step)
    new_strip_columns = int(
        Fraction(radius_new - radius_old, 1) / half_step
    )

    inherited_area = (low_top - t_old) * radius_old
    compact_collar_area = delta * compact_radius
    outer_collar_area = delta * outer_width
    new_strip_area = (low_top - t_new) * (radius_new - radius_old)
    low_total_area = (low_top - t_new) * radius_new
    decomposed_area = (
        inherited_area
        + compact_collar_area
        + outer_collar_area
        + new_strip_area
    )
    if decomposed_area != low_total_area:
        raise RuntimeError("Q209 low-shell area decomposition failed")
    if delta != Fraction(1, 217360):
        raise RuntimeError("Q208-to-Q209 time step failed")
    if outer_cells != 416:
        raise RuntimeError("Q209 outer-collar phase-cell count failed")
    if new_strip_columns != 2:
        raise RuntimeError("Q209 right-strip half-cell count failed")

    vertices = [
        (Fraction(0), t_new),
        (Fraction(radius_new), t_new),
        (Fraction(radius_new), full_top),
        (Fraction(0), full_top),
    ]
    twice_signed_area = sum(
        x_0 * t_1 - t_0 * x_1
        for (x_0, t_0), (x_1, t_1) in zip(
            vertices, vertices[1:] + vertices[:1], strict=True
        )
    )
    expected_twice_area = 2 * radius_new * (full_top - t_new)
    if twice_signed_area != expected_twice_area:
        raise RuntimeError("Q209 boundary orientation audit failed")
    if twice_signed_area <= 0:
        raise RuntimeError("Q209 boundary is not counterclockwise")

    return {
        "indices": {
            "proved_base": old_index,
            "first_unresolved_linear_stage": new_index,
        },
        "linear_schedule": {
            "formula": "t_j=1/(5*j), R_j=j+38",
            "t_208": fraction_text(t_old),
            "t_209": fraction_text(t_new),
            "R_208": str(radius_old),
            "R_209": str(radius_new),
            "delta_208": fraction_text(delta),
        },
        "rectangles": {
            "Q_208": "[1/1040,1/4]x[0,246]",
            "Q_209": "[1/1045,1/4]x[0,247]",
            "P_208": "[1/1040,1/5]x[0,246]",
            "P_209": "[1/1045,1/5]x[0,247]",
            "C_208": "[1/1045,1/1040]x[0,246]",
            "C_208_compact": "[1/1045,1/1040]x[0,38]",
            "C_208_outer": "[1/1045,1/1040]x[38,246]",
            "S_208": "[1/1045,1/5]x[246,247]",
        },
        "low_shell_decomposition": {
            "identity": "P_209=P_208 union C_208 union S_208",
            "refined_identity": (
                "P_209=P_208 union C_208_compact union "
                "C_208_outer union S_208"
            ),
            "interiors_are_disjoint": True,
            "area_P_208": fraction_text(inherited_area),
            "area_compact_collar": fraction_text(compact_collar_area),
            "area_outer_collar": fraction_text(outer_collar_area),
            "area_right_strip": fraction_text(new_strip_area),
            "area_P_209": fraction_text(low_total_area),
            "area_sum_verified": decomposed_area == low_total_area,
        },
        "phase_cell_counts": {
            "q208_bottom_half_unit_cells": 492,
            "outer_old_half_unit_cells": outer_cells,
            "new_right_strip_half_unit_columns": new_strip_columns,
        },
        "boundary_orientation": {
            "coordinate_order": "(x,t)",
            "orientation": "counterclockwise",
            "bottom": "gamma_B(x)=(x,1/1045), 0<=x<=247",
            "right": "gamma_R(t)=(247,t), 1/1045<=t<=1/4",
            "top": "gamma_T(x)=(x,1/4), 247>=x>=0",
            "left": "gamma_L(t)=(0,t), 1/4>=t>=1/1045",
            "twice_signed_area": fraction_text(twice_signed_area),
            "positive_orientation_verified": twice_signed_area > 0,
        },
        "conditional_successor": {
            "outer_transport": (
                "For every one of the 416 old outer half-unit cells "
                "I_(208,k), prove delta_208*M_(208,k)<d_(208,k) "
                "using enclosures valid on "
                "[1/1045,1/1040]xI_(208,k)."
            ),
            "new_strip_half_plane": (
                "Construct one fixed q_208!=0 and eta_208>0 with "
                "q_208 dot V>=eta_208 throughout "
                "[1/1045,1/5]x[246,247]."
            ),
            "conclusion": (
                "Those two antecedents, the compact-core theorem, and the "
                "proved Q208 base imply that P209 is contact-free with "
                "zero first-jet degree. The published Lambda<=1/5 theorem "
                "then covers 1/5<t<=1/4; equality on the new x interval "
                "must come from the new-strip antecedent. Only then is "
                "Q209 certified."
            ),
        },
        "noninheritance_guards": {
            "time_enclosures": (
                "The Q207-to-Q208 derivative enclosures were proved only "
                "on [1/1040,1/1035]. The smaller scalar delta_208 does not "
                "extend them to the disjoint lower-time interval "
                "[1/1045,1/1040]. New interval enclosures are required."
            ),
            "space_translation": (
                "The proved right strip is [245,246]. Oscillatory Xi phase "
                "has no established monotone x-translation theorem, so its "
                "half-plane certificate cannot be shifted to [246,247]."
            ),
            "finite_nonpromotion": (
                "A finite Q209 computation may calibrate or falsify a "
                "candidate estimate, but it is not the parameter-uniform "
                "successor theorem."
            ),
        },
        "ray_aligned_route": {
            "schedule": (
                "c=25, a=100, t_j=25/(100+j), "
                "L_j=101+j, R_j=4*pi*exp(L_j)"
            ),
            "new_strip_identity": (
                "t_(j+1)*L(x)>=t_(j+1)*L_j=25 on every new strip"
            ),
            "new_strip_status": (
                "closed by the existing dominant-saddle tL>=25 theorem"
            ),
            "only_open_successor_region": (
                "C_j^out=[t_(j+1),t_j]x[38,R_j]"
            ),
            "scale": (
                "s_pf=(L^2+(2t)^(-1))^(-1/2)"
                "=sqrt(2t)/sqrt(1+2tL^2)"
            ),
            "open_xi_antecedent": (
                "Construct an a priori source-level kappa_j(t) such that "
                "||partial_t(H,s_pf H_x)||<=kappa_j(t)||(H,s_pf H_x)|| "
                "on every outer descendant collar, with finite one-step "
                "integral. Defining kappa as the quotient of these two "
                "unknown norms is forbidden because it assumes nonvanishing."
            ),
            "pi_provenance": (
                "Here pi is the usual circle constant already present in "
                "the established Xi frequency coordinate "
                "L(x)=log(x/(4*pi)); R_j=4*pi*exp(L_j) is exactly the "
                "inverse coordinate change. No polygon, curvature image, "
                "or arbitrarily chosen circle introduces this pi."
            ),
        },
    }


def build_payload() -> dict:
    sources = load_sources()
    audit = source_audit(sources)
    exact = build_exact()
    coverage = [
        CoverageRow(
            id="q209_cov_01_inherited_p208",
            region=exact["rectangles"]["P_208"],
            status="covered_exact",
            source="Q208 closed-boundary/no-contact certificate",
            obligation="none",
            proof_boundary="Finite inherited low rectangle only.",
        ),
        CoverageRow(
            id="q209_cov_02_compact_collar",
            region=exact["rectangles"]["C_208_compact"],
            status="covered_exact",
            source="compact transversality theorem",
            obligation="none",
            proof_boundary="Uses the union of the origin overlap and 1/4<=x<=38 theorem.",
        ),
        CoverageRow(
            id="q209_cov_03_outer_collar",
            region=exact["rectangles"]["C_208_outer"],
            status="open_xi_antecedent",
            source="adiabatic phase-cell successor lemma",
            obligation=exact["conditional_successor"]["outer_transport"],
            proof_boundary="No interval enclosure valid on the Q209 collar is currently stored.",
        ),
        CoverageRow(
            id="q209_cov_04_new_right_strip",
            region=exact["rectangles"]["S_208"],
            status="open_xi_antecedent",
            source="adiabatic phase-cell successor lemma",
            obligation=exact["conditional_successor"]["new_strip_half_plane"],
            proof_boundary="The Q208 [245,246] half-plane theorem cannot be translated in x.",
        ),
        CoverageRow(
            id="q209_cov_05_high_time",
            region="(1/5,1/4]x[0,247]",
            status="covered_exact",
            source="published Lambda<=1/5 composition",
            obligation="none for t>1/5",
            proof_boundary="At t=1/5 and 246<=x<=247 the new-strip antecedent is still required.",
        ),
    ]
    boundary = [
        BoundaryRow(
            id="q209_edge_01_bottom_compact",
            edge="bottom",
            parameterization="x increases at t=1/1045",
            arc="0<=x<=38",
            status="covered_exact",
            reason="compact transversality theorem",
        ),
        BoundaryRow(
            id="q209_edge_02_bottom_outer",
            edge="bottom",
            parameterization="x increases at t=1/1045",
            arc="38<=x<=246",
            status="open_xi_antecedent",
            reason="outer-collar transport is unproved",
        ),
        BoundaryRow(
            id="q209_edge_03_bottom_new",
            edge="bottom",
            parameterization="x increases at t=1/1045",
            arc="246<=x<=247",
            status="open_xi_antecedent",
            reason="new right-strip half-plane is unproved",
        ),
        BoundaryRow(
            id="q209_edge_04_right_low",
            edge="right",
            parameterization="t increases at x=247",
            arc="1/1045<=t<=1/5",
            status="open_xi_antecedent",
            reason="new right-strip half-plane is unproved",
        ),
        BoundaryRow(
            id="q209_edge_05_right_high",
            edge="right",
            parameterization="t increases at x=247",
            arc="1/5<t<=1/4",
            status="covered_exact",
            reason="published Lambda<=1/5 composition",
        ),
        BoundaryRow(
            id="q209_edge_06_top",
            edge="top",
            parameterization="x decreases at t=1/4",
            arc="247>=x>=0",
            status="covered_exact",
            reason="t=1/4 is strictly above 1/5",
        ),
        BoundaryRow(
            id="q209_edge_07_left",
            edge="left",
            parameterization="t decreases at x=0",
            arc="1/4>=t>=1/1045",
            status="covered_exact",
            reason="high-time theorem plus compact origin theorem",
        ),
    ]
    rows = [
        GateRow(
            id="q209_gate_01_handoff_correction",
            role="source_correction",
            readiness="available_exact",
            claim="The rigorous Q208 theorem supersedes the older Q207 finite-base handoff.",
            formula="proved base Q208; first unresolved linear stage Q209",
            proof_boundary="The historical cofinal artifact is retained and not rewritten.",
        ),
        GateRow(
            id="q209_gate_02_q208_base",
            role="exact_composition",
            readiness="available_exact",
            claim="Q208 is contact-free and has zero first-jet winding.",
            formula=exact["rectangles"]["Q_208"],
            proof_boundary="Finite Q208 only.",
        ),
        GateRow(
            id="q209_gate_03_time_step",
            role="exact_geometry",
            readiness="proved_exact",
            claim="The Q208-to-Q209 linear time step is exact.",
            formula="delta_208=1/1040-1/1045=1/217360",
            proof_boundary="Rational schedule arithmetic only.",
        ),
        GateRow(
            id="q209_gate_04_shell_decomposition",
            role="exact_geometry",
            readiness="proved_exact",
            claim="The complete Q209 low rectangle splits into inherited base, compact collar, outer collar, and new strip.",
            formula=exact["low_shell_decomposition"]["refined_identity"],
            proof_boundary="Set identity with boundary overlaps only.",
        ),
        GateRow(
            id="q209_gate_05_area_audit",
            role="exact_geometry",
            readiness="proved_exact",
            claim="The four low-shell pieces have disjoint interiors and exactly exhaust P209 by rational area.",
            formula=(
                f"area(P209)={exact['low_shell_decomposition']['area_P_209']}"
            ),
            proof_boundary="Area equality supplements, but does not replace, the endpoint identity.",
        ),
        GateRow(
            id="q209_gate_06_orientation",
            role="exact_topology",
            readiness="proved_exact",
            claim="The four Q209 edge parameterizations have the standard counterclockwise (x,t) orientation.",
            formula=(
                "twice signed area="
                + exact["boundary_orientation"]["twice_signed_area"]
                + ">0"
            ),
            proof_boundary="Orientation audit only; no winding is computed for Q209.",
        ),
        GateRow(
            id="q209_gate_07_inherited_coverage",
            role="exact_composition",
            readiness="available_exact",
            claim="P208 and the compact part of the new collar are already contact-free.",
            formula="P208 union C_208_compact",
            proof_boundary="Does not cover x>38 below t=1/1040.",
        ),
        GateRow(
            id="q209_gate_08_outer_gap",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The outer old descendant collar is the first open Q209 region.",
            formula=exact["rectangles"]["C_208_outer"],
            proof_boundary="All 416 new-time transport inequalities remain unproved.",
        ),
        GateRow(
            id="q209_gate_09_strip_gap",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The new right strip is the second open Q209 region.",
            formula=exact["rectangles"]["S_208"],
            proof_boundary="No fixed-half-plane certificate exists on [246,247].",
        ),
        GateRow(
            id="q209_gate_10_high_time",
            role="exact_composition",
            readiness="available_exact",
            claim="The high-time band t>1/5 is already contact-free at every finite x.",
            formula="(1/5,1/4]x[0,247]",
            proof_boundary="Equality t=1/5 on the new strip is not imported.",
        ),
        GateRow(
            id="q209_gate_11_time_noninheritance",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="A smaller successor time step does not extend old derivative enclosures to a disjoint time interval.",
            formula="[1/1040,1/1035] does not certify [1/1045,1/1040]",
            proof_boundary="Requires new interval arithmetic or an analytic uniform theorem.",
        ),
        GateRow(
            id="q209_gate_12_space_noninheritance",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The Q208 right-strip half-plane theorem has no certified one-unit spatial translation.",
            formula="[245,246] does not certify [246,247]",
            proof_boundary="No monotonic spatial translation theorem is available.",
        ),
        GateRow(
            id="q209_gate_13_conditional_q209",
            role="conditional_theorem",
            readiness="open_antecedent",
            claim="Closing exactly the two listed regions would certify Q209 with zero first-jet degree.",
            formula="outer transport plus new-strip half-plane implies Q209",
            proof_boundary="Both Xi antecedents are presently open.",
        ),
        GateRow(
            id="q209_gate_14_finite_nonpromotion",
            role="route_decision",
            readiness="guard_validated",
            claim="A finite Q209 computation is calibration or falsification, not the all-stage theorem.",
            formula="finite Q209 != uniform successor induction",
            proof_boundary="No cofinal conclusion follows from one more rectangle.",
        ),
        GateRow(
            id="q209_gate_15_ray_route",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The ray-aligned schedule closes every new strip and isolates one uniform outer-descendant Xi estimate.",
            formula=exact["ray_aligned_route"]["only_open_successor_region"],
            proof_boundary="The a priori relative Xi estimate remains RH-level and unproved.",
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact Q208-base/Q209 successor-shell coverage gate "
            "with two explicit open regions"
        ),
        "proof_boundary": (
            "This artifact corrects the finite handoff, proves the exact "
            "Q209 geometry, orientation, area decomposition, and coverage "
            "classification, and composes only already certified source "
            "domains. It does not certify either open Q209 region, compute "
            "a Q209 winding, prove a uniform successor theorem, prove "
            "Lambda<=0, prove PF-infinity, prove RH, or establish any "
            "Clay-prize conclusion."
        ),
        "source_audit": audit,
        "exact": exact,
        "coverage": [asdict(row) for row in coverage],
        "boundary": [asdict(row) for row in boundary],
        "rows": [asdict(row) for row in rows],
        "summary": {
            "proved_base_index": 208,
            "first_unresolved_linear_stage": 209,
            "coverage_regions": len(coverage),
            "covered_regions": sum(
                row.status == "covered_exact" for row in coverage
            ),
            "open_shell_regions": sum(
                row.status == "open_xi_antecedent" for row in coverage
            ),
            "boundary_arcs": len(boundary),
            "open_boundary_arcs": sum(
                row.status == "open_xi_antecedent" for row in boundary
            ),
            "outer_old_half_unit_cells": 416,
            "new_right_strip_half_unit_columns": 2,
            "q209_certified": False,
            "ray_aligned_new_strips_open": 0,
            "ray_aligned_old_collar_antecedents": 1,
            "gate_rows": len(rows),
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    shell = exact["low_shell_decomposition"]
    orientation = exact["boundary_orientation"]
    conditional = exact["conditional_successor"]
    guards = exact["noninheritance_guards"]
    ray = exact["ray_aligned_route"]
    coverage_lines = [
        (
            f"{row['id']}: {row['status']}; {row['region']}; "
            f"{row['obligation']}"
        )
        for row in payload["coverage"]
    ]
    boundary_lines = [
        (
            f"{row['id']}: {row['edge']}; {row['arc']}; "
            f"{row['status']}; {row['reason']}"
        )
        for row in payload["boundary"]
    ]
    return "\n".join(
        [
            "# Q208 Base and Q209 Successor-Shell Coverage Gate",
            "",
            "Date: 2026-07-28",
            "",
            "Status: exact finite-base correction and Q209 shell coverage",
            "classification. This is not a proof of Q209, the uniform",
            "successor theorem, Lambda<=0, or RH.",
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
            "validated Q208-base/Q209 shell coverage gate: 15 rows, "
            "5 coverage regions, 2 open shell regions, 7 oriented boundary "
            "arcs, 416 outer old half-cells, 2 new strip half-columns, "
            "Q209 not certified",
            "```",
            "",
            "## Corrected Finite Base",
            "",
            payload["source_audit"]["historical_handoff_correction"],
            "",
            "The current theorem-level finite base is",
            "",
            "```text",
            exact["rectangles"]["Q_208"],
            "Q208 is contact-free and has zero first-jet winding.",
            "```",
            "",
            "The first unresolved rectangle in the linear schedule is",
            "",
            "```text",
            exact["rectangles"]["Q_209"],
            "```",
            "",
            "## Exact Q209 Geometry",
            "",
            "```text",
            exact["linear_schedule"]["formula"],
            (
                "t_208="
                + exact["linear_schedule"]["t_208"]
                + ", t_209="
                + exact["linear_schedule"]["t_209"]
            ),
            "delta_208=" + exact["linear_schedule"]["delta_208"],
            (
                "R_208="
                + exact["linear_schedule"]["R_208"]
                + ", R_209="
                + exact["linear_schedule"]["R_209"]
            ),
            shell["identity"],
            shell["refined_identity"],
            (
                "area(P209)="
                + shell["area_P_209"]
                + " and the four rational piece areas sum exactly to it"
            ),
            "```",
            "",
            "The standard `(x,t)` boundary orientation is",
            "",
            "```text",
            orientation["bottom"],
            orientation["right"],
            orientation["top"],
            orientation["left"],
            "twice signed area=" + orientation["twice_signed_area"] + ">0",
            "```",
            "",
            "## Coverage Table",
            "",
            "```text",
            *coverage_lines,
            "```",
            "",
            "Exactly two low-time shell regions remain open: the outer",
            "old descendant collar and the new right strip.",
            "",
            "## Oriented Boundary Audit",
            "",
            "```text",
            *boundary_lines,
            "```",
            "",
            "The three open boundary rows are traces of the same two open",
            "two-dimensional shell regions; they are not three independent",
            "analytic obligations.",
            "",
            "## Conditional Q209 Theorem",
            "",
            "```text",
            conditional["outer_transport"],
            conditional["new_strip_half_plane"],
            "```",
            "",
            conditional["conclusion"],
            "",
            "Neither antecedent is currently proved.",
            "",
            "## Noninheritance Guards",
            "",
            guards["time_enclosures"],
            "",
            guards["space_translation"],
            "",
            guards["finite_nonpromotion"],
            "",
            "## Proof-Level Route",
            "",
            "The linear Q209 instance is useful only as a bounded calibration.",
            "The preferred cofinal geometry is",
            "",
            "```text",
            ray["schedule"],
            ray["new_strip_identity"],
            ray["new_strip_status"],
            ray["only_open_successor_region"],
            ray["scale"],
            "```",
            "",
            ray["open_xi_antecedent"],
            "",
            "Pi provenance:",
            "",
            ray["pi_provenance"],
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
        "wrote Q208-base/Q209 shell coverage gate: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
