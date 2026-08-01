#!/usr/bin/env python3
"""Localize a hypothetical positive Newman boundary away from t=0."""

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
    "jensen_window_pf_newman_positive_boundary_delta_localization_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_newman_positive_boundary_delta_localization_gate.md"
)
BOUNDARY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
)
DOMINANT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_polymath15_"
    "dominant_saddle_global_ray_certificate.json"
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


def build_exact() -> dict:
    x, time, scale = sp.symbols("x time scale", real=True, positive=True)
    model = x**2 - 2 * time
    if sp.simplify(sp.diff(model, time) + sp.diff(model, x, 2)) != 0:
        raise RuntimeError("quadratic backward-heat identity failed")
    first = sp.diff(model, x)
    laguerre = sp.expand(first**2 - model * sp.diff(model, x, 2))
    if laguerre != 2 * x**2 + 4 * time:
        raise RuntimeError("quadratic Laguerre identity failed")
    first_jet = sp.expand(model**2 + (first / scale) ** 2)
    origin_jet = sp.simplify(first_jet.subs(x, 0))
    if origin_jet != 4 * time**2:
        raise RuntimeError("quadratic origin first-jet identity failed")
    root_jet = sp.expand((first / scale) ** 2).subs(x**2, 2 * time)
    if sp.simplify(root_jet - 8 * time / scale**2) != 0:
        raise RuntimeError("quadratic root first-jet identity failed")

    return {
        "boundary_attainment": (
            "If Lambda>0, then 0<Lambda<=1/5 and H_Lambda has a finite "
            "real multiple zero"
        ),
        "positive_time_equivalence": (
            "Lambda<=0 iff H_t has only simple zeros for every 0<t<=1/5"
        ),
        "delta_contradiction": (
            "To contradict Lambda>0 it is enough, for every rational "
            "delta in (0,1/5], to exclude common zeros of H_t and H_t' "
            "on delta<=t<=1/5; choose any rational delta<Lambda"
        ),
        "dominant_localization": (
            "For fixed delta>0, t>=delta and "
            "L=log(|x|/(4*pi))>=25/delta imply tL>=25 and L>=50, "
            "so the dominant-saddle theorem gives "
            "H_t'(x)^2-H_t(x)H_t''(x)>0"
        ),
        "compact_rectangle": (
            "Every multiple zero with delta<=t<=1/5 lies in "
            "|x|<R_delta, R_delta=4*pi*exp(25/delta)"
        ),
        "cofinal_sequence": (
            "For delta_j=1/(5j), j>=1, Lambda<=0 iff "
            "(H_t,H_t')!=(0,0) on "
            "[delta_j,1/5]x[-4*pi*exp(125j),4*pi*exp(125j)] "
            "for every j"
        ),
        "compact_minimum": (
            "For each fixed j, no common zero on the compact rectangle is "
            "equivalent to min(H_t(x)^2+H_t'(x)^2)>0 there"
        ),
        "corrected_delta_target": (
            "A sufficient nonuniform route is: for each fixed delta>0, "
            "prove the corrected C1 finite-main inequality only on "
            "delta<=t<=1/5 and 38<|x|<4*pi*exp(25/delta), with constants "
            "allowed to deteriorate as delta decreases"
        ),
        "quadratic_model": (
            "G_t(x)=x^2-2t solves partial_t G=-partial_x^2 G; "
            "for t>0 it has the simple real zeros +/-sqrt(2t), at t=0 "
            "it has a double zero, and for t<0 its zeros are nonreal"
        ),
        "quadratic_margins": (
            "L[G_t]=G_t'^2-G_t G_t''=2x^2+4t>0 for t>0, while "
            "G_t(0)^2+(G_t'(0)/A)^2=4t^2 and at either root the same "
            "first-jet quantity is 8t/A^2; both margins tend to zero as "
            "t decreases to zero"
        ),
        "route_reclassification": (
            "A t-independent positive first-jet floor down to t=0 is a "
            "valid stronger sufficient target, but it is not logically "
            "necessary for Lambda<=0 or positive-time simplicity"
        ),
        "proof_handoff": (
            "Use delta-dependent compact or arithmetic transversality "
            "estimates, or keep the uniform scaled small-ball theorem as an "
            "optional stronger route; do not require endpoint simplicity "
            "without stating that extra burden"
        ),
    }


def source_audit() -> dict:
    boundary = json.loads(BOUNDARY_SOURCE.read_text(encoding="utf-8"))
    dominant = json.loads(DOMINANT_SOURCE.read_text(encoding="utf-8"))
    boundary_text = json.dumps(boundary)
    dominant_text = json.dumps(dominant)
    for marker in (
        "If Lambda>0",
        "H_Lambda has a finite real multiple zero",
        "Lambda<=0 if and only if H_t has only simple zeros",
    ):
        if marker not in boundary_text:
            raise RuntimeError(f"boundary source marker missing: {marker}")
    for marker in (
        "t*log(x/(4*pi))>=25",
        "L_t(x)>0",
        "0<t<=1/2",
    ):
        if marker not in dominant_text:
            raise RuntimeError(f"dominant source marker missing: {marker}")
    return {
        "boundary_kind": boundary["kind"],
        "boundary_date": boundary["date"],
        "dominant_kind": dominant["kind"],
        "dominant_date": dominant["date"],
        "dominant_tL_min": dominant["parameters"]["tL_min"],
    }


def cofinal_examples() -> list[dict]:
    return [
        {
            "j": index,
            "delta_j": f"1/{5 * index}",
            "L_cutoff": str(125 * index),
            "radius": f"4*pi*exp({125 * index})",
        }
        for index in range(1, 5)
    ]


def build_rows(exact: dict, audit: dict) -> list[GateRow]:
    return [
        GateRow(
            "npbdlg_01_boundary_attainment",
            "published_input",
            "ready_to_apply",
            "A hypothetical positive Newman boundary is attained at a finite real multiple zero.",
            exact["boundary_attainment"],
            "Uses the separately validated Polymath-15 compactness handoff.",
            audit,
        ),
        GateRow(
            "npbdlg_02_positive_time_equivalence",
            "exact_equivalence",
            "ready_to_apply",
            "The Newman direction is exactly a positive-time simplicity problem.",
            exact["positive_time_equivalence"],
            "Uses the published bound Lambda<=1/5 and simplicity above Lambda.",
        ),
        GateRow(
            "npbdlg_03_delta_contradiction",
            "exact_reduction",
            "ready_to_apply",
            "Uniformity as t tends to zero is not needed in the contradiction argument.",
            exact["delta_contradiction"],
            "The family must hold for a cofinal set of positive deltas.",
        ),
        GateRow(
            "npbdlg_04_dominant_localization",
            "exact_composition",
            "ready_to_apply",
            "For each fixed positive lower time, the dominant-saddle theorem makes the remaining frequency range compact.",
            f"{exact['dominant_localization']}; {exact['compact_rectangle']}",
            "The compact radius grows exponentially as delta decreases.",
        ),
        GateRow(
            "npbdlg_05_cofinal_sequence",
            "exact_equivalence",
            "ready_to_apply",
            "One explicit countable sequence of compact strips is equivalent to the Newman direction.",
            exact["cofinal_sequence"],
            "This is a logical equivalence, not a certificate for any strip.",
            cofinal_examples(),
        ),
        GateRow(
            "npbdlg_06_compact_minimum",
            "exact_reduction",
            "ready_to_apply",
            "Each fixed strip has an ordinary compact nonvanishing formulation.",
            exact["compact_minimum"],
            "No effective positive minimum is supplied.",
        ),
        GateRow(
            "npbdlg_07_delta_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "A delta-dependent corrected first-jet certificate would suffice without endpoint-uniform margins.",
            exact["corrected_delta_target"],
            "Open for every cofinal delta; the displayed radii are far beyond current compact certificates.",
        ),
        GateRow(
            "npbdlg_08_quadratic_calibration",
            "exact_countermodel",
            "guard_validated",
            "The backward heat equation permits Lambda=0 with positive-time simplicity and a double endpoint zero.",
            exact["quadratic_model"],
            "Generic heat-flow calibration, not a model of Xi arithmetic.",
        ),
        GateRow(
            "npbdlg_09_uniform_floor_guard",
            "nonpromotion_gate",
            "guard_validated",
            "A positive first-jet floor uniform down to t=0 is stronger than the required positive-time simplicity statement.",
            f"{exact['quadratic_margins']}; {exact['route_reclassification']}",
            "Does not show that a stronger Xi-specific uniform theorem is false.",
        ),
        GateRow(
            "npbdlg_10_proof_handoff",
            "proof_guard",
            "active",
            "The theorem search may use nonuniform positive-time strips or a clearly labelled stronger endpoint-uniform route.",
            exact["proof_handoff"],
            "No compact strip beyond |x|<=38 is closed here.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    audit = source_audit()
    rows = build_rows(exact, audit)
    return {
        "kind": (
            "jensen_window_pf_newman_positive_boundary_"
            "delta_localization_gate"
        ),
        "date": "2026-07-24",
        "status": (
            "exact delta-localized positive-boundary reduction and "
            "endpoint-uniformity route guard; the compact strip family "
            "remains open and this is not a proof of Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "The gate proves a logical localization and an exact quadratic "
            "countermodel only. It does not certify any new Xi compact "
            "rectangle, positive-time simplicity, Lambda<=0, or RH."
        ),
        "sources": [
            str(BOUNDARY_SOURCE.relative_to(REPO_ROOT)).replace("\\", "/"),
            str(DOMINANT_SOURCE.relative_to(REPO_ROOT)).replace("\\", "/"),
        ],
        "source_audit": audit,
        "exact": exact,
        "cofinal_examples": cofinal_examples(),
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    examples = artifact["cofinal_examples"]
    lines = [
        "# Newman Positive-Boundary Delta Localization Gate",
        "",
        "Date: 2026-07-24",
        "",
        "Status: exact positive-boundary localization and route guard. The",
        "delta-dependent compact-strip family remains open; this is not a proof of `Lambda <= 0` or RH.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_positive_boundary_delta_localization_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_positive_boundary_delta_localization_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_positive_boundary_delta_localization_gate.py",
        "```",
        "",
        "## Exact Localization",
        "",
        "```text",
        exact["boundary_attainment"],
        exact["positive_time_equivalence"],
        exact["delta_contradiction"],
        "```",
        "",
        "For a fixed positive lower time, the complete dominant ray gives",
        "",
        "```text",
        exact["dominant_localization"],
        exact["compact_rectangle"],
        "```",
        "",
        "Therefore one explicit cofinal family is",
        "",
        "```text",
        exact["cofinal_sequence"],
        "```",
        "",
        "| j | delta_j | L cutoff | compact radius |",
        "|---:|---:|---:|---:|",
    ]
    for row in examples:
        lines.append(
            f"| {row['j']} | {row['delta_j']} | {row['L_cutoff']} | "
            f"{row['radius']} |"
        )
    lines.extend(
        [
            "",
            "These rectangles grow rapidly and are not currently certified.",
            "The reduction removes the need for one margin uniform all the way",
            "to `t=0`; it does not make the compact arithmetic problem easy.",
            "",
            "## Endpoint-Uniformity Guard",
            "",
            "```text",
            exact["quadratic_model"],
            exact["quadratic_margins"],
            "```",
            "",
            "This model has the correct backward-heat collision geometry:",
            "positive-time zeros are real and simple, but their separation and",
            "first-jet margins vanish at the endpoint. Hence",
            "",
            "```text",
            exact["route_reclassification"],
            "```",
            "",
            "The current uniform corrected small-ball target remains a valid",
            "stronger route. It must not be presented as logically necessary",
            "unless endpoint simplicity is explicitly added to the burden.",
            "",
            "## Live Handoff",
            "",
            "```text",
            exact["corrected_delta_target"],
            exact["proof_handoff"],
            "```",
            "",
        ]
    )
    return "\n".join(lines)


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
        "built Newman positive-boundary delta-localization gate: "
        f"{len(artifact['rows'])} rows, 5 exact reductions, "
        "1 cofinal-strip criterion, 1 quadratic endpoint-uniformity "
        "countermodel, 1 delta-dependent open target"
    )


if __name__ == "__main__":
    main()
