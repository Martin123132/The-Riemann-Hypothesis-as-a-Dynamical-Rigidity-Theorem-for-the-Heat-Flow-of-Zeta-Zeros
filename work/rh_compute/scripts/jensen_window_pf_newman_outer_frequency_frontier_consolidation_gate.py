#!/usr/bin/env python3
"""Consolidate the proved compact base and the live outer-frequency frontier."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_outer_frequency_frontier_consolidation_gate"
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
DATE = "2026-08-04"

SOURCES = {
    "compact_transversality": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_theta_"
        "compact_transversality_interval_certificate.json"
    ),
    "q208_closed_boundary": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_q208_closed_boundary_winding_certificate.json"
    ),
    "q209_shell_coverage": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_q208_base_q209_shell_coverage_gate.json"
    ),
    "one_sided_phase_bridge": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.json"
    ),
    "complete_theta_spatial_band": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_theta_"
        "complete_spatial_tile_degree_certificate.json"
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_path(key: str) -> Path:
    return REPO_ROOT / SOURCES[key]


def load_source(key: str) -> dict:
    path = source_path(key)
    return json.loads(path.read_text(encoding="utf-8"))


def find_formula(payload: dict, needle: str) -> str:
    for row in payload.get("rows", []):
        formula = str(row.get("formula", ""))
        if needle in formula:
            return formula
    raise RuntimeError(f"source formula missing: {needle}")


def build_payload() -> dict:
    compact = load_source("compact_transversality")
    q208 = load_source("q208_closed_boundary")
    q209 = load_source("q209_shell_coverage")
    phase = load_source("one_sided_phase_bridge")
    band = load_source("complete_theta_spatial_band")

    require(
        compact["status"].startswith(
            "rigorous Arb/Taylor contact exclusion for |x|<=38"
        ),
        "compact theorem status drifted",
    )
    require(
        compact["exact"]["proved_rectangle"]
        == "0<=t<=1/5 and 1/4<=x<=38",
        "compact rectangle drifted",
    )
    require(
        "|x|<=1/4" in compact["exact"]["origin_overlap"],
        "compact origin composition drifted",
    )
    require(
        compact["certificate"]["certified_boxes"] == 1900
        and compact["certificate"]["unresolved_boxes"] == 0,
        "compact interval partition is incomplete",
    )
    require(
        q208["conclusion"]
        == "Q_208 is contact-free and has zero first-jet winding; "
        "therefore Q_1 through Q_208 are certified by containment.",
        "Q208 conclusion drifted",
    )
    require(
        q209["summary"]["first_unresolved_linear_stage"] == 209,
        "first unresolved linear stage drifted",
    )
    require(
        q209["summary"]["open_shell_regions"] == 2
        and q209["summary"]["q209_certified"] is False,
        "Q209 open-shell boundary drifted",
    )
    phase_target = find_formula(
        phase,
        "(B_t(x),B_t'(x))!=(0,0) for every x>38",
    )
    require(
        band["combined_band"]["time"] == ["0", "0.45"]
        and band["combined_band"]["x"] == ["135", "136.5"],
        "complete-theta spatial band drifted",
    )
    require(
        "no common zero" in band["theorem"],
        "complete-theta band theorem drifted",
    )

    critical_time_high = Fraction(1, 5)
    band_time_high = Fraction(9, 20)
    compact_x_high = Fraction(38)
    band_x_low = Fraction(135)
    require(
        critical_time_high <= band_time_high,
        "critical time window is not contained in the local band time window",
    )
    require(
        compact_x_high < band_x_low,
        "local band no longer lies strictly beyond the compact core",
    )

    source_audit = {}
    payloads = {
        "compact_transversality": compact,
        "q208_closed_boundary": q208,
        "q209_shell_coverage": q209,
        "one_sided_phase_bridge": phase,
        "complete_theta_spatial_band": band,
    }
    for key, payload in payloads.items():
        path = source_path(key)
        source_audit[key] = {
            "path": SOURCES[key],
            "sha256": sha256_path(path),
            "kind": payload["kind"],
            "status": payload["status"],
        }

    rows = [
        {
            "id": "ofc_01_compact_core_closed",
            "readiness": "ready_to_apply",
            "claim": (
                "The complete Newman theta first jet has no common zero on "
                "0<=t<=1/5 and |x|<=38. The origin collar is exact and the "
                "remaining compact rectangle has 1900 certified Arb/Taylor boxes."
            ),
        },
        {
            "id": "ofc_02_q208_finite_base",
            "readiness": "ready_to_apply",
            "claim": (
                "Q208=[1/1040,1/4]x[0,246] is contact-free and has zero "
                "first-jet winding; Q1 through Q208 follow by containment."
            ),
        },
        {
            "id": "ofc_03_q209_exact_frontier",
            "readiness": "ready_to_apply",
            "claim": (
                "Q209 is the first unresolved linear stage. Its low shell has "
                "exactly two open regions: the outer old collar and the new right strip."
            ),
        },
        {
            "id": "ofc_04_direct_outer_reduction",
            "readiness": "ready_to_apply",
            "claim": phase_target,
        },
        {
            "id": "ofc_05_local_band_placement",
            "readiness": "ready_to_apply",
            "claim": (
                "The complete-theta theorem on 0<=t<=0.45 and "
                "135<=x<=136.5 lies inside the direct outer-frequency target "
                "during 0<=t<=1/5, but it covers only one finite spatial band."
            ),
        },
        {
            "id": "ofc_06_compact_rescout_redundant",
            "readiness": "ready_to_apply",
            "claim": (
                "A new boundary scout for [0,1/5]x[0,38] is logically "
                "subsumed by the existing full interior transversality theorem."
            ),
        },
        {
            "id": "ofc_07_cofinal_successor_target",
            "readiness": "open",
            "claim": q209["exact"]["ray_aligned_route"]["open_xi_antecedent"],
        },
        {
            "id": "ofc_08_global_outer_target",
            "readiness": "open",
            "claim": (
                "Prove complete first-jet nonvanishing for every x>38 and "
                "0<t<=1/5, either by the cofinal boundary-successor theorem "
                "or by the Xi-specific phase-critical small-ball theorem."
            ),
        },
    ]

    return {
        "kind": STEM,
        "schema_version": 1,
        "date": DATE,
        "status": (
            "exact compact-base and outer-frequency frontier consolidation; "
            "the outer Xi joint-avoidance theorem remains open"
        ),
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "source_audit": source_audit,
        "exact": {
            "compact_core": {
                "domain": "0<=t<=1/5 and |x|<=38",
                "certified_boxes": 1900,
                "unresolved_boxes": 0,
                "contact_free": True,
                "new_boundary_scout_required": False,
            },
            "finite_cofinal_base": {
                "proved_through": 208,
                "first_unresolved_linear_stage": 209,
                "q208": "[1/1040,1/4]x[0,246]",
                "q209": "[1/1045,1/4]x[0,247]",
                "q209_open_shell_regions": 2,
            },
            "local_band": {
                "domain": "0<=t<=0.45 and 135<=x<=136.5",
                "critical_time_window_contained": True,
                "strictly_outside_compact_core": True,
                "global_outer_theorem": False,
            },
            "direct_remaining_domain": "0<t<=1/5 and x>38",
            "direct_phase_target": phase_target,
            "route_split": {
                "minimal_topological_route": (
                    "prove a source-level relative first-jet transport theorem "
                    "on every ray-aligned outer descendant collar"
                ),
                "stronger_direct_route": (
                    "prove Xi-specific phase-critical joint avoidance for "
                    "B_t and B_t' on the complete outer domain"
                ),
                "retired_next_step": (
                    "rescouting the already certified compact rectangle or "
                    "continuing an indefinite fixed-width tile ladder"
                ),
            },
        },
        "rows": rows,
        "summary": {
            "rows": len(rows),
            "compact_certified_boxes": 1900,
            "compact_rescout_required": False,
            "proved_cofinal_base": 208,
            "first_unresolved_linear_stage": 209,
            "q209_open_shell_regions": 2,
            "ready_rows": sum(row["readiness"] == "ready_to_apply" for row in rows),
            "open_rows": sum(row["readiness"] == "open" for row in rows),
        },
        "next_action": (
            "Do not build the proposed [0,1/5]x[0,38] boundary scout. "
            "Use the replayed compact theorem as a closed parent. Before another "
            "large expansion, derive or falsify an a priori source-level relative "
            "first-jet transport inequality on the ray-aligned outer collar; keep "
            "the direct phase-critical-value theorem as the independent alternative."
        ),
        "proof_boundary": (
            "This gate consolidates already proved domains and identifies the exact "
            "remaining outer-frequency obligations. It proves no new outer contact "
            "exclusion, no Q209 shell antecedent, no cofinal successor theorem, no "
            "phase-critical small-ball inequality, Lambda<=0, RH, or prize-level result."
        ),
    }


def success_line(payload: dict) -> str:
    summary = payload["summary"]
    return (
        "validated Newman outer-frequency frontier consolidation gate: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['compact_certified_boxes']} compact boxes, Q208 base, "
        f"Q{summary['first_unresolved_linear_stage']} first unresolved, "
        f"{summary['q209_open_shell_regions']} open shell regions, "
        "no compact rescout"
    )


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# Newman Outer-Frequency Frontier Consolidation Gate",
        "",
        f"Date: {DATE}",
        "",
        "Status: exact dependency and route consolidation. The outer Xi theorem remains open; this is not a proof of RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Closed Compact Base",
        "",
        "The full interval theorem already proves",
        "",
        "```text",
        "(H_t(x),H_t'(x))!=(0,0) for 0<=t<=1/5 and |x|<=38.",
        "1900 certified Arb/Taylor boxes, 0 unresolved boxes.",
        "```",
        "",
        "A new boundary scout on this rectangle would be weaker and redundant.",
        "",
        "## Finite Cofinal Frontier",
        "",
        "```text",
        f"proved finite base: Q_{exact['finite_cofinal_base']['proved_through']}",
        f"first unresolved linear stage: Q_{exact['finite_cofinal_base']['first_unresolved_linear_stage']}",
        f"open Q209 shell regions: {exact['finite_cofinal_base']['q209_open_shell_regions']}",
        "```",
        "",
        "Those regions are the old outer descendant collar and the new right strip.",
        "",
        "## Local Band",
        "",
        "The new complete-theta band `0<=t<=0.45`, `135<=x<=136.5` is a rigorous",
        "outer-frequency calibration. It is not a global outer theorem and should not",
        "be extended indefinitely by fixed-width tiling.",
        "",
        "## Live Outer Target",
        "",
        "```text",
        exact["direct_remaining_domain"],
        exact["direct_phase_target"],
        "```",
        "",
        "The minimal topological route is a source-level relative first-jet transport",
        "theorem on every ray-aligned outer descendant collar. The stronger alternative",
        "is the direct Xi-specific phase-critical-value theorem.",
        "",
        "## Proof Boundary",
        "",
        payload["proof_boundary"],
        "",
        success_line(payload),
        "",
    ]
    return "\n".join(lines)


def write_text_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def write_json_atomic(path: Path, payload: dict) -> None:
    write_text_atomic(path, json.dumps(payload, indent=2) + "\n")


def main() -> int:
    payload = build_payload()
    write_json_atomic(RESULT, payload)
    write_text_atomic(NOTE, render_note(payload))
    print(success_line(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
