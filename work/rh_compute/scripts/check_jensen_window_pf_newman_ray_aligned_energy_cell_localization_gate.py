#!/usr/bin/env python3
"""Independently validate the ray-aligned energy localization gate."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_ray_aligned_energy_cell_localization_gate"
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
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
EXPECTED_IDS = [
    "raecl_01_partition",
    "raecl_02_local_balance",
    "raecl_03_cell_nonnegativity",
    "raecl_04_two_sided_contact_floor",
    "raecl_05_certified_threshold",
    "raecl_06_one_sided_threshold",
    "raecl_07_local_no_contact",
    "raecl_08_flux_telescoping",
    "raecl_09_global_minimum_criterion",
    "raecl_10_sum_threshold_guard",
    "raecl_11_anchor_vertex_cover",
    "raecl_12_outer_anchor_obstruction",
    "raecl_13_flux_localization_tradeoff",
    "raecl_14_compact_ownership",
    "raecl_15_dominant_ownership",
    "raecl_16_oscillatory_effectivity",
    "raecl_17_critical_c1_target",
    "raecl_18_finite_phase_ownership",
    "raecl_19_successor_and_snapshot_guard",
    "raecl_20_route_decision",
]


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, label: str, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing {label}: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid {label}: {exc}")
        return {}


def validate_contact_minimization(issues: list[str]) -> None:
    e_left, e_right, h, u = sp.symbols(
        "e_left e_right h u", positive=True
    )
    lhs = e_left**2 / u + e_right**2 / (h - u)
    floor = (e_left + e_right) ** 2 / h
    expected = (
        (e_left * (h - u) - e_right * u) ** 2
        / (h * u * (h - u))
    )
    if sp.simplify(lhs - floor - expected) != 0:
        issues.append("two-sided contact-floor identity failed")
    minimizer = sp.simplify(h * e_left / (e_left + e_right))
    if sp.simplify(sp.diff(lhs, u).subs(u, minimizer)) != 0:
        issues.append("contact-floor minimizer failed")
    if sp.simplify(lhs.subs(u, minimizer) - floor) != 0:
        issues.append("contact-floor minimum failed")


def validate_telescoping(issues: list[str]) -> None:
    currents = sp.symbols("J0:8")
    telescoped = sum(
        currents[index] - currents[index - 1]
        for index in range(1, len(currents))
    )
    if sp.expand(telescoped) != currents[-1] - currents[0]:
        issues.append("internal-current telescoping failed")


def validate_budget_logic(issues: list[str]) -> None:
    thresholds = [1, 100]
    local_budgets = [1, 0]
    total = sum(local_budgets)
    if not total < sum(thresholds):
        issues.append("sum-threshold counterexample inequality failed")
    if not local_budgets[0] >= thresholds[0]:
        issues.append("sum-threshold counterexample is not contact-compatible")
    if not total >= min(thresholds):
        issues.append("minimum-threshold criterion logic failed")

    owned_two_cells = {0, 2}
    edges_two = [{0, 1}, {1, 2}]
    if not all(edge & owned_two_cells for edge in edges_two):
        issues.append("two-cell outer-anchor cover failed")
    owned_three_cells = {0, 3}
    edges_three = [{0, 1}, {1, 2}, {2, 3}]
    if all(edge & owned_three_cells for edge in edges_three):
        issues.append("three-cell outer-anchor obstruction failed")


def validate_sources(
    stored: dict, source_payloads: dict[str, dict], issues: list[str]
) -> None:
    stored_hashes = (
        stored.get("source_audit", {}).get("source_sha256", {})
    )
    expected_hashes = {
        key: file_hash(path) for key, path in SOURCES.items()
        if path.is_file()
    }
    if stored_hashes != expected_hashes:
        issues.append("source hashes drifted")

    energy = source_payloads["energy_current"]
    if energy.get("exact", {}).get("fixed_slice_balance") != (
        "d/dt integral_a^b E_pf dx+2[J_pf(t,b)-J_pf(t,a)]"
        "=2 integral_a^b B_pf dx"
    ):
        issues.append("energy source balance drifted")
    if energy.get("summary", {}).get("pointwise_contact_exclusion"):
        issues.append("energy source was unexpectedly promoted")

    endpoint = source_payloads["endpoint_margin"]
    if endpoint.get("summary", {}).get("left_raw_margins") != 1:
        issues.append("left outer-anchor count drifted")
    if endpoint.get("summary", {}).get("right_normalized_margins") != 1:
        issues.append("right outer-anchor count drifted")

    dominant = source_payloads["dominant_global"]
    if "t*L>=25" not in dominant.get("exact", {}).get("region", ""):
        issues.append("dominant source region drifted")

    oscillatory = source_payloads["oscillatory_handoff"]
    if "there exists L_epsilon" not in (
        oscillatory.get("exact", {}).get("theorem", "")
    ):
        issues.append("oscillatory existential threshold drifted")
    if "does not provide a practical L_epsilon" not in (
        oscillatory.get("proof_boundary", "")
    ):
        issues.append("oscillatory effectivity guard drifted")

    critical = source_payloads["critical_c1"]
    if critical.get("exact", {}).get("remaining_target") != (
        "Prove T_L[J]>32000000*exp(-3L/2) "
        "for the corrected finite main"
    ):
        issues.append("critical finite-main target drifted")

    successor = source_payloads["successor_lemma"]
    if (
        successor.get("exact", {})
        .get("transport", {})
        .get("strict_gate")
        != "delta_j*M_(j,k)<d_(j,k)"
    ):
        issues.append("all-j successor strict gate drifted")
    if "all-j antecedent open" not in successor.get("status", ""):
        issues.append("all-j successor boundary drifted")

    finite = source_payloads["finite_transport"].get("summary", {})
    if finite.get("finite_successor_instances") != 1:
        issues.append("finite successor count drifted")
    if finite.get("certified_transport_panels") != 525:
        issues.append("finite successor panel count drifted")

    bottom = source_payloads["q208_bottom"]
    if bottom.get("contract", {}).get("time") != "1/1040":
        issues.append("Q208 bottom time drifted")
    if bottom.get("summary", {}).get("certified_leaf_cells") != 492:
        issues.append("Q208 bottom cell count drifted")
    if bottom.get("summary", {}).get("unresolved_boxes") != 0:
        issues.append("Q208 bottom is incomplete")

    top = source_payloads["q208_top"]
    if top.get("contract", {}).get("time") != "1/5":
        issues.append("Q208 top time drifted")
    if top.get("summary", {}).get("certified_leaf_cells") != 492:
        issues.append("Q208 top cell count drifted")
    if top.get("summary", {}).get("unresolved_boxes") != 0:
        issues.append("Q208 top is incomplete")


def validate() -> list[str]:
    issues: list[str] = []
    stored = load_json(RESULT, "localization result", issues)
    source_payloads = {
        key: load_json(path, f"{key} source", issues)
        for key, path in SOURCES.items()
    }
    if issues:
        return issues

    if stored.get("kind") != STEM:
        issues.append("result kind drifted")
    if stored.get("date") != "2026-07-29":
        issues.append("result date drifted")

    rows = stored.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or order drifted")
    expected_summary = {
        "rows": 20,
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
    }
    if stored.get("summary") != expected_summary:
        issues.append("summary drifted")

    exact = stored.get("exact", {})
    if exact.get("local_balance", {}).get("identity") != (
        "D_k(t)=2 integral_(I_k) B_pf(t,x) dx>=0"
    ):
        issues.append("stored local balance drifted")
    if exact.get("cell_contact_floor", {}).get("certified") != (
        "V_pf(t,c)=0 in I_k => D_k(t)>=T_k(t), "
        "T_k=2[m_(k-1)(t)+m_k(t)]^2/h_k"
    ):
        issues.append("stored certified threshold drifted")
    if "min_(1<=k<=N) T_k(t)" not in (
        exact.get("global_minimum_criterion", {}).get("statement", "")
    ):
        issues.append("stored global minimum criterion drifted")
    if exact.get("sum_threshold_guard", {}).get(
        "exact_budget_counterexample"
    ) != (
        "T_1=1, T_2=100, D_1=1, D_2=0 gives "
        "D_total=1<T_1+T_2=101 while D_1>=T_1."
    ):
        issues.append("stored sum-threshold guard drifted")
    if exact.get("anchor_graph", {}).get("current_uniform_set") != (
        "{x_0=38,x_N=R_j}"
    ):
        issues.append("stored outer-anchor set drifted")
    if "T_L[J]>32000000 exp(-3L/2)" not in (
        exact.get("open_handoff", "")
    ):
        issues.append("stored corrected-main handoff drifted")

    boundary = stored.get("proof_boundary", "")
    for required in (
        "does not prove a cofinal interior anchor registry",
        "corrected finite-main C1 lower bound",
        "Lambda<=0",
        "RH",
    ):
        if required not in boundary:
            issues.append(f"proof boundary missing: {required}")

    validate_contact_minimization(issues)
    validate_telescoping(issues)
    validate_budget_logic(issues)
    validate_sources(stored, source_payloads, issues)
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated ray-aligned energy cell-localization gate: 20 rows, "
        "2 exact local balances, 2 exact contact floors, "
        "2 conditional no-contact criteria, 2 outer anchors, "
        "0 cofinal internal anchor registries, "
        "1 finite successor instance, 2 nonpromotion guards, "
        "0 pointwise contact exclusions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
