#!/usr/bin/env python3
"""Reconcile the closed C1 remainder with the live outer arithmetic target."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_remainder_supersession_frontier_gate"
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
DATE = "2026-08-04"

SOURCES = {
    "old_peeling_contract": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_RS_C1_"
        "endpoint_peeling_contract.json"
    ),
    "global_first_order_remainder": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_first_order_"
        "global_remainder_certificate.json"
    ),
    "target_reconciliation": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_"
        "target_reconciliation_gate.json"
    ),
    "signed_contact": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_first_order_"
        "signed_contact_reduction.json"
    ),
    "residual_handoff": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "residual_placement_c1_cumulative_handoff_gate.json"
    ),
    "symmetric_outer_pairing": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_full_support_"
        "ideal_cubic_symmetric_outer_pairing_gate.json"
    ),
    "ordinary_sqrt_transfer": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_theta_forward_sqrt_to_corrected_"
        "rs_C1_transfer_gate.json"
    ),
    "outer_frontier": (
        "work/rh_compute/results/"
        "jensen_window_pf_newman_outer_frequency_frontier_"
        "consolidation_gate.json"
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
    return json.loads(source_path(key).read_text(encoding="utf-8"))


def contains(payload: dict, needle: str) -> bool:
    return needle in json.dumps(payload, sort_keys=True)


def build_payload() -> dict:
    old = load_source("old_peeling_contract")
    remainder = load_source("global_first_order_remainder")
    reconciliation = load_source("target_reconciliation")
    contact = load_source("signed_contact")
    handoff = load_source("residual_handoff")
    pairing = load_source("symmetric_outer_pairing")
    ordinary = load_source("ordinary_sqrt_transfer")
    frontier = load_source("outer_frontier")

    require(
        "uniform second-order heat remainder remains open" in old["status"],
        "old peeling-contract boundary drifted",
    )
    require(
        remainder["summary"]["eta_0"] == 100000
        and remainder["summary"]["eta_1"] == 200000,
        "first-order remainder constants drifted",
    )
    require(
        contains(
            remainder,
            "For every critical disk, including cutoff crossings, "
            "|r_[1](x)|<100000*exp(-5L/4)",
        ),
        "cutoff-uniform value theorem drifted",
    )
    require(
        contains(
            remainder,
            "|partial_x r_[1](x)|<200000*L*exp(-5L/4)",
        ),
        "cutoff-uniform derivative theorem drifted",
    )
    require(
        reconciliation["summary"]["superseded_handoffs"] == 1
        and reconciliation["summary"]["cutoff_uniform_first_order_remainders"]
        == 1,
        "first-order target reconciliation drifted",
    )
    require(
        contact["summary"]["eta_0"] == 100000
        and contact["summary"]["eta_1"] == 200000
        and contact["summary"]["open_frequency_targets"] == 1,
        "signed-contact boundary drifted",
    )
    require(
        contains(
            contact,
            "|X_[1]|<=50000*exp(-5L/4) implies "
            "|U_[1]|>100000*L*exp(-5L/4)",
        ),
        "box-optimal signed target drifted",
    )
    require(
        handoff["counts"]["required_c2_for_boundary_transfer"] == 0
        and handoff["counts"]["whole_jet_c1_homotopies"] == 1
        and handoff["counts"]["open_cross_current_or_abel_targets"] == 1,
        "C1 residual handoff drifted",
    )
    require(
        "signed carrier-plus-near estimate open" in pairing["status"]
        and pairing["counts"]["signed_complete_ideal_bounds"] == 0,
        "Section 11.193 arithmetic frontier drifted",
    )
    require(
        contains(pairing, "mathcal J_H^0=mathscr M_N[P_H^0]+mathscr N_N[P_H^0]+mathscr C_N[P_H^0]"),
        "Hermitian complete-join target drifted",
    )
    require(
        contains(
            ordinary,
            "T_L[O_h]>2*exp(-6h)/(625*x^(23/2))",
        ),
        "ordinary saddle-scale target drifted",
    )
    require(
        frontier["exact"]["direct_remaining_domain"] == "0<t<=1/5 and x>38",
        "outer-frequency domain drifted",
    )

    payloads = {
        "old_peeling_contract": old,
        "global_first_order_remainder": remainder,
        "target_reconciliation": reconciliation,
        "signed_contact": contact,
        "residual_handoff": handoff,
        "symmetric_outer_pairing": pairing,
        "ordinary_sqrt_transfer": ordinary,
        "outer_frontier": frontier,
    }
    source_audit = {}
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
            "id": "c1rsf_01_old_target_identified",
            "readiness": "ready_to_apply",
            "claim": (
                "The July 26 endpoint-peeling contract left the uniform "
                "second-order heat-integrated value/derivative remainder open."
            ),
        },
        {
            "id": "c1rsf_02_remainder_superseded",
            "readiness": "ready_to_apply",
            "claim": (
                "The later global theorem closes that obligation on L>=50 and "
                "0<tL<=25, including every fixed-cutoff cell and adjacent cutoff."
            ),
        },
        {
            "id": "c1rsf_03_global_c1_budget",
            "readiness": "ready_to_apply",
            "claim": (
                "Uniformly, |r_[1]|<100000 exp(-5L/4) and "
                "|r_[1],x|<200000 L exp(-5L/4)."
            ),
        },
        {
            "id": "c1rsf_04_contact_box",
            "readiness": "ready_to_apply",
            "claim": (
                "A contact lies in |X_[1]|<50000 exp(-5L/4), "
                "|U_[1]|<100000 L exp(-5L/4); strict avoidance of this "
                "rectangle is sufficient."
            ),
        },
        {
            "id": "c1rsf_05_c1_boundary_transfer",
            "readiness": "ready_to_apply",
            "claim": (
                "The whole-jet residual homotopy transfers a retained-main "
                "boundary theorem to Xi using C1 bounds only; no C2 residual "
                "bound is required for that transfer."
            ),
        },
        {
            "id": "c1rsf_06_live_c1_frontier",
            "readiness": "open",
            "claim": (
                "The live C1 arithmetic target is the source-specific signed "
                "carrier-plus-near package of Section 11.193, followed by its "
                "canonically paired remote correction and transpose terminal term."
            ),
        },
        {
            "id": "c1rsf_07_remote_pair_boundary",
            "readiness": "open",
            "claim": (
                "Remote positive/negative Fourier modes may be bounded only "
                "after symmetric pairing. A physical derivative-size estimate "
                "is still needed before the C2 pair bound enters the reserve."
            ),
        },
        {
            "id": "c1rsf_08_ordinary_alternative",
            "readiness": "open",
            "claim": (
                "The independent ordinary-theta route requires the direct "
                "saddle-scale first-jet lower bound for T_L[O_h]; its much "
                "smaller threshold cannot be inherited from the corrected main."
            ),
        },
        {
            "id": "c1rsf_09_outer_domain",
            "readiness": "ready_to_apply",
            "claim": "The direct unresolved Newman domain remains 0<t<=1/5 and x>38.",
        },
        {
            "id": "c1rsf_10_route_guard",
            "readiness": "ready_to_apply",
            "claim": (
                "Do not rederive the closed C1 remainder, split the conditionally "
                "convergent outer modes, transfer the ordinary tiny threshold "
                "through a coarser corrected-main envelope, or replace analytic "
                "termination by an indefinite finite tile ladder."
            ),
        },
    ]

    return {
        "kind": STEM,
        "schema_version": 1,
        "date": DATE,
        "status": (
            "cutoff-uniform C1 remainder supersession proved and live outer "
            "arithmetic frontier reconciled; the signed Xi theorem remains open"
        ),
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "source_audit": source_audit,
        "exact": {
            "superseded_target": (
                "uniform second-order heat-integrated value/derivative remainder "
                "after the explicit C1 endpoint peel"
            ),
            "closed_domain": "L>=50 and 0<tL<=25, including adjacent cutoffs",
            "global_remainder": {
                "eta_0": 100000,
                "eta_1": 200000,
                "value": "|r_[1]|<100000*exp(-5L/4)",
                "derivative": "|r_[1],x|<200000*L*exp(-5L/4)",
            },
            "contact_box": {
                "value": "|X_[1]|<50000*exp(-5L/4)",
                "derivative": "|U_[1]|<100000*L*exp(-5L/4)",
                "box_optimal_target": (
                    "|X_[1]|<=50000*exp(-5L/4) implies "
                    "|U_[1]|>100000*L*exp(-5L/4)"
                ),
            },
            "boundary_transfer": {
                "residual_coordinates": 2,
                "whole_jet_c1_homotopies": 1,
                "required_c2_residual_bounds": 0,
            },
            "live_c1_target": {
                "hermitian": (
                    "mathcal J_H^0=mathscr M_N[P_H^0]+mathscr N_N[P_H^0]"
                    "+mathscr C_N[P_H^0]"
                ),
                "transpose": (
                    "mathcal J_T^0=mathscr M_N[P_T^0]+mathscr N_N[P_T^0]"
                    "+mathscr C_N[P_T^0]+2iu_Nmathcal T_N^[N][B_(T,N)^0]"
                ),
                "open_obligation": (
                    "prove or rigorously obstruct the physical signed "
                    "carrier-plus-near estimate with endpoint half-weights, then "
                    "attach the symmetric remote pair and transpose terminal terms"
                ),
            },
            "ordinary_alternative": (
                "T_L[O_h]>2*exp(-6h)/(625*x^(23/2))"
            ),
            "direct_remaining_domain": "0<t<=1/5 and x>38",
        },
        "rows": rows,
        "summary": {
            "rows": len(rows),
            "superseded_remainder_targets": 1,
            "cutoff_uniform_c1_remainders": 1,
            "eta_0": 100000,
            "eta_1": 200000,
            "required_c2_boundary_bounds": 0,
            "ready_rows": sum(row["readiness"] == "ready_to_apply" for row in rows),
            "open_rows": sum(row["readiness"] == "open" for row in rows),
        },
        "next_action": (
            "Resume at Formal Core Section 11.193, not at the old peeling "
            "contract. Keep the endpoint half-weights and K_near intact. Derive "
            "the first nonvanishing normalized physical expansion of the joined "
            "Hermitian and transpose carrier-plus-near packages before taking "
            "absolute values; either prove a signed reserve or produce an "
            "admissible physical obstruction. Attach the symmetric remote C2 "
            "correction only after its physical derivative numerator is bounded."
        ),
        "proof_boundary": (
            "This gate proves that the old C1 remainder target is superseded and "
            "identifies the current arithmetic frontier. It proves no signed "
            "carrier-plus-near estimate, physical remote derivative bound, "
            "ordinary retained-main separation, Q209 shell, cofinal successor, "
            "Lambda<=0, RH, or prize-level conclusion."
        ),
    }


def success_line(payload: dict) -> str:
    summary = payload["summary"]
    return (
        "validated Newman C1 remainder supersession frontier gate: "
        f"{summary['rows']} rows, 0 issues, "
        f"eta_0={summary['eta_0']}, eta_1={summary['eta_1']}, "
        "0 C2 boundary requirements, 1 live carrier-plus-near target"
    )


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    lines = [
        "# Newman C1 Remainder Supersession and Live Arithmetic Frontier",
        "",
        f"Date: {DATE}",
        "",
        "Status: hygiene gate proving that the cutoff-uniform C1 remainder is closed and locating the live signed Xi arithmetic target; this is not a proof of RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Superseded Target",
        "",
        "The old endpoint-peeling contract correctly left a uniform second-order",
        "heat-integrated value/derivative theorem open. The later global theorem",
        "closes exactly that obligation, including adjacent cutoffs:",
        "",
        "```text",
        exact["global_remainder"]["value"],
        exact["global_remainder"]["derivative"],
        "L>=50, 0<tL<=25.",
        "```",
        "",
        "## Contact Box",
        "",
        "A hypothetical contact must satisfy",
        "",
        "```text",
        exact["contact_box"]["value"],
        exact["contact_box"]["derivative"],
        "```",
        "",
        "The whole-jet boundary homotopy uses only these C1 residual coordinates;",
        "it requires no C2 residual estimate.",
        "",
        "## Live C1 Target",
        "",
        "The unresolved arithmetic object is now the Section 11.193 package",
        "",
        "```text",
        exact["live_c1_target"]["hermitian"],
        exact["live_c1_target"]["transpose"],
        "```",
        "",
        "The finite carrier and near kernel must retain their endpoint half-weights.",
        "The positive and negative remote modes must remain symmetrically paired.",
        "",
        "## Independent Alternative",
        "",
        "The saddle-scale ordinary-theta route remains",
        "",
        "```text",
        exact["ordinary_alternative"],
        "```",
        "",
        "Its smaller threshold is not inherited from the corrected-main envelope.",
        "",
        "## Next Stage",
        "",
        payload["next_action"],
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


def main() -> int:
    payload = build_payload()
    write_text_atomic(RESULT, json.dumps(payload, indent=2) + "\n")
    write_text_atomic(NOTE, render_note(payload))
    print(success_line(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
