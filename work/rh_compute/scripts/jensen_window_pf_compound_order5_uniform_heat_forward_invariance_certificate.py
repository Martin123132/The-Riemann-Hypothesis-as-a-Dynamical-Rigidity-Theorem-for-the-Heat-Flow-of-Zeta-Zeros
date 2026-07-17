#!/usr/bin/env python3
"""Propagate contiguous and arbitrary-column order five through lambda zero."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.md"
ENTRY_SOURCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order5_m100_entry_certificate.json"
FLOW_SOURCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.json"
LOWER_SOURCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.json"
TRANSFER_SOURCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.json"


@dataclass(frozen=True)
class ForwardRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_sources() -> dict:
    entry = load_json(ENTRY_SOURCE)
    flow = load_json(FLOW_SOURCE)
    lower = load_json(LOWER_SOURCE)
    transfer = load_json(TRANSFER_SOURCE)
    if entry.get("exact", {}).get("all_shift_entry") != "H_(5,n)(-100)>0 for every integer n>=0":
        raise RuntimeError("order-five endpoint entry source changed")
    if flow.get("summary", {}).get("uniform_eventual_tail_theorems") != 1:
        raise RuntimeError("order-five uniform tail source is not closed")
    if flow.get("summary", {}).get("cooperative_flow_theorems") != 1:
        raise RuntimeError("order-five cooperative flow source is not closed")
    if lower.get("summary", {}).get("lambda_zero_all_shift_theorems") != 1:
        raise RuntimeError("lower order-four interval theorem changed")
    fixed_row = next(
        (row for row in transfer.get("rows", []) if row.get("id") == "o4ntp_09_fixed_order_transfer"),
        None,
    )
    if fixed_row is None or fixed_row.get("readiness") != "ready_to_apply":
        raise RuntimeError("fixed-order arbitrary-column transfer is unavailable")
    return {
        "endpoint_entry": entry["exact"]["all_shift_entry"],
        "uniform_tail": flow["conclusions"]["uniform_eventual_tail"],
        "cooperative_flow": flow["exact_flow"]["cooperative_flow"],
        "variation_of_constants": flow["conclusions"]["variation_of_constants"],
        "conditional_forward": flow["conclusions"]["conditional_forward"],
        "lower_layers": "epsilon_k*H_(k,n)(lambda)>0 for k=1,2,3,4, all n and -100<=lambda<=0",
        "fixed_order_transfer": fixed_row["formula"],
    }


def conclusions() -> dict:
    return {
        "contiguous_order_five": "H_(5,n)(lambda)>0 for every n>=0 and -100<=lambda<=0",
        "signed_layers_through_five": "epsilon_k*H_(k,n)(lambda)>0 for 1<=k<=5, every n>=0 and -100<=lambda<=0",
        "arbitrary_order_five": "R_(5,n)(j_1,j_2,j_3,j_4,j_5;lambda)>0 for every n>=0, 0<=j_1<j_2<j_3<j_4<j_5, and -100<=lambda<=0",
        "fixed_order_five": "every consecutive-row arbitrary-column reshaped-Hankel minor through order five has its required sign on -100<=lambda<=0",
        "next_handoff": "derive a cancellation-preserving contiguous order-six entry coordinate and flow reduction",
    }


def build_artifact() -> dict:
    sources = validate_sources()
    theorem = conclusions()
    rows = [
        ForwardRow("co5uhfic_01_endpoint", "theorem_input", "ready_to_apply", "Every order-five endpoint coordinate is strictly positive at lambda=-100.", sources["endpoint_entry"], "Completed zeta-specific endpoint theorem."),
        ForwardRow("co5uhfic_02_uniform_tail", "theorem_input", "ready_to_apply", "Order-five positivity is uniform at sufficiently large shift throughout the compact heat interval.", sources["uniform_tail"], "Non-effective but compact-uniform eventual theorem."),
        ForwardRow("co5uhfic_03_cooperative_flow", "exact_identity", "ready_to_apply", "The normalized order-five layer obeys a nearest-neighbor cooperative system with positive forward coefficient.", sources["cooperative_flow"], "Exact Plucker and heat-flow algebra."),
        ForwardRow("co5uhfic_04_finite_confinement", "exact_maximum_principle", "ready_to_apply", "The uniform tail confines any first loss to finitely many shifts, where variation of constants preserves positivity.", sources["variation_of_constants"], "Infinite-index maximum principle with compact-uniform tail."),
        ForwardRow("co5uhfic_05_contiguous_theorem", "exact_forward_theorem", "ready_to_apply", "The endpoint, uniform tail, and cooperative system prove contiguous order five throughout the interval.", theorem["contiguous_order_five"], "All shifts and all finite heat parameters from -100 through zero."),
        ForwardRow("co5uhfic_06_lower_layers", "theorem_input", "ready_to_apply", "The completed lower contiguous layers supply every initial-minor sign through order four.", sources["lower_layers"], "Previously proved lower-layer interval theorems."),
        ForwardRow("co5uhfic_07_arbitrary_columns", "published_theorem_composition", "ready_to_apply", "The fixed-order initial-minor criterion transfers signed contiguous order five to all arbitrary column sets.", theorem["arbitrary_order_five"], "Finite-dimensional transfer applied pointwise in lambda."),
        ForwardRow("co5uhfic_08_order_six_handoff", "open_theorem_target", "not_ready_to_apply", "The fixed-order programme now hands off to a genuinely new contiguous order-six layer.", theorem["next_handoff"], "Open order-six analytic problem; not PF-infinity, RH, or Lambda<=0."),
    ]
    return {
        "kind": "jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate",
        "date": "2026-07-13",
        "status": "contiguous and arbitrary-column order-five theorem on -100<=lambda<=0",
        "proof_boundary": "This artifact proves the required signed Hankel layer through fixed order five only. It does not prove order six or higher, PF-infinity, the all-degree Jensen bridge, RH, or Lambda<=0.",
        "sources": [
            "outputs/jensen_window_pf_compound_order5_m100_entry_certificate.md",
            "outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md",
            "outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md",
            "outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md",
            "outputs/formal_core.md",
        ],
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.py",
        "source_contract": sources,
        "conclusions": theorem,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": sum(row.readiness == "ready_to_apply" for row in rows),
            "open_rows": sum(row.readiness == "not_ready_to_apply" for row in rows),
            "uniform_tail_theorems": 1,
            "cooperative_flow_theorems": 1,
            "contiguous_order_five_theorems": 1,
            "arbitrary_column_order_five_theorems": 1,
            "open_order_six_handoffs": 1,
        },
    }


def write_note(path: Path, artifact: dict) -> None:
    source = artifact["source_contract"]
    theorem = artifact["conclusions"]
    lines = [
        "# Jensen-Window PF Compound Order-Five Uniform-Heat Forward Invariance Certificate",
        "", "Date: 2026-07-13", "",
        "Status: contiguous and arbitrary-column order-five theorem on `-100<=lambda<=0`.",
        "This is not a proof of the all-degree Jensen bridge, RH, or `Lambda <= 0`; it proves one fixed compound layer, not PF-infinity.",
        "", "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.py",
        "```", "", "## Cooperative Propagation", "",
        "The endpoint theorem now supplies", "", "```text", source["endpoint_entry"], "```", "",
        "The compact-uniform eventual theorem and exact normalized flow are", "", "```text",
        source["uniform_tail"], source["cooperative_flow"], "```", "",
        "The uniform tail confines a hypothetical first loss to finitely many shifts. Variation of constants on that finite set gives", "", "```text",
        source["variation_of_constants"], "```", "",
        "whose integrand is nonnegative under backward induction from the positive tail. Therefore", "", "```text",
        theorem["contiguous_order_five"], "```", "", "## Arbitrary Columns", "",
        "The signed contiguous layers through order five satisfy the initial-minor hypotheses pointwise in `lambda`. The fixed-order Gasca-Pena transfer therefore gives", "", "```text",
        theorem["arbitrary_order_five"], "```", "",
        "Thus contiguous and arbitrary-column order five are both closed on the full heat interval. The next fixed-order task is", "", "```text",
        theorem["next_handoff"], "```", "",
        "```text",
        "outputs/jensen_window_pf_compound_order5_m100_entry_certificate.md",
        "outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md",
        "outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md",
        "outputs/signed_hankel_jensen_bridge_target.md", "outputs/formal_core.md", "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print("wrote order-five uniform-heat forward invariance: "
          f"{summary['rows']} rows, {summary['contiguous_order_five_theorems']} contiguous theorem, "
          f"{summary['arbitrary_column_order_five_theorems']} arbitrary-column theorem, "
          f"{summary['open_order_six_handoffs']} open order-six handoff")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
