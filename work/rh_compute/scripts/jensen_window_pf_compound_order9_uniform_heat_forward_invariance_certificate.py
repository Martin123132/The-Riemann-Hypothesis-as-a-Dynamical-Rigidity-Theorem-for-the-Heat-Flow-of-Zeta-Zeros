#!/usr/bin/env python3
"""Propagate signed contiguous and arbitrary-column order nine to lambda zero."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate.md"
)
ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_m100_entry_certificate.json"
)
FLOW_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_uniform_tail_flow_reduction.json"
)
LOWER_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate.json"
)
TRANSFER_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_order4_noncontiguous_total_positivity_transfer.json"
)

ENDPOINT_ENTRY = (
    "Q_(9,n)(-100)=H_(9,n)(-100)>0 for every integer n>=0"
)


@dataclass(frozen=True)
class ForwardRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_record(path: Path, artifact: dict) -> dict:
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "sha256": sha256(path),
        "kind": artifact["kind"],
        "status": artifact["status"],
    }


def validate_sources() -> dict:
    entry = load_json(ENTRY_SOURCE)
    flow = load_json(FLOW_SOURCE)
    lower = load_json(LOWER_SOURCE)
    transfer = load_json(TRANSFER_SOURCE)
    if entry.get("exact", {}).get("all_shift_entry") != ENDPOINT_ENTRY:
        raise RuntimeError("order-nine endpoint entry source changed")
    if entry.get("summary", {}).get("all_shift_m100_entry_theorems") != 1:
        raise RuntimeError("order-nine endpoint entry source is not closed")
    if flow.get("summary", {}).get("universal_tail_specializations") != 1:
        raise RuntimeError("order-nine uniform tail source is not closed")
    if flow.get("summary", {}).get("cooperative_flow_recursions") != 1:
        raise RuntimeError("order-nine cooperative flow source is not closed")
    if flow.get("summary", {}).get("conditional_forward_theorems") != 1:
        raise RuntimeError("order-nine conditional forward theorem changed")
    if flow.get("summary", {}).get(
        "conditional_arbitrary_column_theorems"
    ) != 1:
        raise RuntimeError("order-nine conditional column theorem changed")
    if lower.get("summary", {}).get("contiguous_order_eight_theorems") != 1:
        raise RuntimeError("lower order-eight interval theorem changed")
    if lower.get("summary", {}).get(
        "arbitrary_column_order_eight_theorems"
    ) != 1:
        raise RuntimeError("lower arbitrary-column theorem changed")
    fixed_row = next(
        (
            row
            for row in transfer.get("rows", [])
            if row.get("id") == "o4ntp_09_fixed_order_transfer"
        ),
        None,
    )
    if fixed_row is None or fixed_row.get("readiness") != "ready_to_apply":
        raise RuntimeError("fixed-order arbitrary-column transfer is unavailable")

    exact = flow.get("exact", {})
    expected_flow = {
        "uniform_eventual_tail": (
            "there exists N_9 such that Q_(9,n)(lambda)>0 for every n>=N_9 "
            "and -100<=lambda<=0"
        ),
        "cooperative_flow": (
            "Q_n'=a_n*Q_(n+1)+b_n*Q_n, "
            "a_n=(4*n+66)*Q_(8,n)/Q_(8,n+1)>0, "
            "b_n=((4*n+66)/(4*n+62))*(log Q_(8,n+1))'"
        ),
        "conditional_forward": (
            "[Q_(9,n)(-100)>0 for all n] => [Q_(9,n)(lambda)>0 for all n "
            "and -100<=lambda<=0]"
        ),
        "conditional_arbitrary_columns": (
            "completed signed contiguous layers through order nine imply every "
            "arbitrary-column signed layer through order nine"
        ),
    }
    for key, expected in expected_flow.items():
        if exact.get(key) != expected:
            raise RuntimeError(f"order-nine {key} source changed")
    return {
        "sources": [
            source_record(ENTRY_SOURCE, entry),
            source_record(FLOW_SOURCE, flow),
            source_record(LOWER_SOURCE, lower),
            source_record(TRANSFER_SOURCE, transfer),
        ],
        "endpoint_entry": ENDPOINT_ENTRY,
        "uniform_tail": exact["uniform_eventual_tail"],
        "cooperative_flow": exact["cooperative_flow"],
        "variation_of_constants": exact["variation_of_constants"],
        "conditional_forward": exact["conditional_forward"],
        "conditional_arbitrary_columns": exact[
            "conditional_arbitrary_columns"
        ],
        "lower_layers": (
            "epsilon_k*H_(k,n)(lambda)>0 for k=1,2,3,4,5,6,7,8, all n "
            "and -100<=lambda<=0"
        ),
        "fixed_order_transfer": fixed_row["formula"],
    }


def conclusions() -> dict:
    return {
        "contiguous_order_nine": (
            "Q_(9,n)(lambda)=H_(9,n)(lambda)>0 for every n>=0 and "
            "-100<=lambda<=0"
        ),
        "signed_layers_through_nine": (
            "epsilon_k*H_(k,n)(lambda)>0 for 1<=k<=9, every n>=0 and "
            "-100<=lambda<=0"
        ),
        "arbitrary_order_nine": (
            "epsilon_9*R_(9,n)(j_1,j_2,j_3,j_4,j_5,j_6,j_7,j_8,j_9;lambda)>0 "
            "for every n>=0, 0<=j_1<j_2<j_3<j_4<j_5<j_6<j_7<j_8<j_9, "
            "and -100<=lambda<=0"
        ),
        "fixed_order_nine": (
            "every consecutive-row arbitrary-column reshaped-Hankel minor "
            "through order nine has its required sign on -100<=lambda<=0"
        ),
        "next_handoff": (
            "derive a uniform-in-order endpoint and forward mechanism; fixed "
            "order nine alone is not PF-infinity"
        ),
    }


def build_artifact() -> dict:
    sources = validate_sources()
    theorem = conclusions()
    rows = [
        ForwardRow(
            "co9uhfic_01_endpoint",
            "theorem_input",
            "ready_to_apply",
            "Every signed order-nine endpoint coordinate is strictly positive at lambda=-100.",
            sources["endpoint_entry"],
            "Completed zeta-specific endpoint theorem.",
        ),
        ForwardRow(
            "co9uhfic_02_uniform_tail",
            "theorem_input",
            "ready_to_apply",
            "Signed order-nine positivity is uniform at sufficiently large shift throughout the compact heat interval.",
            sources["uniform_tail"],
            "Non-effective but compact-uniform eventual theorem.",
        ),
        ForwardRow(
            "co9uhfic_03_cooperative_flow",
            "exact_identity",
            "ready_to_apply",
            "The signed order-nine layer obeys a nearest-neighbor cooperative system with positive forward coefficient.",
            sources["cooperative_flow"],
            "Exact Plucker and heat-flow algebra inside the completed order-eight cone.",
        ),
        ForwardRow(
            "co9uhfic_04_finite_confinement",
            "exact_maximum_principle",
            "ready_to_apply",
            "The uniform tail confines any first loss to finitely many shifts, where variation of constants preserves positivity.",
            sources["variation_of_constants"],
            "Infinite-index maximum principle with compact-uniform tail.",
        ),
        ForwardRow(
            "co9uhfic_05_contiguous_theorem",
            "exact_forward_theorem",
            "ready_to_apply",
            "Endpoint entry, the uniform tail, and the cooperative system prove signed contiguous order nine throughout the interval.",
            theorem["contiguous_order_nine"],
            "All shifts and all finite heat parameters from -100 through zero.",
        ),
        ForwardRow(
            "co9uhfic_06_lower_layers",
            "theorem_input",
            "ready_to_apply",
            "The completed lower contiguous layers supply every initial-minor sign through order eight.",
            sources["lower_layers"],
            "Previously proved lower-layer interval theorems.",
        ),
        ForwardRow(
            "co9uhfic_07_arbitrary_columns",
            "published_theorem_composition",
            "ready_to_apply",
            "The fixed-order initial-minor criterion transfers signed contiguous order nine to all arbitrary column sets.",
            theorem["arbitrary_order_nine"],
            "Finite-dimensional transfer applied pointwise in lambda.",
        ),
        ForwardRow(
            "co9uhfic_08_all_order_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "Fixed order nine does not supply the all-order theorem required by the Jensen bridge.",
            theorem["next_handoff"],
            "Open uniform-in-order problem; not PF-infinity, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate",
        "date": "2026-07-14",
        "status": (
            "signed contiguous and arbitrary-column order-nine theorem on "
            "-100<=lambda<=0"
        ),
        "proof_boundary": (
            "This artifact proves the required signed Hankel layers through fixed "
            "order nine only. It does not prove order ten or higher, PF-infinity, "
            "the all-degree Jensen bridge, RH, or Lambda<=0."
        ),
        "source_contract": sources,
        "conclusions": theorem,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": sum(row.readiness == "ready_to_apply" for row in rows),
            "open_rows": sum(row.readiness == "not_ready_to_apply" for row in rows),
            "uniform_tail_theorems": 1,
            "cooperative_flow_theorems": 1,
            "contiguous_order_nine_theorems": 1,
            "arbitrary_column_order_nine_theorems": 1,
            "open_all_order_handoffs": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    source = artifact["source_contract"]
    theorem = artifact["conclusions"]
    lines = [
        "# Order-nine uniform-heat forward-invariance certificate",
        "",
        "Date: 2026-07-14",
        "",
        "Status: signed contiguous and arbitrary-column order-nine theorem on `-100<=lambda<=0`; This is not a proof of the all-degree Jensen bridge, PF-infinity, RH,",
        "or `Lambda <= 0`; it proves one fixed compound layer.",
        "",
        "## Cooperative propagation",
        "",
        f"Endpoint entry: `{source['endpoint_entry']}`.",
        "",
        f"Uniform tail: `{source['uniform_tail']}`.",
        "",
        f"Cooperative flow: `{source['cooperative_flow']}`.",
        "",
        "Variation of constants and finite confinement then prove",
        "",
        f"`{theorem['contiguous_order_nine']}`.",
        "",
        "## Arbitrary columns",
        "",
        "The completed signed contiguous layers through order nine satisfy the",
        "fixed-order initial-minor hypotheses pointwise in `lambda`, giving",
        "",
        f"`{theorem['arbitrary_order_nine']}`.",
        "",
        "## Boundary",
        "",
        f"`{theorem['next_handoff']}`.",
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
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(
        "wrote order-nine uniform-heat forward invariance: "
        "one contiguous and one arbitrary-column theorem"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
