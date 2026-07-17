#!/usr/bin/env python3
"""Compose the delayed heat ray and finite prefix at lambda zero."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS = REPO_ROOT / "work/rh_compute/results"
DELAYED_ENTRY_SOURCE = RESULTS / "jensen_window_pf_compound_order10_m100_delayed_entry_certificate.json"
HEAT_LEMMA_SOURCE = RESULTS / "jensen_window_pf_delayed_cooperative_heat_tail_lemma.json"
PREFIX_SOURCE = RESULTS / "jensen_window_pf_compound_order10_lambda0_prefix_certificate.json"
ALL_ORDER_HEAT_SOURCE = RESULTS / "jensen_window_pf_all_order_endpoint_heat_reduction.json"
TRANSFER_SOURCE = RESULTS / "jensen_window_pf_order4_noncontiguous_total_positivity_transfer.json"
DEFAULT_OUT = RESULTS / "jensen_window_pf_compound_order10_lambda0_completion_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order10_lambda0_completion_certificate.md"

DELAYED_ENDPOINT = "Q_(10,n)(-100)>0 for every integer n>=4"
HEAT_RAY = "Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"
PREFIX = "Q_(10,n)(0)>0 for every integer 0<=n<=3"
ALL_SHIFT_Q10 = "Q_(10,n)(0)>0 for every integer n>=0"
LOWER_BASE = "Q_(m,n)(lambda)>0 for 1<=m<=9, every n>=0, and -100<=lambda<=0"
CONTIGUOUS_THROUGH_10 = "Q_(m,n)(0)>0 for every integer 1<=m<=10 and n>=0"
TRANSFER_THEOREM = "[epsilon_k H_(k,s)>0 for 1<=k<=m, all s] => [epsilon_k R_(k,n)(j_1,...,j_k)>0 for 1<=k<=m]"
ARBITRARY_THROUGH_10 = "epsilon_k*R_(k,n)(j_1,...,j_k)(0)>0 for 1<=k<=10, n>=0, and 0<=j_1<...<j_k"


@dataclass(frozen=True)
class CompletionRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_record(path: Path, payload: dict) -> dict:
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "sha256": sha256(path),
        "kind": payload.get("kind"),
        "status": payload.get("status"),
    }


def validate_sources() -> list[dict]:
    delayed = load_json(DELAYED_ENTRY_SOURCE)
    lemma = load_json(HEAT_LEMMA_SOURCE)
    prefix = load_json(PREFIX_SOURCE)
    heat = load_json(ALL_ORDER_HEAT_SOURCE)
    transfer = load_json(TRANSFER_SOURCE)

    if delayed.get("exact", {}).get("delayed_entry") != DELAYED_ENDPOINT:
        raise RuntimeError("delayed endpoint theorem changed")
    if delayed.get("summary", {}).get("delayed_m100_entry_theorems") != 1:
        raise RuntimeError("delayed endpoint theorem is not closed")
    if lemma.get("exact", {}).get("order10_delayed_handoff") != (
        f"[{DELAYED_ENDPOINT}] implies {HEAT_RAY}"
    ):
        raise RuntimeError("delayed heat handoff changed")
    if lemma.get("summary", {}).get("order10_n4_specializations") != 1:
        raise RuntimeError("delayed heat lemma is not closed")
    if prefix.get("finite", {}).get("theorem") != PREFIX:
        raise RuntimeError("lambda-zero prefix theorem changed")
    if prefix.get("summary", {}).get("positive_Q10_rows") != 4:
        raise RuntimeError("lambda-zero prefix certificate is not closed")
    if heat.get("exact", {}).get("known_base") != LOWER_BASE:
        raise RuntimeError("lower heat base changed")
    if heat.get("summary", {}).get("completed_base_order") != 9:
        raise RuntimeError("lower heat base is not closed through order nine")

    transfer_rows = {row.get("id"): row for row in transfer.get("rows", [])}
    fixed_transfer = transfer_rows.get("o4ntp_09_fixed_order_transfer", {})
    if fixed_transfer.get("formula") != TRANSFER_THEOREM:
        raise RuntimeError("fixed-order arbitrary-column transfer changed")
    if transfer.get("summary", {}).get("fixed_order_transfer_theorems") != 1:
        raise RuntimeError("fixed-order arbitrary-column transfer is not closed")

    payloads = (delayed, lemma, prefix, heat, transfer)
    paths = (
        DELAYED_ENTRY_SOURCE,
        HEAT_LEMMA_SOURCE,
        PREFIX_SOURCE,
        ALL_ORDER_HEAT_SOURCE,
        TRANSFER_SOURCE,
    )
    return [source_record(path, payload) for path, payload in zip(paths, payloads)]


def build_artifact() -> dict:
    sources = validate_sources()
    rows = [
        CompletionRow(
            "co10l0c_01_delayed_endpoint",
            "theorem_input",
            "ready_to_apply",
            "The endpoint analysis supplies strict order-ten positivity on the ray beginning at shift four.",
            DELAYED_ENDPOINT,
            "The four smaller endpoint shifts remain negative.",
        ),
        CompletionRow(
            "co10l0c_02_heat_ray",
            "theorem_application",
            "ready_to_apply",
            "The shifted cooperative lemma propagates the delayed endpoint ray through the heat interval.",
            HEAT_RAY,
            "No sign at shifts below four is used.",
        ),
        CompletionRow(
            "co10l0c_03_lambda0_prefix",
            "finite_theorem_input",
            "ready_to_apply",
            "Direct Arb determinants supply the four missing rows at lambda zero.",
            PREFIX,
            "Fixed lambda zero and shifts zero through three only.",
        ),
        CompletionRow(
            "co10l0c_04_all_shift_order10",
            "exact_theorem_composition",
            "ready_to_apply",
            "The heat ray and finite prefix partition every nonnegative shift.",
            ALL_SHIFT_Q10,
            "Contiguous signed Hankel order ten at lambda zero.",
        ),
        CompletionRow(
            "co10l0c_05_contiguous_through10",
            "exact_theorem_composition",
            "ready_to_apply",
            "The completed lower base and order-ten theorem give all contiguous layers through ten.",
            CONTIGUOUS_THROUGH_10,
            "Fixed lambda zero; no order above ten is claimed.",
        ),
        CompletionRow(
            "co10l0c_06_fixed_order_transfer",
            "published_theorem_application",
            "ready_to_apply",
            "The finite initial-minor transfer promotes completed contiguous layers through ten to arbitrary columns.",
            ARBITRARY_THROUGH_10,
            "Consecutive rows and arbitrary increasing columns, orders at most ten.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order10_lambda0_completion_certificate",
        "date": "2026-07-16",
        "status": "rigorous all-shift signed Hankel order ten and fixed-order arbitrary-column completion at lambda zero",
        "proof_boundary": (
            "This proves contiguous signed Hankel positivity at lambda zero through "
            "order ten and the corresponding consecutive-row arbitrary-column signs "
            "through order ten. It proves no order above ten and is not PF-infinity, "
            "RH, or Lambda<=0."
        ),
        "sources": sources,
        "exact": {
            "delayed_endpoint": DELAYED_ENDPOINT,
            "delayed_heat_ray": HEAT_RAY,
            "lambda0_prefix": PREFIX,
            "all_shift_order10_lambda0": ALL_SHIFT_Q10,
            "lower_base": LOWER_BASE,
            "contiguous_through_order10_lambda0": CONTIGUOUS_THROUGH_10,
            "arbitrary_column_definition": "R_(k,n)(j_1,...,j_k)(0)=det[A_(n+i+j_l)(0)]_(0<=i<k,1<=l<=k)",
            "fixed_order_transfer": TRANSFER_THEOREM,
            "arbitrary_columns_through_order10_lambda0": ARBITRARY_THROUGH_10,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "all_shift_order10_lambda0_theorems": 1,
            "contiguous_through_order10_lambda0_theorems": 1,
            "arbitrary_column_through_order10_lambda0_theorems": 1,
            "orders_above_10": 0,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order10_lambda0_completion_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order10_lambda0_completion_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Order-Ten Lambda-Zero Completion Certificate",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous all-shift signed Hankel order-ten completion at lambda zero. This certificate is not a proof of PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "## Contiguous Layer",
        "",
        "The delayed endpoint theorem and shifted cooperative heat lemma give",
        "",
        "```text",
        exact["delayed_endpoint"],
        exact["delayed_heat_ray"],
        "```",
        "",
        "At lambda zero the direct Arb prefix supplies shifts `0,1,2,3`.",
        "Together these disjoint pieces prove",
        "",
        "```text",
        exact["all_shift_order10_lambda0"],
        exact["contiguous_through_order10_lambda0"],
        "```",
        "",
        "## Arbitrary Columns",
        "",
        "The fixed-order initial-minor theorem now applies with `m=10`:",
        "",
        "```text",
        exact["arbitrary_column_definition"],
        exact["arbitrary_columns_through_order10_lambda0"],
        "```",
        "",
        "This is a genuine new fixed-order theorem, but the boundary matters:",
        "nothing here proves order eleven or any higher order. It is not",
        "PF-infinity, RH, or `Lambda<=0`.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print("wrote all-shift order-ten lambda-zero completion and arbitrary-column transfer")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
