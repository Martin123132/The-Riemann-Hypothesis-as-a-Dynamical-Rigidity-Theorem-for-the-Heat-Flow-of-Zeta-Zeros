#!/usr/bin/env python3
"""Compose fixed order-eleven signed-Hankel completion at lambda zero."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS = REPO_ROOT / "work/rh_compute/results"
ENDPOINT_SOURCE = RESULTS / "jensen_window_pf_compound_order11_m100_entry_certificate.json"
DELAYED_HEAT_SOURCE = RESULTS / "jensen_window_pf_delayed_cooperative_heat_tail_lemma.json"
ORDER10_COMPLETION = RESULTS / "jensen_window_pf_compound_order10_lambda0_completion_certificate.json"
PREFIX_SOURCE = RESULTS / "jensen_window_pf_compound_order11_lambda0_prefix_certificate.json"
TRANSFER_SOURCE = RESULTS / "jensen_window_pf_order4_noncontiguous_total_positivity_transfer.json"
DEFAULT_OUT = RESULTS / "jensen_window_pf_compound_order11_lambda0_completion_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order11_lambda0_completion_certificate.md"
GENERATOR_PATH = (
    "work/rh_compute/scripts/"
    "jensen_window_pf_compound_order11_lambda0_completion_certificate.py"
)
CHECKER_PATH = (
    "work/rh_compute/scripts/"
    "check_jensen_window_pf_compound_order11_lambda0_completion_certificate.py"
)

ENDPOINT = "Q_(11,n)(-100)>0 for every integer n>=0"
ORDER10_HEAT = "Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"
DELAYED_HANDOFF = (
    "[Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0 and "
    "Q_(11,n)(-100)>0 for every n>=4] implies "
    "Q_(11,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"
)
ORDER11_HEAT = "Q_(11,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"
PREFIX = "Q_(11,n)(0)>0 for every integer 0<=n<=3"
ALL_SHIFT_Q11 = "Q_(11,n)(0)>0 for every integer n>=0"
CONTIGUOUS_THROUGH_10 = (
    "Q_(m,n)(0)>0 for every integer 1<=m<=10 and n>=0"
)
CONTIGUOUS_THROUGH_11 = (
    "Q_(m,n)(0)>0 for every integer 1<=m<=11 and n>=0"
)
TRANSFER_THEOREM = (
    "[epsilon_k H_(k,s)>0 for 1<=k<=m, all s] => "
    "[epsilon_k R_(k,n)(j_1,...,j_k)>0 for 1<=k<=m]"
)
ARBITRARY_THROUGH_11 = (
    "epsilon_k*R_(k,n)(j_1,...,j_k)(0)>0 for 1<=k<=11, n>=0, "
    "and 0<=j_1<...<j_k"
)
NONPROMOTION_GUARD = (
    "contiguous order four alone does not imply arbitrary-column order four; "
    "the lower initial-minor signs are essential"
)


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


def source_record(path: Path, artifact: dict) -> dict:
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "sha256": sha256(path),
        "kind": artifact["kind"],
        "status": artifact["status"],
    }


def validate_sources() -> list[dict]:
    endpoint = load_json(ENDPOINT_SOURCE)
    delayed = load_json(DELAYED_HEAT_SOURCE)
    order10 = load_json(ORDER10_COMPLETION)
    prefix = load_json(PREFIX_SOURCE)
    transfer = load_json(TRANSFER_SOURCE)
    if (
        endpoint.get("kind")
        != "jensen_window_pf_compound_order11_m100_entry_certificate"
        or endpoint.get("status")
        != "rigorous all-shift signed order-eleven entry at lambda=-100"
        or endpoint.get("exact", {}).get("global_endpoint") != ENDPOINT
        or endpoint.get("summary", {}).get(
            "global_m100_order11_entry_theorems"
        )
        != 1
        or endpoint.get("summary", {}).get("open_rows") != 0
        or endpoint.get("summary", {}).get("heat_interval_theorems") != 0
        or endpoint.get("summary", {}).get("orders_above_11") != 0
        or endpoint.get("summary", {}).get("rh_claims") != 0
    ):
        raise RuntimeError("order-eleven endpoint source is not closed")
    if (
        delayed.get("kind") != "jensen_window_pf_delayed_cooperative_heat_tail_lemma"
        or delayed.get("status") != "exact shifted cooperative heat-flow descent lemma"
        or delayed.get("exact", {}).get("order11_delayed_handoff")
        != DELAYED_HANDOFF
        or delayed.get("summary", {}).get("order11_n4_specializations") != 1
        or delayed.get("summary", {}).get("open_rows") != 0
    ):
        raise RuntimeError("order-eleven delayed heat source changed")
    exact10 = order10.get("exact", {})
    if (
        order10.get("kind")
        != "jensen_window_pf_compound_order10_lambda0_completion_certificate"
        or order10.get("status")
        != (
            "rigorous all-shift signed Hankel order ten and fixed-order "
            "arbitrary-column completion at lambda zero"
        )
        or exact10.get("delayed_heat_ray") != ORDER10_HEAT
        or exact10.get("contiguous_through_order10_lambda0")
        != CONTIGUOUS_THROUGH_10
        or order10.get("summary", {}).get(
            "contiguous_through_order10_lambda0_theorems"
        )
        != 1
        or order10.get("summary", {}).get("open_rows") != 0
        or order10.get("summary", {}).get("orders_above_10") != 0
    ):
        raise RuntimeError("order-ten lambda-zero source is not closed")
    if (
        prefix.get("kind")
        != "jensen_window_pf_compound_order11_lambda0_prefix_certificate"
        or prefix.get("status")
        != "rigorous lambda-zero order-eleven prefix certificate for n=0,1,2,3"
        or prefix.get("finite", {}).get("theorem") != PREFIX
        or prefix.get("finite", {}).get("all_Q11_positive") is not True
        or prefix.get("summary", {}).get("positive_Q11_rows") != 4
        or prefix.get("summary", {}).get("inconclusive_rows") != 0
    ):
        raise RuntimeError("order-eleven lambda-zero prefix is not closed")
    transfer_rows = {row.get("id"): row for row in transfer.get("rows", [])}
    fixed_transfer = transfer_rows.get("o4ntp_09_fixed_order_transfer", {})
    countermodel = transfer.get("exact", {}).get("nonpromotion_countermodel", {})
    if (
        fixed_transfer.get("readiness") != "ready_to_apply"
        or fixed_transfer.get("formula") != TRANSFER_THEOREM
        or countermodel.get("conclusion") != NONPROMOTION_GUARD
        or transfer.get("summary", {}).get("fixed_order_transfer_theorems") != 1
        or transfer.get("summary", {}).get("nonpromotion_guards") != 1
    ):
        raise RuntimeError("fixed-order arbitrary-column transfer changed")
    return [
        source_record(path, artifact)
        for path, artifact in (
            (ENDPOINT_SOURCE, endpoint),
            (DELAYED_HEAT_SOURCE, delayed),
            (ORDER10_COMPLETION, order10),
            (PREFIX_SOURCE, prefix),
            (TRANSFER_SOURCE, transfer),
        )
    ]


def build_artifact() -> dict:
    sources = validate_sources()
    rows = [
        CompletionRow(
            "co11l0c_01_endpoint",
            "theorem_input",
            "ready_to_apply",
            "The completed curvature bridge gives every order-eleven endpoint shift at lambda=-100.",
            ENDPOINT,
            "Fixed order eleven and lambda=-100.",
        ),
        CompletionRow(
            "co11l0c_02_order10_heat_base",
            "theorem_input",
            "ready_to_apply",
            "The completed order-ten layer supplies the positive coefficient in shifted cooperative descent.",
            ORDER10_HEAT,
            "Shifts n>=4 on the heat interval.",
        ),
        CompletionRow(
            "co11l0c_03_delayed_heat",
            "theorem_application",
            "ready_to_apply",
            "The delayed cooperative lemma propagates the order-eleven endpoint ray through the heat interval.",
            ORDER11_HEAT,
            "Uses no order-eleven shift below four.",
        ),
        CompletionRow(
            "co11l0c_04_lambda0_prefix",
            "finite_theorem_input",
            "ready_to_apply",
            "Direct Arb determinants supply the four complementary lambda-zero shifts.",
            PREFIX,
            "Lambda zero and shifts zero through three only.",
        ),
        CompletionRow(
            "co11l0c_05_all_shift_order11",
            "exact_theorem_composition",
            "ready_to_apply",
            "The delayed heat ray and direct prefix partition every nonnegative shift at lambda zero.",
            ALL_SHIFT_Q11,
            "Contiguous signed Hankel order eleven at lambda zero.",
        ),
        CompletionRow(
            "co11l0c_06_contiguous_through11",
            "exact_theorem_composition",
            "ready_to_apply",
            "The completed lower layers and order-eleven theorem give every contiguous layer through eleven.",
            CONTIGUOUS_THROUGH_11,
            "Fixed lambda zero; no order above eleven is claimed.",
        ),
        CompletionRow(
            "co11l0c_07_nonpromotion_guard",
            "countermodel_guard",
            "ready_to_apply",
            "An explicit countermodel forbids promotion from one contiguous layer without all lower initial-minor signs.",
            NONPROMOTION_GUARD,
            "This composition supplies layers one through eleven and makes no inference to order twelve.",
        ),
        CompletionRow(
            "co11l0c_08_fixed_order_transfer",
            "published_theorem_application",
            "ready_to_apply",
            "The initial-minor transfer promotes the completed contiguous layers to arbitrary increasing columns.",
            ARBITRARY_THROUGH_11,
            "Consecutive rows, arbitrary columns, and orders at most eleven.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order11_lambda0_completion_certificate",
        "date": "2026-07-22",
        "status": (
            "rigorous all-shift signed Hankel order eleven and fixed-order "
            "arbitrary-column completion at lambda zero"
        ),
        "proof_boundary": (
            "This proves contiguous signed Hankel positivity at lambda zero "
            "through order eleven and the corresponding consecutive-row "
            "arbitrary-column signs through order eleven. It proves no order "
            "above eleven and is not PF-infinity, Lambda<=0, or RH."
        ),
        "sources": sources,
        "exact": {
            "endpoint": ENDPOINT,
            "order10_heat_base": ORDER10_HEAT,
            "delayed_handoff": DELAYED_HANDOFF,
            "delayed_order11_heat_ray": ORDER11_HEAT,
            "lambda0_prefix": PREFIX,
            "all_shift_order11_lambda0": ALL_SHIFT_Q11,
            "contiguous_through_order10_lambda0": CONTIGUOUS_THROUGH_10,
            "contiguous_through_order11_lambda0": CONTIGUOUS_THROUGH_11,
            "fixed_order_transfer": TRANSFER_THEOREM,
            "nonpromotion_guard": NONPROMOTION_GUARD,
            "arbitrary_columns_through_order11_lambda0": ARBITRARY_THROUGH_11,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "delayed_order11_heat_theorems": 1,
            "all_shift_order11_lambda0_theorems": 1,
            "contiguous_through_order11_lambda0_theorems": 1,
            "arbitrary_column_through_order11_lambda0_theorems": 1,
            "countermodel_nonpromotion_guards": 1,
            "orders_above_11": 0,
            "pf_infinity_theorems": 0,
            "rh_claims": 0,
        },
        "generator": GENERATOR_PATH,
        "checker": CHECKER_PATH,
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Order-Eleven Lambda-Zero Completion Certificate",
        "",
        "Date: 2026-07-22",
        "",
        f"Status: **{artifact['status']}**. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "## Contiguous Layer",
        "",
        f"`{exact['endpoint']}`",
        f"`{exact['delayed_order11_heat_ray']}`",
        f"`{exact['lambda0_prefix']}`",
        "",
        f"Together these prove `{exact['all_shift_order11_lambda0']}` and",
        f"`{exact['contiguous_through_order11_lambda0']}`.",
        "",
        "## Arbitrary Columns",
        "",
        f"`{exact['arbitrary_columns_through_order11_lambda0']}`.",
        "",
        "## Boundary",
        "",
        artifact["proof_boundary"],
        "",
        "This is not a proof of RH.",
        "",
        "## Reproduce",
        "",
        "```powershell",
        f"python {GENERATOR_PATH}",
        f"python {CHECKER_PATH}",
        "```",
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
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(
        "wrote fixed order-eleven lambda-zero completion: "
        f"{ALL_SHIFT_Q11}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
