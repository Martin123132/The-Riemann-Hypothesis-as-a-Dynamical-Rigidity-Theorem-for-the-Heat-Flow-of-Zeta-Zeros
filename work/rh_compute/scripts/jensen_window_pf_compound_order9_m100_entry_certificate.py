#!/usr/bin/env python3
"""Compose all-shift signed order-nine entry at lambda=-100."""

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
    "jensen_window_pf_compound_order9_m100_entry_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order9_m100_entry_certificate.md"
)
CURVATURE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_first_summand_curvature_certificate.json"
)
TRANSFER_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_first_summand_curvature_bridge.json"
)
FINITE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_m100_finite_splice_certificate.json"
)
FLOW_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_uniform_tail_flow_reduction.json"
)

CONTINUOUS_THEOREM = "w_1''(t)<=4200/t^2 for every real t>=1250"
FIRST_DISCRETE_THEOREM = (
    "w_1''(t)<=4200/t^2 => Y_k^(1)<=4200*[-log(1-1/k^2)]"
    "<4201/k^2, k>=1251"
)
FULL_TRANSFER_THEOREM = (
    "|Y_k-Y_k^(1)|<=F_(k-1)+2*F_k+F_(k+1)<550/k^2, k>=1251"
)
FULL_CEILING_THEOREM = (
    "w_1''(t)<=4200/t^2 on t>=1250 => "
    "Y_k<4201/k^2+550/k^2=4751/k^2<4900/k^2, k>=1251"
)
TAIL_THEOREM = "Q_(9,n)(-100)>0 for every n>=1243"
FINITE_THEOREM = "Q_(9,n)(-100)>0 for every 0<=n<=1242"
ENTRY_THEOREM = (
    "Q_(9,n)(-100)=H_(9,n)(-100)>0 for every integer n>=0"
)


@dataclass(frozen=True)
class EntryRow:
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
    curvature = load_json(CURVATURE_SOURCE)
    transfer = load_json(TRANSFER_SOURCE)
    finite = load_json(FINITE_SOURCE)
    flow = load_json(FLOW_SOURCE)
    if curvature.get("theorem") != CONTINUOUS_THEOREM:
        raise RuntimeError("global first-summand curvature theorem changed")
    if curvature.get("summary", {}).get(
        "global_first_summand_curvature_theorems"
    ) != 1:
        raise RuntimeError("global first-summand theorem is not closed")
    exact_transfer = transfer.get("exact", {})
    expected_transfer = {
        "continuous_target": CONTINUOUS_THEOREM,
        "tent_transfer": FIRST_DISCRETE_THEOREM,
        "full_transfer": FULL_TRANSFER_THEOREM,
        "conditional_full_ceiling": FULL_CEILING_THEOREM,
        "conditional_endpoint_tail": (
            "w_1''(t)<=4200/t^2 on t>=1250 => " + TAIL_THEOREM
        ),
    }
    for key, expected in expected_transfer.items():
        if exact_transfer.get(key) != expected:
            raise RuntimeError(f"order-nine transfer {key} changed")
    if transfer.get("summary", {}).get("full_kernel_transfer_theorems") != 1:
        raise RuntimeError("full-kernel transfer is not closed")
    if finite.get("exact", {}).get("combined_prefix") != FINITE_THEOREM:
        raise RuntimeError("finite endpoint prefix changed")
    if finite.get("summary", {}).get("finite_splice_theorems") != 1:
        raise RuntimeError("finite endpoint splice is not closed")
    exact_flow = flow.get("exact", {})
    if exact_flow.get("order9_orientation") != "epsilon_9=1, Q_(9,n)=H_(9,n)":
        raise RuntimeError("order-nine sign orientation changed")
    if exact_flow.get("open_entry") != (
        "Q_(9,n)(-100)>0 for every n>=0, equivalently strict "
        "log-concavity of Q_(8,n)(-100)"
    ):
        raise RuntimeError("order-nine endpoint target changed")
    return {
        "sources": [
            source_record(CURVATURE_SOURCE, curvature),
            source_record(TRANSFER_SOURCE, transfer),
            source_record(FINITE_SOURCE, finite),
            source_record(FLOW_SOURCE, flow),
        ],
        "continuous": curvature["theorem"],
        "tent_transfer": exact_transfer["tent_transfer"],
        "full_transfer": exact_transfer["full_transfer"],
        "full_ceiling": exact_transfer["conditional_full_ceiling"],
        "tail": TAIL_THEOREM,
        "finite": finite["exact"]["combined_prefix"],
        "orientation": exact_flow["order9_orientation"],
    }


def build_artifact() -> dict:
    sources = validate_sources()
    rows = [
        EntryRow(
            "co9m100ec_01_continuous",
            "theorem_input",
            "ready_to_apply",
            "The global first-summand theorem supplies the sole open analytic premise.",
            sources["continuous"],
            "Continuous first summand at lambda=-100.",
        ),
        EntryRow(
            "co9m100ec_02_tent",
            "exact_theorem_application",
            "ready_to_apply",
            "The tent identity converts the continuous ceiling to a discrete first-summand ceiling.",
            sources["tent_transfer"],
            "Integer k>=1251 only.",
        ),
        EntryRow(
            "co9m100ec_03_full_kernel",
            "exact_theorem_application",
            "ready_to_apply",
            "The tenth-power remainder transfer keeps the complete-kernel curvature below its endpoint target.",
            sources["full_transfer"] + "; " + sources["full_ceiling"],
            "Complete theta kernel at lambda=-100, k>=1251.",
        ),
        EntryRow(
            "co9m100ec_04_analytic_tail",
            "analytic_theorem",
            "ready_to_apply",
            "The full curvature ceiling proves every signed endpoint determinant in the analytic tail.",
            sources["tail"],
            "Endpoint shifts n>=1243.",
        ),
        EntryRow(
            "co9m100ec_05_finite_prefix",
            "interval_theorem",
            "ready_to_apply",
            "The finite Arb certificate covers every endpoint shift before the analytic tail.",
            sources["finite"],
            "Endpoint shifts 0<=n<=1242.",
        ),
        EntryRow(
            "co9m100ec_06_all_shift_entry",
            "exact_theorem_composition",
            "ready_to_apply",
            "The finite prefix and analytic tail meet without an index gap.",
            ENTRY_THEOREM,
            "All shifts at lambda=-100 only; no heat-forward claim here.",
            {"orientation": sources["orientation"]},
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order9_m100_entry_certificate",
        "date": "2026-07-14",
        "status": "rigorous all-shift signed order-nine entry at lambda=-100",
        "proof_boundary": (
            "This artifact proves fixed order-nine endpoint entry at lambda=-100. "
            "It does not by itself prove heat-forward invariance, order ten or "
            "higher, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": sources,
        "exact": {
            "continuous_first_summand": CONTINUOUS_THEOREM,
            "first_discrete_ceiling": FIRST_DISCRETE_THEOREM,
            "full_kernel_transfer": FULL_TRANSFER_THEOREM,
            "full_curvature_ceiling": FULL_CEILING_THEOREM,
            "analytic_tail": TAIL_THEOREM,
            "finite_prefix": FINITE_THEOREM,
            "all_shift_entry": ENTRY_THEOREM,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "continuous_curvature_inputs": 1,
            "full_kernel_transfer_applications": 1,
            "analytic_tail_theorems": 1,
            "finite_prefix_theorems": 1,
            "all_shift_m100_entry_theorems": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_m100_entry_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order9_m100_entry_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Order-nine lambda=-100 entry certificate",
        "",
        "Date: 2026-07-14",
        "",
        "Status: rigorous all-shift signed order-nine entry at `lambda=-100`; This is not a proof of PF-infinity, RH, or `Lambda <= 0`; it closes",
        "one fixed compound order at one heat endpoint.",
        "",
        "## Analytic tail",
        "",
        f"`{exact['continuous_first_summand']}`",
        "",
        "feeds the tent and full-kernel bounds",
        "",
        f"`{exact['first_discrete_ceiling']}`",
        "",
        f"`{exact['full_kernel_transfer']}`",
        "",
        f"`{exact['full_curvature_ceiling']}`",
        "",
        "and therefore proves",
        "",
        f"`{exact['analytic_tail']}`.",
        "",
        "## Exact splice",
        "",
        f"The finite interval certificate proves `{exact['finite_prefix']}`.",
        "The two ranges meet exactly between `n=1242` and `n=1243`, giving",
        "",
        f"`{exact['all_shift_entry']}`.",
        "",
        "The existing cooperative-flow theorem can now propagate this endpoint",
        "sign through `-100<=lambda<=0`; that composition remains a separate",
        "fixed-order-nine artifact.",
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
    print("wrote all-shift signed order-nine entry at lambda=-100")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
