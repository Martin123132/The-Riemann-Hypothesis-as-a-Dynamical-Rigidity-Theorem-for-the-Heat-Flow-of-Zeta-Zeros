#!/usr/bin/env python3
"""Compose delayed signed order-ten entry at lambda=-100."""

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
    "jensen_window_pf_compound_order10_m100_delayed_entry_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order10_m100_delayed_entry_certificate.md"
)
CURVATURE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_first_summand_curvature_certificate.json"
)
TRANSFER_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_first_summand_curvature_bridge.json"
)
FINITE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_m100_finite_splice_certificate.json"
)
TAIL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_m100_tail_curvature_reduction.json"
)

CONTINUOUS_THEOREM = "z_1''(t)<=4200/t^2 for every real t>=1251"
FIRST_DISCRETE_THEOREM = (
    "z_1''(t)<=4200/t^2 => Z_k^(1)<=4200*[-log(1-1/k^2)]"
    "<4201/k^2, k>=1252"
)
FULL_TRANSFER_THEOREM = (
    "|Z_k-Z_k^(1)|<10/k^2 for every integer k>=1252"
)
FULL_CEILING_THEOREM = (
    "z_1''(t)<=4200/t^2 on t>=1251 => "
    "Z_k<4201/k^2+10/k^2=4211/k^2<5500/k^2, k>=1252"
)
TAIL_THEOREM = "Q_(10,n)(-100)>0 for every n>=1243"
FINITE_POSITIVE_THEOREM = "Q_(10,n)(-100)>0 for every 4<=n<=1242"
NEGATIVE_PREFIX_THEOREM = "Q_(10,n)(-100)<0 for n=0,1,2,3"
DELAYED_ENTRY_THEOREM = "Q_(10,n)(-100)>0 for every integer n>=4"
SIGN_CHART_THEOREM = (
    "Q_(10,n)(-100)<0 for 0<=n<=3 and "
    "Q_(10,n)(-100)>0 for every n>=4"
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
    tail = load_json(TAIL_SOURCE)

    if curvature.get("theorem") != CONTINUOUS_THEOREM:
        raise RuntimeError("global order-ten first-summand theorem changed")
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
            "z_1''(t)<=4200/t^2 on t>=1251 => " + TAIL_THEOREM
        ),
        "preserved_negative_prefix": NEGATIVE_PREFIX_THEOREM,
    }
    for key, expected in expected_transfer.items():
        if exact_transfer.get(key) != expected:
            raise RuntimeError(f"order-ten transfer {key} changed")
    if transfer.get("summary", {}).get("full_kernel_transfer_theorems") != 1:
        raise RuntimeError("order-ten full-kernel transfer is not closed")

    finite_exact = finite.get("exact", {})
    if finite_exact.get("combined_positive_block") != FINITE_POSITIVE_THEOREM:
        raise RuntimeError("finite order-ten positive block changed")
    if finite_exact.get("preserved_negative_prefix") != NEGATIVE_PREFIX_THEOREM:
        raise RuntimeError("finite order-ten negative prefix changed")
    if finite.get("summary", {}).get("open_rows") != 0:
        raise RuntimeError("finite order-ten source has open rows")

    tail_exact = tail.get("exact", {})
    if tail_exact.get("sufficient_ceiling") != (
        "Z_k<=5500/k^2 for every integer k>=1252"
    ):
        raise RuntimeError("order-ten tail ceiling changed")
    if tail_exact.get("sign_equivalence") != (
        "Q_(10,n)>0 iff L_n>0 iff E_n<0"
    ):
        raise RuntimeError("order-ten tail sign equivalence changed")
    if tail_exact.get("finite_splice") != (
        "[Q_(10,n)(-100)>0 for 4<=n<=1242] plus "
        "[Q_(10,n)(-100)>0 for n>=1243] implies "
        "Q_(10,n)(-100)>0 for every n>=4"
    ):
        raise RuntimeError("order-ten finite splice changed")

    return {
        "sources": [
            source_record(CURVATURE_SOURCE, curvature),
            source_record(TRANSFER_SOURCE, transfer),
            source_record(FINITE_SOURCE, finite),
            source_record(TAIL_SOURCE, tail),
        ],
        "continuous": curvature["theorem"],
        "tent_transfer": exact_transfer["tent_transfer"],
        "full_transfer": exact_transfer["full_transfer"],
        "full_ceiling": exact_transfer["conditional_full_ceiling"],
        "tail": TAIL_THEOREM,
        "finite_positive": finite_exact["combined_positive_block"],
        "negative_prefix": finite_exact["preserved_negative_prefix"],
        "orientation": finite_exact["orientation"],
    }


def build_artifact() -> dict:
    sources = validate_sources()
    rows = [
        EntryRow(
            "co10m100dec_01_continuous",
            "theorem_input",
            "ready_to_apply",
            "The global first-summand theorem closes the sole analytic premise in the transfer bridge.",
            sources["continuous"],
            "First Newman summand at lambda=-100.",
        ),
        EntryRow(
            "co10m100dec_02_tent",
            "exact_theorem_application",
            "ready_to_apply",
            "The tent identity converts the continuous ceiling to a discrete first-summand ceiling.",
            sources["tent_transfer"],
            "Integer k>=1252 only.",
        ),
        EntryRow(
            "co10m100dec_03_full_kernel",
            "exact_theorem_application",
            "ready_to_apply",
            "The twelfth-power perturbation envelope keeps the full curvature far below its endpoint target.",
            sources["full_transfer"] + "; " + sources["full_ceiling"],
            "Complete Newman kernel at lambda=-100, k>=1252.",
        ),
        EntryRow(
            "co10m100dec_04_analytic_tail",
            "analytic_theorem",
            "ready_to_apply",
            "The full curvature ceiling proves every order-ten endpoint sign in the analytic tail.",
            sources["tail"],
            "Endpoint shifts n>=1243.",
        ),
        EntryRow(
            "co10m100dec_05_finite_chart",
            "interval_theorem",
            "ready_to_apply",
            "The finite Arb certificate supplies the positive block and the four negative prefix signs.",
            sources["finite_positive"] + "; " + sources["negative_prefix"],
            "Endpoint shifts 0<=n<=1242.",
        ),
        EntryRow(
            "co10m100dec_06_delayed_entry",
            "exact_theorem_composition",
            "ready_to_apply",
            "The finite positive block and analytic tail meet without an index gap.",
            DELAYED_ENTRY_THEOREM,
            "Delayed positivity begins at n=4; this is not all-shift positivity.",
        ),
        EntryRow(
            "co10m100dec_07_sign_chart",
            "exact_theorem_composition",
            "ready_to_apply",
            "The endpoint signs are completely classified without erasing the failed prefix.",
            SIGN_CHART_THEOREM,
            "Fixed lambda=-100 endpoint only.",
            {"orientation": sources["orientation"]},
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order10_m100_delayed_entry_certificate",
        "date": "2026-07-16",
        "status": "rigorous delayed signed order-ten entry and complete endpoint sign chart at lambda=-100",
        "proof_boundary": (
            "This artifact proves the fixed order-ten sign chart at lambda=-100: "
            "four negative prefix rows and positivity from n=4 onward. It does not "
            "prove all-shift order-ten positivity, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": sources,
        "exact": {
            "continuous_first_summand": CONTINUOUS_THEOREM,
            "first_discrete_ceiling": FIRST_DISCRETE_THEOREM,
            "full_kernel_transfer": FULL_TRANSFER_THEOREM,
            "full_curvature_ceiling": FULL_CEILING_THEOREM,
            "analytic_tail": TAIL_THEOREM,
            "finite_positive_block": FINITE_POSITIVE_THEOREM,
            "negative_prefix": NEGATIVE_PREFIX_THEOREM,
            "delayed_entry": DELAYED_ENTRY_THEOREM,
            "complete_sign_chart": SIGN_CHART_THEOREM,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "continuous_curvature_inputs": 1,
            "full_kernel_transfer_applications": 1,
            "analytic_tail_theorems": 1,
            "finite_positive_block_theorems": 1,
            "negative_prefix_theorems": 1,
            "delayed_m100_entry_theorems": 1,
            "complete_endpoint_sign_charts": 1,
            "all_shift_positive_entry_theorems": 0,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order10_m100_delayed_entry_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order10_m100_delayed_entry_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Order-Ten Lambda=-100 Delayed Entry Certificate",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous delayed signed order-ten entry and complete endpoint",
        "sign chart at `lambda=-100`. This is not all-shift positivity,",
        "and this certificate is not a proof of PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        exact["continuous_first_summand"],
        exact["first_discrete_ceiling"],
        exact["full_kernel_transfer"],
        exact["full_curvature_ceiling"],
        exact["analytic_tail"],
        "```",
        "",
        "The analytic tail joins the finite positive block exactly after",
        "`n=1242`, proving",
        "",
        "```text",
        exact["delayed_entry"],
        exact["complete_sign_chart"],
        "```",
        "",
        "The rows `n=0,1,2,3` remain negative and are part of the theorem, not",
        "a numerical exception to be discarded.",
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
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print("wrote delayed signed order-ten entry at lambda=-100")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
