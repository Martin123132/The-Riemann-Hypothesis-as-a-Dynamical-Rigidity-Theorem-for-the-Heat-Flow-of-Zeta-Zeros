#!/usr/bin/env python3
"""Certify and transport all direct Hardy kernels in t=1e10 blocks 20--35."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import jensen_window_pf_newman_c1_hardy_block20_direct_kernel_gate as block20_direct
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate as transport
from flint import arb, ctx


TELEMETRY = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_crossblock_probe/enabled/chain_telemetry.jsonl"
ATLAS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate.json"
BLOCK20_DIRECT = block20_direct.RESULT
WEIGHT_FIXTURE = transport.WEIGHT_FIXTURE
WEIGHTS = transport.WEIGHTS
ENDPOINTS = transport.ENDPOINTS
RECURSIVE_OUTPUT = transport.RESULT
BLOCK20_TRANSPORT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_direct_kernel_output_gate_cache.jsonl"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_direct_kernel_output_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_direct_kernel_output_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_crossblock_direct_kernel_output_gate.py"
PRECISIONS = (70, 100)
EXPECTED_DIRECT = 5370
EXPECTED_NEW = 5320
REQUESTED_ERROR_SCALE = Fraction(5, 1000)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def magnitude_upper(value: Any) -> Fraction:
    return cells.bound_fraction(abs(value).upper())


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def load_direct() -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain = int(record.get("chain", 0))
        if record["type"] == "chain":
            chains[chain] = {"header": record, "levels": {}, "final": None}
        elif record["type"] == "level":
            chains[chain]["levels"][int(record["level"])] = record
        elif record["type"] == "chain_end":
            chains[chain]["final"] = record
    direct = {
        chain: data
        for chain, data in chains.items()
        if int(data["header"]["mit"]) == 1 and 20 <= int(data["header"]["block"]) <= 35
    }
    require(len(direct) == EXPECTED_DIRECT, "cross-block direct roster drift")
    return direct


def cache_fingerprint() -> str:
    material = {
        "telemetry": file_hash(TELEMETRY),
        "atlas": file_hash(ATLAS),
        "block20_direct": file_hash(BLOCK20_DIRECT),
        "point_ball_builder": file_hash(Path(point_balls.__file__).resolve()),
        "block20_direct_builder": file_hash(Path(block20_direct.__file__).resolve()),
        "builder": file_hash(BUILDER),
        "checker": file_hash(CHECKER),
        "precisions": list(PRECISIONS),
    }
    return hashlib.sha256(json.dumps(material, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_cache(fingerprint: str) -> dict[int, dict[str, Any]]:
    if not CACHE.exists():
        CACHE.write_text(
            json.dumps({"kind": "t1e10_crossblock_direct_kernel_cache", "fingerprint": fingerprint}, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return {}
    lines = [line for line in CACHE.read_text(encoding="utf-8").splitlines() if line.strip()]
    require(lines, "cross-block direct cache empty")
    header = json.loads(lines[0])
    require(header.get("fingerprint") == fingerprint, "cross-block direct cache fingerprint drift")
    rows: dict[int, dict[str, Any]] = {}
    for line in lines[1:]:
        row = json.loads(line)
        chain = int(row["chain"])
        require(chain not in rows, f"duplicate cached direct chain {chain}")
        rows[chain] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    with CACHE.open("a", encoding="utf-8", buffering=1) as handle:
        handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def evaluate(data: dict[str, Any], dps: int) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    header = data["header"]
    level = data["levels"][1]
    require(set(data["levels"]) == {1}, "MIT1 level shape drift")
    require(int(level["degree"]) == 3, "direct kernel degree drift")
    require(level["subtract_one"] is False, "unexpected MIT1 subtract-one route")
    source_tpm = point_balls.binary128_ball(header["tpm_hex"])
    source_kernel = point_balls.adapt(level, point_balls.direct_sum(level, source_tpm))
    logged_kernel = block20_direct.printed_complex_ball(data["final"]["final_state"])
    source_roundoff_gap = source_kernel - logged_kernel
    delta_tpm = -2 * arb.pi() - source_tpm
    phase_lipschitz = arb(0)
    for index in range(int(level["length"]) + 1):
        phase_lipschitz += abs(delta_tpm * block20_direct.polynomial(level, index))
    return {
        "source_kernel": source_kernel,
        "source_roundoff_gap": source_roundoff_gap,
        "phase_lipschitz": phase_lipschitz,
    }


def evaluate_row(chain: int, data: dict[str, Any]) -> dict[str, Any]:
    low = evaluate(data, PRECISIONS[0])
    high = evaluate(data, PRECISIONS[1])
    for name in ("source_kernel", "source_roundoff_gap", "phase_lipschitz"):
        require(low[name].overlaps(high[name]), f"chain {chain} {name} precision nonoverlap")
    source_upper = magnitude_upper(high["source_roundoff_gap"])
    phase_upper = magnitude_upper(high["phase_lipschitz"])
    total_upper = source_upper + phase_upper
    header = data["header"]
    level = data["levels"][1]
    return {
        "chain": chain,
        "block": int(header["block"]),
        "sum_index": int(header["sum_index"]),
        "branch": int(header["branch"]),
        "length_upper_index": int(level["length"]),
        "term_count": int(level["length"]) + 1,
        "conjugated": bool(level["conjugate"]),
        "source_roundoff_gap": point_balls.acb_record(high["source_roundoff_gap"], 55),
        "source_roundoff_gap_abs_upper": decimal(source_upper),
        "phase_normalization_lipschitz_upper": decimal(phase_upper),
        "total_true_gap_majorant_upper": decimal(total_upper),
        "precision_overlap": True,
    }


def admitted_block20_rows() -> dict[int, dict[str, Any]]:
    artifact = json.loads(BLOCK20_DIRECT.read_text(encoding="utf-8"))
    rows: dict[int, dict[str, Any]] = {}
    for source in artifact["rows"]:
        chain = int(source["chain"])
        rows[chain] = {
            "chain": chain,
            "block": 20,
            "sum_index": int(source["sum_index"]),
            "branch": int(source["branch"]),
            "length_upper_index": int(source["length_upper_index"]),
            "term_count": int(source["term_count"]),
            "conjugated": bool(source["conjugated"]),
            "source_roundoff_gap": source["source_roundoff_gap"],
            "source_roundoff_gap_abs_upper": source["source_roundoff_gap_abs_upper"],
            "phase_normalization_lipschitz_upper": source["phase_lipschitz_upper"],
            "total_true_gap_majorant_upper": source["total_true_gap_abs_upper"],
            "precision_overlap": True,
            "reused_admitted_block20_row": True,
        }
    require(len(rows) == 50, "admitted block-20 direct roster drift")
    return rows


def transport_outputs(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    ctx.dps = transport.PRECISION_DIGITS
    ctx.threads = 1
    factors = transport.load_weight_factors()
    by_block_output = {
        (block, output): Fraction(0)
        for block in range(20, 36)
        for output in range(1, 16)
    }
    for row in rows:
        block = int(row["block"])
        local = Fraction(row["total_true_gap_majorant_upper"])
        for output in range(1, 16):
            by_block_output[(block, output)] += factors[
                (block, int(row["sum_index"]), output, int(row["branch"]))
            ] * local

    block_rows: list[dict[str, Any]] = []
    for block in range(20, 36):
        selected = [row for row in rows if int(row["block"]) == block]
        values = [(by_block_output[(block, output)], output) for output in range(1, 16)]
        maximum = max(values)
        block_rows.append(
            {
                "block": block,
                "direct_call_count": len(selected),
                "upper_index": max(int(row["length_upper_index"]) for row in selected),
                "maximum_transported_direct_majorant_upper": decimal(maximum[0]),
                "maximum_witness_output": maximum[1],
            }
        )

    endpoint_artifact = json.loads(ENDPOINTS.read_text(encoding="utf-8"))
    recursive_by_output = {output: Fraction(0) for output in range(1, 16)}
    for row in endpoint_artifact["rows"]:
        block = int(row["block"])
        local = Fraction(row["complete_majorant_upper"])
        for output in range(1, 16):
            recursive_by_output[output] += factors[
                (block, int(row["sum_index"]), output, int(row["branch"]))
            ] * local

    admitted_recursive = json.loads(RECURSIVE_OUTPUT.read_text(encoding="utf-8"))
    admitted_recursive_rows = {
        int(row["output_index"]): row["all_recursive_blocks_majorant_upper"]
        for row in admitted_recursive["output_rows"]
    }
    for output in range(1, 16):
        require(
            decimal(recursive_by_output[output]) == admitted_recursive_rows[output],
            f"recursive output anchor mismatch at {output}",
        )

    block20_transport = json.loads(BLOCK20_TRANSPORT.read_text(encoding="utf-8"))
    block20_expected = {
        int(row["output_index"]): row["transported_direct_majorant_upper"]
        for row in block20_transport["rows"]
    }
    for output in range(1, 16):
        require(
            decimal(by_block_output[(20, output)]) == block20_expected[output],
            f"block-20 direct output anchor mismatch at {output}",
        )

    output_rows: list[dict[str, Any]] = []
    for output in range(1, 16):
        direct_total = sum((by_block_output[(block, output)] for block in range(20, 36)), Fraction(0))
        combined = recursive_by_output[output] + direct_total
        output_rows.append(
            {
                "output_index": output,
                "output_label": transport.output_label(output),
                "recursive_endpoint_majorant_upper": decimal(recursive_by_output[output]),
                "direct_kernel_majorant_upper": decimal(direct_total),
                "combined_exact_point_majorant_upper": decimal(combined),
                "combined_to_requested_scale_ratio": decimal(combined / REQUESTED_ERROR_SCALE),
                "within_requested_scale": combined < REQUESTED_ERROR_SCALE,
            }
        )

    max_direct = max(
        (Fraction(row["direct_kernel_majorant_upper"]), int(row["output_index"]))
        for row in output_rows
    )
    max_combined = max(
        (Fraction(row["combined_exact_point_majorant_upper"]), int(row["output_index"]))
        for row in output_rows
    )
    max_ratio = max(
        (Fraction(row["combined_to_requested_scale_ratio"]), int(row["output_index"]))
        for row in output_rows
    )
    aggregate = {
        "direct_call_count": len(rows),
        "recursive_call_count": len(endpoint_artifact["rows"]),
        "total_call_count": len(rows) + len(endpoint_artifact["rows"]),
        "output_count": 15,
        "requested_error_scale": decimal(REQUESTED_ERROR_SCALE),
        "recursive_output_anchor_exact": True,
        "block20_direct_output_anchor_exact": True,
        "maximum_transported_direct_majorant_upper": decimal(max_direct[0]),
        "maximum_transported_direct_majorant_witness_output": max_direct[1],
        "maximum_combined_exact_point_majorant_upper": decimal(max_combined[0]),
        "maximum_combined_exact_point_witness_output": max_combined[1],
        "maximum_combined_to_requested_scale_ratio": decimal(max_ratio[0]),
        "maximum_combined_to_requested_scale_ratio_witness_output": max_ratio[1],
        "outputs_within_requested_scale": sum(bool(row["within_requested_scale"]) for row in output_rows),
    }
    return block_rows, output_rows, aggregate


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    blocks = "\n".join(
        f"{row['block']:>5} {row['direct_call_count']:>6} {row['upper_index']:>5}  "
        f"{row['maximum_transported_direct_majorant_upper']}"
        for row in artifact["block_rows"]
    )
    return f"""# Hardy t=1e10 cross-block direct-kernel output gate

Date: 2026-08-09
Status: rigorous finite exact-point direct-kernel and combined output budget; not a proof of a height-uniform theorem or RH

Every direct kernel in blocks 20--35 is enclosed from its exact binary128
coefficients.  The exact source-phase finite sum is compared with the logged
40-decimal state at two Arb precisions.  The source-to-mathematical `2*pi`
shift is then bounded by `sum_n |Delta tpm * phi(n)|` using
`|exp(iu)-exp(iv)|<=|u-v)|`.

```text
block  calls  upper index  maximum transported direct bound
{blocks}
```

The complete 5370-call direct column is below
`{aggregate['maximum_transported_direct_majorant_upper']}`.  Adding it to the
independently anchored 1414-call recursive endpoint column gives
`{aggregate['maximum_combined_exact_point_majorant_upper']}`, or
`{aggregate['maximum_combined_to_requested_scale_ratio']}` times the nominal
0.005 scale.  Exactly `{aggregate['outputs_within_requested_scale']}/15`
outputs close under this unsigned exact-point budget.

This excludes coefficient-cell inflation outside block 20, source recurrence
arithmetic outside the admitted correction model, the outer Hardy remainder,
height uniformity, Lambda<=0, PF-infinity, RH, and a prize theorem.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-limit-seconds", type=float, default=10800.0)
    args = parser.parse_args()
    dependencies = (
        TELEMETRY,
        ATLAS,
        BLOCK20_DIRECT,
        WEIGHT_FIXTURE,
        WEIGHTS,
        ENDPOINTS,
        RECURSIVE_OUTPUT,
        BLOCK20_TRANSPORT,
        CHECKER,
    )
    for path in dependencies:
        require(path.is_file(), f"missing cross-block direct dependency: {path}")
    priority = set_low_priority()
    ctx.threads = 1
    started = time.monotonic()
    direct = load_direct()
    new_direct = {chain: data for chain, data in direct.items() if int(data["header"]["block"]) >= 21}
    require(len(new_direct) == EXPECTED_NEW, "new direct roster drift")
    fingerprint = cache_fingerprint()
    cached = load_cache(fingerprint)
    require(set(cached).issubset(new_direct), "direct cache contains an unexpected chain")
    for chain in sorted(new_direct):
        if chain in cached:
            continue
        row = evaluate_row(chain, new_direct[chain])
        append_cache(row)
        cached[chain] = row
        if time.monotonic() - started >= args.runtime_limit_seconds:
            print(f"parked t=1e10 cross-block direct replay after {len(cached)}/{EXPECTED_NEW} new calls")
            return 2

    require(len(cached) == EXPECTED_NEW, "cross-block direct cache incomplete")
    all_rows = admitted_block20_rows()
    all_rows.update(cached)
    require(len(all_rows) == EXPECTED_DIRECT, "combined direct roster drift")
    rows = [all_rows[chain] for chain in sorted(all_rows)]
    block_rows, output_rows, aggregate = transport_outputs(rows)

    local_max = max((Fraction(row["total_true_gap_majorant_upper"]), int(row["chain"])) for row in rows)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_crossblock_direct_kernel_output_gate",
        "status": "rigorous_t1e10_all_direct_kernels_and_complete_exact_point_output_budget",
        "scope": {
            "first_block": 20,
            "last_block": 35,
            "direct_calls": EXPECTED_DIRECT,
            "newly_evaluated_calls": EXPECTED_NEW,
            "reused_block20_calls": 50,
            "recursive_calls": 1414,
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "logged_state_contract": "Each ES50.40E4 component is enclosed by one full last-decimal unit.",
        },
        "local_aggregate": {
            "precision_overlap_count": sum(bool(row["precision_overlap"]) for row in rows),
            "maximum_total_true_gap_majorant_upper": decimal(local_max[0]),
            "maximum_total_true_gap_witness": local_max[1],
            "length_histogram": dict(sorted(Counter(int(row["length_upper_index"]) for row in rows).items())),
        },
        "block_rows": block_rows,
        "output_rows": output_rows,
        "aggregate": aggregate,
        "rows": rows,
        "runtime": {
            "priority": priority,
            "active_workers": 1,
            "elapsed_seconds": time.monotonic() - started,
            "cache_rows_available": len(cached),
        },
        "dependencies": {
            path.stem: {"path": relative(path), "sha256": file_hash(path)}
            for path in dependencies[:-1]
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "cache": {"path": relative(CACHE), "sha256": file_hash(CACHE)},
        },
        "next_obligation": (
            "Tighten the recursive endpoint column enough to close the complete exact-point output budget, "
            "then extend coefficient-cell and source-arithmetic transport through blocks 21--35."
        ),
        "proof_boundary": (
            "Rigorous finite exact-point direct-kernel and recursive-endpoint output budget at t=10^10 only. "
            "It excludes later coefficient-cell and source-arithmetic transport, the outer Hardy remainder, "
            "height uniformity, Lambda <= 0, PF-infinity, RH, and a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built t=1e10 cross-block direct-kernel output gate: "
        f"{EXPECTED_DIRECT} direct, combined max {aggregate['maximum_combined_exact_point_majorant_upper']}, "
        f"{aggregate['outputs_within_requested_scale']}/15 below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
