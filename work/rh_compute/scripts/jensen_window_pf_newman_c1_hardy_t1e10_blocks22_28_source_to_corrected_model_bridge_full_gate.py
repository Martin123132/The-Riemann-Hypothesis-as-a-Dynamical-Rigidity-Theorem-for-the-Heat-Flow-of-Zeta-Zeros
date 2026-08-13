#!/usr/bin/env python3
"""Promote the signed exact-point bridge to all recursive blocks 22--28 calls."""

from __future__ import annotations

import argparse
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

import jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_full_gate as block21_full
import jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_pilot_gate as bridge
import jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_pilot_gate as geometry_pilot
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate as transport
from flint import ctx


TELEMETRY = bridge.TELEMETRY
ENDPOINTS = bridge.ENDPOINTS
WEIGHTS = bridge.WEIGHTS
PILOT_RESULT = geometry_pilot.RESULT
BLOCK21_RESULT = block21_full.RESULT
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_full_gate_cache.jsonl"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_full_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_full_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_full_gate.py"
EXPECTED_CALLS = 733
REQUESTED_ERROR_SCALE = Fraction(5, 1000)
CPU_PARK_THRESHOLD = 75.0
ROW_CACHE_SCHEMA = "blocks22_28_bridge_rows_v1"
LEGACY_CACHE_FINGERPRINTS = {
    "8172987c1c2ddc1ac3643c527963bd879443e030d925c5a22a596c1c7efcdd62",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return bridge.decimal(value)


def load_later_chains() -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    active: set[int] = set()
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain = int(record.get("chain", 0))
        kind = record["type"]
        if kind == "chain":
            if 22 <= int(record["block"]) <= 28 and int(record["mit"]) == 2:
                active.add(chain)
                chains[chain] = {"header": record, "levels": {}, "q_terms": {}, "steps": {}}
            continue
        if chain not in active:
            continue
        if kind == "level":
            chains[chain]["levels"][int(record["level"])] = record
        elif kind == "q_terms":
            chains[chain]["q_terms"][int(record["nit"])] = record
        elif kind == "recurrence":
            chains[chain]["steps"][int(record["nit"])] = record
    require(len(chains) == EXPECTED_CALLS, f"blocks22--28 recursive roster drift: {len(chains)}")
    for chain, data in chains.items():
        require(set(data["levels"]) >= {1, 2}, f"chain {chain} level payload missing")
        require(1 in data["q_terms"] and 1 in data["steps"], f"chain {chain} recurrence payload missing")
        require(int(data["levels"][2]["length"]) == 1, f"chain {chain} child length drift")
    return chains


def load_endpoint_rows() -> dict[int, dict[str, Any]]:
    artifact = json.loads(ENDPOINTS.read_text(encoding="utf-8"))
    rows = {int(row["chain"]): row for row in artifact["rows"] if 22 <= int(row["block"]) <= 28}
    require(len(rows) == EXPECTED_CALLS, "blocks22--28 endpoint roster drift")
    return rows


def cache_fingerprint() -> str:
    payload = {
        "row_cache_schema": ROW_CACHE_SCHEMA,
        "bridge_evaluator": file_hash(bridge.BUILDER),
        "telemetry": file_hash(TELEMETRY),
        "endpoints": file_hash(ENDPOINTS),
        "weights": file_hash(WEIGHTS),
        "expected_calls": EXPECTED_CALLS,
        "precisions": list(bridge.PRECISIONS),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def write_header(fingerprint: str) -> None:
    CACHE.write_text(json.dumps({"kind": "blocks22_28_bridge_full_cache", "fingerprint": fingerprint}, sort_keys=True) + "\n", encoding="utf-8")


def append_cache(row: dict[str, Any]) -> None:
    with CACHE.open("a", encoding="utf-8", buffering=1) as handle:
        handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def seed_from_pilot() -> dict[int, dict[str, Any]]:
    artifact = json.loads(PILOT_RESULT.read_text(encoding="utf-8"))
    require(bool(artifact["passed"]), "cannot seed from failed later-block geometry pilot")
    rows: dict[int, dict[str, Any]] = {}
    for source in artifact["rows"]:
        row = dict(source)
        row["seeded_from_geometry_pilot"] = True
        append_cache(row)
        rows[int(row["chain"])] = row
    require(set(rows) == set(geometry_pilot.WITNESS_CHAINS), "later-block pilot seed roster drift")
    return rows


def load_cache(fingerprint: str) -> dict[int, dict[str, Any]]:
    if not CACHE.exists():
        write_header(fingerprint)
        return seed_from_pilot()
    lines = [line for line in CACHE.read_text(encoding="utf-8").splitlines() if line.strip()]
    require(lines, "blocks22--28 full cache empty")
    header = json.loads(lines[0])
    if header.get("fingerprint") != fingerprint:
        if header.get("fingerprint") in LEGACY_CACHE_FINGERPRINTS:
            migrated = CACHE.with_name(f"{CACHE.stem}.migrating{CACHE.suffix}")
            require(not migrated.exists(), f"cache migration file already exists: {migrated}")
            migrated.write_text(
                json.dumps({"kind": "blocks22_28_bridge_full_cache", "fingerprint": fingerprint}, sort_keys=True)
                + "\n"
                + "\n".join(lines[1:])
                + ("\n" if len(lines) > 1 else ""),
                encoding="utf-8",
            )
            os.replace(migrated, CACHE)
        else:
            stale = CACHE.with_name(f"{CACHE.stem}.stale.{header.get('fingerprint', 'unknown')[:12]}{CACHE.suffix}")
            require(not stale.exists(), f"stale cache destination already exists: {stale}")
            CACHE.replace(stale)
            write_header(fingerprint)
            return seed_from_pilot()
    rows: dict[int, dict[str, Any]] = {}
    for line in lines[1:]:
        row = json.loads(line)
        chain = int(row["chain"])
        require(chain not in rows, f"duplicate later-block full cache chain {chain}")
        rows[chain] = row
    return rows


def cpu_percent() -> float | None:
    try:
        import psutil

        return float(psutil.cpu_percent(interval=1.0))
    except Exception:
        return None


def promoted_outputs(
    rows: list[dict[str, Any]],
    endpoints: dict[int, dict[str, Any]],
    factors: dict[tuple[int, int, int, int], Fraction],
) -> list[dict[str, Any]]:
    block21 = json.loads(BLOCK21_RESULT.read_text(encoding="utf-8"))
    prior = {int(row["output_index"]): row for row in block21["output_rows"]}
    result: list[dict[str, Any]] = []
    for output in range(1, 16):
        retained = Fraction(0)
        exact = Fraction(0)
        block_exact = {block: Fraction(0) for block in range(22, 29)}
        for row in rows:
            chain = int(row["chain"])
            source = endpoints[chain]
            block = int(source["block"])
            factor = factors[(block, int(source["sum_index"]), output, int(source["branch"]))]
            retained += factor * Fraction(source["complete_majorant_upper"])
            contribution = factor * Fraction(row["exact_correction_magnitude_upper"])
            exact += contribution
            block_exact[block] += contribution
        old = prior[output]
        old_complete = Fraction(old["promoted_hybrid_complete_majorant_upper"])
        final_complete = old_complete - retained + exact
        require(final_complete >= 0, f"output {output} all-recursive exact budget became negative")
        result.append(
            {
                "output_index": output,
                "output_label": old["output_label"],
                "blocks22_28_retained_majorant_upper": decimal(retained),
                "blocks22_28_exact_correction_majorant_upper": decimal(exact),
                "blocks22_28_absolute_improvement": decimal(retained - exact),
                "blocks22_28_improvement_factor": decimal(retained / exact),
                "exact_correction_by_block": {str(block): decimal(block_exact[block]) for block in range(22, 29)},
                "prior_block21_promoted_complete_upper": decimal(old_complete),
                "all_recursive_exact_complete_majorant_upper": decimal(final_complete),
                "margin_below_requested_scale": decimal(REQUESTED_ERROR_SCALE - final_complete),
                "to_requested_scale_ratio": decimal(final_complete / REQUESTED_ERROR_SCALE),
                "within_requested_scale": final_complete < REQUESTED_ERROR_SCALE,
            }
        )
    return result


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy t=1e10 complete recursive exact-point bridge

Date: 2026-08-09

Status: finite all-recursive exact-point certificate; not a proof of a coefficient neighborhood, height uniformity, or RH

This gate promotes the blocks 22--28 geometry pilot to all
`{aggregate['later_block_call_count']}` remaining recursive calls.  Together
with the admitted 374 block-20 and 307 block-21 gates, the complete recursive
roster is `{aggregate['all_recursive_call_count']}` exact-point calls.

Every later-block signed source bridge closes independently, and every exact
correction lies inside its earlier retained majorant.  Outer transport sums
absolute local magnitudes, so no signed cross-call cancellation is used.

The blocks 22--28 retained outer column is tightened by at least a factor
`{aggregate['minimum_later_block_improvement_factor']}`.  The worst complete
output, including the 5,370 admitted direct calls, is reduced from
`{aggregate['maximum_prior_block21_promoted_complete_upper']}` to
`{aggregate['maximum_all_recursive_exact_complete_majorant_upper']}`.  Its
margin below `0.005` is
`{aggregate['minimum_margin_below_requested_scale']}`.

What this changes: finite exact-point uncertainty across the saved t=1e10
roster is no longer the active obstruction.  What it does not change:
coefficient-cell transport, the outer Hardy input neighborhood, arbitrary
height, and the analytic passage from those uniform bounds to RH remain open.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-limit-seconds", type=float, default=10800.0)
    args = parser.parse_args()
    for path in (TELEMETRY, ENDPOINTS, WEIGHTS, PILOT_RESULT, BLOCK21_RESULT, CHECKER):
        require(path.is_file(), f"missing blocks22--28 full dependency: {path}")
    priority = bridge.set_low_priority()
    ctx.threads = 1
    started = time.monotonic()
    chains = load_later_chains()
    endpoints = load_endpoint_rows()
    require(set(chains) == set(endpoints), "later-block telemetry/endpoint roster mismatch")
    factors = transport.load_weight_factors()
    fingerprint = cache_fingerprint()
    cached = load_cache(fingerprint)
    require(set(cached).issubset(chains), "later-block full cache roster drift")
    high_cpu_streak = 0
    parked_for_cpu = False
    pending = [chain for chain in sorted(chains) if chain not in cached]
    for chain in pending:
        if time.monotonic() - started >= args.runtime_limit_seconds:
            print(f"parked blocks22--28 full bridge after {len(cached)}/{EXPECTED_CALLS} calls (runtime limit)")
            return 75
        source = endpoints[chain]
        factor = factors[(int(source["block"]), int(source["sum_index"]), bridge.OUTPUT_INDEX, int(source["branch"]))]
        row = bridge.build_row(chain, chains[chain], source, factor)
        row["seeded_from_geometry_pilot"] = False
        append_cache(row)
        cached[chain] = row
        if len(cached) % 10 == 0 or len(cached) == EXPECTED_CALLS:
            print(f"closed later-block bridge call {len(cached)}/{EXPECTED_CALLS} (chain {chain})", flush=True)
            sampled = cpu_percent()
            if sampled is not None:
                high_cpu_streak = high_cpu_streak + 1 if sampled > CPU_PARK_THRESHOLD else 0
                if high_cpu_streak >= 2:
                    parked_for_cpu = True
                    break
    if parked_for_cpu:
        print(f"parked blocks22--28 full bridge after {len(cached)}/{EXPECTED_CALLS} calls (sustained CPU)")
        return 75
    require(len(cached) == EXPECTED_CALLS, "blocks22--28 full cache incomplete")
    rows = [cached[chain] for chain in sorted(cached)]
    require(all(bool(row["signed_bridge_closed"]) for row in rows), "later-block signed bridge closure drift")
    require(all(bool(row["exact_correction_within_retained_majorant"]) for row in rows), "later-block retained majorant failure")
    outputs = promoted_outputs(rows, endpoints, factors)
    maximum_complete = max((Fraction(row["all_recursive_exact_complete_majorant_upper"]), int(row["output_index"])) for row in outputs)
    minimum_margin = min((Fraction(row["margin_below_requested_scale"]), int(row["output_index"])) for row in outputs)
    minimum_improvement = min((Fraction(row["blocks22_28_improvement_factor"]), int(row["output_index"])) for row in outputs)
    aggregate = {
        "later_block_call_count": len(rows),
        "all_recursive_call_count": 374 + 307 + len(rows),
        "direct_call_count": 5370,
        "total_call_count": 374 + 307 + len(rows) + 5370,
        "seeded_geometry_pilot_call_count": sum(bool(row["seeded_from_geometry_pilot"]) for row in rows),
        "call_count_by_block": {str(block): sum(int(row["block"]) == block for row in rows) for block in range(22, 29)},
        "signed_bridge_closed_count": sum(bool(row["signed_bridge_closed"]) for row in rows),
        "exact_corrections_within_retained_count": sum(bool(row["exact_correction_within_retained_majorant"]) for row in rows),
        "maximum_exact_to_retained_ratio_upper": decimal(max(Fraction(row["exact_to_retained_ratio_upper"]) for row in rows)),
        "maximum_formula_target_gap_upper": decimal(max(Fraction(row["formula_target_gap_upper"]) for row in rows)),
        "minimum_integral_accuracy_bits": min(int(row["minimum_integral_accuracy_bits"]) for row in rows),
        "output_count": len(outputs),
        "minimum_later_block_improvement_factor": decimal(minimum_improvement[0]),
        "minimum_later_block_improvement_witness_output": minimum_improvement[1],
        "maximum_prior_block21_promoted_complete_upper": decimal(max(Fraction(row["prior_block21_promoted_complete_upper"]) for row in outputs)),
        "maximum_all_recursive_exact_complete_majorant_upper": decimal(maximum_complete[0]),
        "maximum_all_recursive_exact_witness_output": maximum_complete[1],
        "minimum_margin_below_requested_scale": decimal(minimum_margin[0]),
        "minimum_margin_witness_output": minimum_margin[1],
        "outputs_within_requested_scale": sum(bool(row["within_requested_scale"]) for row in outputs),
        "uses_signed_cross_call_cancellation": False,
        "remaining_recursive_exact_point_calls": 0,
    }
    passed = (
        aggregate["signed_bridge_closed_count"] == EXPECTED_CALLS
        and aggregate["exact_corrections_within_retained_count"] == EXPECTED_CALLS
        and aggregate["outputs_within_requested_scale"] == 15
        and aggregate["all_recursive_call_count"] == 1414
    )
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_full_gate",
        "status": "all_1414_recursive_exact_point_bridges_and_outer_transport_closed" if passed else "later_block_full_bridge_falsified",
        "scope": "All 733 recursive calls in blocks 22--28, completing all 1,414 recursive exact-point calls and 15 outer outputs at t=1e10.",
        "precisions_decimal_digits": list(bridge.PRECISIONS),
        "runtime": {
            "priority": priority,
            "flint_threads": 1,
            "worker_count": 1,
            "elapsed_seconds": time.monotonic() - started,
            "resumable_cache": relative(CACHE),
            "cache_fingerprint": fingerprint,
            "cpu_park_threshold_percent": CPU_PARK_THRESHOLD,
        },
        "provenance": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "endpoints": {"path": relative(ENDPOINTS), "sha256": file_hash(ENDPOINTS)},
            "weights": {"path": relative(WEIGHTS), "sha256": file_hash(WEIGHTS)},
            "geometry_pilot": {"path": relative(PILOT_RESULT), "sha256": file_hash(PILOT_RESULT)},
            "block21_full": {"path": relative(BLOCK21_RESULT), "sha256": file_hash(BLOCK21_RESULT)},
            "bridge_evaluator": {"path": relative(bridge.BUILDER), "sha256": file_hash(bridge.BUILDER)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "rows": rows,
        "output_rows": outputs,
        "aggregate": aggregate,
        "passed": passed,
        "boundary": {
            "proved": "Complete finite exact-point signed source bridge for all 1,414 recursive calls, with absolute-magnitude transport through all 15 outer outputs.",
            "not_proved": "No coefficient-neighborhood, outer Hardy coefficient-cell, arbitrary-height, or RH conclusion follows from the finite exact-point roster alone.",
            "next_action": "Use the exact-point corrections as centers for coefficient-cell transport, beginning with the worst correction/retained-ratio witnesses and proving selector-stable ray and special-function enclosures over their local cells.",
        },
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built complete later-block bridge gate: "
        f"passed={passed}, calls={len(rows)}, all_recursive_max={aggregate['maximum_all_recursive_exact_complete_majorant_upper']}"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
