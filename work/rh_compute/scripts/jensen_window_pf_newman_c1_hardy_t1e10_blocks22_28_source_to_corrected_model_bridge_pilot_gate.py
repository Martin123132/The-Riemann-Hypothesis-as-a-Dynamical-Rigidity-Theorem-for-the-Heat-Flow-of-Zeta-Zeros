#!/usr/bin/env python3
"""Test every new block-22--28 bridge geometry before full promotion."""

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

import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
import jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_pilot_gate as bridge
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate as transport
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_retained_endpoint_gate as retained
from flint import ctx


TELEMETRY = bridge.TELEMETRY
ENDPOINTS = bridge.ENDPOINTS
WEIGHTS = bridge.WEIGHTS
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_pilot_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_pilot_gate_cache.jsonl"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_pilot_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_pilot_gate.py"
WITNESS_CHAINS = (
    862, 871, 992, 991,
    1278, 1291, 1306, 1305,
    1712, 1743,
    2186, 2127,
    2562, 2609,
    2990, 3047,
    3406, 3423,
)
EXPECTED_GROUP_COUNT = 18


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return bridge.decimal(value)


def load_selected() -> dict[int, dict[str, Any]]:
    selected = set(WITNESS_CHAINS)
    chains: dict[int, dict[str, Any]] = {}
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain = int(record.get("chain", 0))
        if chain not in selected:
            continue
        kind = record["type"]
        if kind == "chain":
            chains[chain] = {"header": record, "levels": {}, "q_terms": {}, "steps": {}}
        elif kind == "level":
            chains[chain]["levels"][int(record["level"])] = record
        elif kind == "q_terms":
            chains[chain]["q_terms"][int(record["nit"])] = record
        elif kind == "recurrence":
            chains[chain]["steps"][int(record["nit"])] = record
    require(set(chains) == selected, "blocks22--28 witness telemetry missing")
    for chain, data in chains.items():
        require(22 <= int(data["header"]["block"]) <= 28, f"chain {chain} block drift")
        require(int(data["header"]["mit"]) == 2, f"chain {chain} route drift")
        require(set(data["levels"]) >= {1, 2}, f"chain {chain} level payload missing")
        require(1 in data["q_terms"] and 1 in data["steps"], f"chain {chain} recurrence payload missing")
    return chains


def endpoint_rows() -> tuple[dict[int, dict[str, Any]], list[dict[str, Any]]]:
    artifact = json.loads(ENDPOINTS.read_text(encoding="utf-8"))
    later = [row for row in artifact["rows"] if 22 <= int(row["block"]) <= 28]
    require(len(later) == 733, "blocks22--28 endpoint roster drift")
    selected = {int(row["chain"]): row for row in later if int(row["chain"]) in WITNESS_CHAINS}
    require(set(selected) == set(WITNESS_CHAINS), "blocks22--28 endpoint witnesses missing")
    return selected, later


def geometry_key(row: dict[str, Any], levels: dict[int, dict[str, Any]]) -> tuple[Any, ...]:
    parent = levels[int(row["chain"])]["levels"][1]
    phi3 = cells.binary128_fraction(parent["coefficients_hex"][2])
    return (
        int(row["block"]),
        int(row["parent_length"]),
        int(row["child_length"]),
        row["selector"],
        "positive" if phi3 > 0 else "negative",
    )


def verify_preselection(
    later_rows: list[dict[str, Any]],
    factors: dict[tuple[int, int, int, int], Fraction],
) -> list[tuple[Any, ...]]:
    levels = retained.load_recursive()
    best: dict[tuple[Any, ...], tuple[Fraction, int]] = {}
    for row in later_rows:
        key = geometry_key(row, levels)
        factor = factors[(int(row["block"]), int(row["sum_index"]), bridge.OUTPUT_INDEX, int(row["branch"]))]
        score = factor * Fraction(row["complete_majorant_upper"])
        candidate = (score, int(row["chain"]))
        if key not in best or candidate > best[key]:
            best[key] = candidate
    require(len(best) == EXPECTED_GROUP_COUNT, "blocks22--28 geometry group count drift")
    expected = tuple(value[1] for _, value in sorted(best.items()))
    require(expected == WITNESS_CHAINS, f"blocks22--28 preselected roster drift: {expected}")
    return sorted(best)


def cache_fingerprint() -> str:
    payload = {
        "builder": file_hash(BUILDER),
        "checker": file_hash(CHECKER),
        "bridge_evaluator": file_hash(bridge.BUILDER),
        "telemetry": file_hash(TELEMETRY),
        "endpoints": file_hash(ENDPOINTS),
        "weights": file_hash(WEIGHTS),
        "witnesses": list(WITNESS_CHAINS),
        "precisions": list(bridge.PRECISIONS),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def write_header(fingerprint: str) -> None:
    CACHE.write_text(json.dumps({"kind": "blocks22_28_bridge_pilot_cache", "fingerprint": fingerprint}, sort_keys=True) + "\n", encoding="utf-8")


def append_cache(row: dict[str, Any]) -> None:
    with CACHE.open("a", encoding="utf-8", buffering=1) as handle:
        handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def load_cache(fingerprint: str) -> dict[int, dict[str, Any]]:
    if not CACHE.exists():
        write_header(fingerprint)
        return {}
    lines = [line for line in CACHE.read_text(encoding="utf-8").splitlines() if line.strip()]
    require(lines, "blocks22--28 pilot cache empty")
    header = json.loads(lines[0])
    if header.get("fingerprint") != fingerprint:
        stale = CACHE.with_name(f"{CACHE.stem}.stale.{header.get('fingerprint', 'unknown')[:12]}{CACHE.suffix}")
        require(not stale.exists(), f"stale cache destination already exists: {stale}")
        CACHE.replace(stale)
        write_header(fingerprint)
        return {}
    rows: dict[int, dict[str, Any]] = {}
    for line in lines[1:]:
        row = json.loads(line)
        chain = int(row["chain"])
        require(chain not in rows, f"duplicate blocks22--28 pilot chain {chain}")
        rows[chain] = row
    return rows


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy t=1e10 blocks 22--28 signed bridge geometry pilot

Date: 2026-08-09

Status: representative finite length-and-branch gate closed; not a proof and not a full 733-call certificate

The roster was fixed before exact integration by taking the largest retained
output-15 contribution in each occupied
`(block, parent length, child length, W3/W4, sign(Phi3))` class.  It contains
`{aggregate['witness_count']}` witnesses across parent lengths
`{', '.join(str(value) for value in aggregate['parent_lengths'])}`.

All `{aggregate['closed_witness_count']}` signed source-to-corrected-model
bridges close against the independent endpoint contour calculation.  All
`{aggregate['exact_corrections_within_retained_count']}` exact corrections are
inside the earlier retained majorants.  The largest exact/retained ratio is
`{aggregate['maximum_exact_to_retained_ratio_upper']}`.

The maximum zero-containing formula-target gap is
`{aggregate['maximum_formula_target_gap_upper']}`; this includes the deliberate
origin boxes on four rays and is not a fitted numerical residual.

Boundary: the pilot covers every new exact-point geometry present in blocks
22--28, but it does not replace the complete 733-call roster, coefficient-cell
transport, height uniformity, or the remaining analytic proof obligations.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-limit-seconds", type=float, default=10800.0)
    args = parser.parse_args()
    for path in (TELEMETRY, ENDPOINTS, WEIGHTS, CHECKER):
        require(path.is_file(), f"missing blocks22--28 pilot dependency: {path}")
    priority = bridge.set_low_priority()
    ctx.threads = 1
    started = time.monotonic()
    selected = load_selected()
    endpoints, later = endpoint_rows()
    factors = transport.load_weight_factors()
    geometry_groups = verify_preselection(later, factors)
    fingerprint = cache_fingerprint()
    cached = load_cache(fingerprint)
    require(set(cached).issubset(WITNESS_CHAINS), "blocks22--28 pilot cache roster drift")
    for chain in WITNESS_CHAINS:
        if chain in cached:
            continue
        if time.monotonic() - started >= args.runtime_limit_seconds:
            print(f"parked blocks22--28 bridge pilot after {len(cached)}/{len(WITNESS_CHAINS)} witnesses")
            return 75
        source = endpoints[chain]
        factor = factors[(int(source["block"]), int(source["sum_index"]), bridge.OUTPUT_INDEX, int(source["branch"]))]
        row = bridge.build_row(chain, selected[chain], source, factor)
        append_cache(row)
        cached[chain] = row
        print(f"closed later-block bridge witness {len(cached)}/{len(WITNESS_CHAINS)} (chain {chain})", flush=True)
    rows = [cached[chain] for chain in WITNESS_CHAINS]
    require(all(bool(row["signed_bridge_closed"]) for row in rows), "blocks22--28 pilot closure drift")
    require(all(bool(row["exact_correction_within_retained_majorant"]) for row in rows), "blocks22--28 pilot retained majorant failure")
    aggregate = {
        "witness_count": len(rows),
        "geometry_group_count": len(geometry_groups),
        "parent_lengths": sorted({int(row["parent_length"]) for row in rows}),
        "blocks": sorted({int(row["block"]) for row in rows}),
        "closed_witness_count": sum(bool(row["signed_bridge_closed"]) for row in rows),
        "exact_corrections_within_retained_count": sum(bool(row["exact_correction_within_retained_majorant"]) for row in rows),
        "maximum_exact_to_retained_ratio_upper": decimal(max(Fraction(row["exact_to_retained_ratio_upper"]) for row in rows)),
        "maximum_formula_target_gap_upper": decimal(max(Fraction(row["formula_target_gap_upper"]) for row in rows)),
        "minimum_integral_accuracy_bits": min(int(row["minimum_integral_accuracy_bits"]) for row in rows),
        "remaining_full_roster_call_count": 733,
    }
    passed = (
        aggregate["closed_witness_count"] == len(rows)
        and aggregate["exact_corrections_within_retained_count"] == len(rows)
        and aggregate["geometry_group_count"] == EXPECTED_GROUP_COUNT
    )
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_pilot_gate",
        "status": "all_new_later_block_exact_point_geometries_closed" if passed else "later_block_exact_point_geometry_pilot_falsified",
        "scope": "One preselected maximum-weight witness in each of 18 occupied exact-point geometry classes across blocks 22--28.",
        "selection_rule": "Maximum output-15 retained transported contribution in each (block, parent length, child length, selector, Phi3 sign) class.",
        "precisions_decimal_digits": list(bridge.PRECISIONS),
        "runtime": {
            "priority": priority,
            "flint_threads": 1,
            "worker_count": 1,
            "elapsed_seconds": time.monotonic() - started,
            "resumable_cache": relative(CACHE),
            "cache_fingerprint": fingerprint,
        },
        "provenance": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "endpoints": {"path": relative(ENDPOINTS), "sha256": file_hash(ENDPOINTS)},
            "weights": {"path": relative(WEIGHTS), "sha256": file_hash(WEIGHTS)},
            "bridge_evaluator": {"path": relative(bridge.BUILDER), "sha256": file_hash(bridge.BUILDER)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "geometry_groups": [list(group) for group in geometry_groups],
        "rows": rows,
        "aggregate": aggregate,
        "passed": passed,
        "boundary": {
            "proved": "Every distinct exact-point length/selector/cubic-sign geometry in blocks 22--28 has a closed maximum-weight witness.",
            "not_proved": "The other 715 calls, coefficient neighborhoods, arbitrary height, and RH remain open.",
            "next_action": "Run the same evaluator over all 733 blocks22--28 recursive calls, checkpointing each row, then replace their retained outer columns by exact magnitudes.",
        },
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built blocks22--28 bridge geometry pilot: "
        f"passed={passed}, witnesses={len(rows)}, max_ratio={aggregate['maximum_exact_to_retained_ratio_upper']}"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
