#!/usr/bin/env python3
"""Scale the signed corrected-model coefficient cells to all later recursive calls."""

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

import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
import jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_pilot_gate as pilot
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate as transport
from flint import ctx


TELEMETRY = pilot.TELEMETRY
BLOCK21_POINTS = pilot.BLOCK21_POINTS
LATER_POINTS = pilot.LATER_POINTS
WEIGHTS = transport.WEIGHTS
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_full_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_full_gate_cache.jsonl"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_full_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_full_gate.py"
EXPECTED_CALLS = 1040
EXPECTED_BY_BLOCK = {21: 307, 22: 240, 23: 185, 24: 134, 25: 90, 26: 61, 27: 21, 28: 2}
REQUESTED_SCALE = Fraction(5, 1000)
CPU_PARK_THRESHOLD = 75.0


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


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


def cpu_percent() -> float | None:
    try:
        import psutil

        return float(psutil.cpu_percent(interval=1.0))
    except Exception:  # pragma: no cover
        return None


def load_point_artifacts() -> tuple[dict[int, dict[str, Any]], dict[str, Any]]:
    block21 = json.loads(BLOCK21_POINTS.read_text(encoding="utf-8"))
    later = json.loads(LATER_POINTS.read_text(encoding="utf-8"))
    require(bool(block21["passed"]) and bool(later["passed"]), "exact-point dependency failed")
    rows = {int(row["chain"]): row for row in block21["rows"]}
    for row in later["rows"]:
        chain = int(row["chain"])
        require(chain not in rows, f"duplicate exact-point chain {chain}")
        rows[chain] = row
    require(len(rows) == EXPECTED_CALLS, "exact-point later roster drift")
    return rows, later


def load_chains(selected: set[int]) -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    with TELEMETRY.open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            chain = int(record.get("chain", 0))
            if chain not in selected:
                continue
            if record["type"] == "chain":
                chains[chain] = {"header": record, "levels": {}, "q_terms": {}, "steps": {}}
            elif record["type"] == "level":
                chains[chain]["levels"][int(record["level"])] = record
            elif record["type"] == "q_terms":
                chains[chain]["q_terms"][int(record["nit"])] = record
            elif record["type"] == "recurrence":
                chains[chain]["steps"][int(record["nit"])] = record
    require(set(chains) == selected, "later coefficient-cell telemetry roster missing")
    return chains


def cache_fingerprint() -> str:
    payload = {
        "builder": file_hash(BUILDER),
        "checker": file_hash(CHECKER),
        "pilot_builder": file_hash(pilot.BUILDER),
        "telemetry": file_hash(TELEMETRY),
        "block21_points": file_hash(BLOCK21_POINTS),
        "later_points": file_hash(LATER_POINTS),
        "weights": file_hash(WEIGHTS),
        "expected_calls": EXPECTED_CALLS,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_cache(fingerprint: str) -> dict[int, dict[str, Any]]:
    if not CACHE.exists():
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        CACHE.write_text(json.dumps({"kind": "later_cell_bridge_full_cache", "fingerprint": fingerprint}, sort_keys=True) + "\n", encoding="utf-8")
        return {}
    lines = [line for line in CACHE.read_text(encoding="utf-8").splitlines() if line.strip()]
    require(lines, "empty full coefficient-cell cache")
    header = json.loads(lines[0])
    if header.get("fingerprint") != fingerprint:
        stale = CACHE.with_name(f"{CACHE.stem}.stale.{header.get('fingerprint', 'unknown')[:12]}{CACHE.suffix}")
        require(not stale.exists(), f"stale full cache destination exists: {stale}")
        CACHE.replace(stale)
        CACHE.write_text(json.dumps({"kind": "later_cell_bridge_full_cache", "fingerprint": fingerprint}, sort_keys=True) + "\n", encoding="utf-8")
        return {}
    rows: dict[int, dict[str, Any]] = {}
    for line in lines[1:]:
        row = json.loads(line)
        chain = int(row["chain"])
        require(chain not in rows, f"duplicate full cache chain {chain}")
        rows[chain] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    with CACHE.open("a", encoding="utf-8", buffering=1) as handle:
        handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def output_rows(
    rows: list[dict[str, Any]],
    point_rows: dict[int, dict[str, Any]],
    point_artifact: dict[str, Any],
) -> list[dict[str, Any]]:
    factors = transport.load_weight_factors()
    prior = {int(row["output_index"]): row for row in point_artifact["output_rows"]}
    result: list[dict[str, Any]] = []
    for output in range(1, 16):
        point_column = Fraction(0)
        cell_column = Fraction(0)
        by_block: dict[int, Fraction] = {block: Fraction(0) for block in EXPECTED_BY_BLOCK}
        for row in rows:
            chain = int(row["chain"])
            point = point_rows[chain]
            key = (int(row["block"]), int(row["sum_index"]), output, int(row["branch"]))
            factor = factors[key]
            point_column += factor * Fraction(point["exact_correction_magnitude_upper"])
            contribution = factor * Fraction(row["cell_correction_magnitude_upper"])
            cell_column += contribution
            by_block[int(row["block"])] += contribution
        prior_complete = Fraction(prior[output]["all_recursive_exact_complete_majorant_upper"])
        cell_complete = prior_complete - point_column + cell_column
        require(cell_complete >= 0, f"output {output} cell complete budget became negative")
        result.append(
            {
                "output_index": output,
                "output_label": prior[output]["output_label"],
                "later_recursive_point_correction_upper": decimal(point_column),
                "later_recursive_cell_correction_upper": decimal(cell_column),
                "later_recursive_cell_inflation_upper": decimal(cell_column - point_column),
                "cell_correction_by_block": {str(block): decimal(value) for block, value in by_block.items()},
                "prior_all_recursive_exact_complete_majorant_upper": decimal(prior_complete),
                "all_recursive_cell_complete_majorant_upper": decimal(cell_complete),
                "margin_below_requested_scale": decimal(REQUESTED_SCALE - cell_complete),
                "within_requested_scale": cell_complete < REQUESTED_SCALE,
            }
        )
    return result


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy t=1e10 complete later-recursive signed coefficient-cell bridge

Date: 2026-08-09

Status: complete finite corrected-model cell column at one saved height; not a complete recurrence, height-uniform theorem, or RH proof

The two-witness pilot is scaled to all `{aggregate['call_count']}` recursive
calls in blocks 21--28.  Every row has a nonzero-radius selector-stable box,
overlapping 70/110-digit Arb enclosures of the signed analytic correction
`Qexact-Qpaper`, and an enclosure containing the saved exact-point correction.

The calculation keeps

```text
Dcorr = Ppi - I69pi - Qpaper = Qexact - Qpaper
```

as one complex interval before taking its magnitude.  It uses no signed
cancellation between different calls.

Cells closed: `{aggregate['closed_cell_count']}/{aggregate['call_count']}`

Maximum cell/point inflation: `{aggregate['maximum_cell_to_point_inflation_upper']}`

Worst all-recursive cell complete output: `{aggregate['maximum_all_recursive_cell_complete_majorant_upper']}`

Minimum margin below `0.005`: `{aggregate['minimum_margin_below_requested_scale']}`

Boundary: this closes the analytic corrected-model coefficient-cell column at
the finite saved height only.  Source binary128 rounding over boxes, the
complete cubic recurrence tail in blocks 21--28, the outer Hardy
representation/remainder, and height-uniform selector and accumulation
theorems remain open.  It does not prove Lambda<=0, PF-infinity, RH, or a
prize-level conclusion.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-limit-seconds", type=float, default=10800.0)
    args = parser.parse_args()
    for path in (TELEMETRY, BLOCK21_POINTS, LATER_POINTS, WEIGHTS, CHECKER):
        require(path.is_file(), f"missing full coefficient-cell dependency: {path}")
    started = time.monotonic()
    priority = set_low_priority()
    ctx.threads = 1
    points, point_artifact = load_point_artifacts()
    chains = load_chains(set(points))
    fingerprint = cache_fingerprint()
    cached = load_cache(fingerprint)
    require(set(cached).issubset(points), "full cache contains unexpected chain")
    pending = [chain for chain in sorted(points) if chain not in cached]
    high_cpu_streak = 0
    for index, chain in enumerate(pending, 1):
        if time.monotonic() - started >= args.runtime_limit_seconds:
            print(f"parked complete later coefficient-cell bridge after {len(cached)}/{EXPECTED_CALLS} calls")
            return 75
        row = pilot.build_row(chain, chains[chain], points[chain])
        append_cache(row)
        cached[chain] = row
        if len(cached) % 10 == 0 or len(cached) == EXPECTED_CALLS:
            print(f"closed later coefficient cell {len(cached)}/{EXPECTED_CALLS} (chain {chain})", flush=True)
            sampled = cpu_percent()
            if sampled is not None:
                high_cpu_streak = high_cpu_streak + 1 if sampled > CPU_PARK_THRESHOLD else 0
                if high_cpu_streak >= 2 and len(cached) < EXPECTED_CALLS:
                    print(f"parked complete later coefficient-cell bridge after {len(cached)}/{EXPECTED_CALLS} calls (sustained CPU)")
                    return 75
    require(len(cached) == EXPECTED_CALLS, "full coefficient-cell cache incomplete")
    rows = [cached[chain] for chain in sorted(cached)]
    outputs = output_rows(rows, points, point_artifact)
    maximum_inflation = max((Fraction(row["cell_to_point_inflation_upper"]), int(row["chain"])) for row in rows)
    maximum_extra = max((int(row["transport_extra_halvings"]), int(row["chain"])) for row in rows)
    maximum_output = max((Fraction(row["all_recursive_cell_complete_majorant_upper"]), int(row["output_index"])) for row in outputs)
    minimum_margin = min((Fraction(row["margin_below_requested_scale"]), int(row["output_index"])) for row in outputs)
    block_histogram = Counter(int(row["block"]) for row in rows)
    aggregate = {
        "call_count": len(rows),
        "call_count_by_block": {str(block): block_histogram[block] for block in EXPECTED_BY_BLOCK},
        "w3_count": sum(row["selector"] == "W3" for row in rows),
        "w4_count": sum(row["selector"] == "W4" for row in rows),
        "negative_phi3_count": sum(row["phi3_sign"] == "negative" for row in rows),
        "positive_phi3_count": sum(row["phi3_sign"] == "positive" for row in rows),
        "closed_cell_count": sum(bool(row["signed_bridge_closed_on_analytic_cell"]) for row in rows),
        "precision_overlap_count": sum(bool(row["precision_overlap"]) for row in rows),
        "point_overlap_count": sum(bool(row["point_correction_overlap"]) for row in rows),
        "nonzero_radius_cell_count": sum(bool(row["radii_nonzero"]) for row in rows),
        "within_retained_majorant_count": sum(bool(row["cell_correction_within_retained_majorant"]) for row in rows),
        "maximum_cell_to_point_inflation_upper": decimal(maximum_inflation[0]),
        "maximum_cell_to_point_inflation_witness_chain": maximum_inflation[1],
        "maximum_transport_extra_halvings": maximum_extra[0],
        "maximum_transport_extra_halvings_witness_chain": maximum_extra[1],
        "output_count": len(outputs),
        "outputs_within_requested_scale": sum(bool(row["within_requested_scale"]) for row in outputs),
        "maximum_all_recursive_cell_complete_majorant_upper": decimal(maximum_output[0]),
        "maximum_all_recursive_cell_complete_witness_output": maximum_output[1],
        "minimum_margin_below_requested_scale": decimal(minimum_margin[0]),
        "minimum_margin_witness_output": minimum_margin[1],
        "uses_signed_cross_call_cancellation": False,
        "remaining_later_corrected_model_point_only_calls": 0,
    }
    passed = (
        aggregate["call_count_by_block"] == {str(block): count for block, count in EXPECTED_BY_BLOCK.items()}
        and aggregate["closed_cell_count"] == EXPECTED_CALLS
        and aggregate["precision_overlap_count"] == EXPECTED_CALLS
        and aggregate["point_overlap_count"] == EXPECTED_CALLS
        and aggregate["nonzero_radius_cell_count"] == EXPECTED_CALLS
        and aggregate["within_retained_majorant_count"] == EXPECTED_CALLS
        and aggregate["outputs_within_requested_scale"] == 15
    )
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_full_gate",
        "status": "all_1040_later_recursive_signed_coefficient_cells_and_15_outputs_closed" if passed else "later_recursive_signed_coefficient_cell_full_gate_falsified",
        "scope": "All recursive calls in blocks 21 through 28 and all 15 outer outputs at the saved t=1e10 fixture.",
        "passed": passed,
        "precisions_decimal_digits": list(pilot.PRECISIONS),
        "signed_bridge_identity": "Dcorr=Ppi-I69pi-Qpaper=Qexact-Qpaper",
        "rows": rows,
        "output_rows": outputs,
        "aggregate": aggregate,
        "runtime": {
            "priority": priority,
            "flint_threads": 1,
            "worker_count": 1,
            "elapsed_seconds": time.monotonic() - started,
            "resumable_cache": relative(CACHE),
            "cache_fingerprint": fingerprint,
            "cpu_park_threshold_percent": CPU_PARK_THRESHOLD,
        },
        "dependencies": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "block21_points": {"path": relative(BLOCK21_POINTS), "sha256": file_hash(BLOCK21_POINTS)},
            "later_points": {"path": relative(LATER_POINTS), "sha256": file_hash(LATER_POINTS)},
            "weights": {"path": relative(WEIGHTS), "sha256": file_hash(WEIGHTS)},
            "pilot_builder": {"path": relative(pilot.BUILDER), "sha256": file_hash(pilot.BUILDER)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": "Add complete cubic recurrence-tail and source-rounding transport on the 1,040 accepted cells, then certify the outer Hardy representation and height-uniform route theorem.",
        "proof_boundary": "Finite saved-height analytic corrected-model cell column only; no complete later recurrence, outer Hardy remainder, height uniformity, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    require(passed, "complete later coefficient-cell bridge did not close")
    tmp = RESULT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(RESULT)
    note_tmp = NOTE.with_suffix(".md.tmp")
    note_tmp.write_text(render_note(artifact), encoding="utf-8")
    note_tmp.replace(NOTE)
    print(
        "built complete later coefficient-cell signed bridge: "
        f"{aggregate['closed_cell_count']}/{EXPECTED_CALLS} cells, "
        f"15/15 outputs, max {aggregate['maximum_all_recursive_cell_complete_majorant_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
