#!/usr/bin/env python3
"""Promote the signed exact-point bridge to every recursive block-21 call."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_pilot_gate as pilot
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate as transport
from flint import ctx


TELEMETRY = pilot.TELEMETRY
ENDPOINTS = pilot.ENDPOINTS
WEIGHTS = pilot.WEIGHTS
PILOT_RESULT = pilot.RESULT
HYBRID_RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_hybrid_exact_output_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_full_gate_cache.jsonl"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_full_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_full_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_full_gate.py"
EXPECTED_CALLS = 307
REQUESTED_ERROR_SCALE = Fraction(5, 1000)
CPU_PARK_THRESHOLD = 75.0
ROW_CACHE_SCHEMA = "block21_bridge_rows_v1"
LEGACY_CACHE_FINGERPRINTS = {
    "028453772aa8c3790b565472f91c57fc65ceafb1d99612471d9563c746305a15",
    "d7f6bd32275f8f39b96cda8d2c68f23271f3a7544e81f80f0acad50c2cf8afab",
    "5403f64323ee3935ff992d83460226efcbf6b66b2904b3db5e0a861222e1a7c1",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return pilot.decimal(value)


def set_low_priority() -> str:
    return pilot.set_low_priority()


def load_block21() -> dict[int, dict[str, Any]]:
    chains: dict[int, dict[str, Any]] = {}
    active: set[int] = set()
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain = int(record.get("chain", 0))
        kind = record["type"]
        if kind == "chain":
            if int(record["block"]) == 21 and int(record["mit"]) == 2:
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
    require(len(chains) == EXPECTED_CALLS, f"block-21 recursive roster drift: {len(chains)}")
    for chain, data in chains.items():
        require(set(data["levels"]) >= {1, 2}, f"chain {chain} level payload missing")
        require(1 in data["q_terms"] and 1 in data["steps"], f"chain {chain} recurrence payload missing")
        require(int(data["levels"][1]["length"]) == 119, f"chain {chain} parent length drift")
    return chains


def load_endpoint_rows() -> dict[int, dict[str, Any]]:
    artifact = json.loads(ENDPOINTS.read_text(encoding="utf-8"))
    rows = {int(row["chain"]): row for row in artifact["rows"] if int(row["block"]) == 21}
    require(len(rows) == EXPECTED_CALLS, "block-21 endpoint roster drift")
    return rows


def cache_fingerprint() -> str:
    payload = {
        "row_cache_schema": ROW_CACHE_SCHEMA,
        "pilot_builder": file_hash(pilot.BUILDER),
        "telemetry": file_hash(TELEMETRY),
        "endpoints": file_hash(ENDPOINTS),
        "weights": file_hash(WEIGHTS),
        "expected_calls": EXPECTED_CALLS,
        "precisions": list(pilot.PRECISIONS),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def write_cache_header(fingerprint: str) -> None:
    CACHE.write_text(
        json.dumps({"kind": "block21_bridge_full_cache", "fingerprint": fingerprint}, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def append_cache(row: dict[str, Any]) -> None:
    with CACHE.open("a", encoding="utf-8", buffering=1) as handle:
        handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def load_cache(fingerprint: str) -> dict[int, dict[str, Any]]:
    if not CACHE.exists():
        write_cache_header(fingerprint)
        return seed_from_pilot()
    lines = [line for line in CACHE.read_text(encoding="utf-8").splitlines() if line.strip()]
    require(lines, "block-21 full bridge cache empty")
    header = json.loads(lines[0])
    if header.get("fingerprint") != fingerprint:
        if header.get("fingerprint") in LEGACY_CACHE_FINGERPRINTS:
            migrated = CACHE.with_name(f"{CACHE.stem}.migrating{CACHE.suffix}")
            require(not migrated.exists(), f"cache migration file already exists: {migrated}")
            migrated.write_text(
                json.dumps({"kind": "block21_bridge_full_cache", "fingerprint": fingerprint}, sort_keys=True)
                + "\n"
                + "\n".join(lines[1:])
                + ("\n" if len(lines) > 1 else ""),
                encoding="utf-8",
            )
            os.replace(migrated, CACHE)
            header = {"fingerprint": fingerprint}
        else:
            stale = CACHE.with_name(f"{CACHE.stem}.stale.{header.get('fingerprint', 'unknown')[:12]}{CACHE.suffix}")
            require(not stale.exists(), f"stale cache destination already exists: {stale}")
            CACHE.replace(stale)
            write_cache_header(fingerprint)
            return seed_from_pilot()
    rows: dict[int, dict[str, Any]] = {}
    for line in lines[1:]:
        row = json.loads(line)
        chain = int(row["chain"])
        require(chain not in rows, f"duplicate full bridge cache chain {chain}")
        rows[chain] = row
    return rows


def seed_from_pilot() -> dict[int, dict[str, Any]]:
    artifact = json.loads(PILOT_RESULT.read_text(encoding="utf-8"))
    require(bool(artifact["passed"]), "cannot seed from a failed block-21 pilot")
    rows: dict[int, dict[str, Any]] = {}
    for source in artifact["rows"]:
        row = dict(source)
        row["seeded_from_pilot"] = True
        append_cache(row)
        rows[int(row["chain"])] = row
    require(set(rows) == set(pilot.WITNESS_CHAINS), "block-21 pilot seed roster drift")
    return rows


def cpu_percent() -> float | None:
    try:
        import psutil

        return float(psutil.cpu_percent(interval=1.0))
    except Exception:
        return None


def output_rows(
    rows: list[dict[str, Any]],
    endpoint_rows: dict[int, dict[str, Any]],
    factors: dict[tuple[int, int, int, int], Fraction],
) -> list[dict[str, Any]]:
    hybrid = json.loads(HYBRID_RESULT.read_text(encoding="utf-8"))
    prior = {int(row["output_index"]): row for row in hybrid["output_rows"]}
    result: list[dict[str, Any]] = []
    for output in range(1, 16):
        retained = Fraction(0)
        exact = Fraction(0)
        for row in rows:
            chain = int(row["chain"])
            source = endpoint_rows[chain]
            factor = factors[(21, int(source["sum_index"]), output, int(source["branch"]))]
            retained += factor * Fraction(source["complete_majorant_upper"])
            exact += factor * Fraction(row["exact_correction_magnitude_upper"])
        old = prior[output]
        old_hybrid = Fraction(old["hybrid_complete_majorant_upper"])
        old_later = Fraction(old["blocks21_28_retained_majorant_upper"])
        new_later = old_later - retained + exact
        new_hybrid = old_hybrid - retained + exact
        require(new_later >= 0 and new_hybrid >= 0, f"output {output} promoted budget became negative")
        result.append(
            {
                "output_index": output,
                "output_label": old["output_label"],
                "block21_retained_majorant_upper": decimal(retained),
                "block21_exact_correction_majorant_upper": decimal(exact),
                "block21_absolute_improvement": decimal(retained - exact),
                "block21_improvement_factor": decimal(retained / exact),
                "blocks21_28_after_block21_promotion_upper": decimal(new_later),
                "prior_hybrid_complete_majorant_upper": decimal(old_hybrid),
                "promoted_hybrid_complete_majorant_upper": decimal(new_hybrid),
                "margin_below_requested_scale": decimal(REQUESTED_ERROR_SCALE - new_hybrid),
                "to_requested_scale_ratio": decimal(new_hybrid / REQUESTED_ERROR_SCALE),
                "within_requested_scale": new_hybrid < REQUESTED_ERROR_SCALE,
            }
        )
    return result


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy t=1e10 complete block-21 signed bridge gate

Date: 2026-08-09

Status: complete finite exact-point block-21 signed bridge; not a proof of a coefficient neighborhood, arbitrary height, or RH

The representative four-branch pilot is promoted here to every one of the
`{aggregate['call_count']}` recursive calls in block 21.  Each row independently
closes the signed source correction chain against the sign-aware endpoint
contour identity.  No correction is inferred from an outer output and no
cross-call cancellation is used.

```text
Dsrc -> exact W1 source replacement -> source-only t2 deletion
     -> pi transport -> total source/paper W2--W5 transport
     = Qexact - Qpaper.
```

Every exact local correction lies inside the earlier retained-decay majorant:
`{aggregate['exact_corrections_within_retained_count']}/{aggregate['call_count']}`.
The worst exact/retained ratio is
`{aggregate['maximum_exact_to_retained_ratio_upper']}`.

Transporting exact magnitudes through all 15 outer outputs lowers the worst
hybrid complete column from
`{aggregate['maximum_prior_hybrid_complete_majorant_upper']}` to
`{aggregate['maximum_promoted_hybrid_complete_majorant_upper']}`.  The minimum
remaining margin below `0.005` is
`{aggregate['minimum_margin_below_requested_scale']}`.

Boundary: this closes the exact-point source bridge for block 21 only.  The
remaining finite exact-point frontier is the 733 recursive calls in blocks
22--28.  Coefficient-cell transport and height uniformity remain separate
proof obligations after that finite roster.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-limit-seconds", type=float, default=10800.0)
    args = parser.parse_args()
    for path in (TELEMETRY, ENDPOINTS, WEIGHTS, PILOT_RESULT, HYBRID_RESULT, CHECKER):
        require(path.is_file(), f"missing full block-21 bridge dependency: {path}")
    priority = set_low_priority()
    ctx.threads = 1
    started = time.monotonic()
    chains = load_block21()
    endpoints = load_endpoint_rows()
    require(set(chains) == set(endpoints), "block-21 telemetry/endpoint roster mismatch")
    factors = transport.load_weight_factors()
    fingerprint = cache_fingerprint()
    cached = load_cache(fingerprint)
    require(set(cached).issubset(chains), "full bridge cache contains unexpected chain")
    high_cpu_streak = 0
    parked_for_cpu = False
    pending = [chain for chain in sorted(chains) if chain not in cached]
    for index, chain in enumerate(pending, 1):
        if time.monotonic() - started >= args.runtime_limit_seconds:
            print(f"parked full block-21 bridge after {len(cached)}/{EXPECTED_CALLS} calls (runtime limit)")
            return 75
        source = endpoints[chain]
        factor = factors[(21, int(source["sum_index"]), pilot.OUTPUT_INDEX, int(source["branch"]))]
        row = pilot.build_row(chain, chains[chain], source, factor)
        row["seeded_from_pilot"] = False
        append_cache(row)
        cached[chain] = row
        if len(cached) % 10 == 0 or len(cached) == EXPECTED_CALLS:
            print(f"closed block-21 bridge call {len(cached)}/{EXPECTED_CALLS} (chain {chain})", flush=True)
            sampled = cpu_percent()
            if sampled is not None:
                high_cpu_streak = high_cpu_streak + 1 if sampled > CPU_PARK_THRESHOLD else 0
                if high_cpu_streak >= 2:
                    parked_for_cpu = True
                    break
    if parked_for_cpu:
        print(f"parked full block-21 bridge after {len(cached)}/{EXPECTED_CALLS} calls (sustained CPU)")
        return 75
    require(len(cached) == EXPECTED_CALLS, "full block-21 bridge cache incomplete")
    rows = [cached[chain] for chain in sorted(cached)]
    require(all(bool(row["signed_bridge_closed"]) for row in rows), "full block-21 signed bridge closure drift")
    require(all(bool(row["exact_correction_within_retained_majorant"]) for row in rows), "full block-21 retained majorant failure")
    outputs = output_rows(rows, endpoints, factors)
    maximum_promoted = max((Fraction(row["promoted_hybrid_complete_majorant_upper"]), int(row["output_index"])) for row in outputs)
    minimum_margin = min((Fraction(row["margin_below_requested_scale"]), int(row["output_index"])) for row in outputs)
    maximum_ratio = max((Fraction(row["exact_to_retained_ratio_upper"]), int(row["chain"])) for row in rows)
    aggregate = {
        "call_count": len(rows),
        "seeded_pilot_call_count": sum(bool(row["seeded_from_pilot"]) for row in rows),
        "w3_call_count": sum(row["selector"] == "W3" for row in rows),
        "w4_call_count": sum(row["selector"] == "W4" for row in rows),
        "negative_phi3_call_count": sum(row["phi3_sign"] == "negative" for row in rows),
        "positive_phi3_call_count": sum(row["phi3_sign"] == "positive" for row in rows),
        "signed_bridge_closed_count": sum(bool(row["signed_bridge_closed"]) for row in rows),
        "exact_corrections_within_retained_count": sum(bool(row["exact_correction_within_retained_majorant"]) for row in rows),
        "maximum_exact_to_retained_ratio_upper": decimal(maximum_ratio[0]),
        "maximum_exact_to_retained_ratio_witness_chain": maximum_ratio[1],
        "maximum_formula_target_gap_upper": decimal(max(Fraction(row["formula_target_gap_upper"]) for row in rows)),
        "minimum_integral_accuracy_bits": min(int(row["minimum_integral_accuracy_bits"]) for row in rows),
        "output_count": len(outputs),
        "maximum_prior_hybrid_complete_majorant_upper": decimal(max(Fraction(row["prior_hybrid_complete_majorant_upper"]) for row in outputs)),
        "maximum_promoted_hybrid_complete_majorant_upper": decimal(maximum_promoted[0]),
        "maximum_promoted_hybrid_witness_output": maximum_promoted[1],
        "minimum_margin_below_requested_scale": decimal(minimum_margin[0]),
        "minimum_margin_witness_output": minimum_margin[1],
        "outputs_within_requested_scale": sum(bool(row["within_requested_scale"]) for row in outputs),
        "uses_signed_cross_call_cancellation": False,
        "remaining_recursive_exact_point_calls_blocks22_28": 733,
    }
    passed = (
        aggregate["signed_bridge_closed_count"] == EXPECTED_CALLS
        and aggregate["exact_corrections_within_retained_count"] == EXPECTED_CALLS
        and aggregate["outputs_within_requested_scale"] == 15
    )
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_full_gate",
        "status": "complete_block21_exact_point_signed_bridge_and_outer_transport_closed" if passed else "complete_block21_exact_point_signed_bridge_falsified",
        "scope": "All 307 recursive block-21 calls and all 15 outer output points at t=1e10.",
        "precisions_decimal_digits": list(pilot.PRECISIONS),
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
            "pilot_result": {"path": relative(PILOT_RESULT), "sha256": file_hash(PILOT_RESULT)},
            "prior_hybrid": {"path": relative(HYBRID_RESULT), "sha256": file_hash(HYBRID_RESULT)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "rows": rows,
        "output_rows": outputs,
        "aggregate": aggregate,
        "passed": passed,
        "boundary": {
            "proved": "Complete finite exact-point signed bridge for block 21, transported by absolute local magnitudes through all outer weights.",
            "not_proved": "Blocks 22--28, coefficient neighborhoods, outer Hardy coefficient cells, arbitrary height, and RH remain open.",
            "next_action": "Generalize the same length-aware bridge to blocks 22--28, beginning with one witness per distinct parent length and contour branch before running the remaining 733-call roster.",
        },
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built complete block-21 bridge gate: "
        f"passed={passed}, calls={len(rows)}, max={aggregate['maximum_promoted_hybrid_complete_majorant_upper']}"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
