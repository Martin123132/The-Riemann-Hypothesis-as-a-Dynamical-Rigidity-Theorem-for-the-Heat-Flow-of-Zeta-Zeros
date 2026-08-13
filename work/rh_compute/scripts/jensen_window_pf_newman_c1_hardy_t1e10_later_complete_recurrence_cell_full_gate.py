#!/usr/bin/env python3
"""Scale the complete recurrence tail and source roundoff over all later cells."""

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

import jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_complete_cubic_legendre_tail_budget_gate as cell_tail
import jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate as coefficient_transport
import jensen_window_pf_newman_c1_hardy_block20_complete_cubic_legendre_tail_budget_gate as point_tail
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_full_gate as later_cells
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate as transport
import jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_stress_cell_pilot_gate as pilot
from flint import arb, ctx


FULL_CELLS = later_cells.RESULT
TELEMETRY = later_cells.TELEMETRY
WEIGHTS = transport.WEIGHTS
ROUNDING_MODE_PROBE = pilot.ROUNDING_MODE_PROBE
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_cell_full_gate.json"
CACHE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_cell_full_gate_cache.jsonl"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_cell_full_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_cell_full_gate.py"
EXPECTED_CALLS = later_cells.EXPECTED_CALLS
EXPECTED_BY_BLOCK = later_cells.EXPECTED_BY_BLOCK
PRECISIONS = pilot.PRECISIONS
REQUESTED_SCALE = pilot.REQUESTED_SCALE
CPU_PARK_THRESHOLD = 75.0


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


def cpu_percent() -> float | None:
    try:
        import psutil

        return float(psutil.cpu_percent(interval=1.0))
    except Exception:  # pragma: no cover
        return None


def cache_fingerprint() -> str:
    payload = {
        "builder": file_hash(BUILDER),
        "pilot_builder": file_hash(pilot.BUILDER),
        "full_cells": file_hash(FULL_CELLS),
        "telemetry": file_hash(TELEMETRY),
        "weights": file_hash(WEIGHTS),
        "rounding_mode_probe": file_hash(ROUNDING_MODE_PROBE),
        "expected_calls": EXPECTED_CALLS,
        "precisions": list(PRECISIONS),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_cache(fingerprint: str) -> dict[int, dict[str, Any]]:
    if not CACHE.exists():
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        CACHE.write_text(json.dumps({"kind": "later_complete_recurrence_cell_full_cache", "fingerprint": fingerprint}, sort_keys=True) + "\n", encoding="utf-8")
        return {}
    lines = [line for line in CACHE.read_text(encoding="utf-8").splitlines() if line.strip()]
    require(lines, "empty later complete-recurrence cache")
    header = json.loads(lines[0])
    if header.get("fingerprint") != fingerprint:
        stale = CACHE.with_name(f"{CACHE.stem}.stale.{header.get('fingerprint', 'unknown')[:12]}{CACHE.suffix}")
        require(not stale.exists(), f"stale later complete-recurrence cache destination exists: {stale}")
        CACHE.replace(stale)
        CACHE.write_text(json.dumps({"kind": "later_complete_recurrence_cell_full_cache", "fingerprint": fingerprint}, sort_keys=True) + "\n", encoding="utf-8")
        return {}
    rows: dict[int, dict[str, Any]] = {}
    for line in lines[1:]:
        row = json.loads(line)
        chain = int(row["chain"])
        require(chain not in rows, f"duplicate later complete-recurrence cache chain {chain}")
        rows[chain] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    with CACHE.open("a", encoding="utf-8", buffering=1) as handle:
        handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def build_row(chain: int, data: dict[str, Any], source_row: dict[str, Any]) -> dict[str, Any]:
    cell = source_row["refined_cell"]
    atlas_row = {"q_cell": cell}
    box = coefficient_transport.coefficient_box(data, atlas_row)
    low = cell_tail.evaluate_cell(data, atlas_row, PRECISIONS[0])
    high = cell_tail.evaluate_cell(data, atlas_row, PRECISIONS[1])
    require(len(low["term_gaps"]) == len(high["term_gaps"]) == 2, f"chain {chain} tail term roster drift")
    require(all(a.overlaps(b) for a, b in zip(low["term_gaps"], high["term_gaps"])), f"chain {chain} tail precision nonoverlap")
    point = point_tail.evaluate_chain(data, PRECISIONS[1])
    require(point["local_majorant"] > 0, f"chain {chain} zero point tail")
    require(point["local_majorant"] <= high["local_majorant"], f"chain {chain} point tail escaped cell tail")

    structural = coefficient_transport.structural_transport(data, atlas_row, box, PRECISIONS[1])
    parent_rounding = pilot.parent_rounding_guards(data, box)
    child_rounding = [
        pilot.rounding_preimage_inside_arb(payload, enclosure)
        for payload, enclosure in zip(data["levels"][2]["coefficients_hex"], structural["child_coefficients"])
    ]
    require(all(child_rounding), f"chain {chain} child rounding preimage escaped reconstructed cell")
    multiplier_hex = data["steps"][1]["multiplier_hex"]
    multiplier_rounding = [
        pilot.rounding_preimage_inside_arb(multiplier_hex[0], structural["multiplier"].real),
        pilot.rounding_preimage_inside_arb(multiplier_hex[1], structural["multiplier"].imag),
    ]
    require(all(multiplier_rounding), f"chain {chain} multiplier rounding preimage escaped anchored cell")
    state_hex = data["steps"][1]["state_before_hex"]
    state_rounding = [
        pilot.rounding_preimage_inside_arb(state_hex[0], structural["child_kernel"].real),
        pilot.rounding_preimage_inside_arb(state_hex[1], structural["child_kernel"].imag),
    ]
    require(all(state_rounding), f"chain {chain} state rounding preimage escaped child kernel")

    source_q = point_balls.binary128_complex(data["steps"][1]["qq_hex"])
    operation_roundoff = pilot.complex_multiply_add_roundoff_upper(
        structural["multiplier"], structural["child_kernel"], source_q
    )
    exact_operation, emitted_operation = pilot.exact_complex_multiply_add(data)
    observed_gap = abs(emitted_operation[0] - exact_operation[0]) + abs(emitted_operation[1] - exact_operation[1])
    point_roundoff = pilot.complex_multiply_add_roundoff_upper(
        point_balls.binary128_complex(multiplier_hex),
        point_balls.binary128_complex(state_hex),
        source_q,
    )
    require(observed_gap <= point_roundoff <= operation_roundoff, f"chain {chain} multiply-add roundoff drift")

    tpm_rounding = pilot.binary128_rounding_preimage(data["header"]["tpm_hex"])
    ctx.dps = PRECISIONS[1]
    minus_two_pi = -2 * arb.pi()
    tpm_contains = tpm_rounding["lower"] <= pilot.lower(minus_two_pi) and pilot.upper(minus_two_pi) <= tpm_rounding["upper"]
    require(tpm_contains, f"chain {chain} tpm rounding cell missed -2*pi")
    return {
        "chain": chain,
        "block": int(source_row["block"]),
        "sum_index": int(source_row["sum_index"]),
        "branch": int(source_row["branch"]),
        "parent_length": int(source_row["parent_length"]),
        "child_length": int(source_row["child_length"]),
        "selector": source_row["selector"],
        "phi3_sign": source_row["phi3_sign"],
        "precision_overlap": True,
        "point_tail_contained_in_cell_tail": True,
        "maximum_dimensionless_z_abs_upper": point_tail.upper_decimal(high["maximum_z"]),
        "minimum_stationary_discriminant_lower": point_tail.lower_decimal(high["minimum_discriminant"]),
        "maximum_factored_legendre_phase_tail_upper": point_tail.upper_decimal(high["maximum_tail"]),
        "point_complete_tail_majorant_upper": point_tail.upper_decimal(point["local_majorant"]),
        "cell_complete_tail_majorant_upper": point_tail.upper_decimal(high["local_majorant"]),
        "cell_to_point_tail_inflation_upper": point_tail.upper_decimal(high["local_majorant"] / point["local_majorant"]),
        "parent_binary128_rounding_preimages_inside_cell": bool(parent_rounding["all_preimages_inside_cell"]),
        "parent_maximum_rounding_radius_to_cell_radius_ratio": parent_rounding["maximum_rounding_radius_to_cell_radius_ratio"],
        "child_binary128_rounding_preimages_inside_reconstructed_cell": all(child_rounding),
        "multiplier_binary128_rounding_preimages_inside_anchored_cell": all(multiplier_rounding),
        "child_state_binary128_rounding_preimages_inside_cell_kernel": all(state_rounding),
        "tpm_rounding_preimage_contains_minus_two_pi": tpm_contains,
        "source_complex_multiply_add_observed_gap": decimal(observed_gap),
        "source_complex_multiply_add_point_roundoff_upper": decimal(point_roundoff),
        "source_complex_multiply_add_cell_roundoff_upper": decimal(operation_roundoff),
    }


def output_rows(rows: list[dict[str, Any]], full: dict[str, Any]) -> list[dict[str, Any]]:
    factors = transport.load_weight_factors()
    prior_outputs = {int(row["output_index"]): row for row in full["output_rows"]}
    outputs: list[dict[str, Any]] = []
    for output_index in range(1, 16):
        tail = Fraction(0)
        roundoff = Fraction(0)
        by_block = {block: Fraction(0) for block in EXPECTED_BY_BLOCK}
        for row in rows:
            key = (int(row["block"]), int(row["sum_index"]), output_index, int(row["branch"]))
            factor = factors[key]
            contribution = factor * Fraction(row["cell_complete_tail_majorant_upper"])
            tail += contribution
            roundoff += factor * Fraction(row["source_complex_multiply_add_cell_roundoff_upper"])
            by_block[int(row["block"])] += contribution
        prior = Fraction(prior_outputs[output_index]["all_recursive_cell_complete_majorant_upper"])
        complete = prior + tail + roundoff
        outputs.append(
            {
                "output_index": output_index,
                "output_label": prior_outputs[output_index]["output_label"],
                "prior_corrected_model_cell_complete_majorant_upper": decimal(prior),
                "all_later_cell_complete_tail_majorant_upper": decimal(tail),
                "all_later_cell_source_multiply_add_roundoff_upper": decimal(roundoff),
                "tail_by_block": {str(block): decimal(value) for block, value in by_block.items()},
                "later_complete_recurrence_partial_majorant_upper": decimal(complete),
                "margin_below_requested_scale": decimal(REQUESTED_SCALE - complete),
                "within_requested_scale": complete < REQUESTED_SCALE,
            }
        )
    return outputs


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy t=1e10 complete later-cell recurrence-tail transport

Date: 2026-08-09

Status: complete later-cell tail and source multiply-add rounding columns; not a proof of source q, the outer Hardy representation, or RH

The accepted nonzero coefficient cells for all `{aggregate['call_count']}`
recursive calls in blocks 21--28 preserve the exact factored cubic Legendre
remainder at 70 and 110 decimal digits.  Every saved point tail is contained
in its cell enclosure.  The common transformed-child phase is removed at unit
modulus before absolute values.

The pinned source environment uses radix-2, 113-bit, round-to-nearest
ties-to-even arithmetic.  Parent, reconstructed-child, anchored-multiplier,
child-state, and `tpm=-2*pi` rounding preimages are contained for every call,
and every emitted source complex multiply-add lies inside its operation bound.

```text
later cells                                      {aggregate['call_count']}
precision overlaps                               {aggregate['precision_overlap_count']} / {aggregate['call_count']}
maximum local complete tail                      {aggregate['maximum_local_cell_complete_tail_upper']}
maximum cell/point tail inflation                {aggregate['maximum_cell_to_point_tail_inflation_upper']}
maximum transported complete tail                {aggregate['maximum_transported_complete_tail_upper']}
maximum transported source multiply-add rounding {aggregate['maximum_transported_source_roundoff_upper']}
maximum partial complete output                  {aggregate['maximum_later_complete_recurrence_partial_majorant_upper']}
minimum remaining 0.005 margin                   {aggregate['minimum_margin_below_requested_scale']}
outputs below 0.005                              {aggregate['outputs_within_requested_scale']} / 15
```

No cancellation between recursive calls is used.

## Boundary

This closes the analytic complete-tail and source recurrence multiply-add
rounding columns on the finite saved later cells.  It does not yet enclose
source `q` or special-function variation over those cells.  The outer Hardy
representation/remainder, height uniformity, `Lambda<=0`, PF-infinity, RH,
and a prize-level conclusion remain open.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-limit-seconds", type=float, default=10800.0)
    args = parser.parse_args()
    for path in (FULL_CELLS, TELEMETRY, WEIGHTS, ROUNDING_MODE_PROBE, CHECKER):
        require(path.is_file(), f"missing later complete-recurrence dependency: {path}")
    started = time.monotonic()
    priority = set_low_priority()
    ctx.threads = 1
    rounding_mode = json.loads(ROUNDING_MODE_PROBE.read_text(encoding="utf-8"))
    require(bool(rounding_mode["passed"]) and rounding_mode["probe_values"]["rounding_mode"] == "nearest", "rounding mode dependency failed")
    full = json.loads(FULL_CELLS.read_text(encoding="utf-8"))
    require(bool(full["passed"]), "full later-cell dependency failed")
    source_rows = {int(row["chain"]): row for row in full["rows"]}
    require(len(source_rows) == EXPECTED_CALLS, "later source-cell roster drift")
    chains = later_cells.load_chains(set(source_rows))
    fingerprint = cache_fingerprint()
    cached = load_cache(fingerprint)
    require(set(cached).issubset(source_rows), "later complete-recurrence cache roster drift")
    pending = [chain for chain in sorted(source_rows) if chain not in cached]
    high_cpu_streak = 0
    for chain in pending:
        if time.monotonic() - started >= args.runtime_limit_seconds:
            print(f"parked later complete-recurrence cell transport after {len(cached)}/{EXPECTED_CALLS}")
            return 75
        row = build_row(chain, chains[chain], source_rows[chain])
        append_cache(row)
        cached[chain] = row
        if len(cached) % 25 == 0 or len(cached) == EXPECTED_CALLS:
            print(f"closed later complete-recurrence cell {len(cached)}/{EXPECTED_CALLS} (chain {chain})", flush=True)
            sampled = cpu_percent()
            if sampled is not None:
                high_cpu_streak = high_cpu_streak + 1 if sampled > CPU_PARK_THRESHOLD else 0
                if high_cpu_streak >= 2 and len(cached) < EXPECTED_CALLS:
                    print(f"parked later complete-recurrence cell transport after {len(cached)}/{EXPECTED_CALLS} (sustained CPU)")
                    return 75
    require(len(cached) == EXPECTED_CALLS, "later complete-recurrence cache incomplete")
    rows = [cached[chain] for chain in sorted(cached)]
    outputs = output_rows(rows, full)
    block_histogram = Counter(int(row["block"]) for row in rows)

    def maximum(field: str) -> tuple[Fraction, int]:
        return max((Fraction(row[field]), int(row["chain"])) for row in rows)

    def output_maximum(field: str) -> tuple[Fraction, int]:
        return max((Fraction(row[field]), int(row["output_index"])) for row in outputs)

    def output_minimum(field: str) -> tuple[Fraction, int]:
        return min((Fraction(row[field]), int(row["output_index"])) for row in outputs)

    max_local_tail = maximum("cell_complete_tail_majorant_upper")
    max_inflation = maximum("cell_to_point_tail_inflation_upper")
    max_rounding_ratio = maximum("parent_maximum_rounding_radius_to_cell_radius_ratio")
    max_tail_output = output_maximum("all_later_cell_complete_tail_majorant_upper")
    max_roundoff_output = output_maximum("all_later_cell_source_multiply_add_roundoff_upper")
    max_complete = output_maximum("later_complete_recurrence_partial_majorant_upper")
    min_margin = output_minimum("margin_below_requested_scale")
    aggregate = {
        "call_count": len(rows),
        "call_count_by_block": {str(block): block_histogram[block] for block in EXPECTED_BY_BLOCK},
        "precision_overlap_count": sum(bool(row["precision_overlap"]) for row in rows),
        "point_tail_containment_count": sum(bool(row["point_tail_contained_in_cell_tail"]) for row in rows),
        "parent_rounding_preimage_containment_count": sum(bool(row["parent_binary128_rounding_preimages_inside_cell"]) for row in rows),
        "child_rounding_preimage_containment_count": sum(bool(row["child_binary128_rounding_preimages_inside_reconstructed_cell"]) for row in rows),
        "multiplier_rounding_preimage_containment_count": sum(bool(row["multiplier_binary128_rounding_preimages_inside_anchored_cell"]) for row in rows),
        "child_state_rounding_preimage_containment_count": sum(bool(row["child_state_binary128_rounding_preimages_inside_cell_kernel"]) for row in rows),
        "tpm_rounding_containment_count": sum(bool(row["tpm_rounding_preimage_contains_minus_two_pi"]) for row in rows),
        "source_multiply_add_roundoff_containment_count": sum(
            Fraction(row["source_complex_multiply_add_observed_gap"])
            <= Fraction(row["source_complex_multiply_add_point_roundoff_upper"])
            <= Fraction(row["source_complex_multiply_add_cell_roundoff_upper"])
            for row in rows
        ),
        "maximum_local_cell_complete_tail_upper": decimal(max_local_tail[0]),
        "maximum_local_cell_complete_tail_witness": max_local_tail[1],
        "maximum_cell_to_point_tail_inflation_upper": decimal(max_inflation[0]),
        "maximum_cell_to_point_tail_inflation_witness": max_inflation[1],
        "maximum_parent_rounding_radius_to_cell_radius_ratio": decimal(max_rounding_ratio[0]),
        "maximum_parent_rounding_radius_to_cell_radius_witness": max_rounding_ratio[1],
        "maximum_transported_complete_tail_upper": decimal(max_tail_output[0]),
        "maximum_transported_complete_tail_witness_output": max_tail_output[1],
        "maximum_transported_source_roundoff_upper": decimal(max_roundoff_output[0]),
        "maximum_transported_source_roundoff_witness_output": max_roundoff_output[1],
        "maximum_later_complete_recurrence_partial_majorant_upper": decimal(max_complete[0]),
        "maximum_later_complete_recurrence_partial_witness_output": max_complete[1],
        "minimum_margin_below_requested_scale": decimal(min_margin[0]),
        "minimum_margin_witness_output": min_margin[1],
        "outputs_within_requested_scale": sum(bool(row["within_requested_scale"]) for row in outputs),
        "uses_cross_call_cancellation": False,
        "remaining_later_point_only_tail_calls": 0,
    }
    complete_counts = all(
        aggregate[key] == EXPECTED_CALLS
        for key in (
            "precision_overlap_count",
            "point_tail_containment_count",
            "parent_rounding_preimage_containment_count",
            "child_rounding_preimage_containment_count",
            "multiplier_rounding_preimage_containment_count",
            "child_state_rounding_preimage_containment_count",
            "tpm_rounding_containment_count",
            "source_multiply_add_roundoff_containment_count",
        )
    )
    passed = (
        aggregate["call_count_by_block"] == {str(block): count for block, count in EXPECTED_BY_BLOCK.items()}
        and complete_counts
        and aggregate["outputs_within_requested_scale"] == 15
    )
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_cell_full_gate",
        "status": "all_1040_later_cells_close_complete_tail_and_pinned_source_multiply_add_rounding" if passed else "later_complete_recurrence_cell_full_gate_falsified",
        "scope": "All recursive calls in blocks 21 through 28 and all 15 outer outputs at the saved t=1e10 fixture.",
        "passed": passed,
        "precisions_decimal_digits": list(PRECISIONS),
        "theorem": pilot.build()["theorem"],
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
            "full_later_cells": {"path": relative(FULL_CELLS), "sha256": file_hash(FULL_CELLS)},
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "weights": {"path": relative(WEIGHTS), "sha256": file_hash(WEIGHTS)},
            "rounding_mode_probe": {"path": relative(ROUNDING_MODE_PROBE), "sha256": file_hash(ROUNDING_MODE_PROBE)},
            "stress_pilot": {"path": relative(pilot.RESULT), "sha256": file_hash(pilot.RESULT)},
            "cell_tail_builder": {"path": relative(Path(cell_tail.__file__).resolve()), "sha256": file_hash(Path(cell_tail.__file__).resolve())},
            "point_tail_builder": {"path": relative(Path(point_tail.__file__).resolve()), "sha256": file_hash(Path(point_tail.__file__).resolve())},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": "Enclose the source q and special-function variation on the accepted cells, then certify the outer Hardy representation/remainder and replace the finite saved-height atlas by a uniform theorem.",
        "proof_boundary": "Complete finite saved-height later-cell analytic tail and source multiply-add rounding columns only. Source q/special-function cell variation, the outer Hardy remainder, height uniformity, Lambda<=0, PF-infinity, RH, and a prize-level conclusion remain open.",
    }
    require(passed, "complete later recurrence-tail cell transport did not close")
    tmp = RESULT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(RESULT)
    note_tmp = NOTE.with_suffix(".md.tmp")
    note_tmp.write_text(render_note(artifact), encoding="utf-8")
    note_tmp.replace(NOTE)
    print(
        "built complete later recurrence-tail cell transport: "
        f"{aggregate['call_count']}/{EXPECTED_CALLS} cells, "
        f"{aggregate['outputs_within_requested_scale']}/15 outputs, "
        f"max {aggregate['maximum_later_complete_recurrence_partial_majorant_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
