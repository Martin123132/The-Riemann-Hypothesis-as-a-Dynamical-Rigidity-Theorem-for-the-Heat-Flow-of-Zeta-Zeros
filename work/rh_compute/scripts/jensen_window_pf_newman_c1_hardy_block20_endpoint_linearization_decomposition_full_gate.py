#!/usr/bin/env python3
"""Resumable all-call endpoint-linearization decomposition certificate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from time import monotonic
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_pilot_gate as pilot
import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate as log_full
import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate as log_pilot
import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells


CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_full_gate_cache.jsonl"
)
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_full_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_full_gate.md"
PILOT_RESULT = pilot.RESULT
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_full_gate.py"
)
PRECISIONS = pilot.PRECISIONS


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def upper_abs(value: Any) -> Fraction:
    return log_full.upper(abs(value))


def lower_abs(value: Any) -> Fraction:
    return log_full.lower(abs(value))


def load_log_full_cache() -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {}
    with log_full.CACHE.open("r", encoding="utf-8") as handle:
        header = json.loads(next(handle))
        require(header["kind"].endswith("logarithmic_endpoint_ray_arb_full_gate_cache"), "log-ray cache header drift")
        for line in handle:
            row = json.loads(line)
            rows[int(row["chain"])] = row
    return rows


def cache_fingerprint() -> str:
    payload = {
        "dependencies": {
            "log_full_result": file_hash(log_full.RESULT),
            "log_full_cache": file_hash(log_full.CACHE),
            "decomposition_pilot": file_hash(PILOT_RESULT),
            "decomposition_pilot_builder": file_hash(pilot.BUILDER),
            "builder": file_hash(BUILDER),
            "checker": file_hash(CHECKER),
        },
        "precisions": list(PRECISIONS),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_cache(fingerprint: str) -> dict[int, dict[str, Any]]:
    if not CACHE.exists():
        CACHE.write_text(
            json.dumps({"kind": "endpoint_linearization_decomposition_full_gate_cache", "fingerprint": fingerprint}, sort_keys=True)
            + "\n",
            encoding="utf-8",
        )
        return {}
    with CACHE.open("r", encoding="utf-8") as handle:
        lines = [line for line in handle if line.strip()]
    require(lines, "endpoint-linearization full cache empty")
    require(json.loads(lines[0]).get("fingerprint") == fingerprint, "endpoint-linearization cache fingerprint drift")
    rows: dict[int, dict[str, Any]] = {}
    for line in lines[1:]:
        row = json.loads(line)
        chain = int(row["chain"])
        require(chain not in rows, f"duplicate endpoint-linearization chain {chain}")
        rows[chain] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    with CACHE.open("a", encoding="utf-8", buffering=1) as handle:
        handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def evaluate_row(
    chain: int,
    data: dict[str, Any],
    sector_row: dict[str, Any],
    eq69_row: dict[str, Any],
    w2_row: dict[str, Any],
    cutoff: Fraction,
) -> dict[str, Any]:
    low = pilot.decompose(data, sector_row, eq69_row, w2_row, PRECISIONS[0], cutoff)
    high = pilot.decompose(data, sector_row, eq69_row, w2_row, PRECISIONS[1], cutoff)
    for name in high:
        require(low[name].overlaps(high[name]), f"chain {chain} full decomposition {name} precision nonoverlap")
    require(high["generic_w5_gap"].contains(0), f"chain {chain} full generic/W5 identity excludes zero")
    require(
        high["correction_decomposition_gap"].contains(0),
        f"chain {chain} full correction decomposition excludes zero",
    )
    correction_lower = lower_abs(high["exact_correction"])
    correction_upper = upper_abs(high["exact_correction"])
    require(correction_lower > 0, f"chain {chain} exact correction contains zero")
    component_sum = sum(
        upper_abs(high[name])
        for name in ("b_remainder", "c_remainder", "zero_mode", "paper_w2", "paper_w34")
    )
    parent = data["levels"][1]
    phi1 = cells.binary128_fraction(parent["coefficients_hex"][0])
    phi3 = cells.binary128_fraction(parent["coefficients_hex"][2])
    return {
        "chain": chain,
        "sum_index": int(data["header"]["sum_index"]),
        "branch": int(data["header"]["branch"]),
        "phi1_sign": "positive" if phi1 > 0 else "negative",
        "phi3_sign": "positive" if phi3 > 0 else "negative",
        "ray_cutoff": int(cutoff),
        "b_remainder": point_balls.acb_record(high["b_remainder"], 45),
        "c_remainder": point_balls.acb_record(high["c_remainder"], 45),
        "zero_mode": point_balls.acb_record(high["zero_mode"], 45),
        "paper_w2": point_balls.acb_record(high["paper_w2"], 45),
        "paper_w3_or_w4": point_balls.acb_record(high["paper_w34"], 45),
        "correction_prediction": point_balls.acb_record(high["correction_prediction"], 45),
        "exact_correction": point_balls.acb_record(high["exact_correction"], 45),
        "generic_w5_gap_upper": decimal(upper_abs(high["generic_w5_gap"])),
        "correction_decomposition_gap_upper": decimal(upper_abs(high["correction_decomposition_gap"])),
        "exact_correction_magnitude_lower": decimal(correction_lower),
        "exact_correction_magnitude_upper": decimal(correction_upper),
        "component_triangle_sum_upper": decimal(component_sum),
        "component_triangle_to_correction_lower_ratio_upper": decimal(component_sum / correction_lower),
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 endpoint-linearization decomposition full gate

Date: 2026-08-07

Status: rigorous all-374 finite-input decomposition; not a proof of a height-uniform remainder theorem or RH

The endpoint-linearization pilot is replayed at 70 and 110 decimal digits on
all 374 saved recursive calls.  Each chain reuses the certified full-ray
cutoff and is fsynced atomically to an append-only cache.

On every call the exact identities

```text
endpoint_half+a_endpoint+B_linear+C_linear = paper_half+W5,
Q_exact-Q_paper = R_b+R_c+zero_mode-W2-(W3 or W4)
```

enclose zero.  Aggregate bounds are

```text
maximum generic/W5 identity gap       <= {aggregate['maximum_generic_w5_gap_upper']}
maximum correction decomposition gap  <= {aggregate['maximum_correction_decomposition_gap_upper']}
maximum component triangle sum        <= {aggregate['maximum_component_triangle_sum_upper']}
maximum triangle/correction ratio      <= {aggregate['maximum_component_triangle_to_correction_ratio_upper']}
```

The worst triangle ratio is witnessed by chain
`{aggregate['maximum_component_triangle_to_correction_ratio_witness']}`.  It
measures cancellation already lost by a fully componentwise absolute bound;
it is a diagnostic for designing the next analytic majorant, not a fitted
constant and not a recurrence budget.

This gate rigorously localizes the complete finite correction into generic
endpoint-linearization remainders and the exceptional W2--W4 replacements.
It supplies neither uniform-height piece bounds nor recurrence propagation or
outer Hardy control.  It is not a proof of `Lambda<=0`, PF-infinity, RH, or a
prize-level conclusion.
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runtime-limit-seconds", type=float, default=10800.0)
    args = parser.parse_args()
    log_pilot.set_low_priority()
    started = monotonic()

    sector = json.loads(log_pilot.SECTOR.read_text(encoding="utf-8"))
    eq69 = json.loads(log_pilot.EQ69.read_text(encoding="utf-8"))
    w2 = json.loads(log_pilot.W2_W5.read_text(encoding="utf-8"))
    sector_rows = {int(row["chain"]): row for row in sector["rows"]}
    eq69_rows = {int(row["chain"]): row for row in eq69["rows"]}
    w2_rows = {int(row["chain"]): row for row in w2["rows"]}
    recursive = native_q.load_recursive_chains()
    log_rows = load_log_full_cache()
    require(set(recursive) == set(sector_rows) == set(eq69_rows) == set(w2_rows) == set(log_rows), "full decomposition roster drift")

    fingerprint = cache_fingerprint()
    cached = load_cache(fingerprint)
    for chain in sorted(recursive):
        if chain in cached:
            continue
        cutoff = Fraction(int(log_rows[chain]["ray_cutoff"]))
        row = evaluate_row(chain, recursive[chain], sector_rows[chain], eq69_rows[chain], w2_rows[chain], cutoff)
        append_cache(row)
        cached[chain] = row
        if monotonic() - started >= args.runtime_limit_seconds:
            print(f"parked endpoint-linearization decomposition after {len(cached)}/374 calls")
            return 2

    require(len(cached) == 374, "full endpoint-linearization cache incomplete")
    rows = [cached[chain] for chain in sorted(cached)]

    def maximum(name: str) -> tuple[Fraction, int]:
        return max((Fraction(row[name]), int(row["chain"])) for row in rows)

    max_generic, generic_witness = maximum("generic_w5_gap_upper")
    max_gap, gap_witness = maximum("correction_decomposition_gap_upper")
    max_triangle, triangle_witness = maximum("component_triangle_sum_upper")
    max_ratio, ratio_witness = maximum("component_triangle_to_correction_lower_ratio_upper")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_full_gate",
        "status": "rigorous_all_374_endpoint_linearization_decomposition_validated",
        "precisions_decimal_digits": list(PRECISIONS),
        "rows": rows,
        "aggregate": {
            "call_count": len(rows),
            "generic_w5_identity_count": sum(float(row["generic_w5_gap_upper"]) < 1e-40 for row in rows),
            "correction_decomposition_identity_count": sum(
                float(row["correction_decomposition_gap_upper"]) < 1e-10 for row in rows
            ),
            "maximum_generic_w5_gap_upper": decimal(max_generic),
            "maximum_generic_w5_gap_witness": generic_witness,
            "maximum_correction_decomposition_gap_upper": decimal(max_gap),
            "maximum_correction_decomposition_gap_witness": gap_witness,
            "maximum_component_triangle_sum_upper": decimal(max_triangle),
            "maximum_component_triangle_sum_witness": triangle_witness,
            "maximum_component_triangle_to_correction_ratio_upper": decimal(max_ratio),
            "maximum_component_triangle_to_correction_ratio_witness": ratio_witness,
        },
        "dependencies": {
            "log_full_result": {"path": relative(log_full.RESULT), "sha256": file_hash(log_full.RESULT)},
            "log_full_cache": {"path": relative(log_full.CACHE), "sha256": file_hash(log_full.CACHE)},
            "decomposition_pilot": {"path": relative(PILOT_RESULT), "sha256": file_hash(PILOT_RESULT)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "Rigorous all-374 finite-input decomposition only. No height-uniform component majorant, recurrence "
            "accumulation, outer Hardy control, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print("built endpoint-linearization decomposition full gate: 374/374 calls, both exact identities enclose zero")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
