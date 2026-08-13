#!/usr/bin/env python3
"""Resumable all-call exceptional-mode recombination certificate."""

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

import jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_pilot_gate as pilot
import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate as rays
import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells


CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_full_gate_cache.jsonl"
)
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_full_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_full_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_full_gate.py"
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


def cache_fingerprint() -> str:
    payload = {
        "dependencies": {
            "source": file_hash(pilot.SOURCE),
            "sector": file_hash(rays.SECTOR),
            "pilot_result": file_hash(pilot.RESULT),
            "pilot_builder": file_hash(pilot.BUILDER),
            "builder": file_hash(BUILDER),
            "checker": file_hash(CHECKER),
        },
        "precisions": list(PRECISIONS),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_cache(fingerprint: str) -> dict[int, dict[str, Any]]:
    if not CACHE.exists():
        CACHE.write_text(
            json.dumps({"kind": "exceptional_mode_recombination_full_gate_cache", "fingerprint": fingerprint}, sort_keys=True)
            + "\n",
            encoding="utf-8",
        )
        return {}
    with CACHE.open("r", encoding="utf-8") as handle:
        lines = [line for line in handle if line.strip()]
    require(lines, "exceptional-mode full cache empty")
    require(json.loads(lines[0]).get("fingerprint") == fingerprint, "exceptional-mode cache fingerprint drift")
    rows: dict[int, dict[str, Any]] = {}
    for line in lines[1:]:
        row = json.loads(line)
        chain = int(row["chain"])
        require(chain not in rows, f"duplicate exceptional-mode chain {chain}")
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
    source_row: dict[str, Any],
) -> dict[str, Any]:
    low = pilot.recombine(data, sector_row, source_row, PRECISIONS[0])
    high = pilot.recombine(data, sector_row, source_row, PRECISIONS[1])
    for name in (
        "b_exact_mode",
        "b_linear_mode",
        "b_quadratic_mode",
        "w2_model_gap",
        "upper_generic",
        "upper_exceptional",
        "upper_recombined",
        "upper_recombination_gap",
        "w34_model_gap",
        "lower_recombined",
        "lower_recombination_gap",
        "total_recombined",
        "total_gap",
    ):
        require(low[name].overlaps(high[name]), f"chain {chain} full exceptional-mode {name} precision nonoverlap")
    for name in ("w2_model_gap", "w34_model_gap", "upper_recombination_gap", "lower_recombination_gap", "total_gap"):
        require(high[name].contains(0), f"chain {chain} full exceptional-mode {name} excludes zero")
    correction_lower = pilot.lower_abs(high["correction"])
    require(correction_lower > 0, f"chain {chain} exact correction contains zero")
    model_gap = max(pilot.upper_abs(high["w2_model_gap"]), pilot.upper_abs(high["w34_model_gap"]))
    parent = data["levels"][1]
    phi1 = cells.binary128_fraction(parent["coefficients_hex"][0])
    phi3 = cells.binary128_fraction(parent["coefficients_hex"][2])
    return {
        "chain": chain,
        "sum_index": int(data["header"]["sum_index"]),
        "branch": int(data["header"]["branch"]),
        "phi1_sign": "positive" if phi1 > 0 else "negative",
        "phi3_sign": "positive" if phi3 > 0 else "negative",
        "selector": high["selector"],
        "ray_cutoff": int(source_row["ray_cutoff"]),
        "b_exact_mode": point_balls.acb_record(high["b_exact_mode"], 45),
        "b_exact_minus_quadratic": point_balls.acb_record(high["upper_exceptional"], 45),
        "upper_recombined": point_balls.acb_record(high["upper_recombined"], 45),
        "lower_recombined": point_balls.acb_record(high["lower_recombined"], 45),
        "maximum_model_identity_gap_upper": decimal(model_gap),
        "upper_recombination_gap_upper": decimal(pilot.upper_abs(high["upper_recombination_gap"])),
        "lower_recombination_gap_upper": decimal(pilot.upper_abs(high["lower_recombination_gap"])),
        "total_recombination_gap_upper": decimal(pilot.upper_abs(high["total_gap"])),
        "exceptional_mode_tail_upper": decimal(max(high["b_tail_upper"], high["c_tail_upper"])),
        "component_triangle_sum_upper": decimal(high["component_triangle"]),
        "component_triangle_to_correction_ratio_upper": decimal(high["component_triangle"] / correction_lower),
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 exceptional-mode recombination full gate

Date: 2026-08-07

Status: rigorous all-374 finite-input identity and quadrature certificate; not a height-uniform bound and not a proof of RH

The four-branch exceptional-mode construction was replayed serially at 70 and
110 decimal digits on all 374 saved calls.  Every completed call was fsynced
to an append-only cache.  The independent paper-model identities

```text
W2 = M_b,lin(N)-M_b,quad(N),
W3 = conjugate(M_c,lin(0)-M_c,quad(0)),
W4 = Z_quad(0)+i*exp(-2*pi*i*F(N))/(2*pi*F'(N))
```

and the resulting singularity-free recombinations all enclose equality.  The
quadratic modes use the paper's global `Phi2`, while the exact modes retain
the complete cubic phase.  Thus no `1/delta` or `1/eta` term remains in the
objects that now require uniform bounds.

```text
W3 calls                              = {aggregate['w3_call_count']}
W4 calls                              = {aggregate['w4_call_count']}
maximum paper-model identity gap      <= {aggregate['maximum_model_identity_gap_upper']}
maximum complete recombination gap    <= {aggregate['maximum_total_recombination_gap_upper']}
maximum exceptional-mode tail         <= {aggregate['maximum_exceptional_mode_tail_upper']}
maximum four-piece triangle/correction ratio <= {aggregate['maximum_component_triangle_to_correction_ratio_upper']}
```

The worst triangle ratio occurs at chain
`{aggregate['maximum_component_triangle_to_correction_ratio_witness']}`.  It
shows that removing the formal poles does not remove the need to preserve
coupling between the upper and lower endpoint families.  The ratio is a
finite-roster diagnostic, not a fitted theorem constant.

This gate does not provide a height-uniform finite-difference majorant,
recursive accumulation, outer Hardy control, `Lambda<=0`, PF-infinity, RH, or
a prize-level conclusion.
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runtime-limit-seconds", type=float, default=10800.0)
    args = parser.parse_args()
    rays.set_low_priority()
    started = monotonic()

    recursive = native_q.load_recursive_chains()
    sector = json.loads(rays.SECTOR.read_text(encoding="utf-8"))
    sector_rows = {int(row["chain"]): row for row in sector["rows"]}
    source = json.loads(pilot.SOURCE.read_text(encoding="utf-8"))
    source_rows = {int(row["chain"]): row for row in source["rows"]}
    require(set(recursive) == set(sector_rows) == set(source_rows), "exceptional-mode full roster drift")

    fingerprint = cache_fingerprint()
    cached = load_cache(fingerprint)
    for chain in sorted(recursive):
        if chain in cached:
            continue
        row = evaluate_row(chain, recursive[chain], sector_rows[chain], source_rows[chain])
        append_cache(row)
        cached[chain] = row
        if monotonic() - started >= args.runtime_limit_seconds:
            print(f"parked exceptional-mode recombination after {len(cached)}/374 calls")
            return 2

    require(len(cached) == 374, "exceptional-mode full cache incomplete")
    rows = [cached[chain] for chain in sorted(cached)]

    def maximum(name: str) -> tuple[Fraction, int]:
        return max((Fraction(row[name]), int(row["chain"])) for row in rows)

    max_model, model_witness = maximum("maximum_model_identity_gap_upper")
    max_gap, gap_witness = maximum("total_recombination_gap_upper")
    max_tail, tail_witness = maximum("exceptional_mode_tail_upper")
    max_ratio, ratio_witness = maximum("component_triangle_to_correction_ratio_upper")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_full_gate",
        "status": "rigorous_all_374_exceptional_mode_recombination_validated",
        "precisions_decimal_digits": list(PRECISIONS),
        "rows": rows,
        "aggregate": {
            "call_count": len(rows),
            "w3_call_count": sum(row["selector"] == "W3" for row in rows),
            "w4_call_count": sum(row["selector"] == "W4" for row in rows),
            "model_identity_count": sum(float(row["maximum_model_identity_gap_upper"]) < 1e-40 for row in rows),
            "recombination_identity_count": sum(float(row["total_recombination_gap_upper"]) < 1e-10 for row in rows),
            "maximum_model_identity_gap_upper": decimal(max_model),
            "maximum_model_identity_gap_witness": model_witness,
            "maximum_total_recombination_gap_upper": decimal(max_gap),
            "maximum_total_recombination_gap_witness": gap_witness,
            "maximum_exceptional_mode_tail_upper": decimal(max_tail),
            "maximum_exceptional_mode_tail_witness": tail_witness,
            "maximum_component_triangle_to_correction_ratio_upper": decimal(max_ratio),
            "maximum_component_triangle_to_correction_ratio_witness": ratio_witness,
        },
        "dependencies": {
            "source": {"path": relative(pilot.SOURCE), "sha256": file_hash(pilot.SOURCE)},
            "sector": {"path": relative(rays.SECTOR), "sha256": file_hash(rays.SECTOR)},
            "pilot": {"path": relative(pilot.RESULT), "sha256": file_hash(pilot.RESULT)},
            "cache": {"path": relative(CACHE), "sha256": file_hash(CACHE)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "Rigorous all-374 finite-input identity and mode-quadrature certificate only. No height-uniform "
            "finite-difference majorant, recursive accumulation, outer Hardy control, Lambda<=0, PF-infinity, "
            "RH, or prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print("built exceptional-mode recombination full gate: 374/374 calls, reciprocal poles absent from bounded objects")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
