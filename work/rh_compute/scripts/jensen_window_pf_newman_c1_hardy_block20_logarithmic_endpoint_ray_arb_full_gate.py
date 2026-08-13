#!/usr/bin/env python3
"""Rigorous resumable all-call Arb evaluation of the logarithmic endpoint rays."""

from __future__ import annotations

import argparse
from collections import Counter
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

import jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate as pilot
import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import arb, ctx


SECTOR = pilot.SECTOR
LOG_REDUCTION = pilot.LOG_REDUCTION
EQ69 = pilot.EQ69
W2_W5 = pilot.W2_W5
PILOT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate.json"
)
CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate_cache.jsonl"
)
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate.py"
)
PRECISIONS = pilot.PRECISIONS
TAIL_TARGET = Fraction(1, 10**30)
INITIAL_CUTOFF = Fraction(64)
MAXIMUM_CUTOFF = Fraction(65536)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def upper(value: arb) -> Fraction:
    return cells.bound_fraction(value.upper())


def lower(value: arb) -> Fraction:
    return cells.bound_fraction(value.lower())


def upper_abs(value: Any) -> Fraction:
    return upper(abs(value))


def tail_bound(
    data: dict[str, Any],
    sector_row: dict[str, Any],
    endpoint: int,
    family: str,
    cutoff: Fraction,
    dps: int = 70,
) -> Fraction:
    ctx.dps = dps
    ctx.threads = 1
    parent = data["levels"][1]
    phi1, phi2, phi3 = (pilot.exact_ball(value) for value in parent["coefficients_hex"])
    contract = sector_row[f"{family}_contract"]
    ray_direction, sine = pilot.direction(contract["angle"])
    start = int(data["levels"][2]["length"]) + 1 if family == "b" else 1
    endpoint_ball = arb(endpoint)
    endpoint_derivative = phi1 + 2 * phi2 * endpoint_ball + 3 * phi3 * endpoint_ball**2
    endpoint_second = 2 * phi2 + 6 * phi3 * endpoint_ball
    linear_decay = cells.fraction_ball(pilot.fraction_from_record(contract["linear_decay"]))
    quadratic_decay = cells.fraction_ball(pilot.dyadic(contract["quadratic_decay"]["lower_dyadic"]))
    cubic_decay = cells.fraction_ball(pilot.fraction_from_record(contract["cubic_decay"]))
    far = cells.fraction_ball(cutoff)
    pi = arb.pi()
    angular_decay = 2 * pi * cells.fraction_ball(sine)
    decay_rate = 2 * pi * (
        linear_decay + quadratic_decay * far + cubic_decay * far**2
    )
    inverse = 1 / decay_rate
    p0 = abs(endpoint_derivative)
    p1 = abs(endpoint_second)
    p2 = 3 * abs(phi3)
    polynomial_tail = (-decay_rate * far).exp() * (
        p0 * inverse
        + p1 * (far * inverse + inverse**2)
        + p2 * (far**2 * inverse + 2 * far * inverse**2 + 2 * inverse**3)
    )
    denominator = start * (1 - (-angular_decay * far).exp())
    return upper(abs(ray_direction) * polynomial_tail / denominator)


def select_cutoff(
    data: dict[str, Any], sector_row: dict[str, Any]
) -> tuple[Fraction, Fraction]:
    cutoff = INITIAL_CUTOFF
    while True:
        worst = max(
            tail_bound(data, sector_row, endpoint, family, cutoff)
            for family in ("b", "c")
            for endpoint in (0, 104)
        )
        if worst < TAIL_TARGET:
            return cutoff, worst
        cutoff *= 2
        require(cutoff <= MAXIMUM_CUTOFF, "log-ray full cutoff exceeded safety ceiling")


def breaks_for_cutoff(cutoff: Fraction) -> tuple[Fraction, ...]:
    breaks = list(pilot.CENTRAL_BREAKS)
    while breaks[-1] < cutoff:
        breaks.append(breaks[-1] * 2)
    require(breaks[-1] == cutoff, "log-ray full cutoff not dyadic from pilot base")
    return tuple(breaks)


def cache_fingerprint() -> str:
    payload = {
        "dependencies": {
            "sector": file_hash(SECTOR),
            "log_reduction": file_hash(LOG_REDUCTION),
            "equation_69": file_hash(EQ69),
            "w2_w5": file_hash(W2_W5),
            "pilot_result": file_hash(PILOT_RESULT),
            "pilot_builder": file_hash(pilot.BUILDER),
            "builder": file_hash(BUILDER),
            "checker": file_hash(CHECKER),
        },
        "precisions": list(PRECISIONS),
        "tail_target": [TAIL_TARGET.numerator, TAIL_TARGET.denominator],
        "initial_cutoff": [INITIAL_CUTOFF.numerator, INITIAL_CUTOFF.denominator],
        "origin_cutoff": [pilot.ORIGIN_CUTOFF.numerator, pilot.ORIGIN_CUTOFF.denominator],
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def load_cache(fingerprint: str) -> dict[int, dict[str, Any]]:
    if not CACHE.exists():
        header = {
            "kind": "jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate_cache",
            "fingerprint": fingerprint,
        }
        CACHE.write_text(json.dumps(header, sort_keys=True) + "\n", encoding="utf-8")
        return {}
    rows: dict[int, dict[str, Any]] = {}
    with CACHE.open("r", encoding="utf-8") as handle:
        lines = [line for line in handle if line.strip()]
    require(lines, "log-ray full cache is empty")
    header = json.loads(lines[0])
    require(header.get("fingerprint") == fingerprint, "log-ray full cache fingerprint drift")
    for line in lines[1:]:
        row = json.loads(line)
        chain = int(row["chain"])
        require(chain not in rows, f"duplicate log-ray full cache chain: {chain}")
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
) -> dict[str, Any]:
    cutoff, selected_tail = select_cutoff(data, sector_row)
    central_breaks = breaks_for_cutoff(cutoff)
    low = pilot.evaluate(
        data,
        sector_row,
        eq69_row,
        w2_row,
        PRECISIONS[0],
        ray_cutoff=cutoff,
        central_breaks=central_breaks,
    )
    high = pilot.evaluate(
        data,
        sector_row,
        eq69_row,
        w2_row,
        PRECISIONS[1],
        ray_cutoff=cutoff,
        central_breaks=central_breaks,
    )
    for index, (low_ray, high_ray) in enumerate(zip(low["rays"], high["rays"])):
        require(low_ray["total"].overlaps(high_ray["total"]), f"chain {chain} full ray {index} precision nonoverlap")
    for name in (
        "b_source",
        "c_source",
        "endpoint_half",
        "a_endpoint",
        "zero_mode",
        "exact_nonsaddle",
        "independent_target",
        "formula_target_gap",
        "exact_minus_paper",
        "correction_residual_gap",
    ):
        require(low[name].overlaps(high[name]), f"chain {chain} full {name} precision nonoverlap")
    require(high["formula_target_gap"].contains(0), f"chain {chain} full formula/target identity excludes zero")
    require(high["correction_residual_gap"].contains(0), f"chain {chain} full correction/residual identity excludes zero")
    require(high["maximum_tail_upper"] < TAIL_TARGET, f"chain {chain} full evaluated tail target drift")
    exact_correction_abs = abs(high["exact_minus_paper"])
    correction_lower = lower(exact_correction_abs)
    correction_upper = upper(exact_correction_abs)
    require(correction_lower > 0, f"chain {chain} exact nonsaddle correction contains zero")
    parent = data["levels"][1]
    phi1 = cells.binary128_fraction(parent["coefficients_hex"][0])
    phi3 = cells.binary128_fraction(parent["coefficients_hex"][2])
    return {
        "chain": chain,
        "sum_index": int(data["header"]["sum_index"]),
        "branch": int(data["header"]["branch"]),
        "phi1_sign": "positive" if phi1 > 0 else "negative",
        "phi3_sign": "positive" if phi3 > 0 else "negative",
        "b_angle": sector_row["b_contract"]["angle"],
        "c_angle": sector_row["c_contract"]["angle"],
        "ray_cutoff": int(cutoff),
        "selected_tail_bound_upper": decimal(selected_tail),
        "origin_radius_upper": decimal(high["maximum_origin_upper"]),
        "tail_radius_upper": decimal(high["maximum_tail_upper"]),
        "precision_overlap": True,
        "exact_nonsaddle": point_balls.acb_record(high["exact_nonsaddle"], 50),
        "independent_target": point_balls.acb_record(high["independent_target"], 50),
        "formula_target_gap": point_balls.acb_record(high["formula_target_gap"], 30),
        "formula_target_gap_magnitude_upper": decimal(upper_abs(high["formula_target_gap"])),
        "paper_nonsaddle": point_balls.acb_record(high["paper_nonsaddle"], 50),
        "exact_minus_paper": point_balls.acb_record(high["exact_minus_paper"], 50),
        "exact_minus_paper_magnitude_lower": decimal(correction_lower),
        "exact_minus_paper_magnitude_upper": decimal(correction_upper),
        "exact_minus_paper_excludes_zero": True,
        "paper_residual": point_balls.acb_record(high["paper_residual"], 50),
        "correction_residual_gap": point_balls.acb_record(high["correction_residual_gap"], 30),
        "correction_residual_gap_magnitude_upper": decimal(upper_abs(high["correction_residual_gap"])),
        "formula_target_identity_contains_zero": True,
        "correction_residual_identity_contains_zero": True,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Hardy block-20 logarithmic endpoint-ray full Arb gate

Date: 2026-08-07

Status: rigorous all-374 finite-input logarithmic endpoint-ray enclosure; not a proof of a height-uniform recurrence or RH

## Scope

The four-branch pilot is extended to every recursive block-20 call.  Each
chain is evaluated atomically at 70 and 110 decimal digits and fsynced to an
append-only JSONL cache before the next chain begins.  The run remains serial
with one below-normal RH worker.

For each call the ray cutoff is the smallest power-of-two multiple of 64 for
which the exact polynomial-times-exponential tail bound is below `1e-30`:

```text
{json.dumps(aggregate['cutoff_distribution'], sort_keys=True)}
maximum selected cutoff: {aggregate['maximum_ray_cutoff']}
```

No numerical tail at the cutoff is assumed to be zero.

## Exact identities

All {aggregate['call_count']} calls have 70/110-digit overlap for the four
rays and complete nonsaddle value.  On every call both rigorous complex
identities contain zero:

```text
Q_exact^- - (P^- - I69^-) = 0,
(Q_exact^- - Q_paper^-) - R_paper = 0.                 (1)
```

Aggregate enclosure bounds are

```text
maximum omitted-origin radius       <= {aggregate['maximum_origin_radius_upper']}
maximum far-tail radius             <= {aggregate['maximum_tail_radius_upper']}
maximum |formula-target gap|        <= {aggregate['maximum_formula_target_gap_magnitude_upper']}
maximum |correction-residual gap|    <= {aggregate['maximum_correction_residual_gap_magnitude_upper']}
minimum |Q_exact-Q_paper|            >= {aggregate['minimum_exact_correction_magnitude_lower']}
maximum |Q_exact-Q_paper|            <= {aggregate['maximum_exact_correction_magnitude_upper']}
nonzero exact corrections             {aggregate['exact_correction_excluding_zero_count']} / {aggregate['call_count']}
```

Thus the residual left after exact equation-(69) integration and published
W2--W5 evaluation is independently reconstructed as the exact full-cubic
nonsaddle correction on all 374 calls.  It is not a fitted term and is not a
source patch.

## Interpretation

This closes the finite block-20 contour-integration obligation: the exact
nonsaddle value, the direct parent-minus-saddle value, and the previously
recorded published residual agree as rigorous complex balls.  The remaining
proof problem is no longer numerical identification of this finite residual.
It is to derive a height-uniform analytic bound for the full-cubic versus
quadratic endpoint model and propagate that bound through the complete
recurrence and outer Hardy representation.

## Pi provenance and boundary

Pi is inherited from the paper's Fourier exponential and is evaluated as the
mathematical Arb constant.  It controls the exact logarithmic kernel, phase
and tail decay; no fitted normalization is inserted.

This is rigorous finite-input work for one low-height block.  It proves no
height-uniform W2--W4 estimate, recursive/global error accumulation, outer
Hardy bound, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runtime-limit-seconds", type=float, default=10800.0)
    args = parser.parse_args()
    for path in (SECTOR, LOG_REDUCTION, EQ69, W2_W5, PILOT_RESULT, CHECKER):
        require(path.is_file(), f"missing log-ray full dependency: {path}")
    priority = pilot.set_low_priority()
    fingerprint = cache_fingerprint()
    cached = load_cache(fingerprint)
    sector = json.loads(SECTOR.read_text(encoding="utf-8"))
    eq69 = json.loads(EQ69.read_text(encoding="utf-8"))
    w2 = json.loads(W2_W5.read_text(encoding="utf-8"))
    sector_rows = {int(row["chain"]): row for row in sector["rows"]}
    eq69_rows = {int(row["chain"]): row for row in eq69["rows"]}
    w2_rows = {int(row["chain"]): row for row in w2["rows"]}
    recursive = native_q.load_recursive_chains()
    require(set(recursive) == set(sector_rows) == set(eq69_rows) == set(w2_rows), "log-ray full roster drift")
    require(set(cached).issubset(recursive), "log-ray full cache has foreign chain")

    started = monotonic()
    for chain in sorted(recursive):
        if chain in cached:
            continue
        row = evaluate_row(chain, recursive[chain], sector_rows[chain], eq69_rows[chain], w2_rows[chain])
        append_cache(row)
        cached[chain] = row
        if monotonic() - started >= args.runtime_limit_seconds:
            print(
                f"parked logarithmic endpoint-ray full run after {len(cached)}/374 calls; "
                f"resume with the same command"
            )
            return 2

    rows = [cached[chain] for chain in sorted(cached)]
    require(len(rows) == 374, "log-ray full completion count drift")
    cutoff_counts = Counter(int(row["ray_cutoff"]) for row in rows)
    minimum_correction = min(
        (Fraction(row["exact_minus_paper_magnitude_lower"]), int(row["chain"])) for row in rows
    )
    maximum_correction = max(
        (Fraction(row["exact_minus_paper_magnitude_upper"]), int(row["chain"])) for row in rows
    )
    maximum_origin = max(Fraction(row["origin_radius_upper"]) for row in rows)
    maximum_tail = max(Fraction(row["tail_radius_upper"]) for row in rows)
    maximum_formula_gap = max(Fraction(row["formula_target_gap_magnitude_upper"]) for row in rows)
    maximum_correction_gap = max(Fraction(row["correction_residual_gap_magnitude_upper"]) for row in rows)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate",
        "status": "rigorous_all_374_logarithmic_endpoint_ray_nonsaddle_enclosures_validated",
        "scope": {
            "block": 20,
            "recursive_call_count": len(rows),
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "precision_overlap_count": len(rows),
            "worker_count": 1,
            "process_priority": priority,
            "cache_fingerprint": fingerprint,
        },
        "quadrature": {
            "origin_cutoff": decimal(pilot.ORIGIN_CUTOFF),
            "tail_target": decimal(TAIL_TARGET),
            "initial_ray_cutoff": int(INITIAL_CUTOFF),
            "cutoff_rule": "smallest dyadic multiple of 64 whose exact four-ray tail maximum is below 1e-30",
            "near_kernel": "principal-log identity (11.264.1)",
            "far_kernel": "K_m(w)=w^m*2F1(1,m;m+1;w)/m",
        },
        "aggregate": {
            "call_count": len(rows),
            "cutoff_distribution": {str(key): cutoff_counts[key] for key in sorted(cutoff_counts)},
            "maximum_ray_cutoff": max(cutoff_counts),
            "maximum_origin_radius_upper": decimal(maximum_origin),
            "maximum_tail_radius_upper": decimal(maximum_tail),
            "maximum_formula_target_gap_magnitude_upper": decimal(maximum_formula_gap),
            "maximum_correction_residual_gap_magnitude_upper": decimal(maximum_correction_gap),
            "formula_target_identity_count": sum(row["formula_target_identity_contains_zero"] for row in rows),
            "correction_residual_identity_count": sum(row["correction_residual_identity_contains_zero"] for row in rows),
            "exact_correction_excluding_zero_count": sum(row["exact_minus_paper_excludes_zero"] for row in rows),
            "minimum_exact_correction_magnitude_lower": decimal(minimum_correction[0]),
            "minimum_exact_correction_witness": minimum_correction[1],
            "maximum_exact_correction_magnitude_upper": decimal(maximum_correction[0]),
            "maximum_exact_correction_witness": maximum_correction[1],
        },
        "rows": rows,
        "sources": {
            "sector_gate": {"path": relative(SECTOR), "sha256": file_hash(SECTOR)},
            "log_reduction_gate": {"path": relative(LOG_REDUCTION), "sha256": file_hash(LOG_REDUCTION)},
            "equation_69_gate": {"path": relative(EQ69), "sha256": file_hash(EQ69)},
            "w2_w5_gate": {"path": relative(W2_W5), "sha256": file_hash(W2_W5)},
            "pilot_result": {"path": relative(PILOT_RESULT), "sha256": file_hash(PILOT_RESULT)},
            "cache": {"path": relative(CACHE), "sha256": file_hash(CACHE)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Derive a height-uniform analytic majorant for the exact full-cubic endpoint-ray value minus the published "
            "quadratic W2--W4 model, then propagate it through all recurrence levels and the outer Hardy formula."
        ),
        "proof_boundary": (
            "Rigorous all-call finite-input enclosure for one low-height block only. This does not prove a height-uniform "
            "W2--W4 remainder, recursive or outer Hardy error bound, Lambda<=0, PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built logarithmic endpoint-ray full Arb gate: 374/374 calls, "
        "both exact identities enclose zero"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
