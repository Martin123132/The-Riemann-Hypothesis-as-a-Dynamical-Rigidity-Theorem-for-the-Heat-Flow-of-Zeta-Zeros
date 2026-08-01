#!/usr/bin/env python3
"""Test a tuned step-four exact-anchor lattice across the compact interval."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache as source  # noqa: E402
import jensen_window_pf_compound_order11_compact_high_precision_anchor_reach_pilot as reach  # noqa: E402
import jensen_window_pf_compound_order11_compact_profile_geographic_pilot as single_geo  # noqa: E402
import jensen_window_pf_compound_order11_compact_profile_ladder_pilot as ladder  # noqa: E402
import jensen_window_pf_compound_order11_compact_stencil_precision_diagnostic as precision_diagnostic  # noqa: E402
import jensen_window_pf_compound_order11_sparse_h0_h23_propagation_core as propagation  # noqa: E402
from jensen_window_pf_compound_order11_shifted_taylor_model_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    PRECISION_BITS,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_step4_geographic_pilot.json"
)
DEFAULT_SOURCE_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_step4_geographic_h0_h23.jsonl"
)
PROFILE_NAME = "order11_ladder_p896_n96"
PROFILE = ladder.PROFILES[PROFILE_NAME]
CENTERS = single_geo.TARGETS
ANCHOR_SHIFTS = tuple(Fraction(value) for value in range(-12, 13, 4))
TEST_OFFSETS = (Fraction(-2), Fraction(0), Fraction(2))
STENCIL_RADIUS = 9
QUARTER = Fraction(1, 4)
MAXIMUM_PROPAGATION_DISTANCE = Fraction(2)
ANCHORS = tuple(
    center + shift for center in CENTERS for shift in ANCHOR_SHIFTS
)


def exact_task(task: tuple[int, Fraction]) -> dict:
    index, target = task
    source.profiles.PROFILE_SPECS[PROFILE_NAME] = PROFILE
    row = source.exact_task((index, target, PROFILE_NAME))
    return {**row, "kind": "order11_compact_step4_geographic_h0_h23_jet"}


def tasks() -> list[tuple[int, Fraction]]:
    return list(enumerate(ANCHORS))


def validate_source_row(
    record: dict,
    task: tuple[int, Fraction],
) -> None:
    index, target = task
    source.profiles.PROFILE_SPECS[PROFILE_NAME] = PROFILE
    precision_diagnostic.validate_point_record(record, target)
    if (
        record.get("kind") != "order11_compact_step4_geographic_h0_h23_jet"
        or record.get("index") != index
        or record.get("profile") != PROFILE_NAME
        or record.get("contract_id") != source.row_contract(PROFILE_NAME)
    ):
        raise RuntimeError(f"invalid step-four geographic source row {index}")


def load_source_cache(path: Path = DEFAULT_SOURCE_CACHE) -> list[dict]:
    if not path.exists():
        return []
    records = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise RuntimeError(
                    f"invalid step-four geographic JSONL row {line_number}"
                ) from exc
    expected = tasks()
    if len(records) > len(expected):
        raise RuntimeError("step-four geographic source cache has too many rows")
    for record, task in zip(records, expected):
        validate_source_row(record, task)
    return records


def build_source_cache(
    path: Path = DEFAULT_SOURCE_CACHE,
    *,
    workers: int = 4,
) -> list[dict]:
    expected = tasks()
    records = load_source_cache(path)
    remaining = expected[len(records) :]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=workers) as pool, path.open(
        "a", encoding="utf-8", newline="\n"
    ) as handle:
        for record, task in zip(pool.map(exact_task, remaining), remaining):
            validate_source_row(record, task)
            handle.write(json.dumps(record, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
            records.append(record)
    return records


def expansion_row(
    expansion: Fraction,
    anchor_pairs: dict[Fraction, tuple[list[flint.arb], dict]],
    h_rows: list[dict],
) -> dict:
    targets = [expansion + shift for shift in range(-STENCIL_RADIUS, STENCIL_RADIUS + 1)]
    old_limit = propagation.MAXIMUM_PROPAGATION_DISTANCE
    propagation.MAXIMUM_PROPAGATION_DISTANCE = MAXIMUM_PROPAGATION_DISTANCE
    try:
        propagated = propagation.propagated_point_source(
            targets,
            {anchor: pair[0] for anchor, pair in anchor_pairs.items()},
            h_rows,
            anchor_diagnostic_source={
                anchor: pair[1] for anchor, pair in anchor_pairs.items()
            },
        )
        chosen = [
            propagation.nearest_anchor(target, sorted(anchor_pairs))
            for target in targets
        ]
    finally:
        propagation.MAXIMUM_PROPAGATION_DISTANCE = old_limit
    blocks = [
        reach.evaluate_block(
            expansion,
            expansion - QUARTER,
            expansion,
            h_rows,
            propagated,
            "left",
        ),
        reach.evaluate_block(
            expansion,
            expansion,
            expansion + QUARTER,
            h_rows,
            propagated,
            "right",
        ),
    ]
    completed = [block for block in blocks if "scaled_curvature_upper" in block]
    largest = (
        max(completed, key=lambda block: float(block["scaled_curvature_upper"]))
        if completed
        else None
    )
    return {
        "expansion_anchor": str(expansion),
        "chosen_anchor_range": [str(min(chosen)), str(max(chosen))],
        "maximum_actual_propagation_distance": str(
            max(abs(target - anchor) for target, anchor in zip(targets, chosen))
        ),
        "blocks": blocks,
        "largest_scaled_curvature_upper": (
            largest["scaled_curvature_upper"] if largest else None
        ),
        "passed": all(block["passed"] for block in blocks),
    }


def center_evaluation(
    center: Fraction,
    records: dict[Fraction, dict],
    h_rows: list[dict],
) -> dict:
    anchors = tuple(center + shift for shift in ANCHOR_SHIFTS)
    pairs = {
        anchor: precision_diagnostic.point_pair(records[anchor])
        for anchor in anchors
    }
    expansion_rows = [
        expansion_row(center + offset, pairs, h_rows) for offset in TEST_OFFSETS
    ]
    completed = [
        block
        for row in expansion_rows
        for block in row["blocks"]
        if "scaled_curvature_upper" in block
    ]
    largest = (
        max(completed, key=lambda block: float(block["scaled_curvature_upper"]))
        if completed
        else None
    )
    return {
        "center": str(center),
        "anchor_grid": [str(anchor) for anchor in anchors],
        "source_rows": [records[anchor] for anchor in anchors],
        "expansion_rows": expansion_rows,
        "largest_scaled_curvature_upper": (
            largest["scaled_curvature_upper"] if largest else None
        ),
        "passed": all(row["passed"] for row in expansion_rows),
    }


def build_artifact(*, workers: int = 4) -> dict:
    rows = build_source_cache(workers=workers)
    records = {Fraction(row["target_t"]): row for row in rows}
    if set(records) != set(ANCHORS) or len(rows) != len(ANCHORS):
        raise RuntimeError("step-four geographic exact source grid changed")
    flint.ctx.prec = PRECISION_BITS
    h_rows = single_geo.load_h_rows()
    evaluations = [
        center_evaluation(center, records, h_rows) for center in CENTERS
    ]
    passing = [row["center"] for row in evaluations if row["passed"]]
    return {
        "kind": "jensen_window_pf_compound_order11_compact_step4_geographic_pilot",
        "date": "2026-07-18",
        "status": "rigorous representative step-four compact lattice pilot",
        "proof_boundary": (
            "This covers six quarter blocks near each listed center. It tests the "
            "proposed local lattice contract but not every compact-interval block."
        ),
        "parameters": {
            "centers": [str(center) for center in CENTERS],
            "anchor_shifts": [str(shift) for shift in ANCHOR_SHIFTS],
            "test_offsets": [str(offset) for offset in TEST_OFFSETS],
            "maximum_propagation_distance": str(MAXIMUM_PROPAGATION_DISTANCE),
            "profile_name": PROFILE_NAME,
            "profile": PROFILE,
            "curvature_cap": CURVATURE_CONSTANT,
            "model_working_precision_bits": PRECISION_BITS,
        },
        "evaluations": evaluations,
        "summary": {
            "centers": len(evaluations),
            "exact_source_rows": len(rows),
            "passing_centers": passing,
            "all_centers_passed": len(passing) == len(evaluations),
            "scaled_curvature_uppers": {
                row["center"]: row["largest_scaled_curvature_upper"]
                for row in evaluations
            },
        },
        "source_cache": {
            "path": DEFAULT_SOURCE_CACHE.relative_to(REPO_ROOT).as_posix(),
            "sha256": reach.sha256(DEFAULT_SOURCE_CACHE),
            "rows": len(rows),
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order11_compact_step4_geographic_pilot.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order11_compact_step4_geographic_pilot.py"
        ),
    }


def main() -> int:
    artifact = build_artifact()
    DEFAULT_OUT.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUT.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    summary = artifact["summary"]
    print(
        "compact step-four geographic pilot: "
        f"passing={len(summary['passing_centers'])}/{summary['centers']}, "
        f"uppers={summary['scaled_curvature_uppers']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
