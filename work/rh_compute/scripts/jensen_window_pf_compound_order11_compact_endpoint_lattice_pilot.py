#!/usr/bin/env python3
"""Certify local step-four lattice pilots at both compact endpoints."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
import json
import math
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

import jensen_window_pf_compound_order10_compact_h2_h24_unit_cache as h_source  # noqa: E402
import jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache as point_source  # noqa: E402
import jensen_window_pf_compound_order11_compact_h2_h24_boundary_extension as boundary_source  # noqa: E402
import jensen_window_pf_compound_order11_compact_high_precision_anchor_reach_pilot as reach  # noqa: E402
import jensen_window_pf_compound_order11_compact_profile_ladder_pilot as ladder  # noqa: E402
import jensen_window_pf_compound_order11_compact_stencil_precision_diagnostic as precision_diagnostic  # noqa: E402
import jensen_window_pf_compound_order11_compact_upper_profile_lattice_pilot as upper  # noqa: E402
import jensen_window_pf_compound_order11_sparse_h0_h23_propagation_core as propagation  # noqa: E402
from jensen_window_pf_compound_order11_shifted_taylor_model_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    PRECISION_BITS,
)


DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_endpoint_lattice_h0_h23.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_endpoint_lattice_pilot.json"
)
LOWER_PROFILE_NAME = "order11_ladder_p896_n96"
UPPER_PROFILE_NAME = "order11_upper_p896_n108"
PROFILES = {
    LOWER_PROFILE_NAME: ladder.PROFILES[LOWER_PROFILE_NAME],
    UPPER_PROFILE_NAME: upper.PROFILES[UPPER_PROFILE_NAME],
}
LOWER_ANCHORS = tuple(Fraction(value) for value in range(5692, 5713, 4))
UPPER_ANCHORS = tuple(Fraction(value) for value in range(38008, 38029, 4))
LOWER_EXPANSIONS = tuple(Fraction(5700) + Fraction(index, 2) for index in range(5))
UPPER_EXPANSIONS = tuple(Fraction(38018) + Fraction(index, 2) for index in range(5))
STENCIL_RADIUS = 9
QUARTER = Fraction(1, 4)
MAXIMUM_PROPAGATION_DISTANCE = Fraction(2)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def tasks() -> list[tuple[int, str, Fraction]]:
    result = []
    for profile_name, anchors in (
        (LOWER_PROFILE_NAME, LOWER_ANCHORS),
        (UPPER_PROFILE_NAME, UPPER_ANCHORS),
    ):
        for target in anchors:
            result.append((len(result), profile_name, target))
    return result


def exact_task(task: tuple[int, str, Fraction]) -> dict:
    index, profile_name, target = task
    point_source.profiles.PROFILE_SPECS[profile_name] = PROFILES[profile_name]
    row = point_source.exact_task((index, target, profile_name))
    return {**row, "kind": "order11_compact_endpoint_lattice_h0_h23_jet"}


def validate_source_row(record: dict, task: tuple[int, str, Fraction]) -> None:
    index, profile_name, target = task
    point_source.profiles.PROFILE_SPECS[profile_name] = PROFILES[profile_name]
    precision_diagnostic.validate_point_record(record, target)
    if (
        record.get("kind") != "order11_compact_endpoint_lattice_h0_h23_jet"
        or record.get("index") != index
        or record.get("profile") != profile_name
        or record.get("contract_id") != point_source.row_contract(profile_name)
    ):
        raise RuntimeError(f"invalid compact endpoint source row {index}")


def load_source_cache(path: Path = DEFAULT_CACHE) -> list[dict]:
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
                    f"invalid compact endpoint JSONL row {line_number}"
                ) from exc
    expected = tasks()
    if len(records) > len(expected):
        raise RuntimeError("compact endpoint source cache has too many rows")
    for record, task in zip(records, expected):
        validate_source_row(record, task)
    return records


def build_source_cache(*, workers: int, generate: bool) -> list[dict]:
    expected = tasks()
    records = load_source_cache()
    if len(records) == len(expected):
        return records
    if not generate:
        raise RuntimeError(
            f"compact endpoint source cache has {len(records)}/{len(expected)} rows"
        )
    DEFAULT_CACHE.parent.mkdir(parents=True, exist_ok=True)
    remaining = expected[len(records) :]
    with ProcessPoolExecutor(max_workers=workers) as pool, DEFAULT_CACHE.open(
        "a", encoding="utf-8", newline="\n"
    ) as handle:
        for record, task in zip(pool.map(exact_task, remaining), remaining):
            validate_source_row(record, task)
            handle.write(json.dumps(record, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
            records.append(record)
    return records


def parse_h_record(record: dict) -> tuple[Fraction, dict]:
    left = Fraction(record["target_t_left"])
    right = Fraction(record["target_t_right"])
    derivatives = record.get("h_derivatives", {})
    if (
        record.get("contract_id") != h_source.ROW_CONTRACT
        or record.get("passed") is not True
        or right != left + 1
        or set(derivatives) != {str(order) for order in range(2, 25)}
    ):
        raise RuntimeError(f"invalid endpoint H row at t={left}")
    return (
        left,
        {
            "target_t_left": left,
            "target_t_right": right,
            "H": {
                order: precision_diagnostic.compact.interval_from_text(
                    derivatives[str(order)]
                )
                for order in range(2, 25)
            },
        },
    )


def load_h_rows() -> list[dict]:
    main_targets = {
        Fraction(value)
        for value in (*range(5692, 5712), *range(38008, 38028))
    }
    rows: dict[Fraction, dict] = {}
    with h_source.DEFAULT_CACHE.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            left = Fraction(record["target_t_left"])
            if left in main_targets:
                key, parsed = parse_h_record(record)
                rows[key] = parsed
    boundary = json.loads(boundary_source.DEFAULT_OUT.read_text(encoding="utf-8"))
    for record in boundary.get("rows", []):
        key, parsed = parse_h_record(record)
        rows[key] = parsed
    expected = main_targets | {Fraction(5691), Fraction(38028)}
    if set(rows) != expected:
        raise RuntimeError("compact endpoint H collar is incomplete")
    return [rows[target] for target in sorted(rows)]


def expansion_row(
    expansion: Fraction,
    anchor_pairs: dict[Fraction, tuple[list[flint.arb], dict]],
    h_rows: list[dict],
    sides: tuple[str, ...],
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
    blocks = []
    if "left" in sides:
        blocks.append(
            reach.evaluate_block(
                expansion,
                expansion - QUARTER,
                expansion,
                h_rows,
                propagated,
                "left",
            )
        )
    if "right" in sides:
        blocks.append(
            reach.evaluate_block(
                expansion,
                expansion,
                expansion + QUARTER,
                h_rows,
                propagated,
                "right",
            )
        )
    completed = [block for block in blocks if "scaled_curvature_upper" in block]
    largest = (
        max(completed, key=lambda block: float(block["scaled_curvature_upper"]))
        if completed
        else None
    )
    return {
        "expansion_anchor": str(expansion),
        "sides": list(sides),
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


def endpoint_evaluation(
    name: str,
    profile_name: str,
    anchors: tuple[Fraction, ...],
    expansions: tuple[Fraction, ...],
    records: dict[Fraction, dict],
    h_rows: list[dict],
) -> dict:
    pairs = {
        anchor: precision_diagnostic.point_pair(records[anchor])
        for anchor in anchors
    }
    rows = []
    for index, expansion in enumerate(expansions):
        sides = ("left", "right")
        if name == "lower" and index == 0:
            sides = ("right",)
        if name == "upper" and index == len(expansions) - 1:
            sides = ("left",)
        rows.append(expansion_row(expansion, pairs, h_rows, sides))
    completed = [
        block
        for row in rows
        for block in row["blocks"]
        if "scaled_curvature_upper" in block
    ]
    largest = max(completed, key=lambda block: float(block["scaled_curvature_upper"]))
    return {
        "name": name,
        "profile_name": profile_name,
        "profile": PROFILES[profile_name],
        "anchors": [str(anchor) for anchor in anchors],
        "expansion_rows": rows,
        "quarter_blocks": sum(len(row["blocks"]) for row in rows),
        "largest_scaled_curvature_upper": largest["scaled_curvature_upper"],
        "passed": all(row["passed"] for row in rows),
    }


def build_artifact(*, workers: int = 4, generate: bool = True) -> dict:
    source_rows = build_source_cache(workers=workers, generate=generate)
    records = {Fraction(record["target_t"]): record for record in source_rows}
    if set(records) != set(LOWER_ANCHORS + UPPER_ANCHORS):
        raise RuntimeError("compact endpoint exact source grid changed")
    flint.ctx.prec = PRECISION_BITS
    h_rows = load_h_rows()
    endpoints = [
        endpoint_evaluation(
            "lower",
            LOWER_PROFILE_NAME,
            LOWER_ANCHORS,
            LOWER_EXPANSIONS,
            records,
            h_rows,
        ),
        endpoint_evaluation(
            "upper",
            UPPER_PROFILE_NAME,
            UPPER_ANCHORS,
            UPPER_EXPANSIONS,
            records,
            h_rows,
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order11_compact_endpoint_lattice_pilot",
        "date": "2026-07-18",
        "status": "rigorous compact endpoint step-four lattice pilot",
        "proof_boundary": (
            "This covers only the listed endpoint quarter blocks. It validates "
            "the boundary collars but is not contiguous compact coverage."
        ),
        "parameters": {
            "compact_target": ["5700", "38020"],
            "lower_profile": PROFILES[LOWER_PROFILE_NAME],
            "upper_profile": PROFILES[UPPER_PROFILE_NAME],
            "maximum_propagation_distance": str(MAXIMUM_PROPAGATION_DISTANCE),
            "curvature_cap": CURVATURE_CONSTANT,
            "model_working_precision_bits": PRECISION_BITS,
        },
        "endpoints": endpoints,
        "sources": {
            "exact_cache": {
                "path": relative(DEFAULT_CACHE),
                "sha256": sha256(DEFAULT_CACHE),
                "rows": len(source_rows),
            },
            "h_boundary_extension": {
                "path": relative(boundary_source.DEFAULT_OUT),
                "sha256": sha256(boundary_source.DEFAULT_OUT),
            },
        },
        "summary": {
            "endpoint_regions": len(endpoints),
            "quarter_blocks": sum(row["quarter_blocks"] for row in endpoints),
            "all_endpoints_passed": all(row["passed"] for row in endpoints),
            "scaled_curvature_uppers": {
                row["name"]: row["largest_scaled_curvature_upper"]
                for row in endpoints
            },
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order11_compact_endpoint_lattice_pilot.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order11_compact_endpoint_lattice_pilot.py"
        ),
    }


def main() -> int:
    artifact = build_artifact(generate=True)
    DEFAULT_OUT.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUT.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    summary = artifact["summary"]
    print(
        "compact endpoint lattice: "
        f"passed={summary['all_endpoints_passed']}, "
        f"blocks={summary['quarter_blocks']}, "
        f"uppers={summary['scaled_curvature_uppers']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
