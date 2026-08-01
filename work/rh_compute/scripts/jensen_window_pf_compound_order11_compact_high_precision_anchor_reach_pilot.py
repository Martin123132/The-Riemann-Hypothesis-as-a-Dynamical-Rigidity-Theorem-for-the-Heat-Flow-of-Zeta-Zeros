#!/usr/bin/env python3
"""Measure one high-precision H0-H23 anchor's rigorous local reach."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
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
import jensen_window_pf_compound_order11_compact_high_precision_stencil_pilot as point_source  # noqa: E402
import jensen_window_pf_compound_order11_compact_stencil_precision_diagnostic as precision_diagnostic  # noqa: E402
import jensen_window_pf_compound_order11_sparse_h0_h23_propagation_core as propagation  # noqa: E402
from jensen_window_pf_compound_order11_shifted_taylor_model_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    PRECISION_BITS,
    shifted_taylor_model_curvature_row,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_high_precision_anchor_reach_pilot.json"
)
ANCHOR = Fraction(20004)
MINIMUM_OFFSET = Fraction(-16)
MAXIMUM_OFFSET = Fraction(16)
OFFSET_STEP = Fraction(1, 2)
QUARTER = Fraction(1, 4)
STENCIL_RADIUS = 9
H_LEFT = ANCHOR + MINIMUM_OFFSET - STENCIL_RADIUS - 1
H_RIGHT = ANCHOR + MAXIMUM_OFFSET + STENCIL_RADIUS + 1


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def offsets() -> tuple[Fraction, ...]:
    count = int((MAXIMUM_OFFSET - MINIMUM_OFFSET) / OFFSET_STEP)
    return tuple(MINIMUM_OFFSET + index * OFFSET_STEP for index in range(count + 1))


def central_record() -> dict:
    artifact = json.loads(point_source.DEFAULT_OUT.read_text(encoding="utf-8"))
    matches = [
        row for row in artifact.get("rows", []) if row.get("target_t") == str(ANCHOR)
    ]
    if len(matches) != 1:
        raise RuntimeError("high-precision source has no unique central anchor")
    precision_diagnostic.validate_point_record(matches[0], ANCHOR)
    return matches[0]


def load_h_rows() -> list[dict]:
    selected = {Fraction(value) for value in range(int(H_LEFT), int(H_RIGHT))}
    derivative_keys = {str(order) for order in range(2, 25)}
    rows: dict[Fraction, dict] = {}
    with h_source.DEFAULT_CACHE.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            left = Fraction(record["target_t_left"])
            if left not in selected:
                continue
            right = Fraction(record["target_t_right"])
            derivatives = record.get("h_derivatives", {})
            if (
                record.get("kind") != "order10_compact_h2_h24_unit_tile"
                or record.get("contract_id") != h_source.ROW_CONTRACT
                or record.get("passed") is not True
                or right != left + 1
                or set(derivatives) != derivative_keys
            ):
                raise RuntimeError(f"invalid compact H row at t={left}")
            rows[left] = {
                "target_t_left": left,
                "target_t_right": right,
                "H": {
                    order: precision_diagnostic.compact.interval_from_text(
                        derivatives[str(order)]
                    )
                    for order in range(2, 25)
                },
            }
    if set(rows) != selected:
        raise RuntimeError("compact H source misses the anchor-reach collar")
    return [rows[target] for target in sorted(rows)]


def compact_block(row: dict, side: str) -> dict:
    return {
        "side": side,
        "anchor": row["anchor"],
        "expansion_anchor": row["expansion_anchor"],
        "right": row["right"],
        "scaled_curvature_upper": row["scaled_curvature_upper"],
        "point_scaled_curvature": row["point_scaled_curvature"],
        "curvature_margin_lower": row["curvature_margin_lower"],
        "final_uniform_remainder_upper": row["final_uniform_remainder_upper"],
        "passed": row["passed"],
    }


def evaluate_block(
    expansion: Fraction,
    block_left: Fraction,
    block_right: Fraction,
    h_rows: list[dict],
    propagated: dict,
    side: str,
) -> dict:
    try:
        row = shifted_taylor_model_curvature_row(
            expansion,
            block_left,
            block_right,
            h_rows,
            point_h_source=propagated,
            require_pass=False,
        )
    except (ValueError, RuntimeError) as exc:
        return {
            "side": side,
            "anchor": str(block_left),
            "expansion_anchor": str(expansion),
            "right": str(block_right),
            "passed": False,
            "failure": f"{type(exc).__name__}: {exc}",
        }
    return compact_block(row, side)


def expansion_row(
    expansion: Fraction,
    anchor_pair: tuple[list[flint.arb], dict],
    h_rows: list[dict],
    *,
    exact_anchor: Fraction = ANCHOR,
) -> dict:
    targets = [expansion + shift for shift in range(-STENCIL_RADIUS, STENCIL_RADIUS + 1)]
    maximum_distance = max(abs(target - exact_anchor) for target in targets)
    old_limit = propagation.MAXIMUM_PROPAGATION_DISTANCE
    propagation.MAXIMUM_PROPAGATION_DISTANCE = maximum_distance
    try:
        propagated = propagation.propagated_point_source(
            targets,
            {exact_anchor: anchor_pair[0]},
            h_rows,
            anchor_diagnostic_source={exact_anchor: anchor_pair[1]},
        )
    finally:
        propagation.MAXIMUM_PROPAGATION_DISTANCE = old_limit
    left = evaluate_block(
        expansion,
        expansion - QUARTER,
        expansion,
        h_rows,
        propagated,
        "left",
    )
    right = evaluate_block(
        expansion,
        expansion,
        expansion + QUARTER,
        h_rows,
        propagated,
        "right",
    )
    blocks = [left, right]
    completed = [block for block in blocks if "scaled_curvature_upper" in block]
    largest = (
        max(completed, key=lambda block: float(block["scaled_curvature_upper"]))
        if completed
        else None
    )
    return {
        "offset": str(expansion - exact_anchor),
        "expansion_anchor": str(expansion),
        "maximum_propagation_distance": str(maximum_distance),
        "blocks": blocks,
        "largest_scaled_curvature_upper": (
            largest["scaled_curvature_upper"] if largest else None
        ),
        "passed": all(block["passed"] for block in blocks),
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    record = central_record()
    anchor_pair = precision_diagnostic.point_pair(record)
    h_rows = load_h_rows()
    rows = [expansion_row(ANCHOR + offset, anchor_pair, h_rows) for offset in offsets()]
    completed = [row for row in rows if row["largest_scaled_curvature_upper"]]
    largest = max(
        completed,
        key=lambda row: float(row["largest_scaled_curvature_upper"]),
    )
    passing_offsets = [Fraction(row["offset"]) for row in rows if row["passed"]]
    return {
        "kind": "jensen_window_pf_compound_order11_compact_high_precision_anchor_reach_pilot",
        "date": "2026-07-18",
        "status": "rigorous single-anchor compact reach pilot",
        "proof_boundary": (
            "This covers only quarter blocks around half-grid expansions within "
            "sixteen units of t=20004; it does not certify the full compact interval."
        ),
        "theorem": (
            "the listed quarter blocks are rigorously below the order-eleven "
            f"curvature cap {CURVATURE_CONSTANT} exactly where passed=true"
        ),
        "parameters": {
            "anchor": str(ANCHOR),
            "offset_range": [str(MINIMUM_OFFSET), str(MAXIMUM_OFFSET)],
            "offset_step": str(OFFSET_STEP),
            "quarter_block_width": str(QUARTER),
            "stencil_radius": STENCIL_RADIUS,
            "working_precision_bits": PRECISION_BITS,
            "curvature_cap": CURVATURE_CONSTANT,
            "exact_anchor_profile": point_source.PROFILE,
        },
        "sources": {
            "exact_anchor_stencil": {
                "path": relative(point_source.DEFAULT_OUT),
                "sha256": sha256(point_source.DEFAULT_OUT),
                "selected_target": str(ANCHOR),
            },
            "h2_h24_unit_cache": {
                "path": relative(h_source.DEFAULT_CACHE),
                "sha256": sha256(h_source.DEFAULT_CACHE),
                "selected_tile_range": [str(H_LEFT), str(H_RIGHT - 1)],
            },
        },
        "rows": rows,
        "summary": {
            "expansion_rows": len(rows),
            "quarter_blocks": 2 * len(rows),
            "passing_expansion_rows": len(passing_offsets),
            "all_blocks_passed": all(row["passed"] for row in rows),
            "passing_offset_range": (
                [str(min(passing_offsets)), str(max(passing_offsets))]
                if passing_offsets
                else None
            ),
            "largest_scaled_curvature_upper": largest[
                "largest_scaled_curvature_upper"
            ],
            "largest_scaled_curvature_offset": largest["offset"],
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order11_compact_high_precision_anchor_reach_pilot.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order11_compact_high_precision_anchor_reach_pilot.py"
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
        "compact high-precision anchor reach: "
        f"{summary['passing_expansion_rows']}/{summary['expansion_rows']} expansions, "
        f"largest={summary['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
