#!/usr/bin/env python3
"""Compare ordinary and high-precision compact order-eleven point stencils."""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
import jensen_window_pf_compound_order10_compact_h2_h24_unit_cache as h_source  # noqa: E402
import jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache as sparse_source  # noqa: E402
import jensen_window_pf_compound_order11_compact_high_precision_stencil_pilot as high_source  # noqa: E402
import jensen_window_pf_compound_order11_compact_step4_pilot_anchors as step4_source  # noqa: E402
import jensen_window_pf_compound_order11_compact_unit_pilot_anchors as unit_source  # noqa: E402
import jensen_window_pf_compound_order11_sparse_h0_h23_propagation_core as propagation  # noqa: E402
from jensen_window_pf_compound_order11_shifted_taylor_model_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    PRECISION_BITS,
    shifted_taylor_model_curvature_row,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_stencil_precision_diagnostic.json"
)
EXPANSION_ANCHOR = Fraction(20004)
BLOCK_LEFT = EXPANSION_ANCHOR
BLOCK_RIGHT = EXPANSION_ANCHOR + Fraction(1, 4)
POINT_TARGETS = tuple(EXPANSION_ANCHOR + shift for shift in range(-9, 10))
H_TARGETS = POINT_TARGETS
PROPAGATED_ANCHOR_GRIDS = {
    "step2": tuple(Fraction(value) for value in range(19996, 20014, 2)),
    "step4": tuple(Fraction(value) for value in range(19996, 20013, 4)),
    "step8": (Fraction(19996), Fraction(20004), Fraction(20012)),
    "central": (EXPANSION_ANCHOR,),
}
PROPAGATION_LIMITS = {
    "step2": Fraction(1),
    "step4": Fraction(2),
    "step8": Fraction(4),
    "central": Fraction(9),
}


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_selected_h_rows() -> list[dict]:
    selected = set(H_TARGETS)
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
                    order: compact.interval_from_text(derivatives[str(order)])
                    for order in range(2, 25)
                },
            }
    if set(rows) != selected:
        missing = sorted(selected - set(rows))
        raise RuntimeError(f"compact H source misses selected rows: {missing}")
    return [rows[target] for target in H_TARGETS]


def validate_point_record(record: dict, target: Fraction) -> None:
    if (
        record.get("target_t") != str(target)
        or record.get("passed") is not True
        or set(record.get("h_derivatives", {}))
        != {str(order) for order in range(24)}
    ):
        raise RuntimeError(f"invalid exact H0-H23 point row at t={target}")


def point_pair(record: dict) -> tuple[list[flint.arb], dict]:
    target = Fraction(record["target_t"])
    validate_point_record(record, target)
    derivatives = record["h_derivatives"]
    return (
        [
            compact.interval_from_text(derivatives[str(order)])
            / math.factorial(order)
            for order in range(24)
        ],
        {
            "target_t": str(target),
            "mode_bracket": [record["mode_left"], record["mode_right"]],
            "maximum_panel_error_upper": record["maximum_panel_error_upper"],
            "maximum_tail_moment_upper": record["maximum_tail_moment_upper"],
            "minimum_tail_slope_lower": record["minimum_tail_slope_lower"],
        },
    )


def rows_by_target(rows: list[dict]) -> dict[Fraction, dict]:
    result: dict[Fraction, dict] = {}
    for record in rows:
        target = Fraction(record["target_t"])
        validate_point_record(record, target)
        if target in result:
            raise RuntimeError(f"duplicate exact point row at t={target}")
        result[target] = record
    return result


def load_sparse_selected() -> dict[Fraction, dict]:
    selected = set(POINT_TARGETS)
    result: dict[Fraction, dict] = {}
    with sparse_source.DEFAULT_CACHE.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            target = Fraction(record["target_t"])
            if target in selected:
                validate_point_record(record, target)
                result[target] = record
    return result


def ordinary_records() -> dict[Fraction, dict]:
    result = load_sparse_selected()
    for artifact_path in (step4_source.DEFAULT_OUT, unit_source.DEFAULT_OUT):
        artifact = load_json(artifact_path)
        for target, record in rows_by_target(artifact.get("rows", [])).items():
            if target not in POINT_TARGETS:
                continue
            if target in result:
                raise RuntimeError(f"ordinary stencil duplicates t={target}")
            result[target] = record
    if set(result) != set(POINT_TARGETS):
        missing = sorted(set(POINT_TARGETS) - set(result))
        extra = sorted(set(result) - set(POINT_TARGETS))
        raise RuntimeError(
            f"ordinary exact stencil mismatch: missing={missing}, extra={extra}"
        )
    return result


def high_precision_records() -> dict[Fraction, dict]:
    result = rows_by_target(load_json(high_source.DEFAULT_OUT).get("rows", []))
    if set(result) != set(POINT_TARGETS):
        raise RuntimeError("high-precision exact stencil is not the required grid")
    return result


def point_source(records: dict[Fraction, dict]) -> dict[Fraction, tuple[list, dict]]:
    return {target: point_pair(records[target]) for target in POINT_TARGETS}


def model_row(
    records: dict[Fraction, dict],
    h_rows: list[dict],
) -> dict:
    return shifted_taylor_model_curvature_row(
        EXPANSION_ANCHOR,
        BLOCK_LEFT,
        BLOCK_RIGHT,
        h_rows,
        point_h_source=point_source(records),
        require_pass=False,
    )


def propagated_model_row(
    records: dict[Fraction, dict],
    h_rows: list[dict],
    grid_name: str,
) -> dict:
    anchor_targets = PROPAGATED_ANCHOR_GRIDS[grid_name]
    pairs = {target: point_pair(records[target]) for target in anchor_targets}
    anchor_source = {target: pair[0] for target, pair in pairs.items()}
    diagnostics = {target: pair[1] for target, pair in pairs.items()}
    old_limit = propagation.MAXIMUM_PROPAGATION_DISTANCE
    propagation.MAXIMUM_PROPAGATION_DISTANCE = PROPAGATION_LIMITS[grid_name]
    try:
        propagated = propagation.propagated_point_source(
            list(POINT_TARGETS),
            anchor_source,
            h_rows,
            anchor_diagnostic_source=diagnostics,
        )
    finally:
        propagation.MAXIMUM_PROPAGATION_DISTANCE = old_limit
    row = shifted_taylor_model_curvature_row(
        EXPANSION_ANCHOR,
        BLOCK_LEFT,
        BLOCK_RIGHT,
        h_rows,
        point_h_source=propagated,
        require_pass=False,
    )
    return {
        "anchor_grid": [str(target) for target in anchor_targets],
        "anchor_count": len(anchor_targets),
        "maximum_propagation_distance": str(PROPAGATION_LIMITS[grid_name]),
        "row": row,
    }


def source_descriptor(path: Path) -> dict:
    return {"path": relative(path), "sha256": sha256(path)}


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    h_rows = load_selected_h_rows()
    ordinary = model_row(ordinary_records(), h_rows)
    high_records = high_precision_records()
    high_precision = model_row(high_records, h_rows)
    propagated_rows = {
        name: propagated_model_row(high_records, h_rows, name)
        for name in PROPAGATED_ANCHOR_GRIDS
    }
    with localcontext() as context:
        context.prec = 80
        ordinary_upper = Decimal(ordinary["scaled_curvature_upper"])
        high_upper = Decimal(high_precision["scaled_curvature_upper"])
        reduction = ordinary_upper / high_upper
    high_passed = high_precision.get("passed") is True
    passing_propagated_grids = [
        name
        for name, payload in propagated_rows.items()
        if payload["row"].get("passed") is True
    ]
    return {
        "kind": "jensen_window_pf_compound_order11_compact_stencil_precision_diagnostic",
        "date": "2026-07-18",
        "status": "rigorous local source-precision diagnostic",
        "proof_boundary": (
            "This compares two exact point-source profiles on one quarter block; "
            "it is not a compact-interval certificate."
        ),
        "domain": {
            "expansion_anchor": str(EXPANSION_ANCHOR),
            "block": [str(BLOCK_LEFT), str(BLOCK_RIGHT)],
            "point_stencil": [str(POINT_TARGETS[0]), str(POINT_TARGETS[-1])],
        },
        "parameters": {
            "working_precision_bits": PRECISION_BITS,
            "curvature_cap": CURVATURE_CONSTANT,
            "ordinary_profiles": ["medium"],
            "high_precision_profile": high_source.PROFILE,
        },
        "sources": {
            "h2_h24_unit_cache": source_descriptor(h_source.DEFAULT_CACHE),
            "ordinary_step8_cache": source_descriptor(sparse_source.DEFAULT_CACHE),
            "ordinary_step4_pilot": source_descriptor(step4_source.DEFAULT_OUT),
            "ordinary_unit_pilot": source_descriptor(unit_source.DEFAULT_OUT),
            "high_precision_stencil": source_descriptor(high_source.DEFAULT_OUT),
        },
        "ordinary_unit_stencil": ordinary,
        "high_precision_stencil": high_precision,
        "high_precision_propagated_stencils": propagated_rows,
        "comparison": {
            "ordinary_scaled_curvature_upper": ordinary["scaled_curvature_upper"],
            "high_precision_scaled_curvature_upper": high_precision[
                "scaled_curvature_upper"
            ],
            "high_precision_point_scaled_curvature": high_precision[
                "point_scaled_curvature"
            ],
            "ordinary_to_high_precision_reduction_factor": format(reduction, "E"),
            "high_precision_passed_cap": high_passed,
            "propagated_scaled_curvature_uppers": {
                name: payload["row"]["scaled_curvature_upper"]
                for name, payload in propagated_rows.items()
            },
            "passing_propagated_grids": passing_propagated_grids,
            "diagnosis": (
                "source precision closes this pilot quarter block"
                if high_passed
                else "source precision alone does not close this pilot quarter block"
            ),
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order11_compact_stencil_precision_diagnostic.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order11_compact_stencil_precision_diagnostic.py"
        ),
    }


def main() -> int:
    artifact = build_artifact()
    DEFAULT_OUT.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUT.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    comparison = artifact["comparison"]
    print(
        "compact precision diagnostic: "
        f"ordinary={comparison['ordinary_scaled_curvature_upper']}, "
        f"high={comparison['high_precision_scaled_curvature_upper']}, "
        f"passed={comparison['high_precision_passed_cap']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
