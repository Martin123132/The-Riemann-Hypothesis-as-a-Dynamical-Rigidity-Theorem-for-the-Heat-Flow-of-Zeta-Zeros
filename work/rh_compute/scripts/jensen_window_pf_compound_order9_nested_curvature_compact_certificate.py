#!/usr/bin/env python3
"""Certify order-nine first-summand curvature from t=5700 to mode u=2."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as order4  # noqa: E402
import jensen_window_pf_compound_order6_nested_curvature_compact_certificate as order6  # noqa: E402
import jensen_window_pf_compound_order7_nested_curvature_compact_certificate as order7  # noqa: E402
import jensen_window_pf_compound_order8_nested_curvature_compact_h13_h14_cache as order8  # noqa: E402
import jensen_window_pf_compound_order9_nested_curvature_compact_h15_h16_cache as order9  # noqa: E402
from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    integrate_h_derivatives,
)
from jensen_window_pf_compound_order9_nested_curvature_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    evaluate_nested_curvature_from_h_cover,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


DEFAULT_RIGHT_COLLAR = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_nested_curvature_compact_right_collar.json"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_nested_curvature_compact_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order9_nested_curvature_compact_certificate.md"
)
TARGET_T = 5700
COLLAR_T = 7
INITIAL_CENTRAL_TILE_COUNT = 1000
RIGHT_COLLAR_TILE_COUNT = 3


@dataclass(frozen=True)
class CompactRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_tasks() -> list[tuple[int, Fraction, Fraction]]:
    return order9.deterministic_tasks()


def build_right_collar(path: Path, *, overwrite: bool = False) -> dict:
    if path.exists() and not overwrite:
        return json.loads(path.read_text(encoding="utf-8"))
    tasks = deterministic_tasks()
    flint.ctx.prec = order4.DEFAULT_PRECISION_BITS
    records = []
    for offset in range(RIGHT_COLLAR_TILE_COUNT):
        index = len(tasks) + offset
        left = order4.OUTER_END + offset * order4.TILE_WIDTH
        right = left + order4.TILE_WIDTH
        result = integrate_h_derivatives(
            left,
            right,
            order4.PANELS,
            window_y=order4.WINDOW_Y,
            eighth_envelope_bound=order4.EIGHTH_ENVELOPE,
            max_moment=16,
        )
        if not result.get("passed"):
            raise RuntimeError(f"order-nine right collar failed: {result}")
        records.append(
            {
                "kind": "order9_compact_full_right_collar_tile",
                "index": index,
                "mode_left": str(left),
                "mode_right": str(right),
                "passed": True,
                "t_left": order4.interval_text(
                    potential_jet_arb(arb_rational(left), 1)[1]
                ),
                "t_right": order4.interval_text(
                    potential_jet_arb(arb_rational(right), 1)[1]
                ),
                "h_derivatives": {
                    str(degree): order4.interval_text(
                        result["h_derivatives"][degree]
                    )
                    for degree in range(2, 17)
                },
                "normalizer_lower": result["normalizer_lower"],
                "minimum_tail_slope_lower": result[
                    "minimum_tail_slope_lower"
                ],
                "maximum_tail_moment_upper": result[
                    "maximum_tail_moment_upper"
                ],
                "maximum_simpson_error_upper": result[
                    "maximum_simpson_error_upper"
                ],
            }
        )
    artifact = {
        "kind": "order9_compact_full_right_collar",
        "date": "2026-07-14",
        "tile_count": len(records),
        "all_rows_passed": True,
        "rows": records,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return artifact


def append_right_collar(
    caches: list[list[dict]],
    collar: dict,
) -> list[list[dict]]:
    expected = len(caches[0])
    records = collar.get("rows", [])
    if (
        collar.get("tile_count") != RIGHT_COLLAR_TILE_COUNT
        or collar.get("all_rows_passed") is not True
        or len(records) != RIGHT_COLLAR_TILE_COUNT
    ):
        raise RuntimeError("order-nine right collar is incomplete")
    order_groups = (tuple(range(2, 9)), (9, 10), (11, 12), (13, 14), (15, 16))
    result = [list(records) for records in caches]
    for offset, record in enumerate(records):
        if record.get("index") != expected + offset or record.get("passed") is not True:
            raise RuntimeError("order-nine right collar is not cache aligned")
        for cache, degrees in zip(result, order_groups):
            cache.append(
                {
                    **{
                        key: value
                        for key, value in record.items()
                        if key != "h_derivatives"
                    },
                    "h_derivatives": {
                        str(degree): record["h_derivatives"][str(degree)]
                        for degree in degrees
                    },
                }
            )
    return result


def first_target_tile(records: list[dict]) -> int:
    target = flint.arb(TARGET_T)
    for index, record in enumerate(records):
        left = order4.interval_from_text(record["t_left"])
        right = order4.interval_from_text(record["t_right"])
        if bool(right > target):
            if not bool(left < target):
                raise RuntimeError("compact target does not lie inside its first tile")
            return index
    raise RuntimeError("compact cache does not reach t=5700")


def collar_indices(
    records: list[dict],
    central_left: int,
    central_right: int,
) -> tuple[int, int, dict]:
    lower_t = order4.interval_from_text(records[central_left]["t_left"])
    upper_t = order4.interval_from_text(records[central_right - 1]["t_right"])
    left = central_left
    while left > 0:
        candidate = order4.interval_from_text(records[left]["t_left"])
        if bool(candidate < lower_t - COLLAR_T):
            break
        left -= 1
    right = central_right - 1
    while right + 1 < len(records):
        candidate = order4.interval_from_text(records[right]["t_right"])
        if bool(candidate > upper_t + COLLAR_T):
            break
        right += 1
    outer_lower = order4.interval_from_text(records[left]["t_left"])
    outer_upper = order4.interval_from_text(records[right]["t_right"])
    if not bool(
        outer_lower < lower_t - COLLAR_T
        and outer_upper > upper_t + COLLAR_T
    ):
        raise RuntimeError("compact cache misses the required t+-7 collar")
    return left, right, {
        "left_t_collar_ball": (lower_t - outer_lower)
        .str(40)
        .replace("e", "E"),
        "right_t_collar_ball": (outer_upper - upper_t)
        .str(40)
        .replace("e", "E"),
    }


def derivative_cover(
    caches: list[list[dict]],
    left: int,
    right: int,
) -> tuple[dict[int, flint.arb], dict]:
    derivatives = {}
    order_groups = (tuple(range(2, 9)), (9, 10), (11, 12), (13, 14), (15, 16))
    for records, degrees in zip(caches, order_groups):
        selected = records[left : right + 1]
        for degree in degrees:
            derivatives[degree] = order4.hull(
                [
                    order4.interval_from_text(
                        record["h_derivatives"][str(degree)]
                    )
                    for record in selected
                ]
            )
    return derivatives, {
        "cover_tile_index_first": left,
        "cover_tile_index_last": right,
        "cover_tile_count": right - left + 1,
        "derivative_orders": [2, 16],
    }


def compact_block(
    caches: list[list[dict]],
    central_left: int,
    central_right: int,
) -> dict:
    cover_left, cover_right, collar = collar_indices(
        caches[0], central_left, central_right
    )
    derivatives, diagnostics = derivative_cover(caches, cover_left, cover_right)
    diagnostics.update(collar)
    result = evaluate_nested_curvature_from_h_cover(
        Fraction(caches[0][central_left]["mode_left"]),
        Fraction(caches[0][central_right - 1]["mode_right"]),
        derivatives,
        cover_diagnostics=diagnostics,
    )
    result["central_tile_index_first"] = central_left
    result["central_tile_index_last"] = central_right - 1
    result["central_tile_count"] = central_right - central_left
    return result


def certify_adaptive_block(
    caches: list[list[dict]],
    left: int,
    right: int,
) -> list[dict]:
    result = compact_block(caches, left, right)
    if result.get("passed"):
        return [result]
    if right - left <= 1:
        raise RuntimeError(
            f"order-nine compact failed on irreducible tile {left}: "
            f"{result.get('failure')} {result.get('detail', '')}"
        )
    midpoint = (left + right) // 2
    return certify_adaptive_block(caches, left, midpoint) + certify_adaptive_block(
        caches, midpoint, right
    )


def compact_certificate(caches: list[list[dict]]) -> dict:
    flint.ctx.prec = order4.DEFAULT_PRECISION_BITS
    central_left = first_target_tile(caches[0])
    central_right = order4.mode_index(Fraction(2))
    accepted = []
    cursor = central_left
    while cursor < central_right:
        endpoint = min(cursor + INITIAL_CENTRAL_TILE_COUNT, central_right)
        accepted.extend(certify_adaptive_block(caches, cursor, endpoint))
        cursor = endpoint
    for previous, current in zip(accepted, accepted[1:]):
        if previous["central_mode"][1] != current["central_mode"][0]:
            raise RuntimeError("order-nine compact blocks have a mode gap")
    if accepted[-1]["central_mode"][1] != "2":
        raise RuntimeError("order-nine compact blocks do not reach mode two")
    largest = max(accepted, key=lambda row: Decimal(row["scaled_curvature_upper"]))
    weakest_v = min(accepted, key=lambda row: Decimal(row["V_lower"]))
    start_mode = Fraction(caches[0][central_left]["mode_left"])
    start_t = potential_jet_arb(arb_rational(start_mode), 1)[1]
    return {
        "target_t": TARGET_T,
        "mode_range": [str(start_mode), "2"],
        "start_t_ball": start_t.str(40).replace("e", "E"),
        "end_t_ball": potential_jet_arb(flint.arb(2), 1)[1]
        .str(40)
        .replace("e", "E"),
        "blocks": accepted,
        "block_count": len(accepted),
        "all_blocks_passed": True,
        "largest_scaled_curvature_upper": largest["scaled_curvature_upper"],
        "largest_scaled_mode": largest["central_mode"],
        "weakest_V_lower": weakest_v["V_lower"],
        "weakest_V_mode": weakest_v["central_mode"],
        "theorem": "w_1''(t)<=4200/t^2 for every real 5700<=t<=V'(2)",
    }


def source_contract(collar_path: Path) -> dict:
    paths = (
        order4.DEFAULT_CACHE,
        order6.DEFAULT_EXTENSION_CACHE,
        order7.DEFAULT_EXTENSION_CACHE,
        order8.DEFAULT_CACHE,
        order9.DEFAULT_CACHE,
    )
    return {
        "aligned_caches": [
            {
                "path": path.relative_to(REPO_ROOT).as_posix(),
                "sha256": sha256(path),
                "rows": 107452,
            }
            for path in paths
        ],
        "right_collar": collar_path.relative_to(REPO_ROOT).as_posix(),
        "right_collar_sha256": sha256(collar_path),
    }


def build_artifact(caches: list[list[dict]], collar_path: Path) -> dict:
    compact = compact_certificate(caches)
    rows = [
        CompactRow(
            "co9ncc_01_aligned_h_cache",
            "interval_input",
            "ready_to_apply",
            "Aligned quadrature caches enclose H derivatives through order sixteen.",
            "H^(2),...,H^(16) on strict t+-7 collars",
            "First Newman summand at lambda=-100 only.",
        ),
        CompactRow(
            "co9ncc_02_common_nested_core",
            "exact_interval_algebra",
            "ready_to_apply",
            "Six stable logarithms preserve positive coordinates and certify the curvature ceiling.",
            "B,J,R,S,T,U,V>0 and t^2*w_1''(t)<4200",
            "Outward-rounded common-collar arithmetic.",
        ),
        CompactRow(
            "co9ncc_03_compact_theorem",
            "interval_theorem",
            "ready_to_apply",
            "Adaptive rational mode blocks cover t=5700 through saddle mode two.",
            compact["theorem"],
            "Compact saddle range only.",
            {
                "blocks": compact["block_count"],
                "largest_scaled_upper": compact[
                    "largest_scaled_curvature_upper"
                ],
                "weakest_V_lower": compact["weakest_V_lower"],
            },
        ),
        CompactRow(
            "co9ncc_04_open_handoffs",
            "open_handoff",
            "not_ready_to_apply",
            "Join the localized lower bridge and certify the saddle rays beyond mode two.",
            "cover 1250<=t<=5700 and u>=2",
            "Those ranges are not proved by this artifact.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order9_nested_curvature_compact_certificate",
        "date": "2026-07-14",
        "status": "rigorous order-nine curvature theorem on 5700<=t<=V'(2)",
        "proof_boundary": (
            "This artifact proves only the displayed compact first-summand "
            "curvature range. It does not prove the lower bridge, later saddle "
            "rays, full-kernel transfer, order-nine entry, PF-infinity, or RH."
        ),
        "source_contract": source_contract(collar_path),
        "parameters": {
            "target_t": TARGET_T,
            "collar_t": COLLAR_T,
            "initial_central_tile_count": INITIAL_CENTRAL_TILE_COUNT,
            "precision_bits": order4.DEFAULT_PRECISION_BITS,
            "curvature_constant": CURVATURE_CONSTANT,
        },
        "compact": compact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 3,
            "open_rows": 1,
            "compact_blocks": compact["block_count"],
            "compact_curvature_theorems": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_nested_curvature_compact_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order9_nested_curvature_compact_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    compact = artifact["compact"]
    lines = [
        "# Jensen-Window PF Compound Order-Nine Nested-Curvature Compact Certificate",
        "",
        "Date: 2026-07-14",
        "",
        "Status: rigorous first-summand curvature theorem on `5700<=t<=V'(2)`.",
        "This is not a proof of the lower bridge, later saddle rays, PF-infinity, or RH.",
        "",
        "Aligned H2-H16 interval caches provide every strict `t+-7` collar.",
        "The common nested core verifies `B,J,R,S,T,U,V>0` on every block.",
        "",
        "```text",
        compact["theorem"],
        f"adaptive blocks={compact['block_count']}",
        "largest scaled upper=" + compact["largest_scaled_curvature_upper"] + "<4200",
        "weakest V lower=" + compact["weakest_V_lower"] + ">0",
        "```",
        "",
        "The first tile overlaps the localized target at t=5700 and the final",
        "tile reaches saddle mode u=2 without a mode gap.",
        "The localized `1250<=t<=5700` bridge and the `u>=2` rays remain",
        "separate handoffs.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--right-collar", type=Path, default=DEFAULT_RIGHT_COLLAR)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--overwrite-collar", action="store_true")
    args = parser.parse_args()

    tasks = deterministic_tasks()
    caches = [
        order4.load_cache(order4.DEFAULT_CACHE, tasks),
        order6.load_extension_cache(order6.DEFAULT_EXTENSION_CACHE, tasks),
        order7.load_extension_cache(order7.DEFAULT_EXTENSION_CACHE, tasks),
        order8.load_cache(order8.DEFAULT_CACHE, tasks),
        order9.load_cache(order9.DEFAULT_CACHE, tasks),
    ]
    if any(len(records) != len(tasks) for records in caches):
        raise RuntimeError("order-nine compact source caches are incomplete")
    collar = build_right_collar(
        args.right_collar,
        overwrite=args.overwrite_collar,
    )
    caches = append_right_collar(caches, collar)
    artifact = build_artifact(caches, args.right_collar)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(
        "wrote order-nine nested compact certificate: "
        f"{artifact['summary']['compact_blocks']} blocks, scaled upper "
        f"{artifact['compact']['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
