#!/usr/bin/env python3
"""Certify order-eight first-summand curvature from t=999 to saddle mode u=2."""

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

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as order4_compact  # noqa: E402
import jensen_window_pf_compound_order6_nested_curvature_compact_certificate as order6_compact  # noqa: E402
import jensen_window_pf_compound_order7_nested_curvature_compact_certificate as order7_compact  # noqa: E402
import jensen_window_pf_compound_order8_nested_curvature_compact_h13_h14_cache as h13_h14_cache  # noqa: E402
from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    integrate_h_derivatives,
)
from jensen_window_pf_compound_order8_nested_curvature_interval_core import (  # noqa: E402
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
    "jensen_window_pf_compound_order8_nested_curvature_compact_right_collar.json"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_nested_curvature_compact_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order8_nested_curvature_compact_certificate.md"
)
SOURCE_BRIDGE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate.json"
)
TARGET_T = 999
COLLAR_T = 6
INITIAL_CENTRAL_TILE_COUNT = 1000
RIGHT_COLLAR_TILE_COUNT = 2


@dataclass(frozen=True)
class CompactRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_tasks() -> list[tuple[int, Fraction, Fraction]]:
    return order7_compact.deterministic_tasks()


def build_right_collar(path: Path, *, overwrite: bool = False) -> dict:
    if path.exists() and not overwrite:
        return load_json(path)
    tasks = deterministic_tasks()
    flint.ctx.prec = order4_compact.DEFAULT_PRECISION_BITS
    records = []
    for offset in range(RIGHT_COLLAR_TILE_COUNT):
        index = len(tasks) + offset
        left = order4_compact.OUTER_END + offset * order4_compact.TILE_WIDTH
        right = left + order4_compact.TILE_WIDTH
        result = integrate_h_derivatives(
            left,
            right,
            order4_compact.PANELS,
            window_y=order4_compact.WINDOW_Y,
            eighth_envelope_bound=order4_compact.EIGHTH_ENVELOPE,
            max_moment=14,
        )
        if not result.get("passed"):
            raise RuntimeError(f"order-eight right collar failed: {result}")
        records.append(
            {
                "kind": "order8_compact_full_right_collar_tile",
                "index": index,
                "mode_left": str(left),
                "mode_right": str(right),
                "passed": True,
                "t_left": order4_compact.interval_text(
                    potential_jet_arb(arb_rational(left), 1)[1]
                ),
                "t_right": order4_compact.interval_text(
                    potential_jet_arb(arb_rational(right), 1)[1]
                ),
                "h_derivatives": {
                    str(order): order4_compact.interval_text(
                        result["h_derivatives"][order]
                    )
                    for order in range(2, 15)
                },
                "normalizer_lower": result["normalizer_lower"],
                "minimum_tail_slope_lower": result["minimum_tail_slope_lower"],
                "maximum_tail_moment_upper": result["maximum_tail_moment_upper"],
                "maximum_simpson_error_upper": result["maximum_simpson_error_upper"],
            }
        )
    artifact = {
        "kind": "order8_compact_full_right_collar",
        "date": "2026-07-13",
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
    base: list[dict],
    high: list[dict],
    top: list[dict],
    ultra: list[dict],
    collar: dict,
) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    expected = len(base)
    records = collar.get("rows", [])
    if (
        collar.get("tile_count") != RIGHT_COLLAR_TILE_COUNT
        or collar.get("all_rows_passed") is not True
        or len(records) != RIGHT_COLLAR_TILE_COUNT
    ):
        raise RuntimeError("order-eight right collar is not cache-aligned")

    def derivative_row(orders: tuple[int, ...], kind: str) -> dict:
        row = {
            key: value
            for key, value in collar.items()
            if key not in {"kind", "h_derivatives"}
        }
        row["kind"] = kind
        row["h_derivatives"] = {
            str(order): collar["h_derivatives"][str(order)]
            for order in orders
        }
        return row

    result = (list(base), list(high), list(top), list(ultra))
    for offset, record in enumerate(records):
        if record.get("index") != expected + offset or record.get("passed") is not True:
            raise RuntimeError("order-eight right collar tile is not aligned")
        collar = record
        result[0].append(derivative_row(tuple(range(2, 9)), "order4_compact_h_tile"))
        result[1].append(derivative_row((9, 10), "order6_compact_h9_h10_tile"))
        result[2].append(derivative_row((11, 12), "order7_compact_h11_h12_tile"))
        result[3].append(derivative_row((13, 14), "order8_compact_h13_h14_tile"))
    return result


def first_target_tile(records: list[dict]) -> int:
    target = flint.arb(TARGET_T)
    for index, record in enumerate(records):
        left = order4_compact.interval_from_text(record["t_left"])
        right = order4_compact.interval_from_text(record["t_right"])
        if bool(right > target):
            if not bool(left < target):
                raise RuntimeError("compact target tile does not straddle t=999")
            return index
    raise RuntimeError("compact H cache does not reach t=999")


def collar_indices(
    records: list[dict],
    central_left: int,
    central_right: int,
) -> tuple[int, int, dict]:
    lower_t = order4_compact.interval_from_text(records[central_left]["t_left"])
    upper_t = order4_compact.interval_from_text(
        records[central_right - 1]["t_right"]
    )
    left_index = central_left
    while left_index > 0:
        candidate = order4_compact.interval_from_text(records[left_index]["t_left"])
        if bool(candidate < lower_t - COLLAR_T):
            break
        left_index -= 1
    right_index = central_right - 1
    while right_index + 1 < len(records):
        candidate = order4_compact.interval_from_text(records[right_index]["t_right"])
        if bool(candidate > upper_t + COLLAR_T):
            break
        right_index += 1
    outer_lower = order4_compact.interval_from_text(records[left_index]["t_left"])
    outer_upper = order4_compact.interval_from_text(records[right_index]["t_right"])
    if not bool(
        outer_lower < lower_t - COLLAR_T
        and outer_upper > upper_t + COLLAR_T
    ):
        raise RuntimeError("cached H tiles do not contain the required t+-6 collar")
    return left_index, right_index, {
        "left_t_collar_ball": (lower_t - outer_lower).str(40).replace("e", "E"),
        "right_t_collar_ball": (outer_upper - upper_t).str(40).replace("e", "E"),
    }


def derivative_cover(
    base: list[dict],
    high: list[dict],
    top: list[dict],
    ultra: list[dict],
    left: int,
    right: int,
) -> tuple[dict[int, flint.arb], dict]:
    derivatives, diagnostics = order4_compact.derivative_cover(base, left, right)
    for records, orders in (
        (high, (9, 10)),
        (top, (11, 12)),
        (ultra, (13, 14)),
    ):
        selected = records[left : right + 1]
        for order in orders:
            derivatives[order] = order4_compact.hull(
                [
                    order4_compact.interval_from_text(
                        row["h_derivatives"][str(order)]
                    )
                    for row in selected
                ]
            )
    diagnostics["extended_derivative_orders"] = [9, 10, 11, 12, 13, 14]
    return derivatives, diagnostics


def compact_block(
    base: list[dict],
    high: list[dict],
    top: list[dict],
    ultra: list[dict],
    central_left: int,
    central_right: int,
) -> dict:
    cover_left, cover_right, collar = collar_indices(
        base, central_left, central_right
    )
    derivatives, diagnostics = derivative_cover(
        base, high, top, ultra, cover_left, cover_right
    )
    diagnostics.update(collar)
    left = Fraction(base[central_left]["mode_left"])
    right = Fraction(base[central_right - 1]["mode_right"])
    result = evaluate_nested_curvature_from_h_cover(
        left,
        right,
        derivatives,
        cover_diagnostics=diagnostics,
    )
    result["central_tile_index_first"] = central_left
    result["central_tile_index_last"] = central_right - 1
    result["central_tile_count"] = central_right - central_left
    result["cover_tile_count"] = cover_right - cover_left + 1
    return result


def certify_adaptive_block(
    base: list[dict],
    high: list[dict],
    top: list[dict],
    ultra: list[dict],
    left: int,
    right: int,
) -> list[dict]:
    result = compact_block(base, high, top, ultra, left, right)
    if result.get("passed"):
        return [result]
    if right - left <= 1:
        raise RuntimeError(
            f"order-eight compact failed on irreducible tile {left}: "
            f"{result.get('failure')} {result.get('detail', '')}"
        )
    midpoint = (left + right) // 2
    return certify_adaptive_block(
        base, high, top, ultra, left, midpoint
    ) + certify_adaptive_block(base, high, top, ultra, midpoint, right)


def source_contract(extension_path: Path, collar_path: Path) -> dict:
    bridge = load_json(SOURCE_BRIDGE)
    if bridge.get("summary", {}).get("passed_blocks") != 185:
        raise RuntimeError("order-eight shifted-jet bridge is not certified")
    manifest = load_json(h13_h14_cache.DEFAULT_MANIFEST)
    if manifest.get("cache", {}).get("row_count") != 107452:
        raise RuntimeError("order-eight H13/H14 manifest is incomplete")
    return {
        "base_cache": order4_compact.DEFAULT_CACHE.relative_to(REPO_ROOT).as_posix(),
        "base_cache_sha256": sha256(order4_compact.DEFAULT_CACHE),
        "base_cache_rows": 107452,
        "h9_h10_cache": order6_compact.DEFAULT_EXTENSION_CACHE.relative_to(REPO_ROOT).as_posix(),
        "h9_h10_cache_sha256": sha256(order6_compact.DEFAULT_EXTENSION_CACHE),
        "h9_h10_cache_rows": 107452,
        "h11_h12_cache": order7_compact.DEFAULT_EXTENSION_CACHE.relative_to(REPO_ROOT).as_posix(),
        "h11_h12_cache_sha256": sha256(order7_compact.DEFAULT_EXTENSION_CACHE),
        "h11_h12_cache_rows": 107452,
        "h13_h14_cache": extension_path.relative_to(REPO_ROOT).as_posix(),
        "h13_h14_cache_sha256": sha256(extension_path),
        "h13_h14_cache_rows": 107452,
        "right_collar": collar_path.relative_to(REPO_ROOT).as_posix(),
        "right_collar_sha256": sha256(collar_path),
        "shifted_bridge": SOURCE_BRIDGE.relative_to(REPO_ROOT).as_posix(),
        "shifted_bridge_sha256": sha256(SOURCE_BRIDGE),
    }


def compact_certificate(
    base: list[dict],
    high: list[dict],
    top: list[dict],
    ultra: list[dict],
) -> dict:
    flint.ctx.prec = order4_compact.DEFAULT_PRECISION_BITS
    central_left = first_target_tile(base)
    central_right = order4_compact.mode_index(Fraction(2))
    accepted = []
    cursor = central_left
    while cursor < central_right:
        endpoint = min(cursor + INITIAL_CENTRAL_TILE_COUNT, central_right)
        accepted.extend(
            certify_adaptive_block(base, high, top, ultra, cursor, endpoint)
        )
        cursor = endpoint
    for previous, current in zip(accepted, accepted[1:]):
        if previous["central_mode"][1] != current["central_mode"][0]:
            raise RuntimeError("order-eight compact cover has a mode gap")
    if accepted[-1]["central_mode"][1] != "2":
        raise RuntimeError("order-eight compact cover does not reach mode u=2")

    start_mode = Fraction(base[central_left]["mode_left"])
    start_t = potential_jet_arb(arb_rational(start_mode), 1)[1]
    first_right_t = order4_compact.interval_from_text(base[central_left]["t_right"])
    if not bool(start_t < TARGET_T < first_right_t):
        raise RuntimeError("first compact block does not cover t=999")

    largest = max(accepted, key=lambda row: Decimal(row["scaled_curvature_upper"]))
    weakest_u = min(accepted, key=lambda row: Decimal(row["U_lower"]))
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
        "weakest_U_lower": weakest_u["U_lower"],
        "weakest_U_mode": weakest_u["central_mode"],
        "theorem": "s_1''(t)<=4000/t^2 for every real 999<=t<=V'(2)",
    }


def build_artifact(
    base: list[dict],
    high: list[dict],
    top: list[dict],
    ultra: list[dict],
    extension_path: Path,
    collar_path: Path,
) -> dict:
    compact = compact_certificate(base, high, top, ultra)
    rows = [
        CompactRow(
            "co8ncc_01_aligned_h_cache",
            "interval_input",
            "ready_to_apply",
            "Aligned quadrature caches enclose H derivatives through order fourteen on the compact mode range.",
            "H^(2),...,H^(14) on a strict t+-6 collar",
            "First Newman summand at lambda=-100 only.",
            {"base_rows": 107452, "extension_rows": 107452},
        ),
        CompactRow(
            "co8ncc_02_common_nested_core",
            "exact_interval_algebra",
            "ready_to_apply",
            "Five cancellation-preserving stable logarithms convert every common derivative cover into B,J,R,S,T,U and s'' intervals.",
            "B,J,R,S,T,U>0 and t^2 s_1''(t)<4000",
            "Outward-rounded common-collar arithmetic; no point sampling.",
        ),
        CompactRow(
            "co8ncc_03_compact_theorem",
            "interval_theorem",
            "ready_to_apply",
            "Adaptive rational mode blocks prove the order-eight first-summand curvature ceiling from t=999 through saddle mode two.",
            compact["theorem"],
            "Compact saddle range only.",
            {
                "blocks": compact["block_count"],
                "largest_scaled_upper": compact["largest_scaled_curvature_upper"],
                "weakest_U_lower": compact["weakest_U_lower"],
            },
        ),
        CompactRow(
            "co8ncc_04_finite_ray_handoff",
            "open_handoff",
            "not_ready_to_apply",
            "Continue the same order-eight curvature theorem on the finite saddle ray beyond mode two.",
            "s_1''(t)<=4000/t^2 for 2<=u<=20",
            "The finite and asymptotic saddle rays remain separate tasks.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order8_nested_curvature_compact_certificate",
        "date": "2026-07-13",
        "status": "rigorous order-eight nested first-summand curvature theorem on 999<=t<=V'(2)",
        "proof_boundary": (
            "This artifact proves the first-summand curvature only from t=999 "
            "through saddle mode u=2. It does not prove the finite or asymptotic "
            "saddle rays, order-eight entry, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": source_contract(extension_path, collar_path),
        "parameters": {
            "target_t": TARGET_T,
            "collar_t": COLLAR_T,
            "initial_central_tile_count": INITIAL_CENTRAL_TILE_COUNT,
            "precision_bits": order4_compact.DEFAULT_PRECISION_BITS,
            "panels": order4_compact.PANELS,
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
            "open_finite_rays": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order8_nested_curvature_compact_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order8_nested_curvature_compact_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    compact = artifact["compact"]
    lines = [
        "# Jensen-Window PF Compound Order-Eight Nested-Curvature Compact Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous first-summand curvature theorem on `999<=t<=V'(2)`.",
        "This is not a proof of the finite/asymptotic saddle rays, order-eight",
        "entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order8_nested_curvature_compact_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order8_nested_curvature_compact_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_nested_curvature_compact_certificate.py",
        "```",
        "",
        "## Common-Collar Core",
        "",
        "Four aligned 107,452-row caches and one right collar enclose H2-H14",
        "on a strict t+-6 neighborhood. Five stable logarithms are evaluated",
        "with outward-rounded common-collar arithmetic, and every accepted",
        "block proves B,J,R,S,T,U>0 before bounding s_1''.",
        "",
        "## Compact Theorem",
        "",
        "```text",
        compact["theorem"] + ",",
        f"all {compact['block_count']} adaptive blocks pass,",
        "largest scaled upper=" + compact["largest_scaled_curvature_upper"] + "<4000,",
        "weakest U lower=" + compact["weakest_U_lower"] + ">0.",
        "```",
        "",
        "The first block overlaps the shifted-jet theorem at t=999. The final",
        "block reaches saddle mode u=2 with no mode gap.",
        "",
        "## Remaining Ray",
        "",
        "```text",
        "prove s_1''(t)<=4000/t^2 for every saddle mode 2<=u<=20,",
        "then prove the asymptotic ray u>=20.",
        "```",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate.md",
        "outputs/jensen_window_pf_compound_order8_first_summand_curvature_bridge.md",
        "outputs/formal_core.md",
        "```",
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
    base = order4_compact.load_cache(order4_compact.DEFAULT_CACHE, tasks)
    high = order6_compact.load_extension_cache(
        order6_compact.DEFAULT_EXTENSION_CACHE, tasks
    )
    top = order7_compact.load_extension_cache(
        order7_compact.DEFAULT_EXTENSION_CACHE, tasks
    )
    ultra = h13_h14_cache.load_cache(h13_h14_cache.DEFAULT_CACHE, tasks)
    if any(len(records) != len(tasks) for records in (base, high, top, ultra)):
        raise RuntimeError("compact source caches are incomplete")
    collar = build_right_collar(
        args.right_collar, overwrite=args.overwrite_collar
    )
    base, high, top, ultra = append_right_collar(
        base, high, top, ultra, collar
    )
    artifact = build_artifact(
        base,
        high,
        top,
        ultra,
        h13_h14_cache.DEFAULT_CACHE,
        args.right_collar,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "wrote order-eight nested compact certificate: "
        f"{artifact['summary']['compact_blocks']} blocks, scaled upper "
        f"{artifact['compact']['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
