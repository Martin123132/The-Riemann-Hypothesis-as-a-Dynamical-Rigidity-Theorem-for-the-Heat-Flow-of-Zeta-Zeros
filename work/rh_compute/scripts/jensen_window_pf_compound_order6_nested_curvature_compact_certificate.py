#!/usr/bin/env python3
"""Certify the order-six nested curvature from t=321 to saddle mode u=2."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from time import perf_counter


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as order4_compact  # noqa: E402
from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    integrate_h_derivatives,
)
from jensen_window_pf_compound_order6_nested_curvature_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    evaluate_nested_curvature_from_h_cover,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


DEFAULT_EXTENSION_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_nested_curvature_compact_h9_h10_tiles.jsonl"
)
DEFAULT_RIGHT_COLLAR = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_nested_curvature_compact_right_collar.json"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_nested_curvature_compact_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.md"
)
SOURCE_COMPACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_localized_curvature_compact_certificate.json"
)
SOURCE_BRIDGE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_first_summand_curvature_bridge.json"
)
TARGET_T = 321
COLLAR_T = 4
INITIAL_CENTRAL_TILE_COUNT = 3000
DEFAULT_WORKERS = max(1, min(12, (os.cpu_count() or 4) - 1))


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
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_tasks() -> list[tuple[int, Fraction, Fraction]]:
    return [
        (index, left, right)
        for index, (left, right) in enumerate(
            order4_compact.fraction_grid(
                order4_compact.OUTER_START,
                order4_compact.OUTER_END,
                order4_compact.TILE_WIDTH,
            )
        )
    ]


def initialize_worker() -> None:
    flint.ctx.prec = order4_compact.DEFAULT_PRECISION_BITS


def extension_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, left, right = task
    flint.ctx.prec = order4_compact.DEFAULT_PRECISION_BITS
    result = integrate_h_derivatives(
        left,
        right,
        order4_compact.PANELS,
        window_y=order4_compact.WINDOW_Y,
        eighth_envelope_bound=order4_compact.EIGHTH_ENVELOPE,
        max_moment=10,
    )
    if not result.get("passed"):
        return {
            "index": index,
            "mode_left": str(left),
            "mode_right": str(right),
            "passed": False,
            "failure": result.get("failure"),
        }
    return {
        "kind": "order6_compact_h9_h10_tile",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        "passed": True,
        "h_derivatives": {
            "9": order4_compact.interval_text(result["h_derivatives"][9]),
            "10": order4_compact.interval_text(result["h_derivatives"][10]),
        },
    }


def build_right_collar(path: Path, *, overwrite: bool = False) -> dict:
    if path.exists() and not overwrite:
        return load_json(path)
    index = len(deterministic_tasks())
    left = order4_compact.OUTER_END
    right = left + order4_compact.TILE_WIDTH
    flint.ctx.prec = order4_compact.DEFAULT_PRECISION_BITS
    result = integrate_h_derivatives(
        left,
        right,
        order4_compact.PANELS,
        window_y=order4_compact.WINDOW_Y,
        eighth_envelope_bound=order4_compact.EIGHTH_ENVELOPE,
        max_moment=10,
    )
    if not result.get("passed"):
        raise RuntimeError(f"right-collar derivative tile failed: {result}")
    t_left = potential_jet_arb(arb_rational(left), 1)[1]
    t_right = potential_jet_arb(arb_rational(right), 1)[1]
    record = {
        "kind": "order6_compact_full_right_collar_tile",
        "index": index,
        "mode_left": str(left),
        "mode_right": str(right),
        "passed": True,
        "t_left": order4_compact.interval_text(t_left),
        "t_right": order4_compact.interval_text(t_right),
        "h_derivatives": {
            str(order): order4_compact.interval_text(result["h_derivatives"][order])
            for order in range(2, 11)
        },
        "normalizer_lower": result["normalizer_lower"],
        "minimum_tail_slope_lower": result["minimum_tail_slope_lower"],
        "maximum_tail_moment_upper": result["maximum_tail_moment_upper"],
        "maximum_simpson_error_upper": result["maximum_simpson_error_upper"],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return record


def append_right_collar(
    base_records: list[dict], extension_records: list[dict], collar: dict
) -> tuple[list[dict], list[dict]]:
    expected_index = len(base_records)
    if collar.get("index") != expected_index or collar.get("passed") is not True:
        raise RuntimeError("right-collar tile is not aligned with the compact caches")
    base_row = {
        key: value
        for key, value in collar.items()
        if key != "kind"
    }
    base_row["h_derivatives"] = {
        str(order): collar["h_derivatives"][str(order)] for order in range(2, 9)
    }
    extension_row = {
        "kind": "order6_compact_h9_h10_tile",
        "index": collar["index"],
        "mode_left": collar["mode_left"],
        "mode_right": collar["mode_right"],
        "passed": True,
        "h_derivatives": {
            str(order): collar["h_derivatives"][str(order)] for order in (9, 10)
        },
    }
    return [*base_records, base_row], [*extension_records, extension_row]


def load_extension_cache(
    path: Path, tasks: list[tuple[int, Fraction, Fraction]]
) -> list[dict]:
    if not path.exists():
        return []
    records = [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(records) > len(tasks):
        raise RuntimeError("order-six extension cache has too many rows")
    for record, (index, left, right) in zip(records, tasks):
        if (
            record.get("index") != index
            or record.get("mode_left") != str(left)
            or record.get("mode_right") != str(right)
            or record.get("passed") is not True
        ):
            raise RuntimeError(f"order-six extension cache mismatch at tile {index}")
    return records


def build_extension_cache(
    path: Path,
    tasks: list[tuple[int, Fraction, Fraction]],
    *,
    workers: int,
    overwrite: bool,
    max_tiles: int | None,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_extension_cache(path, tasks)
    stop = len(tasks) if max_tiles is None else min(len(tasks), max_tiles)
    remaining = tasks[len(records) : stop]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    if workers == 1:
        iterator = map(extension_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers, initializer=initialize_worker
        )
        iterator = executor.map(extension_task, remaining, chunksize=16)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        f"H9/H10 tile {record.get('index')} failed: "
                        f"{record.get('failure')}"
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 1000 == 0:
                    handle.flush()
                    print(
                        f"order-six compact H9/H10 tiles: {len(records)}/{stop} "
                        f"({perf_counter() - start:.1f}s)"
                    )
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def source_contract(
    extension_path: Path,
    extension_rows: int,
    right_collar_path: Path,
) -> dict:
    compact = load_json(SOURCE_COMPACT)
    bridge = load_json(SOURCE_BRIDGE)
    if compact.get("summary", {}).get("all_blocks_passed") is not True:
        raise RuntimeError("order-four compact H cache is not certified")
    if bridge.get("exact", {}).get("continuous_target") != (
        "p_1''(t)<=200/t^2 for every real t>=321"
    ):
        raise RuntimeError("order-six continuous target changed")
    base_hash = sha256(order4_compact.DEFAULT_CACHE)
    if base_hash != compact.get("cache", {}).get("sha256"):
        raise RuntimeError("order-four compact H cache hash changed")
    return {
        "base_cache": order4_compact.DEFAULT_CACHE.relative_to(REPO_ROOT).as_posix(),
        "base_cache_sha256": base_hash,
        "base_cache_rows": compact["cache"]["row_count"],
        "extension_cache": extension_path.relative_to(REPO_ROOT).as_posix(),
        "extension_cache_sha256": sha256(extension_path),
        "extension_cache_rows": extension_rows,
        "right_collar": right_collar_path.relative_to(REPO_ROOT).as_posix(),
        "right_collar_sha256": sha256(right_collar_path),
        "continuous_target": bridge["exact"]["continuous_target"],
    }


def first_target_tile(records: list[dict]) -> int:
    target = flint.arb(TARGET_T)
    for index, record in enumerate(records):
        left = order4_compact.interval_from_text(record["t_left"])
        right = order4_compact.interval_from_text(record["t_right"])
        if bool(right > target):
            if not bool(left < target):
                raise RuntimeError("compact target tile does not straddle t=321")
            return index
    raise RuntimeError("compact H cache does not reach t=321")


def collar_indices(
    records: list[dict], central_left: int, central_right: int
) -> tuple[int, int, dict]:
    lower_t = order4_compact.interval_from_text(records[central_left]["t_left"])
    upper_t = order4_compact.interval_from_text(records[central_right - 1]["t_right"])
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
        raise RuntimeError("cached H tiles do not contain the required t+-4 collar")
    return left_index, right_index, {
        "left_t_collar_ball": (lower_t - outer_lower).str(40).replace("e", "E"),
        "right_t_collar_ball": (outer_upper - upper_t).str(40).replace("e", "E"),
    }


def derivative_cover(
    base_records: list[dict], extension_records: list[dict], left: int, right: int
) -> tuple[dict[int, flint.arb], dict]:
    derivatives, diagnostics = order4_compact.derivative_cover(
        base_records, left, right
    )
    selected = extension_records[left : right + 1]
    for order in (9, 10):
        derivatives[order] = order4_compact.hull(
            [
                order4_compact.interval_from_text(row["h_derivatives"][str(order)])
                for row in selected
            ]
        )
    diagnostics["extended_derivative_orders"] = [9, 10]
    return derivatives, diagnostics


def compact_block(
    base_records: list[dict],
    extension_records: list[dict],
    central_left: int,
    central_right: int,
) -> dict:
    cover_left, cover_right, collar = collar_indices(
        base_records, central_left, central_right
    )
    derivatives, diagnostics = derivative_cover(
        base_records, extension_records, cover_left, cover_right
    )
    diagnostics.update(collar)
    left = Fraction(base_records[central_left]["mode_left"])
    right = Fraction(base_records[central_right - 1]["mode_right"])
    result = evaluate_nested_curvature_from_h_cover(
        left, right, derivatives, cover_diagnostics=diagnostics
    )
    result["central_tile_index_first"] = central_left
    result["central_tile_index_last"] = central_right - 1
    result["central_tile_count"] = central_right - central_left
    result["cover_tile_count"] = cover_right - cover_left + 1
    return result


def certify_adaptive_block(
    base_records: list[dict], extension_records: list[dict], left: int, right: int
) -> list[dict]:
    result = compact_block(base_records, extension_records, left, right)
    if result.get("passed"):
        return [result]
    if right - left <= 1:
        raise RuntimeError(
            f"order-six nested curvature failed on irreducible tile {left}: "
            f"{result.get('failure')} {result.get('detail', '')}"
        )
    midpoint = (left + right) // 2
    return certify_adaptive_block(
        base_records, extension_records, left, midpoint
    ) + certify_adaptive_block(base_records, extension_records, midpoint, right)


def compact_certificate(
    base_records: list[dict], extension_records: list[dict]
) -> dict:
    flint.ctx.prec = order4_compact.DEFAULT_PRECISION_BITS
    tasks = deterministic_tasks()
    if (
        len(base_records) != len(tasks) + 1
        or len(extension_records) != len(tasks) + 1
    ):
        raise RuntimeError("compact assembly requires complete aligned caches")
    central_left = first_target_tile(base_records)
    central_right = order4_compact.mode_index(Fraction(2))
    accepted: list[dict] = []
    cursor = central_left
    while cursor < central_right:
        endpoint = min(cursor + INITIAL_CENTRAL_TILE_COUNT, central_right)
        accepted.extend(
            certify_adaptive_block(
                base_records, extension_records, cursor, endpoint
            )
        )
        cursor = endpoint

    for previous, current in zip(accepted, accepted[1:]):
        if previous["central_mode"][1] != current["central_mode"][0]:
            raise RuntimeError("order-six compact cover has a mode gap")
    if accepted[-1]["central_mode"][1] != "2":
        raise RuntimeError("order-six compact cover does not reach mode u=2")

    start_mode = Fraction(base_records[central_left]["mode_left"])
    start_t = potential_jet_arb(arb_rational(start_mode), 1)[1]
    end_t = potential_jet_arb(flint.arb(2), 1)[1]
    first_right_t = order4_compact.interval_from_text(
        base_records[central_left]["t_right"]
    )
    if not bool(start_t < TARGET_T < first_right_t):
        raise RuntimeError("first compact block does not cover t=321")

    worst_margin = min(
        range(len(accepted)), key=lambda i: Decimal(accepted[i]["margin_lower"])
    )
    largest_scaled = max(
        range(len(accepted)),
        key=lambda i: Decimal(accepted[i]["scaled_curvature_upper"]),
    )
    weakest_s = min(
        range(len(accepted)), key=lambda i: Decimal(accepted[i]["S_lower"])
    )
    selected_indices = sorted(
        {
            0,
            len(accepted) // 4,
            len(accepted) // 2,
            3 * len(accepted) // 4,
            len(accepted) - 1,
            worst_margin,
            largest_scaled,
            weakest_s,
        }
    )

    def selected_row(index: int) -> dict:
        row = accepted[index]
        return {
            "index": index,
            "central_mode": row["central_mode"],
            "central_t_ball": row["central_t_ball"],
            "central_tile_count": row["central_tile_count"],
            "cover_tile_count": row["cover_tile_count"],
            "J_lower": row["J_lower"],
            "R_lower": row["R_lower"],
            "S_lower": row["S_lower"],
            "p_second_upper": row["p_second_upper"],
            "margin_lower": row["margin_lower"],
            "scaled_curvature_upper": row["scaled_curvature_upper"],
        }

    return {
        "central_start_mode": str(start_mode),
        "central_end_mode": "2",
        "central_start_t_ball": start_t.str(50).replace("e", "E"),
        "central_end_t_ball": end_t.str(50).replace("e", "E"),
        "target_t": TARGET_T,
        "certified_t_range": "321<=t<=V'(2)",
        "accepted": accepted,
        "selected": [selected_row(index) for index in selected_indices],
        "accepted_blocks": len(accepted),
        "maximum_central_tile_count": max(
            row["central_tile_count"] for row in accepted
        ),
        "largest_scaled_index": largest_scaled,
        "largest_scaled_curvature_upper": accepted[largest_scaled][
            "scaled_curvature_upper"
        ],
        "worst_margin_index": worst_margin,
        "worst_margin_lower": accepted[worst_margin]["margin_lower"],
        "weakest_S_index": weakest_s,
        "weakest_S_lower": accepted[weakest_s]["S_lower"],
        "all_blocks_passed": True,
    }


def build_artifact(
    base_records: list[dict],
    extension_records: list[dict],
    extension_path: Path,
    right_collar_path: Path,
) -> dict:
    contract = source_contract(
        extension_path, len(extension_records) - 1, right_collar_path
    )
    compact = compact_certificate(base_records, extension_records)
    rows = [
        CompactRow(
            "co6ncc_01_cache_extension",
            "interval_input",
            "ready_to_apply",
            "The certified compact H cache is extended by aligned ninth- and tenth-derivative enclosures.",
            "H^(r) enclosed tilewise for r=2,...,10",
            "First Newman summand at lambda=-100 only.",
            {
                "base_rows": contract["base_cache_rows"],
                "extension_rows": contract["extension_cache_rows"],
            },
        ),
        CompactRow(
            "co6ncc_02_stable_layers",
            "interval_certificate",
            "ready_to_apply",
            "Every compact block keeps all three nested stable coordinates strictly positive.",
            "J_1(t)>0, R_1(t)>0, S_1(t)>0",
            "Outward-rounded common-derivative hulls with a t+-4 collar.",
            {"weakest_S_lower": compact["weakest_S_lower"]},
        ),
        CompactRow(
            "co6ncc_03_curvature",
            "interval_theorem",
            "ready_to_apply",
            "The order-six first-summand curvature ceiling holds from t=321 to saddle mode two.",
            "p_1''(t)<=200/t^2 for 321<=t<=V'(2)",
            "Finite adaptive Arb cover only.",
            {
                "largest_scaled_upper": compact["largest_scaled_curvature_upper"],
                "worst_margin_lower": compact["worst_margin_lower"],
            },
        ),
        CompactRow(
            "co6ncc_04_ray_handoff",
            "open_handoff",
            "not_ready_to_apply",
            "Continue the same nested curvature theorem from mode u=2 along the unbounded saddle ray.",
            "p_1''(t)<=200/t^2 for u>=2",
            "The finite and asymptotic ray are separate certificates.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order6_nested_curvature_compact_certificate",
        "date": "2026-07-13",
        "status": "rigorous order-six nested first-summand curvature theorem on the compact t>=321 range",
        "proof_boundary": (
            "This artifact proves the continuous first-summand curvature only on "
            "321<=t<=V'(2). It does not prove the u>=2 ray, full-kernel transfer "
            "by itself, order-six entry, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": contract,
        "compact": compact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 3,
            "open_rows": 1,
            "accepted_blocks": compact["accepted_blocks"],
            "all_blocks_passed": True,
            "compact_curvature_theorems": 1,
            "open_ray_targets": 1,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order6_nested_curvature_compact_certificate.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    compact = artifact["compact"]
    contract = artifact["source_contract"]
    lines = [
        "# Jensen-Window PF Compound Order-Six Nested Curvature Compact Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous first-summand order-six curvature theorem on",
        "`321<=t<=V'(2)`. This is not a proof of the remaining saddle ray,",
        "order-six entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_nested_curvature_compact_certificate.py",
        "```",
        "",
        "## Extended Derivative Cache",
        "",
        "The existing 107452-tile compact cache supplies `H^(2),...,H^(8)`.",
        "An aligned extension supplies only the two genuinely new columns",
        "`H^(9),H^(10)`. Both files are hashed in the JSON artifact.",
        "",
        "```text",
        contract["base_cache"],
        contract["extension_cache"],
        contract["right_collar"],
        f"aligned rows={contract['extension_cache_rows']}",
        "```",
        "",
        "## Certified Cover",
        "",
        "Every adaptive block uses a common `t+-4` derivative collar and",
        "outward-rounded stable power-series arithmetic. It proves",
        "",
        "```text",
        "J_1(t)>0, R_1(t)>0, S_1(t)>0,",
        "p_1''(t)<=200/t^2 for 321<=t<=V'(2).",
        "```",
        "",
        f"Accepted blocks: `{compact['accepted_blocks']}`.",
        f"Largest scaled curvature upper: `{compact['largest_scaled_curvature_upper']}`.",
        f"Weakest curvature margin lower: `{compact['worst_margin_lower']}`.",
        f"Weakest third-gap lower: `{compact['weakest_S_lower']}`.",
        "",
        "## Remaining Ray",
        "",
        "The unresolved continuous region begins at saddle mode `u=2` and is",
        "handled separately by finite-ray and asymptotic-ray certificates.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order6_first_summand_curvature_bridge.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extension-cache", type=Path, default=DEFAULT_EXTENSION_CACHE)
    parser.add_argument("--right-collar", type=Path, default=DEFAULT_RIGHT_COLLAR)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite-cache", action="store_true")
    parser.add_argument("--overwrite-right-collar", action="store_true")
    parser.add_argument("--max-tiles", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()

    tasks = deterministic_tasks()
    extension = build_extension_cache(
        args.extension_cache,
        tasks,
        workers=max(1, args.workers),
        overwrite=args.overwrite_cache,
        max_tiles=args.max_tiles,
    )
    print(f"order-six compact H9/H10 cache rows: {len(extension)}/{len(tasks)}")
    if args.cache_only:
        return 0
    if len(extension) != len(tasks):
        raise RuntimeError("complete the extension cache before certification")
    base = order4_compact.load_cache(order4_compact.DEFAULT_CACHE, tasks)
    collar = build_right_collar(
        args.right_collar, overwrite=args.overwrite_right_collar
    )
    base, extension = append_right_collar(base, extension, collar)
    artifact = build_artifact(
        base,
        extension,
        args.extension_cache,
        args.right_collar,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-six nested compact certificate: "
        f"{summary['accepted_blocks']} blocks, "
        f"{summary['compact_curvature_theorems']} compact theorem, "
        f"{summary['open_ray_targets']} open ray"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
