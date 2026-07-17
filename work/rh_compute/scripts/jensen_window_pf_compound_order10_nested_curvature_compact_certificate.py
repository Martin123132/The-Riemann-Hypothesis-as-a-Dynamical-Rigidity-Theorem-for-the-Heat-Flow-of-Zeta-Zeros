#!/usr/bin/env python3
"""Certify the compact order-ten first-summand curvature bridge."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
import math
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

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
import jensen_window_pf_compound_order10_compact_h2_h24_unit_cache as h_source  # noqa: E402
import jensen_window_pf_compound_order10_compact_propagated_point_h0_h7_cache as propagated_source  # noqa: E402
from jensen_window_pf_compound_order10_localized_final_gap_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    PRECISION_BITS,
    localized_seventh_formula_continuation_row,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
    arb_upper_text,
)


H_CACHE = h_source.DEFAULT_CACHE
H_MANIFEST = h_source.DEFAULT_MANIFEST
POINT_CACHE = propagated_source.DEFAULT_CACHE
POINT_MANIFEST = propagated_source.DEFAULT_MANIFEST
ORDER9_THEOREM = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_first_summand_curvature_certificate.json"
)
INTERVAL_CORE = (
    SCRIPT_DIR
    / "jensen_window_pf_compound_order10_localized_final_gap_interval_core.py"
)
DEFAULT_BLOCK_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_nested_curvature_compact_blocks.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_nested_curvature_compact_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order10_nested_curvature_compact_certificate.md"
)

COMPACT_START = Fraction(5700)
UNIT_END = Fraction(10000)
COMPACT_END = Fraction(38020)
UNIT_WIDTH = Fraction(1)
DOUBLE_WIDTH = Fraction(2)
SADDLE_RAY_START = Fraction(2001, 1000)
THEOREM_CURVATURE_CONSTANT = 4200
DEFAULT_WORKERS = max(1, min(4, os.cpu_count() or 1))
GENERATOR_PATH = (
    "work/rh_compute/scripts/"
    "jensen_window_pf_compound_order10_nested_curvature_compact_certificate.py"
)
CHECKER_PATH = (
    "work/rh_compute/scripts/"
    "check_jensen_window_pf_compound_order10_nested_curvature_compact_certificate.py"
)


_H_MAP: dict[Fraction, dict] = {}
_POINT_MAP: dict[Fraction, tuple[list, dict]] = {}
_POINT_PROFILES: dict[Fraction, str] = {}
_SOURCE_CONTRACT_ID = ""


@dataclass(frozen=True)
class CompactLedgerRow:
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


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_point_manifest() -> dict:
    live_contract = propagated_source.validate_sources()
    manifest = load_json(POINT_MANIFEST)
    cache = manifest.get("cache", {})
    expected_parameters = {
        "target_start_t": str(propagated_source.TARGET_START),
        "target_end_t": str(propagated_source.TARGET_END),
        "target_step_t": str(propagated_source.TARGET_STEP),
        "output_orders": [0, propagated_source.OUTPUT_MAX_MOMENT],
        "remainder_order": propagated_source.REMAINDER_MOMENT,
    }
    if (
        manifest.get("kind")
        != "jensen_window_pf_compound_order10_compact_propagated_point_h0_h7_cache"
        or manifest.get("parameters") != expected_parameters
        or manifest.get("source_contract") != live_contract
        or cache.get("path") != relative(POINT_CACHE)
        or cache.get("sha256") != sha256(POINT_CACHE)
        or cache.get("row_count") != 32335
        or cache.get("all_rows_passed") is not True
        or cache.get("coefficient_orders") != [0, 7]
    ):
        raise RuntimeError("invalid propagated H0-H7 point source manifest")
    return {
        "kind": manifest["kind"],
        "domain": "5692<=t<=38026",
        "path": relative(POINT_CACHE),
        "sha256": sha256(POINT_CACHE),
        "manifest_path": relative(POINT_MANIFEST),
        "manifest_sha256": sha256(POINT_MANIFEST),
        "rows": 32335,
        "propagation_source_contract_id": live_contract["id"],
    }


def validate_h_manifest() -> dict:
    manifest = load_json(H_MANIFEST)
    parameters = manifest.get("parameters", {})
    cache = manifest.get("cache", {})
    expected_parameters = {
        "start_t": str(h_source.DEFAULT_START_T),
        "end_t": str(h_source.DEFAULT_END_T),
        "tile_width_t": str(h_source.DEFAULT_TILE_WIDTH_T),
        "initial_mode_bracket": [str(h_source.MODE_LEFT), str(h_source.MODE_RIGHT)],
        "mode_bisections": h_source.MODE_BISECTIONS,
        "max_moment": h_source.MAX_MOMENT,
        "precision_bits": h_source.PRECISION_BITS,
        "panels": h_source.PANELS,
        "window_y": h_source.WINDOW_Y,
        "eighth_envelope": str(compact.EIGHTH_ENVELOPE),
        "row_contract": h_source.ROW_CONTRACT,
    }
    if (
        manifest.get("kind")
        != "jensen_window_pf_compound_order10_compact_h2_h24_unit_cache"
        or parameters != expected_parameters
        or cache.get("path") != relative(H_CACHE)
        or cache.get("sha256") != sha256(H_CACHE)
        or cache.get("row_count") != 32336
        or cache.get("all_rows_passed") is not True
        or cache.get("h_derivative_orders") != [2, 24]
    ):
        raise RuntimeError("invalid compact H2-H24 source manifest")
    return {
        "kind": manifest["kind"],
        "domain": "5692<=t<=38028",
        "path": relative(H_CACHE),
        "sha256": sha256(H_CACHE),
        "manifest_path": relative(H_MANIFEST),
        "manifest_sha256": sha256(H_MANIFEST),
        "rows": 32336,
    }


def validate_order9_theorem() -> dict:
    artifact = load_json(ORDER9_THEOREM)
    if (
        artifact.get("kind")
        != "jensen_window_pf_compound_order9_first_summand_curvature_certificate"
        or artifact.get("status")
        != "rigorous global order-nine first-summand curvature theorem on t>=1250"
        or artifact.get("theorem")
        != "w_1''(t)<=4200/t^2 for every real t>=1250"
    ):
        raise RuntimeError("invalid inherited order-nine curvature theorem")
    return {
        "kind": artifact["kind"],
        "path": relative(ORDER9_THEOREM),
        "sha256": sha256(ORDER9_THEOREM),
        "theorem": artifact["theorem"],
    }


def source_contract() -> dict:
    sources = [
        validate_h_manifest(),
        validate_point_manifest(),
        validate_order9_theorem(),
        {
            "kind": "order10_localized_interval_core",
            "path": relative(INTERVAL_CORE),
            "sha256": sha256(INTERVAL_CORE),
        },
    ]
    geometry = {
        "compact_domain": [str(COMPACT_START), str(COMPACT_END)],
        "unit_blocks": [str(COMPACT_START), str(UNIT_END), str(UNIT_WIDTH)],
        "double_blocks": [str(UNIT_END), str(COMPACT_END), str(DOUBLE_WIDTH)],
        "point_order": 7,
        "remainder_order": 8,
        "curvature_cap": CURVATURE_CONSTANT,
    }
    fingerprint = hashlib.sha256(
        json.dumps(
            {"sources": sources, "geometry": geometry},
            sort_keys=True,
            separators=(",", ":"),
        ).encode("ascii")
    ).hexdigest()
    return {
        "id": fingerprint,
        "precision_load_invariant": (
            f"set flint.ctx.prec={PRECISION_BITS} before parsing every "
            "serialized Arb interval"
        ),
        "geometry": geometry,
        "sources": sources,
    }


def deterministic_blocks() -> list[tuple[int, Fraction, Fraction]]:
    blocks = []
    index = 0
    anchor = COMPACT_START
    while anchor < UNIT_END:
        blocks.append((index, anchor, anchor + UNIT_WIDTH))
        anchor += UNIT_WIDTH
        index += 1
    anchor = UNIT_END
    while anchor < COMPACT_END:
        blocks.append((index, anchor, anchor + DOUBLE_WIDTH))
        anchor += DOUBLE_WIDTH
        index += 1
    if not blocks or blocks[-1][2] != COMPACT_END:
        raise RuntimeError("compact block partition does not close")
    return blocks


def load_h_map(path: Path) -> dict[Fraction, dict]:
    result = {}
    derivative_keys = {str(order) for order in range(2, 25)}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            left = Fraction(record["target_t_left"])
            right = Fraction(record["target_t_right"])
            derivatives = record.get("h_derivatives", {})
            if (
                record.get("kind") != "order10_compact_h2_h24_unit_tile"
                or record.get("contract_id") != h_source.ROW_CONTRACT
                or record.get("passed") is not True
                or right != left + 1
                or set(derivatives) != derivative_keys
            ):
                raise RuntimeError(f"invalid H source row at t={left}")
            result[left] = {
                "target_t_left": left,
                "target_t_right": right,
                "H": {
                    order: compact.interval_from_text(derivatives[str(order)])
                    for order in range(2, 25)
                },
            }
    expected = {Fraction(value) for value in range(5692, 38028)}
    if set(result) != expected:
        raise RuntimeError("compact H source is not the exact unit grid")
    return result


def load_point_map(
    path: Path,
) -> tuple[dict[Fraction, tuple[list, dict]], dict[Fraction, str]]:
    result = {}
    profiles = {}
    keys = {str(order) for order in range(8)}
    propagated_contract_id = load_json(POINT_MANIFEST)["source_contract"]["id"]
    seen = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            target = Fraction(record["target_t"])
            coefficients = record.get("h_taylor_coefficients", {})
            profile_name = record.get("base_profile")
            if (
                record.get("kind")
                != "order10_compact_propagated_point_h0_h7_jet"
                or record.get("source_contract_id") != propagated_contract_id
                or profile_name not in ("medium", "far")
                or record.get("passed") is not True
                or set(coefficients) != keys
            ):
                raise RuntimeError(f"invalid propagated point row at t={target}")
            result[target] = (
                [
                    compact.interval_from_text(coefficients[str(order)])
                    for order in range(8)
                ],
                {
                    "target_t": str(target),
                    "mode_bracket": [
                        record["base_mode_left"],
                        record["base_mode_right"],
                    ],
                    "maximum_panel_error_upper": record[
                        "maximum_panel_error_upper"
                    ],
                    "maximum_tail_moment_upper": record[
                        "maximum_tail_moment_upper"
                    ],
                    "minimum_tail_slope_lower": record[
                        "minimum_tail_slope_lower"
                    ],
                },
            )
            profiles[target] = profile_name
            seen.append(target)
    if (
        len(seen) != 32335
        or seen[0] != Fraction(5692)
        or seen[-1] != Fraction(38026)
        or any(right != left + 1 for left, right in zip(seen, seen[1:]))
    ):
        raise RuntimeError("propagated point source is not contiguous")
    expected = {Fraction(value) for value in range(5692, 38027)}
    if set(result) != expected:
        raise RuntimeError("propagated point source does not cover 5692..38026")
    return result, profiles


def initialize_worker(
    h_path: str,
    point_path: str,
    contract_id: str,
) -> None:
    global _H_MAP, _POINT_MAP, _POINT_PROFILES, _SOURCE_CONTRACT_ID
    flint.ctx.prec = PRECISION_BITS
    _H_MAP = load_h_map(Path(h_path))
    _POINT_MAP, _POINT_PROFILES = load_point_map(Path(point_path))
    _SOURCE_CONTRACT_ID = contract_id


def block_task(task: tuple[int, Fraction, Fraction]) -> dict:
    index, anchor, right = task
    h_rows = [
        _H_MAP[Fraction(target)]
        for target in range(int(anchor - 8), int(right + 8))
    ]
    profiles = sorted(
        {
            _POINT_PROFILES[anchor + shift]
            for shift in range(-8, 9)
        }
    )
    try:
        row = localized_seventh_formula_continuation_row(
            anchor,
            right,
            h_rows,
            point_order=7,
            remainder_order=8,
            point_h_source=_POINT_MAP,
            require_pass=False,
        )
    except Exception as exc:  # pragma: no cover - retained in certificate cache
        return {
            "kind": "order10_compact_curvature_block",
            "source_contract_id": _SOURCE_CONTRACT_ID,
            "index": index,
            "anchor": str(anchor),
            "right": str(right),
            "point_profiles": profiles,
            "passed": False,
            "failure": f"{type(exc).__name__}: {exc}",
        }
    return {
        "kind": "order10_compact_curvature_block",
        "source_contract_id": _SOURCE_CONTRACT_ID,
        "index": index,
        "point_profiles": profiles,
        **row,
    }


def load_block_cache(
    path: Path,
    expected: list[tuple[int, Fraction, Fraction]],
    contract_id: str,
) -> list[dict]:
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
                    f"invalid compact block JSONL row {line_number}"
                ) from exc
    if len(records) > len(expected):
        raise RuntimeError("compact block cache has too many rows")
    for record, task in zip(records, expected):
        index, anchor, right = task
        if (
            record.get("kind") != "order10_compact_curvature_block"
            or record.get("source_contract_id") != contract_id
            or record.get("index") != index
            or record.get("anchor") != str(anchor)
            or record.get("right") != str(right)
            or record.get("passed") is not True
        ):
            raise RuntimeError(f"compact block cache mismatch at row {index}")
    return records


def build_block_cache(
    path: Path,
    expected: list[tuple[int, Fraction, Fraction]],
    contract_id: str,
    *,
    workers: int,
    overwrite: bool,
    max_blocks: int | None,
) -> list[dict]:
    if overwrite and path.exists():
        path.unlink()
    records = load_block_cache(path, expected, contract_id)
    stop = len(expected) if max_blocks is None else min(len(expected), max_blocks)
    remaining = expected[len(records) : stop]
    if not remaining:
        return records
    path.parent.mkdir(parents=True, exist_ok=True)
    initargs = (
        str(H_CACHE),
        str(POINT_CACHE),
        contract_id,
    )
    started = perf_counter()
    if workers == 1:
        initialize_worker(*initargs)
        iterator = map(block_task, remaining)
        executor = None
    else:
        executor = ProcessPoolExecutor(
            max_workers=workers,
            initializer=initialize_worker,
            initargs=initargs,
        )
        iterator = executor.map(block_task, remaining, chunksize=2)
    try:
        with path.open("a", encoding="utf-8") as handle:
            for completed, record in enumerate(iterator, start=1):
                if record.get("passed") is not True:
                    raise RuntimeError(
                        "compact block failed: "
                        + json.dumps(record, sort_keys=True)
                    )
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                records.append(record)
                if completed % 20 == 0:
                    handle.flush()
                if completed % 100 == 0:
                    elapsed = perf_counter() - started
                    rate = completed / elapsed if elapsed else 0.0
                    eta = (len(remaining) - completed) / rate if rate else math.inf
                    print(
                        "order-ten compact curvature blocks: "
                        f"{len(records)}/{stop} ({elapsed:.1f}s; "
                        f"ETA {eta / 3600:.2f}h)",
                        flush=True,
                    )
            handle.flush()
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    return records


def build_artifact(
    records: list[dict],
    block_cache: Path,
    contract: dict,
) -> dict:
    flint.ctx.prec = PRECISION_BITS
    largest = max(records, key=lambda row: float(row["scaled_curvature_upper"]))
    weakest_w = min(records, key=lambda row: float(row["W_lower"]))
    smallest_margin = min(
        records,
        key=lambda row: float(row["curvature_margin_lower"]),
    )
    theorem_margin = Decimal(THEOREM_CURVATURE_CONSTANT) - Decimal(
        largest["scaled_curvature_upper"]
    )
    if theorem_margin <= 0:
        raise RuntimeError(
            "compact scaled curvature does not fit below the 4200 transfer cap"
        )
    transition = potential_jet_arb(arb_rational(SADDLE_RAY_START), 1)[1]
    if not bool(transition < arb_rational(COMPACT_END)):
        raise RuntimeError("u=2.001 transition is not below compact endpoint")
    theorem = "z_1''(t)<=4200/t^2 for every real 5700<=t<=38020"
    rows = [
        CompactLedgerRow(
            "co10nccc_01_sources",
            "interval_input",
            "ready_to_apply",
            "Exact H2-H24 unit tiles and rigorously propagated H0-H7 point jets cover every compact block collar.",
            "H: 5692..38028; propagated points: 5692..38026; point order 7; block remainder order 8",
            "First Newman summand at lambda=-100 only.",
            {
                "source_contract_id": contract["id"],
                "point_source": "step-eight exact H0-H23 jets plus H24 Taylor remainders",
            },
        ),
        CompactLedgerRow(
            "co10nccc_02_partition",
            "interval_partition",
            "ready_to_apply",
            "A mixed exact partition covers the compact real-t interval without gaps or extrapolation.",
            "unit blocks on 5700..10000; width-two blocks on 10000..38020",
            "The partition is specific to the localized seventh-layer certificate.",
            {
                "block_count": len(records),
                "unit_block_count": 4300,
                "double_block_count": 14010,
            },
        ),
        CompactLedgerRow(
            "co10nccc_03_formula",
            "exact_inequality",
            "ready_to_apply",
            "Every block uses the cancellation-preserving stable-log second-derivative inequality.",
            records[0]["proof_formula"],
            "Drops only the nonpositive -chi(W)*(W')^2 term.",
        ),
        CompactLedgerRow(
            "co10nccc_04_blocks",
            "interval_certificate",
            "ready_to_apply",
            "Every compact block has positive W and lies strictly below the scaled-curvature cap.",
            theorem,
            "Continuous real-t compact interval only.",
            {
                "largest_scaled_curvature_upper": largest[
                    "scaled_curvature_upper"
                ],
                "largest_block": [largest["anchor"], largest["right"]],
                "smallest_margin_lower": smallest_margin[
                    "curvature_margin_lower"
                ],
                "theorem_margin_lower": str(theorem_margin),
                "weakest_W_lower": weakest_w["W_lower"],
                "weakest_W_block": [weakest_w["anchor"], weakest_w["right"]],
            },
        ),
        CompactLedgerRow(
            "co10nccc_05_transition",
            "exact_overlap",
            "ready_to_apply",
            "The compact endpoint reaches beyond the proved finite saddle ray start u=2.001.",
            "V'(2.001)<38020",
            "This is an overlap statement, not a full-half-line composition.",
            {"V_prime_2_001_upper": arb_upper_text(transition)},
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order10_nested_curvature_compact_certificate",
        "date": "2026-07-16",
        "status": "rigorous order-ten first-summand curvature theorem on 5700<=t<=38020",
        "proof_boundary": (
            "This certificate concerns the first Newman summand at lambda=-100. "
            "It does not yet transfer the theorem to the full Newman kernel and "
            "does not prove RH, the Jensen hierarchy, or Lambda<=0."
        ),
        "theorem": theorem,
        "source_contract": contract,
        "block_cache": {
            "path": relative(block_cache),
            "sha256": sha256(block_cache),
            "row_count": len(records),
            "all_rows_passed": True,
            "source_contract_id": contract["id"],
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "blocks": len(records),
            "unit_blocks": 4300,
            "double_blocks": 14010,
            "largest_scaled_curvature_upper": largest[
                "scaled_curvature_upper"
            ],
            "smallest_curvature_margin_lower": smallest_margin[
                "curvature_margin_lower"
            ],
            "theorem_margin_lower": str(theorem_margin),
            "weakest_W_lower": weakest_w["W_lower"],
            "saddle_transition_upper": arb_upper_text(transition),
            "compact_first_summand_theorems": 1,
            "full_half_line_theorems": 0,
            "full_kernel_theorems": 0,
            "rh_claims": 0,
        },
        "generator": GENERATOR_PATH,
        "checker": CHECKER_PATH,
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    lines = [
        "# Order-ten compact first-summand curvature certificate",
        "",
        "Date: 2026-07-17",
        "",
        f"Status: **{artifact['status']}**. This certificate is not a proof of PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "## Theorem",
        "",
        f"`{artifact['theorem']}`.",
        "",
        "The exact partition uses unit blocks from `5700` through `10000` and width-two blocks from `10000` through `38020`.",
        f"All `{summary['blocks']}` blocks pass; the largest scaled upper bound is `{summary['largest_scaled_curvature_upper']}`.",
        f"The weakest proved `W` floor is `{summary['weakest_W_lower']}`.",
        "",
        "## Transition",
        "",
        f"The rigorous upper enclosure `V'(2.001) <= {summary['saddle_transition_upper']}` is below `38020`, so this compact theorem overlaps the finite saddle ray.",
        "",
        "## Boundary",
        "",
        artifact["proof_boundary"],
        "",
        "## Reproduce",
        "",
        "```powershell",
        f"python {GENERATOR_PATH}",
        f"python {CHECKER_PATH}",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--block-cache", type=Path, default=DEFAULT_BLOCK_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-blocks", type=int)
    parser.add_argument("--cache-only", action="store_true")
    args = parser.parse_args()

    contract = source_contract()
    expected = deterministic_blocks()
    records = build_block_cache(
        args.block_cache,
        expected,
        contract["id"],
        workers=max(1, args.workers),
        overwrite=args.overwrite,
        max_blocks=args.max_blocks,
    )
    print(f"order-ten compact curvature blocks: {len(records)}/{len(expected)}")
    if args.cache_only:
        return 0
    if len(records) != len(expected):
        raise RuntimeError("complete the compact block cache before promotion")
    artifact = build_artifact(records, args.block_cache, contract)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    print(
        "wrote order-ten compact curvature theorem: "
        f"{len(records)} blocks, largest scaled upper "
        f"{artifact['summary']['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
