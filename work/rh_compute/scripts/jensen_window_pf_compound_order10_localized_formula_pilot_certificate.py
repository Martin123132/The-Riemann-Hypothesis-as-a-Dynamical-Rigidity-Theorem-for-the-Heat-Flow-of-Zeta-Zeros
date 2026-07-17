#!/usr/bin/env python3
"""Certify a local order-ten first-summand curvature interval."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
from jensen_window_pf_compound_order10_localized_final_gap_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    PRECISION_BITS,
    localized_seventh_formula_continuation_row,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_localized_formula_pilot_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order10_localized_formula_pilot_certificate.md"
)
H_SOURCES = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_shifted_jet_pilot_h2_h24_hundredth_tiles.jsonl",
)
H_MANIFESTS = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_shifted_jet_pilot_h2_h24_hundredth_cache.json",
)
POINT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_shifted_point_h0_h8_eighth_grid.jsonl"
)
POINT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_shifted_point_h0_h8_eighth_grid_cache.json"
)
START_T = Fraction(1251)
END_T = Fraction(2503, 2)
BLOCK_WIDTH = Fraction(1, 4)
H_START_T = Fraction(1243)
H_END_T = Fraction(1260)
H_STEP_T = Fraction(1, 100)
H_ORDERS = tuple(range(2, 25))
POINT_START_T = Fraction(1243)
POINT_END_T = Fraction(2519, 2)
POINT_STEP_T = Fraction(1, 8)
POINT_ORDERS = tuple(range(9))
PILOT_POINT_ORDER = 7
PILOT_REMAINDER_ORDER = 8


@dataclass(frozen=True)
class PilotRow:
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
    return path.relative_to(REPO_ROOT).as_posix()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_source_manifests() -> dict:
    expected_h = (("1243", "1260", 1700),)
    files = []
    for cache, manifest_path, expected in zip(H_SOURCES, H_MANIFESTS, expected_h):
        manifest = load_json(manifest_path)
        parameters = manifest.get("parameters", {})
        cache_record = manifest.get("cache", {})
        if (
            parameters.get("start_t") != expected[0]
            or parameters.get("end_t") != expected[1]
            or parameters.get("tile_width_t") != "1/100"
            or parameters.get("max_moment") != 24
            or parameters.get("precision_bits") != 256
            or cache_record.get("row_count") != expected[2]
            or cache_record.get("all_rows_passed") is not True
            or cache_record.get("sha256") != sha256(cache)
        ):
            raise RuntimeError(f"invalid hundredth H source: {manifest_path}")
        files.extend(
            (
                {"path": relative(cache), "sha256": sha256(cache)},
                {"path": relative(manifest_path), "sha256": sha256(manifest_path)},
            )
        )
    point_manifest = load_json(POINT_MANIFEST)
    parameters = point_manifest.get("parameters", {})
    cache_record = point_manifest.get("cache", {})
    if (
        parameters.get("start_t") != str(POINT_START_T)
        or parameters.get("end_t") != str(POINT_END_T)
        or parameters.get("step_t") != str(POINT_STEP_T)
        or parameters.get("precision_bits") != PRECISION_BITS
        or cache_record.get("row_count") != 133
        or cache_record.get("all_rows_passed") is not True
        or cache_record.get("sha256") != sha256(POINT_SOURCE)
    ):
        raise RuntimeError("invalid eighth-grid exact-point source")
    files.extend(
        (
            {"path": relative(POINT_SOURCE), "sha256": sha256(POINT_SOURCE)},
            {"path": relative(POINT_MANIFEST), "sha256": sha256(POINT_MANIFEST)},
        )
    )
    return {
        "precision_load_invariant": (
            "set flint.ctx.prec=512 before parsing every serialized Arb ball"
        ),
        "H_domain": "1243<=t<=1260 on exact 1/100 tiles",
        "H_derivative_orders": [2, 24],
        "H_rows": 1700,
        "point_domain": "1243<=t<=1259.5 on the exact 1/8 grid",
        "point_derivative_orders": [0, 8],
        "point_rows": 133,
        "files": files,
    }


def load_h_rows() -> list[dict]:
    rows = []
    for path in H_SOURCES:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                record = json.loads(line)
                derivatives = record.get("h_derivatives", {})
                if (
                    record.get("passed") is not True
                    or set(derivatives) != {str(order) for order in H_ORDERS}
                ):
                    raise RuntimeError(f"invalid H row in {path}")
                rows.append(
                    {
                        "target_t_left": Fraction(record["target_t_left"]),
                        "target_t_right": Fraction(record["target_t_right"]),
                        "H": {
                            order: compact.interval_from_text(
                                derivatives[str(order)]
                            )
                            for order in H_ORDERS
                        },
                    }
                )
    rows.sort(key=lambda row: row["target_t_left"])
    if len(rows) != 1700:
        raise RuntimeError(f"unexpected H row count: {len(rows)}")
    if rows[0]["target_t_left"] != H_START_T or rows[-1]["target_t_right"] != H_END_T:
        raise RuntimeError("H source endpoints changed")
    for previous, current in zip(rows, rows[1:]):
        if (
            previous["target_t_right"] != current["target_t_left"]
            or previous["target_t_right"] - previous["target_t_left"] != H_STEP_T
        ):
            raise RuntimeError("H source is not a contiguous hundredth grid")
    return rows


def load_point_source() -> dict[Fraction, tuple[list, dict]]:
    points = {}
    with POINT_SOURCE.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            target = Fraction(record["target_t"])
            derivatives = record.get("h_derivatives", {})
            if (
                record.get("passed") is not True
                or set(derivatives) != {str(order) for order in POINT_ORDERS}
            ):
                raise RuntimeError(f"invalid point row at t={target}")
            points[target] = (
                [
                    compact.interval_from_text(derivatives[str(order)])
                    / math.factorial(order)
                    for order in POINT_ORDERS
                ],
                {
                    "target_t": str(target),
                    "mode_bracket": [record["mode_left"], record["mode_right"]],
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
    expected = {
        POINT_START_T + index * POINT_STEP_T
        for index in range(133)
    }
    if set(points) != expected:
        raise RuntimeError("point source is not the expected eighth grid")
    return points


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    source_contract = validate_source_manifests()
    h_rows = load_h_rows()
    point_source = load_point_source()
    blocks = []
    block_specs = (
        (START_T, START_T, START_T + BLOCK_WIDTH),
        (END_T, END_T - BLOCK_WIDTH, END_T),
    )
    for expansion_anchor, left, right in block_specs:
        blocks.append(
            localized_seventh_formula_continuation_row(
                expansion_anchor,
                right,
                h_rows,
                point_order=PILOT_POINT_ORDER,
                remainder_order=PILOT_REMAINDER_ORDER,
                point_h_source=point_source,
                block_left=left,
            )
        )
    if len(blocks) != 2 or any(row.get("passed") is not True for row in blocks):
        raise RuntimeError("order-ten pilot did not produce two passing blocks")
    maximum = max(blocks, key=lambda row: float(row["scaled_curvature_upper"]))
    minimum_margin = min(
        blocks,
        key=lambda row: float(row["curvature_margin_lower"]),
    )
    theorem = (
        "z_1''(t)<=5500/t^2 for every real 1251<=t<=1251.5"
    )
    rows = [
        PilotRow(
            "co10lfpc_01_sources",
            "interval_input",
            "ready_to_apply",
            "Exact hundredth H2-H24 tiles and exact eighth-grid H0-H8 point jets cover the pilot interval.",
            "1700 H tiles; 133 exact point jets; 512-bit deserialization",
            "First Newman summand near t=1251 only.",
            source_contract,
        ),
        PilotRow(
            "co10lfpc_02_formula",
            "exact_inequality",
            "ready_to_apply",
            "The stable-log second derivative is bounded componentwise without composing a pessimistic sixth-order logarithmic remainder.",
            blocks[0]["proof_formula"],
            "Drops only the nonpositive -chi(W)*(W')^2 term.",
        ),
        PilotRow(
            "co10lfpc_03_blocks",
            "interval_certificate",
            "ready_to_apply",
            "Two endpoint-anchored Taylor/remainder blocks pass the 5500 scaled-curvature cap.",
            "[1251,1251.5]=[1251,1251.25] union [1251.25,1251.5]",
            "Continuous real-t pilot interval only.",
            {"blocks": blocks},
        ),
        PilotRow(
            "co10lfpc_04_theorem",
            "interval_theorem",
            "ready_to_apply",
            "The guarded componentwise enclosures prove the displayed local first-summand curvature theorem.",
            theorem,
            "Does not extend beyond t=1251.5 and is not a full-kernel theorem.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order10_localized_formula_pilot_certificate",
        "date": "2026-07-16",
        "status": "rigorous local order-ten first-summand curvature theorem on 1251<=t<=1251.5",
        "proof_boundary": (
            "This artifact proves only the displayed first-summand real-t pilot "
            "interval. It does not prove the complete t>=1251 curvature bound, "
            "the full Newman kernel transfer, the order-ten endpoint tail, "
            "delayed entry, PF-infinity, RH, or Lambda<=0."
        ),
        "theorem": theorem,
        "source_contract": source_contract,
        "parameters": {
            "start_t": str(START_T),
            "end_t": str(END_T),
            "block_width": str(BLOCK_WIDTH),
            "point_order": PILOT_POINT_ORDER,
            "remainder_order": PILOT_REMAINDER_ORDER,
            "curvature_constant": CURVATURE_CONSTANT,
            "precision_bits": PRECISION_BITS,
        },
        "blocks": blocks,
        "diagnostics": {
            "largest_scaled_curvature_anchor": maximum["anchor"],
            "largest_scaled_curvature_upper": maximum[
                "scaled_curvature_upper"
            ],
            "smallest_margin_anchor": minimum_margin["anchor"],
            "smallest_margin_lower": minimum_margin[
                "curvature_margin_lower"
            ],
            "point_scaled_curvature_range": [
                min(row["point_scaled_curvature"] for row in blocks),
                max(row["point_scaled_curvature"] for row in blocks),
            ],
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "accepted_blocks": len(blocks),
            "first_summand_local_curvature_theorems": 1,
            "full_half_line_theorems": 0,
            "full_kernel_theorems": 0,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order10_localized_formula_pilot_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order10_localized_formula_pilot_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    diagnostics = artifact["diagnostics"]
    lines = [
        "# Order-Ten Localized Formula Pilot Certificate",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous first-summand local curvature theorem. This is not a proof",
        "of RH and does not cover the full half-line, full Newman kernel, or order-ten endpoint tail.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order10_localized_formula_pilot_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order10_localized_formula_pilot_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_localized_formula_pilot_certificate.py",
        "```",
        "",
        "## Theorem",
        "",
        "```text",
        artifact["theorem"],
        "```",
        "",
        "The interval is split into two exact width-`1/4` blocks, expanded",
        "rightward from `1251` and leftward from `1251.5`. Every block",
        "uses a hundredth-tile `H2-H24` common remainder and an exact",
        "eighth-grid `H0-H8` point jet parsed at 512-bit precision.",
        "",
        "```text",
        artifact["blocks"][0]["proof_formula"],
        "largest scaled curvature upper: "
        + diagnostics["largest_scaled_curvature_upper"],
        "smallest cap margin lower: " + diagnostics["smallest_margin_lower"],
        "```",
        "",
        "The monolithic sixth-order composed-log remainder is deliberately not",
        "used. Keeping `s`, `w`, and `W` separate preserves the cancellation that",
        "is visible in the exact point jets while retaining a real-interval proof.",
        "",
        "The next task is to extend this economical higher-order component bound",
        "through the adaptive lower regimes and into the saddle handoff.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-ten localized formula pilot: "
        f"{summary['accepted_blocks']} blocks, "
        f"largest scaled upper "
        f"{artifact['diagnostics']['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
