#!/usr/bin/env python3
"""Scout higher-order order-ten curvature blocks at selected real anchors."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
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

import jensen_window_pf_compound_order4_localized_curvature_compact_certificate as compact  # noqa: E402
import jensen_window_pf_compound_order9_localized_lower_bridge_certificate as point_cache  # noqa: E402
import jensen_window_pf_compound_order10_selected_h2_h24_tenth_cache as selected  # noqa: E402
from jensen_window_pf_compound_order10_localized_final_gap_interval_core import (  # noqa: E402
    PRECISION_BITS,
    localized_seventh_formula_continuation_row,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_higher_order_scaling_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order10_higher_order_scaling_scout.md"
)
H_CACHE = selected.DEFAULT_CACHE
H_MANIFEST = selected.DEFAULT_MANIFEST
POINT_CACHE = point_cache.DEFAULT_POINT_CACHE
POINT_MANIFEST = point_cache.DEFAULT_POINT_MANIFEST
WIDTHS = (Fraction(1, 8), Fraction(1, 4), Fraction(1, 2))
POINT_ORDER = 7
REMAINDER_ORDER = 8


@dataclass(frozen=True)
class ScoutRow:
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


def validate_sources() -> dict:
    h_manifest = json.loads(H_MANIFEST.read_text(encoding="utf-8"))
    h_parameters = h_manifest.get("parameters", {})
    h_record = h_manifest.get("cache", {})
    if (
        h_parameters.get("anchors") != [str(value) for value in selected.ANCHORS]
        or h_parameters.get("tile_width") != "1/10"
        or h_parameters.get("rows_per_anchor") != selected.ROWS_PER_ANCHOR
        or h_parameters.get("max_moment") != 24
        or h_record.get("row_count") != len(selected.ANCHORS) * selected.ROWS_PER_ANCHOR
        or h_record.get("all_rows_passed") is not True
        or h_record.get("sha256") != sha256(H_CACHE)
    ):
        raise RuntimeError("selected H2-H24 source changed")
    point_manifest = json.loads(POINT_MANIFEST.read_text(encoding="utf-8"))
    point_record = point_manifest.get("cache", {})
    if (
        point_record.get("row_count") != point_cache.POINT_ROW_COUNT
        or point_record.get("all_rows_passed") is not True
        or point_record.get("sha256") != sha256(POINT_CACHE)
    ):
        raise RuntimeError("exact half-grid point source changed")
    paths = (H_CACHE, H_MANIFEST, POINT_CACHE, POINT_MANIFEST)
    return {
        "selected_anchors": [str(value) for value in selected.ANCHORS],
        "H_rows": h_record["row_count"],
        "H_orders": [2, 24],
        "point_rows": point_record["row_count"],
        "point_orders": [0, 8],
        "precision_load_invariant": (
            "all serialized Arb balls are parsed after flint.ctx.prec=512"
        ),
        "files": [
            {
                "path": path.relative_to(REPO_ROOT).as_posix(),
                "sha256": sha256(path),
            }
            for path in paths
        ],
    }


def load_h_groups() -> dict[Fraction, list[dict]]:
    groups = {anchor: [] for anchor in selected.ANCHORS}
    with H_CACHE.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            index = int(record["index"])
            anchor = selected.ANCHORS[index // selected.ROWS_PER_ANCHOR]
            derivatives = record.get("h_derivatives", {})
            if (
                record.get("passed") is not True
                or set(derivatives) != {str(order) for order in range(2, 25)}
            ):
                raise RuntimeError(f"invalid selected H row {index}")
            groups[anchor].append(
                {
                    "target_t_left": Fraction(record["target_t_left"]),
                    "target_t_right": Fraction(record["target_t_right"]),
                    "H": {
                        order: compact.interval_from_text(derivatives[str(order)])
                        for order in range(2, 25)
                    },
                }
            )
    for anchor, rows in groups.items():
        if len(rows) != selected.ROWS_PER_ANCHOR:
            raise RuntimeError(f"selected collar at {anchor} is incomplete")
        rows.sort(key=lambda row: row["target_t_left"])
        for previous, current in zip(rows, rows[1:]):
            if previous["target_t_right"] != current["target_t_left"]:
                raise RuntimeError(f"selected collar at {anchor} has a gap")
    return groups


def point_source_for(anchor: Fraction) -> dict:
    source = {}
    for shift in range(-8, 9):
        target = anchor + shift
        index = point_cache._aligned_index(  # noqa: SLF001
            target,
            point_cache.POINT_START_T,
            point_cache.POINT_STEP_T,
        )
        loaded, value = point_cache._load_point_row(index)  # noqa: SLF001
        if loaded != target:
            raise RuntimeError(f"point source mismatch at t={target}")
        source[target] = value
    return source


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    sources = validate_sources()
    groups = load_h_groups()
    point_cache.initialize_worker(str(point_cache.DEFAULT_H_CACHE), str(POINT_CACHE))
    blocks = []
    for anchor in selected.ANCHORS:
        points = point_source_for(anchor)
        for width in WIDTHS:
            block = localized_seventh_formula_continuation_row(
                anchor,
                anchor + width,
                groups[anchor],
                point_order=POINT_ORDER,
                remainder_order=REMAINDER_ORDER,
                point_h_source=points,
                require_pass=False,
            )
            blocks.append({**block, "scout_width": str(width)})

    transitions = {}
    for width in WIDTHS:
        matching = [
            row
            for row in blocks
            if row["scout_width"] == str(width) and row["passed"]
        ]
        transitions[str(width)] = matching[0]["anchor"] if matching else None
    passing = [row for row in blocks if row["passed"]]
    failing = [row for row in blocks if not row["passed"]]
    rows = [
        ScoutRow(
            "co10hos_01_sources",
            "interval_input",
            "ready_to_apply",
            "Selected guarded tenth-tile H2-H24 collars and exact half-grid point jets support the scaling test.",
            "15 selected anchors; Taylor order 7; common remainder order 8",
            "Disjoint collars only.",
            sources,
        ),
        ScoutRow(
            "co10hos_02_selected_blocks",
            "rigorous_interval_scout",
            "ready_to_apply",
            "Every reported pass or failure is a rigorous real-interval enclosure on its selected block.",
            "z_1''(t)<=5500/t^2 tested at widths 1/8,1/4,1/2",
            "No interval between selected anchors is covered.",
            {"blocks": blocks},
        ),
        ScoutRow(
            "co10hos_03_transitions",
            "strategy_diagnostic",
            "ready_to_apply",
            "The first selected passing anchors identify economical candidate regime changes.",
            str(transitions),
            "Selected-anchor evidence, not threshold theorems.",
        ),
        ScoutRow(
            "co10hos_04_boundary",
            "proof_boundary",
            "ready_to_apply",
            "The scout does not compose its disjoint blocks into a continuum bridge.",
            "selected blocks != [1251,5700]",
            "A complete adaptive cache or analytic envelope is still required.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order10_higher_order_scaling_scout",
        "date": "2026-07-16",
        "status": "rigorous selected-block order-ten higher-order scaling scout",
        "proof_boundary": (
            "This artifact certifies only the listed disjoint first-summand real-t "
            "blocks. It does not prove a continuous lower bridge, a full-kernel "
            "curvature theorem, the endpoint tail, delayed entry, RH, or Lambda<=0."
        ),
        "sources": sources,
        "parameters": {
            "point_order": POINT_ORDER,
            "remainder_order": REMAINDER_ORDER,
            "widths": [str(width) for width in WIDTHS],
            "curvature_constant": 5500,
            "precision_bits": PRECISION_BITS,
        },
        "blocks": blocks,
        "transitions": transitions,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "selected_anchors": len(selected.ANCHORS),
            "selected_blocks": len(blocks),
            "passing_blocks": len(passing),
            "failing_blocks": len(failing),
            "continuum_bridge_theorems": 0,
            "full_kernel_theorems": 0,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order10_higher_order_scaling_scout.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order10_higher_order_scaling_scout.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Order-Ten Higher-Order Scaling Scout",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous selected real-interval blocks. This is not a proof;",
        "it is not a",
        "continuous bridge, full-kernel theorem, endpoint tail, or RH proof.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order10_higher_order_scaling_scout.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order10_higher_order_scaling_scout.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_higher_order_scaling_scout.py",
        "```",
        "",
        "## Selected Transitions",
        "",
        "```text",
    ]
    for width, anchor in artifact["transitions"].items():
        lines.append(f"width {width}: first selected passing anchor {anchor}")
    lines.extend(
        [
            "```",
            "",
            "These are strategy transitions, not proved threshold locations.",
            "The blocks between selected anchors remain uncovered.",
            "",
            "The scout supports an adaptive design: short blocks near the lower",
            "edge, quarter blocks in the middle transition, and half-unit blocks",
            "from the far lower bridge toward the existing saddle handoff.",
        ]
    )
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
        "wrote order-ten higher-order scaling scout: "
        f"{summary['selected_blocks']} blocks, "
        f"{summary['passing_blocks']} passed, "
        f"{summary['failing_blocks']} failed"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
