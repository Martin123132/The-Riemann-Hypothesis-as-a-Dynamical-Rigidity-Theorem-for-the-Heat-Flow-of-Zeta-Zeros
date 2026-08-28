#!/usr/bin/env python3
"""Assemble the cached upper-event-cell value and derivative landscape.

This is a route-selection diagnostic.  It records rigorous local boxes from
the append-only caches, but does not infer interval sign from sampled values.
"""

from __future__ import annotations

import hashlib
import json
from decimal import Decimal
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_local_complete_K_T_scout as complete
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_point_scout as point


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "upper_event_cell_phase_landscape_scout"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
RADIUS = "0.0000001"
PRECISION_BITS = 384
VARIANT = "production"
INTEGER_OFFSETS = tuple(str(index) for index in range(1, 21))
CANDIDATE_TROUGHS = tuple(
    format((Decimal("1.5") * index).normalize(), "f") for index in range(1, 14)
)
CANDIDATE_PEAKS = tuple(
    format(Decimal("0.75") + Decimal("1.5") * index, "f") for index in range(13)
)
DERIVATIVE_OFFSETS = ("5", "10", "15", "20")
DERIVATIVE_RADII = {"5": RADIUS, "10": RADIUS, "15": RADIUS, "20": "0.005"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_hash(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def point_record(rows: dict[str, dict[str, Any]], offset: str) -> dict[str, Any]:
    center = point.center_from_offset(offset)
    key = point.cache_key(center, RADIUS, PRECISION_BITS, VARIANT)
    require(key in rows, f"missing point row: {key}")
    row = rows[key]
    require(row.get("passed") is True, f"failed point row: {key}")
    record = row["record"]
    q_value = record["projection"]["Q_K_minus_T_height_box"]
    require(record["projection"]["strictly_negative_on_box"] is True, f"nonnegative point box: {key}")
    return {
        "offset": offset,
        "center_height": center,
        "radius": RADIUS,
        "Q_K_minus_T": q_value,
        "strictly_negative_on_local_box": True,
        "cache_key": key,
        "canonical_row_sha256": canonical_hash(row),
    }


def derivative_record(rows: dict[str, dict[str, Any]], offset: str) -> dict[str, Any]:
    center = complete.center_from_offset(offset)
    radius = DERIVATIVE_RADII[offset]
    assembly = complete.assemble(rows, center, radius, PRECISION_BITS, VARIANT)
    derivative = assembly["Q_K_minus_T_derivative_component_sum"]
    return {
        "offset": offset,
        "center_height": center,
        "radius": radius,
        "Q_K_minus_T_derivative": derivative,
        "strictly_positive": assembly["strictly_positive_derivative"],
        "strictly_negative": assembly["strictly_negative_derivative"],
        "component_row_hashes": assembly["component_row_hashes"],
    }


def selected_rows(rows: dict[str, dict[str, Any]], offsets: tuple[str, ...]) -> list[dict[str, Any]]:
    return [point_record(rows, offset) for offset in offsets]


def build_artifact() -> dict[str, Any]:
    ctx.prec = 512
    point_rows = point.load_cache()
    component_rows = complete.load_cache()
    integers = selected_rows(point_rows, INTEGER_OFFSETS)
    troughs = selected_rows(point_rows, CANDIDATE_TROUGHS)
    peaks = selected_rows(point_rows, CANDIDATE_PEAKS)
    derivatives = [derivative_record(component_rows, offset) for offset in DERIVATIVE_OFFSETS]

    require(derivatives[0]["strictly_positive"], "offset-5 derivative is not positive")
    require(derivatives[1]["strictly_negative"], "offset-10 derivative is not negative")
    require(derivatives[2]["strictly_negative"], "offset-15 derivative is not negative")
    require(derivatives[3]["strictly_positive"], "offset-20 derivative is not positive")

    worst_peak = max(peaks, key=lambda row: arb(row["Q_K_minus_T"]["upper"]))
    worst_peak_upper = arb(worst_peak["Q_K_minus_T"]["upper"])
    require(worst_peak_upper < 0, "candidate-peak upper envelope reached zero")

    cells = []
    for center_text in CANDIDATE_PEAKS:
        center = Decimal(center_text)
        cells.append(
            {
                "center_offset": center_text,
                "left_offset": format(center - Decimal("0.75"), "f"),
                "right_offset": format(center + Decimal("0.75"), "f"),
            }
        )

    return {
        "kind": "rh_diagnostic_upper_event_cell_phase_landscape",
        "date": "2026-08-28",
        "status": "cache_backed_phase_landscape_and_route_selection_not_interval_theorem",
        "passed": True,
        "precision_bits": PRECISION_BITS,
        "variant": VARIANT,
        "local_box_radius": RADIUS,
        "resource_policy": "cache-only assembly; no new numerical worker",
        "value_landscape": {
            "integer_grid": integers,
            "candidate_trough_grid": troughs,
            "candidate_peak_grid": peaks,
            "all_integer_boxes_strictly_negative": True,
            "all_candidate_trough_boxes_strictly_negative": True,
            "all_candidate_peak_boxes_strictly_negative": True,
            "closest_sampled_candidate_peak_to_zero": worst_peak,
            "candidate_peak_upper_envelope": worst_peak["Q_K_minus_T"]["upper"],
            "observed_dominant_spacing": "approximately 1.5 height units",
            "spacing_is_a_diagnostic_observation_not_a_proved_period": True,
        },
        "derivative_landscape": {
            "rows": derivatives,
            "sign_sequence": ["positive", "negative", "negative", "positive"],
            "single_positive_derivative_chain_across_offsets_5_to_20_rejected": True,
            "continuity_consequences": [
                "at least one stationary point lies strictly between offsets 5 and 10",
                "at least one stationary point lies strictly between offsets 15 and 20",
            ],
            "stationary_points_are_not_claimed_unique": True,
        },
        "route_decision": {
            "selected_next_mechanism": (
                "certified phase-cell envelope with interval root isolation or a local Taylor/Chebyshev model"
            ),
            "candidate_full_cells": cells,
            "candidate_full_cell_union": ["0", "19.5"],
            "remaining_tail_before_existing_offset20_sign_island": ["19.5", "19.93"],
            "rejected_mechanism": "one globally positive derivative chain",
            "why": "the complete derivative is rigorously negative at offsets 10 and 15",
        },
        "sources": {
            point.CACHE.relative_to(ROOT).as_posix(): file_hash(point.CACHE),
            complete.CACHE.relative_to(ROOT).as_posix(): file_hash(complete.CACHE),
            Path(point.__file__).resolve().relative_to(ROOT).as_posix(): file_hash(Path(point.__file__).resolve()),
            Path(complete.__file__).resolve().relative_to(ROOT).as_posix(): file_hash(Path(complete.__file__).resolve()),
            Path(__file__).resolve().relative_to(ROOT).as_posix(): file_hash(Path(__file__).resolve()),
        },
        "proof_boundary": (
            "Every listed local Arb box is rigorous, and the derivative signs rule out a single "
            "positive-derivative bridge.  The sampled candidate peaks are not certified extrema; "
            "their negativity does not fill any interval, count all stationary points, or imply RH."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    values = artifact["value_landscape"]
    derivatives = artifact["derivative_landscape"]["rows"]
    worst = values["closest_sampled_candidate_peak_to_zero"]
    lines = [
        "# Upper-event-cell phase landscape scout",
        "",
        "Date: 2026-08-28",
        "",
        "Status: cache-backed route diagnostic; not an interval theorem or RH proof",
        "",
        "## What the certified local boxes say",
        "",
        "The complete derivative has the sign sequence",
        "",
        "```text",
    ]
    for row in derivatives:
        d = row["Q_K_minus_T_derivative"]
        sign = "positive" if row["strictly_positive"] else "negative"
        lines.append(f"offset {row['offset']}: {sign}, [{d['lower']}, {d['upper']}]")
    lines.extend(
        [
            "```",
            "",
            "Thus a single positive-derivative bridge from offset 20 toward offset 5 is impossible.",
            "Continuity gives at least one stationary point in `(5,10)` and at least one in",
            "`(15,20)`, without proving uniqueness.",
            "",
            "All 20 integer-height boxes, all 13 candidate-trough boxes, and all 13",
            "phase-predicted candidate-peak boxes are strictly negative.  The closest sampled",
            "candidate peak to zero is",
            "",
            "```text",
            f"offset {worst['offset']}: upper <= {worst['Q_K_minus_T']['upper']} < 0.",
            "```",
            "",
            "The observed spacing is approximately `1.5` height units.  This is route-sizing",
            "evidence, not a proved period or a certified list of extrema.",
            "",
            "## Next certificate",
            "",
            "Build a phase-normalized interval model on the 13 candidate cells covering offsets",
            "`[0,19.5]`, isolate every derivative zero or certify a Taylor/Chebyshev upper",
            "envelope, and handle the remaining `[19.5,19.93]` tail before the existing",
            "offset-20 sign island.",
            "",
            "## Proof boundary",
            "",
            artifact["proof_boundary"],
            "",
        ]
    )
    return "\n".join(lines)


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def main() -> int:
    artifact = build_artifact()
    atomic_write(RESULT, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE, render_note(artifact))
    print(
        "completed upper-event-cell phase landscape: "
        f"integer={len(INTEGER_OFFSETS)} troughs={len(CANDIDATE_TROUGHS)} "
        f"peaks={len(CANDIDATE_PEAKS)} worst_peak={artifact['value_landscape']['candidate_peak_upper_envelope']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
