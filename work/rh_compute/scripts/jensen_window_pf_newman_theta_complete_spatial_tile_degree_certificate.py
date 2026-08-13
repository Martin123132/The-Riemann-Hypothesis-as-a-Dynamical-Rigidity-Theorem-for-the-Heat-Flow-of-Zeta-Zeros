#!/usr/bin/env python3
"""Certify two adjacent complete-theta positive-time spatial tiles."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from time import perf_counter

for variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(variable, "1")

import psutil


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import jensen_window_pf_newman_theta_complete_spatial_tile_boundary_scout as tile_scout  # noqa: E402
import jensen_window_pf_newman_theta_fourth_summand_compact_contact_degree_certificate as compact  # noqa: E402


continuation = tile_scout.continuation
fifth = tile_scout.fifth
STEM = "jensen_window_pf_newman_theta_complete_spatial_tile_degree_certificate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
DEFAULT_PARTIAL = REPO_ROOT / "work/rh_compute/results" / f"{STEM}_partial.json"
SCOUT_PARENT = REPO_ROOT / "work/rh_compute/results" / f"{tile_scout.STEM}.json"
CENTRAL_PARENT = tile_scout.TAIL_PARENT
INDEX_PARENT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_first_jet_winding_gate.json"
)

DATE = "2026-08-04"
SCHEMA_VERSION = 1
TIME_LOW = Fraction(0)
TIME_HIGH = Fraction(9, 20)
TILES = {
    "left": (Fraction(135), Fraction(271, 2)),
    "right": (Fraction(136), Fraction(273, 2)),
}
PRECISION_BITS = 192
INITIAL_BOUNDARY_SEGMENTS = 8
MAX_SUBDIVISION_DEPTH = 10
CHECKPOINT_EVERY = 8


@dataclass(frozen=True)
class BoundaryCell:
    tile_id: str
    edge: str
    coordinate_low: Fraction
    coordinate_high: Fraction
    direction: int = 1
    depth: int = 0

    @property
    def x_bounds(self) -> tuple[Fraction, Fraction]:
        return TILES[self.tile_id]

    @property
    def coordinate_center(self) -> Fraction:
        return (self.coordinate_low + self.coordinate_high) / 2

    @property
    def coordinate_radius(self) -> Fraction:
        return (self.coordinate_high - self.coordinate_low) / 2

    def parameter_box(self) -> tuple[Fraction, Fraction, Fraction, Fraction]:
        x_low, x_high = self.x_bounds
        if self.edge == "x_low":
            return self.coordinate_center, self.coordinate_radius, x_low, Fraction(0)
        if self.edge == "x_high":
            return self.coordinate_center, self.coordinate_radius, x_high, Fraction(0)
        if self.edge == "t_low":
            return TIME_LOW, Fraction(0), self.coordinate_center, self.coordinate_radius
        if self.edge == "t_high":
            return TIME_HIGH, Fraction(0), self.coordinate_center, self.coordinate_radius
        raise RuntimeError(f"unknown boundary edge: {self.edge}")

    def oriented_start(self) -> Fraction:
        return self.coordinate_low if self.direction > 0 else self.coordinate_high

    def oriented_end(self) -> Fraction:
        return self.coordinate_high if self.direction > 0 else self.coordinate_low

    def point(self, coordinate: Fraction) -> tuple[Fraction, Fraction]:
        x_low, x_high = self.x_bounds
        if self.edge == "x_low":
            return coordinate, x_low
        if self.edge == "x_high":
            return coordinate, x_high
        if self.edge == "t_low":
            return TIME_LOW, coordinate
        if self.edge == "t_high":
            return TIME_HIGH, coordinate
        raise RuntimeError(f"unknown boundary edge: {self.edge}")

    def serialized(self) -> dict:
        return {
            "tile_id": self.tile_id,
            "edge": self.edge,
            "coordinate_low": fraction_text(self.coordinate_low),
            "coordinate_high": fraction_text(self.coordinate_high),
            "direction": self.direction,
            "depth": self.depth,
        }


class ResourcePark(RuntimeError):
    """Raised between completed interval cells after sustained high CPU."""


class CpuMonitor:
    def __init__(self) -> None:
        self.samples: list[float] = []
        self.consecutive_high = 0

    def sample(self) -> None:
        value = float(psutil.cpu_percent(interval=0.5))
        self.samples.append(value)
        self.consecutive_high = self.consecutive_high + 1 if value > 75.0 else 0
        if self.consecutive_high >= 2:
            raise ResourcePark("two consecutive daytime CPU samples exceeded 75%")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fraction_text(value: Fraction) -> str:
    return continuation.fraction_text(value)


def fraction_grid(low: Fraction, high: Fraction, segments: int) -> list[Fraction]:
    return [low + (high - low) * index / segments for index in range(segments + 1)]


def load_tail_radii() -> list:
    payload = json.loads(CENTRAL_PARENT.read_text(encoding="utf-8"))
    require(
        payload["status"]
        == "rigorous complete-theta local positive-time exclusion complete",
        "complete-tail parent is incomplete",
    )
    values = [
        continuation.arb(text).upper()
        for text in payload["positive_moment_budget"]["tail_derivative_bounds"]
    ]
    require(len(values) == 2, "complete-tail derivative budget is incomplete")
    require(all(value.lower() > 0 for value in values), "tail radius is not positive")
    return values


def point_field(time: Fraction, x: Fraction) -> list:
    family = fifth.FifthTransformFamily(
        Fraction(1),
        time,
        time,
        x,
        weighted=True,
    )
    return [family.point_jet(order) for order in range(2)]


def interval_field(cell: BoundaryCell, tail_radii: list) -> tuple[list, dict]:
    time_center, time_radius, x_center, x_radius = cell.parameter_box()
    family = fifth.FifthTransformFamily(
        Fraction(1),
        time_center,
        time_center + time_radius,
        x_center,
        weighted=True,
    )
    complete_boxes: list = []
    five_boxes: list = []
    remainders: list[dict] = []
    for derivative in range(2):
        five_box, remainder = family.transform_box(
            derivative,
            time_radius,
            x_radius,
        )
        complete_box = five_box + continuation.arb(0, tail_radii[derivative])
        five_boxes.append(five_box)
        complete_boxes.append(complete_box)
        remainders.append(remainder)
    return complete_boxes, {
        "five_term_boxes": five_boxes,
        "tail_radii": tail_radii,
        "taylor_remainders": remainders,
    }


def certify_boundary_cell(cell: BoundaryCell, tail_radii: list) -> dict:
    started = perf_counter()
    fields, details = interval_field(cell, tail_radii)
    separations = [value.abs_lower() for value in fields]
    selected = 0 if separations[0] >= separations[1] else 1
    passed = separations[selected] > 0
    tail_upper = details["tail_radii"][selected]
    perturbation_ratio = (
        tail_upper / separations[selected] if passed else continuation.arb(0)
    )
    return {
        "cell": cell.serialized(),
        "passed": bool(passed),
        "separated_by": ("F" if selected == 0 else "F_x") if passed else None,
        "separation_lower": continuation.arb_text(separations[selected]),
        "selected_tail_abs_upper": continuation.arb_text(tail_upper),
        "selected_tail_to_margin_upper": continuation.arb_text(
            perturbation_ratio.upper()
        ),
        "field_boxes": [continuation.ball_endpoints(value) for value in fields],
        "five_term_boxes": [
            continuation.ball_endpoints(value)
            for value in details["five_term_boxes"]
        ],
        "tail_radii": [continuation.arb_text(value) for value in details["tail_radii"]],
        "taylor_remainders": details["taylor_remainders"],
        "elapsed_seconds": perf_counter() - started,
    }


def cell_from_serialized(payload: dict) -> BoundaryCell:
    return BoundaryCell(
        tile_id=payload["tile_id"],
        edge=payload["edge"],
        coordinate_low=Fraction(payload["coordinate_low"]),
        coordinate_high=Fraction(payload["coordinate_high"]),
        direction=int(payload["direction"]),
        depth=int(payload["depth"]),
    )


def edge_coordinate_bounds(tile_id: str, edge: str) -> tuple[Fraction, Fraction, int]:
    x_low, x_high = TILES[tile_id]
    if edge == "x_low":
        return TIME_LOW, TIME_HIGH, 1
    if edge == "t_high":
        return x_low, x_high, 1
    if edge == "x_high":
        return TIME_LOW, TIME_HIGH, -1
    if edge == "t_low":
        return x_low, x_high, -1
    raise RuntimeError(f"unknown edge: {edge}")


def initial_boundary_cells(segments: int) -> list[BoundaryCell]:
    cells: list[BoundaryCell] = []
    for tile_id in TILES:
        for edge in ("x_low", "t_high", "x_high", "t_low"):
            low, high, direction = edge_coordinate_bounds(tile_id, edge)
            grid = fraction_grid(low, high, segments)
            for index in range(segments):
                cells.append(
                    BoundaryCell(
                        tile_id=tile_id,
                        edge=edge,
                        coordinate_low=grid[index],
                        coordinate_high=grid[index + 1],
                        direction=direction,
                    )
                )
    return cells


def split_cell(cell: BoundaryCell) -> list[BoundaryCell]:
    middle = (cell.coordinate_low + cell.coordinate_high) / 2
    return [
        BoundaryCell(
            cell.tile_id,
            cell.edge,
            cell.coordinate_low,
            middle,
            cell.direction,
            cell.depth + 1,
        ),
        BoundaryCell(
            cell.tile_id,
            cell.edge,
            middle,
            cell.coordinate_high,
            cell.direction,
            cell.depth + 1,
        ),
    ]


def partial_payload(
    queue: list[BoundaryCell],
    certified: list[dict],
    evaluated: int,
    segments: int,
) -> dict:
    return {
        "kind": f"{STEM}_partial",
        "schema_version": SCHEMA_VERSION,
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "configuration": {
            "boundary_segments": segments,
            "tiles": {
                tile_id: [fraction_text(bounds[0]), fraction_text(bounds[1])]
                for tile_id, bounds in TILES.items()
            },
            "scout_parent_sha256": sha256_path(SCOUT_PARENT),
            "central_parent_sha256": sha256_path(CENTRAL_PARENT),
        },
        "evaluated_cells": evaluated,
        "queue": [cell.serialized() for cell in queue],
        "certified": certified,
    }


def load_partial(
    path: Path,
    segments: int,
) -> tuple[list[BoundaryCell], list[dict], int]:
    require(path.exists(), "--resume requested but partial checkpoint is missing")
    payload = json.loads(path.read_text(encoding="utf-8"))
    require(payload["kind"] == f"{STEM}_partial", "partial kind mismatch")
    require(payload["schema_version"] == SCHEMA_VERSION, "partial schema mismatch")
    require(
        payload["source_sha256"] == sha256_path(Path(__file__).resolve()),
        "partial source hash mismatch",
    )
    expected = partial_payload([], [], 0, segments)["configuration"]
    require(payload["configuration"] == expected, "partial configuration drifted")
    return (
        [cell_from_serialized(row) for row in payload["queue"]],
        list(payload["certified"]),
        int(payload["evaluated_cells"]),
    )


def minimum_record(records: list[dict]) -> dict:
    require(records, "cannot minimize an empty record list")
    return min(
        records,
        key=lambda row: float(continuation.arb(row["separation_lower"]).lower()),
    )


def maximum_ratio_record(records: list[dict]) -> dict:
    require(records, "cannot maximize an empty record list")
    return max(
        records,
        key=lambda row: float(
            continuation.arb(row["selected_tail_to_margin_upper"]).upper()
        ),
    )


def certify_boundaries(
    segments: int,
    tail_radii: list,
    monitor: CpuMonitor,
    partial_path: Path,
    resume: bool,
    progress: bool,
) -> tuple[dict, list[dict], bool]:
    if resume:
        queue, certified, evaluated = load_partial(partial_path, segments)
    else:
        queue = initial_boundary_cells(segments)
        certified = []
        evaluated = 0
    initial_count = len(TILES) * 4 * segments
    parked = False
    cursor = 0
    while cursor < len(queue):
        cell = queue[cursor]
        cursor += 1
        record = certify_boundary_cell(cell, tail_radii)
        evaluated += 1
        if record["passed"]:
            certified.append(record)
        elif cell.depth < MAX_SUBDIVISION_DEPTH:
            queue.extend(split_cell(cell))
        else:
            raise RuntimeError(
                f"boundary cell unresolved at depth {cell.depth}: {cell.serialized()}"
            )
        if evaluated % CHECKPOINT_EVERY == 0:
            remaining = queue[cursor:]
            continuation.write_json_atomic(
                partial_path,
                partial_payload(remaining, certified, evaluated, segments),
            )
            if progress:
                print(
                    "complete-theta tile boundary: "
                    f"evaluated={evaluated}, certified={len(certified)}, "
                    f"queued={len(remaining)}",
                    flush=True,
                )
            try:
                monitor.sample()
            except ResourcePark:
                parked = True
                queue = remaining
                break
    if not parked:
        queue = []
    continuation.write_json_atomic(
        partial_path,
        partial_payload(queue, certified, evaluated, segments),
    )

    coverage: dict[str, dict] = {}
    for tile_id in TILES:
        coverage[tile_id] = {}
        for edge in ("x_low", "t_high", "x_high", "t_low"):
            rows = [
                row
                for row in certified
                if row["cell"]["tile_id"] == tile_id
                and row["cell"]["edge"] == edge
            ]
            area = sum(
                Fraction(row["cell"]["coordinate_high"])
                - Fraction(row["cell"]["coordinate_low"])
                for row in rows
            )
            low, high, _ = edge_coordinate_bounds(tile_id, edge)
            expected = high - low
            if not parked:
                require(area == expected, f"boundary coverage failed on {tile_id}:{edge}")
            coverage[tile_id][edge] = {
                "certified_length": fraction_text(area),
                "expected_length": fraction_text(expected),
                "passed": area == expected,
            }

    minimum = minimum_record(certified)
    maximum = maximum_ratio_record(certified)
    return (
        {
            "initial_cells": initial_count,
            "evaluated_cells": evaluated,
            "certified_leaf_cells": len(certified),
            "queued_cells": len(queue),
            "resource_parked": parked,
            "maximum_depth": max(int(row["cell"]["depth"]) for row in certified),
            "tile_edge_counts": {
                tile_id: {
                    edge: sum(
                        row["cell"]["tile_id"] == tile_id
                        and row["cell"]["edge"] == edge
                        for row in certified
                    )
                    for edge in ("x_low", "t_high", "x_high", "t_low")
                }
                for tile_id in TILES
            },
            "coverage_audit": coverage,
            "minimum_separation_lower": minimum["separation_lower"],
            "minimum_separation_cell": minimum["cell"],
            "maximum_tail_to_margin_upper": maximum[
                "selected_tail_to_margin_upper"
            ],
            "maximum_ratio_cell": maximum["cell"],
            "all_boundary_cells_nonzero": not parked and not queue,
        },
        certified,
        not parked,
    )


def ordered_tile_records(tile_id: str, records: list[dict]) -> list[dict]:
    ordered: list[dict] = []
    for edge in ("x_low", "t_high", "x_high", "t_low"):
        rows = [
            row
            for row in records
            if row["cell"]["tile_id"] == tile_id and row["cell"]["edge"] == edge
        ]
        direction = edge_coordinate_bounds(tile_id, edge)[2]
        rows.sort(
            key=lambda row: Fraction(row["cell"]["coordinate_low"]),
            reverse=direction < 0,
        )
        low, high, _ = edge_coordinate_bounds(tile_id, edge)
        covered = sum(
            Fraction(row["cell"]["coordinate_high"])
            - Fraction(row["cell"]["coordinate_low"])
            for row in rows
        )
        require(covered == high - low, f"tile edge coverage failed on {tile_id}:{edge}")
        ordered.extend(rows)
    require(ordered, f"tile {tile_id} boundary is empty")
    return ordered


def certify_tile_winding(tile_id: str, records: list[dict]) -> dict:
    ordered = ordered_tile_records(tile_id, records)
    vertices: list[dict] = []
    point_cache: dict[tuple[Fraction, Fraction], list] = {}
    previous_end: tuple[Fraction, Fraction] | None = None
    for index, record in enumerate(ordered):
        cell = cell_from_serialized(record["cell"])
        start = cell.point(cell.oriented_start())
        end = cell.point(cell.oriented_end())
        if previous_end is not None:
            require(start == previous_end, f"boundary discontinuity before {tile_id}:{index}")
        previous_end = end
        if start not in point_cache:
            point_cache[start] = point_field(start[0], start[1])
        if end not in point_cache:
            point_cache[end] = point_field(end[0], end[1])
        field_box = [continuation.arb(item["ball"]) for item in record["field_boxes"]]
        for endpoint in (start, end):
            endpoint_fields = point_cache[endpoint]
            require(
                all(field_box[q].contains(endpoint_fields[q]) for q in range(2)),
                f"complete segment box omits finite endpoint on {tile_id}:{cell.edge}",
            )
        vertices.append(
            {
                "location": {
                    "edge": cell.edge,
                    "time": fraction_text(start[0]),
                    "x": fraction_text(start[1]),
                },
                "fields": point_cache[start],
                "field_balls": [
                    continuation.ball_endpoints(value) for value in point_cache[start]
                ],
                "segment_field_boxes": record["field_boxes"],
                "segment_separated_by": record["separated_by"],
            }
        )
    first = cell_from_serialized(ordered[0]["cell"])
    require(
        previous_end == first.point(first.oriented_start()),
        f"tile {tile_id} boundary does not close",
    )
    winding = compact.classify_polygon_winding(vertices)
    require(winding["winding_t_x"] == 0, f"tile {tile_id} winding is not zero")
    return {
        "tile_id": tile_id,
        "x": [fraction_text(TILES[tile_id][0]), fraction_text(TILES[tile_id][1])],
        "orientation": (
            "positive boundary orientation in coordinate order (t,x): "
            "x=x_low upward in t, t=t_high rightward in x, "
            "x=x_high downward in t, t=t_low leftward in x"
        ),
        "certified_segments": len(ordered),
        "reference_polygon": "finite five-term endpoint values",
        "convex_segment_homotopy": (
            "Each true complete-theta segment image, its complete endpoint values, "
            "the finite five-term endpoint values, and the finite endpoint chord lie "
            "in the same convex complete-field interval rectangle. The selected "
            "component excludes zero, so this is a nonvanishing homotopy."
        ),
        "polygon": winding,
        "vertices": [
            {
                "location": row["location"],
                "field_balls": row["field_balls"],
                "segment_field_boxes": row["segment_field_boxes"],
                "segment_separated_by": row["segment_separated_by"],
            }
            for row in vertices
        ],
    }


def exact_degree_composition(boundary: dict, windings: dict) -> dict:
    index_payload = json.loads(INDEX_PARENT.read_text(encoding="utf-8"))
    central_payload = json.loads(CENTRAL_PARENT.read_text(encoding="utf-8"))
    require("+floor(m/2)" in json.dumps(index_payload, sort_keys=True), "index parent drifted")
    require(boundary["all_boundary_cells_nonzero"] is True, "boundary is incomplete")
    require(
        all(row["polygon"]["winding_t_x"] == 0 for row in windings.values()),
        "an adjacent tile has nonzero winding",
    )
    require(
        "have no common zero" in central_payload["theorem"],
        "central complete-theta theorem is missing",
    )
    return {
        "nontriviality": (
            "All theta components are positive at x=0, so the complete heat solution "
            "is nontrivial."
        ),
        "imported_local_index": (
            "Every multiplicity-m real heat contact has local index -floor(m/2) "
            "in coordinate order (t,x), including degenerate contacts."
        ),
        "tile_exclusion": (
            "Each adjacent tile has nonvanishing boundary and degree zero. Any contact "
            "would contribute a strictly negative integer, so neither tile contains one."
        ),
        "central_composition": (
            "The hash-pinned complete-tail parent excludes contacts on "
            "[0,0.45]x[135.5,136]. The two adjacent tile theorems meet it on their "
            "closed vertical faces, giving the contiguous band [135,136.5]."
        ),
    }


def render_note(payload: dict) -> str:
    boundary = payload["boundary_certificate"]
    lines = [
        "# Complete-Theta Positive-Time Spatial Tile Degree Certificate",
        "",
        f"Date: {DATE}",
        "",
        "Status: rigorous local complete-theta contact exclusion. This is not a proof of RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py --progress",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Theorem",
        "",
        "The complete theta-kernel heat transform and its first x derivative have no",
        "common zero in",
        "",
        "```text",
        "0<=t<=0.45,  135<=x<=136.5.",
        "```",
        "",
        f"The two new boundary covers use `{boundary['certified_leaf_cells']}` cells",
        f"and have minimum separation `{boundary['minimum_separation_lower']}`.",
        "Both endpoint polygons have winding zero. The all-multiplicity same-sign",
        "index theorem excludes every interior contact, and the central parent fills",
        "the interval between the two new tiles.",
        "",
        "## Proof Boundary",
        "",
        payload["proof_boundary"],
        "",
    ]
    return "\n".join(lines)


def build_payload(
    args: argparse.Namespace,
    baseline: list[float],
    priority_lowered: bool,
    monitor: CpuMonitor,
) -> tuple[dict, bool]:
    started = perf_counter()
    tail_radii = load_tail_radii()
    boundary, records, complete = certify_boundaries(
        args.boundary_segments,
        tail_radii,
        monitor,
        args.partial,
        args.resume,
        args.progress,
    )
    if not complete:
        return (
            {
                "kind": STEM,
                "schema_version": SCHEMA_VERSION,
                "date": DATE,
                "status": "parked after resource threshold",
                "boundary_certificate": boundary,
                "resource_policy": {
                    "mode": "daytime",
                    "active_compute_workers": 1,
                    "thread_caps": 1,
                    "below_normal_priority_applied": priority_lowered,
                    "baseline_cpu_percent": baseline,
                    "runtime_cpu_percent": monitor.samples,
                    "resource_parked": True,
                    "elapsed_seconds": perf_counter() - started,
                },
            },
            False,
        )

    windings = {
        tile_id: certify_tile_winding(tile_id, records) for tile_id in TILES
    }
    degree = exact_degree_composition(boundary, windings)
    payload = {
        "kind": STEM,
        "schema_version": SCHEMA_VERSION,
        "date": DATE,
        "status": "rigorous complete-theta spatial tile exclusion certificate complete",
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "parents": {
            "point_scout": {
                "path": str(SCOUT_PARENT.relative_to(REPO_ROOT)).replace("\\", "/"),
                "sha256": sha256_path(SCOUT_PARENT),
            },
            "central_complete_theta": {
                "path": str(CENTRAL_PARENT.relative_to(REPO_ROOT)).replace("\\", "/"),
                "sha256": sha256_path(CENTRAL_PARENT),
            },
            "all_multiplicity_index": {
                "path": str(INDEX_PARENT.relative_to(REPO_ROOT)).replace("\\", "/"),
                "sha256": sha256_path(INDEX_PARENT),
            },
        },
        "resource_policy": {
            "mode": "daytime",
            "active_compute_workers": 1,
            "thread_caps": 1,
            "below_normal_priority_applied": priority_lowered,
            "baseline_cpu_percent": baseline,
            "baseline_mean_percent": sum(baseline) / len(baseline),
            "runtime_cpu_percent": monitor.samples,
            "resource_parked": False,
            "elapsed_seconds": perf_counter() - started,
        },
        "tiles": {
            tile_id: {
                "time": [fraction_text(TIME_LOW), fraction_text(TIME_HIGH)],
                "x": [fraction_text(bounds[0]), fraction_text(bounds[1])],
                "coordinate_orientation": "(t,x)",
            }
            for tile_id, bounds in TILES.items()
        },
        "combined_band": {
            "time": [fraction_text(TIME_LOW), fraction_text(TIME_HIGH)],
            "x": ["135", "136.5"],
            "coordinate_orientation": "(t,x)",
        },
        "family": "complete theta-kernel heat transform",
        "interval_configuration": {
            "precision_bits": PRECISION_BITS,
            "absolute_integration_tolerance": continuation.ABS_TOL,
            "taylor_time_order": continuation.TAYLOR_TIME_ORDER,
            "taylor_x_order": continuation.TAYLOR_X_ORDER,
            "boundary_segments_per_edge": args.boundary_segments,
            "maximum_subdivision_depth": MAX_SUBDIVISION_DEPTH,
            "cutoff_tail_radius_per_five_term_transform": continuation.TAIL_RADIUS,
            "complete_tail_derivative_radii": [
                continuation.arb_text(value) for value in tail_radii
            ],
        },
        "boundary_certificate": boundary,
        "boundary_cells": records,
        "tile_windings": windings,
        "degree_composition": degree,
        "theorem": (
            "The complete theta-kernel heat transform and its first x derivative "
            "have no common zero in [0,0.45]x[135,136.5]."
        ),
        "scaling_frontier": (
            "The reusable tile cover now extends one half-unit on each side of the "
            "central parent. Continue in contiguous rational tiles while developing "
            "an analytic outer-frequency no-contact threshold."
        ),
        "proof_boundary": (
            "This is a rigorous complete-theta contact exclusion only on the declared "
            "positive-time band 0<=t<=0.45 and 135<=x<=136.5. It does not control "
            "other frequencies, negative heat time, all real x, the global "
            "de Bruijn-Newman constant, a degree-uniform Jensen remainder, Lambda<=0, "
            "RH, or a prize-level conclusion."
        ),
    }
    return payload, True


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--partial", type=Path, default=DEFAULT_PARTIAL)
    parser.add_argument("--boundary-segments", type=int, default=INITIAL_BOUNDARY_SEGMENTS)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--baseline-seconds", type=int, default=5)
    parser.add_argument("--baseline-max-percent", type=float, default=60.0)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    require(args.boundary_segments >= 4, "at least four boundary segments are required")
    priority_lowered = continuation.request_below_normal_priority()
    baseline = [float(psutil.cpu_percent(interval=1.0)) for _ in range(args.baseline_seconds)]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "complete-theta tile degree baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, "
        f"mean={baseline_mean:.1f}%",
        flush=True,
    )
    if baseline_mean > args.baseline_max_percent:
        print("complete-theta tile degree run deferred: daytime baseline is already busy")
        return 2
    require(priority_lowered, "could not apply below-normal process priority")
    continuation.flint.ctx.prec = PRECISION_BITS
    monitor = CpuMonitor()
    payload, complete = build_payload(args, baseline, priority_lowered, monitor)
    continuation.write_json_atomic(args.out, payload)
    if not complete:
        print(
            "complete-theta tile degree certificate parked: "
            f"certified={payload['boundary_certificate']['certified_leaf_cells']}, "
            f"queued={payload['boundary_certificate']['queued_cells']}"
        )
        return 3
    continuation.write_text_atomic(args.note, render_note(payload))
    print(
        "built complete-theta spatial tile degree certificate: "
        f"boundary_cells={payload['boundary_certificate']['certified_leaf_cells']}, "
        f"max_depth={payload['boundary_certificate']['maximum_depth']}, "
        f"windings={[row['polygon']['winding_t_x'] for row in payload['tile_windings'].values()]}, "
        "combined_band=[135,136.5], no_contacts=True"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
