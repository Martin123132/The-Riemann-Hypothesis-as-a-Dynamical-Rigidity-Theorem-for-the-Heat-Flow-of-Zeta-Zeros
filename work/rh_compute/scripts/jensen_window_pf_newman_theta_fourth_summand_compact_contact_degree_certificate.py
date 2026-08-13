#!/usr/bin/env python3
"""Certify a compact degree count for the four-summand contact branch."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
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

import jensen_window_pf_newman_theta_fourth_summand_tangency_interval_continuation_certificate as continuation  # noqa: E402


STEM = "jensen_window_pf_newman_theta_fourth_summand_compact_contact_degree_certificate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
DEFAULT_PARTIAL = REPO_ROOT / "work/rh_compute/results" / f"{STEM}_partial.json"
INTERVAL_PARENT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_fourth_summand_"
    "tangency_interval_continuation_certificate.json"
)
INDEX_PARENT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_first_jet_winding_gate.json"
)

SCHEMA_VERSION = 1
DATE = "2026-08-03"

LAMBDA_LOW = Fraction(0)
LAMBDA_HIGH = Fraction(11, 50)
TIME_LOW = Fraction(-3, 50)
TIME_HIGH = Fraction(9, 20)
X_LOW = Fraction(271, 2)
X_HIGH = Fraction(136)
PRECISION_BITS = 192
INITIAL_LAMBDA_STEP = Fraction(1, 50)
INITIAL_BOUNDARY_SEGMENTS = 8
MAX_SUBDIVISION_DEPTH = 8
CHECKPOINT_EVERY = 20


@dataclass(frozen=True)
class BoundaryCell:
    edge: str
    lambda_low: Fraction
    lambda_high: Fraction
    coordinate_low: Fraction
    coordinate_high: Fraction
    direction: int = 1
    depth: int = 0

    @property
    def lambda_center(self) -> Fraction:
        return (self.lambda_low + self.lambda_high) / 2

    @property
    def lambda_radius(self) -> Fraction:
        return (self.lambda_high - self.lambda_low) / 2

    @property
    def coordinate_center(self) -> Fraction:
        return (self.coordinate_low + self.coordinate_high) / 2

    @property
    def coordinate_radius(self) -> Fraction:
        return (self.coordinate_high - self.coordinate_low) / 2

    def parameter_box(self) -> tuple[Fraction, Fraction, Fraction, Fraction]:
        if self.edge == "x_low":
            return self.coordinate_center, self.coordinate_radius, X_LOW, Fraction(0)
        if self.edge == "x_high":
            return self.coordinate_center, self.coordinate_radius, X_HIGH, Fraction(0)
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
        if self.edge == "x_low":
            return coordinate, X_LOW
        if self.edge == "x_high":
            return coordinate, X_HIGH
        if self.edge == "t_low":
            return TIME_LOW, coordinate
        if self.edge == "t_high":
            return TIME_HIGH, coordinate
        raise RuntimeError(f"unknown boundary edge: {self.edge}")

    def serialized(self) -> dict:
        return {
            "edge": self.edge,
            "lambda_low": fraction_text(self.lambda_low),
            "lambda_high": fraction_text(self.lambda_high),
            "coordinate_low": fraction_text(self.coordinate_low),
            "coordinate_high": fraction_text(self.coordinate_high),
            "direction": self.direction,
            "depth": self.depth,
        }


class ResourcePark(RuntimeError):
    """Raised between completed point rows or interval cells."""


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


def point_field(lambda_value: Fraction, time: Fraction, x: Fraction) -> list:
    family = continuation.TransformFamily(
        lambda_value,
        time,
        time,
        x,
        weighted=True,
    )
    return [family.point_jet(order) for order in range(3)]


def interval_field(cell: BoundaryCell) -> tuple[list, dict]:
    time_center, time_radius, x_center, x_radius = cell.parameter_box()
    time_high = time_center + time_radius
    weighted = continuation.TransformFamily(
        cell.lambda_center,
        time_center,
        time_high,
        x_center,
        weighted=True,
    )
    fourth = continuation.TransformFamily(
        cell.lambda_center,
        time_center,
        time_high,
        x_center,
        weighted=False,
    )
    delta_lambda = continuation.symmetric_ball(cell.lambda_radius)
    fields: list = []
    remainders: list[dict] = []
    for derivative in range(2):
        weighted_box, weighted_remainder = weighted.transform_box(
            derivative, time_radius, x_radius
        )
        fourth_box, fourth_remainder = fourth.transform_box(
            derivative, time_radius, x_radius
        )
        fields.append(weighted_box + delta_lambda * fourth_box)
        remainders.append(
            {"weighted": weighted_remainder, "fourth": fourth_remainder}
        )
    return fields, {"taylor_remainders": remainders}


def certify_boundary_cell(cell: BoundaryCell) -> dict:
    started = perf_counter()
    fields, details = interval_field(cell)
    separations = [value.abs_lower() for value in fields]
    selected = 0 if separations[0] >= separations[1] else 1
    passed = separations[selected] > 0
    return {
        "cell": cell.serialized(),
        "passed": bool(passed),
        "separated_by": ("F" if selected == 0 else "F_x") if passed else None,
        "separation_lower": continuation.arb_text(separations[selected]),
        "field_boxes": [continuation.ball_endpoints(value) for value in fields],
        "taylor_remainders": details["taylor_remainders"],
        "elapsed_seconds": perf_counter() - started,
    }


def initial_boundary_cell(
    edge: str,
    lambda_index: int,
    coordinate_index: int,
    segments: int,
) -> BoundaryCell:
    require(0 <= lambda_index < 11, "lambda index must be in 0..10")
    require(0 <= coordinate_index < segments, "coordinate index is out of range")
    lambda_low = Fraction(lambda_index, 50)
    lambda_high = Fraction(lambda_index + 1, 50)
    if edge in ("x_low", "x_high"):
        coordinate_grid = fraction_grid(TIME_LOW, TIME_HIGH, segments)
        direction = 1 if edge == "x_low" else -1
    elif edge in ("t_low", "t_high"):
        coordinate_grid = fraction_grid(X_LOW, X_HIGH, segments)
        direction = -1 if edge == "t_low" else 1
    else:
        raise RuntimeError(f"unknown pilot edge: {edge}")
    return BoundaryCell(
        edge=edge,
        lambda_low=lambda_low,
        lambda_high=lambda_high,
        coordinate_low=coordinate_grid[coordinate_index],
        coordinate_high=coordinate_grid[coordinate_index + 1],
        direction=direction,
    )


def cell_from_serialized(payload: dict) -> BoundaryCell:
    return BoundaryCell(
        edge=payload["edge"],
        lambda_low=Fraction(payload["lambda_low"]),
        lambda_high=Fraction(payload["lambda_high"]),
        coordinate_low=Fraction(payload["coordinate_low"]),
        coordinate_high=Fraction(payload["coordinate_high"]),
        direction=int(payload["direction"]),
        depth=int(payload["depth"]),
    )


def edge_coordinate_bounds(edge: str) -> tuple[Fraction, Fraction, int]:
    if edge == "x_low":
        return TIME_LOW, TIME_HIGH, 1
    if edge == "t_high":
        return X_LOW, X_HIGH, 1
    if edge == "x_high":
        return TIME_LOW, TIME_HIGH, -1
    if edge == "t_low":
        return X_LOW, X_HIGH, -1
    raise RuntimeError(f"unknown edge: {edge}")


def initial_boundary_cells(segments: int) -> list[BoundaryCell]:
    cells: list[BoundaryCell] = []
    for edge in ("x_low", "t_high", "x_high", "t_low"):
        coordinate_low, coordinate_high, direction = edge_coordinate_bounds(edge)
        grid = fraction_grid(coordinate_low, coordinate_high, segments)
        for lambda_index in range(11):
            for coordinate_index in range(segments):
                cells.append(
                    BoundaryCell(
                        edge=edge,
                        lambda_low=Fraction(lambda_index, 50),
                        lambda_high=Fraction(lambda_index + 1, 50),
                        coordinate_low=grid[coordinate_index],
                        coordinate_high=grid[coordinate_index + 1],
                        direction=direction,
                    )
                )
    return cells


def split_cell(cell: BoundaryCell, segments: int) -> list[BoundaryCell]:
    coordinate_low, coordinate_high, _ = edge_coordinate_bounds(cell.edge)
    initial_coordinate_width = (coordinate_high - coordinate_low) / segments
    normalized_lambda = (cell.lambda_high - cell.lambda_low) / INITIAL_LAMBDA_STEP
    normalized_coordinate = (
        (cell.coordinate_high - cell.coordinate_low) / initial_coordinate_width
    )
    if normalized_lambda >= normalized_coordinate and cell.lambda_high > cell.lambda_low:
        middle = (cell.lambda_low + cell.lambda_high) / 2
        return [
            BoundaryCell(
                cell.edge,
                cell.lambda_low,
                middle,
                cell.coordinate_low,
                cell.coordinate_high,
                cell.direction,
                cell.depth + 1,
            ),
            BoundaryCell(
                cell.edge,
                middle,
                cell.lambda_high,
                cell.coordinate_low,
                cell.coordinate_high,
                cell.direction,
                cell.depth + 1,
            ),
        ]
    middle = (cell.coordinate_low + cell.coordinate_high) / 2
    return [
        BoundaryCell(
            cell.edge,
            cell.lambda_low,
            cell.lambda_high,
            cell.coordinate_low,
            middle,
            cell.direction,
            cell.depth + 1,
        ),
        BoundaryCell(
            cell.edge,
            cell.lambda_low,
            cell.lambda_high,
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
            "maximum_depth": MAX_SUBDIVISION_DEPTH,
            "precision_bits": PRECISION_BITS,
        },
        "evaluated_cells": evaluated,
        "queue": [cell.serialized() for cell in queue],
        "certified": certified,
    }


def load_partial(path: Path, segments: int) -> tuple[list[BoundaryCell], list[dict], int]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    require(payload["kind"] == f"{STEM}_partial", "partial checkpoint kind mismatch")
    require(
        payload["source_sha256"] == sha256_path(Path(__file__).resolve()),
        "partial checkpoint source hash mismatch",
    )
    require(
        payload["configuration"]["boundary_segments"] == segments,
        "partial checkpoint segment count mismatch",
    )
    return (
        [cell_from_serialized(row) for row in payload["queue"]],
        list(payload["certified"]),
        int(payload["evaluated_cells"]),
    )


def minimum_record(records: list[dict], key: str) -> dict:
    require(records, "cannot minimize an empty record list")
    return min(records, key=lambda row: float(continuation.arb(row[key]).lower()))


def certify_uniform_boundary(
    segments: int,
    monitor: CpuMonitor,
    partial_path: Path,
    resume: bool,
    progress: bool,
) -> tuple[dict, list[dict], bool]:
    if resume:
        require(partial_path.exists(), "--resume requested but partial checkpoint is missing")
        queue, certified, evaluated = load_partial(partial_path, segments)
    else:
        queue = initial_boundary_cells(segments)
        certified = []
        evaluated = 0
    initial_count = 4 * 11 * segments
    parked = False
    cursor = 0
    while cursor < len(queue):
        cell = queue[cursor]
        cursor += 1
        record = certify_boundary_cell(cell)
        evaluated += 1
        if record["passed"]:
            certified.append(record)
        elif cell.depth < MAX_SUBDIVISION_DEPTH:
            queue.extend(split_cell(cell, segments))
        else:
            raise RuntimeError(
                f"boundary cell remained unresolved at depth {cell.depth}: {cell.serialized()}"
            )
        if evaluated % CHECKPOINT_EVERY == 0:
            remaining = queue[cursor:]
            continuation.write_json_atomic(
                partial_path,
                partial_payload(remaining, certified, evaluated, segments),
            )
            if progress:
                print(
                    "compact boundary: "
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

    edge_areas: dict[str, Fraction] = {}
    expected_areas: dict[str, Fraction] = {}
    for edge in ("x_low", "t_high", "x_high", "t_low"):
        edge_records = [row for row in certified if row["cell"]["edge"] == edge]
        area = sum(
            (
                Fraction(row["cell"]["lambda_high"])
                - Fraction(row["cell"]["lambda_low"])
            )
            * (
                Fraction(row["cell"]["coordinate_high"])
                - Fraction(row["cell"]["coordinate_low"])
            )
            for row in edge_records
        )
        low, high, _ = edge_coordinate_bounds(edge)
        expected = (LAMBDA_HIGH - LAMBDA_LOW) * (high - low)
        edge_areas[edge] = area
        expected_areas[edge] = expected
        if not parked:
            require(area == expected, f"boundary area audit failed on {edge}")

    minimum = minimum_record(certified, "separation_lower")
    summary = {
        "initial_cells": initial_count,
        "evaluated_cells": evaluated,
        "certified_leaf_cells": len(certified),
        "queued_cells": len(queue),
        "resource_parked": parked,
        "maximum_depth": max(int(row["cell"]["depth"]) for row in certified),
        "edge_counts": {
            edge: sum(row["cell"]["edge"] == edge for row in certified)
            for edge in ("x_low", "t_high", "x_high", "t_low")
        },
        "edge_area_audit": {
            edge: {
                "certified_area": fraction_text(edge_areas[edge]),
                "expected_area": fraction_text(expected_areas[edge]),
                "passed": edge_areas[edge] == expected_areas[edge],
            }
            for edge in edge_areas
        },
        "minimum_separation_lower": minimum["separation_lower"],
        "minimum_separation_cell": minimum["cell"],
        "all_boundary_cells_nonzero": not parked and not queue,
    }
    return summary, certified, not parked


def ordered_lambda_zero_records(records: list[dict]) -> list[dict]:
    selected = [
        row for row in records if Fraction(row["cell"]["lambda_low"]) == LAMBDA_LOW
    ]
    ordered: list[dict] = []
    for edge in ("x_low", "t_high", "x_high", "t_low"):
        edge_rows = [row for row in selected if row["cell"]["edge"] == edge]
        direction = edge_coordinate_bounds(edge)[2]
        edge_rows.sort(
            key=lambda row: Fraction(row["cell"]["coordinate_low"]),
            reverse=direction < 0,
        )
        coordinate_low, coordinate_high, _ = edge_coordinate_bounds(edge)
        covered = sum(
            Fraction(row["cell"]["coordinate_high"])
            - Fraction(row["cell"]["coordinate_low"])
            for row in edge_rows
        )
        require(
            covered == coordinate_high - coordinate_low,
            f"lambda=0 edge coverage failed on {edge}",
        )
        ordered.extend(edge_rows)
    require(ordered, "lambda=0 boundary partition is empty")
    return ordered


def transformed_coordinates(fields: list, ray: tuple[int, int]) -> tuple:
    a, b = ray
    along = a * fields[0] + b * fields[1]
    transverse = -b * fields[0] + a * fields[1]
    return along, transverse


def strict_sign(value) -> int | None:
    if value.upper() < 0:
        return -1
    if value.lower() > 0:
        return 1
    return None


def classify_polygon_winding(
    vertices: list[dict],
) -> dict:
    candidates = (
        (1, 0),
        (1, 1),
        (2, 1),
        (1, 2),
        (3, 1),
        (1, 3),
        (-1, 2),
        (2, -1),
    )
    for ray in candidates:
        transformed = [transformed_coordinates(row["fields"], ray) for row in vertices]
        signs = [strict_sign(pair[1]) for pair in transformed]
        if any(sign is None for sign in signs):
            continue
        crossings: list[dict] = []
        winding = 0
        failed = False
        for index, ((u_left, v_left), left_sign) in enumerate(
            zip(transformed, signs, strict=True)
        ):
            right_index = (index + 1) % len(vertices)
            u_right, v_right = transformed[right_index]
            right_sign = signs[right_index]
            if left_sign == right_sign:
                continue
            denominator = v_right - v_left
            u_cross = (u_left * v_right - u_right * v_left) / denominator
            crossing_sign = strict_sign(u_cross)
            if crossing_sign is None:
                failed = True
                break
            if crossing_sign > 0:
                contribution = 1 if left_sign < right_sign else -1
                winding += contribution
                crossings.append(
                    {
                        "edge_index": index,
                        "from_vertex": vertices[index]["location"],
                        "to_vertex": vertices[right_index]["location"],
                        "orientation": "up" if contribution > 0 else "down",
                        "contribution": contribution,
                        "ray_coordinate": continuation.ball_endpoints(u_cross),
                    }
                )
        if failed:
            continue
        return {
            "ray": list(ray),
            "vertex_transverse_signs_strict": True,
            "crossings": crossings,
            "winding_t_x": winding,
        }
    raise RuntimeError("no rational ray gave a strict polygon crossing classification")


def certify_lambda_zero_winding(records: list[dict]) -> dict:
    ordered = ordered_lambda_zero_records(records)
    vertices: list[dict] = []
    point_cache: dict[tuple[Fraction, Fraction], list] = {}
    previous_end: tuple[Fraction, Fraction] | None = None
    for index, record in enumerate(ordered):
        cell = cell_from_serialized(record["cell"])
        start = cell.point(cell.oriented_start())
        end = cell.point(cell.oriented_end())
        if previous_end is not None:
            require(start == previous_end, f"oriented boundary discontinuity before cell {index}")
        previous_end = end
        if start not in point_cache:
            point_cache[start] = point_field(LAMBDA_LOW, start[0], start[1])[:2]
        if end not in point_cache:
            point_cache[end] = point_field(LAMBDA_LOW, end[0], end[1])[:2]
        field_box = [
            continuation.arb(item["ball"]) for item in record["field_boxes"]
        ]
        for endpoint in (start, end):
            endpoint_fields = point_cache[endpoint]
            require(
                all(field_box[q].contains(endpoint_fields[q]) for q in range(2)),
                f"segment field box does not contain endpoint on {cell.edge}",
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
    first_cell = cell_from_serialized(ordered[0]["cell"])
    require(
        previous_end == first_cell.point(first_cell.oriented_start()),
        "oriented boundary does not close",
    )
    winding = classify_polygon_winding(vertices)
    require(winding["winding_t_x"] == -1, "lambda=0 boundary winding is not -1")
    return {
        "lambda": "0",
        "orientation": (
            "positive boundary orientation in coordinate order (t,x): "
            "x=x_low upward in t, t=t_high rightward in x, "
            "x=x_high downward in t, t=t_low leftward in x"
        ),
        "certified_segments": len(ordered),
        "convex_segment_homotopy": (
            "Each true segment image and its endpoint chord lie in the same "
            "convex interval rectangle, which excludes the origin by one component."
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


def certify_branch_containment() -> dict:
    payload = json.loads(INTERVAL_PARENT.read_text(encoding="utf-8"))
    charts = payload["continuation"]["charts"]
    records: list[dict] = []
    for chart in charts:
        lambda_low = Fraction(chart["lambda"]["low"])
        lambda_high = Fraction(chart["lambda"]["high"])
        time_center = Fraction(chart["root_box"]["time_center"])
        time_radius = Fraction(chart["root_box"]["time_radius"])
        x_center = Fraction(chart["root_box"]["x_center"])
        x_radius = Fraction(chart["root_box"]["x_radius"])
        margins = {
            "lambda_low": lambda_low - LAMBDA_LOW,
            "lambda_high": LAMBDA_HIGH - lambda_high,
            "time_low": time_center - time_radius - TIME_LOW,
            "time_high": TIME_HIGH - time_center - time_radius,
            "x_low": x_center - x_radius - X_LOW,
            "x_high": X_HIGH - x_center - x_radius,
        }
        require(all(value >= 0 for value in margins.values()), "branch chart left cylinder")
        require(
            chart["signs"]["F_xx_strictly_negative"] is True,
            "branch chart lost regularity",
        )
        records.append(
            {
                "chart": chart["index"],
                "margins": {key: fraction_text(value) for key, value in margins.items()},
                "F_xx": chart["signs"]["F_xx"],
            }
        )
    spatial_margins = [
        Fraction(value)
        for row in records
        for key, value in row["margins"].items()
        if key in ("time_low", "time_high", "x_low", "x_high")
    ]
    minimum = min(spatial_margins)
    require(minimum > 0, "branch chart touches the spatial cylinder boundary")
    return {
        "charts": len(records),
        "all_charts_inside_cylinder": True,
        "all_branch_contacts_regular": True,
        "minimum_spatial_margin": fraction_text(minimum),
        "records": records,
    }


def exact_degree_composition(boundary: dict, winding: dict, branch: dict) -> dict:
    index_payload = json.loads(INDEX_PARENT.read_text(encoding="utf-8"))
    rows_text = json.dumps(index_payload, sort_keys=True)
    require("+floor(m/2)" in rows_text, "index parent lacks all-multiplicity charge")
    require(boundary["all_boundary_cells_nonzero"] is True, "boundary is not complete")
    require(winding["polygon"]["winding_t_x"] == -1, "wrong base winding")
    require(branch["all_charts_inside_cylinder"] is True, "branch is not inside cylinder")
    return {
        "nontriviality": (
            "For every lambda in [0,0.22], phi_n(u)>0 and the coefficients are "
            "nonnegative with the first three positive, so F_lambda(t,0)>0; "
            "the heat solution is not identically zero."
        ),
        "imported_local_index": (
            "A multiplicity-m real heat contact has local index +floor(m/2) "
            "in (x,t), hence -floor(m/2) in (t,x), including degenerate contacts."
        ),
        "boundary_homotopy": (
            "Uniform boundary nonvanishing makes the boundary winding independent "
            "of lambda; the certified lambda=0 value is -1 in (t,x)."
        ),
        "tracked_charge": (
            "The interval branch supplies one regular interior contact for every "
            "lambda and F_xx<0, so it contributes local index -1."
        ),
        "complete_count": (
            "Every additional contact would contribute a further strictly negative "
            "integer, contradicting total degree -1. Therefore the tracked contact "
            "is the unique contact in the compact cylinder for every lambda."
        ),
        "triple_contact_note": (
            "A separate exclusion of F=F_x=F_xx=0 is unnecessary because the "
            "all-multiplicity index theorem already assigns every triple or higher "
            "contact a nonzero charge of the same sign."
        ),
    }


def render_note(payload: dict) -> str:
    boundary = payload["boundary_certificate"]
    winding = payload["lambda_zero_winding"]
    lines = [
        "# Four-Summand Compact Contact Degree Certificate",
        "",
        f"Date: {DATE}",
        "",
        "Status: rigorous finite-homotopy compact contact count. This is not a proof of RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py --progress",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Compact Theorem",
        "",
        "For every `0<=lambda<=0.22`, the tracked four-summand branch contact is",
        "the unique common zero of `(F_lambda,F_lambda,x)` in",
        "",
        "```text",
        "-0.06<=t<=0.45,  135.5<=x<=136.",
        "```",
        "",
        f"The uniform boundary cover has `{boundary['certified_leaf_cells']}` interval cells;",
        f"its minimum one-component separation is `{boundary['minimum_separation_lower']}`.",
        f"At lambda zero, `{winding['certified_segments']}` convex image segments give",
        f"the exact `(t,x)` boundary winding `{winding['polygon']['winding_t_x']}`.",
        "",
        "The all-multiplicity heat-contact theorem assigns local index",
        "`-floor(m/2)` in `(t,x)`. The certified regular branch already contributes",
        "`-1`; any hidden contact would make the total degree more negative.",
        "",
        "## Why No Triple-Contact Search Is Needed",
        "",
        payload["degree_composition"]["triple_contact_note"],
        "",
        "## Proof Boundary",
        "",
        payload["proof_boundary"],
        "",
    ]
    return "\n".join(lines)


def build_full_payload(
    args: argparse.Namespace,
    baseline: list[float],
    priority_lowered: bool,
    monitor: CpuMonitor,
) -> tuple[dict, bool]:
    started = perf_counter()
    boundary, records, complete = certify_uniform_boundary(
        args.boundary_segments,
        monitor,
        args.partial,
        args.resume,
        args.progress,
    )
    if not complete:
        payload = {
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
        }
        return payload, False

    winding = certify_lambda_zero_winding(records)
    branch = certify_branch_containment()
    degree = exact_degree_composition(boundary, winding, branch)
    payload = {
        "kind": STEM,
        "schema_version": SCHEMA_VERSION,
        "date": DATE,
        "status": "rigorous compact contact degree certificate complete",
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "parents": {
            "interval_continuation": {
                "path": str(INTERVAL_PARENT.relative_to(REPO_ROOT)).replace("\\", "/"),
                "sha256": sha256_path(INTERVAL_PARENT),
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
        "cylinder": {
            "lambda": [fraction_text(LAMBDA_LOW), fraction_text(LAMBDA_HIGH)],
            "time": [fraction_text(TIME_LOW), fraction_text(TIME_HIGH)],
            "x": [fraction_text(X_LOW), fraction_text(X_HIGH)],
            "coordinate_orientation": "(t,x)",
        },
        "interval_configuration": {
            "precision_bits": PRECISION_BITS,
            "absolute_integration_tolerance": continuation.ABS_TOL,
            "taylor_time_order": continuation.TAYLOR_TIME_ORDER,
            "taylor_x_order": continuation.TAYLOR_X_ORDER,
            "boundary_segments_per_edge_per_lambda_slab": args.boundary_segments,
            "initial_lambda_step": fraction_text(INITIAL_LAMBDA_STEP),
            "maximum_subdivision_depth": MAX_SUBDIVISION_DEPTH,
            "cutoff_tail_radius": continuation.TAIL_RADIUS,
        },
        "boundary_certificate": boundary,
        "boundary_cells": records,
        "lambda_zero_winding": winding,
        "branch_containment": branch,
        "degree_composition": degree,
        "theorem": (
            "For every 0<=lambda<=0.22, the certified four-summand branch contact "
            "is the unique common zero of (F_lambda,partial_x F_lambda) in "
            "[-0.06,0.45]x[135.5,136]."
        ),
        "proof_boundary": (
            "This counts contacts only for the finite homotopy H_1+H_2+H_3+lambda H_4 "
            "inside the declared compact cylinder and only for 0<=lambda<=0.22. It does "
            "not control contacts outside that cylinder, lambda>0.22, summands n>=5, the "
            "complete theta/Xi transform, an Xi collision, a degree-uniform Jensen remainder, "
            "Lambda<=0, RH, or a prize-level conclusion."
        ),
    }
    return payload, True


def boundary_vertices(segments_per_edge: int) -> list[tuple[str, Fraction, Fraction]]:
    vertices: list[tuple[str, Fraction, Fraction]] = []
    # Positive orientation in the (t,x) coordinate order.
    for time in fraction_grid(TIME_LOW, TIME_HIGH, segments_per_edge)[:-1]:
        vertices.append(("x_low", time, X_LOW))
    for x in fraction_grid(X_LOW, X_HIGH, segments_per_edge)[:-1]:
        vertices.append(("t_high", TIME_HIGH, x))
    for time in reversed(fraction_grid(TIME_LOW, TIME_HIGH, segments_per_edge)[1:]):
        vertices.append(("x_high", time, X_HIGH))
    for x in reversed(fraction_grid(X_LOW, X_HIGH, segments_per_edge)[1:]):
        vertices.append(("t_low", TIME_LOW, x))
    return vertices


def midpoint_float(value) -> float:
    return float(value.mid())


def polygon_winding(points: list[complex]) -> tuple[int, float]:
    total = 0.0
    for index, left in enumerate(points):
        right = points[(index + 1) % len(points)]
        delta = math.atan2(right.imag, right.real) - math.atan2(left.imag, left.real)
        while delta <= -math.pi:
            delta += 2 * math.pi
        while delta > math.pi:
            delta -= 2 * math.pi
        total += delta
    winding = round(total / (2 * math.pi))
    return winding, total / (2 * math.pi)


def point_scout(segments_per_edge: int, monitor: CpuMonitor) -> dict:
    started = perf_counter()
    lambdas = [LAMBDA_LOW, (LAMBDA_LOW + LAMBDA_HIGH) / 2, LAMBDA_HIGH]
    vertices = boundary_vertices(segments_per_edge)
    rows: list[dict] = []
    for lambda_value in lambdas:
        points: list[complex] = []
        weakest: dict | None = None
        curvature_min = math.inf
        for edge, time, x in vertices:
            jets = point_field(lambda_value, time, x)
            real = midpoint_float(jets[0])
            imag = midpoint_float(jets[1])
            curvature = midpoint_float(jets[2])
            points.append(complex(real, imag))
            norm = math.hypot(real, imag)
            if weakest is None or norm < weakest["norm"]:
                weakest = {
                    "edge": edge,
                    "time": fraction_text(time),
                    "x": fraction_text(x),
                    "F": f"{real:.17e}",
                    "F_x": f"{imag:.17e}",
                    "norm": norm,
                }
            curvature_min = min(curvature_min, abs(curvature))
        winding, raw_winding = polygon_winding(points)
        rows.append(
            {
                "lambda": fraction_text(lambda_value),
                "sampled_boundary_winding_t_x": winding,
                "raw_winding": raw_winding,
                "weakest_point": weakest,
                "minimum_sampled_abs_F_xx": curvature_min,
            }
        )
        monitor.sample()
    return {
        "status": "point scout only; not an interval certificate",
        "rectangle": {
            "lambda": [fraction_text(LAMBDA_LOW), fraction_text(LAMBDA_HIGH)],
            "time": [fraction_text(TIME_LOW), fraction_text(TIME_HIGH)],
            "x": [fraction_text(X_LOW), fraction_text(X_HIGH)],
        },
        "segments_per_edge": segments_per_edge,
        "rows": rows,
        "elapsed_seconds": perf_counter() - started,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--point-scout", action="store_true")
    parser.add_argument("--segments-per-edge", type=int, default=16)
    parser.add_argument("--pilot-edge", choices=("x_low", "t_high", "x_high", "t_low"))
    parser.add_argument("--pilot-lambda-index", type=int, default=0)
    parser.add_argument("--pilot-coordinate-index", type=int, default=0)
    parser.add_argument("--boundary-segments", type=int, default=INITIAL_BOUNDARY_SEGMENTS)
    parser.add_argument("--partial", type=Path, default=DEFAULT_PARTIAL)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--baseline-seconds", type=int, default=5)
    parser.add_argument("--baseline-max-percent", type=float, default=60.0)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    priority_lowered = continuation.request_below_normal_priority()
    baseline = [float(psutil.cpu_percent(interval=1.0)) for _ in range(args.baseline_seconds)]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "compact degree baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, mean={baseline_mean:.1f}%",
        flush=True,
    )
    if baseline_mean > args.baseline_max_percent:
        print("compact degree run deferred: daytime baseline is already busy")
        return 2
    require(priority_lowered, "could not apply below-normal process priority")
    continuation.flint.ctx.prec = PRECISION_BITS
    monitor = CpuMonitor()
    if args.point_scout:
        payload = point_scout(args.segments_per_edge, monitor)
        payload["resource_policy"] = {
            "active_compute_workers": 1,
            "thread_caps": 1,
            "below_normal_priority_applied": priority_lowered,
            "baseline_cpu_percent": baseline,
            "runtime_cpu_percent": monitor.samples,
        }
        print(json.dumps(payload, indent=2))
        return 0
    if args.pilot_edge is not None:
        cell = initial_boundary_cell(
            args.pilot_edge,
            args.pilot_lambda_index,
            args.pilot_coordinate_index,
            args.boundary_segments,
        )
        record = certify_boundary_cell(cell)
        record["resource_policy"] = {
            "active_compute_workers": 1,
            "thread_caps": 1,
            "below_normal_priority_applied": priority_lowered,
            "baseline_cpu_percent": baseline,
        }
        print(json.dumps(record, indent=2))
        return 0 if record["passed"] else 4
    payload, complete = build_full_payload(
        args,
        baseline,
        priority_lowered,
        monitor,
    )
    if not complete:
        continuation.write_json_atomic(args.out, payload)
        print(
            "compact degree certificate parked: "
            f"certified={payload['boundary_certificate']['certified_leaf_cells']}, "
            f"queued={payload['boundary_certificate']['queued_cells']}"
        )
        return 3
    continuation.write_json_atomic(args.out, payload)
    continuation.write_text_atomic(args.note, render_note(payload))
    if args.partial.exists():
        args.partial.unlink()
    print(
        "built compact contact degree certificate: "
        f"boundary_cells={payload['boundary_certificate']['certified_leaf_cells']}, "
        f"winding={payload['lambda_zero_winding']['polygon']['winding_t_x']}, "
        f"unique_contact=True, elapsed={payload['resource_policy']['elapsed_seconds']:.3f}s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
