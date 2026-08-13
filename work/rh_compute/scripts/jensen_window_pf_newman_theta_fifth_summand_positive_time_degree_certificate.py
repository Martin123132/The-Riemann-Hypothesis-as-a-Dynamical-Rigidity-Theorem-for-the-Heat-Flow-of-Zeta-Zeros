#!/usr/bin/env python3
"""Certify the fifth-summand homotopy on the local positive-time rectangle."""

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

import jensen_window_pf_newman_theta_fifth_summand_positive_time_boundary_scout as scout  # noqa: E402
import jensen_window_pf_newman_theta_fourth_summand_compact_contact_degree_certificate as compact  # noqa: E402


continuation = scout.continuation
STEM = "jensen_window_pf_newman_theta_fifth_summand_positive_time_degree_certificate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
DEFAULT_PARTIAL = REPO_ROOT / "work/rh_compute/results" / f"{STEM}_partial.json"
SCOUT_PARENT = REPO_ROOT / "work/rh_compute/results" / f"{scout.STEM}.json"
FOUR_TERM_PARENT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_fourth_summand_"
    "positive_time_degree_certificate.json"
)
INDEX_PARENT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_first_jet_winding_gate.json"
)

SCHEMA_VERSION = 1
DATE = "2026-08-04"
MU_LOW = Fraction(0)
MU_HIGH = Fraction(1)
TIME_LOW = Fraction(0)
TIME_HIGH = Fraction(9, 20)
X_LOW = Fraction(271, 2)
X_HIGH = Fraction(136)
PRECISION_BITS = 192
INITIAL_MU_STEP = Fraction(1)
INITIAL_BOUNDARY_SEGMENTS = 8
MAX_SUBDIVISION_DEPTH = 10
CHECKPOINT_EVERY = 8


@dataclass(frozen=True)
class BoundaryCell:
    edge: str
    mu_low: Fraction
    mu_high: Fraction
    coordinate_low: Fraction
    coordinate_high: Fraction
    direction: int = 1
    depth: int = 0

    @property
    def mu_center(self) -> Fraction:
        return (self.mu_low + self.mu_high) / 2

    @property
    def mu_radius(self) -> Fraction:
        return (self.mu_high - self.mu_low) / 2

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
            "mu_low": fraction_text(self.mu_low),
            "mu_high": fraction_text(self.mu_high),
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


def point_field(mu: Fraction, time: Fraction, x: Fraction) -> list:
    family = scout.FifthTransformFamily(mu, time, time, x, weighted=True)
    return [family.point_jet(order) for order in range(2)]


def interval_field(cell: BoundaryCell) -> tuple[list, dict]:
    time_center, time_radius, x_center, x_radius = cell.parameter_box()
    time_high = time_center + time_radius
    base = scout.FifthTransformFamily(
        Fraction(0),
        time_center,
        time_high,
        x_center,
        weighted=True,
    )
    fifth = scout.FifthTransformFamily(
        Fraction(0),
        time_center,
        time_high,
        x_center,
        weighted=False,
    )
    mu_box = (
        continuation.arb_fraction(cell.mu_center)
        + continuation.symmetric_ball(cell.mu_radius)
    )
    fields: list = []
    base_boxes: list = []
    fifth_boxes: list = []
    remainders: list[dict] = []
    for derivative in range(2):
        base_box, base_remainder = base.transform_box(
            derivative, time_radius, x_radius
        )
        fifth_box, fifth_remainder = fifth.transform_box(
            derivative, time_radius, x_radius
        )
        base_boxes.append(base_box)
        fifth_boxes.append(fifth_box)
        fields.append(base_box + mu_box * fifth_box)
        remainders.append({"base": base_remainder, "fifth": fifth_remainder})
    return fields, {
        "base_boxes": base_boxes,
        "fifth_boxes": fifth_boxes,
        "taylor_remainders": remainders,
    }


def certify_boundary_cell(cell: BoundaryCell) -> dict:
    started = perf_counter()
    fields, details = interval_field(cell)
    separations = [value.abs_lower() for value in fields]
    selected = 0 if separations[0] >= separations[1] else 1
    passed = separations[selected] > 0
    fifth_upper = details["fifth_boxes"][selected].abs_upper()
    perturbation_ratio = (
        fifth_upper / separations[selected] if passed else continuation.arb(0)
    )
    return {
        "cell": cell.serialized(),
        "passed": bool(passed),
        "separated_by": ("F" if selected == 0 else "F_x") if passed else None,
        "separation_lower": continuation.arb_text(separations[selected]),
        "selected_fifth_abs_upper": continuation.arb_text(fifth_upper),
        "selected_fifth_to_margin_upper": continuation.arb_text(
            perturbation_ratio.upper()
        ),
        "field_boxes": [continuation.ball_endpoints(value) for value in fields],
        "base_boxes": [
            continuation.ball_endpoints(value) for value in details["base_boxes"]
        ],
        "fifth_boxes": [
            continuation.ball_endpoints(value) for value in details["fifth_boxes"]
        ],
        "taylor_remainders": details["taylor_remainders"],
        "elapsed_seconds": perf_counter() - started,
    }


def cell_from_serialized(payload: dict) -> BoundaryCell:
    return BoundaryCell(
        edge=payload["edge"],
        mu_low=Fraction(payload["mu_low"]),
        mu_high=Fraction(payload["mu_high"]),
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
        for coordinate_index in range(segments):
            cells.append(
                BoundaryCell(
                    edge=edge,
                    mu_low=MU_LOW,
                    mu_high=MU_HIGH,
                    coordinate_low=grid[coordinate_index],
                    coordinate_high=grid[coordinate_index + 1],
                    direction=direction,
                )
            )
    return cells


def split_cell(cell: BoundaryCell, segments: int) -> list[BoundaryCell]:
    coordinate_low, coordinate_high, _ = edge_coordinate_bounds(cell.edge)
    coordinate_width = (coordinate_high - coordinate_low) / segments
    normalized_mu = (cell.mu_high - cell.mu_low) / INITIAL_MU_STEP
    normalized_coordinate = (
        (cell.coordinate_high - cell.coordinate_low) / coordinate_width
    )
    if normalized_mu >= normalized_coordinate:
        middle = (cell.mu_low + cell.mu_high) / 2
        return [
            BoundaryCell(
                cell.edge,
                cell.mu_low,
                middle,
                cell.coordinate_low,
                cell.coordinate_high,
                cell.direction,
                cell.depth + 1,
            ),
            BoundaryCell(
                cell.edge,
                middle,
                cell.mu_high,
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
            cell.mu_low,
            cell.mu_high,
            cell.coordinate_low,
            middle,
            cell.direction,
            cell.depth + 1,
        ),
        BoundaryCell(
            cell.edge,
            cell.mu_low,
            cell.mu_high,
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
            continuation.arb(row["selected_fifth_to_margin_upper"]).upper()
        ),
    )


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
    initial_count = 4 * segments
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
                    "fifth-summand boundary: "
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
                Fraction(row["cell"]["mu_high"])
                - Fraction(row["cell"]["mu_low"])
            )
            * (
                Fraction(row["cell"]["coordinate_high"])
                - Fraction(row["cell"]["coordinate_low"])
            )
            for row in edge_records
        )
        low, high, _ = edge_coordinate_bounds(edge)
        expected = (MU_HIGH - MU_LOW) * (high - low)
        edge_areas[edge] = area
        expected_areas[edge] = expected
        if not parked:
            require(area == expected, f"boundary area audit failed on {edge}")

    minimum = minimum_record(certified)
    maximum_ratio = maximum_ratio_record(certified)
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
        "maximum_selected_fifth_to_margin_upper": maximum_ratio[
            "selected_fifth_to_margin_upper"
        ],
        "maximum_ratio_cell": maximum_ratio["cell"],
        "all_boundary_cells_nonzero": not parked and not queue,
    }
    return summary, certified, not parked


def ordered_base_records(records: list[dict]) -> list[dict]:
    selected = [row for row in records if Fraction(row["cell"]["mu_low"]) == MU_LOW]
    ordered: list[dict] = []
    for edge in ("x_low", "t_high", "x_high", "t_low"):
        edge_rows = [row for row in selected if row["cell"]["edge"] == edge]
        direction = edge_coordinate_bounds(edge)[2]
        edge_rows.sort(
            key=lambda row: Fraction(row["cell"]["coordinate_low"]),
            reverse=direction < 0,
        )
        low, high, _ = edge_coordinate_bounds(edge)
        covered = sum(
            Fraction(row["cell"]["coordinate_high"])
            - Fraction(row["cell"]["coordinate_low"])
            for row in edge_rows
        )
        require(covered == high - low, f"mu=0 edge coverage failed on {edge}")
        ordered.extend(edge_rows)
    require(ordered, "mu=0 boundary partition is empty")
    return ordered


def certify_base_winding(records: list[dict]) -> dict:
    ordered = ordered_base_records(records)
    vertices: list[dict] = []
    point_cache: dict[tuple[Fraction, Fraction], list] = {}
    previous_end: tuple[Fraction, Fraction] | None = None
    for index, record in enumerate(ordered):
        cell = cell_from_serialized(record["cell"])
        start = cell.point(cell.oriented_start())
        end = cell.point(cell.oriented_end())
        if previous_end is not None:
            require(start == previous_end, f"oriented boundary discontinuity before {index}")
        previous_end = end
        if start not in point_cache:
            point_cache[start] = point_field(MU_LOW, start[0], start[1])
        if end not in point_cache:
            point_cache[end] = point_field(MU_LOW, end[0], end[1])
        field_box = [continuation.arb(item["ball"]) for item in record["field_boxes"]]
        for endpoint in (start, end):
            endpoint_fields = point_cache[endpoint]
            require(
                all(field_box[q].contains(endpoint_fields[q]) for q in range(2)),
                f"segment field box omits endpoint on {cell.edge}",
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
    winding = compact.classify_polygon_winding(vertices)
    require(winding["winding_t_x"] == 0, "mu=0 boundary winding is not zero")
    return {
        "mu": "0",
        "orientation": (
            "positive boundary orientation in coordinate order (t,x): "
            "x=x_low upward in t, t=t_high rightward in x, "
            "x=x_high downward in t, t=t_low leftward in x"
        ),
        "certified_segments": len(ordered),
        "convex_segment_homotopy": (
            "Each true mu=0 segment image and its endpoint chord lie in the same "
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


def exact_degree_composition(boundary: dict, winding: dict) -> dict:
    index_payload = json.loads(INDEX_PARENT.read_text(encoding="utf-8"))
    rows_text = json.dumps(index_payload, sort_keys=True)
    require("+floor(m/2)" in rows_text, "index parent lacks all-multiplicity charge")
    require(boundary["all_boundary_cells_nonzero"] is True, "boundary is incomplete")
    require(winding["polygon"]["winding_t_x"] == 0, "base winding is not zero")
    return {
        "nontriviality": (
            "For 0<=mu<=1 all five theta coefficients are nonnegative and the first "
            "four are positive, so F_mu(t,0)>0 and the heat solution is nontrivial."
        ),
        "imported_local_index": (
            "Every multiplicity-m real heat contact has local index +floor(m/2) "
            "in (x,t), hence -floor(m/2) in (t,x), including degenerate contacts."
        ),
        "boundary_homotopy": (
            "Uniform lateral-boundary nonvanishing preserves the mu=0 boundary "
            "winding zero throughout 0<=mu<=1."
        ),
        "zero_degree_exclusion": (
            "Any contact would contribute a strictly negative integer to total degree. "
            "Since the degree is zero, the rectangle contains no contacts."
        ),
    }


def render_note(payload: dict) -> str:
    boundary = payload["boundary_certificate"]
    winding = payload["mu_zero_winding"]
    lines = [
        "# Fifth-Summand Positive-Time Degree Certificate",
        "",
        f"Date: {DATE}",
        "",
        "Status: rigorous finite-homotopy contact-exclusion theorem. This is not a proof of RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py --progress",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Compact Theorem",
        "",
        "For every `0<=mu<=1`, the finite transform",
        "`H_1+H_2+H_3+H_4+mu H_5` has no common zero of `(F,F_x)` in",
        "",
        "```text",
        "0<=t<=0.45,  135.5<=x<=136.",
        "```",
        "",
        f"The lateral-boundary cover has `{boundary['certified_leaf_cells']}` cells;",
        f"its minimum separation is `{boundary['minimum_separation_lower']}`.",
        f"At mu zero, `{winding['certified_segments']}` convex image segments give",
        "winding zero. Every possible contact has negative local index, so degree",
        "zero excludes regular and degenerate contacts alike.",
        "",
        "## Scaling Signal",
        "",
        "The initial coefficient width is the whole interval `[0,1]`. The stored",
        "fifth-component boxes and perturbation-to-margin ratios provide the first",
        "calibration row for a uniform later-summand tail theorem.",
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
    tail_audit = scout.extended_tail_audit()
    boundary, records, complete = certify_uniform_boundary(
        args.boundary_segments,
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

    winding = certify_base_winding(records)
    degree = exact_degree_composition(boundary, winding)
    payload = {
        "kind": STEM,
        "schema_version": SCHEMA_VERSION,
        "date": DATE,
        "status": "rigorous fifth-summand positive-time exclusion certificate complete",
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "parents": {
            "point_scout": {
                "path": str(SCOUT_PARENT.relative_to(REPO_ROOT)).replace("\\", "/"),
                "sha256": sha256_path(SCOUT_PARENT),
            },
            "four_term_degree": {
                "path": str(FOUR_TERM_PARENT.relative_to(REPO_ROOT)).replace("\\", "/"),
                "sha256": sha256_path(FOUR_TERM_PARENT),
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
            "mu": ["0", "1"],
            "time": [fraction_text(TIME_LOW), fraction_text(TIME_HIGH)],
            "x": [fraction_text(X_LOW), fraction_text(X_HIGH)],
            "coordinate_orientation": "(t,x)",
        },
        "family": "H_1+H_2+H_3+H_4+mu H_5",
        "tail_audit": tail_audit,
        "interval_configuration": {
            "precision_bits": PRECISION_BITS,
            "absolute_integration_tolerance": continuation.ABS_TOL,
            "taylor_time_order": continuation.TAYLOR_TIME_ORDER,
            "taylor_x_order": continuation.TAYLOR_X_ORDER,
            "boundary_segments_per_edge": args.boundary_segments,
            "initial_mu_step": fraction_text(INITIAL_MU_STEP),
            "maximum_subdivision_depth": MAX_SUBDIVISION_DEPTH,
            "cutoff_tail_radius_per_transform": continuation.TAIL_RADIUS,
        },
        "boundary_certificate": boundary,
        "boundary_cells": records,
        "mu_zero_winding": winding,
        "degree_composition": degree,
        "theorem": (
            "For every 0<=mu<=1, H_1+H_2+H_3+H_4+mu H_5 and its first x "
            "derivative have no common zero in [0,0.45]x[135.5,136]."
        ),
        "scaling_frontier": (
            "The full fifth coefficient is handled in one initial mu slab if "
            "maximum_depth=0. Use the stored H_5 derivative boxes to seek a summable "
            "all-n boundary perturbation bound before adding further terms individually."
        ),
        "proof_boundary": (
            "This excludes contacts only for the finite fifth-summand homotopy inside "
            "the declared local positive-time rectangle. It does not control other "
            "frequency rectangles, negative heat time, summands n>=6, the complete "
            "theta/Xi transform, an Xi collision, a degree-uniform Jensen remainder, "
            "Lambda<=0, RH, or a prize-level conclusion."
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
    priority_lowered = continuation.request_below_normal_priority()
    baseline = [float(psutil.cpu_percent(interval=1.0)) for _ in range(args.baseline_seconds)]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "fifth-summand degree baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, "
        f"mean={baseline_mean:.1f}%",
        flush=True,
    )
    if baseline_mean > args.baseline_max_percent:
        print("fifth-summand degree run deferred: daytime baseline is already busy")
        return 2
    require(priority_lowered, "could not apply below-normal process priority")
    continuation.flint.ctx.prec = PRECISION_BITS
    monitor = CpuMonitor()
    payload, complete = build_full_payload(args, baseline, priority_lowered, monitor)
    if not complete:
        continuation.write_json_atomic(args.out, payload)
        print(
            "fifth-summand degree certificate parked: "
            f"certified={payload['boundary_certificate']['certified_leaf_cells']}, "
            f"queued={payload['boundary_certificate']['queued_cells']}"
        )
        return 3
    continuation.write_json_atomic(args.out, payload)
    continuation.write_text_atomic(args.note, render_note(payload))
    if args.partial.exists():
        args.partial.unlink()
    print(
        "built fifth-summand degree certificate: "
        f"boundary_cells={payload['boundary_certificate']['certified_leaf_cells']}, "
        f"max_depth={payload['boundary_certificate']['maximum_depth']}, "
        f"winding={payload['mu_zero_winding']['polygon']['winding_t_x']}, "
        f"no_contacts=True, elapsed={payload['resource_policy']['elapsed_seconds']:.3f}s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
