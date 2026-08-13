#!/usr/bin/env python3
"""Check the complete-theta positive-time spatial tile degree certificate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import importlib
import json
import os
from pathlib import Path
import sys

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

from flint import arb  # noqa: E402


STEM = "jensen_window_pf_newman_theta_complete_spatial_tile_degree_certificate"
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
SOURCE = SCRIPT_DIR / f"{STEM}.py"


class CheckFailure(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CheckFailure(message)


def strict_sign(value: arb) -> int | None:
    if value.upper() < 0:
        return -1
    if value.lower() > 0:
        return 1
    return None


def ray_winding(vertices: list[list[arb]], ray: tuple[int, int]) -> tuple[int, int]:
    a, b = ray
    transformed = [
        (a * fields[0] + b * fields[1], -b * fields[0] + a * fields[1])
        for fields in vertices
    ]
    signs = [strict_sign(pair[1]) for pair in transformed]
    require(all(sign is not None for sign in signs), f"ray {ray} hits a vertex ball")
    winding = 0
    positive_crossings = 0
    for index, ((u_left, v_left), left_sign) in enumerate(
        zip(transformed, signs, strict=True)
    ):
        right_index = (index + 1) % len(vertices)
        u_right, v_right = transformed[right_index]
        right_sign = signs[right_index]
        if left_sign == right_sign:
            continue
        u_cross = (u_left * v_right - u_right * v_left) / (v_right - v_left)
        crossing_sign = strict_sign(u_cross)
        require(crossing_sign is not None, f"ray {ray} crossing is not strict")
        if crossing_sign > 0:
            winding += 1 if left_sign < right_sign else -1
            positive_crossings += 1
    return winding, positive_crossings


def static_audit(payload: dict, core) -> dict:
    require(payload["kind"] == STEM, "result kind mismatch")
    require(payload["schema_version"] == 1, "schema version mismatch")
    require(
        payload["status"]
        == "rigorous complete-theta spatial tile exclusion certificate complete",
        "certificate is not complete",
    )
    require(payload["source_sha256"] == core.sha256_path(SOURCE), "source hash mismatch")
    for parent in payload["parents"].values():
        path = REPO_ROOT / parent["path"]
        require(parent["sha256"] == core.sha256_path(path), f"parent hash mismatch: {path}")

    resource = payload["resource_policy"]
    require(resource["active_compute_workers"] == 1, "worker count is not one")
    require(resource["thread_caps"] == 1, "thread cap is not one")
    require(resource["below_normal_priority_applied"] is True, "priority gate missing")
    require(resource["resource_parked"] is False, "production certificate parked")

    expected_tiles = {
        "left": {"time": ["0", "0.45"], "x": ["135", "135.5"], "coordinate_orientation": "(t,x)"},
        "right": {"time": ["0", "0.45"], "x": ["136", "136.5"], "coordinate_orientation": "(t,x)"},
    }
    require(payload["tiles"] == expected_tiles, "tile geometry mismatch")
    require(payload["combined_band"]["x"] == ["135", "136.5"], "combined band mismatch")
    require(payload["combined_band"]["time"] == ["0", "0.45"], "combined time mismatch")

    config = payload["interval_configuration"]
    require(config["precision_bits"] == 192, "production precision mismatch")
    require(config["absolute_integration_tolerance"] == "1e-34", "production tolerance mismatch")
    require(config["taylor_time_order"] == 10, "time Taylor order mismatch")
    require(config["taylor_x_order"] == 28, "x Taylor order mismatch")
    require(config["boundary_segments_per_edge"] == 8, "segment count mismatch")
    tail_radii = [arb(value) for value in config["complete_tail_derivative_radii"]]
    rebuilt_tail = core.load_tail_radii()
    require(len(tail_radii) == len(rebuilt_tail) == 2, "tail radius count mismatch")
    for stored, rebuilt in zip(tail_radii, rebuilt_tail, strict=True):
        require(stored.overlaps(rebuilt), "tail radius drifted")

    boundary = payload["boundary_certificate"]
    records = payload["boundary_cells"]
    require(boundary["initial_cells"] == 64, "initial cell count mismatch")
    require(boundary["evaluated_cells"] == 64, "unexpected adaptive evaluations")
    require(boundary["certified_leaf_cells"] == 64, "leaf count mismatch")
    require(boundary["queued_cells"] == 0, "production queue is not empty")
    require(boundary["maximum_depth"] == 0, "an initial cell was subdivided")
    require(boundary["all_boundary_cells_nonzero"] is True, "boundary flag failed")
    require(len(records) == 64, "stored cell count mismatch")

    expected_keys = {
        (tile_id, edge, index)
        for tile_id in core.TILES
        for edge in ("x_low", "t_high", "x_high", "t_low")
        for index in range(8)
    }
    actual_keys: set[tuple[str, str, int]] = set()
    minimum: arb | None = None
    maximum_ratio: arb | None = None
    for record in records:
        cell = core.cell_from_serialized(record["cell"])
        require(record["passed"] is True, "stored boundary cell failed")
        require(cell.depth == 0, "stored cell is subdivided")
        require(record["separated_by"] in ("F", "F_x"), "separator missing")
        separation = arb(record["separation_lower"])
        require(separation.lower() > 0, "nonpositive separation")
        fields = [arb(item["ball"]) for item in record["field_boxes"]]
        five = [arb(item["ball"]) for item in record["five_term_boxes"]]
        selected = 0 if record["separated_by"] == "F" else 1
        require(not fields[selected].contains(0), "selected complete box contains zero")
        require(fields[selected].contains(five[selected]), "complete box omits five-term box")
        tail = arb(record["tail_radii"][selected])
        stored_ratio = arb(record["selected_tail_to_margin_upper"])
        recomputed_ratio = tail.upper() / separation
        relative_drift = abs(float(stored_ratio.mid()) - float(recomputed_ratio.mid())) / float(
            recomputed_ratio.mid()
        )
        require(relative_drift < 1e-5, "stored tail ratio drifted")
        require(stored_ratio.upper() < arb("4e-21"), "tail-to-margin ratio is too large")
        minimum = separation if minimum is None or separation < minimum else minimum
        maximum_ratio = (
            stored_ratio
            if maximum_ratio is None or stored_ratio > maximum_ratio
            else maximum_ratio
        )
        low, high, _ = core.edge_coordinate_bounds(cell.tile_id, cell.edge)
        width = (high - low) / 8
        index = int((cell.coordinate_low - low) / width)
        actual_keys.add((cell.tile_id, cell.edge, index))
    require(actual_keys == expected_keys, "boundary cells do not form the exact grids")
    require(minimum is not None and maximum_ratio is not None, "boundary extrema missing")
    require(
        arb(boundary["minimum_separation_lower"]).overlaps(minimum),
        "minimum separation summary mismatch",
    )
    require(
        arb(boundary["maximum_tail_to_margin_upper"]).overlaps(maximum_ratio),
        "maximum ratio summary mismatch",
    )
    for tile_audit in boundary["coverage_audit"].values():
        for edge_audit in tile_audit.values():
            require(edge_audit["passed"] is True, "edge coverage audit failed")
            require(
                edge_audit["certified_length"] == edge_audit["expected_length"],
                "edge lengths differ",
            )

    ray_checks: dict[str, dict] = {}
    for tile_id, winding in payload["tile_windings"].items():
        require(tile_id in core.TILES, f"unknown winding tile: {tile_id}")
        require(winding["certified_segments"] == 32, "winding segment count mismatch")
        require(winding["polygon"]["winding_t_x"] == 0, "stored winding is not zero")
        vertices = winding["vertices"]
        require(len(vertices) == 32, "winding vertex count mismatch")
        vertex_fields = [
            [arb(item["ball"]) for item in vertex["field_balls"]]
            for vertex in vertices
        ]
        for index, vertex in enumerate(vertices):
            segment = [arb(item["ball"]) for item in vertex["segment_field_boxes"]]
            selected = 0 if vertex["segment_separated_by"] == "F" else 1
            require(not segment[selected].contains(0), "winding segment contains origin")
            require(
                all(segment[q].contains(vertex_fields[index][q]) for q in range(2)),
                "segment omits starting finite endpoint",
            )
            require(
                all(
                    segment[q].contains(vertex_fields[(index + 1) % len(vertices)][q])
                    for q in range(2)
                ),
                "segment omits ending finite endpoint",
            )
        tile_rays = {}
        for ray in ((1, 0), (1, 1), (2, 1)):
            value, crossings = ray_winding(vertex_fields, ray)
            require(value == 0, f"tile {tile_id} has nonzero winding on ray {ray}")
            tile_rays[str(ray)] = {
                "winding": value,
                "positive_ray_crossings": crossings,
            }
        ray_checks[tile_id] = tile_rays

    require("[0,0.45]x[135,136.5]" in payload["theorem"], "combined theorem missing")
    degree_text = json.dumps(payload["degree_composition"], sort_keys=True)
    for phrase in ("-floor(m/2)", "degree zero", "[135,136.5]"):
        require(phrase in degree_text, f"degree composition lost phrase: {phrase}")
    for phrase in ("other frequencies", "all real x", "Lambda<=0", "RH", "prize-level"):
        require(phrase in payload["proof_boundary"], f"proof boundary lost phrase: {phrase}")
    return {
        "boundary_cells": len(records),
        "minimum_separation": str(minimum),
        "maximum_tail_to_margin": str(maximum_ratio),
        "tile_ray_checks": ray_checks,
    }


def independent_replay(payload: dict, core) -> dict:
    core.continuation.ABS_TOL = "1e-38"
    core.continuation.TAYLOR_TIME_ORDER = 11
    core.continuation.TAYLOR_X_ORDER = 30
    core.continuation.flint.ctx.prec = 224
    tail_radii = core.load_tail_radii()
    monitor = core.CpuMonitor()
    anchors: list[dict] = []
    for tile_id in core.TILES:
        for edge in ("x_low", "t_high", "x_high", "t_low"):
            rows = [
                row
                for row in payload["boundary_cells"]
                if row["cell"]["tile_id"] == tile_id
                and row["cell"]["edge"] == edge
            ]
            source = min(
                rows,
                key=lambda row: float(arb(row["separation_lower"]).lower()),
            )
            cell = core.cell_from_serialized(source["cell"])
            rebuilt = core.certify_boundary_cell(cell, tail_radii)
            require(rebuilt["passed"] is True, f"independent {tile_id}:{edge} anchor failed")
            require(
                arb(rebuilt["selected_tail_to_margin_upper"]).upper() < arb("4e-21"),
                f"independent {tile_id}:{edge} tail ratio is too large",
            )
            anchors.append(
                {
                    "tile_id": tile_id,
                    "edge": edge,
                    "cell": rebuilt["cell"],
                    "separated_by": rebuilt["separated_by"],
                    "separation_lower": rebuilt["separation_lower"],
                    "tail_to_margin_upper": rebuilt["selected_tail_to_margin_upper"],
                }
            )
        monitor.sample()

    point_cache: dict[tuple[Fraction, Fraction], list[arb]] = {}
    independent_rays: dict[str, dict] = {}
    for tile_id, winding in payload["tile_windings"].items():
        vertices: list[list[arb]] = []
        for vertex in winding["vertices"]:
            point = (
                Fraction(vertex["location"]["time"]),
                Fraction(vertex["location"]["x"]),
            )
            if point not in point_cache:
                point_cache[point] = core.point_field(point[0], point[1])
            vertices.append(point_cache[point])
        tile_rays = {}
        for ray in ((1, 0), (1, 1), (2, 1)):
            value, crossings = ray_winding(vertices, ray)
            require(value == 0, f"independent tile {tile_id} winding failed on {ray}")
            tile_rays[str(ray)] = {
                "winding": value,
                "positive_ray_crossings": crossings,
            }
        independent_rays[tile_id] = tile_rays
        monitor.sample()
    return {
        "precision_bits": 224,
        "absolute_integration_tolerance": "1e-38",
        "taylor_time_order": 11,
        "taylor_x_order": 30,
        "boundary_anchors": anchors,
        "point_evaluations": len(point_cache),
        "tile_ray_checks": independent_rays,
        "runtime_cpu_percent": monitor.samples,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=RESULT)
    parser.add_argument("--baseline-seconds", type=int, default=5)
    parser.add_argument("--baseline-max-percent", type=float, default=60.0)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    core = importlib.import_module(STEM)
    core.continuation.flint.ctx.prec = 192
    priority_lowered = core.continuation.request_below_normal_priority()
    baseline = [float(psutil.cpu_percent(interval=1.0)) for _ in range(args.baseline_seconds)]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "complete-theta tile checker baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, "
        f"mean={baseline_mean:.1f}%",
        flush=True,
    )
    if baseline_mean > args.baseline_max_percent:
        print("complete-theta tile checker deferred: daytime baseline is already busy")
        return 2
    require(priority_lowered, "checker could not apply below-normal priority")
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    static = static_audit(payload, core)
    independent = independent_replay(payload, core)
    print(
        "checked complete-theta spatial tile degree certificate: "
        f"boundary_cells={static['boundary_cells']}, tile_windings=[0, 0], "
        f"independent_anchors={len(independent['boundary_anchors'])}, "
        "combined_band=[135,136.5], no_contacts=True, 0 issues"
    )
    print(json.dumps({"static": static, "independent": independent}, indent=2))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except CheckFailure as error:
        print(f"CHECK FAILED: {error}", file=sys.stderr)
        sys.exit(1)
