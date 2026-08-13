#!/usr/bin/env python3
"""Check the four-summand compact contact degree certificate."""

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


STEM = "jensen_window_pf_newman_theta_fourth_summand_compact_contact_degree_certificate"
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
    transformed: list[tuple[arb, arb]] = []
    for fields in vertices:
        transformed.append(
            (a * fields[0] + b * fields[1], -b * fields[0] + a * fields[1])
        )
    signs = [strict_sign(pair[1]) for pair in transformed]
    require(all(sign is not None for sign in signs), f"ray {ray} hits a vertex enclosure")
    winding = 0
    positive_crossings = 0
    for index, ((u_left, v_left), left_sign) in enumerate(zip(transformed, signs, strict=True)):
        right_index = (index + 1) % len(vertices)
        u_right, v_right = transformed[right_index]
        right_sign = signs[right_index]
        if left_sign == right_sign:
            continue
        u_cross = (
            (u_left * v_right - u_right * v_left) / (v_right - v_left)
        )
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
        payload["status"] == "rigorous compact contact degree certificate complete",
        "certificate is not complete",
    )
    require(payload["source_sha256"] == core.sha256_path(SOURCE), "source hash mismatch")
    for parent in payload["parents"].values():
        path = REPO_ROOT / parent["path"]
        require(parent["sha256"] == core.sha256_path(path), f"parent hash mismatch: {path}")

    resource = payload["resource_policy"]
    require(resource["active_compute_workers"] == 1, "worker count is not one")
    require(resource["thread_caps"] == 1, "thread cap is not one")
    require(resource["below_normal_priority_applied"] is True, "below-normal priority missing")
    require(resource["resource_parked"] is False, "production certificate parked")

    require(payload["cylinder"]["lambda"] == ["0", "0.22"], "lambda cylinder mismatch")
    require(payload["cylinder"]["time"] == ["-0.06", "0.45"], "time cylinder mismatch")
    require(payload["cylinder"]["x"] == ["135.5", "136"], "x cylinder mismatch")
    require(payload["cylinder"]["coordinate_orientation"] == "(t,x)", "orientation mismatch")

    boundary = payload["boundary_certificate"]
    records = payload["boundary_cells"]
    require(boundary["initial_cells"] == 352, "initial boundary count mismatch")
    require(boundary["evaluated_cells"] == 352, "unexpected adaptive evaluations")
    require(boundary["certified_leaf_cells"] == 352, "leaf count mismatch")
    require(boundary["maximum_depth"] == 0, "unexpected boundary subdivision")
    require(boundary["all_boundary_cells_nonzero"] is True, "boundary nonzero flag failed")
    require(len(records) == 352, "stored boundary cell count mismatch")
    for audit in boundary["edge_area_audit"].values():
        require(audit["passed"] is True, "edge area audit failed")
        require(audit["certified_area"] == audit["expected_area"], "edge areas differ")

    expected_keys = {
        (edge, Fraction(k, 50), Fraction(k + 1, 50), coordinate_index)
        for edge in ("x_low", "t_high", "x_high", "t_low")
        for k in range(11)
        for coordinate_index in range(8)
    }
    actual_keys: set[tuple] = set()
    minimum_separation: arb | None = None
    edge_bounds = {
        "x_low": (Fraction(-3, 50), Fraction(9, 20)),
        "x_high": (Fraction(-3, 50), Fraction(9, 20)),
        "t_low": (Fraction(271, 2), Fraction(136)),
        "t_high": (Fraction(271, 2), Fraction(136)),
    }
    for record in records:
        cell = record["cell"]
        require(record["passed"] is True, "stored boundary cell failed")
        require(record["separated_by"] in ("F", "F_x"), "missing separating component")
        separation = arb(record["separation_lower"])
        require(separation.lower() > 0, "nonpositive boundary separation")
        if minimum_separation is None or separation < minimum_separation:
            minimum_separation = separation
        field_boxes = [arb(item["ball"]) for item in record["field_boxes"]]
        selected = 0 if record["separated_by"] == "F" else 1
        require(not field_boxes[selected].contains(0), "selected field box contains zero")
        low, high = edge_bounds[cell["edge"]]
        width = (high - low) / 8
        coordinate_low = Fraction(cell["coordinate_low"])
        coordinate_index = int((coordinate_low - low) / width)
        actual_keys.add(
            (
                cell["edge"],
                Fraction(cell["lambda_low"]),
                Fraction(cell["lambda_high"]),
                coordinate_index,
            )
        )
    require(actual_keys == expected_keys, "boundary cells do not form the exact initial grid")
    require(
        arb(boundary["minimum_separation_lower"]).overlaps(minimum_separation),
        "minimum separation summary mismatch",
    )

    winding = payload["lambda_zero_winding"]
    require(winding["lambda"] == "0", "winding lambda mismatch")
    require(winding["certified_segments"] == 32, "winding segment count mismatch")
    require(winding["polygon"]["winding_t_x"] == -1, "stored winding is not -1")
    vertices = winding["vertices"]
    require(len(vertices) == 32, "winding vertex count mismatch")
    vertex_fields = [
        [arb(item["ball"]) for item in vertex["field_balls"]]
        for vertex in vertices
    ]
    for index, vertex in enumerate(vertices):
        segment_fields = [arb(item["ball"]) for item in vertex["segment_field_boxes"]]
        selected = 0 if vertex["segment_separated_by"] == "F" else 1
        require(not segment_fields[selected].contains(0), "winding segment contains origin")
        next_fields = vertex_fields[(index + 1) % len(vertices)]
        require(
            all(segment_fields[q].contains(vertex_fields[index][q]) for q in range(2)),
            "segment omits its starting endpoint",
        )
        require(
            all(segment_fields[q].contains(next_fields[q]) for q in range(2)),
            "segment omits its ending endpoint",
        )
    ray_checks = {}
    for ray in ((1, 0), (1, 1), (2, 1)):
        value, crossings = ray_winding(vertex_fields, ray)
        require(value == -1, f"stored polygon has wrong winding on ray {ray}")
        ray_checks[str(ray)] = {"winding": value, "positive_ray_crossings": crossings}

    branch = payload["branch_containment"]
    require(branch["charts"] == 11, "branch chart count mismatch")
    require(branch["all_charts_inside_cylinder"] is True, "branch containment failed")
    require(branch["all_branch_contacts_regular"] is True, "branch regularity failed")
    require(Fraction(branch["minimum_spatial_margin"]) > 0, "branch has no spatial margin")
    require("unique common zero" in payload["theorem"], "complete-count theorem missing")
    for phrase in ("outside that cylinder", "n>=5", "RH", "prize-level"):
        require(phrase in payload["proof_boundary"], f"proof boundary lost phrase: {phrase}")
    return {
        "boundary_cells": len(records),
        "minimum_separation": str(minimum_separation),
        "winding_segments": len(vertices),
        "ray_checks": ray_checks,
        "branch_charts": branch["charts"],
    }


def independent_replay(payload: dict, core) -> dict:
    core.continuation.ABS_TOL = "1e-38"
    core.continuation.TAYLOR_TIME_ORDER = 11
    core.continuation.TAYLOR_X_ORDER = 30
    core.continuation.flint.ctx.prec = 224
    monitor = core.CpuMonitor()
    anchors: list[dict] = []
    for edge in ("x_low", "t_high", "x_high", "t_low"):
        edge_rows = [row for row in payload["boundary_cells"] if row["cell"]["edge"] == edge]
        source = min(edge_rows, key=lambda row: float(arb(row["separation_lower"]).lower()))
        cell = core.cell_from_serialized(source["cell"])
        rebuilt = core.certify_boundary_cell(cell)
        require(rebuilt["passed"] is True, f"independent {edge} anchor failed")
        anchors.append(
            {
                "edge": edge,
                "cell": rebuilt["cell"],
                "separated_by": rebuilt["separated_by"],
                "separation_lower": rebuilt["separation_lower"],
            }
        )
        monitor.sample()

    point_vertices: list[list[arb]] = []
    for vertex in payload["lambda_zero_winding"]["vertices"]:
        time = Fraction(vertex["location"]["time"])
        x = Fraction(vertex["location"]["x"])
        point_vertices.append(core.point_field(Fraction(0), time, x)[:2])
    alternate_rays = {}
    for ray in ((1, 0), (1, 1), (2, 1)):
        winding, crossings = ray_winding(point_vertices, ray)
        require(winding == -1, f"independent point winding failed on ray {ray}")
        alternate_rays[str(ray)] = {
            "winding": winding,
            "positive_ray_crossings": crossings,
        }
    monitor.sample()
    return {
        "precision_bits": 224,
        "absolute_integration_tolerance": "1e-38",
        "taylor_time_order": 11,
        "taylor_x_order": 30,
        "boundary_anchors": anchors,
        "point_polygon_ray_checks": alternate_rays,
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
    priority_lowered = core.continuation.request_below_normal_priority()
    baseline = [float(psutil.cpu_percent(interval=1.0)) for _ in range(args.baseline_seconds)]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "compact degree checker baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, mean={baseline_mean:.1f}%",
        flush=True,
    )
    if baseline_mean > args.baseline_max_percent:
        print("compact degree checker deferred: daytime baseline is already busy")
        return 2
    require(priority_lowered, "checker could not apply below-normal priority")
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    static = static_audit(payload, core)
    independent = independent_replay(payload, core)
    print(
        "checked compact contact degree certificate: "
        f"boundary_cells={static['boundary_cells']}, winding=-1, "
        f"independent_anchors={len(independent['boundary_anchors'])}, "
        "unique_contact=True, 0 issues"
    )
    print(json.dumps({"static": static, "independent": independent}, indent=2))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except CheckFailure as error:
        print(f"CHECK FAILED: {error}", file=sys.stderr)
        sys.exit(1)
