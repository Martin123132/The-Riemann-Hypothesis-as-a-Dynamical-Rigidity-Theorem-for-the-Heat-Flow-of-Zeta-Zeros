#!/usr/bin/env python3
"""Scout the four-term positive-time boundary from lambda=0.22 to 1."""

from __future__ import annotations

import argparse
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


STEM = "jensen_window_pf_newman_theta_fourth_summand_positive_time_boundary_scout"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
DATE = "2026-08-03"

LAMBDA_LOW = Fraction(11, 50)
LAMBDA_HIGH = Fraction(1)
LAMBDA_STEP = Fraction(1, 50)
TIME_LOW = Fraction(0)
TIME_HIGH = Fraction(9, 20)
X_LOW = Fraction(271, 2)
X_HIGH = Fraction(136)
PRECISION_BITS = 192


class ResourcePark(RuntimeError):
    """Raised between completed boundary vertices after sustained high CPU."""


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


def boundary_vertices(segments: int) -> list[tuple[str, Fraction, Fraction]]:
    vertices: list[tuple[str, Fraction, Fraction]] = []
    for time in fraction_grid(TIME_LOW, TIME_HIGH, segments)[:-1]:
        vertices.append(("x_low", time, X_LOW))
    for x in fraction_grid(X_LOW, X_HIGH, segments)[:-1]:
        vertices.append(("t_high", TIME_HIGH, x))
    for time in reversed(fraction_grid(TIME_LOW, TIME_HIGH, segments)[1:]):
        vertices.append(("x_high", time, X_HIGH))
    for x in reversed(fraction_grid(X_LOW, X_HIGH, segments)[1:]):
        vertices.append(("t_low", TIME_LOW, x))
    return vertices


def affine_point_jets(time: Fraction, x: Fraction) -> tuple[list, list]:
    base = continuation.TransformFamily(
        Fraction(0),
        time,
        time,
        x,
        weighted=True,
    )
    fourth = continuation.TransformFamily(
        Fraction(0),
        time,
        time,
        x,
        weighted=False,
    )
    return (
        [base.point_jet(order) for order in range(2)],
        [fourth.point_jet(order) for order in range(2)],
    )


def midpoint(value) -> float:
    return float(value.mid())


def radius(value) -> float:
    return float(value.rad())


def affine_value(row: dict, lambda_value: Fraction) -> complex:
    lam = float(lambda_value)
    return complex(
        row["base_midpoint"][0] + lam * row["fourth_midpoint"][0],
        row["base_midpoint"][1] + lam * row["fourth_midpoint"][1],
    )


def polygon_winding(points: list[complex]) -> tuple[int, float, float]:
    total = 0.0
    maximum_increment = 0.0
    for index, left in enumerate(points):
        right = points[(index + 1) % len(points)]
        delta = math.atan2(right.imag, right.real) - math.atan2(left.imag, left.real)
        while delta <= -math.pi:
            delta += 2 * math.pi
        while delta > math.pi:
            delta -= 2 * math.pi
        maximum_increment = max(maximum_increment, abs(delta))
        total += delta
    return round(total / (2 * math.pi)), total / (2 * math.pi), maximum_increment


def affine_vertex_minimum(row: dict) -> dict:
    a0, a1 = row["base_midpoint"]
    b0, b1 = row["fourth_midpoint"]
    denominator = b0 * b0 + b1 * b1
    if denominator == 0:
        lam = float(LAMBDA_LOW)
    else:
        lam = -(a0 * b0 + a1 * b1) / denominator
        lam = min(float(LAMBDA_HIGH), max(float(LAMBDA_LOW), lam))
    value = complex(a0 + lam * b0, a1 + lam * b1)
    return {
        "lambda": f"{lam:.17g}",
        "norm": abs(value),
        "F": f"{value.real:.17e}",
        "F_x": f"{value.imag:.17e}",
    }


def render_note(payload: dict) -> str:
    summary = payload["summary"]
    lines = [
        "# Four-Summand Positive-Time Boundary Scout",
        "",
        f"Date: {DATE}",
        "",
        "Status: finite high-precision diagnostic, not an interval certificate and not a proof of RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py --progress",
        "```",
        "",
        "## Result",
        "",
        f"The sampled winding values are `{summary['distinct_windings']}` on all",
        f"`{summary['lambda_nodes']}` coefficient nodes from `0.22` through `1`.",
        f"The weakest sampled boundary norm is `{summary['minimum_sampled_norm']:.17e}`.",
        f"The largest polygon angle increment is `{summary['maximum_angle_increment']:.8f}` radians.",
        "",
        "This scout decides whether the declared positive-time rectangle is a viable",
        "target for a parameter-uniform ACB/Taylor boundary certificate. It does not",
        "exclude an event between sampled vertices or coefficient nodes.",
        "",
        "## Proof Boundary",
        "",
        payload["proof_boundary"],
        "",
    ]
    return "\n".join(lines)


def build(args: argparse.Namespace, baseline: list[float], monitor: CpuMonitor) -> dict:
    started = perf_counter()
    vertices = boundary_vertices(args.segments_per_edge)
    rows: list[dict] = []
    for index, (edge, time, x) in enumerate(vertices):
        base, fourth = affine_point_jets(time, x)
        row = {
            "index": index,
            "edge": edge,
            "time": fraction_text(time),
            "x": fraction_text(x),
            "base_midpoint": [midpoint(value) for value in base],
            "base_radius": [radius(value) for value in base],
            "fourth_midpoint": [midpoint(value) for value in fourth],
            "fourth_radius": [radius(value) for value in fourth],
        }
        row["affine_lambda_minimum"] = affine_vertex_minimum(row)
        rows.append(row)
        if (index + 1) % args.checkpoint_every == 0:
            partial = {
                "kind": f"{STEM}_partial",
                "source_sha256": sha256_path(Path(__file__).resolve()),
                "segments_per_edge": args.segments_per_edge,
                "completed_vertices": len(rows),
                "vertices": rows,
            }
            continuation.write_json_atomic(args.partial, partial)
            if args.progress:
                print(
                    f"positive-time scout: vertices={len(rows)}/{len(vertices)}",
                    flush=True,
                )
            monitor.sample()

    lambda_nodes = [Fraction(index, 50) for index in range(11, 51)]
    lambda_rows: list[dict] = []
    global_minimum: dict | None = None
    maximum_increment = 0.0
    for lambda_value in lambda_nodes:
        points = [affine_value(row, lambda_value) for row in rows]
        winding, raw_winding, angle_increment = polygon_winding(points)
        weakest_index = min(range(len(points)), key=lambda index: abs(points[index]))
        weakest_value = points[weakest_index]
        weakest_row = rows[weakest_index]
        candidate = {
            "lambda": fraction_text(lambda_value),
            "vertex_index": weakest_index,
            "edge": weakest_row["edge"],
            "time": weakest_row["time"],
            "x": weakest_row["x"],
            "F": f"{weakest_value.real:.17e}",
            "F_x": f"{weakest_value.imag:.17e}",
            "norm": abs(weakest_value),
        }
        if global_minimum is None or candidate["norm"] < global_minimum["norm"]:
            global_minimum = candidate
        maximum_increment = max(maximum_increment, angle_increment)
        lambda_rows.append(
            {
                "lambda": fraction_text(lambda_value),
                "winding_t_x": winding,
                "raw_winding": raw_winding,
                "maximum_angle_increment": angle_increment,
                "weakest_vertex": candidate,
            }
        )

    affine_minimum_index = min(
        range(len(rows)),
        key=lambda index: rows[index]["affine_lambda_minimum"]["norm"],
    )
    affine_minimum = {
        "vertex_index": affine_minimum_index,
        "edge": rows[affine_minimum_index]["edge"],
        "time": rows[affine_minimum_index]["time"],
        "x": rows[affine_minimum_index]["x"],
        **rows[affine_minimum_index]["affine_lambda_minimum"],
    }
    distinct_windings = sorted({row["winding_t_x"] for row in lambda_rows})
    require(global_minimum is not None, "sampled boundary is empty")
    return {
        "kind": STEM,
        "schema_version": 1,
        "date": DATE,
        "status": "finite positive-time boundary scout complete",
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "resource_policy": {
            "mode": "daytime",
            "active_compute_workers": 1,
            "thread_caps": 1,
            "below_normal_priority_applied": True,
            "baseline_cpu_percent": baseline,
            "baseline_mean_percent": sum(baseline) / len(baseline),
            "runtime_cpu_percent": monitor.samples,
            "elapsed_seconds": perf_counter() - started,
        },
        "rectangle": {
            "lambda": [fraction_text(LAMBDA_LOW), fraction_text(LAMBDA_HIGH)],
            "time": [fraction_text(TIME_LOW), fraction_text(TIME_HIGH)],
            "x": [fraction_text(X_LOW), fraction_text(X_HIGH)],
            "coordinate_orientation": "(t,x)",
        },
        "configuration": {
            "precision_bits": PRECISION_BITS,
            "absolute_integration_tolerance": continuation.ABS_TOL,
            "segments_per_edge": args.segments_per_edge,
            "boundary_vertices": len(rows),
            "lambda_step": fraction_text(LAMBDA_STEP),
        },
        "summary": {
            "lambda_nodes": len(lambda_rows),
            "distinct_windings": distinct_windings,
            "winding_constant_zero": distinct_windings == [0],
            "minimum_sampled_norm": global_minimum["norm"],
            "minimum_sampled_point": global_minimum,
            "affine_vertex_minimum": affine_minimum,
            "maximum_angle_increment": maximum_increment,
            "angle_resolution_below_pi_over_two": maximum_increment < math.pi / 2,
        },
        "lambda_rows": lambda_rows,
        "affine_boundary_vertices": rows,
        "route_decision": (
            "If every sampled winding is zero, the angle increments are resolved, and "
            "the boundary norm remains separated, promote this rectangle to an adaptive "
            "parameter-uniform ACB/Taylor lateral-boundary certificate. Otherwise isolate "
            "the indicated edge event before changing the rectangle."
        ),
        "proof_boundary": (
            "Finite 192-bit point diagnostic on a sampled boundary and forty coefficient "
            "nodes only. It is not interval coverage between vertices or lambda nodes, not "
            "a Brouwer-degree theorem, not control outside the declared rectangle, not a "
            "fifth-summand or complete-Xi result, and not a proof of Lambda<=0 or RH."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--partial", type=Path, default=REPO_ROOT / "work/rh_compute/results" / f"{STEM}_partial.json")
    parser.add_argument("--segments-per-edge", type=int, default=32)
    parser.add_argument("--checkpoint-every", type=int, default=8)
    parser.add_argument("--baseline-seconds", type=int, default=5)
    parser.add_argument("--baseline-max-percent", type=float, default=60.0)
    parser.add_argument("--progress", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    require(args.segments_per_edge >= 8, "at least eight segments per edge are required")
    require(args.checkpoint_every >= 1, "checkpoint interval must be positive")
    priority_lowered = continuation.request_below_normal_priority()
    baseline = [float(psutil.cpu_percent(interval=1.0)) for _ in range(args.baseline_seconds)]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "positive-time boundary scout baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, "
        f"mean={baseline_mean:.1f}%",
        flush=True,
    )
    if baseline_mean > args.baseline_max_percent:
        print("positive-time boundary scout deferred: daytime baseline is already busy")
        return 2
    require(priority_lowered, "could not apply below-normal process priority")
    continuation.flint.ctx.prec = PRECISION_BITS
    monitor = CpuMonitor()
    try:
        payload = build(args, baseline, monitor)
    except ResourcePark as exc:
        print(f"positive-time boundary scout parked: {exc}")
        return 3
    continuation.write_json_atomic(args.out, payload)
    continuation.write_text_atomic(args.note, render_note(payload))
    if args.partial.exists():
        args.partial.unlink()
    print(
        "built positive-time boundary scout: "
        f"vertices={payload['configuration']['boundary_vertices']}, "
        f"lambda_nodes={payload['summary']['lambda_nodes']}, "
        f"windings={payload['summary']['distinct_windings']}, "
        f"minimum_norm={payload['summary']['minimum_sampled_norm']:.6e}, "
        f"elapsed={payload['resource_policy']['elapsed_seconds']:.3f}s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
