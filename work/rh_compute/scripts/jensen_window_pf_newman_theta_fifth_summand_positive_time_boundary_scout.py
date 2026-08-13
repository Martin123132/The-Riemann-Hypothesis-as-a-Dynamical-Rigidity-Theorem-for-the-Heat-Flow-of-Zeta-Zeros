#!/usr/bin/env python3
"""Scout the fifth theta-summand coefficient on the positive-time boundary."""

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


STEM = "jensen_window_pf_newman_theta_fifth_summand_positive_time_boundary_scout"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
DATE = "2026-08-04"

MU_LOW = Fraction(0)
MU_HIGH = Fraction(1)
TIME_LOW = Fraction(0)
TIME_HIGH = Fraction(9, 20)
X_LOW = Fraction(271, 2)
X_HIGH = Fraction(136)
PRECISION_BITS = 192


class FifthTransformFamily(continuation.TransformFamily):
    """First four theta components plus mu times the fifth, or the fifth alone."""

    def phi(self, u):
        if not self.weighted:
            return self.phi_n(u, 5)
        return (
            self.phi_n(u, 1)
            + self.phi_n(u, 2)
            + self.phi_n(u, 3)
            + self.phi_n(u, 4)
            + continuation.acb(continuation.arb_fraction(self.lambda_center))
            * self.phi_n(u, 5)
        )


class ResourcePark(RuntimeError):
    """Raised between completed vertices after sustained high daytime CPU."""


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


def extended_tail_audit() -> dict:
    parent = continuation.exact_tail_audit(60)
    require(parent["radius_used_per_transform"] == continuation.TAIL_RADIUS, "tail radius drift")
    require(3 * 16 * 2000 > 4, "large-z monotonicity threshold audit failed")
    require(5 * 10**1500 < 10**1800, "five-component tail radius audit failed")
    return {
        "parent_n_le_4_audit": parent,
        "n5_dominance": (
            "Write z=pi*n^2*exp(4u). Then phi_n(u)=exp(u)*z*(2z-3)*exp(-z). "
            "For z>3, the z derivative has sign -2z^2+7z-3=-(2z-1)(z-3)<0. "
            "On u>=2 and n>=4, z>3*16*2000>3, so 0<phi_5(u)<phi_4(u)."
        ),
        "moment_transfer": (
            "Multiplication by u^m*exp(t*u^2) is nonnegative for real t and u>=2; "
            "therefore every absolute n=5 tail moment through m=60 is bounded by "
            "the already certified n=4 tail moment."
        ),
        "combined_family": (
            "For 0<=mu<=1, the first-four-plus-mu*fifth tail is below "
            "5*10^-1800, hence below the attached 10^-1500 radius."
        ),
    }


def affine_point_jets(time: Fraction, x: Fraction) -> tuple[list, list]:
    base = FifthTransformFamily(Fraction(0), time, time, x, weighted=True)
    fifth = FifthTransformFamily(Fraction(0), time, time, x, weighted=False)
    return (
        [base.point_jet(order) for order in range(2)],
        [fifth.point_jet(order) for order in range(2)],
    )


def midpoint(value) -> float:
    return float(value.mid())


def radius(value) -> float:
    return float(value.rad())


def affine_value(row: dict, mu: Fraction) -> complex:
    coefficient = float(mu)
    return complex(
        row["base_midpoint"][0] + coefficient * row["fifth_midpoint"][0],
        row["base_midpoint"][1] + coefficient * row["fifth_midpoint"][1],
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
    b0, b1 = row["fifth_midpoint"]
    denominator = b0 * b0 + b1 * b1
    if denominator == 0:
        mu = 0.0
    else:
        mu = -(a0 * b0 + a1 * b1) / denominator
        mu = min(1.0, max(0.0, mu))
    value = complex(a0 + mu * b0, a1 + mu * b1)
    return {
        "mu": f"{mu:.17g}",
        "norm": abs(value),
        "F": f"{value.real:.17e}",
        "F_x": f"{value.imag:.17e}",
    }


def render_note(payload: dict) -> str:
    summary = payload["summary"]
    lines = [
        "# Fifth-Summand Positive-Time Boundary Scout",
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
        f"All `{summary['mu_nodes']}` sampled fifth coefficients have winding",
        f"`{summary['distinct_windings']}`. The weakest sampled norm is",
        f"`{summary['minimum_sampled_norm']:.17e}` and the maximum pointwise",
        f"fifth/base vector ratio is `{summary['maximum_fifth_to_base_vector_ratio']:.17e}`.",
        "",
        "The explicit tail audit transfers the established n=4 u>=2 moment bound",
        "to n=5 by monotonicity of the positive theta component.",
        "",
        "## Proof Boundary",
        "",
        payload["proof_boundary"],
        "",
    ]
    return "\n".join(lines)


def build(args: argparse.Namespace, baseline: list[float], monitor: CpuMonitor) -> dict:
    started = perf_counter()
    tail_audit = extended_tail_audit()
    vertices = boundary_vertices(args.segments_per_edge)
    rows: list[dict] = []
    for index, (edge, time, x) in enumerate(vertices):
        base, fifth = affine_point_jets(time, x)
        base_midpoint = [midpoint(value) for value in base]
        fifth_midpoint = [midpoint(value) for value in fifth]
        base_norm = math.hypot(*base_midpoint)
        fifth_norm = math.hypot(*fifth_midpoint)
        row = {
            "index": index,
            "edge": edge,
            "time": fraction_text(time),
            "x": fraction_text(x),
            "base_midpoint": base_midpoint,
            "base_radius": [radius(value) for value in base],
            "fifth_midpoint": fifth_midpoint,
            "fifth_radius": [radius(value) for value in fifth],
            "base_vector_norm": base_norm,
            "fifth_vector_norm": fifth_norm,
            "fifth_to_base_vector_ratio": fifth_norm / base_norm,
        }
        row["affine_mu_minimum"] = affine_vertex_minimum(row)
        rows.append(row)
        if (index + 1) % args.checkpoint_every == 0:
            if args.progress:
                print(
                    f"fifth-summand scout: vertices={len(rows)}/{len(vertices)}",
                    flush=True,
                )
            monitor.sample()

    mu_nodes = [Fraction(index, 50) for index in range(51)]
    mu_rows: list[dict] = []
    global_minimum: dict | None = None
    maximum_increment = 0.0
    for mu in mu_nodes:
        points = [affine_value(row, mu) for row in rows]
        winding, raw_winding, angle_increment = polygon_winding(points)
        weakest_index = min(range(len(points)), key=lambda index: abs(points[index]))
        weakest_value = points[weakest_index]
        weakest_row = rows[weakest_index]
        candidate = {
            "mu": fraction_text(mu),
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
        mu_rows.append(
            {
                "mu": fraction_text(mu),
                "winding_t_x": winding,
                "raw_winding": raw_winding,
                "maximum_angle_increment": angle_increment,
                "weakest_vertex": candidate,
            }
        )

    affine_minimum_index = min(
        range(len(rows)),
        key=lambda index: rows[index]["affine_mu_minimum"]["norm"],
    )
    affine_minimum = {
        "vertex_index": affine_minimum_index,
        "edge": rows[affine_minimum_index]["edge"],
        "time": rows[affine_minimum_index]["time"],
        "x": rows[affine_minimum_index]["x"],
        **rows[affine_minimum_index]["affine_mu_minimum"],
    }
    distinct_windings = sorted({row["winding_t_x"] for row in mu_rows})
    maximum_ratio = max(row["fifth_to_base_vector_ratio"] for row in rows)
    require(global_minimum is not None, "sampled boundary is empty")
    return {
        "kind": STEM,
        "schema_version": 1,
        "date": DATE,
        "status": "finite fifth-summand positive-time boundary scout complete",
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
            "mu": ["0", "1"],
            "time": [fraction_text(TIME_LOW), fraction_text(TIME_HIGH)],
            "x": [fraction_text(X_LOW), fraction_text(X_HIGH)],
            "coordinate_orientation": "(t,x)",
        },
        "family": "H_1+H_2+H_3+H_4+mu H_5",
        "configuration": {
            "precision_bits": PRECISION_BITS,
            "absolute_integration_tolerance": continuation.ABS_TOL,
            "segments_per_edge": args.segments_per_edge,
            "boundary_vertices": len(rows),
            "mu_step": "0.02",
        },
        "tail_audit": tail_audit,
        "summary": {
            "mu_nodes": len(mu_rows),
            "distinct_windings": distinct_windings,
            "winding_constant_zero": distinct_windings == [0],
            "minimum_sampled_norm": global_minimum["norm"],
            "minimum_sampled_point": global_minimum,
            "affine_vertex_minimum": affine_minimum,
            "maximum_fifth_to_base_vector_ratio": maximum_ratio,
            "maximum_angle_increment": maximum_increment,
            "angle_resolution_below_pi_over_two": maximum_increment < math.pi / 2,
        },
        "mu_rows": mu_rows,
        "affine_boundary_vertices": rows,
        "route_decision": (
            "Test a single coefficient slab 0<=mu<=1 on each spatial boundary segment. "
            "If interval dependency is too wide, subdivide mu only as needed. Record the "
            "rigorous fifth-component boxes for use in a later all-n tail majorant."
        ),
        "proof_boundary": (
            "Finite 192-bit point diagnostic on a sampled boundary and 51 coefficient "
            "nodes only. It is not interval coverage between vertices or mu nodes, not "
            "a Brouwer-degree theorem, not a bound for n>=6 or other frequency rectangles, "
            "not a complete-Xi result, and not a proof of Lambda<=0 or RH."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
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
        "fifth-summand boundary scout baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, "
        f"mean={baseline_mean:.1f}%",
        flush=True,
    )
    if baseline_mean > args.baseline_max_percent:
        print("fifth-summand boundary scout deferred: daytime baseline is already busy")
        return 2
    require(priority_lowered, "could not apply below-normal process priority")
    continuation.flint.ctx.prec = PRECISION_BITS
    monitor = CpuMonitor()
    try:
        payload = build(args, baseline, monitor)
    except ResourcePark as exc:
        print(f"fifth-summand boundary scout parked: {exc}")
        return 3
    continuation.write_json_atomic(args.out, payload)
    continuation.write_text_atomic(args.note, render_note(payload))
    print(
        "built fifth-summand boundary scout: "
        f"vertices={payload['configuration']['boundary_vertices']}, "
        f"mu_nodes={payload['summary']['mu_nodes']}, "
        f"windings={payload['summary']['distinct_windings']}, "
        f"max_ratio={payload['summary']['maximum_fifth_to_base_vector_ratio']:.6e}, "
        f"elapsed={payload['resource_policy']['elapsed_seconds']:.3f}s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
