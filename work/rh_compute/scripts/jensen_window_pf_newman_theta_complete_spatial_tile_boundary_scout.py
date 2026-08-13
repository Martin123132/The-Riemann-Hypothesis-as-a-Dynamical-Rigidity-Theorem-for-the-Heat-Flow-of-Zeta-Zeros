#!/usr/bin/env python3
"""Scout complete-theta contact degree on reusable positive-time x tiles."""

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

import jensen_window_pf_newman_theta_fifth_summand_positive_time_boundary_scout as fifth  # noqa: E402


continuation = fifth.continuation
STEM = "jensen_window_pf_newman_theta_complete_spatial_tile_boundary_scout"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
DEFAULT_PARTIAL = REPO_ROOT / "work/rh_compute/results" / f"{STEM}_partial.json"
FIFTH_PARENT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_fifth_summand_"
    "positive_time_degree_certificate.json"
)
TAIL_PARENT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_complete_tail_"
    "positive_time_degree_certificate.json"
)

DATE = "2026-08-04"
TIME_LOW = Fraction(0)
TIME_HIGH = Fraction(9, 20)
PRECISION_BITS = 192
DEFAULT_TILES = ("135:135.5", "135.5:136", "136:136.5")


class ResourcePark(RuntimeError):
    """Raised between completed point evaluations after sustained high CPU."""


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


def fraction_key(time: Fraction, x: Fraction) -> str:
    return f"{time.numerator}/{time.denominator}:{x.numerator}/{x.denominator}"


def fraction_grid(low: Fraction, high: Fraction, segments: int) -> list[Fraction]:
    return [low + (high - low) * index / segments for index in range(segments + 1)]


def parse_tile(text: str) -> tuple[Fraction, Fraction]:
    parts = text.split(":")
    require(len(parts) == 2, f"tile must have X_LOW:X_HIGH form: {text!r}")
    low, high = (Fraction(part.strip()) for part in parts)
    require(low < high, f"tile has nonpositive width: {text!r}")
    return low, high


def boundary_vertices(
    x_low: Fraction,
    x_high: Fraction,
    segments: int,
) -> list[tuple[str, Fraction, Fraction]]:
    times = fraction_grid(TIME_LOW, TIME_HIGH, segments)
    xs = fraction_grid(x_low, x_high, segments)
    vertices: list[tuple[str, Fraction, Fraction]] = []
    for time in times[:-1]:
        vertices.append(("x_low", time, x_low))
    for x in xs[:-1]:
        vertices.append(("t_high", TIME_HIGH, x))
    for time in reversed(times[1:]):
        vertices.append(("x_high", time, x_high))
    for x in reversed(xs[1:]):
        vertices.append(("t_low", TIME_LOW, x))
    return vertices


def point_jets(time: Fraction, x: Fraction) -> list:
    family = fifth.FifthTransformFamily(
        Fraction(1),
        time,
        time,
        x,
        weighted=True,
    )
    return [family.point_jet(order) for order in range(2)]


def midpoint(value) -> float:
    return float(value.mid())


def radius(value) -> float:
    return float(value.rad())


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


def component_abs_lower(midpoint_value: float, total_radius: float) -> float:
    return max(0.0, abs(midpoint_value) - total_radius)


def load_parents() -> tuple[dict, dict, list[float]]:
    finite = json.loads(FIFTH_PARENT.read_text(encoding="utf-8"))
    tail = json.loads(TAIL_PARENT.read_text(encoding="utf-8"))
    require(
        finite["status"]
        == "rigorous fifth-summand positive-time exclusion certificate complete",
        "fifth-summand parent is incomplete",
    )
    require(
        tail["status"]
        == "rigorous complete-theta local positive-time exclusion complete",
        "complete-tail parent is incomplete",
    )
    require(tail["rectangle"]["time"] == ["0", "0.45"], "tail time range drifted")
    tail_bounds = [
        math.nextafter(float(continuation.arb(value).upper()), math.inf)
        for value in tail["positive_moment_budget"]["tail_derivative_bounds"]
    ]
    require(len(tail_bounds) == 2, "complete-tail derivative budget is incomplete")
    require(all(value > 0 for value in tail_bounds), "complete-tail budget is not positive")
    return finite, tail, tail_bounds


def partial_payload(
    tiles: list[tuple[Fraction, Fraction]],
    segments: int,
    point_rows: dict[str, dict],
) -> dict:
    return {
        "kind": f"{STEM}_partial",
        "schema_version": 1,
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "configuration": {
            "tiles": [
                [fraction_text(x_low), fraction_text(x_high)]
                for x_low, x_high in tiles
            ],
            "segments_per_edge": segments,
            "precision_bits": PRECISION_BITS,
            "absolute_integration_tolerance": continuation.ABS_TOL,
            "fifth_parent_sha256": sha256_path(FIFTH_PARENT),
            "tail_parent_sha256": sha256_path(TAIL_PARENT),
        },
        "evaluated_points": len(point_rows),
        "point_cache": point_rows,
    }


def load_partial(
    path: Path,
    tiles: list[tuple[Fraction, Fraction]],
    segments: int,
    ordered_points: dict[str, tuple[Fraction, Fraction]],
) -> dict[str, dict]:
    require(path.exists(), "--resume requested but partial point cache is missing")
    payload = json.loads(path.read_text(encoding="utf-8"))
    require(payload["kind"] == f"{STEM}_partial", "partial cache kind mismatch")
    require(payload["schema_version"] == 1, "partial cache schema mismatch")
    require(
        payload["source_sha256"] == sha256_path(Path(__file__).resolve()),
        "partial cache source hash mismatch",
    )
    expected = partial_payload(tiles, segments, {})["configuration"]
    require(payload["configuration"] == expected, "partial cache configuration drifted")
    rows = dict(payload["point_cache"])
    require(payload["evaluated_points"] == len(rows), "partial point count drifted")
    require(set(rows).issubset(ordered_points), "partial cache contains an unknown point")
    for key, row in rows.items():
        time, x = ordered_points[key]
        require(row["time"] == fraction_text(time), f"partial time drifted at {key}")
        require(row["x"] == fraction_text(x), f"partial x drifted at {key}")
    return rows


def render_note(payload: dict) -> str:
    lines = [
        "# Complete-Theta Positive-Time Spatial Tile Boundary Scout",
        "",
        f"Date: {DATE}",
        "",
        "Status: finite high-precision boundary diagnostic, not an interval tile certificate and not a proof of RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py --progress",
        "```",
        "",
        "## Result",
        "",
        f"The run requested `{payload['cache']['requested_vertex_references']}` boundary",
        f"vertices and evaluated `{payload['cache']['unique_point_evaluations']}` unique",
        f"points, reusing `{payload['cache']['cache_hits']}` shared-edge references.",
        "",
    ]
    for tile in payload["tiles"]:
        summary = tile["summary"]
        lines.extend(
            [
                f"- `{tile['id']}` on x=`{tile['x'][0]}` to `{tile['x'][1]}`: "
                f"winding `{summary['winding_t_x']}`, minimum sampled complete norm "
                f"lower `{summary['minimum_complete_norm_lower']:.17e}`, maximum angle "
                f"increment `{summary['maximum_angle_increment']:.6f}`.",
            ]
        )
    lines.extend(
        [
            "",
            "The rigorous n>=6 derivative budget is applied at every sampled point,",
            "but no claim is made between adjacent boundary vertices.",
            "",
            "## Proof Boundary",
            "",
            payload["proof_boundary"],
            "",
        ]
    )
    return "\n".join(lines)


def build(args: argparse.Namespace, baseline: list[float], monitor: CpuMonitor) -> dict:
    started = perf_counter()
    finite_parent, tail_parent, tail_bounds = load_parents()
    tiles = [parse_tile(text) for text in (args.tile or DEFAULT_TILES)]
    require(len(set(tiles)) == len(tiles), "duplicate tile requested")
    tiles.sort()

    tile_vertices: list[list[tuple[str, Fraction, Fraction]]] = []
    ordered_points: dict[str, tuple[Fraction, Fraction]] = {}
    requested_references = 0
    for x_low, x_high in tiles:
        vertices = boundary_vertices(x_low, x_high, args.segments_per_edge)
        tile_vertices.append(vertices)
        requested_references += len(vertices)
        for _, time, x in vertices:
            ordered_points.setdefault(fraction_key(time, x), (time, x))

    point_rows = (
        load_partial(
            args.partial,
            tiles,
            args.segments_per_edge,
            ordered_points,
        )
        if args.resume
        else {}
    )
    tail_vector_upper = math.hypot(*tail_bounds)
    for key, (time, x) in ordered_points.items():
        if key in point_rows:
            continue
        jets = point_jets(time, x)
        mids = [midpoint(value) for value in jets]
        finite_radii = [radius(value) for value in jets]
        complete_radii = [
            finite_radius + tail_bound
            for finite_radius, tail_bound in zip(finite_radii, tail_bounds, strict=True)
        ]
        component_lowers = [
            component_abs_lower(mid, total_radius)
            for mid, total_radius in zip(mids, complete_radii, strict=True)
        ]
        five_norm = math.hypot(*mids)
        require(five_norm > 0, f"five-term point vector vanished at {key}")
        complete_norm_lower = math.hypot(*component_lowers)
        point_rows[key] = {
            "time": fraction_text(time),
            "x": fraction_text(x),
            "five_term_midpoint": mids,
            "five_term_radius": finite_radii,
            "complete_tail_abs_upper": tail_bounds,
            "complete_component_radius": complete_radii,
            "complete_component_abs_lower": component_lowers,
            "five_term_vector_norm": five_norm,
            "complete_vector_norm_lower": complete_norm_lower,
            "tail_vector_to_five_term_norm_upper": tail_vector_upper / five_norm,
            "complete_vertex_excludes_origin": complete_norm_lower > 0,
        }
        if len(point_rows) % args.checkpoint_every == 0:
            continuation.write_json_atomic(
                args.partial,
                partial_payload(tiles, args.segments_per_edge, point_rows),
            )
            if args.progress:
                print(
                    "complete-theta tile scout: "
                    f"unique_points={len(point_rows)}/{len(ordered_points)}",
                    flush=True,
                )
            monitor.sample()

    continuation.write_json_atomic(
        args.partial,
        partial_payload(tiles, args.segments_per_edge, point_rows),
    )

    tile_rows: list[dict] = []
    for tile_index, ((x_low, x_high), vertices) in enumerate(
        zip(tiles, tile_vertices, strict=True)
    ):
        keys = [fraction_key(time, x) for _, time, x in vertices]
        points = [complex(*point_rows[key]["five_term_midpoint"]) for key in keys]
        winding, raw_winding, maximum_increment = polygon_winding(points)
        weakest_key = min(
            keys,
            key=lambda key: point_rows[key]["complete_vector_norm_lower"],
        )
        weakest = point_rows[weakest_key]
        maximum_tail_ratio = max(
            point_rows[key]["tail_vector_to_five_term_norm_upper"] for key in keys
        )
        tile_rows.append(
            {
                "id": f"tile_{tile_index}",
                "x": [fraction_text(x_low), fraction_text(x_high)],
                "time": [fraction_text(TIME_LOW), fraction_text(TIME_HIGH)],
                "coordinate_orientation": "(t,x)",
                "vertex_keys": keys,
                "summary": {
                    "boundary_vertices": len(keys),
                    "winding_t_x": winding,
                    "raw_winding": raw_winding,
                    "maximum_angle_increment": maximum_increment,
                    "angle_resolution_below_pi_over_two": maximum_increment < math.pi / 2,
                    "all_sampled_complete_vertices_exclude_origin": all(
                        point_rows[key]["complete_vertex_excludes_origin"] for key in keys
                    ),
                    "minimum_complete_norm_lower": weakest["complete_vector_norm_lower"],
                    "minimum_complete_norm_point": {
                        "key": weakest_key,
                        "time": weakest["time"],
                        "x": weakest["x"],
                        "five_term_midpoint": weakest["five_term_midpoint"],
                    },
                    "maximum_complete_tail_to_five_term_vector_ratio": maximum_tail_ratio,
                },
            }
        )

    shared_edges: list[dict] = []
    for left_index in range(len(tiles) - 1):
        left = tiles[left_index]
        right = tiles[left_index + 1]
        if left[1] != right[0]:
            continue
        shared_x = left[1]
        shared_keys = [
            fraction_key(time, shared_x)
            for time in fraction_grid(TIME_LOW, TIME_HIGH, args.segments_per_edge)
        ]
        require(all(key in point_rows for key in shared_keys), "shared edge cache is incomplete")
        shared_edges.append(
            {
                "left_tile": f"tile_{left_index}",
                "right_tile": f"tile_{left_index + 1}",
                "x": fraction_text(shared_x),
                "cached_points": len(shared_keys),
            }
        )

    all_promising = all(
        tile["summary"]["winding_t_x"] == 0
        and tile["summary"]["angle_resolution_below_pi_over_two"]
        and tile["summary"]["all_sampled_complete_vertices_exclude_origin"]
        for tile in tile_rows
    )
    return {
        "kind": STEM,
        "schema_version": 1,
        "date": DATE,
        "status": "finite complete-theta positive-time spatial tile scout complete",
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "parents": {
            "fifth_summand_degree": {
                "path": str(FIFTH_PARENT.relative_to(REPO_ROOT)).replace("\\", "/"),
                "sha256": sha256_path(FIFTH_PARENT),
            },
            "complete_tail": {
                "path": str(TAIL_PARENT.relative_to(REPO_ROOT)).replace("\\", "/"),
                "sha256": sha256_path(TAIL_PARENT),
            },
        },
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
        "configuration": {
            "precision_bits": PRECISION_BITS,
            "absolute_integration_tolerance": continuation.ABS_TOL,
            "segments_per_edge": args.segments_per_edge,
            "time": [fraction_text(TIME_LOW), fraction_text(TIME_HIGH)],
            "tail_derivative_bounds": [
                value
                for value in tail_parent["positive_moment_budget"][
                    "tail_derivative_bounds"
                ]
            ],
            "tail_bound_source_safety_factor": tail_parent["positive_moment_budget"][
                "safety_factor"
            ],
            "finite_parent_minimum_separation": finite_parent["boundary_certificate"][
                "minimum_separation_lower"
            ],
        },
        "cache": {
            "requested_vertex_references": requested_references,
            "unique_point_evaluations": len(point_rows),
            "cache_hits": requested_references - len(point_rows),
            "shared_vertical_edges": shared_edges,
        },
        "tiles": tile_rows,
        "point_cache": point_rows,
        "summary": {
            "tiles": len(tile_rows),
            "distinct_windings": sorted(
                {tile["summary"]["winding_t_x"] for tile in tile_rows}
            ),
            "all_tiles_promising_for_interval_followup": all_promising,
            "minimum_complete_norm_lower": min(
                tile["summary"]["minimum_complete_norm_lower"] for tile in tile_rows
            ),
            "maximum_angle_increment": max(
                tile["summary"]["maximum_angle_increment"] for tile in tile_rows
            ),
            "maximum_complete_tail_to_five_term_vector_ratio": max(
                tile["summary"]["maximum_complete_tail_to_five_term_vector_ratio"]
                for tile in tile_rows
            ),
        },
        "route_decision": (
            "If every tile has sampled winding zero, origin-excluding complete-theta "
            "vertices, and angle increments below pi/2, promote the two adjacent tiles "
            "to an adaptive ACB/Taylor boundary cover. Reuse each shared vertical edge "
            "once, then absorb n>=6 with the uniform derivative budget. Any failed tile "
            "must be isolated before further spatial expansion."
        ),
        "proof_boundary": (
            "Finite 192-bit point diagnostic on the sampled boundaries of the declared "
            "tiles. The n>=6 derivative budget is rigorous at each sampled point, but "
            "there is no interval coverage between vertices, no certified segment "
            "homotopy, and no Brouwer-degree theorem for the new tiles. It does not "
            "control other frequencies, negative heat time, all real x, the global "
            "de Bruijn-Newman constant, Lambda<=0, RH, or a prize-level conclusion."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tile",
        action="append",
        help="Rational or decimal X_LOW:X_HIGH tile; repeat for multiple tiles.",
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--partial", type=Path, default=DEFAULT_PARTIAL)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--segments-per-edge", type=int, default=32)
    parser.add_argument("--checkpoint-every", type=int, default=16)
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
        "complete-theta tile scout baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, "
        f"mean={baseline_mean:.1f}%",
        flush=True,
    )
    if baseline_mean > args.baseline_max_percent:
        print("complete-theta tile scout deferred: daytime baseline is already busy")
        return 2
    require(priority_lowered, "could not apply below-normal process priority")
    continuation.flint.ctx.prec = PRECISION_BITS
    monitor = CpuMonitor()
    try:
        payload = build(args, baseline, monitor)
    except ResourcePark as exc:
        print(f"complete-theta tile scout parked: {exc}")
        return 3
    continuation.write_json_atomic(args.out, payload)
    continuation.write_text_atomic(args.note, render_note(payload))
    print(
        "built complete-theta spatial tile scout: "
        f"tiles={payload['summary']['tiles']}, "
        f"unique_points={payload['cache']['unique_point_evaluations']}, "
        f"cache_hits={payload['cache']['cache_hits']}, "
        f"windings={payload['summary']['distinct_windings']}, "
        f"promising={payload['summary']['all_tiles_promising_for_interval_followup']}, "
        f"elapsed={payload['resource_policy']['elapsed_seconds']:.3f}s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
