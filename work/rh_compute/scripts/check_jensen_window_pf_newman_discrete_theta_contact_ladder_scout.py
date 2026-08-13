#!/usr/bin/env python3
"""Independently check the discrete-theta m=2 contact-ladder scout."""

from __future__ import annotations

import argparse
import functools
import hashlib
import json
import math
import os
from pathlib import Path
import sys

import mpmath as mp
import psutil


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_discrete_theta_contact_ladder_scout"
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
BUILDER = REPO_ROOT / "work/rh_compute/scripts" / f"{STEM}.py"


class ResourcePark(RuntimeError):
    """Raised between complete roster rows after sustained high CPU."""


class CpuMonitor:
    def __init__(self) -> None:
        self.samples: list[float] = []
        self.consecutive_high = 0

    def sample(self) -> None:
        value = float(psutil.cpu_percent(interval=0.25))
        self.samples.append(value)
        if value > 75.0:
            self.consecutive_high += 1
        else:
            self.consecutive_high = 0
        if self.consecutive_high >= 2:
            raise ResourcePark("two consecutive checker CPU samples exceeded 75%")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def integration_breaks(index: mp.mpf) -> list[mp.mpf]:
    scale = 1 / (4 * mp.pi * index * index)
    candidates = [
        mp.mpf("0"),
        scale,
        4 * scale,
        16 * scale,
        mp.mpf("0.02"),
        mp.mpf("0.1"),
        mp.mpf("0.5"),
        mp.mpf("1"),
        mp.mpf("2"),
    ]
    result: list[mp.mpf] = []
    for value in sorted(set(candidates)):
        if value < 0 or value > 2:
            continue
        if not result or value > result[-1]:
            result.append(value)
    require(result[0] == 0 and result[-1] == 2, "invalid integration chart")
    return result


@functools.lru_cache(maxsize=None)
def normalized_jets(
    index_value: float, time_value: float, x_value: float
) -> tuple[mp.mpf, ...]:
    index = mp.mpf(str(index_value))
    time = mp.mpf(str(time_value))
    x = mp.mpf(str(x_value))
    pi_n2 = mp.pi * index * index

    def kernel(u: mp.mpf) -> mp.mpf:
        exp4u = mp.exp(4 * u)
        return (
            pi_n2
            * mp.exp(5 * u + time * u * u)
            * (2 * pi_n2 * exp4u - 3)
            * mp.exp(-pi_n2 * exp4u)
        )

    breaks = integration_breaks(index)
    mass = mp.quad(kernel, breaks)
    require(mass > 0, f"nonpositive high-precision mass at index={index}")
    moments = [
        mp.quad(
            lambda u, order=order: u**order * kernel(u) * mp.e ** (1j * x * u),
            breaks,
        )
        / mass
        for order in range(4)
    ]
    return tuple(mp.re((1j**order) * moments[order]) for order in range(4))


def request_below_normal_priority() -> bool:
    try:
        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
        else:
            process.nice(5)
        return True
    except (psutil.Error, PermissionError, OSError):
        return False


def solve_field_witness(
    indices: list[float], time: float, x: float, target: float
) -> tuple[list[mp.mpf], list[mp.mpf], mp.mpf]:
    jets = [normalized_jets(index, time, x) for index in indices]
    x_mp = mp.mpf(str(x))
    target_mp = mp.mpf(str(target))
    matrix = mp.matrix(
        [
            [mp.mpf(1) for _ in indices],
            [row[0] for row in jets],
            [row[1] for row in jets],
            [x_mp * row[3] - 6 * (target_mp + mp.mpf("0.5")) * row[2] for row in jets],
        ]
    )
    weights = mp.lu_solve(matrix, mp.matrix([1, 0, 0, 0]))
    mixed = [sum(weights[k] * jets[k][j] for k in range(4)) for j in range(4)]
    require(abs(mixed[2]) > mp.mpf("1e-40"), "high-precision witness has zero second jet")
    computed = x_mp * mixed[3] / (6 * mixed[2]) - mp.mpf("0.5")
    return [weights[k] for k in range(4)], mixed, computed


def solve_contact_vertex(
    indices: list[float], time: float, x: float
) -> tuple[list[mp.mpf], list[mp.mpf]]:
    jets = [normalized_jets(index, time, x) for index in indices]
    matrix = mp.matrix(
        [
            [mp.mpf(1) for _ in indices],
            [row[0] for row in jets],
            [row[1] for row in jets],
        ]
    )
    weights = mp.lu_solve(matrix, mp.matrix([1, 0, 0]))
    mixed = [sum(weights[k] * jets[k][j] for k in range(3)) for j in range(4)]
    return [weights[k] for k in range(3)], mixed


def check_target_witnesses(
    payload: dict, monitor: CpuMonitor
) -> tuple[int, mp.mpf, mp.mpf]:
    checked = 0
    minimum_weight = mp.inf
    maximum_field_error = mp.mpf(0)
    for roster in payload["roster_results"]:
        for key, witness in roster["target_witnesses"].items():
            target = float(key)
            weights, mixed, computed = solve_field_witness(
                witness["support_indices"],
                witness["time"],
                witness["x"],
                target,
            )
            local_minimum = min(weights)
            positive_support = sum(
                weight > mp.mpf("1e-20") for weight in weights
            )
            field_error = abs(computed - mp.mpf(str(target)))
            require(local_minimum > mp.mpf("1e-20"), f"nonpositive recomputed witness {roster['name']} ell={key}")
            require(positive_support >= 3, f"degenerate recomputed witness {roster['name']} ell={key}")
            require(abs(mixed[0]) < mp.mpf("1e-55"), f"value contact residual {roster['name']} ell={key}")
            require(abs(mixed[1]) < mp.mpf("1e-55"), f"derivative contact residual {roster['name']} ell={key}")
            require(field_error < mp.mpf("1e-50"), f"field residual {roster['name']} ell={key}")
            minimum_weight = min(minimum_weight, local_minimum)
            maximum_field_error = max(maximum_field_error, field_error)
            checked += 1
        monitor.sample()
    return checked, minimum_weight, maximum_field_error


def check_unbounded_guards(
    payload: dict, monitor: CpuMonitor
) -> tuple[int, mp.mpf]:
    checked = 0
    minimum_third_jet = mp.inf
    for roster in payload["roster_results"]:
        guard = roster["unbounded_adjacent_ratio_guard"]
        if not guard or not guard["detected"]:
            continue
        witness = guard["witness"]
        left_weights, left_mixed = solve_contact_vertex(
            witness["left"]["support_indices"], guard["time"], guard["x"]
        )
        right_weights, right_mixed = solve_contact_vertex(
            witness["right"]["support_indices"], guard["time"], guard["x"]
        )
        require(min(left_weights) >= -mp.mpf("1e-40"), f"negative left vertex {roster['name']}")
        require(min(right_weights) >= -mp.mpf("1e-40"), f"negative right vertex {roster['name']}")
        require(left_mixed[2] * right_mixed[2] < 0, f"second jets do not straddle zero {roster['name']}")
        alpha = -left_mixed[2] / (right_mixed[2] - left_mixed[2])
        third = (1 - alpha) * left_mixed[3] + alpha * right_mixed[3]
        require(0 < alpha < 1, f"invalid segment parameter {roster['name']}")
        require(abs(third) > mp.mpf("1e-14"), f"third jet also vanishes {roster['name']}")
        minimum_third_jet = min(minimum_third_jet, abs(third))
        checked += 1
        monitor.sample()
    return checked, minimum_third_jet


def check_targeted_refinements(payload: dict, monitor: CpuMonitor) -> int:
    refinements = payload.get("targeted_integer_refinement")
    require(isinstance(refinements, list) and len(refinements) == 2, "missing targeted refinements")
    checked = 0
    for refinement in refinements:
        require(
            refinement["maximum_coarse_fine_jet_difference"] < 1e-10,
            "targeted refinement quadrature drift",
        )
        require(refinement["nearest_definite_interval"] is not None, "missing nearest interval")
        witness = refinement["best_target_witness"]
        if witness is not None:
            target = refinement["configuration"]["target"]
            weights, mixed, computed = solve_field_witness(
                witness["support_indices"], witness["time"], witness["x"], target
            )
            require(min(weights) > mp.mpf("1e-20"), "nonpositive targeted witness")
            require(abs(mixed[0]) < mp.mpf("1e-55"), "targeted value residual")
            require(abs(mixed[1]) < mp.mpf("1e-55"), "targeted derivative residual")
            require(
                abs(computed - mp.mpf(str(target))) < mp.mpf("1e-50"),
                "targeted field residual",
            )
            checked += 1
        monitor.sample()
    return checked


def check_structure(payload: dict) -> None:
    require(payload.get("kind") == STEM, "wrong result kind")
    require(payload.get("schema_version") == 1, "wrong schema version")
    require(payload.get("source_sha256") == sha256_path(BUILDER), "builder hash drift")
    resource = payload["resource_policy"]
    require(resource["active_compute_workers"] == 1, "worker policy drift")
    require(resource["baseline_mean_percent"] <= 60.0, "busy baseline should have deferred")
    require(not resource.get("resource_parked", False), "result is a parked partial run")
    require(payload["fixed_theta_diagnostic"] is not None, "missing fixed-weight rung")
    require(payload.get("targeted_integer_refinement") is not None, "missing refinement rung")
    for roster in payload["roster_results"]:
        require(roster["target_fields_realized"] == len(roster["target_witnesses"]), "target count drift")
        require(roster["maximum_coarse_fine_jet_difference"] < 1e-10, "quadrature drift too large")
    fixed = payload["fixed_theta_diagnostic"]
    require(fixed["maximum_coarse_fine_difference"] < 1e-10, "fixed-weight quadrature drift")
    for value in fixed["component_index_tail_moment_bounds"].values():
        require(math.isfinite(value) and value >= 0.0, "invalid index-tail bound")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--baseline-seconds", type=int, default=5)
    parser.add_argument("--baseline-max-percent", type=float, default=60.0)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    priority_lowered = request_below_normal_priority()
    baseline = [
        float(psutil.cpu_percent(interval=1.0))
        for _ in range(args.baseline_seconds)
    ]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "theta-contact checker baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, "
        f"mean={baseline_mean:.1f}%"
    )
    if baseline_mean > args.baseline_max_percent:
        print("theta-contact checker deferred: daytime baseline is already busy")
        return 2
    mp.mp.dps = args.dps
    with args.result.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    check_structure(payload)
    monitor = CpuMonitor()
    try:
        witness_count, minimum_weight, maximum_field_error = check_target_witnesses(
            payload, monitor
        )
        unbounded_count, minimum_third_jet = check_unbounded_guards(
            payload, monitor
        )
        targeted_count = check_targeted_refinements(payload, monitor)
    except ResourcePark as exc:
        print(
            "theta-contact checker parked after a complete roster row: "
            f"{exc}; samples={monitor.samples}"
        )
        return 2
    print(
        "validated discrete theta contact ladder scout: "
        f"{len(payload['roster_results'])} rosters, "
        f"{witness_count} high-precision field witnesses, "
        f"{unbounded_count} unbounded-ratio guards, "
        f"{targeted_count} refined target witnesses, "
        f"minimum positive weight={mp.nstr(minimum_weight, 8)}, "
        f"maximum field error={mp.nstr(maximum_field_error, 5)}, "
        f"minimum crossing third jet={mp.nstr(minimum_third_jet, 8)}, "
        f"runtime CPU samples={monitor.samples}, "
        f"below-normal={priority_lowered}, "
        "0 issues"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
