#!/usr/bin/env python3
"""Independently audit the fixed-theta mass-profile primal-dual scout."""

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
STEM = "jensen_window_pf_newman_fixed_theta_mass_profile_dual_scout"
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
BUILDER = REPO_ROOT / "work/rh_compute/scripts" / f"{STEM}.py"
PARENT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_discrete_theta_contact_ladder_scout.json"
)


class ResourcePark(RuntimeError):
    """Raised between complete high-precision projection rows."""


class CpuMonitor:
    def __init__(self) -> None:
        self.samples: list[float] = []
        self.consecutive_high = 0

    def sample(self) -> None:
        value = float(psutil.cpu_percent(interval=0.5))
        self.samples.append(value)
        self.consecutive_high = self.consecutive_high + 1 if value > 75.0 else 0
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


def integration_breaks(index: mp.mpf) -> list[mp.mpf]:
    scale = 1 / (4 * mp.pi * index * index)
    candidates = [
        mp.mpf("0"),
        scale,
        4 * scale,
        16 * scale,
        mp.mpf("0.02"),
        mp.mpf("0.05"),
        mp.mpf("0.1"),
        mp.mpf("0.2"),
        mp.mpf("0.5"),
        mp.mpf("1"),
        mp.mpf("2"),
    ]
    result: list[mp.mpf] = []
    for value in sorted(set(candidates)):
        if 0 <= value <= 2 and (not result or value > result[-1]):
            result.append(value)
    require(result[0] == 0 and result[-1] == 2, "invalid integration chart")
    return result


@functools.lru_cache(maxsize=None)
def mass_and_jets(
    index_value: int, time_value: float, x_value: float
) -> tuple[mp.mpf, tuple[mp.mpf, ...]]:
    index = mp.mpf(index_value)
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
    require(mass > 0, f"nonpositive high-precision mass at n={index_value}")
    derivatives = (
        lambda u: mp.cos(x * u),
        lambda u: -u * mp.sin(x * u),
        lambda u: -(u**2) * mp.cos(x * u),
        lambda u: (u**3) * mp.sin(x * u),
    )
    jets = tuple(
        mp.quad(lambda u, derivative=derivative: kernel(u) * derivative(u), breaks)
        / mass
        for derivative in derivatives
    )
    return mass, jets


def matrix_from_jets(jets: list[tuple[mp.mpf, ...]]) -> mp.matrix:
    return mp.matrix(
        [
            [mp.mpf(1) for _ in jets],
            [row[0] for row in jets],
            [row[1] for row in jets],
        ]
    )


def profile_at(
    depth: int, time: float, x: float
) -> tuple[list[mp.mpf], list[tuple[mp.mpf, ...]]]:
    rows = [mass_and_jets(index, time, x) for index in range(1, depth + 1)]
    masses = [row[0] for row in rows]
    total = sum(masses)
    return [mass / total for mass in masses], [row[1] for row in rows]


def relative_error(left: mp.mpf, right: mp.mpf) -> mp.mpf:
    return abs(left - right) / max(abs(left), abs(right), mp.mpf("1e-100"))


def solve_face_projection(
    profile: list[mp.mpf],
    jets: list[tuple[mp.mpf, ...]],
    active: list[int],
    zeroed: list[int],
) -> dict:
    matrix = matrix_from_jets(jets)
    projected = profile.copy()
    for position in active + zeroed:
        projected[position] = mp.mpf(0)
    target = mp.matrix([1, 0, 0]) - matrix * mp.matrix(projected)
    active_matrix = mp.matrix(
        [[matrix[row, position] for position in active] for row in range(3)]
    )
    solution = mp.lu_solve(active_matrix, target)
    for local, position in enumerate(active):
        projected[position] = solution[local]
    perturbation = [projected[index] - profile[index] for index in range(len(profile))]
    residual = matrix * mp.matrix(projected) - mp.matrix([1, 0, 0])
    require(max(abs(value) for value in residual) < mp.mpf("1e-60"), "contact solve residual")
    require(
        min(projected) >= -mp.mpf("1e-60"),
        "negative high-precision projection: "
        + ",".join(mp.nstr(value, 12) for value in projected),
    )

    signs = mp.matrix(
        [
            -mp.mpf("0.5") if perturbation[position] > 0 else mp.mpf("0.5")
            for position in active
        ]
    )
    dual = mp.lu_solve(active_matrix.T, signs)
    observable = [sum(matrix[row, index] * dual[row] for row in range(3)) for index in range(len(profile))]
    minimum_slack = mp.inf
    for index, value in enumerate(observable):
        if index in active:
            expected = -mp.mpf("0.5") if perturbation[index] > 0 else mp.mpf("0.5")
            require(abs(value - expected) < mp.mpf("1e-60"), "active KKT equality")
            continue
        if index in zeroed:
            slack = value - mp.mpf("0.5")
        elif profile[index] == 0:
            require(projected[index] == 0, "inactive zero-profile component moved")
            slack = value + mp.mpf("0.5")
        else:
            require(perturbation[index] == 0, "unexpected inactive perturbation")
            slack = mp.mpf("0.5") - abs(value)
        require(slack > -mp.mpf("1e-45"), f"KKT dual inequality at component {index + 1}")
        minimum_slack = min(minimum_slack, slack)
    television = mp.mpf("0.5") * sum(abs(value) for value in perturbation)
    return {
        "projected": projected,
        "perturbation": perturbation,
        "television": television,
        "dual": dual,
        "dual_observable": observable,
        "minimum_dual_slack": minimum_slack,
    }


def fisher_projection(
    profile: list[mp.mpf], jets: list[tuple[mp.mpf, ...]]
) -> dict:
    mean = [
        sum(profile[index] * jets[index][order] for index in range(len(profile)))
        for order in range(2)
    ]
    centered = [
        (jets[index][0] - mean[0], jets[index][1] - mean[1])
        for index in range(len(profile))
    ]
    covariance = mp.matrix(
        [
            [
                sum(
                    profile[index] * centered[index][left] * centered[index][right]
                    for index in range(len(profile))
                )
                for right in range(2)
            ]
            for left in range(2)
        ]
    )
    alpha = mp.lu_solve(covariance, mp.matrix(mean))
    observable = [
        centered[index][0] * alpha[0] + centered[index][1] * alpha[1]
        for index in range(len(profile))
    ]
    projected = [
        profile[index] * (1 - observable[index])
        for index in range(len(profile))
    ]
    distance_squared = mean[0] * alpha[0] + mean[1] * alpha[1]
    require(min(projected) > 0, "high-precision Fisher projection is not positive")
    return {
        "alpha": alpha,
        "observable": observable,
        "projected": projected,
        "distance_squared": distance_squared,
    }


def check_structure(payload: dict) -> None:
    require(payload.get("kind") == STEM, "wrong result kind")
    require(payload.get("schema_version") == 1, "wrong schema version")
    require(payload.get("source_sha256") == sha256_path(BUILDER), "builder hash drift")
    require(
        payload["parent"]["result_sha256"] == sha256_path(PARENT_RESULT),
        "parent result hash drift",
    )
    resource = payload["resource_policy"]
    require(resource["active_compute_workers"] == 1, "worker policy drift")
    require(resource["thread_caps"] == 1, "native thread cap drift")
    require(resource["baseline_mean_percent"] <= 60.0, "busy baseline should defer")
    require(not resource["resource_parked"], "cannot certify a parked partial result")
    require(payload["grid"]["depths"] == [4, 8, 12], "depth ladder drift")
    require(len(payload["depth_results"]) == 3, "depth result count drift")
    require(payload["maximum_coarse_fine_difference"] < 1e-10, "quadrature drift")
    require(
        not payload["numerical_resolution_gate"]["distance_values_resolved"],
        "selected locator objective should be marked unresolved",
    )
    objectives = [
        row["nearest_total_variation_contact"]["solver_objective"]
        for row in payload["depth_results"]
    ]
    require(objectives[1] <= objectives[0], "depth-8 contact set did not enlarge")
    require(objectives[2] <= objectives[1], "depth-12 contact set did not enlarge")
    closure = payload["infinite_endpoint_closure"]["nearest_total_variation_contact"]
    require(closure["solver_objective"] <= objectives[2], "endpoint closure did not enlarge")
    for row in payload["depth_results"]:
        require(row["tv_to_tail_ratio"] > 1e12, "distance not separated from tail budget")


def check_projections(payload: dict, monitor: CpuMonitor) -> dict:
    cases = (
        (4, [0, 1, 2], [3]),
        (8, [0, 1, 2], list(range(3, 8))),
        (12, [0, 1, 2], list(range(3, 12))),
    )
    maximum_jet_error = mp.mpf(0)
    maximum_profile_error = mp.mpf(0)
    maximum_tv_relative_error = mp.mpf(0)
    maximum_fisher_relative_error = mp.mpf(0)
    minimum_dual_slack = mp.inf
    fisher_rows: dict[int, dict] = {}

    for depth, active, zeroed in cases:
        stored = payload["depth_results"][[4, 8, 12].index(depth)][
            "nearest_total_variation_contact"
        ]
        time = float(stored["time"])
        x = float(stored["x"])
        require(time == 0.2 and x == 133.675, "selected finite minimum drift")
        profile, jets = profile_at(depth, time, x)
        stored_profile = [mp.mpf(str(value)) for value in stored["retained_profile_weights"]]
        stored_jets = stored["component_jets"]
        maximum_profile_error = max(
            maximum_profile_error,
            max(relative_error(profile[index], stored_profile[index]) for index in range(depth)),
        )
        maximum_jet_error = max(
            maximum_jet_error,
            max(
                abs(jets[index][order] - mp.mpf(str(stored_jets[index][order])))
                for index in range(depth)
                for order in range(4)
            ),
        )
        projection = solve_face_projection(profile, jets, active, zeroed)
        stored_tv = mp.mpf(str(stored["solver_objective"]))
        maximum_tv_relative_error = max(
            maximum_tv_relative_error,
            relative_error(projection["television"], stored_tv),
        )
        minimum_dual_slack = min(minimum_dual_slack, projection["minimum_dual_slack"])

        fisher = fisher_projection(profile, jets)
        fisher_rows[depth] = fisher
        stored_fisher = stored["fisher_projection_at_same_point"]
        stored_distance = mp.mpf(str(stored_fisher["chi_square_distance_squared"]))
        maximum_fisher_relative_error = max(
            maximum_fisher_relative_error,
            relative_error(fisher["distance_squared"], stored_distance),
            relative_error(fisher["alpha"][0], mp.mpf(str(stored_fisher["feature_dual_original"][0]))),
            relative_error(fisher["alpha"][1], mp.mpf(str(stored_fisher["feature_dual_original"][1]))),
        )
        monitor.sample()

    closure = payload["infinite_endpoint_closure"]["nearest_total_variation_contact"]
    profile, jets = profile_at(12, float(closure["time"]), float(closure["x"]))
    profile.append(mp.mpf(0))
    jets.append((mp.mpf(1), mp.mpf(0), mp.mpf(0), mp.mpf(0)))
    closure_projection = solve_face_projection(
        profile, jets, [0, 1, 2], list(range(3, 12))
    )
    maximum_tv_relative_error = max(
        maximum_tv_relative_error,
        relative_error(
            closure_projection["television"], mp.mpf(str(closure["solver_objective"]))
        ),
    )
    minimum_dual_slack = min(
        minimum_dual_slack, closure_projection["minimum_dual_slack"]
    )
    monitor.sample()

    observable = fisher_rows[12]["observable"]
    require(all(value > 0 for value in observable[2:]), "Fisher tail sign failed")
    require(
        all(observable[index + 1] > observable[index] for index in range(2, len(observable) - 1)),
        "Fisher observable tail is not strictly increasing in n",
    )
    require(
        maximum_jet_error < mp.mpf("1e-8"),
        f"high-precision jet mismatch {mp.nstr(maximum_jet_error, 8)}",
    )
    require(
        maximum_profile_error < mp.mpf("1e-6"),
        f"high-precision profile mismatch {mp.nstr(maximum_profile_error, 8)}",
    )
    stored_closure = mp.mpf(str(closure["solver_objective"]))
    require(
        closure_projection["television"] < stored_closure / 50,
        "resolution gate did not expose the expected locator discrepancy",
    )
    require(
        maximum_fisher_relative_error > mp.mpf("0.5"),
        "resolution stress no longer exposes the double-precision Fisher dual",
    )
    return {
        "maximum_jet_error": maximum_jet_error,
        "maximum_profile_error": maximum_profile_error,
        "maximum_tv_relative_error": maximum_tv_relative_error,
        "maximum_fisher_relative_error": maximum_fisher_relative_error,
        "minimum_dual_slack": minimum_dual_slack,
        "closure_tv": closure_projection["television"],
        "fisher_alpha": fisher_rows[12]["alpha"],
        "fisher_observable": observable,
    }


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
    baseline = [float(psutil.cpu_percent(interval=1.0)) for _ in range(args.baseline_seconds)]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "fixed-theta dual checker baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, mean={baseline_mean:.1f}%"
    )
    if baseline_mean > args.baseline_max_percent:
        print("fixed-theta dual checker deferred: daytime baseline is already busy")
        return 2
    mp.mp.dps = args.dps
    with args.result.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    check_structure(payload)
    monitor = CpuMonitor()
    try:
        audit = check_projections(payload, monitor)
    except ResourcePark as exc:
        print(f"fixed-theta dual checker parked: {exc}; samples={monitor.samples}")
        return 2
    print(
        "validated fixed-theta mass-profile dual scout: "
        f"max jet error={mp.nstr(audit['maximum_jet_error'], 8)}, "
        f"max profile relative error={mp.nstr(audit['maximum_profile_error'], 8)}, "
        f"max TV relative error={mp.nstr(audit['maximum_tv_relative_error'], 8)}, "
        f"max Fisher relative error={mp.nstr(audit['maximum_fisher_relative_error'], 8)}, "
        f"minimum KKT slack={mp.nstr(audit['minimum_dual_slack'], 8)}, "
        f"closure TV={mp.nstr(audit['closure_tv'], 12)}, "
        f"Fisher alpha=({mp.nstr(audit['fisher_alpha'][0], 12)},"
        f"{mp.nstr(audit['fisher_alpha'][1], 12)}), "
        f"runtime CPU samples={monitor.samples}, below-normal={priority_lowered}, 0 issues"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
