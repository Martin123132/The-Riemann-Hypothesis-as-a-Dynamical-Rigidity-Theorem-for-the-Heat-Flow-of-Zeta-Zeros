#!/usr/bin/env python3
"""Independently check the fourth theta-summand tangency homotopy."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys

import mpmath as mp
import psutil


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_theta_fourth_summand_tangency_homotopy_scout"
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
BUILDER = REPO_ROOT / "work/rh_compute/scripts" / f"{STEM}.py"
PARENT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_fixed_theta_mass_profile_dual_scout.json"
)


class ResourcePark(RuntimeError):
    """Raised between complete independently integrated anchor rows."""


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


def relative_error(left: mp.mpf, right: mp.mpf) -> mp.mpf:
    return abs(left - right) / max(abs(left), abs(right), mp.mpf("1e-100"))


def independent_breaks(index: int, x: mp.mpf) -> list[mp.mpf]:
    scale = 1 / (4 * mp.pi * index * index)
    values = {
        mp.mpf(0),
        scale,
        4 * scale,
        16 * scale,
        mp.mpf("0.02"),
        mp.mpf("0.05"),
        mp.mpf("0.1"),
        mp.mpf("0.2"),
        mp.mpf("0.5"),
        mp.mpf("0.75"),
        mp.mpf(1),
        mp.mpf(2),
    }
    phase_step = mp.pi / x
    phase_index = 1
    while phase_index * phase_step < mp.mpf("0.75"):
        values.add(phase_index * phase_step)
        phase_index += 1
    return sorted(value for value in values if 0 <= value <= 2)


def component_jets(
    index: int, time: mp.mpf, x: mp.mpf
) -> list[mp.mpf]:
    n = mp.mpf(index)
    pi_n2 = mp.pi * n * n

    def kernel(u: mp.mpf) -> mp.mpf:
        exp4u = mp.exp(4 * u)
        return (
            pi_n2
            * mp.exp(5 * u + time * u * u)
            * (2 * pi_n2 * exp4u - 3)
            * mp.exp(-pi_n2 * exp4u)
        )

    derivatives = (
        lambda u: mp.cos(x * u),
        lambda u: -u * mp.sin(x * u),
        lambda u: -(u**2) * mp.cos(x * u),
        lambda u: (u**3) * mp.sin(x * u),
    )
    breaks = independent_breaks(index, x)
    return [
        mp.quad(lambda u, derivative=derivative: kernel(u) * derivative(u), breaks)
        for derivative in derivatives
    ]


def field_at(
    coefficient: mp.mpf, time: mp.mpf, x: mp.mpf
) -> tuple[list[mp.mpf], list[mp.mpf]]:
    rows = [component_jets(index, time, x) for index in range(1, 5)]
    field = [
        sum(rows[index][order] for index in range(3))
        + coefficient * rows[3][order]
        for order in range(4)
    ]
    return field, rows[3]


def tangent(field: list[mp.mpf], fourth: list[mp.mpf]) -> mp.matrix:
    jacobian = mp.matrix(
        [[-field[2], field[1]], [-field[3], field[2]]]
    )
    return mp.lu_solve(jacobian, mp.matrix([-fourth[0], -fourth[1]]))


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
    require(resource["thread_caps"] == 1, "thread cap drift")
    require(resource["baseline_mean_percent"] <= 60.0, "busy baseline should defer")
    require(not resource["resource_parked"], "cannot check a parked result")
    require(payload["precision"]["decimal_digits"] >= 50, "precision drift")
    require(not payload["precision"]["interval_certified"], "proof status drift")
    require("Jacobi-theta" in payload["exact_setup"]["pi_provenance"], "pi provenance missing")

    rows = payload["continuation"]["rows"]
    require(len(rows) == 12, "continuation row count drift")
    expected = [mp.mpf(index) / 50 for index in range(12)]
    actual = [mp.mpf(row["lambda"]) for row in rows]
    require(actual == expected, "lambda grid drift")
    times = [mp.mpf(row["time"]) for row in rows]
    xs = [mp.mpf(row["x"]) for row in rows]
    require(all(times[index + 1] < times[index] for index in range(11)), "time monotonicity")
    require(all(xs[index + 1] > xs[index] for index in range(11)), "x monotonicity")
    require(all(mp.mpf(row["dt_dlambda"]) < 0 for row in rows), "response sign drift")
    require(all(mp.mpf(row["dx_dlambda"]) > 0 for row in rows), "x response sign drift")
    require(
        min(abs(mp.mpf(row["field_jets_0_to_3"][2])) for row in rows)
        > mp.mpf("1e-22"),
        "sampled Jacobian degeneracy",
    )
    require(
        max(
            max(abs(mp.mpf(row["field_jets_0_to_3"][order])) for order in range(2))
            for row in rows
        )
        < mp.mpf("1e-40"),
        "stored contact residual drift",
    )
    crossing = payload["boundary_crossing"]
    crossing_lambda = mp.mpf(crossing["lambda"])
    require(mp.mpf("0.2") < crossing_lambda < mp.mpf("0.22"), "crossing bracket drift")
    require(crossing_lambda < 1, "crossing no longer precedes the true coefficient")
    require(mp.mpf(crossing["dt_dlambda"]) < -2, "crossing is not transverse")
    for component in payload["precision"]["cutoff_tail_log10_bounds"].values():
        for value in component.values():
            require(math.isfinite(value) and value < -1000, "cutoff tail weakened")


def check_anchor(row: dict) -> dict:
    coefficient = mp.mpf(row["lambda"])
    time = mp.mpf(row["time"])
    x = mp.mpf(row["x"])
    field, fourth = field_at(coefficient, time, x)
    response = tangent(field, fourth)
    stored_field = [mp.mpf(value) for value in row["field_jets_0_to_3"]]
    maximum_jet_relative_error = max(
        relative_error(field[order], stored_field[order]) for order in (2, 3)
    )
    maximum_response_relative_error = max(
        relative_error(response[0], mp.mpf(row["dt_dlambda"])),
        relative_error(response[1], mp.mpf(row["dx_dlambda"])),
    )
    residual = max(abs(field[0]), abs(field[1]))
    determinant = -(field[2] ** 2)
    require(residual < mp.mpf("1e-38"), "independent anchor contact residual")
    require(field[2] < -mp.mpf("1e-22"), "independent anchor lost regularity")
    require(response[0] < 0 and response[1] > 0, "independent response sign")
    require(
        relative_error(determinant, mp.mpf(row["jacobian_determinant"]))
        < mp.mpf("1e-20"),
        "Jacobian determinant mismatch",
    )
    return {
        "residual": residual,
        "jet_relative_error": maximum_jet_relative_error,
        "response_relative_error": maximum_response_relative_error,
    }


def check_anchors(payload: dict, monitor: CpuMonitor) -> dict:
    rows = payload["continuation"]["rows"]
    anchors = [rows[0], rows[10], rows[11], payload["boundary_crossing"]]
    audits = []
    for row in anchors:
        audits.append(check_anchor(row))
        monitor.sample()
    maximum_residual = max(audit["residual"] for audit in audits)
    maximum_jet_error = max(audit["jet_relative_error"] for audit in audits)
    maximum_response_error = max(
        audit["response_relative_error"] for audit in audits
    )
    require(maximum_jet_error < mp.mpf("1e-18"), "independent jet mismatch")
    require(maximum_response_error < mp.mpf("1e-18"), "independent response mismatch")
    return {
        "anchor_count": len(anchors),
        "maximum_residual": maximum_residual,
        "maximum_jet_relative_error": maximum_jet_error,
        "maximum_response_relative_error": maximum_response_error,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--dps", type=int, default=75)
    parser.add_argument("--baseline-seconds", type=int, default=5)
    parser.add_argument("--baseline-max-percent", type=float, default=60.0)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    priority_lowered = request_below_normal_priority()
    baseline = [float(psutil.cpu_percent(interval=1.0)) for _ in range(args.baseline_seconds)]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "theta homotopy checker baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, mean={baseline_mean:.1f}%"
    )
    if baseline_mean > args.baseline_max_percent:
        print("theta homotopy checker deferred: daytime baseline is already busy")
        return 2
    mp.mp.dps = args.dps
    with args.result.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    check_structure(payload)
    monitor = CpuMonitor()
    try:
        audit = check_anchors(payload, monitor)
    except ResourcePark as exc:
        print(f"theta homotopy checker parked: {exc}; samples={monitor.samples}")
        return 2
    print(
        "validated fourth-summand tangency homotopy: "
        f"{audit['anchor_count']} independent anchors, "
        f"maximum contact residual={mp.nstr(audit['maximum_residual'], 8)}, "
        f"maximum jet relative error={mp.nstr(audit['maximum_jet_relative_error'], 8)}, "
        f"maximum response relative error={mp.nstr(audit['maximum_response_relative_error'], 8)}, "
        f"runtime CPU samples={monitor.samples}, below-normal={priority_lowered}, 0 issues"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
