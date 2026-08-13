#!/usr/bin/env python3
"""Scout which parts of the discrete theta roster constrain m=2 contacts."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import psutil
from scipy.special import roots_legendre


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_discrete_theta_contact_ladder_scout"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SCHEMA_VERSION = 1
DATE = "2026-08-03"

TARGET_FIELDS = (-20.0, -5.0, -1.0, 0.0, 1.0, 5.0, 20.0)
TIMES = (0.0, 0.05, 0.1, 0.15, 0.2)
X_MIN = 2.0
X_MAX = 160.0
FIXED_WEIGHT_X_MAX = 80.0
X_STEP = 0.5
CUTOFF = 2.0
COARSE_ORDER = 600
FINE_ORDER = 900
INTEGER_TAIL_START = 9


@dataclass(frozen=True)
class Roster:
    name: str
    role: str
    indices: tuple[float, ...]


ROSTERS = (
    Roster(
        "continuum_log_mesh",
        "continuous-shift relaxation",
        tuple(float(value) for value in np.geomspace(1.0, 8.0, 17)),
    ),
    Roster(
        "geometric_sqrt2",
        "geometric q-lattice relaxation",
        tuple(float(2.0 ** (index / 2.0)) for index in range(7)),
    ),
    Roster(
        "integer_1_8",
        "true integer shifts with free positive mixture weights",
        tuple(float(index) for index in range(1, 9)),
    ),
    Roster(
        "integer_1_12",
        "extended true integer shifts with free positive mixture weights",
        tuple(float(index) for index in range(1, 13)),
    ),
)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def write_text_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


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


def sample_baseline(seconds: int) -> list[float]:
    return [float(psutil.cpu_percent(interval=1.0)) for _ in range(seconds)]


def component_kernel(index: float, time: float, u: np.ndarray) -> np.ndarray:
    """Return exp(t u^2) phi_index(u) on the positive half-line."""
    pi_n2 = math.pi * index * index
    exp4u = np.exp(4.0 * u)
    return (
        pi_n2
        * np.exp(5.0 * u + time * u * u)
        * (2.0 * pi_n2 * exp4u - 3.0)
        * np.exp(-pi_n2 * exp4u)
    )


def normalized_component_jets(
    indices: tuple[float, ...],
    time: float,
    xs: np.ndarray,
    order: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return normalized jets j=0..3 and unnormalized component masses."""
    nodes, weights = roots_legendre(order)
    u = (CUTOFF / 2.0) * (nodes + 1.0)
    quadrature_weights = (CUTOFF / 2.0) * weights
    xu = np.outer(xs, u)
    cosine = np.cos(xu)
    sine = np.sin(xu)

    jets = np.empty((len(indices), 4, len(xs)), dtype=np.float64)
    masses = np.empty(len(indices), dtype=np.float64)
    for component, index in enumerate(indices):
        weighted = component_kernel(index, time, u) * quadrature_weights
        mass = float(np.sum(weighted))
        if not math.isfinite(mass) or mass <= 0.0:
            raise RuntimeError(f"nonpositive component mass for index={index}")
        normalized = weighted / mass
        masses[component] = mass
        jets[component, 0] = cosine @ normalized
        jets[component, 1] = -(sine @ (u * normalized))
        jets[component, 2] = -(cosine @ (u * u * normalized))
        jets[component, 3] = sine @ (u * u * u * normalized)
    return jets, masses


def solve_target_witness(
    jets_at_x: np.ndarray,
    x: float,
    target: float,
) -> dict | None:
    """Find a positive four-component division-free contact witness."""
    supports = np.asarray(
        list(itertools.combinations(range(jets_at_x.shape[0]), 4)),
        dtype=np.int64,
    )
    if len(supports) == 0:
        return None
    selected = jets_at_x[supports]
    field_rows = x * selected[:, :, 3] - 6.0 * (target + 0.5) * selected[:, :, 2]
    matrices = np.stack(
        (
            np.ones((len(supports), 4)),
            selected[:, :, 0],
            selected[:, :, 1],
            field_rows,
        ),
        axis=1,
    )
    conditions = np.linalg.cond(matrices)
    valid = np.isfinite(conditions) & (conditions <= 1.0e12)
    if not np.any(valid):
        return None
    valid_positions = np.flatnonzero(valid)
    valid_matrices = matrices[valid]
    right_hand_sides = np.zeros((len(valid_matrices), 4, 1), dtype=np.float64)
    right_hand_sides[:, 0, 0] = 1.0
    try:
        weights = np.linalg.solve(valid_matrices, right_hand_sides)[:, :, 0]
    except np.linalg.LinAlgError:
        return None
    selected = selected[valid]
    conditions = conditions[valid]
    mixed = np.einsum("ki,kij->kj", weights, selected)
    minimum_weights = np.min(weights, axis=1)
    third_largest_weights = np.sort(weights, axis=1)[:, -3]
    admissible = (
        (minimum_weights > 1.0e-12)
        & (third_largest_weights > 1.0e-10)
        & (np.abs(mixed[:, 2]) > 1.0e-12)
    )
    if not np.any(admissible):
        return None
    weights = weights[admissible]
    mixed = mixed[admissible]
    minimum_weights = minimum_weights[admissible]
    conditions = conditions[admissible]
    valid_positions = valid_positions[admissible]
    scores = np.minimum.reduce(
        (
            minimum_weights,
            np.abs(mixed[:, 2]),
            1.0 / np.maximum(conditions, 1.0),
        )
    )
    winner = int(np.argmax(scores))
    winner_matrix = matrices[valid_positions[winner]]
    winner_weights = weights[winner]
    winner_mixed = mixed[winner]
    computed = x * float(winner_mixed[3]) / (6.0 * float(winner_mixed[2])) - 0.5
    equation_residual = float(
        np.max(
            np.abs(
                winner_matrix @ winner_weights
                - np.array([1.0, 0.0, 0.0, 0.0])
            )
        )
    )
    support = supports[valid_positions[winner]]
    return {
        "support_positions": [int(value) for value in support],
        "weights": [float(value) for value in winner_weights],
        "mixed_jets": [float(value) for value in winner_mixed],
        "computed_field": computed,
        "target_field": target,
        "field_error": computed - target,
        "linear_equation_residual": equation_residual,
        "condition_number": float(conditions[winner]),
        "minimum_weight": float(minimum_weights[winner]),
        "third_largest_weight": float(
            third_largest_weights[admissible][winner]
        ),
        "score": float(scores[winner]),
    }


def sparse_contact_vertices(jets_at_x: np.ndarray) -> list[dict]:
    """Enumerate positive three-component vertices of the contact slice."""
    supports = np.asarray(
        list(itertools.combinations(range(jets_at_x.shape[0]), 3)),
        dtype=np.int64,
    )
    if len(supports) == 0:
        return []
    selected = jets_at_x[supports]
    matrices = np.stack(
        (
            np.ones((len(supports), 3)),
            selected[:, :, 0],
            selected[:, :, 1],
        ),
        axis=1,
    )
    conditions = np.linalg.cond(matrices)
    valid = np.isfinite(conditions) & (conditions <= 1.0e12)
    if not np.any(valid):
        return []
    valid_positions = np.flatnonzero(valid)
    matrices = matrices[valid]
    selected = selected[valid]
    right_hand_sides = np.zeros((len(matrices), 3, 1), dtype=np.float64)
    right_hand_sides[:, 0, 0] = 1.0
    try:
        weights = np.linalg.solve(matrices, right_hand_sides)[:, :, 0]
    except np.linalg.LinAlgError:
        return []
    admissible = np.min(weights, axis=1) >= -1.0e-10
    weights = np.maximum(weights[admissible], 0.0)
    selected = selected[admissible]
    valid_positions = valid_positions[admissible]
    if len(weights) == 0:
        return []
    weights /= np.sum(weights, axis=1)[:, None]
    mixed = np.einsum("ki,kij->kj", weights, selected)
    return [
        {
            "support_positions": [int(value) for value in supports[position]],
            "weights": [float(value) for value in weights[row]],
            "h2": float(mixed[row, 2]),
            "h3": float(mixed[row, 3]),
        }
        for row, position in enumerate(valid_positions)
    ]


def unbounded_ratio_guard(vertices: list[dict]) -> dict:
    """Detect a contact segment crossing H''=0 with nonzero H'''."""
    best: dict | None = None
    for left, right in itertools.combinations(vertices, 2):
        if left["h2"] * right["h2"] >= 0.0:
            continue
        alpha = -left["h2"] / (right["h2"] - left["h2"])
        h3 = (1.0 - alpha) * left["h3"] + alpha * right["h3"]
        margin = abs(h3)
        if best is None or margin > best["h3_at_h2_zero_abs"]:
            best = {
                "left": left,
                "right": right,
                "segment_parameter": alpha,
                "h3_at_h2_zero": h3,
                "h3_at_h2_zero_abs": margin,
            }
    return {
        "detected": best is not None and best["h3_at_h2_zero_abs"] > 1.0e-12,
        "witness": best,
    }


def cutoff_tail_log10_bound(
    index: float, derivative_order: int, time_max: float = 0.2
) -> float:
    """Return log10 of a rigorous envelope for the omitted u>=2 integral."""
    epsilon = time_max / (4.0 * math.e * math.e)
    alpha = math.pi * index * index - epsilon
    y0 = math.exp(4.0 * CUTOFF)
    denominator = alpha - 1.5 / y0
    log_bound = (
        2.0 * math.log(math.pi)
        + 4.0 * math.log(index)
        + math.lgamma(derivative_order + 1.0)
        - math.log(2.0)
        + 1.5 * math.log(y0)
        - alpha * y0
        - math.log(denominator)
    )
    return log_bound / math.log(10.0)


def component_index_tail_bound(start: int, derivative_order: int) -> float:
    """Bound sum_{n>=start} integral_0^inf u^j e^(tu^2) phi_n for t<=1/5."""
    epsilon = 0.2 / (4.0 * math.e * math.e)
    total = 0.0
    for index in range(start, start + 200):
        alpha = math.pi * index * index - epsilon
        term = (
            math.pi**2
            * index**4
            * math.factorial(derivative_order)
            / 2.0
            * math.exp(-alpha)
            / (alpha - 1.5)
        )
        total += term
        if term == 0.0:
            break
    return total


def scan_roster(roster: Roster, xs: np.ndarray) -> dict:
    target_witnesses: dict[str, dict] = {}
    unbounded_witness: dict | None = None
    quadrature_max_difference = 0.0
    evaluated_points = 0
    contact_vertex_points = 0
    definite_denominator_points = 0
    sampled_vertex_field_min = math.inf
    sampled_vertex_field_max = -math.inf
    interval_target_distances = {
        f"{target:g}": math.inf for target in TARGET_FIELDS
    }
    interval_target_nearest: dict[str, dict | None] = {
        f"{target:g}": None for target in TARGET_FIELDS
    }

    for time in TIMES:
        coarse, _ = normalized_component_jets(
            roster.indices, time, xs, COARSE_ORDER
        )
        fine, _ = normalized_component_jets(
            roster.indices, time, xs, FINE_ORDER
        )
        quadrature_max_difference = max(
            quadrature_max_difference,
            float(np.max(np.abs(fine - coarse))),
        )
        for x_position, x in enumerate(xs):
            jets_at_x = fine[:, :, x_position]
            evaluated_points += 1
            for target in TARGET_FIELDS:
                key = f"{target:g}"
                if key in target_witnesses:
                    continue
                witness = solve_target_witness(jets_at_x, float(x), target)
                if witness is not None:
                    witness["time"] = time
                    witness["x"] = float(x)
                    witness["support_indices"] = [
                        roster.indices[position]
                        for position in witness["support_positions"]
                    ]
                    target_witnesses[key] = witness
            vertices = sparse_contact_vertices(jets_at_x)
            if vertices:
                contact_vertex_points += 1
                second_jets = np.asarray([vertex["h2"] for vertex in vertices])
                if np.all(second_jets > 1.0e-12) or np.all(second_jets < -1.0e-12):
                    definite_denominator_points += 1
                    fields = np.asarray(
                        [
                            float(x) * vertex["h3"] / (6.0 * vertex["h2"]) - 0.5
                            for vertex in vertices
                        ]
                    )
                    local_minimum = float(np.min(fields))
                    local_maximum = float(np.max(fields))
                    sampled_vertex_field_min = min(
                        sampled_vertex_field_min, local_minimum
                    )
                    sampled_vertex_field_max = max(
                        sampled_vertex_field_max, local_maximum
                    )
                    for target in TARGET_FIELDS:
                        if local_minimum <= target <= local_maximum:
                            distance = 0.0
                        else:
                            distance = min(
                                abs(target - local_minimum),
                                abs(target - local_maximum),
                            )
                        key = f"{target:g}"
                        if distance < interval_target_distances[key]:
                            interval_target_distances[key] = distance
                            interval_target_nearest[key] = {
                                "distance": distance,
                                "time": time,
                                "x": float(x),
                                "local_field_min": local_minimum,
                                "local_field_max": local_maximum,
                            }
            if unbounded_witness is None:
                guard = unbounded_ratio_guard(vertices)
                if guard["detected"]:
                    guard["time"] = time
                    guard["x"] = float(x)
                    for side in ("left", "right"):
                        guard["witness"][side]["support_indices"] = [
                            roster.indices[position]
                            for position in guard["witness"][side]["support_positions"]
                        ]
                    unbounded_witness = guard
            if (
                len(target_witnesses) == len(TARGET_FIELDS)
                and unbounded_witness is not None
            ):
                break
        if (
            len(target_witnesses) == len(TARGET_FIELDS)
            and unbounded_witness is not None
        ):
            break

    return {
        "name": roster.name,
        "role": roster.role,
        "indices": list(roster.indices),
        "evaluated_time_x_points": evaluated_points,
        "target_fields": list(TARGET_FIELDS),
        "target_witnesses": target_witnesses,
        "target_fields_realized": len(target_witnesses),
        "all_target_fields_realized": len(target_witnesses) == len(TARGET_FIELDS),
        "unbounded_adjacent_ratio_guard": unbounded_witness,
        "sampled_contact_polytope_geometry": {
            "contact_vertex_points": contact_vertex_points,
            "definite_denominator_points": definite_denominator_points,
            "sampled_vertex_field_min": (
                sampled_vertex_field_min
                if math.isfinite(sampled_vertex_field_min)
                else None
            ),
            "sampled_vertex_field_max": (
                sampled_vertex_field_max
                if math.isfinite(sampled_vertex_field_max)
                else None
            ),
            "target_distances_to_definite_field_intervals": {
                key: value if math.isfinite(value) else None
                for key, value in interval_target_distances.items()
            },
            "nearest_definite_field_intervals": interval_target_nearest,
        },
        "maximum_coarse_fine_jet_difference": quadrature_max_difference,
    }


def fixed_theta_diagnostic(xs: np.ndarray) -> dict:
    indices = tuple(float(index) for index in range(1, INTEGER_TAIL_START))
    best: dict | None = None
    quadrature_max_difference = 0.0
    sparse_contact_upper_bound: dict | None = None

    for time in TIMES:
        coarse, coarse_masses = normalized_component_jets(
            indices, time, xs, COARSE_ORDER
        )
        fine, fine_masses = normalized_component_jets(
            indices, time, xs, FINE_ORDER
        )
        quadrature_max_difference = max(
            quadrature_max_difference,
            float(np.max(np.abs(fine - coarse))),
            float(np.max(np.abs(fine_masses - coarse_masses))),
        )
        actual_weights = fine_masses / float(np.sum(fine_masses))
        mixed = np.tensordot(actual_weights, fine, axes=(0, 0))
        value_scale = np.tensordot(
            actual_weights, np.abs(fine[:, 0, :]), axes=(0, 0)
        )
        derivative_scale = np.tensordot(
            actual_weights, np.abs(fine[:, 1, :]), axes=(0, 0)
        )
        value_scale = np.maximum(value_scale, 1.0e-300)
        derivative_scale = np.maximum(derivative_scale, 1.0e-300)
        residual = np.maximum(
            np.abs(mixed[0]) / value_scale,
            np.abs(mixed[1]) / derivative_scale,
        )
        position = int(np.argmin(residual))
        candidate = {
            "time": time,
            "x": float(xs[position]),
            "normalized_contact_residual": float(residual[position]),
            "mixed_jets": [float(mixed[j, position]) for j in range(4)],
            "retained_component_weights": [float(value) for value in actual_weights],
            "local_value_scale": float(value_scale[position]),
            "local_derivative_scale": float(derivative_scale[position]),
        }
        if best is None or candidate["normalized_contact_residual"] < best["normalized_contact_residual"]:
            best = candidate

        for x_position, x in enumerate(xs):
            vertices = sparse_contact_vertices(fine[:, :, x_position])
            for vertex in vertices:
                sparse = np.zeros(len(indices))
                sparse[np.asarray(vertex["support_positions"])] = vertex["weights"]
                television = 0.5 * float(np.sum(np.abs(sparse - actual_weights)))
                if (
                    sparse_contact_upper_bound is None
                    or television < sparse_contact_upper_bound["total_variation"]
                ):
                    sparse_contact_upper_bound = {
                        "time": time,
                        "x": float(x),
                        "total_variation": television,
                        "support_positions": vertex["support_positions"],
                        "support_indices": [
                            indices[position]
                            for position in vertex["support_positions"]
                        ],
                        "weights": vertex["weights"],
                        "h2": vertex["h2"],
                        "h3": vertex["h3"],
                    }

    tail_bounds = {
        str(order): component_index_tail_bound(INTEGER_TAIL_START, order)
        for order in range(4)
    }
    return {
        "retained_indices": list(indices),
        "omitted_indices": f"n>={INTEGER_TAIL_START}",
        "minimum_sampled_fixed_weight_contact_residual": best,
        "nearest_explicit_sparse_contact_total_variation_upper_bound": sparse_contact_upper_bound,
        "component_index_tail_moment_bounds": tail_bounds,
        "maximum_coarse_fine_difference": quadrature_max_difference,
        "proof_boundary": (
            "Finite-grid diagnostic only. The retained fixed-weight sum and its "
            "tail bounds do not certify a continuum rectangle or exclude a contact."
        ),
    }


def targeted_integer_refinement() -> list[dict]:
    roster = ROSTERS[-1]
    configurations = (
        {
            "target": -1.0,
            "time_min": 0.0,
            "time_max": 0.02,
            "time_step": 0.0025,
            "x_min": 48.0,
            "x_max": 50.0,
            "x_step": 0.025,
        },
        {
            "target": 20.0,
            "time_min": 0.0,
            "time_max": 0.02,
            "time_step": 0.0025,
            "x_min": 78.0,
            "x_max": 82.0,
            "x_step": 0.025,
        },
    )
    results: list[dict] = []
    for config in configurations:
        times = np.arange(
            config["time_min"],
            config["time_max"] + config["time_step"] / 2.0,
            config["time_step"],
        )
        xs = np.arange(
            config["x_min"],
            config["x_max"] + config["x_step"] / 2.0,
            config["x_step"],
        )
        best_interval: dict | None = None
        best_witness: dict | None = None
        strongest_pole: dict | None = None
        maximum_quadrature_drift = 0.0
        for time in times:
            coarse, _ = normalized_component_jets(
                roster.indices, float(time), xs, COARSE_ORDER
            )
            fine, _ = normalized_component_jets(
                roster.indices, float(time), xs, FINE_ORDER
            )
            maximum_quadrature_drift = max(
                maximum_quadrature_drift,
                float(np.max(np.abs(fine - coarse))),
            )
            for position, x in enumerate(xs):
                jets_at_x = fine[:, :, position]
                witness = solve_target_witness(
                    jets_at_x, float(x), config["target"]
                )
                if witness is not None and (
                    best_witness is None or witness["score"] > best_witness["score"]
                ):
                    witness["time"] = float(time)
                    witness["x"] = float(x)
                    witness["support_indices"] = [
                        roster.indices[index]
                        for index in witness["support_positions"]
                    ]
                    best_witness = witness

                vertices = sparse_contact_vertices(jets_at_x)
                if not vertices:
                    continue
                second_jets = np.asarray([vertex["h2"] for vertex in vertices])
                if np.all(second_jets > 1.0e-12) or np.all(second_jets < -1.0e-12):
                    fields = np.asarray(
                        [
                            float(x) * vertex["h3"] / (6.0 * vertex["h2"]) - 0.5
                            for vertex in vertices
                        ]
                    )
                    lower = float(np.min(fields))
                    upper = float(np.max(fields))
                    target = config["target"]
                    distance = (
                        0.0
                        if lower <= target <= upper
                        else min(abs(target - lower), abs(target - upper))
                    )
                    candidate = {
                        "distance": distance,
                        "time": float(time),
                        "x": float(x),
                        "field_min": lower,
                        "field_max": upper,
                    }
                    if best_interval is None or distance < best_interval["distance"]:
                        best_interval = candidate
                else:
                    pole = unbounded_ratio_guard(vertices)
                    if pole["witness"] is not None and (
                        strongest_pole is None
                        or pole["witness"]["h3_at_h2_zero_abs"]
                        > strongest_pole["witness"]["h3_at_h2_zero_abs"]
                    ):
                        pole["time"] = float(time)
                        pole["x"] = float(x)
                        for side in ("left", "right"):
                            pole["witness"][side]["support_indices"] = [
                                roster.indices[index]
                                for index in pole["witness"][side]["support_positions"]
                            ]
                        strongest_pole = pole
        results.append(
            {
                "roster": roster.name,
                "configuration": config,
                "time_rows": len(times),
                "x_rows": len(xs),
                "evaluated_points": len(times) * len(xs),
                "best_target_witness": best_witness,
                "nearest_definite_interval": best_interval,
                "strongest_pole_candidate": strongest_pole,
                "maximum_coarse_fine_jet_difference": maximum_quadrature_drift,
                "proof_boundary": (
                    "A finite local refinement can falsify the coarse-grid gap but "
                    "cannot certify its persistence on a continuum domain."
                ),
            }
        )
    return results


def build_payload(baseline: list[float], priority_lowered: bool) -> dict:
    xs = np.arange(X_MIN, X_MAX + X_STEP / 2.0, X_STEP)
    started = perf_counter()
    roster_results: list[dict] = []
    runtime_cpu_samples: list[float] = []
    consecutive_high_samples = 0
    resource_parked = False
    for roster in ROSTERS:
        roster_results.append(scan_roster(roster, xs))
        sample = float(psutil.cpu_percent(interval=0.25))
        runtime_cpu_samples.append(sample)
        if sample > 75.0:
            consecutive_high_samples += 1
        else:
            consecutive_high_samples = 0
        if consecutive_high_samples >= 2:
            resource_parked = True
            break
    fixed_xs = xs[xs <= FIXED_WEIGHT_X_MAX]
    targeted_refinement = None
    fixed = None
    if not resource_parked:
        targeted_refinement = targeted_integer_refinement()
        fixed = fixed_theta_diagnostic(fixed_xs)
    cutoff_bounds = {
        roster.name: {
            str(order): {
                "log10_upper": max(
                    cutoff_tail_log10_bound(index, order)
                    for index in roster.indices
                )
            }
            for order in range(4)
        }
        for roster in ROSTERS
    }
    elapsed = perf_counter() - started
    return {
        "kind": STEM,
        "schema_version": SCHEMA_VERSION,
        "date": DATE,
        "status": "finite adversarial m=2 theta-contact ladder scout",
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "resource_policy": {
            "mode": "daytime",
            "active_compute_workers": 1,
            "below_normal_priority_applied": priority_lowered,
            "baseline_cpu_percent": baseline,
            "baseline_mean_percent": sum(baseline) / len(baseline),
            "runtime_cpu_percent": runtime_cpu_samples,
            "resource_parked": resource_parked,
            "elapsed_seconds": elapsed,
        },
        "exact_setup": {
            "component": (
                "phi_a(u)=pi*a^2*exp(5u)*(2*pi*a^2*exp(4u)-3)"
                "*exp(-pi*a^2*exp(4u))"
            ),
            "translate_law": "phi_a(u)=a^(-1/2)*phi_1(u+(log a)/2)",
            "normalized_jets": (
                "h_a^(j)(c)=M_a^(-1)*integral_0^infinity "
                "exp(tu^2)phi_a(u)*partial_c^j cos(cu)du"
            ),
            "m2_contact": "sum p_a h_a=0 and sum p_a h_a'=0",
            "field": "ell=c*(sum p_a h_a''')/(6*sum p_a h_a'')-1/2",
            "division_free_target": (
                "c*sum p_a h_a'''-6*(ell+1/2)*sum p_a h_a''=0"
            ),
            "vertex_fact": (
                "With normalization and two contact equations, every vertex "
                "has support at most three; a prescribed field adds one linear "
                "equation and needs support at most four."
            ),
        },
        "grid": {
            "times": list(TIMES),
            "x_min": X_MIN,
            "x_max": X_MAX,
            "fixed_weight_x_max": FIXED_WEIGHT_X_MAX,
            "x_step": X_STEP,
            "x_rows": len(xs),
            "cutoff": CUTOFF,
            "coarse_order": COARSE_ORDER,
            "fine_order": FINE_ORDER,
        },
        "analytic_truncation_envelopes": {
            "cutoff_tail_moment_bounds": cutoff_bounds,
            "derivation": (
                "For t<=1/5, t*u^2<=exp(4u)/(20e^2), u^j<=j!*exp(u), "
                "and phi_a<=2*pi^2*a^4*exp(9u)*exp(-pi*a^2*exp(4u)). "
                "After y=exp(4u), log(y/y0)<=(y-y0)/y0 gives the stored "
                "u>=2 envelope. The same argument from y=1 gives the integer tail."
            ),
        },
        "roster_results": roster_results,
        "targeted_integer_refinement": targeted_refinement,
        "fixed_theta_diagnostic": fixed,
        "conclusion": {
            "free_roster_test": (
                "A roster passes this diagnostic only when every stored target field "
                "has a positive division-free four-component witness."
            ),
            "fixed_weight_test": (
                "The exact theta weights are not varied. The sampled residual and "
                "sparse-contact distance only measure the finite scout geometry."
            ),
            "proof_boundary": (
                "Positive point witnesses are high-accuracy quadrature diagnostics, "
                "not interval certificates. No fixed-weight contact exclusion, Xi "
                "field bound, Lambda<=0, RH, or prize-level conclusion is proved."
            ),
        },
    }


def render_note(payload: dict) -> str:
    lines = [
        "# Discrete Theta Contact Ladder Scout",
        "",
        f"Date: {DATE}",
        "",
        "Status: finite adversarial diagnostic. This is not a proof of `Lambda<=0` or RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Division-Free Test",
        "",
        "At multiplicity two, normalized positive component weights satisfy",
        "",
        "```text",
        "sum p_a h_a=0,  sum p_a h_a'=0,",
        "c sum p_a h_a'''-6(ell+1/2)sum p_a h_a''=0.",
        "```",
        "",
        "Together with `sum p_a=1`, this is a four-equation linear witness problem.",
        "The adjacent-jet ratio is never divided during the search.",
        "",
        "## Ladder Results",
        "",
        "| roster | targets realized | unbounded guard | sampled definite field span | max jet drift |",
        "|:---|---:|:---:|:---|---:|",
    ]
    for result in payload["roster_results"]:
        unbounded = result["unbounded_adjacent_ratio_guard"]
        geometry = result["sampled_contact_polytope_geometry"]
        field_span = (
            f"[{geometry['sampled_vertex_field_min']:.3g},"
            f" {geometry['sampled_vertex_field_max']:.3g}]"
            if geometry["sampled_vertex_field_min"] is not None
            else "none"
        )
        lines.append(
            f"| {result['name']} | "
            f"{result['target_fields_realized']}/{len(result['target_fields'])} | "
            f"{bool(unbounded and unbounded['detected'])} | "
            f"{field_span} | "
            f"{result['maximum_coarse_fine_jet_difference']:.3e} |"
        )
    lines.extend(
        [
            "",
            "The target set is `-20,-5,-1,0,1,5,20`. Each realized value has an explicit",
            "positive four-component witness in the JSON result.",
            "",
            "## Integer Gap Refinement",
            "",
        ]
    )
    refinements = payload["targeted_integer_refinement"]
    if refinements is not None:
        for refinement in refinements:
            target = refinement["configuration"]["target"]
            nearest_interval = refinement["nearest_definite_interval"]
            witness_status = refinement["best_target_witness"] is not None
            pole = refinement["strongest_pole_candidate"]
            lines.append(
                f"- `ell={target:g}`: target witness `{witness_status}`; "
                f"nearest definite interval distance "
                f"`{nearest_interval['distance']:.6e}`; pole candidate "
                f"`{bool(pole and pole['detected'])}` over "
                f"`{refinement['evaluated_points']}` refined points."
            )
        lines.extend(
            [
                "",
                "These are local finite refinements. Failure to find a witness is not",
                "a continuum exclusion theorem.",
                "",
            ]
        )
    lines.extend(
        [
            "## Fixed Theta Rung",
            "",
        ]
    )
    fixed = payload["fixed_theta_diagnostic"]
    if fixed is None:
        lines.extend(
            [
                "The run parked after two consecutive high-CPU samples. The fixed",
                "theta rung remains unfinished and must be resumed in a fresh cycle.",
                "",
                "## Proof Boundary",
                "",
                payload["conclusion"]["proof_boundary"],
                "",
            ]
        )
        return "\n".join(lines)
    nearest = fixed["minimum_sampled_fixed_weight_contact_residual"]
    sparse = fixed["nearest_explicit_sparse_contact_total_variation_upper_bound"]
    lines.extend(
        [
            f"The smallest sampled normalized fixed-weight contact residual was "
            f"`{nearest['normalized_contact_residual']:.6e}` at "
            f"`t={nearest['time']}, x={nearest['x']}`.",
            "",
            f"The nearest stored sparse free-weight contact has total-variation "
            f"distance at most `{sparse['total_variation']:.6e}` from the retained "
            "fixed mixture. This is an upper bound from explicit contact vertices,",
            "not the exact distance to the full contact polytope.",
            "",
            "The result stores analytic envelopes for both the omitted `u>2` integral",
            "and the fixed-weight integer tail `n>=9`. Quadrature agreement and finite",
            "grid coverage remain diagnostic rather than interval-certified.",
            "",
            "## Decision Gate",
            "",
            "After local refinement, the free integer roster realizes all seven tested",
            "fields. Thus the sampled evidence rejects the integer shift locations alone",
            "as the missing rigidity. This does not prove realization of every real field.",
            "The next local theorem must use the exact unit theta coefficients, or the",
            "equivalent time-dependent normalized mass profile, jointly with the contact",
            "equations and the infinite roster. A global zero-flow bypass remains separate.",
            "",
            "## Proof Boundary",
            "",
            payload["conclusion"]["proof_boundary"],
            "",
        ]
    )
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--baseline-seconds", type=int, default=5)
    parser.add_argument("--baseline-max-percent", type=float, default=60.0)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    priority_lowered = request_below_normal_priority()
    baseline = sample_baseline(args.baseline_seconds)
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "theta-contact ladder baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, "
        f"mean={baseline_mean:.1f}%"
    )
    if baseline_mean > args.baseline_max_percent:
        print("theta-contact ladder deferred: daytime baseline is already busy")
        return 2
    payload = build_payload(baseline, priority_lowered)
    write_json_atomic(args.out, payload)
    write_text_atomic(args.note, render_note(payload))
    realized = [
        result["target_fields_realized"] for result in payload["roster_results"]
    ]
    print(
        "built discrete theta contact ladder scout: "
        f"free target counts={realized}, "
        f"elapsed={payload['resource_policy']['elapsed_seconds']:.3f}s"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
