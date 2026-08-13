#!/usr/bin/env python3
"""Project the genuine theta mass profile onto the m=2 contact slice."""

from __future__ import annotations

import argparse
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

import numpy as np
import psutil
from scipy.optimize import linprog

import jensen_window_pf_newman_discrete_theta_contact_ladder_scout as ladder


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_fixed_theta_mass_profile_dual_scout"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
PARENT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_discrete_theta_contact_ladder_scout.json"
)
SCHEMA_VERSION = 1
DATE = "2026-08-03"

DEPTHS = (4, 8, 12)
TIMES = (0.0, 0.05, 0.1, 0.15, 0.2)
X_MIN = 2.0
X_MAX = 160.0
X_STEP = 0.5
REFINE_TIME_RADIUS = 0.02
REFINE_TIME_STEP = 0.0025
REFINE_X_RADIUS = 1.0
REFINE_X_STEP = 0.025
TOTAL_MASS_LOWER_BOUND = 9.0 / 2900.0


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


def sample_cpu(seconds: int) -> list[float]:
    return [float(psutil.cpu_percent(interval=1.0)) for _ in range(seconds)]


def contact_matrix(
    weights: np.ndarray, jets_at_x: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    absolute_features = np.abs(jets_at_x[:, :2])
    maximum_features = np.max(absolute_features, axis=0)
    profile_residuals = np.abs(weights @ jets_at_x[:, :2])
    # In perturbation coordinates the scaled right-hand side is nonzero.  The
    # cap keeps feature coefficients below 1e12 when a residual is minute.
    scales = np.maximum(profile_residuals, 1.0e-12 * maximum_features)
    scales = np.maximum(scales, 1.0e-300)
    value_scale, derivative_scale = (float(value) for value in scales)
    matrix = np.vstack(
        (
            np.ones(jets_at_x.shape[0]),
            jets_at_x[:, 0] / value_scale,
            jets_at_x[:, 1] / derivative_scale,
        )
    )
    return matrix, np.asarray((value_scale, derivative_scale))


def tv_projection(weights: np.ndarray, jets_at_x: np.ndarray) -> dict | None:
    """Solve the exact finite-dimensional total-variation projection LP."""
    count = len(weights)
    matrix, scales = contact_matrix(weights, jets_at_x)
    right_hand_side = np.asarray((1.0, 0.0, 0.0))

    # Write delta=u-v.  This exposes -A p as the equality right-hand side;
    # formulating directly in q lets a tiny nonzero defect be lost inside an
    # absolute feasibility tolerance when the target right-hand side is zero.
    objective = 0.5 * np.ones(2 * count)
    equality = np.zeros((3, 2 * count))
    equality[:, :count] = matrix
    equality[:, count:] = -matrix
    equality_rhs = np.asarray((0.0, -np.dot(weights, matrix[1]), -np.dot(weights, matrix[2])))
    inequality = np.zeros((count, 2 * count))
    inequality[:, :count] = -np.eye(count)
    inequality[:, count:] = np.eye(count)
    inequality_rhs = weights.copy()

    result = linprog(
        objective,
        A_ub=inequality,
        b_ub=inequality_rhs,
        A_eq=equality,
        b_eq=equality_rhs,
        bounds=(0.0, None),
        method="highs-ds",
        options={
            "primal_feasibility_tolerance": 1.0e-10,
            "dual_feasibility_tolerance": 1.0e-10,
        },
    )
    if not result.success:
        return None

    positive = np.asarray(result.x[:count])
    negative = np.asarray(result.x[count:])
    perturbation = positive - negative
    projected = weights + perturbation
    equality_dual_scaled = np.asarray(result.eqlin.marginals)
    equality_dual_original = equality_dual_scaled.copy()
    equality_dual_original[1:] /= scales
    inequality_dual = np.asarray(result.ineqlin.marginals)
    lower_dual = np.asarray(result.lower.marginals)
    stationarity = (
        objective
        - equality.T @ equality_dual_scaled
        - inequality.T @ inequality_dual
        - lower_dual
    )
    dual_objective = float(
        equality_rhs @ equality_dual_scaled
        + inequality_rhs @ inequality_dual
    )
    actual_tv = 0.5 * float(np.sum(np.abs(projected - weights)))
    return {
        "total_variation": actual_tv,
        "solver_objective": float(result.fun),
        "projected_weights": projected.tolist(),
        "weight_perturbation": perturbation.tolist(),
        "positive_perturbation_variables": positive.tolist(),
        "negative_perturbation_variables": negative.tolist(),
        "contact_row_scales": scales.tolist(),
        "unscaled_profile_contact_residuals": [
            float(np.dot(weights, jets_at_x[:, 0])),
            float(np.dot(weights, jets_at_x[:, 1])),
        ],
        "scaled_perturbation_right_hand_side": equality_rhs.tolist(),
        "scaled_contact_residual_inf": float(
            np.max(np.abs(matrix @ projected - right_hand_side))
        ),
        "unscaled_contact_residuals": [
            float(np.dot(projected, jets_at_x[:, 0])),
            float(np.dot(projected, jets_at_x[:, 1])),
        ],
        "equality_dual_scaled": equality_dual_scaled.tolist(),
        "equality_dual_original": equality_dual_original.tolist(),
        "dual_objective": dual_objective,
        "primal_dual_gap": float(result.fun - dual_objective),
        "stationarity_residual_inf": float(np.max(np.abs(stationarity))),
        "minimum_projected_weight": float(np.min(projected)),
        "solver_status": int(result.status),
    }


def fisher_projection(weights: np.ndarray, jets_at_x: np.ndarray) -> dict:
    """Return the closed-form affine chi-square projection and its dual."""
    matrix, scales = contact_matrix(weights, jets_at_x)
    features = matrix[1:].T
    mean = weights @ features
    centered = features - mean
    covariance = (centered.T * weights) @ centered
    eigenvalues = np.linalg.eigvalsh(covariance)
    condition = (
        float(eigenvalues[-1] / eigenvalues[0])
        if eigenvalues[0] > 0.0
        else math.inf
    )
    try:
        feature_dual_scaled = np.linalg.solve(covariance, mean)
    except np.linalg.LinAlgError:
        feature_dual_scaled = np.linalg.pinv(covariance, rcond=1.0e-15) @ mean
    dual_observable = centered @ feature_dual_scaled
    perturbation = -weights * dual_observable
    projected = weights + perturbation
    distance_squared = float(mean @ feature_dual_scaled)
    feature_dual_original = feature_dual_scaled / scales
    normalization_dual = -float(mean @ feature_dual_scaled)
    residual = matrix @ projected - np.asarray((1.0, 0.0, 0.0))
    return {
        "chi_square_distance_squared": distance_squared,
        "projected_weights": projected.tolist(),
        "weight_perturbation": perturbation.tolist(),
        "contact_row_scales": scales.tolist(),
        "scaled_feature_mean": mean.tolist(),
        "scaled_feature_covariance": covariance.tolist(),
        "covariance_eigenvalues": eigenvalues.tolist(),
        "covariance_condition_number": condition,
        "normalization_dual": normalization_dual,
        "feature_dual_scaled": feature_dual_scaled.tolist(),
        "feature_dual_original": feature_dual_original.tolist(),
        "dual_observable": dual_observable.tolist(),
        "scaled_contact_residual_inf": float(np.max(np.abs(residual))),
        "minimum_projected_weight": float(np.min(projected)),
        "positivity_active": bool(np.min(projected) < -1.0e-11),
    }


def evaluate_depth(
    depth: int,
    time: float,
    xs: np.ndarray,
    jets: np.ndarray,
    masses: np.ndarray,
) -> dict:
    retained_jets = jets[:depth]
    retained_masses = masses[:depth]
    weights = retained_masses / float(np.sum(retained_masses))
    best_tv: dict | None = None
    best_fisher: dict | None = None
    feasible_points = 0
    fisher_positive_points = 0

    for position, x in enumerate(xs):
        local_jets = retained_jets[:, :, position]
        fisher = fisher_projection(weights, local_jets)
        if not fisher["positivity_active"]:
            fisher_positive_points += 1
        fisher_candidate = {
            "depth": depth,
            "time": float(time),
            "x": float(x),
            "retained_profile_weights": weights.tolist(),
            "component_jets": local_jets.tolist(),
            **fisher,
        }
        if (
            best_fisher is None
            or fisher_candidate["chi_square_distance_squared"]
            < best_fisher["chi_square_distance_squared"]
        ):
            best_fisher = fisher_candidate

        television = tv_projection(weights, local_jets)
        if television is None:
            continue
        feasible_points += 1
        tv_candidate = {
            "depth": depth,
            "time": float(time),
            "x": float(x),
            "retained_profile_weights": weights.tolist(),
            "component_jets": local_jets.tolist(),
            **television,
            "fisher_projection_at_same_point": fisher,
        }
        if best_tv is None or tv_candidate["total_variation"] < best_tv["total_variation"]:
            best_tv = tv_candidate

    return {
        "depth": depth,
        "time": float(time),
        "x_points": len(xs),
        "tv_feasible_points": feasible_points,
        "fisher_positive_points": fisher_positive_points,
        "nearest_total_variation_contact": best_tv,
        "nearest_affine_fisher_contact": best_fisher,
    }


def evaluate_infinite_endpoint_closure(
    time: float,
    xs: np.ndarray,
    jets: np.ndarray,
    masses: np.ndarray,
) -> dict:
    """Enlarge the depth-12 contact polytope by the exact n->infinity limit."""
    weights = masses[: max(DEPTHS)] / float(np.sum(masses[: max(DEPTHS)]))
    extended_weights = np.concatenate((weights, np.zeros(1)))
    endpoint = np.zeros((1, 4, len(xs)))
    endpoint[:, 0, :] = 1.0
    extended_jets = np.concatenate((jets[: max(DEPTHS)], endpoint), axis=0)
    best: dict | None = None
    feasible_points = 0
    for position, x in enumerate(xs):
        local_jets = extended_jets[:, :, position]
        television = tv_projection(extended_weights, local_jets)
        if television is None:
            continue
        feasible_points += 1
        candidate = {
            "depth": "12_plus_infinity_endpoint",
            "time": float(time),
            "x": float(x),
            "retained_profile_weights": extended_weights.tolist(),
            "component_jets": local_jets.tolist(),
            "endpoint": {
                "label": "n=infinity",
                "normalized_jets": [1.0, 0.0, 0.0, 0.0],
                "justification": (
                    "The normalized theta component concentrates at u=0 as n tends "
                    "to infinity, so (h,h',h'',h''') tends to (1,0,0,0)."
                ),
            },
            **television,
        }
        if best is None or candidate["total_variation"] < best["total_variation"]:
            best = candidate
    return {
        "time": float(time),
        "x_points": len(xs),
        "tv_feasible_points": feasible_points,
        "nearest_total_variation_contact": best,
    }


def update_best(current: dict | None, candidate: dict | None, key: str) -> dict | None:
    if candidate is None:
        return current
    if current is None or candidate[key] < current[key]:
        return candidate
    return current


def tail_budget(depth: int) -> dict:
    raw = {
        str(order): ladder.component_index_tail_bound(depth + 1, order)
        for order in range(4)
    }
    return {
        "first_omitted_index": depth + 1,
        "raw_moment_bounds": raw,
        "normalized_omitted_mass_upper_bound": raw["0"] / TOTAL_MASS_LOWER_BOUND,
        "total_mass_lower_bound": TOTAL_MASS_LOWER_BOUND,
    }


def rounded_grid(start: float, stop: float, step: float) -> np.ndarray:
    count = int(round((stop - start) / step))
    return np.asarray([round(start + index * step, 12) for index in range(count + 1)])


def build_payload(baseline: list[float], priority_lowered: bool) -> dict:
    started = perf_counter()
    indices = tuple(float(index) for index in range(1, max(DEPTHS) + 1))
    main_xs = rounded_grid(X_MIN, X_MAX, X_STEP)
    coarse_fine_max_difference = 0.0
    main_rows: list[dict] = []
    best_tv_by_depth: dict[int, dict | None] = {depth: None for depth in DEPTHS}
    best_fisher_by_depth: dict[int, dict | None] = {depth: None for depth in DEPTHS}
    closure_best: dict | None = None
    closure_main_rows: list[dict] = []
    runtime_cpu_samples: list[float] = []
    consecutive_high = 0
    parked = False

    for time in TIMES:
        coarse_jets, coarse_masses = ladder.normalized_component_jets(
            indices, time, main_xs, ladder.COARSE_ORDER
        )
        fine_jets, fine_masses = ladder.normalized_component_jets(
            indices, time, main_xs, ladder.FINE_ORDER
        )
        coarse_fine_max_difference = max(
            coarse_fine_max_difference,
            float(np.max(np.abs(coarse_jets - fine_jets))),
            float(np.max(np.abs(coarse_masses - fine_masses))),
        )
        for depth in DEPTHS:
            row = evaluate_depth(depth, time, main_xs, fine_jets, fine_masses)
            main_rows.append(row)
            best_tv_by_depth[depth] = update_best(
                best_tv_by_depth[depth],
                row["nearest_total_variation_contact"],
                "total_variation",
            )
            best_fisher_by_depth[depth] = update_best(
                best_fisher_by_depth[depth],
                row["nearest_affine_fisher_contact"],
                "chi_square_distance_squared",
            )
        closure_row = evaluate_infinite_endpoint_closure(
            time, main_xs, fine_jets, fine_masses
        )
        closure_main_rows.append(closure_row)
        closure_best = update_best(
            closure_best,
            closure_row["nearest_total_variation_contact"],
            "total_variation",
        )
        sample = float(psutil.cpu_percent(interval=1.0))
        runtime_cpu_samples.append(sample)
        consecutive_high = consecutive_high + 1 if sample > 75.0 else 0
        if consecutive_high >= 2:
            parked = True
            break

    refinement_rows: list[dict] = []
    closure_refinement_rows: list[dict] = []
    refinement_anchor = best_tv_by_depth[max(DEPTHS)]
    if not parked and refinement_anchor is not None:
        time_start = max(0.0, refinement_anchor["time"] - REFINE_TIME_RADIUS)
        time_stop = min(0.2, refinement_anchor["time"] + REFINE_TIME_RADIUS)
        x_start = max(X_MIN, refinement_anchor["x"] - REFINE_X_RADIUS)
        x_stop = min(X_MAX, refinement_anchor["x"] + REFINE_X_RADIUS)
        refine_times = rounded_grid(time_start, time_stop, REFINE_TIME_STEP)
        refine_xs = rounded_grid(x_start, x_stop, REFINE_X_STEP)
        for time in refine_times:
            fine_jets, fine_masses = ladder.normalized_component_jets(
                indices, float(time), refine_xs, ladder.FINE_ORDER
            )
            for depth in DEPTHS:
                row = evaluate_depth(
                    depth, float(time), refine_xs, fine_jets, fine_masses
                )
                refinement_rows.append(row)
                best_tv_by_depth[depth] = update_best(
                    best_tv_by_depth[depth],
                    row["nearest_total_variation_contact"],
                    "total_variation",
                )
                best_fisher_by_depth[depth] = update_best(
                    best_fisher_by_depth[depth],
                    row["nearest_affine_fisher_contact"],
                    "chi_square_distance_squared",
                )
            closure_row = evaluate_infinite_endpoint_closure(
                float(time), refine_xs, fine_jets, fine_masses
            )
            closure_refinement_rows.append(closure_row)
            closure_best = update_best(
                closure_best,
                closure_row["nearest_total_variation_contact"],
                "total_variation",
            )
            sample = float(psutil.cpu_percent(interval=1.0))
            runtime_cpu_samples.append(sample)
            consecutive_high = consecutive_high + 1 if sample > 75.0 else 0
            if consecutive_high >= 2:
                parked = True
                break

    budgets = {str(depth): tail_budget(depth) for depth in DEPTHS}
    depth_results = []
    for depth in DEPTHS:
        television = best_tv_by_depth[depth]
        tail = budgets[str(depth)]["normalized_omitted_mass_upper_bound"]
        depth_results.append(
            {
                "depth": depth,
                "nearest_total_variation_contact": television,
                "nearest_affine_fisher_contact": best_fisher_by_depth[depth],
                "tail_budget": budgets[str(depth)],
                "tv_to_tail_ratio": (
                    television["total_variation"] / tail
                    if television is not None and tail > 0.0
                    else None
                ),
            }
        )

    deepest = depth_results[-1]["nearest_total_variation_contact"]
    closure_objective = (
        closure_best["solver_objective"] if closure_best is not None else None
    )
    distance_resolved = bool(
        closure_objective is not None
        and closure_objective > 100.0 * coarse_fine_max_difference
    )
    elapsed = perf_counter() - started
    return {
        "kind": STEM,
        "schema_version": SCHEMA_VERSION,
        "date": DATE,
        "status": (
            "parked after resource threshold"
            if parked
            else "finite fixed-theta mass-profile primal-dual scout complete"
        ),
        "source_sha256": sha256_path(Path(__file__).resolve()),
        "parent": {
            "result": str(PARENT_RESULT.relative_to(REPO_ROOT)).replace("\\", "/"),
            "result_sha256": sha256_path(PARENT_RESULT),
            "builder_sha256": sha256_path(Path(ladder.__file__).resolve()),
        },
        "resource_policy": {
            "mode": "daytime",
            "active_compute_workers": 1,
            "thread_caps": 1,
            "below_normal_priority_applied": priority_lowered,
            "baseline_cpu_percent": baseline,
            "baseline_mean_percent": sum(baseline) / len(baseline),
            "runtime_cpu_percent": runtime_cpu_samples,
            "resource_parked": parked,
            "elapsed_seconds": elapsed,
        },
        "exact_setup": {
            "genuine_profile": "p_(n,t)=M_(n,t)/sum_k M_(k,t)",
            "contact_slice": "sum q_n=1, sum q_n h_n=0, sum q_n h_n'=0, q_n>=0",
            "total_variation_problem": "min_q (1/2)sum_n |q_n-p_n| on the contact slice",
            "fisher_problem": "min_delta sum_n delta_n^2/p_n subject to A(p+delta)=(1,0,0)",
            "fisher_dual": (
                "For scaled features f=(h/s0,h'/s1), mu=E_p f and "
                "C=Cov_p(f), alpha=C^(-1)mu, delta_n=-p_n(f_n-mu).alpha."
            ),
            "coefficient_boundary": (
                "Only p_(n,t) is projected. Contact weights are not refitted when "
                "evaluating the genuine profile; refitting defines the comparison set."
            ),
        },
        "grid": {
            "depths": list(DEPTHS),
            "times": list(TIMES),
            "x_min": X_MIN,
            "x_max": X_MAX,
            "x_step": X_STEP,
            "coarse_order": ladder.COARSE_ORDER,
            "fine_order": ladder.FINE_ORDER,
            "refinement_anchor": refinement_anchor,
            "refinement_time_radius": REFINE_TIME_RADIUS,
            "refinement_time_step": REFINE_TIME_STEP,
            "refinement_x_radius": REFINE_X_RADIUS,
            "refinement_x_step": REFINE_X_STEP,
        },
        "maximum_coarse_fine_difference": coarse_fine_max_difference,
        "numerical_resolution_gate": {
            "distance_values_resolved": distance_resolved,
            "smallest_locator_objective": closure_objective,
            "maximum_coarse_fine_jet_or_mass_difference": coarse_fine_max_difference,
            "required_objective_to_drift_factor": 100.0,
            "role": (
                "When false, the double-precision LP identifies candidate regions "
                "and dual directions only. Its distance value and active face must "
                "be replaced by an independent high-precision face solve."
            ),
        },
        "main_slice_rows": main_rows,
        "refinement_slice_rows": refinement_rows,
        "infinite_endpoint_closure": {
            "main_slice_rows": closure_main_rows,
            "refinement_slice_rows": closure_refinement_rows,
            "nearest_total_variation_contact": closure_best,
            "role": (
                "Enlarged finite contact polytope containing the exact limiting "
                "component direction; a positive distance here cannot be attributed "
                "to the largest retained index alone."
            ),
        },
        "depth_results": depth_results,
        "diagnostic_decision": {
            "deepest_nearest_tv": deepest,
            "question": (
                "Do the nearest distance, perturbation support, and dual direction "
                "stabilize before the analytic omitted-index budget?"
            ),
            "proof_boundary": (
                "The LPs and Fisher projections are finite high-accuracy diagnostics. "
                "The selected LP objective lies below the conservative quadrature "
                "resolution gate and is a locator, not a certified distance. The "
                "stored tail envelopes are analytic, but this scout does not interval-"
                "certify the quadrature, optimize over continuous (t,x), or turn a "
                "sampled dual into an infinite arithmetic separator."
            ),
        },
    }


def render_note(payload: dict) -> str:
    lines = [
        "# Fixed Theta Mass-Profile Dual Scout",
        "",
        f"Date: {DATE}",
        "",
        "Status: finite primal-dual diagnostic. This is not a proof of contact exclusion or RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Fixed-Coefficient Question",
        "",
        "The true unit-coefficient theta sum induces the normalized mass profile",
        "",
        "```text",
        "p_(n,t)=M_(n,t)/sum_k M_(k,t).",
        "```",
        "",
        "At each sampled `(t,x)` this scout solves the finite-dimensional",
        "total-variation projection from `p` to the positive contact slice and the",
        "closed-form affine Fisher projection. The latter returns a two-coordinate",
        "dual observable that can be tested for an arithmetic pattern in `n`.",
        "The smallest double-precision objective is below the conservative",
        "coarse/fine quadrature resolution gate, so the table is a locator rather",
        "than a certified distance report; the independent checker resolves the",
        "selected face at 80 decimal places.",
        "",
        "## Depth Stability",
        "",
        "| depth | locator TV objective | t | x | normalized omitted-mass bound | TV/tail | Fisher d^2 at TV point |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for result in payload["depth_results"]:
        nearest = result["nearest_total_variation_contact"]
        if nearest is None:
            lines.append(f"| {result['depth']} | infeasible | - | - | - | - | - |")
            continue
        tail = result["tail_budget"]["normalized_omitted_mass_upper_bound"]
        fisher = nearest["fisher_projection_at_same_point"]
        ratio = result["tv_to_tail_ratio"]
        lines.append(
            f"| {result['depth']} | {nearest['total_variation']:.6e} | "
            f"{nearest['time']:.6g} | {nearest['x']:.6g} | {tail:.3e} | "
            f"{ratio:.3e} | {fisher['chi_square_distance_squared']:.6e} |"
        )
    closure = payload["infinite_endpoint_closure"]["nearest_total_variation_contact"]
    if closure is not None:
        tail = payload["depth_results"][-1]["tail_budget"][
            "normalized_omitted_mass_upper_bound"
        ]
        lines.append(
            f"| 12 + infinity endpoint | {closure['total_variation']:.6e} | "
            f"{closure['time']:.6g} | {closure['x']:.6g} | {tail:.3e} | "
            f"{closure['total_variation'] / tail:.3e} | n/a |"
        )
    lines.extend(
        [
            "",
            "## Dual Audit",
            "",
        ]
    )
    deepest = payload["depth_results"][-1]["nearest_total_variation_contact"]
    if deepest is not None:
        fisher = deepest["fisher_projection_at_same_point"]
        lines.extend(
            [
                f"At the deepest nearest point, the TV primal-dual gap is "
                f"`{deepest['primal_dual_gap']:.3e}` and the scaled contact residual is "
                f"`{deepest['scaled_contact_residual_inf']:.3e}`.",
                "",
                "The Fisher feature dual in original `(h,h')` coordinates is",
                "",
                "```text",
                f"{fisher['feature_dual_original']}",
                "```",
                "",
                f"Its affine projection has minimum weight "
                f"`{fisher['minimum_projected_weight']:.3e}`; positivity-active is "
                f"`{fisher['positivity_active']}`.",
                "",
            ]
        )
    lines.extend(
        [
            "## Decision Gate",
            "",
            "A depth-stable numerical gap is not itself a theorem. A useful next step",
            "requires the dual direction to simplify into an identity or inequality",
            "that explicitly uses the unit theta coefficients and controls every",
            "omitted component. If the direction wanders with `(t,x)`, this local",
            "projection route should be deprioritized in favor of the global flow and",
            "degree-uniform Jensen branches.",
            "",
            "## Proof Boundary",
            "",
            payload["diagnostic_decision"]["proof_boundary"],
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
    baseline = sample_cpu(args.baseline_seconds)
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "fixed-theta dual baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, "
        f"mean={baseline_mean:.1f}%"
    )
    if baseline_mean > args.baseline_max_percent:
        print("fixed-theta dual scout deferred: daytime baseline is already busy")
        return 2
    payload = build_payload(baseline, priority_lowered)
    write_json_atomic(args.out, payload)
    write_text_atomic(args.note, render_note(payload))
    deepest = payload["depth_results"][-1]["nearest_total_variation_contact"]
    summary = "none" if deepest is None else f"{deepest['total_variation']:.6e}"
    print(
        "built fixed-theta mass-profile dual scout: "
        f"deepest nearest TV={summary}, "
        f"parked={payload['resource_policy']['resource_parked']}, "
        f"elapsed={payload['resource_policy']['elapsed_seconds']:.3f}s"
    )
    return 3 if payload["resource_policy"]["resource_parked"] else 0


if __name__ == "__main__":
    sys.exit(main())
