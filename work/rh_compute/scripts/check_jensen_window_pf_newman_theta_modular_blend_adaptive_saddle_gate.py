#!/usr/bin/env python3
"""Validate the modular-blend frequency-adaptive saddle gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate import (
    DEFAULT_OUT,
    X_VALUES,
    TIMES,
    build_exact,
    compute_adaptive_rows,
)


EXPECTED_IDS = [f"ntmbasg_{index:02d}_{suffix}" for index, suffix in (
    (1, "component_phase"),
    (2, "zero_time_saddle"),
    (3, "transition_scale"),
    (4, "positive_time_crossing"),
    (5, "uniform_scale"),
    (6, "adaptive_scout"),
    (7, "collar_guard"),
    (8, "live_handoff"),
)]


def validate(path: Path, rebuild: bool) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != "jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate":
        issues.append("artifact kind mismatch")
    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")
    build_exact()
    diagnostics = artifact.get("diagnostics", [])
    expected_keys = [(t, x) for t in TIMES for x in X_VALUES]
    actual_keys = [(row.get("t"), row.get("x")) for row in diagnostics]
    if actual_keys != expected_keys:
        issues.append(f"diagnostic key/order mismatch: {actual_keys}")
    expected_first = {
        ("0", 80): 3,
        ("0", 100): 4,
        ("0", 120): 4,
        ("0", 150): 5,
        ("0", 200): 6,
        ("0.5", 80): 3,
        ("0.5", 100): 4,
        ("0.5", 120): 4,
        ("0.5", 150): 5,
        ("0.5", 200): 6,
    }
    for row in diagnostics:
        key = (row["t"], row["x"])
        if row.get("first_n_relative_error_below_half") != expected_first[key]:
            issues.append(f"adaptive first-N drift at {key}")
        if row.get("first_n_margin_relative_error_below_half") != expected_first[key]:
            issues.append(f"adaptive margin first-N drift at {key}")
        if not row.get("adaptive_bound_covers_observed_n"):
            issues.append(f"two-index collar failed at {key}")
        if not row.get("adaptive_bound_covers_observed_margin_n"):
            issues.append(f"two-index margin collar failed at {key}")
        if mp.mpf(row.get("full_monotonicity_margin", "0")) <= 0:
            issues.append(f"full monotonicity margin is not positive at {key}")
        if mp.mpf(row["transition_index_a5"]) <= 0 or mp.mpf(row["transition_index_a9"]) <= 0:
            issues.append(f"nonpositive transition index at {key}")
    convergence = artifact.get("convergence", {})
    if convergence.get("stable_first_n") is not True:
        issues.append("coarse/fine first-N instability")
    if convergence.get("stable_margin_first_n") is not True:
        issues.append("coarse/fine margin first-N instability")
    if mp.mpf(convergence.get("max_relative_near_transition_ratio_delta", "1")) >= mp.mpf("1e-20"):
        issues.append("near-transition quadrature drift")
    if mp.mpf(
        convergence.get("max_relative_near_transition_margin_ratio_delta", "1")
    ) >= mp.mpf("1e-20"):
        issues.append("near-transition margin quadrature drift")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove",
        "Lambda<=0",
        "RH",
        "global sign is false",
        "Laguerre remainder sign",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    if rebuild and not issues:
        params = artifact["parameters"]
        fresh = compute_adaptive_rows(
            params["fine_nodes_per_panel"], params["dps"], params["theta_terms"]
        )
        for stored, rebuilt in zip(diagnostics, fresh, strict=True):
            if stored["first_n_relative_error_below_half"] != rebuilt["first_n_relative_error_below_half"]:
                issues.append(f"rebuild first-N drift at {(stored['t'], stored['x'])}")
            if (
                stored["first_n_margin_relative_error_below_half"]
                != rebuilt["first_n_margin_relative_error_below_half"]
            ):
                issues.append(f"rebuild margin first-N drift at {(stored['t'], stored['x'])}")
            n = stored["first_n_relative_error_below_half"]
            left = mp.mpf(stored["partials"][n - 1]["partial_to_full_ratio"])
            right = mp.mpf(rebuilt["partials"][n - 1]["partial_to_full_ratio"])
            if abs(left - right) / abs(right) >= mp.mpf("1e-35"):
                issues.append(f"rebuild ratio drift at {(stored['t'], stored['x'])}")
            margin_n = stored["first_n_margin_relative_error_below_half"]
            left_margin = mp.mpf(
                stored["partials"][margin_n - 1]["partial_margin_to_full_ratio"]
            )
            right_margin = mp.mpf(
                rebuilt["partials"][margin_n - 1]["partial_margin_to_full_ratio"]
            )
            if abs(left_margin - right_margin) / abs(right_margin) >= mp.mpf("1e-35"):
                issues.append(f"rebuild margin ratio drift at {(stored['t'], stored['x'])}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()
    issues = validate(args.artifact, args.rebuild)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman modular-blend adaptive-saddle gate: "
        f"{len(artifact['rows'])} rows, 0 issues, 2 exact saddle laws, "
        "1 square-root transition theorem, 10 adaptive diagnostics, "
        "10 matched monotonicity diagnostics, 1 collar guard, "
        "1 retired monotonicity branch"
    )


if __name__ == "__main__":
    main()
