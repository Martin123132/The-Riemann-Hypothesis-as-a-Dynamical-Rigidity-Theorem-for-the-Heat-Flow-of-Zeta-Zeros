#!/usr/bin/env python3
"""Validate the finite planar curvature-energy scaling scout."""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_planar_curvature_energy_scout.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_planar_curvature_energy_scout.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_planar_curvature_energy_scout.py"
)


def close(
    left: float,
    right: float,
    tolerance: float = 2e-11,
) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def check_schema(
    payload: dict,
    note_text: str,
    issues: list[str],
) -> None:
    if payload.get("kind") != (
        "jensen_window_pf_mertens_planar_curvature_energy_scout"
    ):
        issues.append("result kind mismatch")
    if payload.get("date") != "2026-07-24":
        issues.append("result date mismatch")
    if payload.get("status") != (
        "finite double-precision scaling diagnostic for the open "
        "planar curvature-energy gate"
    ):
        issues.append("result status mismatch")
    parameters = payload.get("parameters", {})
    if parameters.get("sizes") != [16, 32, 64, 128, 256, 512, 1024]:
        issues.append("size grid mismatch")
    if parameters.get("alphas") != [0.125, 0.25, 0.5, 0.75]:
        issues.append("alpha grid mismatch")
    rows = payload.get("rows", [])
    if len(rows) != 28:
        issues.append(f"expected 28 rows, found {len(rows)}")
    pairs = {(row.get("K"), row.get("alpha")) for row in rows}
    if len(pairs) != len(rows):
        issues.append("duplicate K/alpha rows")
    required_fields = {
        "K",
        "alpha",
        "R",
        "mixed_variation",
        "mixed_variation_times_K",
        "analytic_variation_bound",
        "variation_to_bound_ratio",
        "curvature_energy",
        "energy_per_K",
        "energy_per_K_log_2K",
        "max_planar_prefix",
        "max_prefix_per_K",
        "signed_offdiagonal",
        "double_abel_value",
        "direct_abel_abs_error",
        "elapsed_seconds",
        "energy_doubling_exponent",
        "variation_doubling_exponent",
    }
    for index, row in enumerate(rows):
        missing = required_fields - set(row)
        if missing:
            issues.append(f"row {index} missing fields {sorted(missing)}")
    proof_boundary = str(payload.get("proof_boundary", "")).lower()
    for marker in (
        "finite",
        "non-rigorous",
        "does not prove",
        "lambda <= 0",
    ):
        if marker not in proof_boundary:
            issues.append(f"proof boundary missing {marker}")
    for marker in (
        "finite double-precision scaling diagnostic",
        "This diagnostic does not prove",
        "E_(alpha,K)=O_epsilon(K^(1+epsilon))",
        "compatible with",
        "Finite Observation",
        "does not",
        "Lambda <= 0",
    ):
        if marker not in note_text:
            issues.append(f"note missing marker: {marker}")


def check_internal_metrics(payload: dict, issues: list[str]) -> None:
    rows = payload["rows"]
    audit = payload.get("audit", {})
    expected = {
        "row_count": len(rows),
        "max_K": max(row["K"] for row in rows),
        "alpha_count": len({row["alpha"] for row in rows}),
        "max_direct_abel_abs_error": max(
            row["direct_abel_abs_error"] for row in rows
        ),
        "max_variation_times_K": max(
            row["mixed_variation_times_K"] for row in rows
        ),
        "max_energy_per_K": max(
            row["energy_per_K"] for row in rows
        ),
        "max_prefix_per_K": max(
            row["max_prefix_per_K"] for row in rows
        ),
        "all_analytic_variation_bounds_pass": all(
            row["mixed_variation"] < row["analytic_variation_bound"]
            for row in rows
        ),
        "finite_energy_target_proved": False,
        "asymptotic_scaling_proved": False,
        "rh_proved": False,
        "lambda_le_zero_proved": False,
    }
    for key, value in expected.items():
        actual = audit.get(key)
        if isinstance(value, float):
            if not close(float(actual), value, 2e-12):
                issues.append(f"audit mismatch for {key}")
        elif actual != value:
            issues.append(f"audit mismatch for {key}")
    for row in rows:
        size = row["K"]
        if not close(
            row["mixed_variation_times_K"],
            size * row["mixed_variation"],
        ):
            issues.append(f"K*V mismatch at K={size}, alpha={row['alpha']}")
        if not close(
            row["energy_per_K"],
            row["curvature_energy"] / size,
        ):
            issues.append(f"E/K mismatch at K={size}, alpha={row['alpha']}")
        if not close(
            row["signed_offdiagonal"],
            row["double_abel_value"],
            2e-12,
        ):
            issues.append(
                f"direct/Abel mismatch at K={size}, alpha={row['alpha']}"
            )
        if row["direct_abel_abs_error"] >= 1e-12:
            issues.append(
                f"Abel error exceeds finite tolerance at "
                f"K={size}, alpha={row['alpha']}"
            )
        if row["mixed_variation"] >= row["analytic_variation_bound"]:
            issues.append(
                f"analytic variation bound failed at "
                f"K={size}, alpha={row['alpha']}"
            )
    if audit["max_variation_times_K"] >= 7.0:
        issues.append("recorded finite K*V maximum is not below 7")
    if audit["max_energy_per_K"] >= 1.0:
        issues.append("recorded finite E/K maximum is not below 1")
    if audit["max_prefix_per_K"] >= 3.0:
        issues.append("recorded finite M/K maximum is not below 3")


def check_selected_reproduction(
    payload: dict,
    issues: list[str],
) -> None:
    spec = importlib.util.spec_from_file_location(
        "planar_energy_scout_builder",
        BUILDER,
    )
    if spec is None or spec.loader is None:
        issues.append("could not load builder")
        return
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    selected = module.build_payload(
        max_size=64,
        alphas=(0.125, 0.5),
    )
    cached = {
        (row["K"], row["alpha"]): row for row in payload["rows"]
    }
    compare_fields = (
        "R",
        "mixed_variation",
        "curvature_energy",
        "max_planar_prefix",
        "signed_offdiagonal",
        "double_abel_value",
    )
    for row in selected["rows"]:
        key = (row["K"], row["alpha"])
        expected = cached.get(key)
        if expected is None:
            issues.append(f"selected row missing from cache: {key}")
            continue
        for field in compare_fields:
            left = row[field]
            right = expected[field]
            if isinstance(left, float):
                if not close(left, right, 4e-11):
                    issues.append(
                        f"selected reproduction mismatch {key} {field}"
                    )
            elif left != right:
                issues.append(
                    f"selected reproduction mismatch {key} {field}"
                )
    source = BUILDER.read_text(encoding="utf-8")
    for marker in (
        "OPENBLAS_NUM_THREADS",
        "SetPriorityClass",
        "0x00004000",
    ):
        if marker not in source:
            issues.append(f"resource-policy marker missing: {marker}")


def main() -> int:
    issues: list[str] = []
    for path in (RESULT, NOTE, BUILDER):
        if not path.exists():
            issues.append(f"missing file {path}")
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}", file=sys.stderr)
        return 1
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    note_text = NOTE.read_text(encoding="utf-8")
    check_schema(payload, note_text, issues)
    check_internal_metrics(payload, issues)
    check_selected_reproduction(payload, issues)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}", file=sys.stderr)
        return 1
    print(
        "validated Mertens planar curvature-energy scout: "
        "28 rows, 0 issues, K<=1024, 4 alpha values"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
