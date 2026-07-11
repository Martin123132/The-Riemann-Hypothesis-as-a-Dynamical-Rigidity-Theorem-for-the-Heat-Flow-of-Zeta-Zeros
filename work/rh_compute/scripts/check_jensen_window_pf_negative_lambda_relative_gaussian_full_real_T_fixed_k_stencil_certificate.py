#!/usr/bin/env python3
"""Validate the full-real-T fixed-k stencil certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate import (  # noqa: E402
    DEFAULT_EVENNESS_JSON,
    DEFAULT_FINITE_SEGMENT_JSON,
    DEFAULT_INDICES,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    build_artifact,
    result_line,
)


EXPECTED_KIND = (
    "jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate"
)
EXPECTED_IDS = [
    "nlrgfrtfksc_01_order42_source",
    "nlrgfrtfksc_02_full_disk_majorant",
    "nlrgfrtfksc_03_factored_compact_ray_bound",
    "nlrgfrtfksc_04_full_real_tail_majorant",
    "nlrgfrtfksc_05_scaled_upper_gamma_monotonicity",
    "nlrgfrtfksc_06_ray_residual_budget_certificate",
    "nlrgfrtfksc_07_ray_stencil_certificate",
    "nlrgfrtfksc_08_full_real_T_fixed_k_composition",
    "nlrgfrtfksc_09_all_k_promotion_rejected",
    "nlrgfrtfksc_10_acceptance_gate",
]


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_arb(text: object) -> flint.arb:
    return flint.arb(str(text).replace("E", "e"))


def issue(location: str, code: str, detail: str) -> dict[str, str]:
    return {"location": location, "code": code, "detail": detail}


def validate(
    artifact: dict,
    note_text: str,
    evenness_path: Path,
    finite_segment_path: Path,
) -> list[dict[str, str]]:
    issues: list[dict[str, str]] = []
    if artifact.get("kind") != EXPECTED_KIND:
        issues.append(issue("artifact", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("date") != "2026-07-10":
        issues.append(issue("artifact", "bad-date", repr(artifact.get("date"))))
    if artifact.get("status") != "full real-T fixed-k stencil certificate":
        issues.append(issue("artifact", "bad-status", repr(artifact.get("status"))))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "every t>=1156",
        "does not cover all k",
        "does not by itself prove cone entry",
        "does not prove rh",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("artifact", "weak-proof-boundary", required))

    rows = artifact.get("matrix_rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append(issue("matrix_rows", "bad-row-ids", repr([row.get("id") for row in rows])))
    if any(row.get("readiness") == "ready_to_apply" for row in rows):
        issues.append(issue("matrix_rows", "unexpected-ready-row", "ready_to_apply"))
    if rows and rows[5].get("role") != "real_interval_certificate":
        issues.append(issue("matrix_rows[5]", "missing-ray-certificate", repr(rows[5].get("role"))))
    if rows and rows[7].get("role") != "real_interval_certificate":
        issues.append(issue("matrix_rows[7]", "missing-full-T-certificate", repr(rows[7].get("role"))))
    if rows and rows[8].get("role") != "rejected_route":
        issues.append(issue("matrix_rows[8]", "missing-all-k-rejection", repr(rows[8].get("role"))))

    diagnostics = artifact.get("diagnostics", {})
    params = diagnostics.get("parameters", {})
    expected_params = {
        "full_T_interval": "T>=1156",
        "finite_segment": "1156<=T<=10000",
        "unbounded_ray": "T>=10000",
        "ray_start_T": 10000,
        "indices": list(DEFAULT_INDICES),
        "k": 22,
        "polynomial_M": 20,
        "first_residual_degree": 42,
        "disk_radius": "0.38",
        "split_x": "0.2",
        "split_y_at_ray_start": "400",
        "full_kernel_sum_n": 80,
        "ratio_cutoff_n": 80,
        "precision_bits": 4096,
        "value_budget_A": "0.5",
        "derivative_budget_B": "0.009",
    }
    for key, expected in expected_params.items():
        if params.get(key) != expected:
            issues.append(issue("parameters", f"bad-{key}", f"{params.get(key)!r} != {expected!r}"))
    for key in (
        "source_full_kernel_evenness_certified",
        "source_order42_residual_zero_certified",
        "source_finite_segment_certified",
        "full_disk_admissible",
        "full_real_majorants_decrease_after_split",
        "all_four_ray_residual_budgets_certified",
        "all_upper_gamma_scaled_tails_decrease_on_ray",
        "ray_fixed_k_stencil_system_certified",
        "full_real_T_fixed_k_stencil_system_certified",
    ):
        if diagnostics.get(key) is not True:
            issues.append(issue("diagnostics", f"false-{key}", repr(diagnostics.get(key))))
    if diagnostics.get("remaining_full_T_fixed_k_stencil_sources") != []:
        issues.append(
            issue(
                "diagnostics",
                "open-full-T-sources",
                repr(diagnostics.get("remaining_full_T_fixed_k_stencil_sources")),
            )
        )
    if diagnostics.get("target_closing") is not False:
        issues.append(issue("diagnostics", "target-closing", repr(diagnostics.get("target_closing"))))

    for key in (
        "full_disk_majorant_upper",
        "full_disk_degree4_tail_upper",
        "full_disk_degree2_tail_upper",
        "full_real_value_majorant_upper_at_split",
        "full_real_derivative_majorant_upper_at_split",
        "full_real_value_n_tail_upper",
        "full_real_derivative_n_tail_upper",
        "full_real_value_monotonicity_margin_lower",
        "full_real_derivative_monotonicity_margin_lower",
    ):
        try:
            value = parse_arb(diagnostics[key])
        except (KeyError, ValueError) as exc:
            issues.append(issue("diagnostics", f"unparseable-{key}", str(exc)))
        else:
            if not bool(value > 0):
                issues.append(issue("diagnostics", f"nonpositive-{key}", str(value)))
    for key in ("full_disk_degree4_ratio_upper", "full_disk_degree2_ratio_upper"):
        try:
            ratio = parse_arb(diagnostics[key])
        except (KeyError, ValueError) as exc:
            issues.append(issue("diagnostics", f"unparseable-{key}", str(exc)))
        else:
            if not bool(0 < ratio < 1):
                issues.append(issue("diagnostics", f"bad-{key}", str(ratio)))

    ray_rows = diagnostics.get("ray_residual_rows", [])
    if [row.get("index") for row in ray_rows] != list(DEFAULT_INDICES):
        issues.append(issue("ray_residual_rows", "bad-indices", repr(ray_rows)))
    for position, row in enumerate(ray_rows):
        location = f"ray_residual_rows[{position}]"
        for key in (
            "value_budget_certified_on_ray",
            "derivative_budget_certified_on_ray",
            "upper_gamma_scaled_tails_decrease_on_ray",
            "ray_residual_budget_certified",
        ):
            if row.get(key) is not True:
                issues.append(issue(location, f"false-{key}", repr(row.get(key))))
        if row.get("T_interval") != "T>=10000":
            issues.append(issue(location, "bad-T-interval", repr(row.get("T_interval"))))
        if row.get("compact_value_T_exponent") != -18:
            issues.append(issue(location, "bad-value-exponent", repr(row.get("compact_value_T_exponent"))))
        if row.get("compact_derivative_T_exponent") != -19:
            issues.append(issue(location, "bad-derivative-exponent", repr(row.get("compact_derivative_T_exponent"))))
        if row.get("selected_real_tail_moment_row_count") != 5:
            issues.append(issue(location, "bad-selected-tail-count", repr(row.get("selected_real_tail_moment_row_count"))))
        try:
            value_scaled = parse_arb(row["value_residual_scaled_uniform_upper"])
            derivative_scaled = parse_arb(row["derivative_residual_scaled_uniform_upper"])
            value_fraction = parse_arb(row["value_fraction_of_budget_upper"])
            derivative_fraction = parse_arb(row["derivative_fraction_of_budget_upper"])
            hazard_margin = parse_arb(
                row["upper_gamma_minimum_monotonicity_margin_lower"]
            )
        except (KeyError, ValueError) as exc:
            issues.append(issue(location, "unparseable", str(exc)))
            continue
        if not bool(0 < value_scaled < flint.arb("1e-14")):
            issues.append(issue(location, "large-value-bound", str(value_scaled)))
        if not bool(0 < derivative_scaled < flint.arb("1e-16")):
            issues.append(issue(location, "large-derivative-bound", str(derivative_scaled)))
        if not bool(0 < value_fraction < flint.arb("1e-14")):
            issues.append(issue(location, "large-value-fraction", str(value_fraction)))
        if not bool(0 < derivative_fraction < flint.arb("1e-14")):
            issues.append(issue(location, "large-derivative-fraction", str(derivative_fraction)))
        if not bool(hazard_margin > 300):
            issues.append(issue(location, "small-hazard-margin", str(hazard_margin)))

    perturbation_rows = diagnostics.get("perturbation_rows", [])
    if [row.get("name") for row in perturbation_rows] != [
        "normalizer",
        "B_product",
        "companion_product",
        "weighted_gap_derivative",
    ]:
        issues.append(issue("perturbation_rows", "bad-targets", repr(perturbation_rows)))
    for position, row in enumerate(perturbation_rows):
        if row.get("certified_positive_after_perturbation") is not True:
            issues.append(issue(f"perturbation_rows[{position}]", "uncertified", repr(row)))

    summary = artifact.get("summary", {})
    expected_counts = {
        "matrix_rows": 10,
        "ray_residual_rows": 4,
        "full_kernel_n_tail_channels": 4,
        "positive_retained_perturbation_margins": 4,
        "remaining_full_T_fixed_k_stencil_source_count": 0,
        "ready_to_apply_rows": 0,
    }
    for key, expected in expected_counts.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected!r}"))
    for key in (
        "ray_fixed_k_stencil_system_certified",
        "full_real_T_fixed_k_stencil_system_certified",
    ):
        if summary.get(key) is not True:
            issues.append(issue("summary", f"false-{key}", repr(summary.get(key))))
    if summary.get("target_closing") is not False:
        issues.append(issue("summary", "target-closing", repr(summary.get("target_closing"))))

    repo_root = Path(__file__).resolve().parents[3]
    for key, value in artifact.items():
        if key.startswith("source_") and isinstance(value, str):
            if not (repo_root / value).exists():
                issues.append(issue("sources", "missing-source", value))

    invariant_text = " ".join(str(value) for value in artifact.get("invariants", [])).lower()
    for required in (
        "cover every real t>=1156",
        "full-kernel evenness",
        "all full-kernel n tails",
        "proved decreasing",
        "not promoted to all k",
        "lambda <= 0",
    ):
        if required not in invariant_text:
            issues.append(issue("invariants", "missing-invariant", required))

    for required in (
        "full real-T fixed-k stencil certificate",
        "every real `T>=1156`",
        "order-42 residual zero",
        "T^-18",
        "[1156,10000] union [10000,infinity)",
        "does not prove RH or Lambda <= 0",
        result_line(artifact),
    ):
        if required.lower() not in note_text.lower():
            issues.append(issue("note", "missing-note-phrase", required))

    try:
        reproduced = build_artifact(evenness_path, finite_segment_path)
    except Exception as exc:  # pragma: no cover
        issues.append(issue("reproduction", "rebuild-failed", repr(exc)))
    else:
        for key in ("diagnostics", "matrix_rows", "summary"):
            if artifact.get(key) != reproduced.get(key):
                issues.append(issue("reproduction", f"{key}-mismatch", "static artifact differs from rebuild"))
    return issues


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--evenness-json", type=Path, default=DEFAULT_EVENNESS_JSON)
    parser.add_argument("--finite-segment-json", type=Path, default=DEFAULT_FINITE_SEGMENT_JSON)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    artifact = load_json(args.json)
    note_text = args.note.read_text(encoding="utf-8")
    issues = validate(artifact, note_text, args.evenness_json, args.finite_segment_json)
    if issues:
        print(json.dumps({"issues": issues}, indent=2))
        raise SystemExit(1)
    print(result_line(artifact))


if __name__ == "__main__":
    main()
