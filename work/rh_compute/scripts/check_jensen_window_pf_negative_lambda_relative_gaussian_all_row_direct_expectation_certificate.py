#!/usr/bin/env python3
"""Validate the all-row direct relative-Gaussian expectation certificate."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate import (  # noqa: E402
    DEFAULT_FIRST_OMITTED_JSON,
    DEFAULT_FLOATING_GRID_JSON,
    DEFAULT_INDICES,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_T_GRID,
    build_artifact,
    result_line,
)


EXPECTED_KIND = (
    "jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate"
)
EXPECTED_MATRIX_IDS = [
    "nlrgardec_01_uniform_split_geometry",
    "nlrgardec_02_degree384_exact_moment_core",
    "nlrgardec_03_rowwise_real_tail_certificate",
    "nlrgardec_04_global_n_tail_normalization",
    "nlrgardec_05_complete_recorded_grid_certificate",
    "nlrgardec_06_floating_grid_noninput_check",
    "nlrgardec_07_real_T_promotion_rejected",
    "nlrgardec_08_acceptance_gate",
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
    first_omitted_path: Path,
    floating_grid_path: Path,
) -> list[dict[str, str]]:
    issues: list[dict[str, str]] = []
    if artifact.get("kind") != EXPECTED_KIND:
        issues.append(issue("artifact", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("date") != "2026-07-10":
        issues.append(issue("artifact", "bad-date", repr(artifact.get("date"))))
    if artifact.get("status") != "complete recorded-grid direct expectation certificate":
        issues.append(issue("artifact", "bad-status", repr(artifact.get("status"))))

    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "20 recorded expectations",
        "does not certify any t interval",
        "does not prove rh",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("artifact", "weak-proof-boundary", required))

    rows = artifact.get("matrix_rows", [])
    if [row.get("id") for row in rows] != EXPECTED_MATRIX_IDS:
        issues.append(issue("matrix_rows", "bad-row-ids", repr([row.get("id") for row in rows])))
    if any(row.get("readiness") == "ready_to_apply" for row in rows):
        issues.append(issue("matrix_rows", "unexpected-ready-row", "ready_to_apply"))
    if rows and rows[4].get("readiness") != "available_finite_grid_certificate":
        issues.append(issue("matrix_rows[4]", "missing-grid-certificate", repr(rows[4].get("readiness"))))
    if rows and rows[6].get("role") != "rejected_route":
        issues.append(issue("matrix_rows[6]", "missing-promotion-rejection", repr(rows[6].get("role"))))

    diagnostics = artifact.get("diagnostics", {})
    params = diagnostics.get("parameters", {})
    expected_params = {
        "T_grid": list(DEFAULT_T_GRID),
        "indices": list(DEFAULT_INDICES),
        "grid_row_count": 20,
        "phi_term_count": 30,
        "ratio_cutoff_n": 80,
        "polynomial_M": 20,
        "taylor_degree": 384,
        "complex_disk_radius": "0.38",
        "core_x_cap_squared": "64/625",
        "core_x_cap": "8/25",
        "maximum_split_y": 200,
        "precision_bits": 4096,
    }
    for key, expected in expected_params.items():
        if params.get(key) != expected:
            issues.append(issue("parameters", f"bad-{key}", f"{params.get(key)!r} != {expected!r}"))

    for key in (
        "disk_radius_admissible",
        "global_n_tail_extension_certified",
        "all_value_expectations_certified_negative",
        "all_derivative_expectations_certified_negative",
        "all_value_ratios_below_one",
        "all_derivative_ratios_below_one",
        "all_denominator_crosschecks_pass",
        "all_floating_diagnostics_agree_within_tolerance",
        "complete_recorded_grid_expectations_certified",
    ):
        if diagnostics.get(key) is not True:
            issues.append(issue("diagnostics", f"false-{key}", repr(diagnostics.get(key))))
    if diagnostics.get("quadrature_needed_for_recorded_grid_expectations") is not False:
        issues.append(
            issue(
                "diagnostics",
                "quadrature-still-required",
                repr(diagnostics.get("quadrature_needed_for_recorded_grid_expectations")),
            )
        )
    if diagnostics.get("remaining_recorded_grid_integral_sources") != []:
        issues.append(
            issue(
                "diagnostics",
                "open-grid-integral-sources",
                repr(diagnostics.get("remaining_recorded_grid_integral_sources")),
            )
        )
    if diagnostics.get("target_closing") is not False:
        issues.append(issue("diagnostics", "target-closing", repr(diagnostics.get("target_closing"))))

    if diagnostics.get("value_worst_ratio_location") != {"T": 10000, "index": 21}:
        issues.append(
            issue(
                "diagnostics",
                "bad-value-worst-location",
                repr(diagnostics.get("value_worst_ratio_location")),
            )
        )
    if diagnostics.get("derivative_worst_ratio_location") != {"T": 10000, "index": 21}:
        issues.append(
            issue(
                "diagnostics",
                "bad-derivative-worst-location",
                repr(diagnostics.get("derivative_worst_ratio_location")),
            )
        )

    grid_rows = diagnostics.get("grid_rows", [])
    expected_keys = [(T, index) for T in DEFAULT_T_GRID for index in DEFAULT_INDICES]
    actual_keys = [(row.get("T"), row.get("index")) for row in grid_rows]
    if actual_keys != expected_keys:
        issues.append(issue("grid_rows", "bad-grid-keys", repr(actual_keys)))
    for position, row in enumerate(grid_rows):
        location = f"grid_rows[{position}]"
        for key in (
            "full_value_certified_negative",
            "full_derivative_certified_negative",
            "source_denominator_lowers_contained_by_computed_balls",
            "floating_diagnostic_agrees_within_tolerance",
            "row_direct_expectation_certified",
        ):
            if row.get(key) is not True:
                issues.append(issue(location, f"false-{key}", repr(row.get(key))))
        if row.get("quadrature_used_as_proof_input") is not False:
            issues.append(issue(location, "quadrature-proof-input", repr(row.get("quadrature_used_as_proof_input"))))
        try:
            value_ball = parse_arb(row["full_value_expectation_ball"])
            derivative_ball = parse_arb(row["full_derivative_expectation_ball"])
            value_ratio = parse_arb(row["value_full_ratio_to_first_omitted_upper"])
            derivative_ratio = parse_arb(row["derivative_full_ratio_to_first_omitted_upper"])
            value_margin = parse_arb(row["value_remaining_margin_below_one_lower"])
            derivative_margin = parse_arb(row["derivative_remaining_margin_below_one_lower"])
            q = parse_arb(row["x_to_disk_radius_ratio_upper"])
            first_value = parse_arb(row["value_first_omitted_expectation_lower"])
            first_derivative = parse_arb(row["derivative_first_omitted_expectation_lower"])
            value_radius = parse_arb(row["value_total_added_radius_upper"])
            derivative_radius = parse_arb(row["derivative_total_added_radius_upper"])
        except (KeyError, ValueError) as exc:
            issues.append(issue(location, "unparseable-row", str(exc)))
            continue
        checks = (
            (bool(value_ball < 0), "value-ball-not-negative", value_ball),
            (bool(derivative_ball < 0), "derivative-ball-not-negative", derivative_ball),
            (bool(0 < value_ratio < 1), "value-ratio-not-between-zero-one", value_ratio),
            (bool(0 < derivative_ratio < 1), "derivative-ratio-not-between-zero-one", derivative_ratio),
            (bool(value_margin > 0), "value-margin-not-positive", value_margin),
            (bool(derivative_margin > 0), "derivative-margin-not-positive", derivative_margin),
            (bool(0 < q < 1), "disk-ratio-invalid", q),
            (bool(first_value > 0), "value-first-omitted-not-positive", first_value),
            (bool(first_derivative > 0), "derivative-first-omitted-not-positive", first_derivative),
            (bool(value_radius > 0), "value-radius-not-positive", value_radius),
            (bool(derivative_radius > 0), "derivative-radius-not-positive", derivative_radius),
        )
        for passed, code, value in checks:
            if not passed:
                issues.append(issue(location, code, str(value)))
        if row.get("selected_compact_moment_row_count") != 13:
            issues.append(
                issue(
                    location,
                    "bad-selected-moment-count",
                    repr(row.get("selected_compact_moment_row_count")),
                )
            )
        proof_boundary = str(row.get("proof_boundary", "")).lower()
        if "no uniform collar claim" not in proof_boundary:
            issues.append(issue(location, "weak-row-boundary", proof_boundary))

    summary = artifact.get("summary", {})
    expected_summary_counts = {
        "matrix_rows": 8,
        "grid_rows": 20,
        "certified_value_rows": 20,
        "certified_derivative_rows": 20,
        "below_one_value_rows": 20,
        "below_one_derivative_rows": 20,
        "remaining_recorded_grid_integral_source_count": 0,
        "ready_to_apply_rows": 0,
    }
    for key, expected in expected_summary_counts.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected!r}"))
    if summary.get("complete_recorded_grid_expectations_certified") is not True:
        issues.append(issue("summary", "incomplete-grid", repr(summary.get("complete_recorded_grid_expectations_certified"))))
    if summary.get("quadrature_needed_for_recorded_grid_expectations") is not False:
        issues.append(issue("summary", "quadrature-needed", repr(summary.get("quadrature_needed_for_recorded_grid_expectations"))))
    if summary.get("target_closing") is not False:
        issues.append(issue("summary", "target-closing", repr(summary.get("target_closing"))))

    referenced = []
    for key, value in artifact.items():
        if key.startswith("source_") and isinstance(value, str):
            referenced.append(value)
    repo_root = Path(__file__).resolve().parents[3]
    for relative in referenced:
        if not (repo_root / relative).exists():
            issues.append(issue("sources", "missing-source", relative))

    invariants = " ".join(str(value) for value in artifact.get("invariants", [])).lower()
    for required in (
        "40 certified absolute ratios",
        "not a proof input",
        "no interpolation",
        "signed-stencil aggregation",
        "lambda <= 0",
    ):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant", required))

    expected_note_phrases = (
        "complete recorded-grid direct expectation certificate",
        "y_* = min(200, (64/625) T)",
        "degree-384",
        "remaining recorded-grid integration sources: `0`",
        "not used to derive any enclosure",
        "does not prove RH or Lambda <= 0",
        result_line(artifact),
    )
    note_lower = note_text.lower()
    for phrase in expected_note_phrases:
        if phrase.lower() not in note_lower:
            issues.append(issue("note", "missing-note-phrase", phrase))

    try:
        reproduced = build_artifact(first_omitted_path, floating_grid_path)
    except Exception as exc:  # pragma: no cover - reported as a checker issue
        issues.append(issue("reproduction", "rebuild-failed", repr(exc)))
    else:
        if artifact.get("diagnostics") != reproduced.get("diagnostics"):
            issues.append(issue("reproduction", "diagnostics-mismatch", "static artifact differs from rebuild"))
        if artifact.get("summary") != reproduced.get("summary"):
            issues.append(issue("reproduction", "summary-mismatch", "static artifact differs from rebuild"))
        if artifact.get("matrix_rows") != reproduced.get("matrix_rows"):
            issues.append(issue("reproduction", "matrix-mismatch", "static artifact differs from rebuild"))
    return issues


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--first-omitted-json", type=Path, default=DEFAULT_FIRST_OMITTED_JSON)
    parser.add_argument("--floating-grid-json", type=Path, default=DEFAULT_FLOATING_GRID_JSON)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    artifact = load_json(args.json)
    note_text = args.note.read_text(encoding="utf-8")
    issues = validate(artifact, note_text, args.first_omitted_json, args.floating_grid_json)
    if issues:
        print(json.dumps({"issues": issues}, indent=2))
        raise SystemExit(1)
    print(result_line(artifact))


if __name__ == "__main__":
    main()
