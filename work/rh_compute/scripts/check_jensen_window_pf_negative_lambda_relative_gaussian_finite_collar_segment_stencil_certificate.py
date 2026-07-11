#!/usr/bin/env python3
"""Validate the finite-collar-segment fixed-k stencil certificate."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate import (  # noqa: E402
    DEFAULT_INDICES,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_RECORDED_GRID_JSON,
    build_artifact,
    result_line,
)


EXPECTED_KIND = (
    "jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate"
)
EXPECTED_IDS = [
    "nlrgfcssc_01_two_regime_split",
    "nlrgfcssc_02_compact_polynomial_majorant",
    "nlrgfcssc_03_uniform_cauchy_remainder",
    "nlrgfcssc_04_two_regime_real_tail",
    "nlrgfcssc_05_global_normalization_tail",
    "nlrgfcssc_06_uniform_residual_budget_certificate",
    "nlrgfcssc_07_finite_collar_segment_stencil_certificate",
    "nlrgfcssc_08_T_gt_10000_promotion_rejected",
    "nlrgfcssc_09_all_k_promotion_rejected",
    "nlrgfcssc_10_acceptance_gate",
]


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_arb(text: object) -> flint.arb:
    return flint.arb(str(text).replace("E", "e"))


def issue(location: str, code: str, detail: str) -> dict[str, str]:
    return {"location": location, "code": code, "detail": detail}


def validate(artifact: dict, note_text: str, recorded_grid_path: Path) -> list[dict[str, str]]:
    issues: list[dict[str, str]] = []
    if artifact.get("kind") != EXPECTED_KIND:
        issues.append(issue("artifact", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("date") != "2026-07-10":
        issues.append(issue("artifact", "bad-date", repr(artifact.get("date"))))
    if artifact.get("status") != "bounded real-T fixed-k stencil certificate":
        issues.append(issue("artifact", "bad-status", repr(artifact.get("status"))))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "1156<=t<=10000",
        "does not cover t>10000",
        "does not cover all k",
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
        issues.append(issue("matrix_rows[5]", "missing-residual-certificate", repr(rows[5].get("role"))))
    if rows and rows[6].get("role") != "real_interval_certificate":
        issues.append(issue("matrix_rows[6]", "missing-stencil-certificate", repr(rows[6].get("role"))))
    if rows and (rows[7].get("role"), rows[8].get("role")) != (
        "rejected_route",
        "rejected_route",
    ):
        issues.append(issue("matrix_rows", "missing-promotion-rejections", repr(rows[7:9])))

    diagnostics = artifact.get("diagnostics", {})
    params = diagnostics.get("parameters", {})
    expected_params = {
        "T_interval": "1156<=T<=10000",
        "T_left": 1156,
        "T_transition": "15625/8",
        "T_right": 10000,
        "indices": list(DEFAULT_INDICES),
        "polynomial_M": 20,
        "taylor_degree": 384,
        "phi_term_count": 30,
        "ratio_cutoff_n": 80,
        "disk_radius": "0.38",
        "x_cap": "8/25",
        "fixed_split_y": 200,
        "precision_bits": 4096,
        "value_budget_A": "0.5",
        "derivative_budget_B": "0.009",
    }
    for key, expected in expected_params.items():
        if params.get(key) != expected:
            issues.append(issue("parameters", f"bad-{key}", f"{params.get(key)!r} != {expected!r}"))
    split_rule = str(diagnostics.get("split_rule", ""))
    for required in ("1156<=T<=15625/8", "x_*=8/25", "15625/8<=T<=10000", "y_*=200"):
        if required not in split_rule:
            issues.append(issue("diagnostics", "bad-split-rule", required))
    for key in (
        "disk_radius_admissible",
        "global_n_tail_extension_certified",
        "all_four_value_budgets_certified_uniformly",
        "all_four_derivative_budgets_certified_uniformly",
        "all_four_residual_budgets_certified_uniformly",
        "all_retained_perturbation_margins_positive",
        "finite_collar_segment_stencil_system_certified",
    ):
        if diagnostics.get(key) is not True:
            issues.append(issue("diagnostics", f"false-{key}", repr(diagnostics.get(key))))
    if diagnostics.get("remaining_finite_segment_stencil_sources") != []:
        issues.append(
            issue(
                "diagnostics",
                "open-segment-sources",
                repr(diagnostics.get("remaining_finite_segment_stencil_sources")),
            )
        )
    if diagnostics.get("target_closing") is not False:
        issues.append(issue("diagnostics", "target-closing", repr(diagnostics.get("target_closing"))))

    for key in (
        "cauchy_value_scaled_uniform_upper",
        "cauchy_derivative_scaled_uniform_upper",
        "normalization_value_scaled_uniform_upper",
        "normalization_derivative_scaled_uniform_upper",
    ):
        try:
            value = parse_arb(diagnostics[key])
        except (KeyError, ValueError) as exc:
            issues.append(issue("diagnostics", f"unparseable-{key}", str(exc)))
        else:
            if not bool(value > 0):
                issues.append(issue("diagnostics", f"nonpositive-{key}", str(value)))

    uniform_rows = diagnostics.get("uniform_residual_rows", [])
    if [row.get("index") for row in uniform_rows] != list(DEFAULT_INDICES):
        issues.append(issue("uniform_residual_rows", "bad-indices", repr(uniform_rows)))
    for position, row in enumerate(uniform_rows):
        location = f"uniform_residual_rows[{position}]"
        for key in (
            "value_budget_certified_uniformly",
            "derivative_budget_certified_uniformly",
            "uniform_residual_budget_certified",
        ):
            if row.get(key) is not True:
                issues.append(issue(location, f"false-{key}", repr(row.get(key))))
        if row.get("T_interval") != "1156<=T<=10000":
            issues.append(issue(location, "bad-T-interval", repr(row.get("T_interval"))))
        if row.get("selected_compact_majorant_row_count") != 13:
            issues.append(issue(location, "bad-selected-count", repr(row.get("selected_compact_majorant_row_count"))))
        try:
            value_scaled = parse_arb(row["value_residual_scaled_uniform_upper"])
            derivative_scaled = parse_arb(row["derivative_residual_scaled_uniform_upper"])
            value_fraction = parse_arb(row["value_fraction_of_budget_upper"])
            derivative_fraction = parse_arb(row["derivative_fraction_of_budget_upper"])
        except (KeyError, ValueError) as exc:
            issues.append(issue(location, "unparseable", str(exc)))
            continue
        if not bool(0 < value_scaled < flint.arb("0.5")):
            issues.append(issue(location, "value-budget-failed", str(value_scaled)))
        if not bool(0 < derivative_scaled < flint.arb("0.009")):
            issues.append(issue(location, "derivative-budget-failed", str(derivative_scaled)))
        if not bool(0 < value_fraction < flint.arb("0.002")):
            issues.append(issue(location, "value-fraction-too-large", str(value_fraction)))
        if not bool(0 < derivative_fraction < flint.arb("0.002")):
            issues.append(issue(location, "derivative-fraction-too-large", str(derivative_fraction)))
        boundary = str(row.get("proof_boundary", "")).lower()
        if "not t>10000" not in boundary:
            issues.append(issue(location, "weak-proof-boundary", boundary))

    perturbation_rows = diagnostics.get("perturbation_rows", [])
    if [row.get("name") for row in perturbation_rows] != [
        "normalizer",
        "B_product",
        "companion_product",
        "weighted_gap_derivative",
    ]:
        issues.append(issue("perturbation_rows", "bad-targets", repr(perturbation_rows)))
    for position, row in enumerate(perturbation_rows):
        location = f"perturbation_rows[{position}]"
        if row.get("certified_positive_after_perturbation") is not True:
            issues.append(issue(location, "uncertified", repr(row)))
        try:
            finite_margin = parse_arb(row["finite_margin_lower"])
            perturbation = parse_arb(row["perturbation_bound_upper"])
            retained = parse_arb(row["retained_margin_lower"])
        except (KeyError, ValueError) as exc:
            issues.append(issue(location, "unparseable", str(exc)))
            continue
        if not bool(0 < perturbation < finite_margin and 0 < retained < finite_margin):
            issues.append(issue(location, "invalid-margin-order", repr(row)))

    summary = artifact.get("summary", {})
    expected_counts = {
        "matrix_rows": 10,
        "uniform_residual_rows": 4,
        "T_regimes": 2,
        "positive_retained_perturbation_margins": 4,
        "remaining_finite_segment_stencil_source_count": 0,
        "ready_to_apply_rows": 0,
    }
    for key, expected in expected_counts.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected!r}"))
    if summary.get("finite_collar_segment_stencil_system_certified") is not True:
        issues.append(issue("summary", "segment-not-certified", repr(summary)))
    if summary.get("target_closing") is not False:
        issues.append(issue("summary", "target-closing", repr(summary.get("target_closing"))))
    for key in (
        "maximum_value_fraction_of_budget_upper",
        "maximum_derivative_fraction_of_budget_upper",
    ):
        try:
            fraction = parse_arb(summary[key])
        except (KeyError, ValueError) as exc:
            issues.append(issue("summary", f"unparseable-{key}", str(exc)))
        else:
            if not bool(0 < fraction < flint.arb("0.002")):
                issues.append(issue("summary", f"large-{key}", str(fraction)))

    repo_root = Path(__file__).resolve().parents[3]
    for key, value in artifact.items():
        if key.startswith("source_") and isinstance(value, str):
            if not (repo_root / value).exists():
                issues.append(issue("sources", "missing-source", value))

    invariant_text = " ".join(str(value) for value in artifact.get("invariants", [])).lower()
    for required in (
        "every real t in [1156,10000]",
        "no floating quadrature",
        "does not use the obstructed formal",
        "not promoted to t>10000",
        "not promoted to all k",
        "lambda <= 0",
    ):
        if required not in invariant_text:
            issues.append(issue("invariants", "missing-invariant", required))

    for required in (
        "bounded real-T fixed-k stencil certificate",
        "1156<=T<=10000",
        "Both regimes",
        "All four rows therefore satisfy",
        "next T obligation is the",
        "does not prove RH or Lambda <= 0",
        result_line(artifact),
    ):
        if required.lower() not in note_text.lower():
            issues.append(issue("note", "missing-note-phrase", required))

    try:
        reproduced = build_artifact(recorded_grid_path)
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
    parser.add_argument("--recorded-grid-json", type=Path, default=DEFAULT_RECORDED_GRID_JSON)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    artifact = load_json(args.json)
    note_text = args.note.read_text(encoding="utf-8")
    issues = validate(artifact, note_text, args.recorded_grid_json)
    if issues:
        print(json.dumps({"issues": issues}, indent=2))
        raise SystemExit(1)
    print(result_line(artifact))


if __name__ == "__main__":
    main()
