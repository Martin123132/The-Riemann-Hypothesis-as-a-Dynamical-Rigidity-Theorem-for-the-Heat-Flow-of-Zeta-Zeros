#!/usr/bin/env python3
"""Validate the recorded-grid fixed-k stencil composition certificate."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate import (  # noqa: E402
    DEFAULT_EXPECTATION_JSON,
    DEFAULT_INDICES,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_RESIDUAL_BUDGET_JSON,
    DEFAULT_T_GRID,
    build_artifact,
    result_line,
)


EXPECTED_KIND = (
    "jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate"
)
EXPECTED_IDS = [
    "nlrgrgscc_01_residual_expectation_identity",
    "nlrgrgscc_02_rational_residual_budgets",
    "nlrgrgscc_03_arb_perturbation_ledger",
    "nlrgrgscc_04_all_row_budget_composition",
    "nlrgrgscc_05_recorded_T_stencil_certificate",
    "nlrgrgscc_06_real_T_promotion_rejected",
    "nlrgrgscc_07_all_k_promotion_rejected",
    "nlrgrgscc_08_acceptance_gate",
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
    expectation_path: Path,
    residual_budget_path: Path,
) -> list[dict[str, str]]:
    issues: list[dict[str, str]] = []
    if artifact.get("kind") != EXPECTED_KIND:
        issues.append(issue("artifact", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("date") != "2026-07-10":
        issues.append(issue("artifact", "bad-date", repr(artifact.get("date"))))
    if artifact.get("status") != "recorded-grid fixed-k stencil composition certificate":
        issues.append(issue("artifact", "bad-status", repr(artifact.get("status"))))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "fixed-k=22",
        "does not certify an interval in t",
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
    if rows and rows[2].get("readiness") != "available_interval_certificate":
        issues.append(issue("matrix_rows[2]", "missing-perturbation-certificate", repr(rows[2].get("readiness"))))
    if rows and rows[4].get("readiness") != "available_finite_grid_certificate":
        issues.append(issue("matrix_rows[4]", "missing-stencil-certificate", repr(rows[4].get("readiness"))))
    if rows and (rows[5].get("role"), rows[6].get("role")) != (
        "rejected_route",
        "rejected_route",
    ):
        issues.append(issue("matrix_rows", "missing-rejected-promotions", repr(rows[5:7])))

    diagnostics = artifact.get("diagnostics", {})
    params = diagnostics.get("parameters", {})
    expected_params = {
        "T_grid": list(DEFAULT_T_GRID),
        "indices": list(DEFAULT_INDICES),
        "k": 22,
        "finite_degree": 40,
        "polynomial_M": 20,
        "ratio_cutoff_n": 80,
        "collar_start_T": 1156,
        "precision_bits": 4096,
        "value_budget_A": "1/2",
        "derivative_budget_B": "9/1000",
        "value_residual_identity": "expectation_value=R_i(u)",
        "derivative_residual_identity": "expectation_derivative=u*R_i'(u)",
    }
    for key, expected in expected_params.items():
        if params.get(key) != expected:
            issues.append(issue("parameters", f"bad-{key}", f"{params.get(key)!r} != {expected!r}"))
    if diagnostics.get("source_float_budgets_used_as_proof_inputs") is not False:
        issues.append(
            issue(
                "diagnostics",
                "float-budget-proof-input",
                repr(diagnostics.get("source_float_budgets_used_as_proof_inputs")),
            )
        )
    for key in (
        "all_retained_margins_positive",
        "all_20_value_residual_budgets_certified",
        "all_20_derivative_residual_budgets_certified",
        "all_5_recorded_T_stencil_systems_certified",
    ):
        if diagnostics.get(key) is not True:
            issues.append(issue("diagnostics", f"false-{key}", repr(diagnostics.get(key))))
    if diagnostics.get("remaining_recorded_grid_stencil_sources") != []:
        issues.append(
            issue(
                "diagnostics",
                "open-grid-stencil-sources",
                repr(diagnostics.get("remaining_recorded_grid_stencil_sources")),
            )
        )
    if diagnostics.get("target_closing") is not False:
        issues.append(issue("diagnostics", "target-closing", repr(diagnostics.get("target_closing"))))

    perturbation_rows = diagnostics.get("perturbation_rows", [])
    if [row.get("name") for row in perturbation_rows] != [
        "normalizer",
        "B_product",
        "companion_product",
        "weighted_gap_derivative",
    ]:
        issues.append(
            issue(
                "perturbation_rows",
                "bad-targets",
                repr([row.get("name") for row in perturbation_rows]),
            )
        )
    for position, row in enumerate(perturbation_rows):
        location = f"perturbation_rows[{position}]"
        if row.get("certified_positive_after_perturbation") is not True:
            issues.append(issue(location, "uncertified-perturbation", repr(row)))
        try:
            finite_margin = parse_arb(row["finite_margin_lower"])
            perturbation = parse_arb(row["perturbation_bound_upper"])
            retained = parse_arb(row["retained_margin_lower"])
        except (KeyError, ValueError) as exc:
            issues.append(issue(location, "unparseable", str(exc)))
            continue
        if not bool(finite_margin > 0 and perturbation > 0 and retained > 0):
            issues.append(issue(location, "nonpositive-ledger-entry", repr(row)))
        if not bool(retained < finite_margin and perturbation < finite_margin):
            issues.append(issue(location, "invalid-retained-margin-order", repr(row)))

    residual_rows = diagnostics.get("grid_residual_rows", [])
    expected_keys = [(T, index) for T in DEFAULT_T_GRID for index in DEFAULT_INDICES]
    actual_keys = [(row.get("T"), row.get("index")) for row in residual_rows]
    if actual_keys != expected_keys:
        issues.append(issue("grid_residual_rows", "bad-grid", repr(actual_keys)))
    A = flint.arb(1) / 2
    B = flint.arb(9) / 1000
    for position, row in enumerate(residual_rows):
        location = f"grid_residual_rows[{position}]"
        for key in (
            "value_budget_certified",
            "derivative_budget_certified",
            "row_residual_budget_certified",
        ):
            if row.get(key) is not True:
                issues.append(issue(location, f"false-{key}", repr(row.get(key))))
        try:
            value_scaled = parse_arb(row["value_residual_scaled_upper"])
            derivative_scaled = parse_arb(row["derivative_residual_scaled_upper"])
            value_fraction = parse_arb(row["value_fraction_of_rational_budget_upper"])
            derivative_fraction = parse_arb(row["derivative_fraction_of_rational_budget_upper"])
        except (KeyError, ValueError) as exc:
            issues.append(issue(location, "unparseable", str(exc)))
            continue
        if not bool(0 < value_scaled < A):
            issues.append(issue(location, "value-budget-failed", str(value_scaled)))
        if not bool(0 < derivative_scaled < B):
            issues.append(issue(location, "derivative-budget-failed", str(derivative_scaled)))
        if not bool(0 < value_fraction < 1 and 0 < derivative_fraction < 1):
            issues.append(issue(location, "fraction-not-below-one", repr(row)))
        if row.get("identity_used") != "value expectation=R_i(u); derivative expectation=u*R_i'(u)":
            issues.append(issue(location, "bad-identity", repr(row.get("identity_used"))))

    t_rows = diagnostics.get("recorded_T_stencil_rows", [])
    if [row.get("T") for row in t_rows] != list(DEFAULT_T_GRID):
        issues.append(issue("recorded_T_stencil_rows", "bad-T-grid", repr(t_rows)))
    for position, row in enumerate(t_rows):
        location = f"recorded_T_stencil_rows[{position}]"
        for key in (
            "all_four_normalizers_certified_positive",
            "B_product_certified_positive",
            "companion_product_certified_positive",
            "weighted_gap_derivative_certified_positive",
            "recorded_T_stencil_system_certified",
        ):
            if row.get(key) is not True:
                issues.append(issue(location, f"false-{key}", repr(row.get(key))))
        if row.get("certified_indices") != list(DEFAULT_INDICES):
            issues.append(issue(location, "bad-indices", repr(row.get("certified_indices"))))
        boundary = str(row.get("proof_boundary", "")).lower()
        if "not an interval in t" not in boundary or "all-k" not in boundary:
            issues.append(issue(location, "weak-proof-boundary", boundary))

    summary = artifact.get("summary", {})
    expected_counts = {
        "matrix_rows": 8,
        "grid_residual_rows": 20,
        "recorded_T_stencil_rows": 5,
        "positive_retained_perturbation_margins": 4,
        "certified_recorded_T_stencil_systems": 5,
        "remaining_recorded_grid_stencil_source_count": 0,
        "ready_to_apply_rows": 0,
    }
    for key, expected in expected_counts.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected!r}"))
    if summary.get("all_5_recorded_T_stencil_systems_certified") is not True:
        issues.append(issue("summary", "incomplete-stencil-grid", repr(summary)))
    if summary.get("target_closing") is not False:
        issues.append(issue("summary", "target-closing", repr(summary.get("target_closing"))))
    for key in (
        "maximum_value_fraction_of_rational_budget_upper",
        "maximum_derivative_fraction_of_rational_budget_upper",
    ):
        try:
            fraction = parse_arb(summary[key])
        except (KeyError, ValueError) as exc:
            issues.append(issue("summary", f"unparseable-{key}", str(exc)))
        else:
            if not bool(0 < fraction < flint.arb("0.001")):
                issues.append(issue("summary", f"large-{key}", str(fraction)))

    repo_root = Path(__file__).resolve().parents[3]
    for key, value in artifact.items():
        if key.startswith("source_") and isinstance(value, str):
            if not (repo_root / value).exists():
                issues.append(issue("sources", "missing-source", value))

    invariant_text = " ".join(str(value) for value in artifact.get("invariants", [])).lower()
    for required in (
        "not proof inputs",
        "exact rationals",
        "arb ball arithmetic",
        "not promoted to the real-t collar",
        "not promoted to all k",
        "lambda <= 0",
    ):
        if required not in invariant_text:
            issues.append(issue("invariants", "missing-invariant", required))

    for required in (
        "recorded-grid fixed-k stencil composition certificate",
        "|R_i(u)|  <= (1/2) u^3",
        "not proof inputs",
        "positive perturbation margins",
        "does not interpolate between them",
        "does not prove RH or Lambda <= 0",
        result_line(artifact),
    ):
        if required.lower() not in note_text.lower():
            issues.append(issue("note", "missing-note-phrase", required))

    try:
        reproduced = build_artifact(expectation_path, residual_budget_path)
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
    parser.add_argument("--expectation-json", type=Path, default=DEFAULT_EXPECTATION_JSON)
    parser.add_argument("--residual-budget-json", type=Path, default=DEFAULT_RESIDUAL_BUDGET_JSON)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    artifact = load_json(args.json)
    note_text = args.note.read_text(encoding="utf-8")
    issues = validate(artifact, note_text, args.expectation_json, args.residual_budget_json)
    if issues:
        print(json.dumps({"issues": issues}, indent=2))
        raise SystemExit(1)
    print(result_line(artifact))


if __name__ == "__main__":
    main()
