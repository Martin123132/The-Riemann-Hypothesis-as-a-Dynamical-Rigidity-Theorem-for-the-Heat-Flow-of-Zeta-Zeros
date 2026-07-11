#!/usr/bin/env python3
"""Validate the exact full-kernel evenness and Cauchy lemma."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    build_artifact,
    result_line,
)


EXPECTED_KIND = (
    "jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma"
)
EXPECTED_IDS = [
    "nlrgfkec_01_jacobi_theta_input",
    "nlrgfkec_02_kernel_operator_identity",
    "nlrgfkec_03_operator_covariance",
    "nlrgfkec_04_full_kernel_evenness",
    "nlrgfkec_05_disk_analyticity",
    "nlrgfkec_06_order42_residual_zero",
    "nlrgfkec_07_factored_cauchy_bounds",
    "nlrgfkec_08_finite_truncation_promotion_rejected",
    "nlrgfkec_09_acceptance_gate",
]


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def issue(location: str, code: str, detail: str) -> dict[str, str]:
    return {"location": location, "code": code, "detail": detail}


def validate(artifact: dict, note_text: str) -> list[dict[str, str]]:
    issues: list[dict[str, str]] = []
    if artifact.get("kind") != EXPECTED_KIND:
        issues.append(issue("artifact", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("date") != "2026-07-10":
        issues.append(issue("artifact", "bad-date", repr(artifact.get("date"))))
    if artifact.get("status") != "exact full-kernel evenness and Cauchy lemma":
        issues.append(issue("artifact", "bad-status", repr(artifact.get("status"))))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "exact full-kernel parity",
        "does not bound a gamma expectation",
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
    if rows and rows[3].get("role") != "exact_analytic_lemma":
        issues.append(issue("matrix_rows[3]", "missing-evenness-lemma", repr(rows[3].get("role"))))
    if rows and rows[5].get("role") != "exact_analytic_lemma":
        issues.append(issue("matrix_rows[5]", "missing-order42-lemma", repr(rows[5].get("role"))))
    if rows and rows[7].get("role") != "rejected_route":
        issues.append(issue("matrix_rows[7]", "missing-finite-rejection", repr(rows[7].get("role"))))

    diagnostics = artifact.get("diagnostics", {})
    symbolic = diagnostics.get("symbolic", {})
    for key in (
        "operator_covariance_verified",
        "annihilator_x_minus_half_verified",
        "annihilator_constant_verified",
    ):
        if symbolic.get(key) is not True:
            issues.append(issue("symbolic", f"false-{key}", repr(symbolic.get(key))))
    for key in (
        "all_symbolic_operator_identities_verified",
        "full_kernel_evenness_certified",
        "all_odd_taylor_coefficients_zero",
        "residual_order_42_certified",
    ):
        if diagnostics.get(key) is not True:
            issues.append(issue("diagnostics", f"false-{key}", repr(diagnostics.get(key))))
    if diagnostics.get("full_kernel_evenness_identity") != "Phi(-u)=Phi(u)":
        issues.append(issue("diagnostics", "bad-evenness-identity", repr(diagnostics.get("full_kernel_evenness_identity"))))
    if diagnostics.get("polynomial_M") != 20:
        issues.append(issue("diagnostics", "bad-M", repr(diagnostics.get("polynomial_M"))))
    if diagnostics.get("subtracted_even_degree") != 40:
        issues.append(issue("diagnostics", "bad-subtracted-degree", repr(diagnostics.get("subtracted_even_degree"))))
    if diagnostics.get("first_residual_degree") != 42:
        issues.append(issue("diagnostics", "bad-first-degree", repr(diagnostics.get("first_residual_degree"))))
    if diagnostics.get("residual_zero_order_at_least") != 42:
        issues.append(issue("diagnostics", "bad-zero-order", repr(diagnostics.get("residual_zero_order_at_least"))))
    if diagnostics.get("target_closing") is not False:
        issues.append(issue("diagnostics", "target-closing", repr(diagnostics.get("target_closing"))))
    analytic = diagnostics.get("analytic", {})
    if analytic.get("disk_radius") != "0.38":
        issues.append(issue("analytic", "bad-radius", repr(analytic.get("disk_radius"))))
    if analytic.get("disk_admissible") is not True:
        issues.append(issue("analytic", "inadmissible-disk", repr(analytic)))
    for formula_key in ("cauchy_value_factor_formula", "cauchy_derivative_factor_formula"):
        formula = str(diagnostics.get(formula_key, ""))
        if "x^42" not in formula or "R^42" not in formula:
            issues.append(issue("diagnostics", f"bad-{formula_key}", formula))

    summary = artifact.get("summary", {})
    expected_summary = {
        "matrix_rows": 9,
        "symbolic_operator_identities": 3,
        "symbolic_operator_identities_verified": 3,
        "disk_radius": "0.38",
        "full_kernel_evenness_certified": True,
        "residual_zero_order_at_least": 42,
        "residual_order_42_certified": True,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected!r}"))

    repo_root = Path(__file__).resolve().parents[3]
    for key, value in artifact.items():
        if key.startswith("source_") and isinstance(value, str):
            if not (repo_root / value).exists():
                issues.append(issue("sources", "missing-source", value))

    invariant_text = " ".join(str(value) for value in artifact.get("invariants", [])).lower()
    for required in (
        "full theta kernel",
        "poisson-summation",
        "order-42 zero",
        "not formal full-gamma",
        "lambda <= 0",
    ):
        if required not in invariant_text:
            issues.append(issue("invariants", "missing-invariant", required))

    for required in (
        "exact full-kernel evenness and Cauchy lemma",
        "Phi(-u)=Phi(u)",
        "zero of order at least `42`",
        "it is not inferred",
        "Cauchy Factors",
        "does not prove RH or Lambda <= 0",
        result_line(artifact),
    ):
        if required.lower() not in note_text.lower():
            issues.append(issue("note", "missing-note-phrase", required))

    try:
        reproduced = build_artifact()
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
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    artifact = load_json(args.json)
    note_text = args.note.read_text(encoding="utf-8")
    issues = validate(artifact, note_text)
    if issues:
        print(json.dumps({"issues": issues}, indent=2))
        raise SystemExit(1)
    print(result_line(artifact))


if __name__ == "__main__":
    main()
