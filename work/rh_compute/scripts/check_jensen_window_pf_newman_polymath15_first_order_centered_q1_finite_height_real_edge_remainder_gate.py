#!/usr/bin/env python3
"""Independently validate the q=1 finite-height real-edge remainder gate."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import jensen_window_pf_newman_polymath15_first_order_centered_q1_finite_height_real_edge_remainder_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "q1_finite_height_real_edge_remainder_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    "q1fher_01_domain",
    "q1fher_02_endpoint_source",
    "q1fher_03_terminal_source",
    "q1fher_04_stable_log",
    "q1fher_05_removable_rewrites",
    "q1fher_06_terminal_bounds",
    "q1fher_07_mu_bounds",
    "q1fher_08_c0_cauchy",
    "q1fher_09_c0_direct_cover",
    "q1fher_10_endpoint_jets",
    "q1fher_11_value_jet",
    "q1fher_12_first_jets",
    "q1fher_13_second_jet",
    "q1fher_14_current_budget",
    "q1fher_15_finite_sign",
    "q1fher_16_handoff",
]
EXPECTED_COUNTS = {
    "rows": 16,
    "direct_arb_boxes": 4096,
    "removable_cauchy_charts": 2,
    "source_derivative_bounds": 6,
    "stable_terminal_log_identities": 1,
    "terminal_value_derivative_bounds": 2,
    "endpoint_mu_bounds": 2,
    "finite_edge_jet_defects": 4,
    "uniform_q1_finite_height_edge_signs": 1,
    "uniform_q_gt_1_edge_signs": 0,
    "xi_level_edge_signs": 0,
    "abel_gaps": 0,
    "winding_bounds": 0,
}
SUMMARY = (
    "validated q=1 finite-height real-edge remainder gate: 16 rows, "
    "4096 Arb boxes, 2 Cauchy charts, 4 finite edge jets, "
    "remainder <1/10000000, 1 retained-model finite-height sign, "
    "0 q>1 signs, 0 Xi-level signs"
)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing result: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid result: {exc}")
        return {}


def validate_payload(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or order drifted")

    status = payload.get("status", "")
    for token in ("finite-height", "q=1", "first-order model"):
        if token not in status:
            issues.append(f"status missing token: {token}")

    boundary = payload.get("proof_boundary", "")
    for token in (
        "omitted higher-order Xi approximation error",
        "q>1 edge sign",
        "adjacent cutoffs",
        "Abel-scalar gap",
        "successor winding cap",
        "Lambda<=0",
        "RH",
    ):
        if token not in boundary:
            issues.append(f"proof boundary missing token: {token}")

    audit = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif audit.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    try:
        gate.verify_symbolics()
    except Exception as exc:  # noqa: BLE001
        issues.append(f"symbolic audit failed: {exc}")

    derivatives = payload.get("derivative_certificate", {})
    try:
        if derivatives != gate.c0_derivative_certificate():
            issues.append("C_0 derivative certificate drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"C_0 derivative recomputation failed: {exc}")
    if derivatives.get("direct_boxes") != 4096:
        issues.append("direct C_0 interval count drifted")
    if derivatives.get("global_C0_derivative_bounds") != list(
        gate.F_DERIVATIVE_BOUNDS
    ):
        issues.append("global C_0 derivative bounds drifted")
    if derivatives.get("C1_derivative_bounds") != list(
        gate.C1_DERIVATIVE_BOUNDS
    ):
        issues.append("C_1 derivative bounds drifted")

    budget = payload.get("majorant_certificate", {})
    try:
        if budget != gate.majorant_certificate():
            issues.append("normalized majorant certificate drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"majorant recomputation failed: {exc}")
    chain = budget.get("normalized_majorant_chain", {})
    if len(chain) != 30:
        issues.append("normalized majorant chain length drifted")
    for name, row in chain.items():
        try:
            if not Fraction(row["upper"]) < Fraction(row["cap"]):
                issues.append(f"majorant cap failed: {name}")
        except (KeyError, ValueError, ZeroDivisionError) as exc:
            issues.append(f"invalid majorant row {name}: {exc}")

    if budget.get("linear_error_constant") != 5400:
        issues.append("linear current constant drifted")
    if budget.get("quadratic_error_constant") != 152500:
        issues.append("quadratic current constant drifted")
    if budget.get("uniform_remainder") != (
        "|a^2 J_edge/S_a^2-K_edge(p)|<5500h<1/10000000"
    ):
        issues.append("uniform finite-height remainder drifted")
    if budget.get("finite_height_margin") != (
        "a^2 J_edge/S_a^2<-3749/10000000<0"
    ):
        issues.append("finite-height sign margin drifted")

    exact = payload.get("exact", {})
    stable = exact.get("stable_log", {})
    for key in (
        "coordinates",
        "real_part",
        "imaginary_part",
        "removable_rewrites",
        "x_derivative",
        "branch",
    ):
        if not stable.get(key):
            issues.append(f"stable logarithm field missing: {key}")


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Effective Domain",
        "Stable Terminal Logarithm",
        "Endpoint Derivatives",
        "Four-Jet Budget",
        "Finite-Height Sign",
        "4096",
        "1/10000000",
        "-3749/10000000",
        "omitted higher-order Xi approximation error",
        "does not extend the sign to q>1",
        "No winding cap, contact exclusion, Lambda<=0",
    ):
        if token not in text:
            issues.append(f"note marker missing: {token}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    validate_payload(payload, issues)
    validate_note(issues)
    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        return 1
    print(SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
