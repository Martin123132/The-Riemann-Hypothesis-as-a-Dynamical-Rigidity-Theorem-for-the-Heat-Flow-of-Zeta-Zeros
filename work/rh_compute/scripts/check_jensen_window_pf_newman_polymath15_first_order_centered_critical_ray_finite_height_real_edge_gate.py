#!/usr/bin/env python3
"""Validate the full critical-ray finite-height real-edge sign gate."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import jensen_window_pf_newman_polymath15_first_order_centered_critical_ray_finite_height_real_edge_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "critical_ray_finite_height_real_edge_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    "crfher_01_domain",
    "crfher_02_exact_source",
    "crfher_03_changed_caps",
    "crfher_04_stable_heat",
    "crfher_05_source_bounds",
    "crfher_06_source_derivatives",
    "crfher_07_terminal_jets",
    "crfher_08_four_jets",
    "crfher_09_current_budget",
    "crfher_10_uniform_sign",
    "crfher_11_q_ge_1",
    "crfher_12_nonpromotion",
    "crfher_13_handoff",
]
EXPECTED_COUNTS = {
    "rows": 13,
    "direct_arb_boxes": 4096,
    "removable_cauchy_charts": 2,
    "normalized_majorant_inequalities": 27,
    "finite_edge_jet_defects": 4,
    "uniform_critical_ray_model_signs": 1,
    "uniform_q_ge_1_model_signs": 1,
    "adjacent_cutoff_signed_splices": 0,
    "xi_level_edge_signs": 0,
    "abel_gaps": 0,
    "winding_bounds": 0,
}
SUMMARY = (
    "validated critical-ray finite-height real-edge gate: 13 rows, "
    "4096 Arb boxes, 27 normalized majorants, 4 finite edge jets, "
    "1 full critical-ray model sign, 1 q>=1 model sign, "
    "0 cutoff splices, 0 Xi-level signs"
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


def validate(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    if [row.get("id") for row in payload.get("rows", [])] != EXPECTED_IDS:
        issues.append("row ids or order drifted")

    for token in ("finite-height", "critical-ray", "first-order model"):
        if token not in payload.get("status", ""):
            issues.append(f"status missing token: {token}")
    for token in (
        "adjacent-cutoff signed splice",
        "higher-order Xi/source",
        "Xi-level edge sign",
        "Abel-scalar gap",
        "successor winding cap",
        "Lambda<=0",
        "RH",
    ):
        if token not in payload.get("proof_boundary", ""):
            issues.append(f"proof boundary missing token: {token}")

    audit = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif audit.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    try:
        gate.q1.verify_symbolics()
    except Exception as exc:  # noqa: BLE001
        issues.append(f"symbolic inheritance failed: {exc}")

    derivatives = payload.get("derivative_certificate", {})
    try:
        if derivatives != gate.q1.c0_derivative_certificate():
            issues.append("C_0 derivative certificate drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"C_0 derivative recomputation failed: {exc}")
    if derivatives.get("direct_boxes") != 4096:
        issues.append("C_0 interval count drifted")

    budget = payload.get("majorant_certificate", {})
    try:
        if budget != gate.majorant_certificate():
            issues.append("critical-ray majorant certificate drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"critical-ray majorant recomputation failed: {exc}")
    chain = budget.get("normalized_majorant_chain", {})
    if len(chain) != 27:
        issues.append("critical-ray majorant count drifted")
    for name, row in chain.items():
        try:
            if not Fraction(row["upper"]) < Fraction(row["cap"]):
                issues.append(f"majorant cap failed: {name}")
        except (KeyError, ValueError, ZeroDivisionError) as exc:
            issues.append(f"invalid majorant row {name}: {exc}")

    expected_source = {
        "rho": "rho<h^2/16",
        "log_r": "|log r|<h^2/16",
        "alpha_prime": "|alpha'|<h^2",
        "alpha_second": "|alpha''|<h^4",
        "correction": "|d|<h^2/4, |d_x|<h^4/32",
        "stable_log": "|W|<4h, |W_x|<h^2",
        "terminal": "|Z+Q|<8h, |Z_x-i*h*theta*Q/2|<6h^2",
        "endpoint_rate": "|mu|<h^2, |mu_x|<h^4",
    }
    if budget.get("source_bounds") != expected_source:
        issues.append("critical-ray source bounds drifted")
    if budget.get("uniform_remainder") != (
        "|a^2 J_edge/S_a^2-K_edge(p)|<5500h<1/10000000"
    ):
        issues.append("critical-ray remainder drifted")
    if budget.get("finite_height_margin") != (
        "a^2 J_edge/S_a^2<-3749/10000000<0"
    ):
        issues.append("critical-ray finite sign drifted")

    domain = payload.get("exact", {}).get("domain", {})
    for token in ("0<tL<=25", "0<t<=1/2", "q=2tL^2", "q>=1"):
        if not any(token in str(value) for value in domain.values()):
            issues.append(f"critical domain token missing: {token}")


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Domain Extension",
        "Enlarged Source Bounds",
        "Four-Jet Transfer",
        "q>=1 Consequence",
        "4096",
        "27",
        "-3749/10000000",
        "adjacent-cutoff signed splice",
        "No Xi-level edge sign",
        "not a proof",
    ):
        if token not in text:
            issues.append(f"note marker missing: {token}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    validate(payload, issues)
    validate_note(issues)
    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        return 1
    print(SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
