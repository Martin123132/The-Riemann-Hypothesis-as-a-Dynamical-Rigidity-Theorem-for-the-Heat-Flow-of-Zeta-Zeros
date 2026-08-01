#!/usr/bin/env python3
"""Independently validate the real-edge projective-current gate."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_real_edge_projective_current_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_real_edge_projective_current_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    "repc_01_edge_coordinates",
    "repc_02_division_free_current",
    "repc_03_second_jet",
    "repc_04_q1_terminal_phase",
    "repc_05_phase_guard",
    "repc_06_asymptotic_jets",
    "repc_07_derivative_collapse",
    "repc_08_current_limit",
    "repc_09_zero_fibre",
    "repc_10_real_source",
    "repc_11_removable_charts",
    "repc_12_interval_cover",
    "repc_13_uniform_symbol_sign",
    "repc_14_handoff",
]
EXPECTED_COUNTS = {
    "rows": 14,
    "division_free_edge_currents": 1,
    "asymptotic_edge_jet_limits": 4,
    "arb_intervals": 3072,
    "removable_charts": 2,
    "strict_curvature_margins": 1,
    "uniform_q1_leading_symbol_signs": 1,
    "uniform_finite_height_edge_signs": 0,
    "abel_gaps": 0,
    "winding_bounds": 0,
}
SUMMARY = (
    "validated real-edge projective-current gate: 14 rows, "
    "1 division-free edge current, 4 asymptotic edge jets, "
    "3072 Arb intervals, 2 removable charts, 1 strict curvature margin, "
    "1 uniform q=1 leading sign, 0 finite-height edge signs, "
    "0 Abel gaps, 0 winding bounds"
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


def validate_symbolics(issues: list[str]) -> None:
    c, cx, gr, grx, alpha, alpha_x = sp.symbols(
        "C C_x G_R G_xR alpha alpha_x", real=True
    )
    d = gr - alpha * c
    dx = grx - alpha_x * c - alpha * cx
    if sp.expand(c * dx - d * cx - (c * grx - gr * cx - alpha_x * c**2)) != 0:
        issues.append("division-free centering cancellation failed")

    p = sp.symbols("p", real=True)
    theta = (1 - p) / 2
    psi = sp.pi * (p**2 / 2 + sp.Rational(3, 8))
    delta = sp.pi * (p**2 / 2 - p + sp.Rational(3, 8))
    f = sp.cos(psi) / (2 * sp.cos(sp.pi * p))
    a = f - sp.cos(delta)
    b = -sp.diff(f, p) / (4 * sp.pi) + theta * sp.sin(delta) / 2
    m = (
        sp.diff(f, p, 2) / (16 * sp.pi**2)
        + sp.sin(delta) / (16 * sp.pi)
        + theta**2 * sp.cos(delta) / 4
    )
    if sp.simplify(b + sp.diff(a, p) / (4 * sp.pi)) != 0:
        issues.append("B derivative collapse failed")
    if sp.simplify(m - sp.diff(a, p, 2) / (16 * sp.pi**2)) != 0:
        issues.append("M derivative collapse failed")
    if sp.simplify(
        a * m
        - b**2
        - (a * sp.diff(a, p, 2) - sp.diff(a, p) ** 2)
        / (16 * sp.pi**2)
    ) != 0:
        issues.append("log-curvature current failed")

    q_phase = -sp.pi * (p**2 / 2 - p + sp.Rational(3, 8))
    r_conjugate_phase = -sp.pi * (
        p**2 / 2 + p + sp.Rational(3, 8)
    )
    phase_delta = sp.expand(q_phase - r_conjugate_phase + 2 * sp.pi * (1 - p))
    if sp.simplify(phase_delta) != 2 * sp.pi:
        issues.append("terminal phase guard exponent failed")


def validate_payload(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or order drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")

    status = payload.get("status", "")
    for token in ("division-free", "q=1", "negative leading-symbol"):
        if token not in status:
            issues.append(f"status missing token: {token}")

    boundary = payload.get("proof_boundary", "")
    for token in (
        "does not yet prove a uniform finite-a edge sign",
        "q>1 edge sign",
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

    interval = payload.get("interval", {})
    recomputed = gate.interval_certificate()
    if interval != recomputed:
        issues.append("interval certificate drifted")
    if interval.get("boxes") != 3072:
        issues.append("interval count drifted")
    if interval.get("curvature_margin") != ">3/50":
        issues.append("curvature margin drifted")
    if interval.get("current_symbol_margin") != "K_edge<-3/8000":
        issues.append("current margin drifted")
    if interval.get("removable_centers") != ["-1/2", "1/2"]:
        issues.append("removable centers drifted")

    points = interval.get("points", {})
    try:
        wrong = gate.flint.arb(
            points["phase_guard"]["wrong_terminal_phase_current_ball"]
        )
        correct = gate.flint.arb(points["phase_guard"]["current_symbol_ball"])
        if not wrong > 0:
            issues.append("wrong-phase positive guard drifted")
        if not correct < 0:
            issues.append("corrected phase sign drifted")
    except (KeyError, ValueError) as exc:
        issues.append(f"phase guard diagnostics invalid: {exc}")


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Exact Edge Current",
        "Correct Terminal Phase",
        "Three-Jet Collapse",
        "Whole-Cell Certificate",
        "3/8000",
        "3072",
        "finite-`a` error bound",
        "reserve. The Abel gap, winding cap, contact exclusion, and RH",
    ):
        if token not in text:
            issues.append(f"note marker missing: {token}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    validate_symbolics(issues)
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
