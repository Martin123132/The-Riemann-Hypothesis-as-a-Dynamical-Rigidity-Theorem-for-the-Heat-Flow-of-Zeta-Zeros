#!/usr/bin/env python3
"""Validate the contiguous terminal-tail recurrence/current theorem."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_recurrence_current_gate as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [f"cttr_{index:02d}_{suffix}" for index, suffix in enumerate((
    "path",
    "tail",
    "recurrence",
    "phase",
    "trace",
    "jets",
    "current",
    "base",
    "away",
    "near",
    "domain",
    "compact",
    "theorem",
    "endpoint",
    "pi",
    "handoff",
), start=1)]
EXPECTED_COUNTS = {
    "rows": 16,
    "exact_endpoint_recurrences": 1,
    "exact_tail_collapses": 1,
    "grouped_jet_identities": 3,
    "arb_compact_boxes": 1024,
    "fixed_finite_tail_current_theorems": 1,
    "growing_tail_estimates": 0,
    "complete_aggregate_cross_current_bounds": 0,
    "abel_gaps": 0,
    "contact_exclusions": 0,
}
SUMMARY = (
    "validated contiguous terminal-tail recurrence/current gate: "
    "16 rows, 1 exact C0 tail collapse, 3 grouped jets, 1024 Arb boxes, "
    "1 all-fixed-M clockwise theorem, 0 growing-tail estimates, "
    "0 aggregate closures"
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


def independent_symbolic_audit(issues: list[str]) -> None:
    p, m = sp.symbols("p m", real=True)
    phi = p**2 / 2 - p + sp.Rational(3, 8)
    shifted = p**2 / 2 - (2 * m + 1) * p + sp.Rational(3, 8)
    if sp.expand(-phi + 2 * m * p + shifted) != 0:
        issues.append("independent phase matching failed")

    c, cp, cpp = sp.symbols("C C_p C_pp", real=True)
    pi = sp.pi
    d = -cp / (4 * pi)
    n = cpp / (16 * pi**2)
    if sp.simplify(c * n - d**2 - (c * cpp - cp**2) / (16 * pi**2)) != 0:
        issues.append("independent grouped-current identity failed")

    t, tangent = sp.symbols("t tangent", positive=True, real=True)
    quadratic = pi**2 * t**2 * (1 + tangent**2) + pi * tangent
    minimum = pi**2 * t**2 - 1 / (4 * t**2)
    square = pi**2 * t**2 * (tangent + 1 / (2 * pi * t**2)) ** 2
    if sp.simplify(quadratic - minimum - square) != 0:
        issues.append("independent completed-square identity failed")


def independent_endpoint_audit(issues: list[str]) -> None:
    flint.ctx.prec = gate.PRECISION_BITS
    pi = flint.arb.pi()
    p = gate.aq(1)
    currents = []
    for k in range(1, 7):
        variable = flint.arb_series([p, flint.arb(1)], prec=3)
        phase = variable * variable / 2 - 2 * k * variable + flint.arb(3) / 8
        trace = (-1) ** k * phase.cos_pi() / (2 * variable.cos_pi())
        coefficients = trace.coeffs()
        value = coefficients[0]
        first = coefficients[1]
        second = 2 * coefficients[2]
        currents.append((value * second - first**2) / (16 * pi**2))
    base = currents[0]
    for k, current in enumerate(currents, start=1):
        expected = base - gate.aq(Fraction(k * (k - 1), 16))
        if not (current - expected).contains(0):
            issues.append(f"independent endpoint diagnostic failed at K={k}")


def validate(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    if payload.get("date") != "2026-07-31":
        issues.append("artifact date drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    if [row.get("id") for row in payload.get("rows", [])] != EXPECTED_IDS:
        issues.append("row ids or order drifted")

    for token in ("contiguous-tail", "fixed-M", "no growing-tail", "or RH"):
        if token not in payload.get("status", ""):
            issues.append(f"status missing token: {token}")
    for token in (
        "M grows with N",
        "complete Xi aggregate",
        "Abel gap",
        "contact exclusion",
        "Lambda<=0",
        "RH",
    ):
        if token not in payload.get("proof_boundary", ""):
            issues.append(f"proof boundary missing token: {token}")

    hashes = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    try:
        if payload.get("symbolic_certificate") != gate.symbolic_certificate():
            issues.append("stored symbolic certificate drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"symbolic recomputation failed: {exc}")

    try:
        interval = gate.compact_box_certificate()
        stored = payload.get("interval_certificate", {})
        for key in (
            "domain",
            "boxes",
            "strict_margin",
            "minimum_lower",
            "minimum_box",
            "sinc_tail",
            "e_strip_minima",
        ):
            if stored.get(key) != interval.get(key):
                issues.append(f"compact interval certificate drifted: {key}")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"compact interval recomputation failed: {exc}")

    exact = payload.get("exact", {})
    if exact.get("theorem") != (
        "For every integer M>=0 and every p in [-1,1], D(2M+2-p)>3/50 "
        "and the fixed-M contiguous terminal-tail leading symbol "
        "satisfies P_M(p)<-3/8000."
    ):
        issues.append("theorem statement drifted")
    if exact.get("endpoint_diagnostic") != "P_M(1)=P_0(1)-M*(M+1)/16.":
        issues.append("endpoint diagnostic drifted")

    independent_symbolic_audit(issues)
    independent_endpoint_audit(issues)


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Contiguous Terminal-Tail Recurrence and Current Gate",
        "Exact Collapse",
        "Uniform Curvature",
        "q_z^2-q*q_zz>7/5",
        "P_M(1)=P_0(1)-M*(M+1)/16",
        "`M` fixed before `N -> infinity`",
        "No growing-prefix finite-height bound",
    ):
        if token not in text:
            issues.append(f"note marker missing: {token}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
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
