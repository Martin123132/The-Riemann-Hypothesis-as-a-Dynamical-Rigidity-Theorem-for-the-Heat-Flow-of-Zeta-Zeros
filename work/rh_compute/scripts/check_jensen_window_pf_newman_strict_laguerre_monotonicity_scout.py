#!/usr/bin/env python3
"""Validate the strict-Laguerre monotonicity target and numerical scout."""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import sympy as sp

import jensen_window_pf_newman_strict_laguerre_monotonicity_scout as target


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_strict_laguerre_monotonicity_scout.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_strict_laguerre_monotonicity_scout.md"
EXPECTED_IDS = [f"nslms_{index:02d}_{suffix}" for index, suffix in enumerate(
    (
        "derivative_identity",
        "monotone_sufficiency",
        "correlation_sine_target",
        "tail_cosine_geometry",
        "theta_primitive_target",
        "dense_scout",
        "high_frequency_scout",
        "theta_tail_nonpromotion",
        "xi_lehmer_counterexample",
    ),
    start=1,
)]


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]
    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = target.build_payload()
    if stored != rebuilt:
        issues.append("stored payload differs from reconstruction")
    if stored.get("kind") != "jensen_window_pf_newman_strict_laguerre_monotonicity_scout":
        issues.append("kind drifted")
    rows = stored.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")

    x = sp.symbols("x", real=True)
    h = sp.Function("h")(x)
    laguerre = sp.diff(h, x) ** 2 - h * sp.diff(h, x, 2)
    expected = sp.diff(h, x) * sp.diff(h, x, 2) - h * sp.diff(h, x, 3)
    if sp.simplify(sp.diff(laguerre, x) - expected) != 0:
        issues.append("independent Laguerre derivative identity failed")

    exact = stored.get("exact", {})
    for phrase in ("integral_x^infinity", "Lambda<=0", "M_t(c)=0"):
        joined = " ".join(str(value) for value in exact.values())
        if phrase not in joined:
            issues.append(f"exact target missing {phrase}")
    tail = exact.get("tail_cosine_reduction", {})
    for phrase in (
        "8*x*integral",
        "Q_t''(0)=-K_(1,t)(0)<0",
        "exactly once",
        "Q_t(i*p) real",
        "pi/8",
    ):
        if phrase not in " ".join(str(value) for value in tail.values()):
            issues.append(f"tail-cosine reduction missing {phrase}")

    moderate = stored.get("moderate_scout", {})
    moderate_rows = moderate.get("rows", [])
    if len(moderate_rows) != 6:
        issues.append("expected six dense time rows")
    if any(row.get("x_count") != 2000 for row in moderate_rows):
        issues.append("dense grid size drifted")
    if any(row.get("all_L_prime_negative") is not True for row in moderate_rows):
        issues.append("dense scout contains a nonnegative derivative")
    if float(moderate.get("max_relative_coarse_fine_delta", "inf")) >= 1e-6:
        issues.append("dense coarse/fine convergence failed")

    high = stored.get("high_frequency_scout", {})
    high_rows = high.get("rows", [])
    if len(high_rows) != 20:
        issues.append("expected twenty high-frequency rows")
    if any(mp.mpf(row.get("L_prime", "0")) >= 0 for row in high_rows):
        issues.append("high-frequency scout contains a nonnegative derivative")
    if mp.mpf(high.get("max_relative_coarse_fine_delta", "inf")) >= mp.mpf("1e-25"):
        issues.append("high-frequency coarse/fine convergence failed")
    if mp.mpf(high.get("max_relative_t0_xi_delta", "inf")) >= mp.mpf("1e-25"):
        issues.append("high-frequency Xi cross-check failed")

    guard = stored.get("countermodel_guard", {})
    for phrase in ("theta-type", "L'(y)>0", "21/5"):
        if phrase not in " ".join(str(value) for value in guard.values()):
            issues.append(f"countermodel guard missing {phrase}")

    xi_guard = stored.get("xi_monotonicity_counterexample", {})
    for phrase in ("M_0(x)<0", "sufficiently small positive t", "strict-monotonicity route"):
        if phrase not in " ".join(str(value) for value in xi_guard.values()):
            issues.append(f"Xi monotonicity counterexample missing {phrase}")

    note = NOTE.read_text(encoding="utf-8")
    for phrase in (
        "not a proof of RH",
        "Exact Sufficient Candidate",
        "Dense Moderate Scout",
        "High-Frequency Scout",
        "Nonpromotion Guard",
        "classical convex-kernel positivity theorem",
        "10.1017/S0004972700047511",
        "Exact Xi Counterexample",
        "Dominated convergence",
        "Do not pursue a global positive sine-transform theorem",
        "0<t<=1/5",
        "v*K_(1,t)(v)",
    ):
        if phrase not in note:
            issues.append(f"note missing required phrase: {phrase}")
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Jensen-window PF Newman strict-Laguerre monotonicity scout: "
        "9 rows, 0 issues, 5 exact target identities, "
        "2 exact classical-route obstructions, 6 dense time rows, "
        "20 high-frequency rows, 1 theta-tail nonpromotion guard, "
        "1 Arb Xi monotonicity rejection"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
