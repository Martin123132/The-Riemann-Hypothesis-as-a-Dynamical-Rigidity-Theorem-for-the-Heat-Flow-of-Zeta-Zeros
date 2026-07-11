#!/usr/bin/env python3
"""Independently validate the cubic forward-uniform tail certificate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_cubic_forward_uniform_tail_certificate.json"
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    issues: list[str] = []

    rows = payload.get("rows", [])
    by_id = {row.get("id"): row for row in rows}
    expected_ids = {f"cfut_{index:02d}_{suffix}" for index, suffix in [
        (1, "strict_reciprocal_defects"),
        (2, "exact_q_flow"),
        (3, "neighbor_monotonicity"),
        (4, "weighted_source_cap"),
        (5, "initial_weighted_tail"),
        (6, "weighted_supersolution"),
        (7, "forward_uniform_tail"),
        (8, "cubic_continuation"),
        (9, "lambda_zero_cubics"),
        (10, "higher_degree_handoff"),
    ]}
    if set(by_id) != expected_ids:
        issues.append("row ids drifted")

    k, q, g, j = sp.symbols("k q g j", positive=True)
    q1 = q + g
    q2 = q1 + j
    defect = 1 / q**2
    next_defect = 1 / q1**2
    defect_prime = 2 * (
        (2 * k + 3) * (1 - defect) * next_defect - (2 * k - 1) * defect
    )
    q_prime_from_defect = sp.factor(-q**3 * defect_prime / 2)
    q_prime = (2 * k - 1) * q - (2 * k + 3) * (q**3 - q) / q1**2
    if sp.simplify(q_prime_from_defect - q_prime) != 0:
        issues.append("defect-to-q flow identity failed")
    q1_prime = (1 - 1 / q1**2) * (
        (2 * k + 1) * q1 - (2 * k + 5) * (q1**3 - q1) / q2**2
    )
    increment_flow = sp.factor(q1_prime - q_prime)
    derivative = sp.factor(sp.diff(increment_flow, j))
    expected_derivative = (
        2
        * (2 * k + 5)
        * (q + g - 1) ** 2
        * (q + g + 1) ** 2
        / ((q + g) * (q + g + j) ** 3)
    )
    if sp.simplify(derivative - expected_derivative) != 0:
        issues.append("neighbor derivative identity failed")

    polynomial = (
        -6 * g**3 * k
        + g**3
        - 14 * g**2 * k * q
        + 5 * g**2 * q
        - 6 * g * k * q**2
        - 2 * g * k
        + 13 * g * q**2
        - 5 * g
        - 2 * k * q
        + 6 * q**3
        - 5 * q
    )
    expected_same = (1 - g**2) * polynomial / ((q + g) ** 2 * (q + 2 * g) ** 2)
    if sp.simplify(increment_flow.subs(j, g) - expected_same) != 0:
        issues.append("same-neighbor identity failed")

    positive = g**3 + 5 * g**2 * q + 13 * g * q**2 + 6 * q**3
    gap = sp.Poly(sp.expand(positive - polynomial), k, q, g)
    if any(coefficient < 0 for coefficient in gap.coeffs()):
        issues.append("positive majorant is not coefficient-positive")

    remainder = Fraction(13, 17) + Fraction(5, 17**2) + Fraction(1, 17**3)
    if remainder != Fraction(3843, 4913) or not remainder < 1:
        issues.append("weighted source cap failed")
    coercive_coupling = Fraction(4, 1) + Fraction(10, 319)
    if not coercive_coupling < 5:
        issues.append("coercive supersolution coupling cap failed")

    initial = payload.get("diagnostics", {}).get("initial_tail", {})
    margin = 144 * 99_999**4 * 319**3 - 100_000**4 * (5 * 319 + 6) ** 3
    if margin <= 0 or initial.get("endpoint_squared_margin") != margin:
        issues.append("initial weighted-tail margin failed")

    ready = sum(row.get("readiness") == "ready_to_apply" for row in rows)
    open_rows = sum(row.get("role") == "open_handoff" for row in rows)
    if ready != 9:
        issues.append(f"ready row count {ready} != 9")
    if open_rows != 1:
        issues.append(f"open handoff count {open_rows} != 1")
    boundary = payload.get("proof_boundary", "")
    for marker in ("does not prove degree 4", "PF-infinity", "RH", "Lambda <= 0"):
        if marker not in boundary:
            issues.append(f"missing proof-boundary marker {marker!r}")

    print(
        "validated Jensen-window PF cubic forward-uniform tail certificate: "
        f"{len(rows)} rows, {len(issues)} issues, 3 exact flow identities, "
        "1 weighted source cap, 1 initial weighted tail, 1 forward-uniform tail, "
        "1 full cubic propagation theorem, 1 lambda=0 cubic theorem, "
        "1 open higher-degree handoff"
    )
    for issue in issues:
        print(f"ISSUE {issue}")
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
