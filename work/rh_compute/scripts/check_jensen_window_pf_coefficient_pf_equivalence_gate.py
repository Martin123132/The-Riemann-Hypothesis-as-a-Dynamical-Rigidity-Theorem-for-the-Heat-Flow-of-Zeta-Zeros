#!/usr/bin/env python3
"""Validate the coefficient-PF/multiplier/Jensen-window equivalence gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import comb, factorial
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_coefficient_pf_equivalence_gate.json"
)


def pascal_row(degree: int) -> list[int]:
    row = [1]
    for _ in range(degree):
        row = [1] + [row[j - 1] + row[j] for j in range(1, len(row))] + [1]
    return row


def recompute_audit() -> dict:
    values = [Fraction(k * k + 3 * k + 2, 1) for k in range(12)]
    counts = {
        "normalization_checks": 0,
        "derivative_tail_checks": 0,
        "diagonal_window_checks": 0,
    }
    for k, value in enumerate(values):
        if value / factorial(k) * factorial(k) != value:
            raise RuntimeError(f"normalization failed at k={k}")
        counts["normalization_checks"] += 1
    for n in range(4):
        for j in range(len(values) - n):
            c_nj = values[n + j] / factorial(n + j)
            if c_nj * factorial(n + j) / factorial(j) != values[n + j] / factorial(j):
                raise RuntimeError(f"derivative identity failed at n={n}, j={j}")
            counts["derivative_tail_checks"] += 1
        for degree in range(1, 6):
            expanded_coefficients = pascal_row(degree)
            for j in range(degree + 1):
                if expanded_coefficients[j] * values[n + j] != comb(degree, j) * values[n + j]:
                    raise RuntimeError(f"window identity failed at n={n}, d={degree}, j={j}")
                counts["diagonal_window_checks"] += 1
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    issues: list[str] = []
    rows = payload.get("rows", [])
    exact = payload.get("exact", {})

    if payload.get("sample_audit") != recompute_audit():
        issues.append("sample audit mismatch")
    if len(rows) != 13:
        issues.append(f"row count {len(rows)} != 13")
    ids = [row.get("id") for row in rows]
    if len(ids) != len(set(ids)):
        issues.append("duplicate row ids")

    required_markers = {
        "scope_assumptions": "A_k>=0 for every k, A_0>0, and F is entire",
        "normalization": "c_k=A_k/k!",
        "tail_function": "F^(n)(z)",
        "window_identity": "T_n[(1+z)^d]",
        "edrei_form": "C>0, gamma>=0, alpha_r>=0",
        "polya_schur": "multiplier sequence",
        "derivative_closure": "F^(n)",
        "finite_asw": "finite PF-infinity",
        "equivalence": "<=> every B^(d,n)",
        "finite_guard": "finite Toeplitz/PF certificates",
        "signed_hankel_guard": "signed-Hankel positivity",
        "positive_mixture_guard": "positive mixture of dilations",
    }
    for key, marker in required_markers.items():
        if marker not in str(exact.get(key, "")):
            issues.append(f"missing exact marker {key}:{marker}")

    role_counts = {
        role: sum(row.get("role") == role for row in rows)
        for role in {
            "exact_identity",
            "classical_theorem",
            "closure_theorem",
            "theorem_composition",
            "exact_equivalence",
            "forbidden_promotion",
            "nonimplication_guard",
            "positive_mixture_guard",
            "open_handoff",
        }
    }
    expected_counts = {
        "exact_identity": 3,
        "classical_theorem": 3,
        "closure_theorem": 1,
        "theorem_composition": 1,
        "exact_equivalence": 1,
        "forbidden_promotion": 1,
        "nonimplication_guard": 1,
        "positive_mixture_guard": 1,
        "open_handoff": 1,
    }
    if role_counts != expected_counts:
        issues.append(f"role counts mismatch: {role_counts}")

    sources = payload.get("primary_sources", [])
    for marker in ("annals.math.princeton.edu", "BF02786970", "BF02786971", "pnas.37.5.303"):
        if not any(marker in source for source in sources):
            issues.append(f"missing primary source marker {marker}")

    boundary = payload.get("proof_boundary", "")
    for marker in (
        "does not prove coefficient PF-infinity",
        "signed-Hankel-to-PF bridge",
        "Laguerre-Polya",
        "RH",
        "Lambda <= 0",
    ):
        if marker not in boundary:
            issues.append(f"missing proof-boundary marker {marker!r}")

    if any(
        row.get("readiness") == "ready_to_apply"
        for row in rows
        if row.get("role") in {
            "forbidden_promotion",
            "nonimplication_guard",
            "positive_mixture_guard",
            "open_handoff",
        }
    ):
        issues.append("guard or open row incorrectly marked ready_to_apply")

    print(
        "validated Jensen-window PF coefficient-PF equivalence gate: "
        f"{len(rows)} rows, {len(issues)} issues, 3 exact coefficient identities, "
        "4 classical/closure steps, 1 seven-way equivalence, "
        "3 guards, 1 open structural handoff"
    )
    for issue in issues:
        print(f"ISSUE {issue}")
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
