#!/usr/bin/env python3
"""Validate the Edrei-log/Stieltjes endpoint-equivalence gate independently."""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import permutations
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_edrei_stieltjes_equivalence_gate.json"
)


def parse_fraction(text: str) -> Fraction:
    return Fraction(text)


def permutation_sign(values: tuple[int, ...]) -> int:
    inversions = sum(
        values[i] > values[j]
        for i in range(len(values))
        for j in range(i + 1, len(values))
    )
    return -1 if inversions % 2 else 1


def determinant_by_permutations(matrix: list[list[Fraction]]) -> Fraction:
    size = len(matrix)
    return sum(
        permutation_sign(perm)
        * product(matrix[row][perm[row]] for row in range(size))
        for perm in permutations(range(size))
    )


def product(values) -> Fraction:
    out = Fraction(1)
    for value in values:
        out *= value
    return out


def recompute_audit(payload: dict) -> list[str]:
    issues: list[str] = []
    audit = payload.get("sample_audit", {})
    degree = int(audit.get("degree", 0))
    gamma = parse_fraction(audit.get("gamma", "0"))
    atoms = [
        (parse_fraction(row["beta"]), int(row["multiplicity"]))
        for row in audit.get("atoms", [])
    ]

    moments = [
        (gamma if n == 0 else Fraction(0))
        + sum(multiplicity * beta ** (n + 1) for beta, multiplicity in atoms)
        for n in range(degree)
    ]
    persisted_moments = [
        parse_fraction(value) for value in audit.get("shifted_log_moments", [])
    ]
    if persisted_moments != moments:
        issues.append("shifted log moments mismatch")

    expected_log_derivative = [(-1) ** n * value for n, value in enumerate(moments)]
    persisted_log_derivative = [
        parse_fraction(value)
        for value in audit.get("log_derivative_coefficients", [])
    ]
    if persisted_log_derivative != expected_log_derivative:
        issues.append("log-derivative coefficient mismatch")

    h = [Fraction(0) for _ in range(degree + 1)]
    h[0] = Fraction(1)
    for n in range(1, degree + 1):
        h[n] = sum(
            expected_log_derivative[k - 1] * h[n - k] for k in range(1, n + 1)
        ) / n
    persisted_h = [parse_fraction(value) for value in audit.get("h_coefficients", [])]
    if persisted_h != h:
        issues.append("independent H-series recurrence mismatch")

    h0 = [
        determinant_by_permutations(
            [[moments[i + j] for j in range(size)] for i in range(size)]
        )
        for size in range(1, 5)
    ]
    h1 = [
        determinant_by_permutations(
            [[moments[i + j + 1] for j in range(size)] for i in range(size)]
        )
        for size in range(1, 4)
    ]
    persisted_h0 = [
        parse_fraction(value)
        for value in audit.get("hankel_s1_determinants_sizes_1_to_4", [])
    ]
    persisted_h1 = [
        parse_fraction(value)
        for value in audit.get("hankel_s2_determinants_sizes_1_to_3", [])
    ]
    if persisted_h0 != h0:
        issues.append("s=1 Hankel determinants mismatch")
    if persisted_h1 != h1:
        issues.append("s=2 Hankel determinants mismatch")
    if audit.get("indexing_checks") != degree:
        issues.append("indexing check count mismatch")
    if audit.get("hankel_checks") != 7:
        issues.append("Hankel check count mismatch")
    if not all(value > 0 for value in h0[:3]) or h0[3] != 0:
        issues.append("unexpected s=1 finite-support sign pattern")
    if not all(value > 0 for value in h1[:2]) or h1[2] != 0:
        issues.append("unexpected s=2 finite-support sign pattern")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    issues = recompute_audit(payload)
    rows = payload.get("rows", [])
    exact = payload.get("exact", {})

    if len(rows) != 14:
        issues.append(f"row count {len(rows)} != 14")
    ids = [row.get("id") for row in rows]
    if len(ids) != len(set(ids)):
        issues.append("duplicate row ids")

    required_markers = {
        "corpus_normalization": "H_lambda(0)=1",
        "log_coordinates": "p_n=(-1)^(n-1)*q_n",
        "indexing_identity": "=p_(r+1)",
        "sokal_criterion": "Stieltjes moment sequence",
        "unique_measure": "m_j in Z_(>=1)",
        "stieltjes_transform": "H'(z)/H(z)=integral",
        "hankel_psd": "positive semidefinite",
        "strict_two_column_target": "D_(m,2)",
        "s_fraction": "continued fraction",
        "endpoint_equivalence": "every shifted Jensen window",
        "phi_integral": "Phi(u)*cosh",
        "phi_log_derivative": "sinh(u*sqrt(z))",
        "finite_guard": "finitely many positive",
        "object_separation": "not the original signed-Hankel matrices",
        "mixture_guard": "positive mixture",
    }
    for key, marker in required_markers.items():
        if marker not in str(exact.get(key, "")):
            issues.append(f"missing exact marker {key}:{marker}")

    expected_roles = {
        "exact_identity": 2,
        "classical_theorem": 4,
        "exact_equivalence": 2,
        "sufficient_theorem_target": 1,
        "finite_evidence_map": 1,
        "object_separation_guard": 1,
        "forbidden_promotion": 1,
        "xi_phi_exact_handoff": 1,
        "positive_mixture_guard": 1,
    }
    role_counts = {
        role: sum(row.get("role") == role for row in rows) for role in expected_roles
    }
    if role_counts != expected_roles:
        issues.append(f"role counts mismatch: {role_counts}")

    sources = payload.get("primary_sources", [])
    for marker in ("jmaa.2022.126432", "Sokal_1-s2.0", "BF02786970", "BF02786971"):
        if not any(marker in source for source in sources):
            issues.append(f"missing primary source marker {marker}")

    boundary = payload.get("proof_boundary", "")
    for marker in (
        "does not prove the all-order Stieltjes property",
        "coefficient PF-infinity",
        "Jensen hyperbolicity",
        "RH",
        "Lambda <= 0",
    ):
        if marker not in boundary:
            issues.append(f"missing proof-boundary marker {marker!r}")

    forbidden_ready_roles = {
        "sufficient_theorem_target",
        "finite_evidence_map",
        "object_separation_guard",
        "forbidden_promotion",
        "xi_phi_exact_handoff",
        "positive_mixture_guard",
    }
    if any(
        row.get("readiness") == "ready_to_apply"
        for row in rows
        if row.get("role") in forbidden_ready_roles
    ):
        issues.append("open, finite, or guard row incorrectly marked ready_to_apply")

    print(
        "validated Jensen-window PF Edrei-Stieltjes equivalence gate: "
        f"{len(rows)} rows, {len(issues)} issues, 12 exact indexing checks, "
        "7 exact Hankel checks, 1 unified endpoint, 3 finite/nonpromotion guards, "
        "1 open Xi/Phi handoff"
    )
    for issue in issues:
        print(f"ISSUE: {issue}")
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
