#!/usr/bin/env python3
"""Validate the Burnol tail-discrepancy and dyadic reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.md"
)

REQUIRED_IDS = {
    f"btdr_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "weighted_mobius_data"),
        (2, "limiting_residue"),
        (3, "finite_and_limiting_discrepancies"),
        (4, "tail_divisor_sum"),
        (5, "discrepancy_difference"),
        (6, "omitted_fractional_tail"),
        (7, "abel_tail_operator"),
        (8, "energy_definitions"),
        (9, "hilbert_tail_split"),
        (10, "stable_prefix_energy"),
        (11, "scalar_tail_equivalence"),
        (12, "post_prefix_energy"),
        (13, "three_gate_criterion"),
        (14, "quotient_reindexing"),
        (15, "short_interval_reindexing"),
        (16, "dyadic_block_energy"),
        (17, "dyadic_equivalence"),
        (18, "summable_envelope"),
        (19, "literature_nonpromotion_guards"),
        (20, "open_dyadic_square_function"),
    )
}

REQUIRED_NOTE = (
    "# Jensen-Window PF Burnol Tail-Discrepancy/Dyadic Reduction",
    "D_(omega,N)(k)",
    "T_(omega,N)(k)",
    "g_k(d+1)-g_k(d)",
    "R_(omega,N)=P_(omega,N)+U_(omega,N)",
    "exact three-gate criterion",
    "r_(omega,N)=O(N^(-1/2-omega))",
    "M_alpha(max(N,floor(k/(q+1))))",
    "V_(omega,N)(K)",
    "2^(alpha-2)",
    "2^(-eta*j)",
    "https://arxiv.org/abs/math/0011254",
    "https://arxiv.org/abs/math/0306251",
    "Equation (BTDR.16) is open",
    "does not establish",
)


def mobius(n: int) -> int:
    return int(sp.mobius(n))


def a_value(omega: mp.mpf, d: int) -> mp.mpf:
    return mobius(d) * mp.power(d, -2 * omega)


def partial_a(omega: mp.mpf, n_max: int) -> mp.mpf:
    return mp.fsum(a_value(omega, d) / d for d in range(1, n_max + 1))


def weighted_mertens(omega: mp.mpf, x: int) -> mp.mpf:
    return mp.fsum(a_value(omega, d) for d in range(1, x + 1))


def b_sum(omega: mp.mpf, n_max: int, k: int) -> mp.mpf:
    return mp.fsum(
        a_value(omega, d) * (k // d)
        for d in range(1, min(n_max, k) + 1)
    )


def discrepancy(
    omega: mp.mpf,
    n_max: int,
    k: int,
) -> mp.mpf:
    return partial_a(omega, n_max) * k - b_sum(omega, n_max, k)


def limiting_discrepancy(omega: mp.mpf, k: int) -> mp.mpf:
    return k / mp.zeta(1 + 2 * omega) - b_sum(omega, k, k)


def h_direct(omega: mp.mpf, n_max: int, k: int) -> mp.mpf:
    return mp.fsum(
        a_value(omega, d) * (k // d)
        for d in range(n_max + 1, k + 1)
    )


def h_by_m(omega: mp.mpf, n_max: int, k: int) -> mp.mpf:
    return mp.fsum(
        weighted_mertens(omega, k // m)
        - weighted_mertens(omega, n_max)
        for m in range(1, k // (n_max + 1) + 1)
    )


def h_by_q(omega: mp.mpf, n_max: int, k: int) -> mp.mpf:
    return mp.fsum(
        q
        * (
            weighted_mertens(omega, k // q)
            - weighted_mertens(
                omega,
                max(n_max, k // (q + 1)),
            )
        )
        for q in range(1, k // (n_max + 1) + 1)
    )


def g_value(k: int, d: int) -> mp.mpf:
    return mp.mpf(d) * mp.frac(mp.mpf(k) / d)


def independent_tail_identity_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    omega = mp.mpf(1) / 4
    n_max = 7
    a_infinity = 1 / mp.zeta(1 + 2 * omega)
    r_n = a_infinity - partial_a(omega, n_max)

    for k in (1, 7, 8, 13, 25):
        finite_s = discrepancy(omega, n_max, k)
        infinite_s = limiting_discrepancy(omega, k)
        h_value = h_direct(omega, n_max, k)
        d_value = h_value - k * r_n
        if abs(finite_s - infinite_s - d_value) > mp.mpf("1e-60"):
            issues.append(f"discrepancy split failed at k={k}")

        finite_tail = mp.fsum(
            a_value(omega, d) * mp.frac(mp.mpf(k) / d)
            for d in range(n_max + 1, k + 1)
        )
        finite_tail += k * (
            a_infinity - partial_a(omega, max(n_max, k))
        )
        if abs(finite_tail + d_value) > mp.mpf("1e-60"):
            issues.append(f"fractional tail identity failed at k={k}")

        abel_tail = r_n * g_value(k, n_max + 1)
        abel_tail += mp.fsum(
            (
                a_infinity - partial_a(omega, d)
            )
            * (g_value(k, d + 1) - g_value(k, d))
            for d in range(n_max + 1, k + 1)
        )
        if abs(abel_tail - finite_tail) > mp.mpf("1e-59"):
            issues.append(f"Abel tail identity failed at k={k}")
    return issues


def independent_reindexing_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    omega = mp.mpf(1) / 4
    n_max = 7
    for k in (8, 13, 25, 41):
        direct = h_direct(omega, n_max, k)
        by_m = h_by_m(omega, n_max, k)
        by_q = h_by_q(omega, n_max, k)
        if abs(direct - by_m) > mp.mpf("1e-60"):
            issues.append(f"multiple-count reindexing failed at k={k}")
        if abs(direct - by_q) > mp.mpf("1e-60"):
            issues.append(f"quotient-interval reindexing failed at k={k}")
    return issues


def independent_energy_split_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    omega = mp.mpf(1) / 4
    alpha = 2 * omega
    n_max = 6
    k_max = 48
    a_infinity = 1 / mp.zeta(1 + alpha)
    r_n = a_infinity - partial_a(omega, n_max)

    r_cut = mp.fsum(
        mp.power(k, alpha - 2)
        * (
            discrepancy(omega, n_max, k)
            - limiting_discrepancy(omega, k)
        )
        ** 2
        for k in range(1, k_max + 1)
    )
    prefix = r_n**2 * mp.fsum(
        mp.power(k, alpha) for k in range(1, n_max + 1)
    )
    post = mp.fsum(
        mp.power(k, alpha - 2)
        * (h_direct(omega, n_max, k) - k * r_n) ** 2
        for k in range(n_max + 1, k_max + 1)
    )
    if abs(r_cut - prefix - post) > mp.mpf("1e-58"):
        issues.append("prefix/post-prefix energy split failed")

    dyadic_actual = mp.mpf("0")
    dyadic_unweighted = mp.mpf("0")
    for j in range(3):
        k_left = (2**j) * n_max
        block = mp.fsum(
            (h_direct(omega, n_max, k) - k * r_n) ** 2
            for k in range(k_left + 1, 2 * k_left + 1)
        )
        dyadic_unweighted += mp.power(k_left, alpha - 2) * block
        dyadic_actual += mp.fsum(
            mp.power(k, alpha - 2)
            * (h_direct(omega, n_max, k) - k * r_n) ** 2
            for k in range(k_left + 1, 2 * k_left + 1)
        )
    lower = mp.power(2, alpha - 2) * dyadic_unweighted
    if not (lower <= dyadic_actual <= dyadic_unweighted):
        issues.append("dyadic weight comparison failed")
    return issues


def independent_power_sum_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    alpha = mp.mpf(1) / 2
    for n_max in (1, 7, 30):
        power_sum = mp.fsum(
            mp.power(k, alpha) for k in range(1, n_max + 1)
        )
        lower = mp.power(n_max, 1 + alpha) / (1 + alpha)
        upper = mp.power(n_max, 1 + alpha)
        if not (lower <= power_sum <= upper):
            issues.append(f"power-sum comparison failed at N={n_max}")
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    note = args.note.read_text(encoding="utf-8")
    issues: list[str] = []

    if payload.get("kind") != (
        "jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction"
    ):
        issues.append("bad kind")
    rows = payload.get("rows", [])
    ids = {item.get("id") for item in rows}
    if len(rows) != 20:
        issues.append(f"expected 20 rows, found {len(rows)}")
    if ids != REQUIRED_IDS:
        issues.append(
            f"row ids differ: missing={sorted(REQUIRED_IDS-ids)}, "
            f"extra={sorted(ids-REQUIRED_IDS)}"
        )
    if len(ids) != len(rows):
        issues.append("duplicate row id")
    for item in rows:
        if not item.get("proof_boundary"):
            issues.append(f"{item.get('id')}: missing proof boundary")

    expected_statuses = {
        "btdr_09_hilbert_tail_split": "available_exact",
        "btdr_13_three_gate_criterion": "available_exact",
        "btdr_19_literature_nonpromotion_guards": "source_backed_guard",
        "btdr_20_open_dyadic_square_function": "open_target",
    }
    by_id = {item.get("id"): item for item in rows}
    for row_id, expected in expected_statuses.items():
        if by_id.get(row_id, {}).get("status") != expected:
            issues.append(f"{row_id}: expected status {expected}")

    expected_audit = {
        "row_count": 20,
        "exact_identity_count": 12,
        "equivalence_step_count": 4,
        "literature_guard_count": 2,
        "open_dyadic_square_function_gate_count": 1,
        "uniform_dyadic_bound_proved": False,
        "rh_proved": False,
        "pf_infinity_proved": False,
        "lambda_le_zero_proved": False,
    }
    audit = payload.get("audit", {})
    for key, expected in expected_audit.items():
        if audit.get(key) != expected:
            issues.append(
                f"audit {key}: expected {expected!r}, "
                f"found {audit.get(key)!r}"
            )

    for marker in REQUIRED_NOTE:
        if marker not in note:
            issues.append(f"note missing: {marker}")

    issues.extend(independent_tail_identity_check())
    issues.extend(independent_reindexing_check())
    issues.extend(independent_energy_split_check())
    issues.extend(independent_power_sum_check())

    lowered = note.lower()
    for forbidden in (
        "this proves rh",
        "the dyadic bound is proved",
        "the short-interval gate is closed",
        "therefore lambda <= 0",
    ):
        if forbidden in lowered:
            issues.append(f"forbidden promotion language: {forbidden}")

    for issue in issues:
        print(f"BTDR [{issue}]")
    print(
        "validated Burnol tail-discrepancy/dyadic reduction: "
        f"{len(rows)} rows, {len(issues)} issues, "
        "12 exact identities, 4 equivalence steps, "
        "2 literature guards, 1 open dyadic square-function gate"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
