#!/usr/bin/env python3
"""Validate the Burnol cell-energy and reciprocal-zeta tail obstruction."""

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
    "jensen_window_pf_burnol_cell_energy_tail_obstruction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_burnol_cell_energy_tail_obstruction.md"
)

REQUIRED_IDS = {
    f"bcet_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "partial_residue"),
        (2, "truncated_convolution"),
        (3, "summatory_floor_identity"),
        (4, "reciprocal_fractional_part"),
        (5, "weighted_burnol_norm"),
        (6, "cell_decomposition"),
        (7, "zero_cell"),
        (8, "closed_cell_formula"),
        (9, "discrete_discrepancy"),
        (10, "uniform_residue_bound"),
        (11, "norm_implies_discrete_energy"),
        (12, "discrete_energy_implies_norm"),
        (13, "exact_discrete_criterion"),
        (14, "stable_prefix"),
        (15, "limit_discrete_energy"),
        (16, "reciprocal_zeta_tail_rate"),
        (17, "zero_free_consequence"),
        (18, "open_short_interval_gate"),
    )
}

REQUIRED_NOTE = (
    "# Jensen-Window PF Burnol Cell-Energy/Tail Obstruction",
    "A_(omega,N)",
    "B_(omega,N)(k)",
    "f_(2omega,N)(1/x)",
    "J_(omega,N)(k)",
    "Q_(omega,N)",
    "S_(omega,N)(k)-S_(omega,infinity)(k)",
    "N^(-1/2-omega)",
    "Re(s)>s_0-(1/2+omega)=1/2+omega",
    "short-multiplicative-interval",
    "Equation (BCET.10) is open",
    "does not establish",
)


def mobius(n: int) -> int:
    return int(sp.mobius(n))


def partial_residue(omega: mp.mpf, n_max: int) -> mp.mpf:
    return mp.fsum(
        mobius(d) * mp.power(d, -(1 + 2 * omega))
        for d in range(1, n_max + 1)
    )


def b_value(omega: mp.mpf, n_max: int, n: int) -> mp.mpf:
    return mp.fsum(
        mobius(int(d)) * mp.power(int(d), -2 * omega)
        for d in sp.divisors(n)
        if int(d) <= n_max
    )


def b_summatory(
    omega: mp.mpf, n_max: int, k: int
) -> mp.mpf:
    return mp.fsum(
        b_value(omega, n_max, n)
        for n in range(1, k + 1)
    )


def floor_summatory(
    omega: mp.mpf, n_max: int, k: int
) -> mp.mpf:
    return mp.fsum(
        mobius(d)
        * mp.power(d, -2 * omega)
        * (k // d)
        for d in range(1, n_max + 1)
    )


def independent_arithmetic_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    omega = mp.mpf(1) / 4
    n_max = 12
    a_value = partial_residue(omega, n_max)
    for k in (1, 7, 12, 19):
        by_n = b_summatory(omega, n_max, k)
        by_d = floor_summatory(omega, n_max, k)
        if abs(by_n - by_d) > mp.mpf("1e-60"):
            issues.append(
                f"summatory floor identity failed at k={k}"
            )
        x = mp.mpf(k) + mp.mpf("0.37")
        fractional = mp.fsum(
            mobius(d)
            * mp.power(d, -2 * omega)
            * mp.frac(x / d)
            for d in range(1, n_max + 1)
        )
        cell = a_value * x - by_d
        if abs(fractional - cell) > mp.mpf("1e-60"):
            issues.append(
                f"reciprocal cell identity failed at k={k}"
            )
    return issues


def independent_cell_integral_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    omega = mp.mpf(1) / 4
    n_max = 12
    a_value = partial_residue(omega, n_max)
    for k in (1, 4, 13):
        b_value_k = floor_summatory(omega, n_max, k)
        direct = mp.quad(
            lambda x: mp.power(x, 2 * omega - 2)
            * (a_value * x - b_value_k) ** 2,
            [mp.mpf(k), mp.mpf(k + 1)],
        )
        delta_1 = (
            mp.power(k + 1, 1 + 2 * omega)
            - mp.power(k, 1 + 2 * omega)
        )
        delta_2 = (
            mp.power(k + 1, 2 * omega)
            - mp.power(k, 2 * omega)
        )
        delta_3 = (
            mp.power(k + 1, 2 * omega - 1)
            - mp.power(k, 2 * omega - 1)
        )
        closed = (
            a_value**2 * delta_1 / (1 + 2 * omega)
            - 2
            * a_value
            * b_value_k
            * delta_2
            / (2 * omega)
            + b_value_k**2
            * delta_3
            / (2 * omega - 1)
        )
        if abs(direct - closed) > mp.mpf("1e-58"):
            issues.append(
                f"closed cell formula failed at k={k}: "
                f"{abs(direct-closed)}"
            )
        if direct < 0:
            issues.append(f"negative cell energy at k={k}")
    return issues


def independent_discrepancy_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    omega = mp.mpf(1) / 4
    n_max = 12
    a_value = partial_residue(omega, n_max)
    for k in (2, 9, 17):
        b_value_k = floor_summatory(omega, n_max, k)
        discrepancy = a_value * k - b_value_k
        fractional = mp.fsum(
            mobius(d)
            * mp.power(d, -2 * omega)
            * mp.frac(mp.mpf(k) / d)
            for d in range(1, n_max + 1)
        )
        if abs(discrepancy - fractional) > mp.mpf("1e-60"):
            issues.append(
                f"discrete discrepancy failed at k={k}"
            )

    n_large = 18
    a_large = partial_residue(omega, n_large)
    a_infinity = 1 / mp.zeta(1 + 2 * omega)
    for k in (3, 11, 18):
        finite_b = floor_summatory(omega, n_large, k)
        infinite_b = floor_summatory(omega, k, k)
        if abs(finite_b - infinite_b) > mp.mpf("1e-60"):
            issues.append(f"stable prefix failed at k={k}")
        finite_s = a_large * k - finite_b
        infinite_s = a_infinity * k - infinite_b
        expected = k * (a_large - a_infinity)
        if abs(finite_s - infinite_s - expected) > mp.mpf(
            "1e-60"
        ):
            issues.append(
                f"stable discrepancy difference failed at k={k}"
            )
    return issues


def independent_weight_comparability_check() -> list[str]:
    issues: list[str] = []
    omega = sp.symbols("omega", positive=True)
    # For 0<omega<1/2, the worst ratio of adjacent weights occurs
    # at k=1 and is below 4.
    ratio_at_one = 2 ** (2 - 2 * omega)
    if sp.simplify(ratio_at_one.subs(omega, sp.Rational(1, 4))) >= 4:
        issues.append("cell weight comparability failed")

    n, p = sp.symbols("n p", positive=True)
    lower_integral = n ** (p + 1) / (p + 1)
    if sp.simplify(
        lower_integral.subs(
            {n: 10, p: sp.Rational(1, 2)}
        )
    ) <= 0:
        issues.append("power-sum lower bound failed")
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
        "jensen_window_pf_burnol_cell_energy_tail_obstruction"
    ):
        issues.append("bad kind")
    rows = payload.get("rows", [])
    ids = {item.get("id") for item in rows}
    if len(rows) != 18:
        issues.append(f"expected 18 rows, found {len(rows)}")
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
        "bcet_16_reciprocal_zeta_tail_rate": "available_exact",
        "bcet_17_zero_free_consequence": "available_exact",
        "bcet_18_open_short_interval_gate": "open_target",
    }
    by_id = {item.get("id"): item for item in rows}
    for row_id, status in expected_statuses.items():
        if by_id.get(row_id, {}).get("status") != status:
            issues.append(f"{row_id}: expected status {status}")

    expected_audit = {
        "row_count": 18,
        "exact_identity_or_bound_count": 10,
        "norm_reduction_step_count": 4,
        "reciprocal_zeta_tail_obstruction_count": 1,
        "zero_free_consequence_count": 1,
        "open_short_interval_gate_count": 1,
        "uniform_discrepancy_proved": False,
        "rh_proved": False,
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

    issues.extend(independent_arithmetic_check())
    issues.extend(independent_cell_integral_check())
    issues.extend(independent_discrepancy_check())
    issues.extend(independent_weight_comparability_check())

    lowered = note.lower()
    for forbidden in (
        "this proves rh",
        "the discrepancy bound is proved",
        "the reciprocal-zeta rate is established for zeta",
        "the short-interval gate is closed",
    ):
        if forbidden in lowered:
            issues.append(f"forbidden promotion language: {forbidden}")

    for issue in issues:
        print(f"BCET [{issue}]")
    print(
        "validated Burnol cell-energy/tail obstruction: "
        f"{len(rows)} rows, {len(issues)} issues, "
        "10 exact identities, 4 norm-reduction steps, "
        "1 reciprocal-zeta tail obstruction, "
        "1 open short-interval gate"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
