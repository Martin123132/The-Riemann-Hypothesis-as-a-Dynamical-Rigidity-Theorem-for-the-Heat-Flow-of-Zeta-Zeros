#!/usr/bin/env python3
"""Validate the Jordan-Muntz/Burnol Hardy-intertwiner audit."""

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
    "jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.md"
)

REQUIRED_IDS = {
    f"jmbh_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "fractional_part_kernel"),
        (2, "fractional_power_abel_identity"),
        (3, "hardy_operator_definition"),
        (4, "remainder_intertwining"),
        (5, "dilation_covariance"),
        (6, "finite_jordan_factorization"),
        (7, "burnol_partial_sum"),
        (8, "reciprocal_unitary"),
        (9, "tail_hardy_operator"),
        (10, "exact_burnol_intertwiner"),
        (11, "hardy_norm"),
        (12, "two_sided_norm_equivalence"),
        (13, "invertibility"),
        (14, "mellin_multiplier"),
        (15, "transform_consistency"),
        (16, "cofinal_burnol_corollary"),
        (17, "nonpromotion_guards"),
        (18, "open_natural_mollifier_gate"),
    )
}

REQUIRED_NOTE = (
    "# Jensen-Window PF Jordan-Muntz/Burnol Hardy Intertwiner",
    "R_0(y)=floor(y)-y=-{y}",
    "R_omega=T_omega R_0",
    "T_omega D_d=d^omega D_d T_omega",
    "f_(epsilon,N)(t)",
    "E_(omega,N)(1/t)",
    "I-omega*H_star",
    "||H_star||_(L2(dt)->L2(dt))=2",
    "(1-2omega)||t^(-omega)f_(2omega,N)||_2",
    "M[H_star g](s)=M[g](s)/s",
    "Balazard-Saias",
    "published sequence theorem",
    "compact support",
    "Equation (JMBH.14) is open",
    "remain unproved",
)


def mobius(n: int) -> int:
    return int(sp.mobius(n))


def fractional_remainder(omega: mp.mpf, y: mp.mpf) -> mp.mpf:
    return mp.fsum(
        mp.power(m, omega)
        for m in range(1, int(mp.floor(y)) + 1)
    ) - mp.power(y, 1 + omega) / (1 + omega)


def fractional_part_integral(
    omega: mp.mpf, y: mp.mpf
) -> mp.mpf:
    total = mp.mpf("0")
    whole = int(mp.floor(y))
    for m in range(whole):
        left = mp.mpf(m)
        right = mp.mpf(m + 1)
        total += (
            (mp.power(right, omega + 1)
             - mp.power(left, omega + 1))
            / (omega + 1)
            - m
            * (mp.power(right, omega) - mp.power(left, omega))
            / omega
        )
    if y > whole:
        left = mp.mpf(whole)
        total += (
            (mp.power(y, omega + 1)
             - mp.power(left, omega + 1))
            / (omega + 1)
            - whole
            * (mp.power(y, omega) - mp.power(left, omega))
            / omega
        )
    return total


def independent_remainder_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    omega = mp.mpf(1) / 4
    for y in (
        mp.mpf("0.37"),
        mp.mpf("2.75"),
        mp.mpf("9.2"),
    ):
        direct = fractional_remainder(omega, y)
        abel = (
            -mp.power(y, omega) * mp.frac(y)
            + omega * fractional_part_integral(omega, y)
        )
        if abs(direct - abel) > mp.mpf("1e-60"):
            issues.append(
                f"fractional Abel identity failed at y={y}: "
                f"{abs(direct-abel)}"
            )
    return issues


def independent_dilation_check() -> list[str]:
    issues: list[str] = []
    x, d, omega, p = sp.symbols(
        "x d omega p", positive=True
    )
    t_power = p / (p + omega) * x ** (p + omega)
    left = d ** (-p) * t_power
    right = d**omega * t_power.subs(x, x / d)
    if sp.simplify(left - right) != 0:
        issues.append("dilation covariance failed on power test")
    return issues


def independent_intertwiner_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    omega = mp.mpf(1) / 4
    n_max = 14
    for t in (mp.mpf("0.19"), mp.mpf("0.73"), mp.mpf("1.4")):
        direct = mp.fsum(
            mobius(d)
            * mp.power(d, -omega)
            * fractional_remainder(omega, 1 / (d * t))
            for d in range(1, n_max + 1)
        )
        g_value = mp.fsum(
            mobius(d)
            * mp.power(d, -2 * omega)
            * mp.power(t, -omega)
            * mp.frac(1 / (d * t))
            for d in range(1, n_max + 1)
        )
        hardy_value = mp.fsum(
            mobius(d)
            * mp.power(d, -omega)
            * fractional_part_integral(
                omega, 1 / (d * t)
            )
            for d in range(1, n_max + 1)
        )
        intertwined = -(g_value - omega * hardy_value)
        if abs(direct - intertwined) > mp.mpf("1e-58"):
            issues.append(
                f"finite Hardy intertwiner failed at t={t}: "
                f"{abs(direct-intertwined)}"
            )
    return issues


def independent_hardy_check() -> list[str]:
    issues: list[str] = []
    a, t = sp.symbols("a t", positive=True)
    g = t ** (-a - 1)
    h_star = g / (a + 1)
    if sp.simplify(sp.diff(h_star, t) + g / t) != 0:
        issues.append("tail Hardy power model failed")

    s = sp.symbols("s", positive=True)
    # For g(t)=t^a on (0,1), H_star g=(1-t^a)/a on
    # (0,1). The two elementary Mellin integrals stay exact and
    # avoid asking a CAS to simplify an exponential-integral form.
    mellin_h = (1 / a) * (1 / s - 1 / (s + a))
    mellin_g = 1 / (s + a)
    if sp.simplify(mellin_h - mellin_g / s) != 0:
        issues.append("tail Hardy Mellin multiplier failed")
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
        "jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner"
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

    by_id = {item.get("id"): item for item in rows}
    expected_statuses = {
        "jmbh_11_hardy_norm": "source_backed",
        "jmbh_16_cofinal_burnol_corollary": (
            "source_backed_candidate"
        ),
        "jmbh_17_nonpromotion_guards": "guard_validated",
        "jmbh_18_open_natural_mollifier_gate": "open_target",
    }
    for row_id, status in expected_statuses.items():
        if by_id.get(row_id, {}).get("status") != status:
            issues.append(f"{row_id}: expected status {status}")

    expected_audit = {
        "row_count": 18,
        "exact_identity_or_corollary_count": 13,
        "classical_source_backed_count": 1,
        "source_backed_cofinal_candidate_count": 1,
        "nonpromotion_guard_count": 2,
        "open_natural_mollifier_gate_count": 1,
        "uniform_partial_norm_proved": False,
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

    issues.extend(independent_remainder_check())
    issues.extend(independent_dilation_check())
    issues.extend(independent_intertwiner_check())
    issues.extend(independent_hardy_check())

    lowered = note.lower()
    for forbidden in (
        "this proves rh",
        "the natural mollifier bound is proved",
        "hardy boundedness proves the arithmetic estimate",
        "compact support is unnecessary",
    ):
        if forbidden in lowered:
            issues.append(f"forbidden promotion language: {forbidden}")

    for issue in issues:
        print(f"JMBH [{issue}]")
    print(
        "validated Jordan-Muntz/Burnol Hardy intertwiner: "
        f"{len(rows)} rows, {len(issues)} issues, "
        "11 exact identities, 1 source-backed cofinal theorem candidate, "
        "2 nonpromotion guards, 1 open natural-mollifier gate"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
