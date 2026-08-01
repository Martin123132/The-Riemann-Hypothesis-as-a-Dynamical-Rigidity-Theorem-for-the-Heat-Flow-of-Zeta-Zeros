#!/usr/bin/env python3
"""Validate the OU-tail/Mertens mean-square reduction."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_ou_mertens_mean_square_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_ou_mertens_mean_square_reduction.md"
)

REQUIRED_IDS = {
    f"ommsr_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "definitions"),
        (2, "pnt_tail"),
        (3, "tail_difference"),
        (4, "mertens_from_tail"),
        (5, "tail_from_mertens"),
        (6, "weighted_hardy"),
        (7, "weighted_copson"),
        (8, "tail_mertens_norm_equivalence"),
        (9, "ou_weight_comparison"),
        (10, "infinite_ou_mertens_equivalence"),
        (11, "continuous_mertens_cell_identity"),
        (12, "discrete_continuous_comparison"),
        (13, "mertens_mellin_identity"),
        (14, "mellin_plancherel"),
        (15, "dyadic_block_comparison"),
        (16, "root_growth_equivalence"),
        (17, "dyadic_mean_square_rh_criterion"),
        (18, "prefix_pair_correlation"),
        (19, "cumulative_correlation"),
        (20, "anchored_correlation_criterion"),
        (21, "generic_hardy_barrier"),
        (22, "averaged_chowla_source"),
        (23, "averaged_chowla_mismatch"),
        (24, "endpoint_slice_guard"),
        (25, "open_anchored_mean_square_gate"),
    )
}

ROLE_STATUS = {
    "exact_definition": "available_exact",
    "exact_identity": "available_exact",
    "exact_inequality": "available_exact",
    "exact_equivalence": "available_exact",
    "classical_theorem_step": "source_backed",
    "literature_guard": "source_backed",
    "proof_guard": "guard_validated",
    "open_target": "open_target",
}

REQUIRED_NOTE = (
    "# Jensen-Window PF OU/Mertens Mean-Square Reduction",
    "r_k:=sum_(n>=k)mu(n)/n",
    "M(k)=sum_(n=1)^k r_n-k*r_(k+1)",
    "sqrt(2(3+alpha))/(1+alpha)",
    "sqrt(2(3-alpha))/(1-alpha)",
    "M_alpha",
    "I_alpha",
    "1/zeta(s)",
    "D_j",
    "limsup_(j->infinity)D_j^(1/j)<=1",
    "sum_(K<=k<2K)M(k)^2",
    "C_h(Y)",
    "origin-anchored integrated-correlation estimate",
    "O_alpha(N^(1-alpha))",
    "o(H^k*X)",
    "cubic scale",
    "exceptional origin slice",
    "https://arxiv.org/abs/2208.06141",
    "https://arxiv.org/abs/2508.00388",
    "https://arxiv.org/abs/1503.05121",
    "open RH-equivalent mean-square gate",
    "not a proof of RH",
)


def mobius_values(limit: int) -> list[int]:
    mu = [1] * (limit + 1)
    prime = [True] * (limit + 1)
    mu[0] = 0
    for p in range(2, limit + 1):
        if not prime[p]:
            continue
        for multiple in range(p, limit + 1, p):
            prime[multiple] = False if multiple != p else prime[multiple]
            mu[multiple] *= -1
        square = p * p
        if square <= limit:
            for multiple in range(square, limit + 1, square):
                mu[multiple] = 0
    return mu


def weighted_norm(values: list[float], alpha: float) -> float:
    return sum(
        (index ** (-alpha)) * (value * value)
        for index, value in enumerate(values, start=1)
    )


def check_schema(payload: dict, note: str, issues: list[str]) -> None:
    if payload.get("kind") != "jensen_window_pf_ou_mertens_mean_square_reduction":
        issues.append("bad kind")
    rows = payload.get("rows", [])
    if len(rows) != 25:
        issues.append(f"expected 25 rows, found {len(rows)}")
    ids = {item.get("id") for item in rows}
    if ids != REQUIRED_IDS:
        issues.append(
            "row id mismatch: "
            f"missing={sorted(REQUIRED_IDS-ids)}, extra={sorted(ids-REQUIRED_IDS)}"
        )
    for item in rows:
        role = item.get("role")
        expected_status = ROLE_STATUS.get(role)
        if expected_status is None:
            issues.append(f"{item.get('id')}: bad role {role!r}")
        elif item.get("status") != expected_status:
            issues.append(
                f"{item.get('id')}: role/status mismatch "
                f"{role!r}/{item.get('status')!r}"
            )
        for field in ("statement", "proof_boundary"):
            if not item.get(field):
                issues.append(f"{item.get('id')}: missing {field}")

    audit = payload.get("audit", {})
    expected_audit = {
        "row_count": 25,
        "exact_reduction_count": 16,
        "classical_theorem_step_count": 4,
        "proof_guard_count": 4,
        "literature_guard_count": 1,
        "open_mean_square_gate_count": 1,
        "dyadic_mean_square_rh_equivalence_proved": True,
        "averaged_chowla_closes_gate": False,
        "ou_uniform_bound_proved": False,
        "full_burnol_bound_proved": False,
        "rh_proved": False,
        "pf_infinity_proved": False,
        "lambda_le_zero_proved": False,
    }
    for key, value in expected_audit.items():
        if audit.get(key) != value:
            issues.append(
                f"audit {key}: expected {value!r}, found {audit.get(key)!r}"
            )

    for needle in REQUIRED_NOTE:
        if needle not in note:
            issues.append(f"note missing {needle!r}")


def check_finite_identities(issues: list[str]) -> None:
    limit = 64
    mu = mobius_values(limit)
    mertens = [0] * (limit + 1)
    for n in range(1, limit + 1):
        mertens[n] = mertens[n - 1] + mu[n]

    tails = [Fraction(0) for _ in range(limit + 2)]
    for k in range(limit, 0, -1):
        tails[k] = tails[k + 1] + Fraction(mu[k], k)

    for k in range(1, limit + 1):
        if tails[k] - tails[k + 1] != Fraction(mu[k], k):
            issues.append(f"tail difference failed at k={k}")
        reconstructed = sum(tails[1 : k + 1], Fraction(0)) - k * tails[k + 1]
        if reconstructed != mertens[k]:
            issues.append(f"Mertens-from-tail failed at k={k}")

    for k in range(1, limit):
        finite_abel = (
            Fraction(mertens[limit], limit)
            - Fraction(mertens[k], k + 1)
            + sum(
                (
                    Fraction(mertens[n], n * (n + 1))
                    for n in range(k + 1, limit)
                ),
                Fraction(0),
            )
        )
        if finite_abel != tails[k + 1]:
            issues.append(f"finite tail-from-Mertens failed at k={k}")

    mellin_left = sum(
        (Fraction(mu[n], n * n) for n in range(1, limit + 1)),
        Fraction(0),
    )
    mellin_right = Fraction(mertens[limit], limit * limit) + sum(
        (
            mertens[k]
            * (Fraction(1, k * k) - Fraction(1, (k + 1) * (k + 1)))
            for k in range(1, limit)
        ),
        Fraction(0),
    )
    if mellin_left != mellin_right:
        issues.append("finite Mertens Mellin/Abel identity failed")

    prefix_square = sum(mertens[k] ** 2 for k in range(1, limit + 1))
    diagonal = sum(mu[n] ** 2 * (limit - n + 1) for n in range(1, limit + 1))
    off_diagonal = 2 * sum(
        mu[m] * mu[m + h] * (limit - m - h + 1)
        for h in range(1, limit)
        for m in range(1, limit - h + 1)
    )
    if prefix_square != diagonal + off_diagonal:
        issues.append("prefix pair-correlation expansion failed")

    cumulative = 0
    for h in range(1, limit):
        running = 0
        for y in range(1, limit - h + 1):
            running += mu[y] * mu[y + h]
            cumulative += running
    if off_diagonal != 2 * cumulative:
        issues.append("cumulative correlation reindexing failed")


def check_weighted_inequalities(issues: list[str]) -> None:
    for alpha in (0.2, 0.5, 0.8):
        for k in range(1, 257):
            weight = (k ** (1.0 - alpha)) - ((k - 1) ** (1.0 - alpha))
            lower = (1.0 - alpha) * (k ** (-alpha))
            upper = k ** (-alpha)
            if weight < lower * (1.0 - 2e-13) or weight > upper * (1.0 + 2e-13):
                issues.append(f"OU weight comparison failed at alpha={alpha}, k={k}")
                break

        size = 96
        vectors = [
            [1.0 if n == 1 else 0.0 for n in range(1, size + 1)],
            [((-1.0) ** n) / math.sqrt(n) for n in range(1, size + 1)],
            [((17 * n) % 23 - 11) / 11.0 for n in range(1, size + 1)],
        ]
        hardy_bound = math.sqrt(2.0 * (3.0 + alpha)) / (1.0 + alpha)
        copson_bound = math.sqrt(2.0 * (3.0 - alpha)) / (1.0 - alpha)
        for vector in vectors:
            hardy = []
            running = 0.0
            for index, value in enumerate(vector, start=1):
                running += value
                hardy.append(running / index)
            copson = [0.0] * size
            running = 0.0
            for index in range(size, 0, -1):
                running += vector[index - 1] / index
                copson[index - 1] = running
            source_norm = weighted_norm(vector, alpha)
            if weighted_norm(hardy, alpha) > (
                hardy_bound * hardy_bound * source_norm * (1.0 + 2e-12)
            ):
                issues.append(f"weighted Hardy check failed at alpha={alpha}")
            if weighted_norm(copson, alpha) > (
                copson_bound * copson_bound * source_norm * (1.0 + 2e-12)
            ):
                issues.append(f"weighted Copson check failed at alpha={alpha}")

    mu = mobius_values(64)
    mertens = [0] * 65
    for n in range(1, 65):
        mertens[n] = mertens[n - 1] + mu[n]
    alpha = 0.37
    m_alpha = sum(
        mertens[k] ** 2 * (k ** (-2.0 - alpha))
        for k in range(1, 65)
    )
    dyadic = 0.0
    for j in range(0, 6):
        start = 2**j
        stop = 2 ** (j + 1)
        d_j = (2.0 ** (-2 * j)) * sum(
            mertens[k] ** 2 for k in range(start, stop)
        )
        dyadic += (2.0 ** (-alpha * j)) * d_j
    if m_alpha > dyadic * (1.0 + 2e-13):
        issues.append("dyadic upper comparison failed")
    if m_alpha < (2.0 ** (-2.0 - alpha)) * dyadic * (1.0 - 2e-13):
        issues.append("dyadic lower comparison failed")

    i_alpha = sum(
        (mertens[k] ** 2)
        * (
            (k ** (-1.0 - alpha) - (k + 1) ** (-1.0 - alpha))
            / (1.0 + alpha)
        )
        for k in range(1, 65)
    )
    if i_alpha > m_alpha * (1.0 + 2e-13):
        issues.append("continuous/discrete upper comparison failed")
    if i_alpha < (2.0 ** (-2.0 - alpha)) * m_alpha * (1.0 - 2e-13):
        issues.append("continuous/discrete lower comparison failed")


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
    check_schema(payload, note, issues)
    check_finite_identities(issues)
    check_weighted_inequalities(issues)
    for issue in issues:
        print(f"OU-MERTENS {issue}")
    print(
        "validated OU/Mertens mean-square reduction: "
        f"{len(payload.get('rows', []))} rows, {len(issues)} issues, "
        f"{payload.get('audit', {}).get('exact_reduction_count', 0)} "
        "exact reductions, "
        f"{payload.get('audit', {}).get('proof_guard_count', 0)} "
        "proof guards, "
        f"{payload.get('audit', {}).get('open_mean_square_gate_count', 0)} "
        "open anchored mean-square gate"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
