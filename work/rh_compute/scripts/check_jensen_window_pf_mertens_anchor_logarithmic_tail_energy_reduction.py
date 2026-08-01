#!/usr/bin/env python3
"""Validate the Mertens anchor/logarithmic tail-energy reduction."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.md"
)

REQUIRED_IDS = {
    f"malter_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "definitions"),
        (2, "adjacent_tail_difference"),
        (3, "prefix_from_all_tails"),
        (4, "normalized_hardy_identity"),
        (5, "tail_from_prefix"),
        (6, "weighted_hardy"),
        (7, "weighted_copson"),
        (8, "prefix_tail_norm_equivalence"),
        (9, "stable_prefix_energy"),
        (10, "power_sum_comparison"),
        (11, "logarithmic_stable_prefix_equivalence"),
        (12, "three_way_energy_equivalence"),
        (13, "cofinal_rh_criterion"),
        (14, "all_cutoff_anchor_resolution"),
        (15, "fixed_cutoff_gauge_guard"),
        (16, "critical_pointwise_consequence"),
        (17, "uniform_prefix_one_log_short"),
        (18, "critical_scalar_model"),
        (19, "scalar_coefficient_envelope"),
        (20, "scalar_prefix_energy_bound"),
        (21, "scalar_logarithmic_divergence"),
        (22, "countermodel_scope_guard"),
        (23, "full_q_distinction_guard"),
        (24, "literature_context"),
        (25, "focused_source_search_guard"),
        (26, "open_logarithmic_gain_gate"),
    )
}

ROLE_STATUS = {
    "exact_definition": "available_exact",
    "exact_identity": "available_exact",
    "exact_inequality": "available_exact",
    "exact_equivalence": "available_exact",
    "exact_consequence": "available_exact",
    "classical_theorem_step": "source_backed",
    "literature_guard": "source_backed",
    "proof_guard": "guard_validated",
    "open_target": "open_target",
}

REQUIRED_NOTE = (
    "# Jensen-Window PF Mertens Anchor/Logarithmic Tail-Energy Reduction",
    "r_(N-1)-r_N=a_N/N",
    "sum_(j=0)^(N-1)r_j-N*r_N",
    "Weighted Hardy/Copson Equivalence",
    "E_alpha<infinity",
    "sum_(N>=1)P_(alpha/2,N)/N<infinity",
    "one logarithm weaker",
    "RH-equivalent target",
    "The Exact Logarithmic Barrier",
    "gamma:=(1+alpha)/2",
    "2^(-(1+alpha))/(1+alpha)<=P_N<=1",
    "scalar/operator countermodel",
    "one cutoff:",
    "full Burnol energy `Q`",
    "not a novelty or priority claim",
    "one open RH-equivalent arithmetic gate",
    "not a proof of RH",
)


def close(
    left: complex | float,
    right: complex | float,
    tolerance: float = 3e-11,
) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def mobius_values(limit: int) -> list[int]:
    mu = [0] * (limit + 1)
    composite = [False] * (limit + 1)
    primes: list[int] = []
    mu[1] = 1
    for n in range(2, limit + 1):
        if not composite[n]:
            primes.append(n)
            mu[n] = -1
        for prime in primes:
            product = n * prime
            if product > limit:
                break
            composite[product] = True
            if n % prime == 0:
                mu[product] = 0
                break
            mu[product] = -mu[n]
    return mu


def check_schema(payload: dict, note: str, issues: list[str]) -> None:
    expected_kind = (
        "jensen_window_pf_mertens_anchor_"
        "logarithmic_tail_energy_reduction"
    )
    if payload.get("kind") != expected_kind:
        issues.append("bad kind")
    rows = payload.get("rows", [])
    if len(rows) != 26:
        issues.append(f"expected 26 rows, found {len(rows)}")
    ids = {item.get("id") for item in rows}
    if ids != REQUIRED_IDS:
        issues.append(
            "row id mismatch: "
            f"missing={sorted(REQUIRED_IDS-ids)}, "
            f"extra={sorted(ids-REQUIRED_IDS)}"
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

    expected_audit = {
        "row_count": 26,
        "exact_reduction_count": 17,
        "classical_theorem_step_count": 2,
        "literature_guard_count": 2,
        "proof_guard_count": 4,
        "open_logarithmic_gain_gate_count": 1,
        "prefix_tail_norm_equivalence_proved": True,
        "logarithmic_stable_prefix_equivalence_proved": True,
        "cofinal_rh_equivalence_proved": True,
        "uniform_stable_prefix_closes_log_gate": False,
        "logarithmic_gain_proved": False,
        "full_burnol_bound_proved": False,
        "rh_proved": False,
        "pf_infinity_proved": False,
        "lambda_le_zero_proved": False,
    }
    audit = payload.get("audit", {})
    for key, value in expected_audit.items():
        if audit.get(key) != value:
            issues.append(
                f"audit {key}: expected {value!r}, "
                f"found {audit.get(key)!r}"
            )
    for needle in REQUIRED_NOTE:
        if needle not in note:
            issues.append(f"note missing {needle!r}")


def check_exact_finite_algebra(issues: list[str]) -> None:
    limit = 37
    tails = [
        Fraction(((19 * n + 7) % 31) - 15, (n + 2) * 17)
        for n in range(limit)
    ] + [Fraction(0)]
    coefficients = [Fraction(0)] * (limit + 1)
    prefix = [Fraction(0)] * (limit + 1)

    for n in range(1, limit + 1):
        coefficients[n] = n * (tails[n - 1] - tails[n])
        prefix[n] = prefix[n - 1] + coefficients[n]
        if tails[n - 1] - tails[n] != coefficients[n] / n:
            issues.append(f"rational adjacent-tail identity failed at n={n}")
        reconstructed = sum(tails[:n]) - n * tails[n]
        if prefix[n] != reconstructed:
            issues.append(f"rational prefix reconstruction failed at n={n}")
        normalized = sum(tails[:n], Fraction(0)) / n - tails[n]
        if prefix[n] / n != normalized:
            issues.append(f"rational Hardy identity failed at n={n}")

    for anchor in (1, 4, 13, 29):
        direct_tail = sum(
            coefficients[n] / n
            for n in range(anchor + 1, limit + 1)
        )
        if direct_tail != tails[anchor]:
            issues.append(
                f"rational coefficient-tail reconstruction failed at "
                f"N={anchor}"
            )
        inverse = (
            -prefix[anchor] / (anchor + 1)
            + sum(
                prefix[n] / (n * (n + 1))
                for n in range(anchor + 1, limit)
            )
            + prefix[limit] / limit
        )
        if inverse != tails[anchor]:
            issues.append(f"rational inverse Abel identity failed at N={anchor}")


def check_mobius_tail_and_prefix(issues: list[str]) -> None:
    limit = 320
    mu = mobius_values(limit)
    for alpha in (0.17, 0.48, 0.83):
        coefficients = [0.0] * (limit + 1)
        prefix = [0.0] * (limit + 1)
        tails = [0.0] * (limit + 1)
        for n in range(1, limit + 1):
            coefficients[n] = mu[n] * (n ** (-alpha))
            prefix[n] = prefix[n - 1] + coefficients[n]
        for n in range(limit - 1, -1, -1):
            tails[n] = tails[n + 1] + coefficients[n + 1] / (n + 1)

        for n in range(1, 151):
            if not close(
                tails[n - 1] - tails[n],
                coefficients[n] / n,
            ):
                issues.append(
                    f"Mobius adjacent tail failed at alpha={alpha}, n={n}"
                )
                break
            reconstructed = sum(tails[:n]) - n * tails[n]
            if not close(prefix[n], reconstructed, tolerance=8e-11):
                issues.append(
                    f"Mobius prefix reconstruction failed at "
                    f"alpha={alpha}, n={n}"
                )
                break

        for anchor in (1, 7, 31, 127):
            inverse = (
                -prefix[anchor] / (anchor + 1)
                + sum(
                    prefix[n] / (n * (n + 1))
                    for n in range(anchor + 1, limit)
                )
                + prefix[limit] / limit
            )
            if not close(inverse, tails[anchor], tolerance=8e-11):
                issues.append(
                    f"Mobius finite inverse tail failed at "
                    f"alpha={alpha}, N={anchor}"
                )

        for n in (1, 2, 5, 17, 63, 149):
            power_sum = sum(k**alpha for k in range(1, n + 1))
            p_value = tails[n] ** 2 * power_sum
            if not close(p_value, tails[n] ** 2 * power_sum):
                issues.append("stable-prefix energy arithmetic failed")
            scaled = p_value / n
            lower = (n**alpha) * tails[n] ** 2 / (1.0 + alpha)
            upper = (n**alpha) * tails[n] ** 2
            if scaled < lower * (1.0 - 4e-13):
                issues.append(
                    f"stable-prefix lower comparison failed at "
                    f"alpha={alpha}, N={n}"
                )
            if scaled > upper * (1.0 + 4e-13):
                issues.append(
                    f"stable-prefix upper comparison failed at "
                    f"alpha={alpha}, N={n}"
                )


def check_weighted_operators(issues: list[str]) -> None:
    limit = 600
    for alpha in (0.12, 0.51, 0.86):
        vector = [0.0] + [
            ((((37 * n + 11) % 53) - 26) / 27.0)
            * math.exp(-n / 83.0)
            for n in range(1, limit + 1)
        ]
        source_norm_sq = sum(
            (n**alpha) * vector[n] ** 2 for n in range(1, limit + 1)
        )

        running = 0.0
        hardy_norm_sq = 0.0
        for n in range(1, limit + 1):
            running += vector[n]
            hardy_norm_sq += (n**alpha) * (running / n) ** 2
        h_alpha = math.sqrt(2.0 * (3.0 - alpha)) / (1.0 - alpha)
        if hardy_norm_sq > h_alpha * h_alpha * source_norm_sq * (
            1.0 + 2e-12
        ):
            issues.append(f"finite weighted Hardy test failed at alpha={alpha}")

        running = 0.0
        copson = [0.0] * (limit + 1)
        for n in range(limit, 0, -1):
            running += vector[n] / n
            copson[n] = running
        copson_norm_sq = sum(
            (n**alpha) * copson[n] ** 2 for n in range(1, limit + 1)
        )
        c_alpha = math.sqrt(2.0 * (3.0 + alpha)) / (1.0 + alpha)
        if copson_norm_sq > c_alpha * c_alpha * source_norm_sq * (
            1.0 + 2e-12
        ):
            issues.append(f"finite weighted Copson test failed at alpha={alpha}")

        shifted_tail = [0.0] * (limit + 1)
        running = 0.0
        for n in range(limit, 0, -1):
            shifted_tail[n - 1] = running
            running += vector[n] / (n + 1)
        shifted_norm_sq = sum(
            (n**alpha) * shifted_tail[n] ** 2
            for n in range(1, limit)
        )
        if shifted_norm_sq > c_alpha * c_alpha * source_norm_sq * (
            1.0 + 2e-12
        ):
            issues.append(
                f"shifted Copson domination failed at alpha={alpha}"
            )


def check_critical_scalar_model(issues: list[str]) -> None:
    limit = 65536
    for alpha in (0.19, 0.5, 0.81):
        gamma = (1.0 + alpha) / 2.0
        prefix = 0.0
        harmonic_energy = 0.0
        log_energy = 0.0
        p_at: dict[int, float] = {}
        prefix_at: dict[int, float] = {}

        for n in range(1, limit + 1):
            r_previous = n ** (-gamma)
            r_current = (n + 1) ** (-gamma)
            coefficient = n * (r_previous - r_current)
            if coefficient <= 0.0:
                issues.append(
                    f"scalar coefficient is nonpositive at alpha={alpha}, n={n}"
                )
                break
            if coefficient > gamma * (n ** (-gamma)) * (1.0 + 3e-13):
                issues.append(
                    f"scalar gamma envelope failed at alpha={alpha}, n={n}"
                )
                break
            if coefficient > gamma * (n ** (-alpha)) * (1.0 + 3e-13):
                issues.append(
                    f"scalar alpha envelope failed at alpha={alpha}, n={n}"
                )
                break

            prefix += coefficient
            harmonic_energy += n**alpha
            p_value = r_current * r_current * harmonic_energy
            log_energy += p_value / n
            if n in (256, 4096, limit):
                p_at[n] = p_value
                prefix_at[n] = prefix

            lower = 2.0 ** (-(1.0 + alpha)) / (1.0 + alpha)
            if p_value < lower * (1.0 - 4e-13) or p_value > 1.0 + 4e-13:
                issues.append(
                    f"scalar stable-prefix bound failed at "
                    f"alpha={alpha}, n={n}, P={p_value}"
                )
                break

        asymptotic_scale = (
            gamma / (1.0 - gamma)
        ) * (limit ** (1.0 - gamma))
        ratio = prefix / asymptotic_scale
        early_scale = (
            gamma / (1.0 - gamma)
        ) * (4096 ** (1.0 - gamma))
        early_ratio = prefix_at.get(4096, 0.0) / early_scale
        if not 0.0 < early_ratio < ratio <= 1.05:
            issues.append(
                f"scalar prefix asymptotic test failed at alpha={alpha}: "
                f"early={early_ratio}, late={ratio}"
            )

        last_term_ratio = (
            (limit ** (alpha - 2.0))
            * prefix
            * prefix
            * limit
            / ((gamma / (1.0 - gamma)) ** 2)
        )
        if not close(last_term_ratio, ratio * ratio, tolerance=2e-12):
            issues.append(
                f"scalar harmonic boundary test failed at alpha={alpha}: "
                f"ratio={last_term_ratio}"
            )

        if p_at.get(limit, 0.0) <= 0.0:
            issues.append(f"scalar P limit test missing at alpha={alpha}")
        if log_energy <= p_at.get(256, 0.0) * math.log(limit / 256):
            issues.append(
                f"scalar logarithmic growth test failed at alpha={alpha}"
            )

        tail_anchor = 37
        tail_reconstruction = sum(
            n
            * (n ** (-gamma) - (n + 1) ** (-gamma))
            / n
            for n in range(tail_anchor + 1, limit + 1)
        ) + (limit + 1) ** (-gamma)
        if not close(
            tail_reconstruction,
            (tail_anchor + 1) ** (-gamma),
            tolerance=2e-10,
        ):
            issues.append(
                f"scalar tail telescoping failed at alpha={alpha}"
            )


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
    check_exact_finite_algebra(issues)
    check_mobius_tail_and_prefix(issues)
    check_weighted_operators(issues)
    check_critical_scalar_model(issues)
    for issue in issues:
        print(f"MERTENS-ANCHOR-LOG-TAIL {issue}")
    audit = payload.get("audit", {})
    print(
        "validated Mertens anchor/logarithmic tail-energy reduction: "
        f"{len(payload.get('rows', []))} rows, {len(issues)} issues, "
        f"{audit.get('exact_reduction_count', 0)} exact reductions, "
        f"{audit.get('proof_guard_count', 0)} proof guards, "
        f"{audit.get('open_logarithmic_gain_gate_count', 0)} open "
        "logarithmic gate"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
