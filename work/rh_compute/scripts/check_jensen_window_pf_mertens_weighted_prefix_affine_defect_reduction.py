#!/usr/bin/env python3
"""Validate the Mertens weighted-prefix and affine-defect reduction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.md"
)

REQUIRED_IDS = {
    f"mwpad_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "weighted_prefix_definitions"),
        (2, "forward_abel"),
        (3, "inverse_abel"),
        (4, "logarithmic_coordinates"),
        (5, "forward_volterra"),
        (6, "inverse_volterra"),
        (7, "sharp_norm_sandwich"),
        (8, "weighted_prefix_mellin"),
        (9, "mellin_plancherel"),
        (10, "shared_reciprocal_zeta_spectrum"),
        (11, "weighted_prefix_cell_identity"),
        (12, "weighted_prefix_cell_comparison"),
        (13, "mertens_cell_comparison"),
        (14, "discrete_energy_equivalence"),
        (15, "weighted_prefix_rh_criterion"),
        (16, "burnol_coefficient_match"),
        (17, "first_block_floor_collapse"),
        (18, "first_block_affine_defect"),
        (19, "weighted_tail_abel"),
        (20, "residual_weight_comparison"),
        (21, "affine_gram"),
        (22, "three_term_block_bound"),
        (23, "anchor_gauge_invariance"),
        (24, "constant_tail_countermodel"),
        (25, "tail_rate_does_not_fix_anchor"),
        (26, "correlation_one_dimensional_collapse"),
        (27, "discrete_energy_difference"),
        (28, "vaughan_source"),
        (29, "vaughan_type_i_ii_identity"),
        (30, "endogenous_test_guard"),
        (31, "cauchy_schwarz_recycling_guard"),
        (32, "precollapse_bounded_shift_route"),
        (33, "open_anchor_gate"),
        (34, "open_type_i_ii_gate"),
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
    "# Jensen-Window PF Mertens Weighted-Prefix/Affine-Defect Reduction",
    "A_alpha(x):=sum_(n<=x)mu(n)n^(-alpha)",
    "a=m+alpha*V_lambda*m",
    "((1+alpha)/(1-alpha))^2",
    "1/[s*zeta(s+alpha)]",
    "same reciprocal-zeta boundary data",
    "E_alpha",
    "H_(omega,N)(k)=A_alpha(k)-A_alpha(N)",
    "unobserved anchor mode",
    "constant `C` from `N` onward",
    "sum_(n=2)^X (X-n+1)mu(n)M(n-1)",
    "Green and Tao's Mobius Vaughan identity",
    "f_X(n):=(X-n+1)M(n-1)",
    "feeds the target energy back",
    "f_(h,Y)(m)=mu(m+h)1_(m<=Y)",
    "B_1(X)-TI_X+TII_X",
    "concrete signed estimate",
    "https://doi.org/10.5802/aif.2401",
    "two open arithmetic gates",
    "not a proof of",
)


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


def close(left: complex | float, right: complex | float, tolerance: float = 2e-11) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def check_schema(payload: dict, note: str, issues: list[str]) -> None:
    expected_kind = (
        "jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction"
    )
    if payload.get("kind") != expected_kind:
        issues.append("bad kind")
    rows = payload.get("rows", [])
    if len(rows) != 34:
        issues.append(f"expected 34 rows, found {len(rows)}")
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

    expected_audit = {
        "row_count": 34,
        "exact_reduction_count": 24,
        "classical_theorem_step_count": 1,
        "literature_guard_count": 1,
        "proof_guard_count": 4,
        "open_handoff_gate_count": 2,
        "weighted_prefix_rh_equivalence_proved": True,
        "first_block_anchor_defect_proved": True,
        "vaughan_identity_closes_gate": False,
        "anchor_gate_proved": False,
        "signed_type_i_ii_gate_proved": False,
        "rh_proved": False,
        "pf_infinity_proved": False,
        "lambda_le_zero_proved": False,
    }
    audit = payload.get("audit", {})
    for key, value in expected_audit.items():
        if audit.get(key) != value:
            issues.append(
                f"audit {key}: expected {value!r}, found {audit.get(key)!r}"
            )
    for needle in REQUIRED_NOTE:
        if needle not in note:
            issues.append(f"note missing {needle!r}")


def check_abel_and_energy(issues: list[str]) -> None:
    limit = 72
    mu = mobius_values(limit)
    for alpha in (0.2, 0.47, 0.8):
        mertens = [0.0] * (limit + 1)
        weighted = [0.0] * (limit + 1)
        for n in range(1, limit + 1):
            mertens[n] = mertens[n - 1] + mu[n]
            weighted[n] = weighted[n - 1] + mu[n] * (n ** (-alpha))

        for k in range(1, limit + 1):
            forward = mertens[k] * (k ** (-alpha)) + sum(
                mertens[n]
                * ((n ** (-alpha)) - ((n + 1) ** (-alpha)))
                for n in range(1, k)
            )
            inverse = (k**alpha) * weighted[k] - sum(
                weighted[n] * (((n + 1) ** alpha) - (n**alpha))
                for n in range(1, k)
            )
            if not close(forward, weighted[k]):
                issues.append(f"forward Abel failed at alpha={alpha}, k={k}")
                break
            if not close(inverse, mertens[k]):
                issues.append(f"inverse Abel failed at alpha={alpha}, k={k}")
                break

        i_energy = sum(
            (mertens[k] ** 2)
            * (
                (k ** (-1.0 - alpha) - (k + 1) ** (-1.0 - alpha))
                / (1.0 + alpha)
            )
            for k in range(1, limit)
        ) + (mertens[limit] ** 2) * (
            limit ** (-1.0 - alpha)
        ) / (1.0 + alpha)
        j_energy = sum(
            (weighted[k] ** 2)
            * (
                (k ** (alpha - 1.0) - (k + 1) ** (alpha - 1.0))
                / (1.0 - alpha)
            )
            for k in range(1, limit)
        ) + (weighted[limit] ** 2) * (
            limit ** (alpha - 1.0)
        ) / (1.0 - alpha)
        ratio = ((1.0 + alpha) / (1.0 - alpha)) ** 2
        if j_energy < i_energy * (1.0 - 3e-11):
            issues.append(f"continuous energy lower bound failed at alpha={alpha}")
        if j_energy > ratio * i_energy * (1.0 + 3e-11):
            issues.append(f"continuous energy upper bound failed at alpha={alpha}")

        m_discrete = sum(
            (mertens[k] ** 2) * (k ** (-2.0 - alpha))
            for k in range(1, limit + 1)
        )
        a_discrete = sum(
            (weighted[k] ** 2) * (k ** (alpha - 2.0))
            for k in range(1, limit + 1)
        )
        if i_energy > m_discrete * (1.0 + 3e-11):
            issues.append(f"Mertens cell upper bound failed at alpha={alpha}")
        if j_energy > a_discrete * (1.0 + 3e-11):
            issues.append(f"weighted-prefix cell upper bound failed at alpha={alpha}")

        for spectral_t in (0.0, 0.5, 3.0, 19.0):
            beta = (1.0 + alpha) / 2.0
            lam = (1.0 - alpha) / 2.0
            multiplier_sq = (
                (beta * beta + spectral_t * spectral_t)
                / (lam * lam + spectral_t * spectral_t)
            )
            if multiplier_sq < 1.0 - 2e-14 or multiplier_sq > ratio * (
                1.0 + 2e-14
            ):
                issues.append(
                    f"Hardy multiplier range failed at alpha={alpha}, "
                    f"t={spectral_t}"
                )

        s = complex(1.3, 0.7)
        mellin_cells = sum(
            weighted[k]
            * ((k ** (-s)) - ((k + 1) ** (-s)))
            / s
            for k in range(1, limit)
        ) + weighted[limit] * (limit ** (-s)) / s
        mellin_series = sum(
            mu[n] * (n ** (-alpha - s)) for n in range(1, limit + 1)
        ) / s
        if not close(mellin_cells, mellin_series, tolerance=5e-11):
            issues.append(f"weighted-prefix Mellin identity failed at alpha={alpha}")


def check_first_block(issues: list[str]) -> None:
    alpha = 0.42
    limit = 128
    n_anchor = 24
    mu = mobius_values(limit)
    coefficients = [0.0] * (limit + 1)
    prefix = [0.0] * (limit + 1)
    for n in range(1, limit + 1):
        coefficients[n] = mu[n] * (n ** (-alpha))
        prefix[n] = prefix[n - 1] + coefficients[n]

    r_tail = sum(
        coefficients[n] / n for n in range(n_anchor + 1, limit + 1)
    )
    tail_abel = (
        -prefix[n_anchor] / (n_anchor + 1)
        + sum(
            prefix[n] / (n * (n + 1))
            for n in range(n_anchor + 1, limit)
        )
        + prefix[limit] / limit
    )
    if not close(r_tail, tail_abel):
        issues.append("finite weighted-tail Abel identity failed")

    residuals: list[tuple[int, float]] = []
    for k in range(n_anchor + 1, 2 * n_anchor + 1):
        h_value = sum(
            coefficients[d] * (k // d) for d in range(n_anchor + 1, k + 1)
        )
        increment = prefix[k] - prefix[n_anchor]
        if not close(h_value, increment):
            issues.append(f"first-block floor collapse failed at k={k}")
        d_value = h_value - k * r_tail
        affine_value = prefix[k] - prefix[n_anchor] - k * r_tail
        if not close(d_value, affine_value):
            issues.append(f"first-block affine identity failed at k={k}")
        residuals.append((k, d_value))

    v_value = sum(value * value for _, value in residuals)
    e_value = sum(
        (k ** (alpha - 2.0)) * value * value for k, value in residuals
    )
    lower = (2.0 ** (alpha - 2.0)) * (
        n_anchor ** (alpha - 2.0)
    ) * v_value
    upper = (n_anchor ** (alpha - 2.0)) * v_value
    if e_value < lower * (1.0 - 3e-12) or e_value > upper * (1.0 + 3e-12):
        issues.append("first-block residual weight comparison failed")

    weights = [
        (k, k ** (alpha - 2.0))
        for k in range(n_anchor + 1, 2 * n_anchor + 1)
    ]
    w0 = sum(weight for _, weight in weights)
    w1 = sum(weight * k for k, weight in weights)
    w2 = sum(weight * k * k for k, weight in weights)
    if w0 * w2 - w1 * w1 <= 0.0:
        issues.append("affine Gram determinant is not positive")
    line_direct = sum(
        weight * (prefix[n_anchor] + k * r_tail) ** 2
        for k, weight in weights
    )
    line_gram = (
        prefix[n_anchor] ** 2 * w0
        + 2.0 * prefix[n_anchor] * r_tail * w1
        + r_tail * r_tail * w2
    )
    if not close(line_direct, line_gram):
        issues.append("affine Gram expansion failed")

    f_block = sum(
        weight * prefix[k] ** 2 for k, weight in weights
    )
    coarse_rhs = (
        (n_anchor ** ((alpha - 2.0) / 2.0)) * math.sqrt(v_value)
        + (n_anchor ** ((alpha - 1.0) / 2.0)) * abs(prefix[n_anchor])
        + (2.0 ** (alpha / 2.0))
        * (n_anchor ** ((1.0 + alpha) / 2.0))
        * abs(r_tail)
    )
    if math.sqrt(f_block) > coarse_rhs * (1.0 + 3e-12):
        issues.append("three-term first-block bound failed")

    constant = 7.0
    constant_tail_sum = (
        -constant / (n_anchor + 1)
        + sum(
            constant / (n * (n + 1))
            for n in range(n_anchor + 1, 20 * n_anchor)
        )
        + constant / (20 * n_anchor)
    )
    if not close(constant_tail_sum, 0.0):
        issues.append("constant-tail Abel countermodel failed")
    constant_energy = constant * constant * w0
    if constant_energy <= 0.0:
        issues.append("constant-anchor countermodel has no positive energy")


def check_correlation_and_vaughan(issues: list[str]) -> None:
    x_limit = 80
    mu = mobius_values(2 * x_limit)
    mertens = [0] * (2 * x_limit + 1)
    for n in range(1, 2 * x_limit + 1):
        mertens[n] = mertens[n - 1] + mu[n]

    cumulative = 0
    for h in range(1, x_limit):
        running = 0
        for y in range(1, x_limit - h + 1):
            running += mu[y] * mu[y + h]
            cumulative += running
    collapsed = sum(
        (x_limit - n + 1) * mu[n] * mertens[n - 1]
        for n in range(2, x_limit + 1)
    )
    if cumulative != collapsed:
        issues.append("one-dimensional correlation collapse failed")

    prefix_energy = sum(mertens[k] ** 2 for k in range(1, x_limit + 1))
    energy_difference = sum(
        (x_limit - n + 1)
        * (
            2 * mu[n] * mertens[n - 1]
            + mu[n] * mu[n]
        )
        for n in range(1, x_limit + 1)
    )
    if prefix_energy != energy_difference:
        issues.append("discrete energy-difference identity failed")

    endogenous_norm = sum(
        ((x_limit - n + 1) * mertens[n - 1]) ** 2
        for n in range(1, x_limit + 1)
    )
    if endogenous_norm > (x_limit**2) * prefix_energy:
        issues.append("endogenous test L2 guard failed")

    block_n = 30
    u_cut = 3
    v_cut = 4
    test = [((17 * n + 5) % 23) - 11 for n in range(2 * block_n + 1)]
    lhs = sum(
        mu[n] * test[n] for n in range(block_n + 1, 2 * block_n + 1)
    )

    type_i = 0
    for d in range(1, u_cut * v_cut + 1):
        a_d = sum(
            mu[b] * mu[c]
            for b in range(1, u_cut + 1)
            for c in range(1, v_cut + 1)
            if b * c == d
        )
        inner = sum(
            test[d * w]
            for w in range(block_n // d + 1, (2 * block_n) // d + 1)
        )
        type_i += a_d * inner

    type_ii = 0
    for d in range(v_cut + 1, (2 * block_n) // u_cut + 1):
        b_d = sum(
            mu[c] for c in range(v_cut + 1, d + 1) if d % c == 0
        )
        lower = max(u_cut, block_n // d) + 1
        upper = (2 * block_n) // d
        inner = sum(mu[w] * test[d * w] for w in range(lower, upper + 1))
        type_ii += b_d * inner
    if lhs != -type_i + type_ii:
        issues.append(
            "finite Vaughan Type I/II identity failed: "
            f"lhs={lhs}, rhs={-type_i + type_ii}"
        )

    aggregate_x = 24
    aggregate_target = 0
    for h in range(1, aggregate_x):
        running = 0
        for y in range(1, aggregate_x - h + 1):
            running += mu[y] * mu[y + h]
            aggregate_target += running
    boundary = sum(
        mu[n] * (aggregate_x - n + 1) for n in range(2, aggregate_x + 1)
    )
    aggregate_type_i = 0
    aggregate_type_ii = 0
    dyadic_n = 1
    while dyadic_n < aggregate_x:
        u_cut = math.isqrt(dyadic_n)
        v_cut = u_cut
        a_coefficients = {
            d: sum(
                mu[b] * mu[c]
                for b in range(1, u_cut + 1)
                for c in range(1, v_cut + 1)
                if b * c == d
            )
            for d in range(1, u_cut * v_cut + 1)
        }
        b_coefficients = {
            d: sum(
                mu[c] for c in range(v_cut + 1, d + 1) if d % c == 0
            )
            for d in range(v_cut + 1, (2 * dyadic_n) // u_cut + 1)
        }
        for h in range(1, aggregate_x):
            for y in range(1, aggregate_x - h + 1):
                block_top = min(2 * dyadic_n, y)
                for d, coefficient in a_coefficients.items():
                    lower = dyadic_n // d + 1
                    upper = block_top // d
                    aggregate_type_i += coefficient * sum(
                        mu[d * w + h] for w in range(lower, upper + 1)
                    )
                for d, coefficient in b_coefficients.items():
                    lower = max(u_cut, dyadic_n // d) + 1
                    upper = block_top // d
                    aggregate_type_ii += coefficient * sum(
                        mu[w] * mu[d * w + h]
                        for w in range(lower, upper + 1)
                    )
        dyadic_n *= 2
    if aggregate_target != boundary - aggregate_type_i + aggregate_type_ii:
        issues.append(
            "pre-collapse Vaughan aggregate identity failed: "
            f"target={aggregate_target}, "
            f"rhs={boundary - aggregate_type_i + aggregate_type_ii}"
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
    check_abel_and_energy(issues)
    check_first_block(issues)
    check_correlation_and_vaughan(issues)
    for issue in issues:
        print(f"MERTENS-WEIGHTED-PREFIX {issue}")
    audit = payload.get("audit", {})
    print(
        "validated Mertens weighted-prefix/affine-defect reduction: "
        f"{len(payload.get('rows', []))} rows, {len(issues)} issues, "
        f"{audit.get('exact_reduction_count', 0)} exact reductions, "
        f"{audit.get('proof_guard_count', 0)} proof guards, "
        f"{audit.get('open_handoff_gate_count', 0)} open handoff gates"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
