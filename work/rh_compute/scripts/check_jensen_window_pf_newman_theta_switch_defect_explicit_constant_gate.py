#!/usr/bin/env python3
"""Independently check the explicit arbitrary-N modular-tail constants."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

import flint
from flint import arb
import sympy as sp


STEM = "jensen_window_pf_newman_theta_switch_defect_explicit_constant_gate"
DEFAULT_ARTIFACT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
RESULT_ROOT = REPO_ROOT / "work" / "rh_compute" / "results"
SOURCE_PATHS = {
    "modular_source_sha256": (
        RESULT_ROOT / "jensen_window_pf_newman_theta_modular_blend_gate.json"
    ),
    "envelope_source_sha256": (
        RESULT_ROOT
        / "jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate.json"
    ),
    "arbitrary_n_source_sha256": (
        RESULT_ROOT
        / "jensen_window_pf_newman_theta_arbitrary_n_stable_remainder_gate.json"
    ),
    "outer_source_sha256": (
        RESULT_ROOT
        / "jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate.json"
    ),
    "full_budget_source_sha256": (
        RESULT_ROOT
        / "jensen_window_pf_newman_theta_full_derivative_budget_certificate.json"
    ),
}
PRECISION_BITS = 256
MAX_ORDER = 9
M_ORDER = 9
FIRST_VALUES = (3, 4, 5, 8, 11, 17, 26, 65)
EXPECTED_IDS = [
    "ntsdecg_01_exact_kernel_recurrences",
    "ntsdecg_02_direct_tail_leibniz_envelope",
    "ntsdecg_03_forward_monomial_integral",
    "ntsdecg_04_reflected_compact_integral",
    "ntsdecg_05_reflected_large_integral",
    "ntsdecg_06_arbitrary_n_monomial_theorem",
    "ntsdecg_07_full_m9_d0_d1_compiler",
    "ntsdecg_08_three_quarter_scale_instantiation",
    "ntsdecg_09_cofinal_separation_handoff",
]


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def ball(value: Fraction | sp.Rational | int) -> arb:
    if isinstance(value, int):
        return arb(value)
    numerator = value.p if isinstance(value, sp.Rational) else value.numerator
    denominator = value.q if isinstance(value, sp.Rational) else value.denominator
    return arb(int(numerator)) / int(denominator)


def rational_power(base: arb, exponent: Fraction | sp.Rational) -> arb:
    return (ball(exponent) * base.log()).exp()


def independent_kernel_rows() -> list[dict]:
    variable = sp.symbols("X", nonnegative=True)
    rows: list[dict] = []
    expression = 2 * variable - 3
    for order in range(MAX_ORDER + 1):
        polynomial = sp.Poly(expression, variable)
        rows.append(
            {
                "q": order,
                "P_q": str(polynomial.as_expr()),
                "degree": polynomial.degree(),
                "A_q": int(
                    sum(abs(coefficient) for coefficient in polynomial.all_coeffs())
                ),
            }
        )
        expression = sp.expand(
            5 * expression
            + 4 * variable * sp.diff(expression, variable)
            - 4 * variable * expression
        )
    return rows


def laurent_norm(expression: sp.Expr, variable: sp.Symbol) -> tuple[sp.Rational, int]:
    norm = sp.Rational(0)
    degrees: list[int] = []
    for term in sp.Add.make_args(sp.expand(expression)):
        coefficient, exponent = term.as_coeff_exponent(variable)
        if coefficient.has(variable) or not exponent.is_integer:
            raise RuntimeError("non-Laurent switch recurrence")
        norm += abs(coefficient)
        degrees.append(int(exponent))
    return norm, max(degrees)


def independent_switch_rows() -> list[dict]:
    variable = sp.symbols("y", positive=True)
    switch_argument = sp.Rational(3, 2) * (
        variable - variable**-1
    )
    argument_derivative = 4 * variable * sp.diff(
        switch_argument,
        variable,
    )
    expression = -argument_derivative
    rows = [
        {
            "a": 0,
            "laurent_polynomial": None,
            "coefficient_norm": "1/2",
            "maximum_degree": 0,
            "normalization": "1",
        }
    ]
    for order in range(1, MAX_ORDER + 1):
        norm, degree = laurent_norm(expression, variable)
        rows.append(
            {
                "a": order,
                "laurent_polynomial": str(expression),
                "coefficient_norm": str(norm),
                "maximum_degree": degree,
                "normalization": "1/sqrt(pi)",
            }
        )
        expression = sp.expand(
            4 * variable * sp.diff(expression, variable)
            - 2
            * switch_argument
            * argument_derivative
            * expression
        )
    return rows


def forward_integral(
    moment: int,
    q_power: int,
    y_power: sp.Rational,
    first: int,
) -> arb:
    pi = arb.pi()
    summation_degree = 2 * q_power - 2
    logarithmic_power = max(
        y_power - sp.Rational(3, 4),
        sp.Rational(0),
    )
    shift = logarithmic_power + sp.Rational(1, 80)
    denominator_factor = 1 - ball(shift) / (pi * first**2)
    geometric_ratio = (
        summation_degree / arb(first) - pi * (2 * first + 1)
    ).exp()
    if denominator_factor.lower() <= 0 or geometric_ratio.upper() >= 1:
        raise RuntimeError("invalid independent forward-tail condition")
    leading = (
        factorial(moment)
        * (arb(1) / 80).exp()
        * pi ** (q_power - 1)
        * first**summation_degree
        * (-pi * first**2).exp()
    )
    return leading / (
        4 * denominator_factor * (1 - geometric_ratio)
    )


def compact_reflected_integral(
    moment: int,
    q_power: int,
    y_power: sp.Rational,
    first: int,
) -> arb:
    pi = arb.pi()
    root_two = arb(2).sqrt()
    interval_length = arb(2).log() / 8
    degree = 2 * q_power
    ratio = (
        degree / arb(first)
        - pi * (2 * first + 1) / root_two
    ).exp()
    if ratio.upper() >= 1:
        raise RuntimeError("invalid independent compact-tail ratio")
    power_factor = (ball(y_power / 2) * arb(2).log()).exp()
    return (
        interval_length ** (moment + 1)
        * (interval_length**2 / 5).exp()
        * power_factor
        * pi**q_power
        * first**degree
        * (-pi * first**2 / root_two).exp()
        / (1 - ratio)
    )


def independent_safe_constant() -> arb:
    return (
        arb(3)
        / 16
        * rational_power(arb(3), Fraction(2, 3))
        * rational_power(arb.pi(), Fraction(2, 3))
    )


def independent_discrete_tail(q_power: int, first: int) -> arb:
    degree = 2 * q_power
    decay = independent_safe_constant()
    left = arb(first - 1)
    critical = rational_power(
        3 * degree / (4 * decay),
        Fraction(3, 4),
    )
    if critical.lower() >= left:
        maximum_at = critical
    elif critical.upper() <= left:
        maximum_at = left
    else:
        raise RuntimeError("independent mode split is indeterminate")
    maximum = (
        maximum_at**degree
        * (
            -decay
            * rational_power(maximum_at, Fraction(4, 3))
        ).exp()
    )
    gamma_shape = Fraction(3 * (degree + 1), 4)
    gamma_point = (
        decay * rational_power(left, Fraction(4, 3))
    )
    integral = (
        arb(3)
        * rational_power(decay, -gamma_shape)
        * gamma_point.gamma_upper(ball(gamma_shape))
        / 4
    )
    return maximum + integral


def large_reflected_integral(
    moment: int,
    q_power: int,
    y_power: sp.Rational,
    first: int,
) -> arb:
    residual_decay = sp.Rational(27, 128)
    young_epsilon = sp.Rational(27, 64)
    linear_heat = sp.Rational(1, 80)
    gamma_shape = (
        y_power - sp.Rational(3, 4) + 1
    ) / 2
    if gamma_shape <= 0:
        raise RuntimeError("independent gamma shape is not positive")
    young_factor = ball(
        linear_heat**2 / (2 * young_epsilon)
    ).exp()
    gamma_point = ball(2 * residual_decay)
    continuous_factor = (
        factorial(moment)
        * young_factor
        * rational_power(ball(residual_decay), -gamma_shape)
        * gamma_point.gamma_upper(ball(gamma_shape))
        / 8
    )
    return (
        continuous_factor
        * arb.pi() ** q_power
        * independent_discrete_tail(q_power, first)
    )


def template_integral(
    moment: int,
    derivative: int,
    first: int,
    kernels: list[dict],
    switches: list[dict],
) -> dict[str, arb]:
    answer = {
        "forward": arb(0),
        "reflected_compact": arb(0),
        "reflected_large": arb(0),
    }
    pi = arb.pi()
    for switch_order in range(derivative + 1):
        kernel_order = derivative - switch_order
        q_power = kernel_order + 2
        binomial = comb(derivative, switch_order)
        kernel_norm = kernels[kernel_order]["A_q"]

        if switch_order:
            switch_norm = (
                ball(sp.Rational(switches[switch_order]["coefficient_norm"]))
                / pi.sqrt()
            )
            switch_degree = switches[switch_order]["maximum_degree"]
            coefficient_forward = (
                binomial * kernel_norm * switch_norm
            )
            coefficient_reflected = (
                binomial
                * kernel_norm
                * 2 ** (kernel_order + 1)
                * switch_norm
            )
            forward_power = (
                switch_degree + kernel_order + sp.Rational(9, 4)
            )
            reflected_power = max(
                sp.Rational(switch_degree) - sp.Rational(5, 4),
                sp.Rational(0),
            )
        else:
            coefficient_forward = arb(binomial * kernel_norm)
            coefficient_reflected = ball(
                sp.Rational(
                    binomial
                    * kernel_norm
                    * 2 ** (kernel_order + 1),
                    2,
                )
            )
            forward_power = kernel_order + sp.Rational(9, 4)
            reflected_power = sp.Rational(0)

        answer["forward"] += (
            coefficient_forward
            * forward_integral(
                moment,
                q_power,
                forward_power,
                first,
            )
        )
        compact = compact_reflected_integral(
            moment,
            q_power,
            reflected_power,
            first,
        )
        large = large_reflected_integral(
            moment,
            q_power,
            reflected_power,
            first,
        )
        answer["reflected_compact"] += coefficient_reflected * compact
        answer["reflected_large"] += coefficient_reflected * large
    answer["total"] = (
        answer["forward"]
        + answer["reflected_compact"]
        + answer["reflected_large"]
    )
    return answer


def independent_heat_terms(order: int) -> list[tuple[int, Fraction]]:
    answer: list[tuple[int, Fraction]] = []
    for ell in range(order // 2 + 1):
        monomial = order - 2 * ell
        answer.append(
            (
                monomial,
                Fraction(
                    factorial(order) * 2**monomial,
                    factorial(ell)
                    * factorial(monomial)
                    * 5 ** (order - ell),
                ),
            )
        )
    return answer


def zero_components() -> dict[str, arb]:
    return {
        "forward": arb(0),
        "reflected_compact": arb(0),
        "reflected_large": arb(0),
        "total": arb(0),
    }


def accumulate(
    target: dict[str, arb],
    source: dict[str, arb],
    coefficient: Fraction | int,
) -> None:
    coefficient_ball = ball(coefficient)
    for component in target:
        target[component] += coefficient_ball * source[component]


def independently_compile(
    first: int,
    kernels: list[dict],
    switches: list[dict],
) -> dict[str, dict[str, arb] | int]:
    memo: dict[tuple[int, int], dict[str, arb]] = {}

    def get(moment: int, derivative: int) -> dict[str, arb]:
        index = (moment, derivative)
        if index not in memo:
            memo[index] = template_integral(
                moment,
                derivative,
                first,
                kernels,
                switches,
            )
        return memo[index]

    d0 = zero_components()
    d1_weighted = zero_components()
    d1_product = zero_components()
    for heat_order in range(10):
        derivative = 9 - heat_order
        outer_coefficient = comb(9, heat_order)
        for moment, heat_coefficient in independent_heat_terms(heat_order):
            full_coefficient = outer_coefficient * heat_coefficient
            accumulate(d0, get(moment, derivative), full_coefficient)
            accumulate(
                d1_weighted,
                get(moment + 1, derivative),
                full_coefficient,
            )
    for heat_order in range(9):
        derivative = 8 - heat_order
        outer_coefficient = 9 * comb(8, heat_order)
        for moment, heat_coefficient in independent_heat_terms(heat_order):
            accumulate(
                d1_product,
                get(moment, derivative),
                outer_coefficient * heat_coefficient,
            )
    d1 = {
        component: d1_weighted[component] + d1_product[component]
        for component in d1_weighted
    }
    return {
        "memo_count": len(memo),
        "d0": d0,
        "d1_weighted": d1_weighted,
        "d1_product": d1_product,
        "d1": d1,
    }


def compare_component_set(
    issues: list[str],
    stored: dict,
    rebuilt: dict[str, arb],
    label: str,
) -> None:
    for component, value in rebuilt.items():
        try:
            stored_value = arb(stored[component]["enclosure"])
            stored_log = arb(stored[component]["log10_enclosure"])
        except Exception as exc:
            issues.append(f"{label} {component} parse failed: {exc}")
            continue
        if not stored_value.overlaps(value):
            issues.append(f"{label} {component} enclosure mismatch")
        rebuilt_log = value.log() / arb(10).log()
        if not stored_log.overlaps(rebuilt_log):
            issues.append(f"{label} {component} log10 mismatch")
        if stored_value.rel_accuracy_bits() < 200:
            issues.append(f"{label} {component} has weak stored accuracy")


def verify_symbolic_inequalities(issues: list[str]) -> None:
    u, y, q = sp.symbols("u y q", positive=True)
    if sp.simplify(
        9 * sp.sinh(4 * u) ** 2
        - sp.Rational(9, 4) * (sp.exp(4 * u) - sp.exp(-4 * u)) ** 2
    ) != 0:
        issues.append("switch phase identity failed")
    large_difference = sp.factor(
        sp.Rational(9, 4) * (y - y**-1) ** 2
        - sp.Rational(9, 16) * y**2
    )
    large_factorization = (
        sp.Rational(9, 16)
        * (3 * y**2 - 2)
        * (y**2 - 2)
        / y**2
    )
    if sp.simplify(large_difference - large_factorization) != 0:
        issues.append("large-y switch lower bound factorization failed")
    minimizer = (8 * q / 9) ** sp.Rational(1, 3)
    phase = q / minimizer + sp.Rational(9, 16) * minimizer**2
    expected = (
        sp.Rational(3, 2)
        * (sp.Rational(9, 8)) ** sp.Rational(1, 3)
        * q ** sp.Rational(2, 3)
    )
    if sp.simplify(phase - expected) != 0:
        issues.append("reflected phase minimum failed")
    if arb(2).exp().lower() <= 4:
        issues.append("u^2<=exp(4u)/16 numerical endpoint failed")


def validate(path: Path) -> list[str]:
    flint.ctx.prec = PRECISION_BITS
    issues: list[str] = []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("kind") != STEM:
        issues.append("artifact kind mismatch")
    expected_parameters = {
        "precision_bits": PRECISION_BITS,
        "max_derivative_order": MAX_ORDER,
        "m_order": M_ORDER,
        "t_cap_exact": "1/5",
        "minimum_retained_count": 2,
        "minimum_first_omitted": 3,
        "retained_count_upper_bound": None,
        "witness_K": list(FIRST_VALUES),
    }
    if artifact.get("parameters") != expected_parameters:
        issues.append("parameter mismatch")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if [row.get("readiness") for row in rows] != (
        ["proved"] * 8 + ["not_ready_to_apply"]
    ):
        issues.append("readiness mismatch")
    if len(rows) == len(EXPECTED_IDS) and rows[-1].get("formula") != (
        "|J-J_N|<=16*d0/x^5; "
        "|J'-J_N'|<=64*d0/x^6+16*d1/x^5"
    ):
        issues.append("m=9 Fourier error multiplier mismatch")

    audit = artifact.get("source_audit", {})
    for key, source in SOURCE_PATHS.items():
        if audit.get(key) != digest(source):
            issues.append(f"source hash mismatch: {key}")
    if audit.get("tail_identity") != (
        "r_N(u)=sum_(n>N)b_n(u)=Phi(u)-sum_(n<=N)b_n(u)"
    ):
        issues.append("tail identity mismatch")
    if audit.get("safe_exponent") != "c_*=c_0/4":
        issues.append("safe exponent mismatch")

    kernels = independent_kernel_rows()
    if artifact.get("kernel_polynomials") != kernels:
        issues.append("independent kernel recurrence mismatch")
    switches = independent_switch_rows()
    stored_switches = artifact.get("switch_derivatives", [])
    if len(stored_switches) != len(switches):
        issues.append("switch row count mismatch")
    else:
        for expected, stored in zip(switches, stored_switches, strict=True):
            for key, value in expected.items():
                if key == "laurent_polynomial" and value is not None:
                    observed_expression = sp.sympify(stored.get(key))
                    expected_expression = sp.sympify(value)
                    matches = (
                        sp.simplify(observed_expression - expected_expression)
                        == 0
                    )
                else:
                    matches = stored.get(key) == value
                if not matches:
                    issues.append(
                        f"independent switch recurrence mismatch a={expected['a']} key={key}"
                    )

    witness_rows = artifact.get("budgets", [])
    if [row.get("K") for row in witness_rows] != list(FIRST_VALUES):
        issues.append("witness K sequence mismatch")
    previous_d0: arb | None = None
    previous_d1: arb | None = None
    for row in witness_rows:
        first = int(row["K"])
        rebuilt = independently_compile(first, kernels, switches)
        if row.get("N") != first - 1:
            issues.append(f"retained/omitted index mismatch at K={first}")
        if row.get("monomial_integrals_reused") != rebuilt["memo_count"]:
            issues.append(f"memoized integral count mismatch at K={first}")
        compare_component_set(
            issues,
            row["d0_m9_t1_5"],
            rebuilt["d0"],
            f"K={first} d0",
        )
        compare_component_set(
            issues,
            row["d1_u_m9_t1_5"],
            rebuilt["d1_weighted"],
            f"K={first} d1-u",
        )
        compare_component_set(
            issues,
            row["d1_lower_m9_t1_5"],
            rebuilt["d1_product"],
            f"K={first} d1-lower",
        )
        compare_component_set(
            issues,
            row["d1_m9_t1_5"],
            rebuilt["d1"],
            f"K={first} d1",
        )
        current_d0 = rebuilt["d0"]["total"]
        current_d1 = rebuilt["d1"]["total"]
        if previous_d0 is not None and current_d0.upper() >= previous_d0.lower():
            issues.append(f"d0 monotonicity failed at K={first}")
        if previous_d1 is not None and current_d1.upper() >= previous_d1.lower():
            issues.append(f"d1 monotonicity failed at K={first}")
        previous_d0 = current_d0
        previous_d1 = current_d1

    verify_symbolic_inequalities(issues)
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "every N>=2",
        "retained first-jet",
        "unbounded adaptive transition",
        "strict Laguerre",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary guard missing: {marker}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        return 1
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman theta switch-defect explicit constant gate: "
        f"{len(artifact['rows'])} rows, "
        f"{len(artifact['budgets'])} all-N formula witnesses, 0 issues, "
        "1 open cofinal retained-separation handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
