#!/usr/bin/env python3
"""Scout where adjacent log-concavity stops controlling Jensen-window PF minors.

This extends the low-degree Cauchy-Binet scout from selected formulas to the
first larger contiguous Jensen-window Toeplitz minors that the rational
countermodel makes negative.
"""

from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ALGEBRA = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_obligation_algebra.json"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json"


def expr_str(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(expr))


def rational_str(expr: sp.Expr) -> str:
    value = sp.Rational(expr)
    if value.q == 1:
        return str(value.p)
    return f"{value.p}/{value.q}"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def coefficient_symbols() -> tuple[sp.Symbol, ...]:
    return sp.symbols("a0:5")


def ratio_substitution() -> tuple[dict[sp.Symbol, sp.Expr], tuple[sp.Symbol, ...], tuple[sp.Symbol, ...]]:
    A, rho, x1, x2, x3 = sp.symbols("A rho x1 x2 x3")
    a = coefficient_symbols()
    return (
        {
            a[0]: A,
            a[1]: A * rho,
            a[2]: A * rho**2 * x1,
            a[3]: A * rho**3 * x1**2 * x2,
            a[4]: A * rho**4 * x1**3 * x2**2 * x3,
        },
        (A, rho, x1, x2, x3),
        (x1, x2, x3),
    )


def toeplitz_det(window: list[sp.Expr], size: int) -> sp.Expr:
    matrix = [
        [window[col - row] if 0 <= col - row < len(window) else sp.Integer(0) for col in range(1, size + 1)]
        for row in range(size)
    ]
    return sp.factor(sp.Matrix(matrix).det())


def positive_monomial_factor(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> tuple[sp.Expr, sp.Expr]:
    poly = sp.Poly(sp.expand(expr), *variables)
    min_powers = [
        min(mon[index] for mon, coeff in poly.terms() if coeff != 0)
        for index, _var in enumerate(variables)
    ]
    monomial = sp.prod(var**power for var, power in zip(variables, min_powers))
    normalized = sp.factor(sp.expand(expr) / monomial)
    return monomial, normalized


def bernstein_coefficients(poly: sp.Expr, variables: tuple[sp.Symbol, ...]) -> tuple[list[int], list[dict]]:
    active = tuple(var for var in variables if poly.has(var))
    if not active:
        return [], [{"index": [], "coefficient": expr_str(poly)}]

    power_poly = sp.Poly(sp.expand(poly), *active)
    degrees = [power_poly.degree(var) for var in active]
    indices = list(product(*[range(degree + 1) for degree in degrees]))
    unknowns = sp.symbols(f"b0:{len(indices)}")
    bernstein_sum = sp.Integer(0)
    for unknown, multi_index in zip(unknowns, indices):
        term = unknown
        for var, i, degree in zip(active, multi_index, degrees):
            term *= sp.binomial(degree, i) * var**i * (1 - var) ** (degree - i)
        bernstein_sum += term

    diff = sp.Poly(sp.expand(bernstein_sum - poly), *active)
    solution = sp.solve([sp.Eq(coeff, 0) for coeff in diff.coeffs()], unknowns, dict=True)
    if not solution:
        raise RuntimeError(f"could not solve Bernstein coefficients for {expr_str(poly)}")
    solved = solution[0]
    rows = [
        {"index": list(multi_index), "coefficient": expr_str(solved[unknown])}
        for unknown, multi_index in zip(unknowns, indices)
    ]
    return degrees, rows


def coefficient_values(coefficients: list[dict]) -> list[sp.Rational]:
    return [sp.Rational(row["coefficient"]) for row in coefficients]


def evaluate_contiguous(values: list[sp.Rational], degree: int, size: int) -> sp.Rational:
    window = [sp.binomial(degree, index) * values[index] for index in range(degree + 1)]
    matrix = [
        [window[col - row] if 0 <= col - row < len(window) else sp.Integer(0) for col in range(1, size + 1)]
        for row in range(size)
    ]
    return sp.factor(sp.Matrix(matrix).det())


def analyze_row(degree: int, size: int, values: list[sp.Rational]) -> dict:
    a = coefficient_symbols()
    substitutions, positive_variables, contraction_variables = ratio_substitution()
    window = [sp.binomial(degree, index) * a[index] for index in range(degree + 1)]
    determinant = toeplitz_det(window, size)
    ratio_expr = sp.factor(determinant.subs(substitutions))
    monomial, normalized = positive_monomial_factor(ratio_expr, positive_variables)
    degrees, coefficients = bernstein_coefficients(normalized, contraction_variables)
    coeff_values = coefficient_values(coefficients)
    counter_value = evaluate_contiguous(values, degree, size)
    return {
        "id": f"d{degree}_contiguous_m{size}",
        "degree": degree,
        "minor_size": size,
        "rows": list(range(size)),
        "cols": list(range(1, size + 1)),
        "determinant": expr_str(determinant),
        "ratio_factorization": expr_str(ratio_expr),
        "positive_monomial_factor": expr_str(monomial),
        "normalized_ratio_polynomial": expr_str(normalized),
        "normalized_ratio_terms": len(sp.Poly(sp.expand(normalized), *[var for var in contraction_variables if normalized.has(var)]).terms())
        if any(normalized.has(var) for var in contraction_variables)
        else 1,
        "bernstein_variables": [sp.sstr(var) for var in contraction_variables if normalized.has(var)],
        "bernstein_multidegree": degrees,
        "bernstein_coefficient_count": len(coefficients),
        "bernstein_min_coefficient": rational_str(min(coeff_values)),
        "bernstein_negative_coefficient_count": sum(1 for value in coeff_values if value < 0),
        "bernstein_coefficients_nonnegative": all(value >= 0 for value in coeff_values),
        "countermodel_value": rational_str(counter_value),
        "countermodel_sign": int(sp.sign(counter_value)),
        "countermodel_positive": bool(counter_value > 0),
        "proof_boundary": (
            "Contiguous-minor log-concavity frontier diagnostic only; not a "
            "positive kernel identity and not Jensen-window PF-infinity."
        ),
    }


def first_row(rows: list[dict], key: str, value: object) -> dict | None:
    for row in rows:
        if row.get(key) == value:
            return row
    return None


def build_payload(algebra: dict) -> dict:
    values = [sp.Rational(item) for item in algebra["finite_countermodel"]["sequence_A0_to_A4"]]
    rows_by_degree: dict[str, list[dict]] = {}
    for degree, max_size in ((3, 8), (4, 6)):
        rows_by_degree[str(degree)] = [
            analyze_row(degree, size, values) for size in range(1, max_size + 1)
        ]

    frontiers = {}
    for degree, rows in rows_by_degree.items():
        first_bernstein = first_row(rows, "bernstein_coefficients_nonnegative", False)
        first_counter = first_row(rows, "countermodel_positive", False)
        frontiers[degree] = {
            "first_negative_bernstein_size": None if first_bernstein is None else first_bernstein["minor_size"],
            "first_negative_bernstein_min_coefficient": None if first_bernstein is None else first_bernstein["bernstein_min_coefficient"],
            "first_countermodel_negative_size": None if first_counter is None else first_counter["minor_size"],
            "first_countermodel_negative_value": None if first_counter is None else first_counter["countermodel_value"],
        }

    algebra_counter = algebra["finite_countermodel"]
    return {
        "kind": "jensen_window_pf_log_concavity_frontier_scout",
        "date": "2026-07-06",
        "target_ansatz": "ansatz_02_positive_cauchy_binet_kernel",
        "source_algebra": "work/rh_compute/results/jensen_window_pf_obligation_algebra.json",
        "proof_boundary": (
            "Symbolic frontier diagnostic only; not a proof of a positive "
            "kernel, Cauchy-Binet identity, Jensen-window PF-infinity, RH, or "
            "Lambda <= 0."
        ),
        "ratio_parameterization": {
            "a0": "A",
            "a1": "A*rho",
            "a2": "A*rho**2*x1",
            "a3": "A*rho**3*x1**2*x2",
            "a4": "A*rho**4*x1**3*x2**2*x3",
            "log_concavity_box": "A>0, rho>0, 0<=x1,x2,x3<=1",
        },
        "summary": {
            "degree3_rows": len(rows_by_degree["3"]),
            "degree4_rows": len(rows_by_degree["4"]),
            "total_rows": sum(len(rows) for rows in rows_by_degree.values()),
            "kernel_identity_found": False,
            "target_closing": False,
            "main_finding": (
                "Adjacent log-concavity certifies only an initial contiguous "
                "minor frontier.  The Bernstein certificate first fails at "
                "degree 3 size 6 and degree 4 size 5; the exact log-concave "
                "countermodel first becomes negative at degree 3 size 8 and "
                "degree 4 size 6."
            ),
        },
        "frontiers": frontiers,
        "rows_by_degree": rows_by_degree,
        "algebra_countermodel_check": {
            "d3_expected_first_negative_size": algebra_counter["d3_first_negative_contiguous_toeplitz_minor"]["size"],
            "d3_expected_first_negative_value": algebra_counter["d3_first_negative_contiguous_toeplitz_minor"]["determinant"],
            "d4_expected_first_negative_size": algebra_counter["d4_first_negative_contiguous_toeplitz_minor"]["size"],
            "d4_expected_first_negative_value": algebra_counter["d4_first_negative_contiguous_toeplitz_minor"]["determinant"],
        },
        "interpretation": (
            "Selected low-degree positivity and adjacent log-concavity are "
            "strictly weaker than Jensen-window PF-infinity.  A valid "
            "Cauchy-Binet/kernel theorem must add structure that controls the "
            "larger contiguous minors and then all remaining minor shapes."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--algebra", type=Path, default=DEFAULT_ALGEBRA)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload(load_json(args.algebra))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "wrote Jensen-window PF log-concavity frontier scout: "
            f"{payload['summary']['total_rows']} contiguous rows, "
            "degree 3 first countermodel negative m=8, degree 4 first countermodel negative m=6"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
