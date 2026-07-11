#!/usr/bin/env python3
"""Low-degree symbolic scout for the Jensen-window Cauchy-Binet ansatz.

This is theorem-search algebra, not a proof.  It checks whether the hard
degree-2/3/4 Jensen-window Toeplitz formulas are already certified by a much
weaker adjacent-log-concavity ratio parametrization.  If so, those low-degree
passes cannot be mistaken for a positive kernel/Cauchy-Binet theorem.
"""

from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ALGEBRA = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_obligation_algebra.json"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json"


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
    """Return a_j = A*rho^j*prod x_i^(j-i) for adjacent ratio contractions.

    The variables x1,x2,x3 encode:

      x1 = (a2/a1)/(a1/a0)
      x2 = (a3/a2)/(a2/a1)
      x3 = (a4/a3)/(a3/a2)

    For a positive log-concave sequence these variables lie in [0,1].
    """

    A, rho, x1, x2, x3 = sp.symbols("A rho x1 x2 x3")
    coeffs = coefficient_symbols()
    substitutions = {
        coeffs[0]: A,
        coeffs[1]: A * rho,
        coeffs[2]: A * rho**2 * x1,
        coeffs[3]: A * rho**3 * x1**2 * x2,
        coeffs[4]: A * rho**4 * x1**3 * x2**2 * x3,
    }
    return substitutions, (A, rho, x1, x2, x3), (x1, x2, x3)


def sympify_formula(formula: str) -> sp.Expr:
    coeffs = coefficient_symbols()
    locals_map = {f"a{i}": coeffs[i] for i in range(len(coeffs))}
    return sp.sympify(formula, locals=locals_map)


def positive_monomial_factor(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> tuple[sp.Expr, sp.Expr]:
    expanded = sp.expand(expr)
    poly = sp.Poly(expanded, *variables)
    min_powers = []
    for index, _var in enumerate(variables):
        min_powers.append(min(mon[index] for mon, coeff in poly.terms() if coeff != 0))
    monomial = sp.prod(var**power for var, power in zip(variables, min_powers))
    normalized = sp.factor(expanded / monomial)
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
    equations = [sp.Eq(coeff, 0) for coeff in diff.coeffs()]
    solution = sp.solve(equations, unknowns, dict=True)
    if not solution:
        raise RuntimeError(f"could not solve Bernstein coefficients for {expr_str(poly)}")
    solved = solution[0]
    rows = []
    for unknown, multi_index in zip(unknowns, indices):
        rows.append({"index": list(multi_index), "coefficient": expr_str(solved[unknown])})
    return degrees, rows


def coefficient_values(coefficients: list[dict]) -> list[sp.Rational]:
    return [sp.Rational(row["coefficient"]) for row in coefficients]


def formula_inputs(algebra: dict) -> list[dict]:
    rows: list[dict] = [
        {
            "id": "degree2_jensen_discriminant",
            "degree": 2,
            "rows": [],
            "cols": [],
            "determinant": algebra["degree2"]["jensen_discriminant"],
            "source_section": "degree2",
        }
    ]
    for degree_key in ("degree3", "degree4"):
        degree = int(degree_key[-1])
        for index, row in enumerate(algebra[degree_key]["selected_toeplitz_minors"], start=1):
            rows.append(
                {
                    "id": f"{degree_key}_selected_minor_{index:02d}",
                    "degree": degree,
                    "rows": row["rows"],
                    "cols": row["cols"],
                    "determinant": row["determinant"],
                    "source_section": degree_key,
                }
            )
    return rows


def analyze_formula(row: dict) -> dict:
    substitutions, positive_variables, contraction_variables = ratio_substitution()
    expr = sympify_formula(row["determinant"])
    ratio_expr = sp.factor(expr.subs(substitutions))
    monomial, normalized = positive_monomial_factor(ratio_expr, positive_variables)
    degrees, coefficients = bernstein_coefficients(normalized, contraction_variables)
    coeff_values = coefficient_values(coefficients)
    nonnegative = all(value >= 0 for value in coeff_values)
    positive = all(value > 0 for value in coeff_values)
    min_coeff = min(coeff_values)

    return {
        **row,
        "ratio_substitution": {
            "a0": "A",
            "a1": "A*rho",
            "a2": "A*rho**2*x1",
            "a3": "A*rho**3*x1**2*x2",
            "a4": "A*rho**4*x1**3*x2**2*x3",
            "log_concavity_box": "A>0, rho>0, 0<=x1,x2,x3<=1",
        },
        "ratio_factorization": expr_str(ratio_expr),
        "positive_monomial_factor": expr_str(monomial),
        "normalized_ratio_polynomial": expr_str(normalized),
        "bernstein_variables": [sp.sstr(var) for var in contraction_variables if normalized.has(var)],
        "bernstein_multidegree": degrees,
        "bernstein_coefficients": coefficients,
        "bernstein_min_coefficient": rational_str(min_coeff),
        "bernstein_coefficients_nonnegative": bool(nonnegative),
        "bernstein_coefficients_strictly_positive": bool(positive),
        "low_degree_log_concavity_certificate": bool(nonnegative),
        "proof_boundary": (
            "Low-degree adjacent-log-concavity certificate only; not a "
            "Cauchy-Binet kernel identity and not Jensen-window PF-infinity."
        ),
    }


def countermodel_ratio_check(algebra: dict) -> dict:
    values = [sp.Rational(item) for item in algebra["finite_countermodel"]["sequence_A0_to_A4"]]
    ratios = [sp.factor(values[index] / values[index - 1]) for index in range(1, len(values))]
    contractions = [
        sp.factor(ratios[index] / ratios[index - 1])
        for index in range(1, len(ratios))
    ]
    gaps = [
        sp.factor(values[index] ** 2 - values[index - 1] * values[index + 1])
        for index in range(1, len(values) - 1)
    ]
    counter = algebra["finite_countermodel"]
    return {
        "sequence_A0_to_A4": [rational_str(value) for value in values],
        "adjacent_ratios": [rational_str(value) for value in ratios],
        "ratio_contractions": [rational_str(value) for value in contractions],
        "ratio_contractions_in_unit_interval": all(0 < value <= 1 for value in contractions),
        "adjacent_log_concavity_gaps": [rational_str(value) for value in gaps],
        "adjacent_log_concavity_gaps_nonnegative": all(value >= 0 for value in gaps),
        "selected_low_degree_minors_positive": bool(
            counter["d3_selected_toeplitz_minors_positive"]
            and counter["d4_selected_toeplitz_minors_positive"]
        ),
        "d3_first_negative_contiguous_toeplitz_minor_size": counter["d3_first_negative_contiguous_toeplitz_minor"]["size"],
        "d4_first_negative_contiguous_toeplitz_minor_size": counter["d4_first_negative_contiguous_toeplitz_minor"]["size"],
        "interpretation": (
            "The finite countermodel sits inside the adjacent-log-concavity "
            "ratio box and passes the selected low-degree formulas, but still "
            "has negative larger contiguous Jensen-window Toeplitz minors."
        ),
    }


def build_payload(algebra: dict) -> dict:
    formula_rows = [analyze_formula(row) for row in formula_inputs(algebra)]
    nonnegative_count = sum(row["bernstein_coefficients_nonnegative"] for row in formula_rows)
    strict_count = sum(row["bernstein_coefficients_strictly_positive"] for row in formula_rows)
    countermodel = countermodel_ratio_check(algebra)
    return {
        "kind": "jensen_window_pf_cauchy_binet_low_degree_scout",
        "date": "2026-07-06",
        "target_ansatz": "ansatz_02_positive_cauchy_binet_kernel",
        "target_obligation": "jwpf_06_sign_regular_to_jensen_pf_conversion",
        "source_algebra": "work/rh_compute/results/jensen_window_pf_obligation_algebra.json",
        "proof_boundary": (
            "Symbolic low-degree theorem-search scout only; not a proof of a "
            "positive kernel, Cauchy-Binet identity, Jensen-window PF-infinity, "
            "Jensen hyperbolicity, RH, or Lambda <= 0."
        ),
        "ratio_parameterization": {
            "a_j": "A*rho**j*prod_{i=1}^{j-1} x_i**(j-i)",
            "meaning": "x_i are adjacent ratio contractions; positive log-concavity gives 0<=x_i<=1.",
            "certificate_method": "After extracting a positive monomial, nonnegative Bernstein coefficients on [0,1]^m certify the normalized low-degree formula.",
        },
        "summary": {
            "formula_rows": len(formula_rows),
            "bernstein_nonnegative_rows": nonnegative_count,
            "bernstein_strictly_positive_rows": strict_count,
            "cauchy_binet_identity_found": False,
            "kernel_identity_found": False,
            "target_closing": False,
            "main_finding": (
                "All selected degree-2/3/4 hard formulas are certified by "
                "adjacent log-concavity in the ratio box, so these low-degree "
                "passes are too weak to identify a positive Cauchy-Binet kernel."
            ),
        },
        "formula_rows": formula_rows,
        "countermodel_ratio_check": countermodel,
        "invariants": [
            "Every formula row has nonnegative Bernstein coefficients after positive monomial extraction.",
            "The degree-2 threshold may have a zero Bernstein boundary coefficient.",
            "No positive Cauchy-Binet or kernel identity is claimed.",
            "The log-concave finite countermodel still fails larger contiguous Jensen-window Toeplitz minors.",
        ],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--algebra", type=Path, default=DEFAULT_ALGEBRA)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    algebra = load_json(args.algebra)
    payload = build_payload(algebra)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        summary = payload["summary"]
        print(
            "wrote Jensen-window PF Cauchy-Binet low-degree scout: "
            f"{summary['formula_rows']} formula rows, "
            f"{summary['bernstein_nonnegative_rows']} Bernstein-certified rows, "
            "0 kernel identities found"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
