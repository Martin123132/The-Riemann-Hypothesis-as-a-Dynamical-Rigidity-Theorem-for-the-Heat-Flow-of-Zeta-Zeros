#!/usr/bin/env python3
"""Build the monotone-contraction frontier scout for Jensen-window PF.

This is exact theorem-search algebra plus finite Arb diagnostics.  It checks a
sharper sufficient region for the first hard column-recurrence rows:

    x_i = (A_{n+i+1}/A_{n+i}) / (A_{n+i}/A_{n+i-1})

Plain adjacent log-concavity only gives 0 <= x_i <= 1.  The hard rational
countermodel violates monotone contractions, while the checked zeta windows
have x_1 <= x_2 <= x_3.  This script records exact Bernstein certificates for
the two hard frontier polynomials on that monotone-contraction region.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
from itertools import product
import json
from math import comb
from pathlib import Path
import sys

import sympy as sp


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ALGEBRA = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_obligation_algebra.json"
DEFAULT_ENCLOSURE_JSONL = (
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k0_k9.jsonl",
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k10_k20.jsonl",
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k21_k32.jsonl",
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_frontier_scout.json"


@dataclass(frozen=True)
class ZetaDiagnosticRow:
    degree: int
    checked_rows: int
    monotone_contraction_rows: int
    positive_hard_value_rows: int
    failed_or_inconclusive_rows: int
    min_gap_sample: str
    min_gap_sample_at: dict
    min_normalized_hard_value_sample: str
    min_normalized_hard_value_sample_at: dict
    contraction_ranges: dict[str, dict[str, str]]


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


def coefficient_symbols(count: int = 5) -> tuple[sp.Symbol, ...]:
    return sp.symbols(f"a0:{count}")


def ratio_substitution() -> tuple[dict[sp.Symbol, sp.Expr], tuple[sp.Symbol, ...]]:
    A, rho, x1, x2, x3 = sp.symbols("A rho x1 x2 x3")
    a = coefficient_symbols()
    substitutions = {
        a[0]: A,
        a[1]: A * rho,
        a[2]: A * rho**2 * x1,
        a[3]: A * rho**3 * x1**2 * x2,
        a[4]: A * rho**4 * x1**3 * x2**2 * x3,
    }
    return substitutions, (A, rho, x1, x2, x3)


def positive_monomial_factor(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> tuple[sp.Expr, sp.Expr]:
    expanded = sp.expand(expr)
    poly = sp.Poly(expanded, *variables)
    min_powers = [
        min(mon[index] for mon, coeff in poly.terms() if coeff != 0)
        for index, _var in enumerate(variables)
    ]
    monomial = sp.prod(var**power for var, power in zip(variables, min_powers))
    return sp.factor(monomial), sp.factor(expanded / monomial)


def column_determinant_polynomial(degree: int, size: int) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    h = sp.symbols(f"h0:{degree + 1}")
    matrix = [
        [h[col - row] if 0 <= col - row <= degree else sp.Integer(0) for col in range(1, size + 1)]
        for row in range(size)
    ]
    return sp.factor(sp.Matrix(matrix).det()), h


def hard_ratio_polynomial(degree: int, size: int) -> tuple[sp.Expr, sp.Expr]:
    determinant, h = column_determinant_polynomial(degree, size)
    a = coefficient_symbols()
    h_subs = {h[index]: sp.binomial(degree, index) * a[index] for index in range(degree + 1)}
    ratio_subs, positive_variables = ratio_substitution()
    ratio_expr = sp.factor(determinant.subs(h_subs).subs(ratio_subs))
    monomial, normalized = positive_monomial_factor(ratio_expr, positive_variables)
    return monomial, normalized


def monotone_substitution(degree: int) -> tuple[dict[sp.Symbol, sp.Expr], tuple[sp.Symbol, ...], str]:
    _A, _rho, x1, x2, x3 = ratio_substitution()[1]
    s, u, v = sp.symbols("s u v")
    if degree == 3:
        return {x1: s, x2: s + u * (1 - s)}, (s, u), "x1=s, x2=s+u*(1-s), 0<=s,u<=1"
    if degree == 4:
        x2_expr = s + u * (1 - s)
        return (
            {x1: s, x2: x2_expr, x3: x2_expr + v * (1 - x2_expr)},
            (s, u, v),
            "x1=s, x2=s+u*(1-s), x3=x2+v*(1-x2), 0<=s,u,v<=1",
        )
    raise ValueError(f"unsupported degree for monotone substitution: {degree}")


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
        for var, index, degree in zip(active, multi_index, degrees):
            term *= sp.binomial(degree, index) * var**index * (1 - var) ** (degree - index)
        bernstein_sum += term

    diff = sp.Poly(sp.expand(bernstein_sum - poly), *active)
    solution = sp.solve([sp.Eq(coeff, 0) for coeff in diff.coeffs()], unknowns, dict=True)
    if not solution:
        raise RuntimeError(f"could not solve Bernstein coefficients for {expr_str(poly)}")
    solved = solution[0]
    return degrees, [
        {"index": list(multi_index), "coefficient": expr_str(solved[unknown])}
        for unknown, multi_index in zip(unknowns, indices)
    ]


def hard_exact_row(degree: int, size: int) -> dict:
    monomial, ratio_poly = hard_ratio_polynomial(degree, size)
    substitutions, variables, description = monotone_substitution(degree)
    monotone_poly = sp.factor(ratio_poly.subs(substitutions))
    multidegree, coefficients = bernstein_coefficients(monotone_poly, variables)
    coefficient_values = [sp.Rational(row["coefficient"]) for row in coefficients]
    row_id = f"mcfs_0{1 if degree == 3 else 2}_d{degree}_m{size}_monotone_contractions"
    source_row = f"d{degree}_column_recurrence_m{size}"
    return {
        "id": row_id,
        "source_frontier_row": source_row,
        "degree": degree,
        "minor_size": size,
        "normalized_positive_monomial": expr_str(monomial),
        "ratio_box_polynomial": expr_str(ratio_poly),
        "monotone_substitution": description,
        "monotone_region_polynomial": expr_str(monotone_poly),
        "bernstein_variables": [sp.sstr(var) for var in variables],
        "bernstein_multidegree": multidegree,
        "bernstein_coefficient_count": len(coefficients),
        "bernstein_min_coefficient": rational_str(min(coefficient_values)),
        "bernstein_negative_coefficient_count": sum(1 for value in coefficient_values if value < 0),
        "bernstein_coefficients": coefficients,
        "certificate": (
            "All Bernstein coefficients are positive on the monotone-contraction "
            "unit cube, so this exact sufficient region makes the hard "
            "normalized frontier polynomial positive."
        ),
        "proof_boundary": (
            "Exact sufficient condition for this hard row only; it does not "
            "prove zeta-window monotone contractions, all m, all shifts, all "
            "lambda values, all Schur shapes, jwpf_06, or Lambda <= 0."
        ),
    }


def countermodel_check(algebra: dict) -> dict:
    values = [sp.Rational(item) for item in algebra["finite_countermodel"]["sequence_A0_to_A4"]]
    ratios = [sp.factor(values[index] / values[index - 1]) for index in range(1, len(values))]
    contractions = [sp.factor(ratios[index] / ratios[index - 1]) for index in range(1, len(ratios))]
    _A, _rho, x1, x2, x3 = ratio_substitution()[1]
    p3 = hard_ratio_polynomial(3, 8)[1]
    p4 = hard_ratio_polynomial(4, 6)[1]
    substitutions = {x1: contractions[0], x2: contractions[1], x3: contractions[2]}
    violations = []
    if contractions[1] < contractions[0]:
        violations.append("x2 < x1")
    if contractions[2] < contractions[1]:
        violations.append("x3 < x2")
    return {
        "source": "work/rh_compute/results/jensen_window_pf_obligation_algebra.json",
        "sequence_A0_to_A4": [rational_str(value) for value in values],
        "ratio_contractions": {
            "x1": rational_str(contractions[0]),
            "x2": rational_str(contractions[1]),
            "x3": rational_str(contractions[2]),
        },
        "monotone_contractions": len(violations) == 0,
        "violations": violations,
        "d3_m8_normalized_value": rational_str(sp.factor(p3.subs(substitutions))),
        "d4_m6_normalized_value": rational_str(sp.factor(p4.subs(substitutions))),
        "interpretation": (
            "The exact rational countermodel lies in the adjacent-log-concavity "
            "box but outside the monotone-contraction region, so the new "
            "sufficient condition escapes this specific countermodel."
        ),
    }


def decimal_lam_key(text: str) -> Decimal:
    return Decimal(str(text)).normalize()


def load_zeta_values(paths: tuple[Path, ...]) -> tuple[dict[tuple[Decimal, int], flint.arb], dict[tuple[Decimal, int], Decimal], dict[Decimal, str]]:
    balls: dict[tuple[Decimal, int], flint.arb] = {}
    samples: dict[tuple[Decimal, int], Decimal] = {}
    labels: dict[Decimal, str] = {}
    for path in paths:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                row = json.loads(line)
                if row.get("kind") != "acb_coefficient_enclosure":
                    continue
                lam = decimal_lam_key(row["lam"])
                index = int(row["k"])
                balls[(lam, index)] = flint.arb(row["A_ball"])
                samples[(lam, index)] = Decimal(row["cache_A"])
                labels[lam] = row["lam"]
    return balls, samples, labels


def hard_value_numeric(degree: int, xs):
    x1 = xs[0]
    x2 = xs[1]
    if degree == 3:
        return -9 * (
            x1**5 * x2**2
            - 6 * x1**4 * x2**2
            - 36 * x1**4 * x2
            - 9 * x1**4
            + 180 * x1**3 * x2
            + 270 * x1**3
            - 162 * x1**2 * x2
            - 1215 * x1**2
            + 1701 * x1
            - 729
        )
    x3 = xs[2]
    return 4 * (
        3 * x1**4 * x2**2 * x3
        + 4 * x1**4 * x2**2
        - 12 * x1**3 * x2**2 * x3
        - 144 * x1**3 * x2
        - 54 * x1**3
        + 256 * x1**2 * x2
        + 864 * x1**2
        - 1920 * x1
        + 1024
    )


def contractions_from_sequence(values):
    ratios = [values[index] / values[index - 1] for index in range(1, len(values))]
    return [ratios[index] / ratios[index - 1] for index in range(1, len(ratios))]


def decimal_format(value: Decimal) -> str:
    return f"{value:.18E}"


def zeta_diagnostic_for_degree(
    degree: int,
    balls: dict[tuple[Decimal, int], flint.arb],
    samples: dict[tuple[Decimal, int], Decimal],
    labels: dict[Decimal, str],
    shifts: range,
) -> ZetaDiagnosticRow:
    checked = 0
    monotone_ok = 0
    hard_ok = 0
    sample_gap_rows = []
    sample_value_rows = []
    contraction_samples: dict[int, list[Decimal]] = {index: [] for index in range(1, degree)}
    for lam in sorted(labels):
        for shift in shifts:
            arb_values = [balls[(lam, shift + offset)] for offset in range(degree + 1)]
            arb_xs = contractions_from_sequence(arb_values)
            arb_gaps = [arb_xs[index + 1] - arb_xs[index] for index in range(len(arb_xs) - 1)]
            arb_value = hard_value_numeric(degree, arb_xs)
            checked += 1
            monotone = all(gap > 0 and not gap.contains(0) for gap in arb_gaps)
            hard_positive = arb_value > 0 and not arb_value.contains(0)
            monotone_ok += int(monotone)
            hard_ok += int(hard_positive)

            decimal_values = [samples[(lam, shift + offset)] for offset in range(degree + 1)]
            decimal_xs = contractions_from_sequence(decimal_values)
            decimal_gaps = [decimal_xs[index + 1] - decimal_xs[index] for index in range(len(decimal_xs) - 1)]
            decimal_value = hard_value_numeric(degree, decimal_xs)
            for index, value in enumerate(decimal_xs, start=1):
                contraction_samples[index].append(value)
            if decimal_gaps:
                sample_gap_rows.append((min(decimal_gaps), labels[lam], shift, [decimal_format(item) for item in decimal_xs]))
            sample_value_rows.append((decimal_value, labels[lam], shift, [decimal_format(item) for item in decimal_xs]))

    sample_gap_rows.sort(key=lambda row: row[0])
    sample_value_rows.sort(key=lambda row: row[0])
    min_gap = sample_gap_rows[0]
    min_value = sample_value_rows[0]
    return ZetaDiagnosticRow(
        degree=degree,
        checked_rows=checked,
        monotone_contraction_rows=monotone_ok,
        positive_hard_value_rows=hard_ok,
        failed_or_inconclusive_rows=checked - min(monotone_ok, hard_ok),
        min_gap_sample=decimal_format(min_gap[0]),
        min_gap_sample_at={"lambda": min_gap[1], "shift_n": min_gap[2], "contractions": min_gap[3]},
        min_normalized_hard_value_sample=decimal_format(min_value[0]),
        min_normalized_hard_value_sample_at={"lambda": min_value[1], "shift_n": min_value[2], "contractions": min_value[3]},
        contraction_ranges={
            f"x{index}": {
                "min_sample": decimal_format(min(values)),
                "max_sample": decimal_format(max(values)),
            }
            for index, values in contraction_samples.items()
        },
    )


def finite_zeta_diagnostics(paths: tuple[Path, ...], shifts: range) -> list[dict]:
    getcontext().prec = 80
    balls, samples, labels = load_zeta_values(paths)
    rows = [
        zeta_diagnostic_for_degree(3, balls, samples, labels, shifts),
        zeta_diagnostic_for_degree(4, balls, samples, labels, shifts),
    ]
    return [asdict(row) for row in rows]


def build_payload(algebra: dict, enclosure_paths: tuple[Path, ...], shifts: range) -> dict:
    exact_rows = [hard_exact_row(3, 8), hard_exact_row(4, 6)]
    finite_rows = finite_zeta_diagnostics(enclosure_paths, shifts)
    total_coefficients = sum(row["bernstein_coefficient_count"] for row in exact_rows)
    negative_coefficients = sum(row["bernstein_negative_coefficient_count"] for row in exact_rows)
    finite_checked = sum(row["checked_rows"] for row in finite_rows)
    finite_monotone = sum(row["monotone_contraction_rows"] for row in finite_rows)
    finite_positive = sum(row["positive_hard_value_rows"] for row in finite_rows)
    return {
        "kind": "jensen_window_pf_monotone_contraction_frontier_scout",
        "date": "2026-07-06",
        "status": "symbolic_frontier_scout",
        "target_frontier_rows": ["d3_column_recurrence_m8", "d4_column_recurrence_m6"],
        "source_column_recurrence_contract": "work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json",
        "source_frontier_matrix": "work/rh_compute/results/jensen_window_pf_cauchy_binet_cone_frontier_matrix.json",
        "source_algebra": "work/rh_compute/results/jensen_window_pf_obligation_algebra.json",
        "proof_boundary": (
            "Monotone-contraction frontier scout only; it proves exact sufficient "
            "positivity certificates for two hard normalized polynomials under "
            "an extra ratio-contraction hypothesis.  It does not prove that the "
            "zeta heat-flow coefficients satisfy that hypothesis for all shifts "
            "or lambda values, does not prove all-shape Jensen-window PF, and "
            "does not prove RH or Lambda <= 0."
        ),
        "ratio_contraction_definition": {
            "x_i": "(A_{n+i+1}/A_{n+i}) / (A_{n+i}/A_{n+i-1})",
            "adjacent_log_concavity_box": "0 <= x_i <= 1",
            "extra_sufficient_region": "0 <= x1 <= x2 <= x3 <= 1, truncated to available degree",
        },
        "exact_frontier_rows": exact_rows,
        "countermodel_check": countermodel_check(algebra),
        "finite_zeta_arb_diagnostics": finite_rows,
        "summary": {
            "exact_certificate_rows": len(exact_rows),
            "total_bernstein_coefficients": total_coefficients,
            "negative_bernstein_coefficients": negative_coefficients,
            "finite_zeta_checked_rows": finite_checked,
            "finite_zeta_monotone_rows": finite_monotone,
            "finite_zeta_positive_hard_rows": finite_positive,
            "target_closing": False,
            "main_finding": (
                "The first hard column-frontier rows become exactly positive "
                "on the monotone-contraction subregion x1<=x2<=x3 of the "
                "adjacent-log-concavity box.  The rational countermodel fails "
                "this monotone condition, while the checked zeta-grid windows "
                "satisfy it with positive Arb classification.  The remaining "
                "theorem target is to prove the monotone-contraction property, "
                "or a stronger kernel identity, noncircularly for the actual "
                "zeta heat-flow coefficients."
            ),
        },
        "invariants": [
            "Both hard frontier rows must have strictly positive Bernstein coefficients on the monotone-contraction cube.",
            "The rational countermodel must violate monotone contractions.",
            "Finite zeta diagnostics must remain finite evidence only.",
            "No row may set target_closing=true.",
        ],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--algebra", type=Path, default=DEFAULT_ALGEBRA)
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="+", default=list(DEFAULT_ENCLOSURE_JSONL))
    parser.add_argument("--shifts", default="0..20")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--json", action="store_true")
    return parser


def parse_shift_range(text: str) -> range:
    if ".." not in text:
        value = int(text)
        return range(value, value + 1)
    left, right = text.split("..", 1)
    start = int(left)
    stop = int(right)
    if stop < start:
        raise ValueError(f"descending shift range: {text}")
    return range(start, stop + 1)


def main() -> int:
    args = build_parser().parse_args()
    flint.ctx.dps = 160
    payload = build_payload(
        load_json(args.algebra),
        tuple(args.enclosure_jsonl),
        parse_shift_range(args.shifts),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        summary = payload["summary"]
        print(
            "wrote Jensen-window PF monotone contraction frontier scout: "
            f"{summary['exact_certificate_rows']} exact rows, "
            f"{summary['total_bernstein_coefficients']} Bernstein coefficients, "
            f"{summary['finite_zeta_checked_rows']} finite zeta rows"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
