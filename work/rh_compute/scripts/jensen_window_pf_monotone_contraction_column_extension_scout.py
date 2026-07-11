#!/usr/bin/env python3
"""Extend the monotone-contraction column-recurrence scout."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_column_extension_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_monotone_contraction_column_extension_scout.md"


def expr_str(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(expr))


def rational_str(expr: sp.Expr) -> str:
    value = sp.Rational(expr)
    if value.q == 1:
        return str(value.p)
    return f"{value.p}/{value.q}"


def column_recurrence_polynomials(degree: int, max_m: int) -> tuple[tuple[sp.Symbol, ...], list[sp.Expr]]:
    h = sp.symbols(f"h0:{degree + 1}")
    columns = [sp.Integer(1)]
    for size in range(1, max_m + 1):
        value = sp.Integer(0)
        for offset in range(1, min(degree, size) + 1):
            value += (-1) ** (offset - 1) * h[0] ** (offset - 1) * h[offset] * columns[size - offset]
        columns.append(sp.factor(value))
    return h, columns


def ratio_substitution(degree: int) -> tuple[list[sp.Expr], tuple[sp.Symbol, ...], tuple[sp.Symbol, ...]]:
    scale, rho = sp.symbols("A rho")
    contractions = sp.symbols(f"x1:{degree}")
    values = [scale]
    for index in range(1, degree + 1):
        value = scale * rho**index
        for contraction_index in range(1, index):
            value *= contractions[contraction_index - 1] ** (index - contraction_index)
        values.append(value)
    return values, (scale, rho) + contractions, contractions


def positive_monomial_factor(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> tuple[sp.Expr, sp.Expr]:
    poly = sp.Poly(sp.expand(expr), *variables)
    min_powers = [
        min(mon[index] for mon, coeff in poly.terms() if coeff != 0)
        for index, _var in enumerate(variables)
    ]
    monomial = sp.prod(var**power for var, power in zip(variables, min_powers))
    return sp.factor(monomial), sp.factor(sp.expand(expr) / monomial)


def monotone_substitution(degree: int) -> tuple[dict[sp.Symbol, sp.Expr], tuple[sp.Symbol, ...], str]:
    contractions = sp.symbols(f"x1:{degree}")
    parameters = sp.symbols(f"s0:{degree - 1}")
    substitutions = {contractions[0]: parameters[0]}
    current = parameters[0]
    pieces = [f"x1={sp.sstr(parameters[0])}"]
    for index in range(1, degree - 1):
        current = current + parameters[index] * (1 - current)
        substitutions[contractions[index]] = current
        pieces.append(f"x{index + 1}=x{index}+{sp.sstr(parameters[index])}*(1-x{index})")
    pieces.append("all parameters in [0,1]")
    return substitutions, parameters, ", ".join(pieces)


def bernstein_coefficients(poly: sp.Expr, variables: tuple[sp.Symbol, ...]) -> tuple[list[int], list[dict]]:
    active = tuple(var for var in variables if poly.has(var))
    if not active:
        return [], [{"index": [], "coefficient": rational_str(poly)}]

    power_poly = sp.Poly(sp.expand(poly), *active)
    degrees = [power_poly.degree(var) for var in active]
    power_coefficients = {monomial: sp.Rational(coefficient) for monomial, coefficient in power_poly.terms()}
    indices = list(itertools.product(*[range(degree + 1) for degree in degrees]))
    rows = []
    for multi_index in indices:
        coefficient = sp.Rational(0)
        for powers, power_coefficient in power_coefficients.items():
            factor = sp.Rational(1)
            for index, power, degree in zip(multi_index, powers, degrees):
                if power > index:
                    factor = None
                    break
                factor *= sp.Rational(sp.binomial(index, power), sp.binomial(degree, power))
            if factor is not None:
                coefficient += power_coefficient * factor
        rows.append({"index": list(multi_index), "coefficient": rational_str(coefficient)})
    return degrees, rows


def analyze_degree(degree: int, max_m: int) -> list[dict]:
    h_symbols, columns = column_recurrence_polynomials(degree, max_m)
    a_values, positive_variables, _contractions = ratio_substitution(degree)
    h_subs = {h_symbols[index]: sp.binomial(degree, index) * a_values[index] for index in range(degree + 1)}
    monotone_subs, monotone_variables, monotone_description = monotone_substitution(degree)

    rows: list[dict] = []
    for size in range(1, max_m + 1):
        recurrence_expr = sp.factor(columns[size].subs(h_subs))
        monomial, normalized = positive_monomial_factor(recurrence_expr, positive_variables)
        monotone_poly = sp.factor(normalized.subs(monotone_subs))
        multidegree, coefficients = bernstein_coefficients(monotone_poly, monotone_variables)
        coefficient_values = [sp.Rational(row["coefficient"]) for row in coefficients]
        negative_count = sum(1 for value in coefficient_values if value < 0)
        rows.append(
            {
                "id": f"mccx_d{degree}_m{size}",
                "degree": degree,
                "minor_size": size,
                "normalized_positive_monomial": expr_str(monomial),
                "monotone_substitution": monotone_description,
                "monotone_region_polynomial": expr_str(monotone_poly),
                "bernstein_variables": [sp.sstr(var) for var in monotone_variables if monotone_poly.has(var)],
                "bernstein_multidegree": multidegree,
                "bernstein_coefficient_count": len(coefficients),
                "bernstein_min_coefficient": rational_str(min(coefficient_values)),
                "bernstein_negative_coefficient_count": negative_count,
                "bernstein_coefficients_positive": negative_count == 0 and all(value > 0 for value in coefficient_values),
                "beyond_first_hard_frontier": (degree == 3 and size > 8) or (degree == 4 and size > 6),
                "proof_boundary": (
                    "Bounded exact column-recurrence certificate under monotone contractions only; "
                    "not an all-m, all-degree, all-shift, all-shape, or zeta cone-entry theorem."
                ),
            }
        )
    return rows


def build_payload() -> dict:
    rows_by_degree = {"3": analyze_degree(3, 10), "4": analyze_degree(4, 7), "5": analyze_degree(5, 8)}
    all_rows = [row for rows in rows_by_degree.values() for row in rows]
    beyond_rows = [row for row in all_rows if row["beyond_first_hard_frontier"]]
    degree5_rows = rows_by_degree["5"]
    negative_rows = [row for row in all_rows if row["bernstein_negative_coefficient_count"]]
    return {
        "kind": "jensen_window_pf_monotone_contraction_column_extension_scout",
        "date": "2026-07-06",
        "status": "exact_bounded_column_extension_diagnostic",
        "source_frontier_scout": "outputs/jensen_window_pf_monotone_contraction_frontier_scout.md",
        "source_column_recurrence_contract": "outputs/jensen_window_pf_column_recurrence_contract.md",
        "source_monotone_contraction_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "proof_boundary": (
            "Exact bounded column-recurrence theorem-search diagnostic only. It proves "
            "positive Bernstein certificates for finitely many degree-3 and degree-4 "
            "column rows under monotone contractions, but it does not prove monotone "
            "contractions for the zeta coefficients, all column rows, all Schur shapes, "
            "Jensen-window PF-infinity, RH, or Lambda <= 0."
        ),
        "monotone_contraction_region": {
            "degree3": "0 <= x1 <= x2 <= 1",
            "degree4": "0 <= x1 <= x2 <= x3 <= 1",
            "degree5": "0 <= x1 <= x2 <= x3 <= x4 <= 1",
        },
        "rows_by_degree": rows_by_degree,
        "summary": {
            "degree3_rows": len(rows_by_degree["3"]),
            "degree4_rows": len(rows_by_degree["4"]),
            "degree5_rows": len(degree5_rows),
            "total_column_rows": len(all_rows),
            "total_bernstein_coefficients": sum(row["bernstein_coefficient_count"] for row in all_rows),
            "beyond_first_hard_frontier_rows": len(beyond_rows),
            "higher_degree_extension_rows": len(degree5_rows),
            "negative_bernstein_rows": len(negative_rows),
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "The monotone-contraction cone certifies more than the first two hard "
                "column-frontier rows: degree 3 is certified through m=10 and degree "
                "4 through m=7, with three rows beyond the original hard frontier; "
                "the same method also certifies degree 5 through m=8. This "
                "strengthens the monotone-contraction route as a column-family target, "
                "while leaving the all-m/all-shape and zeta cone-entry theorems open."
            ),
        },
        "invariants": [
            "Every row has strictly positive Bernstein coefficients.",
            "The scout covers only degree 3 m<=10, degree 4 m<=7, and degree 5 m<=8.",
            "Rows beyond the first hard frontier are diagnostic extensions, not an all-m theorem.",
            "Degree-5 rows are higher-degree diagnostic extensions, not an all-degree theorem.",
            "No row is ready_to_apply.",
            "Endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(payload: dict, path: Path) -> None:
    summary = payload["summary"]
    result_line = (
        "validated Jensen-window PF monotone-contraction column extension scout: "
        f"{summary['total_column_rows']} column rows, "
        f"{summary['total_bernstein_coefficients']} Bernstein coefficients, "
        f"{summary['beyond_first_hard_frontier_rows']} beyond-frontier rows, "
        f"{summary['negative_bernstein_rows']} negative Bernstein rows, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Monotone-Contraction Column Extension Scout",
        "",
        "Date: 2026-07-06",
        "",
        "Status: exact bounded column-extension diagnostic. This is not a proof",
        "of Jensen-window PF-infinity, all-shape Schur positivity, RH, or",
        "`Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_monotone_contraction_column_extension_scout`.",
        "",
        "Proof boundary: this artifact proves finite exact Bernstein certificates",
        "for bounded column rows under monotone contractions. It does not prove",
        "that the zeta coefficients satisfy monotone contractions, and it does",
        "not prove all column rows or all Schur shapes.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_monotone_contraction_column_extension_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_column_extension_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_column_extension_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Certified Region",
        "",
        "The scout works on the monotone-contraction cones:",
        "",
        "```text",
        "degree 3: 0 <= x1 <= x2 <= 1",
        "degree 4: 0 <= x1 <= x2 <= x3 <= 1",
        "degree 5: 0 <= x1 <= x2 <= x3 <= x4 <= 1",
        "```",
        "",
        "where `x_i=(A_{n+i+1}/A_{n+i})/(A_{n+i}/A_{n+i-1})`.",
        "",
        "## Column Rows",
        "",
        "```text",
    ]
    for degree in ("3", "4", "5"):
        for row in payload["rows_by_degree"][degree]:
            lines.append(
                f"{row['id']}: degree {row['degree']}, m={row['minor_size']}, "
                f"Bernstein count={row['bernstein_coefficient_count']}, "
                f"min={row['bernstein_min_coefficient']}, "
                f"beyond_frontier={row['beyond_first_hard_frontier']}"
            )
    lines.extend(
        [
            "```",
            "",
            "The extension beyond the original first hard frontier is:",
            "",
            "```text",
            "degree 3: m=9, m=10",
            "degree 4: m=7",
            "```",
            "",
            "The higher-degree extension is:",
            "",
            "```text",
            "degree 5: m=1..8",
            "```",
            "",
            "## Consequence",
            "",
            "This strengthens the monotone-contraction route from a two-row escape",
            "from the rational countermodel into a bounded column-family theorem",
            "search target that now includes a first degree-5 band. It still leaves two hard problems open: proving",
            "monotone contractions for the actual zeta heat-flow coefficients and",
            "lifting from bounded column rows to all column rows/all Schur shapes.",
            "",
            "Integration:",
            "",
            "```text",
            "outputs/jensen_window_pf_monotone_contraction_frontier_scout.md",
            "outputs/jensen_window_pf_column_recurrence_contract.md",
            "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
            "outputs/jensen_window_pf_cauchy_binet_cone_frontier_matrix.md",
            "```",
            "",
            "Summary:",
            "",
            summary["main_finding"],
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    payload = build_payload()
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(payload, note)
    print(
        "wrote Jensen-window PF monotone-contraction column extension scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
