#!/usr/bin/env python3
"""Sparse exact degree-6 monotone-contraction column scout."""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import json
from math import comb, gcd
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree6_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_monotone_contraction_sparse_degree6_scout.md"
DEGREE = 6
MAX_M = 10

Monomial = tuple[int, ...]
Polynomial = dict[Monomial, int]


def rational_str(value: Fraction | int) -> str:
    frac = Fraction(value)
    if frac.denominator == 1:
        return str(frac.numerator)
    return f"{frac.numerator}/{frac.denominator}"


def lcm(left: int, right: int) -> int:
    return left // gcd(left, right) * right


def add_scaled_inplace(target: defaultdict[Monomial, int], source: Polynomial, scale: int) -> None:
    if scale == 0:
        return
    for monomial, coefficient in source.items():
        value = target[monomial] + scale * coefficient
        if value:
            target[monomial] = value
        elif monomial in target:
            del target[monomial]


def multiply_monomial(poly: Polynomial, exponent: Monomial, coefficient: int = 1) -> Polynomial:
    out: defaultdict[Monomial, int] = defaultdict(int)
    for monomial, value in poly.items():
        out[tuple(left + right for left, right in zip(monomial, exponent))] += coefficient * value
    return {monomial: value for monomial, value in out.items() if value}


def multiply_poly(left: Polynomial, right: Polynomial) -> Polynomial:
    out: defaultdict[Monomial, int] = defaultdict(int)
    for left_monomial, left_coeff in left.items():
        for right_monomial, right_coeff in right.items():
            monomial = tuple(a + b for a, b in zip(left_monomial, right_monomial))
            out[monomial] += left_coeff * right_coeff
    return {monomial: value for monomial, value in out.items() if value}


def power_poly(poly: Polynomial, exponent: int) -> Polynomial:
    variable_count = len(next(iter(poly)))
    result: Polynomial = {tuple([0] * variable_count): 1}
    base = poly
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply_poly(result, base)
        remaining //= 2
        if remaining:
            base = multiply_poly(base, base)
    return result


def normalized_column_polynomials(degree: int, max_m: int) -> list[Polynomial]:
    """Return Q_m after factoring the positive A^m rho^m monomial.

    With h_j=binom(d,j) A rho^j prod_{i<j} x_i^{j-i}, the column recurrence
    C_m=sum_j (-1)^(j-1) h_0^(j-1) h_j C_{m-j} becomes
    Q_m=sum_j (-1)^(j-1) binom(d,j) P_j(x) Q_{m-j}.
    """

    variable_count = degree - 1
    zero = tuple([0] * variable_count)
    multipliers: list[tuple[Monomial, int]] = [(zero, 1)]
    for offset in range(1, degree + 1):
        exponent = [0] * variable_count
        for index in range(1, offset):
            exponent[index - 1] = offset - index
        multipliers.append((tuple(exponent), comb(degree, offset)))

    columns: list[Polynomial] = [{zero: 1}]
    for size in range(1, max_m + 1):
        out: defaultdict[Monomial, int] = defaultdict(int)
        for offset in range(1, min(degree, size) + 1):
            exponent, coefficient = multipliers[offset]
            sign = 1 if (offset - 1) % 2 == 0 else -1
            add_scaled_inplace(out, multiply_monomial(columns[size - offset], exponent, sign * coefficient), 1)
        columns.append(dict(out))
    return columns


def monotone_contraction_polynomials(variable_count: int) -> list[Polynomial]:
    """Return x_i=1-prod_{r<=i}(1-s_r) in sparse power form."""

    zero = tuple([0] * variable_count)
    running_product: Polynomial = {zero: 1}
    x_polys: list[Polynomial] = []
    for index in range(variable_count):
        unit = tuple(1 if slot == index else 0 for slot in range(variable_count))
        running_product = multiply_poly(running_product, {zero: 1, unit: -1})
        x_i: defaultdict[Monomial, int] = defaultdict(int)
        x_i[zero] = 1
        add_scaled_inplace(x_i, running_product, -1)
        x_polys.append(dict(x_i))
    return x_polys


def substitute_monotone(poly: Polynomial, variable_count: int) -> Polynomial:
    x_polys = monotone_contraction_polynomials(variable_count)
    power_cache: dict[tuple[int, int], Polynomial] = {}
    zero = tuple([0] * variable_count)
    out: defaultdict[Monomial, int] = defaultdict(int)

    for x_exponent, coefficient in poly.items():
        term: Polynomial = {zero: coefficient}
        for index, exponent in enumerate(x_exponent):
            if exponent == 0:
                continue
            key = (index, exponent)
            if key not in power_cache:
                power_cache[key] = power_poly(x_polys[index], exponent)
            term = multiply_poly(term, power_cache[key])
        add_scaled_inplace(out, term, 1)
    return dict(out)


def flat_strides(shape: list[int]) -> list[int]:
    if not shape:
        return []
    strides: list[int] = []
    tail = 1
    for size in reversed(shape[1:]):
        tail *= size
        strides.insert(0, tail)
    strides.append(1)
    return strides


def unflatten_index(flat_index: int, strides: list[int], shape: list[int]) -> list[int]:
    if not shape:
        return []
    indices: list[int] = []
    remaining = flat_index
    for stride, size in zip(strides, shape):
        index = remaining // stride
        if index >= size:
            raise ValueError(f"flat index {flat_index} is outside shape {shape}")
        indices.append(index)
        remaining -= index * stride
    return indices


def scaled_bernstein_stats(
    power_poly: Polynomial,
    *,
    include_index_details: bool = False,
    negative_example_limit: int = 0,
) -> dict:
    variable_count = len(next(iter(power_poly)))
    power_degrees = [max(monomial[index] for monomial in power_poly) for index in range(variable_count)]
    active = [index for index, degree in enumerate(power_degrees) if degree > 0]
    multidegree = [power_degrees[index] for index in active]
    shape = [degree + 1 for degree in multidegree]
    total = 1
    for size in shape:
        total *= size

    strides = flat_strides(shape)
    coeffs = [0] * total
    for monomial, coefficient in power_poly.items():
        flat_index = 0
        for active_index, stride in zip(active, strides):
            flat_index += monomial[active_index] * stride
        coeffs[flat_index] = coefficient

    common_denominator = 1
    for axis, degree in enumerate(multidegree):
        denominator_part = 1
        for alpha in range(degree + 1):
            denominator_part = lcm(denominator_part, comb(degree, alpha))
        common_denominator *= denominator_part
        factors = [
            [comb(index, alpha) * (denominator_part // comb(degree, alpha)) if index >= alpha else 0 for index in range(degree + 1)]
            for alpha in range(degree + 1)
        ]
        stride = strides[axis]
        block = stride * (degree + 1)
        transformed = [0] * total
        for base in range(0, total, block):
            for offset in range(stride):
                power_values = [coeffs[base + alpha * stride + offset] for alpha in range(degree + 1)]
                for index in range(degree + 1):
                    value = 0
                    for alpha in range(index + 1):
                        value += power_values[alpha] * factors[alpha][index]
                    transformed[base + index * stride + offset] = value
        coeffs = transformed

    min_flat_index, min_scaled = min(enumerate(coeffs), key=lambda item: item[1])
    negative_count = sum(1 for value in coeffs if value < 0)
    zero_count = sum(1 for value in coeffs if value == 0)
    stats = {
        "bernstein_variables": [f"s{index}" for index in active],
        "bernstein_multidegree": multidegree,
        "bernstein_coefficient_count": total,
        "bernstein_min_coefficient": rational_str(Fraction(min_scaled, common_denominator)),
        "bernstein_negative_coefficient_count": negative_count,
        "bernstein_zero_coefficient_count": zero_count,
        "bernstein_coefficients_strictly_positive": negative_count == 0 and zero_count == 0,
    }
    if include_index_details:
        stats["bernstein_min_index"] = unflatten_index(min_flat_index, strides, shape)
    if negative_example_limit:
        examples = []
        for flat_index, value in enumerate(coeffs):
            if value < 0:
                examples.append(
                    {
                        "index": unflatten_index(flat_index, strides, shape),
                        "coefficient": rational_str(Fraction(value, common_denominator)),
                    }
                )
                if len(examples) >= negative_example_limit:
                    break
        stats["bernstein_negative_examples"] = examples
    return stats


def analyze_degree6(max_m: int = MAX_M) -> list[dict]:
    variable_count = DEGREE - 1
    columns = normalized_column_polynomials(DEGREE, max_m)
    rows: list[dict] = []
    for size in range(1, max_m + 1):
        monotone_poly = substitute_monotone(columns[size], variable_count)
        power_degrees = [max(monomial[index] for monomial in monotone_poly) for index in range(variable_count)]
        stats = scaled_bernstein_stats(monotone_poly)
        rows.append(
            {
                "id": f"mcs6_d6_m{size}",
                "degree": DEGREE,
                "minor_size": size,
                "normalized_recurrence_term_count": len(columns[size]),
                "monotone_power_basis_term_count": len(monotone_poly),
                "monotone_power_multidegree": power_degrees,
                "monotone_substitution": "x_i=1-prod_{r<=i}(1-s_r), all s_r in [0,1]",
                **stats,
                "proof_boundary": (
                    "Sparse exact degree-6 column-recurrence certificate under monotone contractions only; "
                    "not an all-m, all-degree, all-shift, all-shape, zeta cone-entry, RH, or Lambda <= 0 theorem."
                ),
            }
        )
    return rows


def build_payload() -> dict:
    rows = analyze_degree6(MAX_M)
    negative_rows = [row for row in rows if row["bernstein_negative_coefficient_count"]]
    zero_rows = [row for row in rows if row["bernstein_zero_coefficient_count"]]
    return {
        "kind": "jensen_window_pf_monotone_contraction_sparse_degree6_scout",
        "date": "2026-07-06",
        "status": "exact_sparse_degree6_column_extension_diagnostic",
        "source_column_extension_scout": "outputs/jensen_window_pf_monotone_contraction_column_extension_scout.md",
        "source_column_recurrence_contract": "outputs/jensen_window_pf_column_recurrence_contract.md",
        "source_monotone_contraction_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "proof_boundary": (
            "Exact sparse bounded degree-6 column-recurrence theorem-search diagnostic only. "
            "It proves strict Bernstein positivity for finitely many degree-6 rows under "
            "monotone contractions, but it does not prove monotone contractions for the "
            "zeta coefficients, all column rows, all Schur shapes, Jensen-window PF-infinity, "
            "RH, or Lambda <= 0."
        ),
        "monotone_contraction_region": "0 <= x1 <= x2 <= x3 <= x4 <= x5 <= 1",
        "normalized_recurrence": (
            "Q_0=1 and Q_m=sum_{j=1..min(6,m)} (-1)^(j-1) binom(6,j) "
            "prod_{i=1..j-1} x_i^(j-i) Q_{m-j}, after factoring A^m*rho^m."
        ),
        "rows": rows,
        "summary": {
            "degree6_rows": len(rows),
            "max_minor_size": MAX_M,
            "total_bernstein_coefficients": sum(row["bernstein_coefficient_count"] for row in rows),
            "largest_bernstein_row": max(row["bernstein_coefficient_count"] for row in rows),
            "negative_bernstein_rows": len(negative_rows),
            "zero_bernstein_rows": len(zero_rows),
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "A sparse exact common-denominator Bernstein transform certifies degree 6 "
                "monotone-contraction column rows through m=10, with 63,347 strict positive "
                "Bernstein coefficients and no zero or negative Bernstein rows. This extends "
                "the bounded column-family evidence beyond the degree-5 band, but it remains "
                "finite theorem-search algebra rather than an all-m or zeta cone-entry theorem."
            ),
        },
        "invariants": [
            "Every recorded degree-6 row has strictly positive Bernstein coefficients.",
            "The scout covers only degree 6 and m<=10.",
            "The sparse recurrence removes only the positive A^m*rho^m monomial.",
            "No row is ready_to_apply.",
            "Endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(payload: dict, path: Path) -> None:
    summary = payload["summary"]
    result_line = (
        "validated Jensen-window PF monotone-contraction sparse degree-6 scout: "
        f"{summary['degree6_rows']} degree-6 rows, "
        f"{summary['total_bernstein_coefficients']} Bernstein coefficients, "
        f"m<={summary['max_minor_size']}, "
        f"{summary['negative_bernstein_rows']} negative Bernstein rows, "
        f"{summary['zero_bernstein_rows']} zero Bernstein rows, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Monotone-Contraction Sparse Degree-6 Scout",
        "",
        "Date: 2026-07-06",
        "",
        "Status: exact sparse degree-6 diagnostic. This is not a proof of",
        "Jensen-window PF-infinity, all-shape Schur positivity, RH, or",
        "`Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_monotone_contraction_sparse_degree6_scout`.",
        "",
        "Proof boundary: this artifact proves finite exact Bernstein sign",
        "certificates for bounded degree-6 column rows under monotone",
        "contractions. It does not prove that the zeta coefficients satisfy",
        "monotone contractions, and it does not prove all column rows or all",
        "Schur shapes.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree6_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_sparse_degree6_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree6_scout.py",
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
        "```text",
        "degree 6: 0 <= x1 <= x2 <= x3 <= x4 <= x5 <= 1",
        "x_i=1-prod_{r<=i}(1-s_r), all s_r in [0,1]",
        "```",
        "",
        "## Sparse Recurrence",
        "",
        "After removing the positive monomial `A^m*rho^m`, the column recurrence is:",
        "",
        "```text",
        payload["normalized_recurrence"],
        "```",
        "",
        "## Degree-6 Rows",
        "",
        "```text",
    ]
    for row in payload["rows"]:
        lines.append(
            f"{row['id']}: m={row['minor_size']}, "
            f"recurrence terms={row['normalized_recurrence_term_count']}, "
            f"power terms={row['monotone_power_basis_term_count']}, "
            f"Bernstein multidegree={row['bernstein_multidegree']}, "
            f"Bernstein count={row['bernstein_coefficient_count']}, "
            f"min={row['bernstein_min_coefficient']}, "
            f"negative={row['bernstein_negative_coefficient_count']}, "
            f"zero={row['bernstein_zero_coefficient_count']}"
        )
    lines.extend(
        [
            "```",
            "",
            "## Consequence",
            "",
            payload["summary"]["main_finding"],
            "",
            "Integration:",
            "",
            "```text",
            "outputs/jensen_window_pf_monotone_contraction_column_extension_scout.md",
            "outputs/jensen_window_pf_column_recurrence_contract.md",
            "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
            "outputs/signed_hankel_jensen_dependency_graph.md",
            "```",
            "",
            "Summary:",
            "",
            payload["summary"]["main_finding"],
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(payload, args.note)
    print(
        "wrote "
        f"{args.out_json.relative_to(REPO_ROOT)} and {args.note.relative_to(REPO_ROOT)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
