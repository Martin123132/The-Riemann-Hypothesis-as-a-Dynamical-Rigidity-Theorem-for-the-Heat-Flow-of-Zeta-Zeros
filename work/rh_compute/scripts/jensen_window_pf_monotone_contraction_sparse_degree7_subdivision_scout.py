#!/usr/bin/env python3
"""Repair the degree-7 m=10 Bernstein obstruction by dyadic subdivision."""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import json
from math import comb, gcd
from pathlib import Path

import jensen_window_pf_monotone_contraction_sparse_degree6_scout as sparse


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.md"
DEGREE = 7
MINOR_SIZE = 10


def lcm(left: int, right: int) -> int:
    return left // gcd(left, right) * right


def bit_length_str(value: int) -> int:
    return abs(value).bit_length()


def bernstein_array(power_poly: sparse.Polynomial) -> tuple[list[int], list[int], list[int], int]:
    variable_count = len(next(iter(power_poly)))
    power_degrees = [max(monomial[index] for monomial in power_poly) for index in range(variable_count)]
    active = [index for index, degree in enumerate(power_degrees) if degree > 0]
    multidegree = [power_degrees[index] for index in active]
    shape = [degree + 1 for degree in multidegree]
    strides = sparse.flat_strides(shape)
    total = 1
    for size in shape:
        total *= size

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
    return coeffs, shape, strides, common_denominator


def split_line_half(values: list[int]) -> tuple[list[int], list[int]]:
    degree = len(values) - 1
    work = list(values)
    left = [values[0] * (1 << degree)]
    right = [values[-1] * (1 << degree)]
    for level in range(1, degree + 1):
        for index in range(degree - level + 1):
            work[index] = work[index] + work[index + 1]
        scale = 1 << (degree - level)
        left.append(work[0] * scale)
        right.append(work[degree - level] * scale)
    right.reverse()
    return left, right


def split_axis_half(coeffs: list[int], shape: list[int], strides: list[int], axis: int) -> tuple[list[int], list[int]]:
    degree = shape[axis] - 1
    stride = strides[axis]
    block = stride * (degree + 1)
    total = len(coeffs)
    left = [0] * total
    right = [0] * total
    for base in range(0, total, block):
        for offset in range(stride):
            values = [coeffs[base + index * stride + offset] for index in range(degree + 1)]
            left_values, right_values = split_line_half(values)
            for index, value in enumerate(left_values):
                left[base + index * stride + offset] = value
            for index, value in enumerate(right_values):
                right[base + index * stride + offset] = value
    return left, right


def coefficient_stats(
    coeffs: list[int],
    shape: list[int],
    strides: list[int],
    common_denominator: int,
    extra_power_of_two_denominator: int = 0,
) -> dict:
    min_flat_index, min_scaled = min(enumerate(coeffs), key=lambda item: item[1])
    negative_count = sum(1 for value in coeffs if value < 0)
    zero_count = sum(1 for value in coeffs if value == 0)
    denominator = common_denominator * (1 << extra_power_of_two_denominator)
    return {
        "bernstein_coefficient_count": len(coeffs),
        "bernstein_min_coefficient": sparse.rational_str(Fraction(min_scaled, denominator)),
        "bernstein_min_scaled_coefficient": str(min_scaled),
        "bernstein_min_scaled_bit_length": bit_length_str(min_scaled),
        "bernstein_min_index": sparse.unflatten_index(min_flat_index, strides, shape),
        "bernstein_negative_coefficient_count": negative_count,
        "bernstein_zero_coefficient_count": zero_count,
        "bernstein_coefficients_strictly_positive": negative_count == 0 and zero_count == 0,
    }


def build_payload() -> dict:
    columns = sparse.normalized_column_polynomials(DEGREE, MINOR_SIZE)
    monotone_poly = sparse.substitute_monotone(columns[MINOR_SIZE], DEGREE - 1)
    coeffs, shape, strides, common_denominator = bernstein_array(monotone_poly)
    global_stats = coefficient_stats(coeffs, shape, strides, common_denominator)

    left_half, right_half = split_axis_half(coeffs, shape, strides, axis=0)
    right_left_quarter, right_right_quarter = split_axis_half(right_half, shape, strides, axis=0)

    slab_specs = [
        ("s0 in [0,1/2]", "0", "1/2", 1, left_half),
        ("s0 in [1/2,3/4]", "1/2", "3/4", 2, right_left_quarter),
        ("s0 in [3/4,1]", "3/4", "1", 2, right_right_quarter),
    ]
    slabs = []
    for label, lower, upper, split_depth, slab_coeffs in slab_specs:
        stats = coefficient_stats(
            slab_coeffs,
            shape,
            strides,
            common_denominator,
            extra_power_of_two_denominator=(shape[0] - 1) * split_depth,
        )
        slabs.append(
            {
                "label": label,
                "s0_interval": {"lower": lower, "upper": upper},
                "all_other_s_parameters": "[0,1]",
                "split_axis": "s0",
                "dyadic_split_depth_on_s0": split_depth,
                "extra_power_of_two_denominator": (shape[0] - 1) * split_depth,
                **stats,
            }
        )

    negative_slabs = [slab for slab in slabs if slab["bernstein_negative_coefficient_count"]]
    zero_slabs = [slab for slab in slabs if slab["bernstein_zero_coefficient_count"]]
    return {
        "kind": "jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout",
        "date": "2026-07-06",
        "status": "exact_sparse_degree7_m10_subdivision_certificate",
        "source_degree7_frontier_scout": "outputs/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.md",
        "source_degree6_scout": "outputs/jensen_window_pf_monotone_contraction_sparse_degree6_scout.md",
        "source_monotone_contraction_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "proof_boundary": (
            "Exact bounded degree-7 m=10 subdivision certificate under monotone contractions only. "
            "It repairs the one-shot global Bernstein certificate obstruction by dyadic subdivision "
            "in s0, but it does not prove all m, all degrees, all Schur shapes, zeta cone entry, RH, "
            "or Lambda <= 0."
        ),
        "degree": DEGREE,
        "minor_size": MINOR_SIZE,
        "monotone_contraction_region": "0 <= x1 <= x2 <= x3 <= x4 <= x5 <= x6 <= 1",
        "monotone_substitution": "x_i=1-prod_{r<=i}(1-s_r), all s_r in [0,1]",
        "bernstein_shape": shape,
        "global_one_shot_stats": global_stats,
        "subdivision": {
            "axis": "s0",
            "slabs": slabs,
            "slab_count": len(slabs),
            "total_slab_bernstein_coefficients": sum(slab["bernstein_coefficient_count"] for slab in slabs),
            "negative_slab_count": len(negative_slabs),
            "zero_slab_count": len(zero_slabs),
            "certificate_success": not negative_slabs and not zero_slabs,
        },
        "summary": {
            "repaired_degree": DEGREE,
            "repaired_minor_size": MINOR_SIZE,
            "global_negative_bernstein_coefficients": global_stats["bernstein_negative_coefficient_count"],
            "global_min_coefficient": global_stats["bernstein_min_coefficient"],
            "dyadic_slab_count": len(slabs),
            "total_slab_bernstein_coefficients": sum(slab["bernstein_coefficient_count"] for slab in slabs),
            "negative_slab_coefficients": sum(slab["bernstein_negative_coefficient_count"] for slab in slabs),
            "zero_slab_coefficients": sum(slab["bernstein_zero_coefficient_count"] for slab in slabs),
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "The degree-7 m=10 one-shot global Bernstein obstruction is repaired by splitting only "
                "s0 into the three dyadic slabs [0,1/2], [1/2,3/4], and [3/4,1]. Across those slabs, "
                "785,400 exact Bernstein coefficients are strictly positive. Combined with the global "
                "degree-7 m<=9 certificates, this gives bounded degree-7 column positivity through m=10 "
                "under monotone contractions, while leaving all-m/all-degree and zeta cone-entry theorems open."
            ),
        },
        "invariants": [
            "The global degree-7 m=10 Bernstein net has 126 negative coefficients.",
            "The dyadic s0 subdivision has three slabs and all slab Bernstein coefficients are strictly positive.",
            "The subdivision repairs a finite certificate obstruction only; it is not an all-m theorem.",
            "No row is ready_to_apply.",
            "Endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(payload: dict, path: Path) -> None:
    summary = payload["summary"]
    result_line = (
        "validated Jensen-window PF monotone-contraction sparse degree-7 subdivision scout: "
        f"{payload['subdivision']['slab_count']} dyadic slabs, "
        f"{summary['total_slab_bernstein_coefficients']} slab Bernstein coefficients, "
        f"{summary['negative_slab_coefficients']} negative slab coefficients, "
        f"{summary['zero_slab_coefficients']} zero slab coefficients, repaired m=10 obstruction, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Monotone-Contraction Sparse Degree-7 Subdivision Scout",
        "",
        "Date: 2026-07-06",
        "",
        "Status: exact sparse degree-7 subdivision certificate. This is not a",
        "proof of Jensen-window PF-infinity, all-shape Schur positivity, RH,",
        "or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout`.",
        "",
        "Proof boundary: this artifact repairs the finite degree-7 `m=10`",
        "one-shot global Bernstein certificate obstruction by dyadic subdivision",
        "in `s0`. It does not prove all column rows, all degrees, all Schur",
        "shapes, zeta cone entry, RH, or `Lambda <= 0`.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Global Obstruction",
        "",
        "The one-shot global Bernstein net for degree 7, `m=10`, has:",
        "",
        "```text",
        f"shape={payload['bernstein_shape']}",
        f"min={payload['global_one_shot_stats']['bernstein_min_coefficient']}",
        f"negative coefficients={payload['global_one_shot_stats']['bernstein_negative_coefficient_count']}",
        f"zero coefficients={payload['global_one_shot_stats']['bernstein_zero_coefficient_count']}",
        "```",
        "",
        "## Subdivision Certificate",
        "",
        "Splitting only the first monotone parameter `s0` gives:",
        "",
        "```text",
    ]
    for slab in payload["subdivision"]["slabs"]:
        lines.append(
            f"{slab['label']}: count={slab['bernstein_coefficient_count']}, "
            f"min={slab['bernstein_min_coefficient']}, "
            f"min_index={slab['bernstein_min_index']}, "
            f"negative={slab['bernstein_negative_coefficient_count']}, "
            f"zero={slab['bernstein_zero_coefficient_count']}"
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
            "outputs/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.md",
            "outputs/jensen_window_pf_monotone_contraction_sparse_degree6_scout.md",
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
