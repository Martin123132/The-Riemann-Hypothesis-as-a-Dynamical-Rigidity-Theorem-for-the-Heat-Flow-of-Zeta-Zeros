#!/usr/bin/env python3
"""Test contraction log-concavity as a Jensen-window PF bridge condition.

The ratio-condition scout left x2**2 >= x1*x3 as the only natural nearby
condition not rejected by the original rational countermodel.  This script
keeps the same bad degree-3 prefix and chooses a smaller positive x3 so that
the contraction log-concavity condition is satisfied exactly.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ALGEBRA = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_obligation_algebra.json"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_contraction_log_concavity_scout.json"


def rational_str(expr: sp.Expr) -> str:
    value = sp.Rational(expr)
    if value.q == 1:
        return str(value.p)
    return f"{value.p}/{value.q}"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def evaluate_contiguous(values: list[sp.Rational], degree: int, size: int) -> sp.Rational:
    window = [sp.binomial(degree, index) * values[index] for index in range(degree + 1)]
    matrix = [
        [window[col - row] if 0 <= col - row < len(window) else sp.Integer(0) for col in range(1, size + 1)]
        for row in range(size)
    ]
    return sp.factor(sp.Matrix(matrix).det())


def constructed_extension() -> dict:
    rho = sp.Rational(33, 40)
    x1 = sp.Rational(65, 66)
    x2 = sp.Rational(44, 65)
    x3 = sp.Rational(1, 3)

    values = [
        sp.Rational(1),
        rho,
        rho**2 * x1,
        rho**3 * x1**2 * x2,
        rho**4 * x1**3 * x2**2 * x3,
    ]
    ratios = [sp.factor(values[index] / values[index - 1]) for index in range(1, len(values))]
    contractions = [
        sp.factor(ratios[index] / ratios[index - 1])
        for index in range(1, len(ratios))
    ]
    condition_values = {
        "1 - x1": 1 - x1,
        "1 - x2": 1 - x2,
        "1 - x3": 1 - x3,
        "x1 - x2": x1 - x2,
        "x2 - x3": x2 - x3,
        "(1 - x2)^2 - x2^2*(1 - x1)*(1 - x3)": (1 - x2) ** 2 - x2**2 * (1 - x1) * (1 - x3),
        "x2^2 - x1*x3": x2**2 - x1 * x3,
    }
    return {
        "rho": rho,
        "x1": x1,
        "x2": x2,
        "x3": x3,
        "values": values,
        "ratios": ratios,
        "contractions": contractions,
        "condition_values": condition_values,
    }


def bool_conditions(condition_values: dict[str, sp.Rational]) -> dict[str, bool]:
    return {
        "adjacent_log_concavity_box": bool(
            all(condition_values[key] >= 0 for key in ("1 - x1", "1 - x2", "1 - x3"))
        ),
        "decreasing_ratio_contractions": bool(
            all(condition_values[key] >= 0 for key in ("x1 - x2", "x2 - x3"))
        ),
        "second_order_log_concavity": bool(condition_values["(1 - x2)^2 - x2^2*(1 - x1)*(1 - x3)"] >= 0),
        "contraction_log_concavity": bool(condition_values["x2^2 - x1*x3"] >= 0),
    }


def frontier_values(values: list[sp.Rational]) -> dict:
    rows = []
    for degree, size in ((3, 8), (4, 6)):
        value = evaluate_contiguous(values, degree, size)
        rows.append(
            {
                "id": f"d{degree}_contiguous_m{size}",
                "degree": degree,
                "minor_size": size,
                "value": rational_str(value),
                "sign": int(sp.sign(value)),
                "positive": bool(value > 0),
            }
        )
    return {
        "rows": rows,
        "negative_rows": sum(1 for row in rows if row["sign"] < 0),
    }


def build_payload(algebra: dict) -> dict:
    point = constructed_extension()
    frontiers = frontier_values(point["values"])
    original_d3 = algebra["finite_countermodel"]["d3_first_negative_contiguous_toeplitz_minor"]["determinant"]
    return {
        "kind": "jensen_window_pf_contraction_log_concavity_scout",
        "date": "2026-07-06",
        "target_ratio_condition": "rc_07_contraction_log_concavity",
        "source_algebra": "work/rh_compute/results/jensen_window_pf_obligation_algebra.json",
        "proof_boundary": (
            "Ratio-condition rejection diagnostic only; not a proof of a positive "
            "kernel, Cauchy-Binet identity, Jensen-window PF-infinity, RH, or "
            "Lambda <= 0."
        ),
        "constructed_extension": {
            "ratio_parameterization": {
                "a0": "A",
                "a1": "A*rho",
                "a2": "A*rho**2*x1",
                "a3": "A*rho**3*x1**2*x2",
                "a4": "A*rho**4*x1**3*x2**2*x3",
            },
            "ratio_point": {
                "rho": rational_str(point["rho"]),
                "x1": rational_str(point["x1"]),
                "x2": rational_str(point["x2"]),
                "x3": rational_str(point["x3"]),
            },
            "sequence_A0_to_A4": [rational_str(value) for value in point["values"]],
            "adjacent_ratios": [rational_str(value) for value in point["ratios"]],
            "contractions": [rational_str(value) for value in point["contractions"]],
            "condition_values": {
                key: rational_str(value) for key, value in point["condition_values"].items()
            },
            "satisfies": bool_conditions(point["condition_values"]),
        },
        "frontier_failures": frontiers,
        "comparison_to_original_countermodel": {
            "same_A0_to_A3_prefix": True,
            "d3_m8_matches_original_countermodel": frontiers["rows"][0]["value"] == original_d3,
            "original_d3_m8_value": original_d3,
        },
        "summary": {
            "candidate_rows_tested": 1,
            "rejected_by_constructed_extension_rows": 1,
            "negative_frontier_rows": frontiers["negative_rows"],
            "target_closing_rows": 0,
            "main_finding": (
                "Contraction log-concavity x2^2 >= x1*x3 is satisfied by a "
                "positive rational extension of the old bad prefix, but the "
                "degree 3 size 8 and degree 4 size 6 contiguous minors are "
                "still negative."
            ),
        },
        "invariants": [
            "The constructed extension must satisfy x2^2 >= x1*x3 exactly.",
            "At least one advertised frontier minor must remain negative.",
            "This diagnostic rejects rc_07 only as a standalone bridge condition.",
            "No row may set target_closing=true.",
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
    payload = build_payload(load_json(args.algebra))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        summary = payload["summary"]
        print(
            "wrote Jensen-window PF contraction-log-concavity scout: "
            f"{summary['rejected_by_constructed_extension_rows']} rejected by construction, "
            f"{summary['negative_frontier_rows']} negative frontier rows"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
