#!/usr/bin/env python3
"""Audit strengthened ratio-condition candidates for the Jensen-window PF bridge.

The frontier scouts show that adjacent log-concavity is too weak.  This script
tests a few natural strengthening attempts against the exact rational
countermodel or constructed extensions, and separates rejected candidates from
tautological restatements.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ALGEBRA = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_obligation_algebra.json"
DEFAULT_LOW_SCOUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json"
DEFAULT_FRONTIER = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json"
DEFAULT_CONTRACTION_SCOUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_contraction_log_concavity_scout.json"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_ratio_condition_scout.json"


def rational_str(expr: sp.Expr) -> str:
    value = sp.Rational(expr)
    if value.q == 1:
        return str(value.p)
    return f"{value.p}/{value.q}"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def countermodel_ratios(algebra: dict) -> dict:
    values = [sp.Rational(item) for item in algebra["finite_countermodel"]["sequence_A0_to_A4"]]
    ratios = [sp.factor(values[index] / values[index - 1]) for index in range(1, len(values))]
    contractions = [
        sp.factor(ratios[index] / ratios[index - 1])
        for index in range(1, len(ratios))
    ]
    x1, x2, x3 = contractions
    return {
        "sequence_A0_to_A4": [rational_str(value) for value in values],
        "adjacent_ratios": [rational_str(value) for value in ratios],
        "x1": rational_str(x1),
        "x2": rational_str(x2),
        "x3": rational_str(x3),
        "condition_values": {
            "1 - x1": rational_str(1 - x1),
            "1 - x2": rational_str(1 - x2),
            "1 - x3": rational_str(1 - x3),
            "x1 - x2": rational_str(x1 - x2),
            "x2 - x3": rational_str(x2 - x3),
            "(1 - x2)^2 - x2^2*(1 - x1)*(1 - x3)": rational_str(
                (1 - x2) ** 2 - x2**2 * (1 - x1) * (1 - x3)
            ),
            "x2^2 - x1*x3": rational_str(x2**2 - x1 * x3),
        },
    }


def frontier_failures(frontier: dict) -> list[dict]:
    d3 = frontier["rows_by_degree"]["3"][-1]
    d4 = frontier["rows_by_degree"]["4"][-1]
    return [
        {
            "id": d3["id"],
            "degree": d3["degree"],
            "minor_size": d3["minor_size"],
            "countermodel_value": d3["countermodel_value"],
        },
        {
            "id": d4["id"],
            "degree": d4["degree"],
            "minor_size": d4["minor_size"],
            "countermodel_value": d4["countermodel_value"],
        },
    ]


def row(
    row_id: str,
    title: str,
    status: str,
    condition: str,
    countermodel_satisfies: bool,
    values: dict,
    conclusion: str,
    next_action: str,
    target_failures: list[dict],
) -> dict:
    return {
        "id": row_id,
        "title": title,
        "status": status,
        "condition": condition,
        "countermodel_satisfies": countermodel_satisfies,
        "condition_values_at_countermodel": values,
        "target_failures_at_countermodel": target_failures if countermodel_satisfies else [],
        "target_closing": False,
        "conclusion": conclusion,
        "next_action": next_action,
        "proof_boundary": (
            "Ratio-condition theorem-search diagnostic only; not a proof of "
            "Jensen-window PF-infinity, RH, or Lambda <= 0."
        ),
    }


def build_rows(algebra: dict, low_scout: dict, frontier: dict, contraction_scout: dict, ratios: dict) -> list[dict]:
    values = ratios["condition_values"]
    failures = frontier_failures(frontier)
    contraction_extension = contraction_scout["constructed_extension"]
    contraction_failures = contraction_scout["frontier_failures"]["rows"]
    return [
        row(
            "rc_01_adjacent_log_concavity",
            "adjacent log-concavity",
            "rejected_by_countermodel",
            "0 <= x1,x2,x3 <= 1",
            True,
            {"1 - x1": values["1 - x1"], "1 - x2": values["1 - x2"], "1 - x3": values["1 - x3"]},
            "Adjacent log-concavity is satisfied by the rational countermodel but does not prevent larger contiguous Jensen-window PF failures.",
            "Do not use adjacent log-concavity alone as a bridge lemma.",
            failures,
        ),
        row(
            "rc_02_decreasing_ratio_contractions",
            "decreasing adjacent ratio contractions",
            "rejected_by_countermodel",
            "x1 >= x2 >= x3",
            True,
            {"x1 - x2": values["x1 - x2"], "x2 - x3": values["x2 - x3"]},
            "Monotone contraction of the adjacent ratios is still satisfied by the rational countermodel.",
            "Require a stronger structural condition than monotone contraction variables.",
            failures,
        ),
        row(
            "rc_03_second_order_log_concavity",
            "second-order log-concavity of the original five-term prefix",
            "rejected_by_countermodel",
            "(1 - x2)^2 >= x2^2*(1 - x1)*(1 - x3)",
            True,
            {
                "(1 - x2)^2 - x2^2*(1 - x1)*(1 - x3)": values[
                    "(1 - x2)^2 - x2^2*(1 - x1)*(1 - x3)"
                ]
            },
            "The natural two-log-concavity inequality is positive at the countermodel but does not prevent Jensen-window PF failure.",
            "Do not promote finite higher log-concavity of A_0..A_4 into the missing all-minor theorem.",
            failures,
        ),
        row(
            "rc_04_selected_low_degree_bernstein_positivity",
            "selected low-degree Bernstein positivity",
            "rejected_by_countermodel",
            "all 15 selected degree-2/3/4 hard formulas have nonnegative Bernstein coefficients",
            True,
            {
                "low_scout_formula_rows": str(low_scout["summary"]["formula_rows"]),
                "low_scout_bernstein_nonnegative_rows": str(low_scout["summary"]["bernstein_nonnegative_rows"]),
                "kernel_identity_found": str(low_scout["summary"]["kernel_identity_found"]).lower(),
            },
            "The selected hard formulas are all certified, but this selected family misses larger contiguous failures.",
            "Use the frontier scout before treating selected low-degree positivity as structural evidence.",
            failures,
        ),
        row(
            "rc_05_degree3_discriminant",
            "degree-3 Jensen discriminant condition",
            "tautological_window_condition",
            "degree-3 Jensen window has nonnegative discriminant",
            False,
            {"countermodel_d3_cubic_discriminant": algebra["finite_countermodel"]["d3_cubic_discriminant"]},
            "This condition would rule out the degree-3 countermodel, but it is the degree-3 hyperbolicity target itself, not a bridge from sign-regularity.",
            "Use only as a finite window consequence after a noncircular theorem proves it.",
            failures,
        ),
        row(
            "rc_06_degree4_discriminant",
            "degree-4 Jensen discriminant condition",
            "tautological_window_condition",
            "degree-4 Jensen window has nonnegative discriminant and real-root certificate",
            False,
            {"countermodel_d4_quartic_discriminant": algebra["finite_countermodel"]["d4_quartic_discriminant"]},
            "This condition would rule out the degree-4 countermodel, but it is a window hyperbolicity condition, not a structural bridge.",
            "Use only as a finite window consequence after a noncircular theorem proves it.",
            failures,
        ),
        row(
            "rc_07_contraction_log_concavity",
            "log-concavity of the contraction sequence",
            "rejected_by_constructed_extension",
            "x2^2 >= x1*x3",
            False,
            {"x2^2 - x1*x3": values["x2^2 - x1*x3"]},
            "A positive rational extension satisfies this stronger contraction condition while preserving negative frontier minors, so the condition is not a standalone bridge.",
            "Do not use contraction log-concavity alone as a bridge lemma; it can only appear inside a stronger noncircular theorem.",
            failures,
        )
        | {
            "constructed_extension_satisfies": True,
            "condition_values_at_constructed_extension": {
                "rho": contraction_extension["ratio_point"]["rho"],
                "x1": contraction_extension["ratio_point"]["x1"],
                "x2": contraction_extension["ratio_point"]["x2"],
                "x3": contraction_extension["ratio_point"]["x3"],
                "x2^2 - x1*x3": contraction_extension["condition_values"]["x2^2 - x1*x3"],
            },
            "target_failures_at_constructed_extension": contraction_failures,
        },
    ]


def build_payload(algebra: dict, low_scout: dict, frontier: dict, contraction_scout: dict) -> dict:
    ratios = countermodel_ratios(algebra)
    rows = build_rows(algebra, low_scout, frontier, contraction_scout, ratios)
    rejected = sum(row["status"] == "rejected_by_countermodel" for row in rows)
    rejected_by_construction = sum(row["status"] == "rejected_by_constructed_extension" for row in rows)
    tautological = sum(row["status"] == "tautological_window_condition" for row in rows)
    open_rows = sum(row["status"] == "open_candidate_not_validated" for row in rows)
    return {
        "kind": "jensen_window_pf_ratio_condition_scout",
        "date": "2026-07-06",
        "target_ansatz": "ansatz_02_positive_cauchy_binet_kernel",
        "source_algebra": "work/rh_compute/results/jensen_window_pf_obligation_algebra.json",
        "source_low_degree_scout": "work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json",
        "source_frontier_scout": "work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json",
        "source_contraction_log_concavity_scout": "work/rh_compute/results/jensen_window_pf_contraction_log_concavity_scout.json",
        "proof_boundary": (
            "Ratio-condition theorem-search diagnostic only; not a proof of a "
            "positive kernel, Cauchy-Binet identity, Jensen-window PF-infinity, "
            "RH, or Lambda <= 0."
        ),
        "countermodel_ratio_point": ratios,
        "summary": {
            "candidate_rows": len(rows),
            "rejected_by_countermodel_rows": rejected,
            "rejected_by_constructed_extension_rows": rejected_by_construction,
            "tautological_window_condition_rows": tautological,
            "open_candidate_not_validated_rows": open_rows,
            "target_closing_rows": 0,
            "main_finding": (
                "Adjacent log-concavity, decreasing contraction variables, "
                "second-order log-concavity, and selected low-degree "
                "Bernstein positivity are all insufficient bridge conditions; "
                "the rational countermodel satisfies them while the frontier "
                "minors d3 m=8 and d4 m=6 are negative.  The remaining "
                "contraction log-concavity candidate is rejected by a positive "
                "rational extension with the same frontier failures."
            ),
        },
        "candidate_rows": rows,
        "invariants": [
            "Rejected rows are satisfied by the countermodel and still have negative target failures.",
            "Rows rejected by constructed extension must include positive condition values and negative target failures at that extension.",
            "Tautological window rows are not allowed to close jwpf_06 as bridge theorems.",
            "Open candidate rows are not promoted without a source theorem and frontier-minor proof.",
            "No row may set target_closing=true.",
        ],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--algebra", type=Path, default=DEFAULT_ALGEBRA)
    parser.add_argument("--low-scout", type=Path, default=DEFAULT_LOW_SCOUT)
    parser.add_argument("--frontier", type=Path, default=DEFAULT_FRONTIER)
    parser.add_argument("--contraction-scout", type=Path, default=DEFAULT_CONTRACTION_SCOUT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload(
        load_json(args.algebra),
        load_json(args.low_scout),
        load_json(args.frontier),
        load_json(args.contraction_scout),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        summary = payload["summary"]
        print(
            "wrote Jensen-window PF ratio-condition scout: "
            f"{summary['candidate_rows']} candidate rows, "
            f"{summary['rejected_by_countermodel_rows']} rejected by countermodel, "
            f"{summary['rejected_by_constructed_extension_rows']} rejected by construction, "
            f"{summary['open_candidate_not_validated_rows']} open"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
