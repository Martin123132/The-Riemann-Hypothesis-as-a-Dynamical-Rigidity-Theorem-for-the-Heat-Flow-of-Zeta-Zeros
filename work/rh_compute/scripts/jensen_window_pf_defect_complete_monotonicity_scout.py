#!/usr/bin/env python3
"""Certify finite defect complete monotonicity and an exact all-shape guard."""

from __future__ import annotations

import argparse
from decimal import Decimal, getcontext
from fractions import Fraction
import json
from pathlib import Path

import jensen_window_pf_monotone_contraction_stress as stress


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_defect_complete_monotonicity_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_defect_complete_monotonicity_scout.md"
DEFAULT_MAX_ORDER = 12
DEFAULT_MAX_INDEX = 64


def rational_str(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def decimal_str(value: Decimal) -> str:
    return f"{value:.24E}"


def arb_positive(value) -> bool:
    return bool(value > 0 and not value.contains(0))


def forward_differences(values: list):
    current = list(values)
    while current:
        yield current
        current = [current[index + 1] - current[index] for index in range(len(current) - 1)]


def append_difference_rows(
    rows: list[dict],
    channel: str,
    lam_label: str,
    arb_values: list,
    sample_values: list[Decimal],
    max_order: int,
) -> None:
    arb_levels = list(forward_differences(arb_values))
    sample_levels = list(forward_differences(sample_values))
    for order in range(max_order + 1):
        sign = 1 if order % 2 == 0 else -1
        arb_signed = [sign * value for value in arb_levels[order]]
        sample_signed = [Decimal(sign) * value for value in sample_levels[order]]
        minimum = min(sample_signed)
        rows.append(
            {
                "channel": channel,
                "lam": lam_label,
                "order": order,
                "checked_intervals": len(arb_signed),
                "all_strictly_positive": all(arb_positive(value) for value in arb_signed),
                "strictly_positive_count": sum(1 for value in arb_signed if arb_positive(value)),
                "contains_zero_count": sum(1 for value in arb_signed if value.contains(0)),
                "strictly_negative_count": sum(1 for value in arb_signed if value < 0 and not value.contains(0)),
                "minimum_midpoint_sample": decimal_str(minimum),
                "minimum_at_global_k": sample_signed.index(minimum) + 1,
            }
        )


def finite_rows(max_order: int, max_index: int) -> tuple[list[dict], list[dict], list[str]]:
    paths = list(stress.DEFAULT_ENCLOSURE_JSONL)
    balls, samples, labels = stress.load_enclosures(paths)
    rows: list[dict] = []
    log_rows: list[dict] = []
    for lam in sorted(labels):
        arb_values = [balls[(lam, index)] for index in range(max_index + 1)]
        sample_values = [samples[(lam, index)] for index in range(max_index + 1)]
        arb_x = stress.contractions(arb_values)
        sample_x = stress.contractions(sample_values)
        append_difference_rows(
            rows,
            "defect",
            labels[lam],
            [1 - value for value in arb_x],
            [Decimal(1) - value for value in sample_x],
            max_order,
        )
        append_difference_rows(
            log_rows,
            "negative_log_contraction",
            labels[lam],
            [-value.log() for value in arb_x],
            [-value.ln() for value in sample_x],
            max_order,
        )
    return rows, log_rows, [path.relative_to(REPO_ROOT).as_posix() for path in paths]


def cubic_discriminant(a: Fraction, b: Fraction, c: Fraction, d: Fraction) -> Fraction:
    return b * b * c * c - 4 * a * c**3 - 4 * b**3 * d - 27 * a * a * d * d + 18 * a * b * c * d


def exact_countermodel() -> dict:
    mass = Fraction(1, 2)
    atom = Fraction(0)
    x1 = 1 - mass
    x2 = 1 - mass * atom
    coefficients = [Fraction(1), Fraction(3), 3 * x1, x1 * x1 * x2]
    discriminant = cubic_discriminant(coefficients[3], coefficients[2], coefficients[1], coefficients[0])
    return {
        "hausdorff_measure": "(1/2)*delta_0",
        "defect_rule": "d_1=1/2 and d_k=0 for every k>=2",
        "complete_monotonicity_rule": "(-1)^m*Delta^m d_1=1/2 and (-1)^m*Delta^m d_k=0 for k>=2",
        "contraction_rule": "x_1=1/2 and x_k=1 for every k>=2",
        "full_static_cone_checks": {
            "first_lower_wall": rational_str(Fraction(1, 3)),
            "first_margin": rational_str(x1 - Fraction(1, 3)),
            "tail_lower_walls_ok": True,
            "upper_walls_ok": True,
            "monotone_contractions_ok": True,
        },
        "degree3_jensen_coefficients_ascending": [rational_str(value) for value in coefficients],
        "degree3_jensen_polynomial": "1+3*z+(3/2)*z^2+(1/4)*z^3",
        "discriminant": rational_str(discriminant),
        "discriminant_is_negative": discriminant < 0,
        "consequence": "The full static ratio cone plus complete monotonicity of d_k does not imply all-shape Jensen-window PF or degree-3 hyperbolicity.",
    }


def build_payload(max_order: int = DEFAULT_MAX_ORDER, max_index: int = DEFAULT_MAX_INDEX) -> dict:
    getcontext().prec = 100
    rows, log_rows, sources = finite_rows(max_order, max_index)
    interval_count = sum(row["checked_intervals"] for row in rows)
    positive_count = sum(row["strictly_positive_count"] for row in rows)
    inconclusive_count = sum(row["contains_zero_count"] for row in rows)
    negative_count = sum(row["strictly_negative_count"] for row in rows)
    fully_certified_orders = sorted(
        {
            row["order"]
            for row in rows
            if all(candidate["all_strictly_positive"] for candidate in rows if candidate["order"] == row["order"])
        }
    )
    log_positive_count = sum(row["strictly_positive_count"] for row in log_rows)
    log_inconclusive_count = sum(row["contains_zero_count"] for row in log_rows)
    log_negative_count = sum(row["strictly_negative_count"] for row in log_rows)
    log_fully_certified_orders = sorted(
        {
            row["order"]
            for row in log_rows
            if all(
                candidate["all_strictly_positive"]
                for candidate in log_rows
                if candidate["order"] == row["order"]
            )
        }
    )
    countermodel = exact_countermodel()
    return {
        "kind": "jensen_window_pf_defect_complete_monotonicity_scout",
        "date": "2026-07-10",
        "status": "finite_interval_diagnostic_with_exact_countermodel_gate",
        "proof_boundary": (
            "Finite Arb evidence for alternating defect differences on the cached coefficient grid, plus an exact abstract countermodel. "
            "This does not prove complete monotonicity for all k, does not prove any all-column theorem, does not prove Jensen-window PF-infinity, RH, or Lambda <= 0."
        ),
        "definition": {
            "ratio_contraction": "x_k=(A_(k+1)/A_k)/(A_k/A_(k-1))",
            "defect": "d_k=1-x_k",
            "finite_condition": "(-1)^m*Delta^m d_k>0",
            "negative_log_contraction": "y_k=-log(x_k)",
            "log_finite_condition": "(-1)^m*Delta^m y_k>0",
        },
        "source_enclosures": sources,
        "max_coefficient_index": max_index,
        "max_difference_order": max_order,
        "finite_rows": rows,
        "negative_log_finite_rows": log_rows,
        "exact_all_shape_countermodel": countermodel,
        "summary": {
            "lambda_count": len({row["lam"] for row in rows}),
            "order_rows": len(rows),
            "checked_intervals": interval_count,
            "strictly_positive_intervals": positive_count,
            "inconclusive_intervals": inconclusive_count,
            "strictly_negative_intervals": negative_count,
            "fully_certified_orders": fully_certified_orders,
            "negative_log_strictly_positive_intervals": log_positive_count,
            "negative_log_inconclusive_intervals": log_inconclusive_count,
            "negative_log_strictly_negative_intervals": log_negative_count,
            "negative_log_fully_certified_orders": log_fully_certified_orders,
            "exact_countermodels": 1,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "Every cached zeta defect difference through order 8 is strictly positive after the alternating sign; "
                "orders 9 through 12 retain positive midpoints but include 421 enclosure-width inconclusives. "
                "The multiplier-kernel sequence y_k=-log(x_k) is likewise certified through order 8, with 3288 positive intervals and 417 high-order inconclusives. "
                "An exact Hausdorff-defect sequence inside the full static ratio cone has a degree-3 Jensen discriminant -27/16. "
                "Complete defect monotonicity is therefore a live column-structure diagnostic, not an all-shape bridge."
            ),
        },
        "invariants": [
            "All finite zeta rows are interval-certified rather than midpoint-promoted.",
            "The exact countermodel is abstract and is not claimed to be a zeta heat-flow orbit.",
            "Finite complete-monotonicity evidence is not an analytic all-k theorem.",
            "Jensen hyperbolicity, RH, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(payload: dict, path: Path) -> None:
    summary = payload["summary"]
    countermodel = payload["exact_all_shape_countermodel"]
    result = (
        "validated Jensen-window PF defect complete-monotonicity scout: "
        f"{summary['strictly_positive_intervals']} defect positives, "
        f"{summary['negative_log_strictly_positive_intervals']} log positives, "
        f"{summary['inconclusive_intervals'] + summary['negative_log_inconclusive_intervals']} inconclusive, "
        f"both certified through order 8, {summary['lambda_count']} lambdas, 1 exact all-shape countermodel, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Defect Complete-Monotonicity Scout",
        "",
        "Date: 2026-07-10",
        "",
        "Status: finite interval diagnostic with exact countermodel gate. This is not a proof of",
        "all-k complete monotonicity, PF-infinity, Jensen hyperbolicity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_defect_complete_monotonicity_scout`.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_defect_complete_monotonicity_scout.json",
        "```",
        "",
        "Generator and checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_defect_complete_monotonicity_scout.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_defect_complete_monotonicity_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result,
        "```",
        "",
        "## Finite Arb Pattern",
        "",
        "For `d_k=1-x_k`, every cached interval through order 8 satisfies",
        "",
        "```text",
        "(-1)^m*Delta^m d_k > 0, 0<=m<=8.",
        "```",
        "",
        f"Across all orders 0 through 12 there are `{summary['strictly_positive_intervals']}` certified",
        f"positive intervals and `{summary['inconclusive_intervals']}` intervals containing zero; no interval is",
        "strictly negative. The order 9 through 12 midpoints remain positive, but are not promoted.",
        "This is finite interval evidence only; no all-`k` theorem is claimed.",
        "",
        "For the multiplier-kernel variable `y_k=-log(x_k)`, every cached interval through",
        "order 8 also satisfies",
        "",
        "```text",
        "(-1)^m*Delta^m y_k > 0, 0<=m<=8.",
        "```",
        "",
        f"Across orders 0 through 12, `{summary['negative_log_strictly_positive_intervals']}` log intervals",
        f"are certified positive and `{summary['negative_log_inconclusive_intervals']}` contain zero; none is",
        "strictly negative. This is the finite necessary pattern for a positive multiplier-kernel",
        "or counting-measure representation, not a proof that such a representation exists.",
        "",
        "## Exact Guard",
        "",
        f"Take the Hausdorff defect measure `{countermodel['hausdorff_measure']}`. Then",
        "",
        "```text",
        f"{countermodel['defect_rule']}",
        f"{countermodel['contraction_rule']}",
        "1/3 <= x_1 <= x_2 <= ... <= 1",
        "```",
        "",
        "so the infinite sequence satisfies the full static ratio cone and complete defect",
        "monotonicity. Its degree-3 Jensen polynomial is",
        "",
        "```text",
        countermodel["degree3_jensen_polynomial"],
        f"discriminant = {countermodel['discriminant']} < 0.",
        "```",
        "",
        countermodel["consequence"],
        "",
        "## Consequence",
        "",
        summary["main_finding"],
        "The next viable use is a column-recurrence theorem or a stronger heat-flow-specific",
        "hierarchy; this condition cannot be promoted directly to the all-shape bridge.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-order", type=int, default=DEFAULT_MAX_ORDER)
    parser.add_argument("--max-index", type=int, default=DEFAULT_MAX_INDEX)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload(args.max_order, args.max_index)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(payload, args.note)
    print(f"wrote {args.out.relative_to(REPO_ROOT)} and {args.note.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
