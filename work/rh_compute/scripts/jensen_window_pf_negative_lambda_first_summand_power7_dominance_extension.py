#!/usr/bin/env python3
"""Strengthen lambda=-100 first-summand dominance to an inverse seventh power."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_negative_lambda_first_summand_dominance_certificate as base  # noqa: E402
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_upper_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.md"
)
BASE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_dominance_certificate.json"
)
START_K = 316
POWER = 7
PRECISION_BITS = 256


@dataclass(frozen=True)
class DominanceRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def arb_text(value: flint.arb) -> str:
    return value.str(70).replace("e", "E")


def source_contract() -> dict:
    source = json.loads(BASE_SOURCE.read_text(encoding="utf-8"))
    if source.get("summary", {}).get("full_tail_power") != 6:
        raise RuntimeError("base first-summand dominance power changed")
    if source.get("summary", {}).get("positive_analytic_gates") != 15:
        raise RuntimeError("base first-summand analytic gate count changed")
    bound = source.get("diagnostics", {}).get("full_tail_relative_bound")
    if bound != (
        "0<=delta_k=(M_k-M_k^(1))/M_k^(1)<=2/k^6 "
        "for every integer k>=300"
    ):
        raise RuntimeError("base first-summand dominance contract changed")
    return {
        "status": source.get("status"),
        "exact_ratio": source["rows"][0]["formula"],
        "ratio_monotonicity": source["rows"][1]["formula"],
        "epsilon_zero_cap": source["diagnostics"]["epsilon_zero"]["cap"],
    }


def endpoint_diagnostics() -> dict:
    flint.ctx.prec = PRECISION_BITS
    k = flint.arb(START_K)
    log_k = k.log()
    sqrt_k = k.sqrt()
    pi = flint.arb.pi()
    epsilon_cap = flint.arb(base.DEFAULT_EPSILON_ZERO_CAP)

    high_endpoint = flint.arb(17).log() - 3 * pi * sqrt_k + POWER * log_k
    high_derivative_margin = 3 * pi * sqrt_k / 2 - POWER

    ga = 16 * k / log_k - 25 * log_k + 5 - 4 * pi * sqrt_k
    ga_sqrt_margin = 8 * (log_k - 1) / log_k**2 - 2 * pi / sqrt_k
    ga_k_margin = 8 * (log_k - 1) / log_k**2 - 25 / k

    offset_c = flint.arb(22) / 25
    exp_c = (flint.arb(11) / 25).exp()
    gc = (
        16 * k / (log_k + offset_c)
        - 25 * (log_k + offset_c)
        + 5
        - 4 * pi * exp_c * sqrt_k
    )
    gc_sqrt_margin = (
        8 * (log_k - flint.arb(3) / 25) / (log_k + offset_c) ** 2
        - 2 * pi * exp_c / sqrt_k
    )
    gc_k_margin = (
        8 * (log_k - flint.arb(3) / 25) / (log_k + offset_c) ** 2
        - 25 / k
    )

    exp_h = (flint.arb(2) / 5).exp()
    coarse_low = (
        -8 * k / (5 * (log_k + flint.arb(4) / 5))
        + 5 * log_k / 2
        + 1
        + pi * (exp_h - 1) * sqrt_k
    )
    low_endpoint = coarse_low + epsilon_cap.log() + POWER * log_k
    half_negative = (
        4
        * (log_k - flint.arb(1) / 5)
        / (5 * (log_k + flint.arb(4) / 5) ** 2)
    )
    low_sqrt_margin = half_negative - pi * (exp_h - 1) / (2 * sqrt_k)
    low_k_margin = half_negative - flint.arb(2 * POWER + 5) / (2 * k)

    gates = {
        "high_endpoint_below_k_minus_7": bool(high_endpoint < 0),
        "high_comparison_decreasing": bool(high_derivative_margin > 0),
        "D_a_lower_above_100": bool(ga > 100),
        "D_a_lower_increasing_sqrt_side": bool(ga_sqrt_margin > 0),
        "D_a_lower_increasing_k_side": bool(ga_k_margin > 0),
        "D_c_lower_positive": bool(gc > 0),
        "D_c_lower_increasing_sqrt_side": bool(gc_sqrt_margin > 0),
        "D_c_lower_increasing_k_side": bool(gc_k_margin > 0),
        "low_endpoint_below_k_minus_7": bool(low_endpoint < 0),
        "low_comparison_decreasing_sqrt_side": bool(low_sqrt_margin > 0),
        "low_comparison_decreasing_k_side": bool(low_k_margin > 0),
    }
    if not all(gates.values()):
        failed = [name for name, passed in gates.items() if not passed]
        raise RuntimeError(f"power-seven endpoint gates failed: {failed}")

    return {
        "parameters": {
            "lambda": "-100",
            "start_k": START_K,
            "power": POWER,
            "precision_bits": PRECISION_BITS,
            "epsilon_zero_cap": base.DEFAULT_EPSILON_ZERO_CAP,
        },
        "high_endpoint_ball": arb_text(high_endpoint),
        "high_endpoint_upper": arb_upper_text(high_endpoint),
        "high_derivative_margin_lower": arb_lower_text(high_derivative_margin),
        "low_endpoint_ball": arb_text(low_endpoint),
        "low_endpoint_upper": arb_upper_text(low_endpoint),
        "low_sqrt_derivative_margin_lower": arb_lower_text(low_sqrt_margin),
        "low_k_derivative_margin_lower": arb_lower_text(low_k_margin),
        "ga_minus_100_lower": arb_lower_text(ga - 100),
        "gc_lower": arb_lower_text(gc),
        "positive_gates": gates,
        "all_positive_gates": True,
        "half_line_rule": (
            "the same positive ratio-log derivatives as the power-six proof "
            "show every split endpoint comparison strengthens for k>=316"
        ),
        "full_tail_relative_bound": (
            "0<=delta_k=(M_k-M_k^(1))/M_k^(1)<=2/k^7 "
            "for every integer k>=316"
        ),
        "moment_log_error": "0<=log(1+delta_k)<=2/k^7, k>=316",
        "adjacent_B_error": (
            "|B_j-B_j^(1)|<=a_j=2*((j-1)^(-7)+2*j^(-7)+(j+1)^(-7)), "
            "j>=317"
        ),
    }


def build_artifact() -> dict:
    source = source_contract()
    diagnostics = endpoint_diagnostics()
    rows = [
        DominanceRow(
            "fsd7_01_ratio_geometry",
            "theorem_input",
            "ready_to_apply",
            "The exact positive summand ratios and their pointwise monotonicity are unchanged.",
            source["exact_ratio"] + "; " + source["ratio_monotonicity"],
            "Inherited exact kernel geometry.",
        ),
        DominanceRow(
            "fsd7_02_high_region",
            "interval_analytic_certificate",
            "ready_to_apply",
            "The high-u tail is below k^-7 at the endpoint and the comparison improves on the half-line.",
            "17*exp(-3*pi*sqrt(k))<k^(-7), k>=316",
            "One Arb endpoint gate plus exact derivative monotonicity.",
            {
                "endpoint_upper": diagnostics["high_endpoint_upper"],
                "derivative_margin_lower": diagnostics["high_derivative_margin_lower"],
            },
        ),
        DominanceRow(
            "fsd7_03_low_region",
            "interval_analytic_certificate",
            "ready_to_apply",
            "The low-u probability contribution is below k^-7 at the endpoint and decreases thereafter.",
            "epsilon(0)*P_k(u<a(k))<k^(-7), k>=316",
            "One Arb endpoint gate and the inherited split derivative argument.",
            {
                "endpoint_upper": diagnostics["low_endpoint_upper"],
                "sqrt_margin_lower": diagnostics["low_sqrt_derivative_margin_lower"],
                "k_margin_lower": diagnostics["low_k_derivative_margin_lower"],
            },
        ),
        DominanceRow(
            "fsd7_04_full_tail",
            "analytic_theorem",
            "ready_to_apply",
            "The two regions compose to strengthen the full moment-tail dominance by one inverse power.",
            diagnostics["full_tail_relative_bound"],
            "Lambda=-100 integer moments only.",
        ),
        DominanceRow(
            "fsd7_05_log_error",
            "exact_consequence",
            "ready_to_apply",
            "The relative tail bound gives a seventh-power logarithmic perturbation bound.",
            diagnostics["moment_log_error"],
            "Uses log(1+z)<=z for z>=0.",
        ),
        DominanceRow(
            "fsd7_06_B_error",
            "exact_consequence",
            "ready_to_apply",
            "The adjacent log wall inherits an explicit seventh-power stencil error.",
            diagnostics["adjacent_B_error"],
            "Exact four-term triangle inequality.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension",
        "date": "2026-07-13",
        "status": "all-k inverse-seventh-power first-summand dominance extension",
        "proof_boundary": (
            "This artifact strengthens first-summand dominance only for lambda=-100 "
            "and k>=316. It does not prove an order-six curvature inequality, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": source,
        "diagnostics": diagnostics,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": len(rows),
            "positive_analytic_gates": len(diagnostics["positive_gates"]),
            "full_tail_power": POWER,
            "tail_start_k": START_K,
            "dominance_theorems": 1,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    d = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF First-Summand Power-Seven Dominance Extension",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous inverse-seventh-power first-summand dominance at",
        "`lambda=-100`. This is not a proof of order six, PF-infinity, RH, or",
        "`Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.json",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.py",
        "```",
        "",
        "## Strengthened Tail",
        "",
        "The exact summand-ratio monotonicity and tilted-integrand geometry from",
        "the power-six proof are unchanged. Re-running only the endpoint power",
        "comparison at `k=316` gives",
        "",
        "```text",
        d["full_tail_relative_bound"],
        d["moment_log_error"],
        "```",
        "",
        "The two decisive outward-rounded endpoint balls are",
        "",
        "```text",
        "high log comparison=" + d["high_endpoint_ball"] + "<0,",
        "low log comparison=" + d["low_endpoint_ball"] + "<0.",
        "```",
        "",
        "All eleven derivative and endpoint gates are strictly positive. The",
        "same ratio-log derivative checks used in the original certificate show",
        "that both comparisons strengthen for every `k>=316`.",
        "",
        "## Wall Error",
        "",
        "For the adjacent logarithmic coordinate this yields",
        "",
        "```text",
        d["adjacent_B_error"],
        "```",
        "",
        "That extra inverse power is what keeps the third nested stable-log",
        "transfer below a fixed multiple of `k^-2`.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
        "outputs/jensen_window_pf_compound_order6_m100_tail_curvature_reduction.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote power-seven first-summand dominance: "
        f"{summary['rows']} rows, {summary['positive_analytic_gates']} gates, "
        f"tail k>={summary['tail_start_k']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
