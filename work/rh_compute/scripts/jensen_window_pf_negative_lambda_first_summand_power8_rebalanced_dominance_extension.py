#!/usr/bin/env python3
"""Prove inverse-eighth-power first-summand dominance with a rebalanced split."""

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
    "jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.md"
)
BASE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_dominance_certificate.json"
)
START_K = 300
POWER = 8
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


def endpoint_diagnostics(
    power: int = POWER,
    start_k: int = START_K,
) -> dict:
    flint.ctx.prec = PRECISION_BITS
    k = flint.arb(start_k)
    log_k = k.log()
    pi = flint.arb.pi()
    exponent = flint.arb(2) / 5
    k_power = k**exponent
    epsilon_cap = flint.arb(base.DEFAULT_EPSILON_ZERO_CAP)

    high_geometric = (flint.arb(3) / 2) ** 4 * (-5 * pi * k_power).exp()
    high_prefactor = (
        16
        * (2 * pi * k_power / (2 * pi * k_power - 3))
        / (1 - high_geometric)
    )
    high_endpoint = flint.arb(17).log() - 3 * pi * k_power + power * log_k
    high_derivative_margin = 6 * pi * k_power / 5 - power

    ga = 20 * k / log_k - 20 * log_k + 5 - 4 * pi * k_power
    ga_power_margin = (
        10 * (log_k - 1) / log_k**2
        - 8 * pi * k ** (-flint.arb(3) / 5) / 5
    )
    ga_log_margin = 10 * (log_k - 1) / log_k**2 - 20 / k

    offset_c = flint.arb(11) / 10
    exp_c = (flint.arb(11) / 25).exp()
    gc = (
        20 * k / (log_k + offset_c)
        - 20 * (log_k + offset_c)
        + 5
        - 4 * pi * exp_c * k_power
    )
    gc_power_margin = (
        10
        * (log_k + flint.arb(1) / 10)
        / (log_k + offset_c) ** 2
        - 8 * pi * exp_c * k ** (-flint.arb(3) / 5) / 5
    )
    gc_log_margin = (
        10
        * (log_k + flint.arb(1) / 10)
        / (log_k + offset_c) ** 2
        - 20 / k
    )

    exp_h = (flint.arb(2) / 5).exp()
    coarse_low = (
        -2 * k / (log_k + 1)
        + 2 * log_k
        + 1
        + pi * (exp_h - 1) * k_power
    )
    low_endpoint = coarse_low + epsilon_cap.log() + power * log_k
    low_negative = 2 * log_k / (log_k + 1) ** 2
    low_power_margin = (
        low_negative
        - 2
        * pi
        * (exp_h - 1)
        * k ** (-flint.arb(3) / 5)
        / 5
    )
    low_k_margin = low_negative - 10 / k

    a = log_k / 10
    b = a + flint.arb(1) / 10
    c = b + flint.arb(1) / 100
    derivative_a = base.log_integrand_derivative(k, a)
    derivative_c = base.log_integrand_derivative(k, c)

    gates = {
        "high_prefactor_below_17": bool(high_prefactor < 17),
        "high_endpoint_below_k_minus_8": bool(high_endpoint < 0),
        "high_comparison_decreasing": bool(high_derivative_margin > 0),
        "D_a_lower_above_100": bool(ga > 100),
        "D_a_lower_increasing_power_side": bool(ga_power_margin > 0),
        "D_a_lower_increasing_k_side": bool(ga_log_margin > 0),
        "D_c_lower_positive": bool(gc > 0),
        "D_c_lower_increasing_power_side": bool(gc_power_margin > 0),
        "D_c_lower_increasing_k_side": bool(gc_log_margin > 0),
        "low_endpoint_below_k_minus_8": bool(low_endpoint < 0),
        "low_comparison_decreasing_power_side": bool(low_power_margin > 0),
        "low_comparison_decreasing_k_side": bool(low_k_margin > 0),
        "exact_D_a_positive": bool(derivative_a > 0),
        "exact_D_c_positive": bool(derivative_c > 0),
    }
    if not all(gates.values()):
        failed = [name for name, passed in gates.items() if not passed]
        raise RuntimeError(f"power-{power} endpoint gates failed: {failed}")

    return {
        "parameters": {
            "lambda": "-100",
            "start_k": start_k,
            "power": power,
            "precision_bits": PRECISION_BITS,
            "adaptive_a": "a(k)=log(k)/10",
            "adaptive_b": "b(k)=a(k)+1/10",
            "adaptive_c": "c(k)=b(k)+1/100",
            "epsilon_zero_cap": base.DEFAULT_EPSILON_ZERO_CAP,
        },
        "high_prefactor_ball": arb_text(high_prefactor),
        "high_prefactor_margin_lower": arb_lower_text(17 - high_prefactor),
        "high_endpoint_ball": arb_text(high_endpoint),
        "high_endpoint_upper": arb_upper_text(high_endpoint),
        "high_derivative_margin_lower": arb_lower_text(high_derivative_margin),
        "low_endpoint_ball": arb_text(low_endpoint),
        "low_endpoint_upper": arb_upper_text(low_endpoint),
        "low_power_derivative_margin_lower": arb_lower_text(low_power_margin),
        "low_k_derivative_margin_lower": arb_lower_text(low_k_margin),
        "ga_minus_100_lower": arb_lower_text(ga - 100),
        "ga_power_margin_lower": arb_lower_text(ga_power_margin),
        "ga_log_margin_lower": arb_lower_text(ga_log_margin),
        "gc_lower": arb_lower_text(gc),
        "gc_power_margin_lower": arb_lower_text(gc_power_margin),
        "gc_log_margin_lower": arb_lower_text(gc_log_margin),
        "exact_D_a_lower": arb_lower_text(derivative_a),
        "exact_D_c_lower": arb_lower_text(derivative_c),
        "positive_gates": gates,
        "all_positive_gates": True,
        "half_line_rule": (
            "the six positive ratio-log derivatives with leading constants "
            "3/5 or 1 show every split endpoint comparison strengthens for k>=300"
        ),
        "ratio_log_derivatives": [
            "3/5+1/(L-1)-2/L>0",
            "1+1/(L-1)-2/L>0",
            "3/5+1/(L+1/10)-2/(L+11/10)>0",
            "1+1/(L+1/10)-2/(L+11/10)>0",
            "3/5+1/L-2/(L+1)>0",
            "1+1/L-2/(L+1)>0",
        ],
        "high_region_bound": (
            f"epsilon(a(k))<=17*exp(-3*pi*k^(2/5))<k^(-{power}), "
            f"k>={start_k}"
        ),
        "low_region_bound": (
            f"epsilon(0)*P_k(u<a(k))<k^(-{power}), k>={start_k}"
        ),
        "full_tail_relative_bound": (
            f"0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^{power} "
            f"for every integer k>={start_k}"
        ),
        "moment_log_error": (
            f"0<=log(1+delta_k)<2/k^{power}, k>={start_k}"
        ),
        "adjacent_B_error": (
            f"|B_j-B_j^(1)|<=a_j=2*((j-1)^(-{power})+"
            f"2*j^(-{power})+(j+1)^(-{power})), j>={start_k + 1}"
        ),
    }


def build_artifact() -> dict:
    source = source_contract()
    diagnostics = endpoint_diagnostics()
    rows = [
        DominanceRow(
            "fsd8r_01_ratio_geometry",
            "theorem_input",
            "ready_to_apply",
            "The exact positive summand ratios, monotonicity, and epsilon(0) cap are inherited unchanged.",
            source["exact_ratio"] + "; " + source["ratio_monotonicity"],
            "Inherited exact kernel geometry and rigorous epsilon(0) bound.",
        ),
        DominanceRow(
            "fsd8r_02_rebalanced_split",
            "exact_reduction",
            "ready_to_apply",
            "Moving the split inward spends the large high-region reserve to strengthen the low-region power.",
            "a(k)=log(k)/10, b(k)=a(k)+1/10, c(k)=b(k)+1/100",
            "Exact split definition only.",
        ),
        DominanceRow(
            "fsd8r_03_high_region",
            "interval_analytic_certificate",
            "ready_to_apply",
            "The complete higher-summand kernel ratio on the high region is below k^-8 on the full half-line.",
            diagnostics["high_region_bound"],
            "Arb endpoint gates plus exact derivative monotonicity.",
        ),
        DominanceRow(
            "fsd8r_04_saddle_geometry",
            "interval_analytic_certificate",
            "ready_to_apply",
            "The log-integrand derivative stays positive across the comparison collar for every k>=300.",
            "S_k'(a(k))>100 and S_k'(c(k))>0",
            "Arb endpoint gates plus six exact ratio-monotonicity checks.",
        ),
        DominanceRow(
            "fsd8r_05_low_region",
            "interval_analytic_certificate",
            "ready_to_apply",
            "The tilted first-summand probability below the new split is below k^-8.",
            diagnostics["low_region_bound"],
            "Arb endpoint gate and exact half-line derivative comparisons.",
        ),
        DominanceRow(
            "fsd8r_06_full_tail",
            "analytic_theorem",
            "ready_to_apply",
            "The two regions compose to an inverse-eighth-power complete-to-first moment defect.",
            diagnostics["full_tail_relative_bound"],
            "Lambda=-100 integer moments only.",
        ),
        DominanceRow(
            "fsd8r_07_wall_error",
            "exact_consequence",
            "ready_to_apply",
            "The adjacent logarithmic wall inherits the inverse-eighth-power stencil bound.",
            diagnostics["moment_log_error"] + "; " + diagnostics["adjacent_B_error"],
            "Uses log(1+z)<=z and the exact centered stencil.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension",
        "date": "2026-07-13",
        "status": "rebalanced all-k inverse-eighth-power first-summand dominance extension",
        "proof_boundary": (
            "This artifact strengthens first-summand dominance only for lambda=-100 "
            "and k>=300. It does not prove the order-seven curvature ceiling, "
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
            "rebalanced_splits": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    d = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF Rebalanced Power-Eight First-Summand Dominance",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous inverse-eighth-power first-summand dominance at",
        "`lambda=-100`. This is not a proof of order seven, PF-infinity, RH,",
        "or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.json",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.py",
        "```",
        "",
        "## Rebalanced Split",
        "",
        "The exact summand-ratio monotonicity and the rigorous",
        "`epsilon(0)<0.0022` theorem are unchanged. Replace the former split by",
        "",
        "```text",
        "a(k)=log(k)/10, b(k)=a(k)+1/10, c(k)=b(k)+1/100.",
        "q(a(k))=pi*k^(2/5).",
        "```",
        "",
        "The high-region exponential reserve remains large while the tilted",
        "low-region probability gains the extra inverse power.",
        "",
        "## Endpoint Gates",
        "",
        "At `k=300`, outward-rounded Arb arithmetic gives",
        "",
        "```text",
        "high log comparison=" + d["high_endpoint_ball"] + "<0,",
        "low log comparison=" + d["low_endpoint_ball"] + "<0,",
        "high prefactor reserve=" + d["high_prefactor_margin_lower"] + ">0,",
        "minimum low power-derivative margin=" + d["low_power_derivative_margin_lower"] + ">0.",
        "```",
        "",
        "All fourteen endpoint and derivative gates are strictly positive.",
        "The six displayed ratio-log derivatives prove that every comparison",
        "strengthens for the complete half-line `k>=300`.",
        "",
        "## Strengthened Tail",
        "",
        "```text",
        d["high_region_bound"],
        d["low_region_bound"],
        d["full_tail_relative_bound"],
        d["moment_log_error"],
        d["adjacent_B_error"],
        "```",
        "",
        "The inverse-eighth-power wall error is the input needed by the fourth",
        "stable logarithm in the order-seven complete-to-first transfer.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
        "outputs/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.md",
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
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote power-eight rebalanced first-summand dominance: "
        f"{summary['rows']} rows, {summary['positive_analytic_gates']} gates, "
        f"tail k>={summary['tail_start_k']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
