#!/usr/bin/env python3
"""Transfer order-nine sixth-nested curvature from first summand to full kernel."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order5_first_summand_curvature_bridge import (  # noqa: E402
    shifted_positive_polynomial,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_first_summand_curvature_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order9_first_summand_curvature_bridge.md"
)
TAIL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_m100_tail_curvature_reduction.json"
)
POWER10_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_power10_rebalanced_dominance_extension.json"
)
SHARP_ORDER8_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate.json"
)
ORDER4_BRIDGE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_first_summand_curvature_bridge.json"
)
ORDER4_ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_m100_entry_certificate.json"
)
TAIL_FIRST_K = 1251
GAP_FLOOR_START = TAIL_FIRST_K - 1
CONTINUOUS_FIRST_T = 1250
FIRST_CONTINUOUS_CONSTANT = 4200
FULL_TRANSFER_CONSTANT = 550
DISCRETE_FIRST_CONSTANT = 4201
FULL_CURVATURE_CONSTANT = 4900


@dataclass(frozen=True)
class BridgeRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_sources() -> dict:
    tail = load_json(TAIL_SOURCE)
    power10 = load_json(POWER10_SOURCE)
    sharp8 = load_json(SHARP_ORDER8_SOURCE)
    order4_bridge = load_json(ORDER4_BRIDGE_SOURCE)
    order4_entry = load_json(ORDER4_ENTRY_SOURCE)
    if tail.get("exact", {}).get("sufficient_ceiling") != (
        "Y_k<=4900/k^2 for every real/integer k>=1249"
    ):
        raise RuntimeError("order-nine tail reduction contract changed")
    if power10.get("summary", {}).get("full_tail_power") != 10:
        raise RuntimeError("power-ten dominance source changed")
    if power10.get("summary", {}).get("tail_start_k") != 300:
        raise RuntimeError("power-ten dominance start changed")
    if sharp8.get("exact", {}).get("global_first_curvature") != (
        "s_1''(t)<=2500/t^2 for every real t>=1249"
    ):
        raise RuntimeError("sharp order-eight first curvature changed")
    if sharp8.get("exact", {}).get("sharp_full_ceiling") != (
        "W_k<2501/k^2+190/k^2=2691/k^2, k>=1250"
    ):
        raise RuntimeError("sharp order-eight full ceiling changed")
    if order4_bridge.get("exact", {}).get("integer_B_bounds") != (
        "1/(2*m+1)<=B_1(m)<=3/(2*m-1), integer m>=319"
    ):
        raise RuntimeError("first-summand B floor changed")
    if order4_entry.get("tail_arithmetic", {}).get("defect_buffer") != (
        "-3*log(x_k)>3*d_k>=753/(250*(2*k+1))"
    ):
        raise RuntimeError("full-kernel defect floor changed")
    return {
        "tail_ceiling": tail["exact"]["sufficient_ceiling"],
        "power10_bound": power10["diagnostics"]["full_tail_relative_bound"],
        "order8_first_curvature": sharp8["exact"]["global_first_curvature"],
        "order8_full_ceiling": sharp8["exact"]["sharp_full_ceiling"],
        "first_B_floor": order4_bridge["exact"]["integer_B_bounds"],
        "full_B_floor": "B(j)>=d_j>=251/(250*(2*j+1)), j>=320",
    }


def rationalized_error_chain() -> dict[str, sp.Expr]:
    j = sp.symbols("j", integer=True, positive=True)

    def shift(expression: sp.Expr, offset: int) -> sp.Expr:
        return expression.subs(j, j + offset)

    def reduce(expression: sp.Expr) -> sp.Expr:
        return sp.cancel(expression)

    a = reduce(2 * ((j - 1) ** -10 + 2 * j**-10 + (j + 1) ** -10))
    ell = reduce(4 * j * a)
    first_gap = reduce(2 * a + shift(ell, -1) + 2 * ell + shift(ell, 1))
    log_first_gap = reduce(2 * ell + 8 * j * first_gap)
    second_gap = reduce(
        3 * a
        + shift(log_first_gap, -1)
        + 2 * log_first_gap
        + shift(log_first_gap, 1)
    )
    order5_coordinate = reduce(
        2 * log_first_gap + sp.Rational(5, 7) * j * second_gap + ell
    )
    third_gap = reduce(
        4 * a
        + shift(order5_coordinate, -1)
        + 2 * order5_coordinate
        + shift(order5_coordinate, 1)
    )
    order6_coordinate = reduce(
        2 * order5_coordinate
        + log_first_gap
        + sp.Rational(2, 3) * j * third_gap
    )
    fourth_gap = reduce(
        5 * a
        + shift(order6_coordinate, -1)
        + 2 * order6_coordinate
        + shift(order6_coordinate, 1)
    )
    order7_coordinate = reduce(
        2 * order6_coordinate
        + order5_coordinate
        + sp.Rational(2, 3) * j * fourth_gap
    )
    fifth_gap = reduce(
        6 * a
        + shift(order7_coordinate, -1)
        + 2 * order7_coordinate
        + shift(order7_coordinate, 1)
    )
    order8_coordinate = reduce(
        2 * order7_coordinate
        + order6_coordinate
        + sp.Rational(2, 3) * j * fifth_gap
    )
    sixth_gap = reduce(
        7 * a
        + shift(order8_coordinate, -1)
        + 2 * order8_coordinate
        + shift(order8_coordinate, 1)
    )
    order9_coordinate = reduce(
        2 * order8_coordinate
        + order7_coordinate
        + sp.Rational(3, 4) * j * sixth_gap
    )
    total_curvature = reduce(
        shift(order9_coordinate, -1)
        + 2 * order9_coordinate
        + shift(order9_coordinate, 1)
    )
    return {
        "symbol": j,
        "a": a,
        "ell": ell,
        "first_gap": first_gap,
        "log_first_gap": log_first_gap,
        "second_gap": second_gap,
        "order5_coordinate": order5_coordinate,
        "third_gap": third_gap,
        "order6_coordinate": order6_coordinate,
        "fourth_gap": fourth_gap,
        "order7_coordinate": order7_coordinate,
        "fifth_gap": fifth_gap,
        "order8_coordinate": order8_coordinate,
        "sixth_gap": sixth_gap,
        "order9_coordinate": order9_coordinate,
        "total_curvature": total_curvature,
    }


def exact_diagnostics() -> dict:
    errors = rationalized_error_chain()
    j = errors["symbol"]
    k = sp.symbols("k", integer=True, positive=True)

    first_v_floor = shifted_positive_polynomial(
        7 / (2 * j + 1) - 2501 / j**2 - sp.Rational(4, 3) / j,
        j,
        GAP_FLOOR_START,
    )
    full_v_floor = shifted_positive_polynomial(
        sp.Rational(1757, 250) / (2 * j + 1)
        - 2691 / j**2
        - sp.Rational(4, 3) / j,
        j,
        GAP_FLOOR_START,
    )

    total_error = errors["total_curvature"].subs(j, k)
    transfer_difference = sp.cancel(
        sp.Rational(FULL_TRANSFER_CONSTANT, 1) / k**2 - total_error
    )
    numerator, denominator = sp.fraction(transfer_difference)
    shift = sp.symbols("m", integer=True, nonnegative=True)
    transfer_polynomial = sp.Poly(
        sp.expand(numerator.subs(k, TAIL_FIRST_K + shift)), shift
    )
    transfer_coefficients = transfer_polynomial.all_coeffs()
    if any(value <= 0 for value in transfer_coefficients):
        raise RuntimeError("order-nine transfer polynomial is not positive")
    endpoint_error = sp.cancel(total_error.subs(k, TAIL_FIRST_K))
    endpoint_scaled = sp.cancel(TAIL_FIRST_K**2 * endpoint_error)

    t = sp.symbols("t", positive=True)
    phi = 1 / (sp.exp(t) - 1)
    chi = sp.exp(t) / (sp.exp(t) - 1) ** 2
    first_residual = sp.simplify(sp.diff(sp.log(1 - sp.exp(-t)), t) - phi)
    second_residual = sp.simplify(
        sp.diff(sp.log(1 - sp.exp(-t)), t, 2) + chi
    )
    if first_residual != 0 or second_residual != 0:
        raise RuntimeError("stable logarithm derivative identity failed")

    return {
        "continuous_coordinates": (
            "V(t)=7*B(t)-s(t-1)+2*s(t)-s(t+1); "
            "w(t)=2*s(t)-r(t)+log(1-exp(-V(t)))"
        ),
        "order9_curvature": (
            "w''=2*s''-r''+phi(V)*V''-chi(V)*(V')^2; "
            "phi(z)=1/(exp(z)-1), chi(z)=exp(z)/(exp(z)-1)^2"
        ),
        "discrete_identity": (
            "Y_k=w(k-1)-2*w(k)+w(k+1)="
            "integral_[-1,1](1-|v|)*w''(k+v) dv"
        ),
        "moment_error": (
            "a_j=2*((j-1)^(-10)+2*j^(-10)+(j+1)^(-10)); "
            "|B_j-B_j^(1)|<=a_j"
        ),
        "log_defect_error": "L_j=4*j*a_j; |ell_j-ell_j^(1)|<=L_j",
        "first_gap_error": (
            "U1_j=2*a_j+L_(j-1)+2*L_j+L_(j+1); "
            "|J_j-J_j^(1)|<=U1_j"
        ),
        "log_first_gap_error": (
            "V1_j=2*L_j+8*j*U1_j; |h_j-h_j^(1)|<=V1_j"
        ),
        "second_gap_error": (
            "W1_j=3*a_j+V1_(j-1)+2*V1_j+V1_(j+1); "
            "|R_j-R_j^(1)|<=W1_j"
        ),
        "order5_coordinate_error": (
            "E_j=2*V1_j+(5*j/7)*W1_j+L_j; |q_j-q_j^(1)|<=E_j"
        ),
        "third_gap_error": (
            "Z_j=4*a_j+E_(j-1)+2*E_j+E_(j+1); "
            "|S_j-S_j^(1)|<=Z_j"
        ),
        "order6_coordinate_error": (
            "Y1_j=2*E_j+V1_j+(2*j/3)*Z_j; |p_j-p_j^(1)|<=Y1_j"
        ),
        "fourth_gap_error": (
            "O_j=5*a_j+Y1_(j-1)+2*Y1_j+Y1_(j+1); "
            "|T_j-T_j^(1)|<=O_j"
        ),
        "order7_coordinate_error": (
            "N_j=2*Y1_j+E_j+(2*j/3)*O_j; |r_j-r_j^(1)|<=N_j"
        ),
        "fifth_gap_error": (
            "P_j=6*a_j+N_(j-1)+2*N_j+N_(j+1); "
            "|U_j-U_j^(1)|<=P_j"
        ),
        "order8_coordinate_error": (
            "C_j=2*N_j+Y1_j+(2*j/3)*P_j; |s_j-s_j^(1)|<=C_j"
        ),
        "sixth_gap_error": (
            "D_j=7*a_j+C_(j-1)+2*C_j+C_(j+1); "
            "|V_j-V_j^(1)|<=D_j"
        ),
        "sixth_gap_floor": (
            "min(V_j,V_j^(1))>=4/(3*j), j>=1250"
        ),
        "stable_log_lipschitz": (
            "phi(min(V_j,V_j^(1)))<=1/min(V_j,V_j^(1))<=3*j/4"
        ),
        "order9_coordinate_error": (
            "F_j=2*C_j+N_j+(3*j/4)*D_j; |w_j-w_j^(1)|<=F_j"
        ),
        "full_transfer": (
            "|Y_k-Y_k^(1)|<=F_(k-1)+2*F_k+F_(k+1)<550/k^2, k>=1251"
        ),
        "continuous_target": "w_1''(t)<=4200/t^2 for every real t>=1250",
        "tent_transfer": (
            "w_1''(t)<=4200/t^2 => Y_k^(1)<=4200*[-log(1-1/k^2)]"
            "<4201/k^2, k>=1251"
        ),
        "conditional_full_ceiling": (
            "w_1''(t)<=4200/t^2 on t>=1250 => "
            "Y_k<4201/k^2+550/k^2=4751/k^2<4900/k^2, k>=1251"
        ),
        "conditional_endpoint_tail": (
            "w_1''(t)<=4200/t^2 on t>=1250 => "
            "Q_(9,n)(-100)>0 for every n>=1243"
        ),
        "finite_splice": "prove Q_(9,n)(-100)>0 for n=1241,1242",
        "first_V_floor_polynomial": first_v_floor,
        "full_V_floor_polynomial": full_v_floor,
        "transfer_polynomial": {
            "start": TAIL_FIRST_K,
            "degree": transfer_polynomial.degree(),
            "coefficient_count": len(transfer_coefficients),
            "leading_coefficient": str(transfer_coefficients[0]),
            "constant_coefficient": str(transfer_coefficients[-1]),
            "minimum_coefficient": str(min(transfer_coefficients)),
            "denominator_degree": int(sp.degree(denominator, k)),
            "shifted_coefficients": [str(value) for value in transfer_coefficients],
        },
        "endpoint_transfer_error": str(endpoint_error),
        "endpoint_scaled_transfer": str(endpoint_scaled),
        "endpoint_scaled_transfer_decimal": str(sp.N(endpoint_scaled, 30)),
        "endpoint_scaled_reserve_below_550": str(
            sp.N(FULL_TRANSFER_CONSTANT - endpoint_scaled, 30)
        ),
        "stable_derivative_residuals": [str(first_residual), str(second_residual)],
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_diagnostics()
    rows = [
        BridgeRow(
            "co9fscb_01_coordinate",
            "exact_identity",
            "ready_to_apply",
            "Order nine adds one sixth nested stable logarithm to the completed hierarchy.",
            exact["continuous_coordinates"] + "; " + exact["order9_curvature"],
            "Exact differentiation only.",
        ),
        BridgeRow(
            "co9fscb_02_power10_input",
            "theorem_input",
            "ready_to_apply",
            "The rebalanced first-summand theorem supplies the perturbation power required by the sixth stable logarithm.",
            sources["power10_bound"],
            "Completed lambda=-100 moment-tail theorem.",
        ),
        BridgeRow(
            "co9fscb_03_sixth_floor",
            "exact_analytic_composition",
            "ready_to_apply",
            "The sharpened order-eight ceiling keeps first and full sixth-gap coordinates uniformly away from zero.",
            exact["sixth_gap_floor"],
            "Composes sharp order-eight first/full curvature and defect floors.",
            {
                "first": exact["first_V_floor_polynomial"],
                "full": exact["full_V_floor_polynomial"],
            },
        ),
        BridgeRow(
            "co9fscb_04_error_chain",
            "exact_error_propagation",
            "ready_to_apply",
            "The tenth-power moment error survives all six stable logarithms as an inverse-square curvature error.",
            exact["moment_error"] + "; " + exact["order9_coordinate_error"],
            "Exact triangle and mean-value inequalities inside proved positive floors.",
        ),
        BridgeRow(
            "co9fscb_05_full_transfer",
            "analytic_theorem",
            "ready_to_apply",
            "The complete theta kernel changes the order-nine centered curvature by fewer than 550 inverse squares.",
            exact["full_transfer"],
            "Exact degree-168 coefficient-positive rational audit.",
            {
                "endpoint_scaled": exact["endpoint_scaled_transfer_decimal"],
                "reserve": exact["endpoint_scaled_reserve_below_550"],
            },
        ),
        BridgeRow(
            "co9fscb_06_tent",
            "conditional_theorem",
            "ready_to_apply",
            "A continuous first-summand ceiling contributes fewer than 4201 inverse squares discretely.",
            exact["tent_transfer"],
            "Conditional only on the displayed continuous curvature theorem.",
        ),
        BridgeRow(
            "co9fscb_07_full_ceiling",
            "conditional_theorem",
            "ready_to_apply",
            "The continuous theorem and full-kernel transfer fit strictly inside the 4900/k^2 endpoint target.",
            exact["conditional_full_ceiling"],
            "Conditional only on the continuous first-summand ceiling.",
        ),
        BridgeRow(
            "co9fscb_08_endpoint_tail",
            "conditional_theorem",
            "ready_to_apply",
            "The continuous theorem would close the order-nine endpoint tail after the two-index splice.",
            exact["conditional_endpoint_tail"],
            "Conditional only on the continuous first-summand ceiling.",
        ),
        BridgeRow(
            "co9fscb_09_open_continuous",
            "open_handoff",
            "not_ready_to_apply",
            "Prove the continuous sixth-nested first-summand curvature theorem.",
            exact["continuous_target"],
            "This continuous inequality is not proved here.",
        ),
        BridgeRow(
            "co9fscb_10_open_splice",
            "open_handoff",
            "not_ready_to_apply",
            "Extend the rigorous endpoint prefix across the two indices before the analytic bridge starts.",
            exact["finite_splice"],
            "These two finite signs are not proved here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order9_first_summand_curvature_bridge",
        "date": "2026-07-13",
        "status": (
            "exact order-nine first/full curvature transfer with one open "
            "continuous theorem and a two-index finite splice"
        ),
        "proof_boundary": (
            "This artifact proves the sixth-gap floor, full-kernel transfer, "
            "and conditional composition. It does not prove w_1''<=4200/t^2, "
            "the two-index splice, order-nine entry, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 8,
            "open_rows": 2,
            "sixth_gap_floor_theorems": 1,
            "full_kernel_transfer_theorems": 1,
            "conditional_rows": 3,
            "open_continuous_targets": 1,
            "open_finite_splices": 1,
            "finite_splice_indices": 2,
            "transfer_polynomial_degree": exact["transfer_polynomial"]["degree"],
            "positive_transfer_coefficients": exact["transfer_polynomial"]["coefficient_count"],
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_first_summand_curvature_bridge.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order9_first_summand_curvature_bridge.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Nine First-Summand Curvature Bridge",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact first/full-kernel transfer with one open continuous",
        "first-summand theorem and a two-index finite splice. This is not a",
        "proof of order-nine entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order9_first_summand_curvature_bridge.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order9_first_summand_curvature_bridge.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order9_first_summand_curvature_bridge.py",
        "```",
        "",
        "## Sixth Stable Layer",
        "",
        "```text",
        exact["continuous_coordinates"],
        exact["order9_curvature"],
        exact["discrete_identity"],
        exact["sixth_gap_floor"],
        exact["stable_log_lipschitz"],
        "first shifted numerator: "
        + exact["first_V_floor_polynomial"]["shifted_polynomial"]
        + ">0,",
        "full shifted numerator: "
        + exact["full_V_floor_polynomial"]["shifted_polynomial"]
        + ">0.",
        "```",
        "",
        "## Full-Kernel Transfer",
        "",
        "The rebalanced `2/k^10` moment-tail theorem propagates through all",
        "six stable logarithms. The final two stages are",
        "",
        "```text",
        exact["sixth_gap_error"],
        exact["order9_coordinate_error"],
        exact["full_transfer"],
        "1251^2*transfer_error=" + exact["endpoint_scaled_transfer_decimal"],
        "reserve below 550=" + exact["endpoint_scaled_reserve_below_550"],
        "```",
        "",
        "The degree-168 shifted numerator has 169 strictly positive",
        "coefficients, so the transfer bound holds on the complete half-line.",
        "",
        "## Remaining Targets",
        "",
        "```text",
        exact["continuous_target"],
        exact["tent_transfer"],
        exact["conditional_full_ceiling"],
        exact["conditional_endpoint_tail"],
        exact["finite_splice"],
        "```",
        "",
        "Thus the continuous theorem would contribute fewer than 4201 inverse",
        "squares and the full-kernel transfer fewer than 550. Their total 4751",
        "fits strictly inside 4900. The continuous theorem and the two finite",
        "signs remain open.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_first_summand_power10_rebalanced_dominance_extension.md",
        "outputs/jensen_window_pf_compound_order8_m100_sharp_tail_curvature_certificate.md",
        "outputs/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.md",
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
        "wrote order-nine first/full curvature bridge: "
        f"{summary['rows']} rows, "
        f"{summary['positive_transfer_coefficients']} positive transfer coefficients, "
        f"{summary['open_continuous_targets']} open continuous target, "
        f"{summary['finite_splice_indices']} finite splice indices"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
