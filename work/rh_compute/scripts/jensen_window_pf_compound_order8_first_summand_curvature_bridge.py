#!/usr/bin/env python3
"""Transfer the order-eight fifth-nested curvature from first summand to full kernel."""

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
    "jensen_window_pf_compound_order8_first_summand_curvature_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order8_first_summand_curvature_bridge.md"
)
TAIL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_m100_tail_curvature_reduction.json"
)
POWER9_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_power9_rebalanced_dominance_extension.json"
)
ORDER7_ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_m100_entry_certificate.json"
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
TAIL_FIRST_K = 1250
GAP_FLOOR_START = TAIL_FIRST_K - 1
CONTINUOUS_FIRST_T = 999
FIRST_CONTINUOUS_CONSTANT = 4000
FULL_TRANSFER_CONSTANT = 190
DISCRETE_FIRST_CONSTANT = 4001
FULL_CURVATURE_CONSTANT = 4300


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
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_sources() -> dict:
    tail = load_json(TAIL_SOURCE)
    power9 = load_json(POWER9_SOURCE)
    order7 = load_json(ORDER7_ENTRY_SOURCE)
    order4_bridge = load_json(ORDER4_BRIDGE_SOURCE)
    order4_entry = load_json(ORDER4_ENTRY_SOURCE)
    if tail.get("exact", {}).get("sufficient_ceiling") != (
        "W_k<=4300/k^2 for every real/integer k>=1250"
    ):
        raise RuntimeError("order-eight tail reduction contract changed")
    if power9.get("summary", {}).get("full_tail_power") != 9:
        raise RuntimeError("power-nine dominance source changed")
    if power9.get("summary", {}).get("tail_start_k") != 300:
        raise RuntimeError("power-nine dominance start changed")
    if order7.get("exact", {}).get("global_first_curvature") != (
        "r_1''(t)<=600/t^2 for every real t>=320"
    ):
        raise RuntimeError("completed order-seven first curvature changed")
    if order7.get("exact", {}).get("full_ceiling") != (
        "R_k<R_k^(1)+|R_k-R_k^(1)|<601/k^2+262/k^2=863/k^2<900/k^2, k>=321"
    ):
        raise RuntimeError("completed order-seven full ceiling changed")
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
        "power9_bound": power9["diagnostics"]["full_tail_relative_bound"],
        "order7_first_curvature": order7["exact"]["global_first_curvature"],
        "order7_full_ceiling": order7["exact"]["full_ceiling"],
        "first_B_floor": order4_bridge["exact"]["integer_B_bounds"],
        "full_B_floor": "B(j)>=d_j>=251/(250*(2*j+1)), j>=320",
    }


def rationalized_error_chain() -> dict[str, sp.Expr]:
    j = sp.symbols("j", integer=True, positive=True)

    def shift(expression: sp.Expr, offset: int) -> sp.Expr:
        return expression.subs(j, j + offset)

    def reduce(expression: sp.Expr) -> sp.Expr:
        return sp.cancel(expression)

    a = reduce(2 * ((j - 1) ** -9 + 2 * j**-9 + (j + 1) ** -9))
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
    total_curvature = reduce(
        shift(order8_coordinate, -1)
        + 2 * order8_coordinate
        + shift(order8_coordinate, 1)
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
        "total_curvature": total_curvature,
    }


def exact_diagnostics() -> dict:
    errors = rationalized_error_chain()
    j = errors["symbol"]
    k = sp.symbols("k", integer=True, positive=True)

    first_u_floor = shifted_positive_polynomial(
        6 / (2 * j + 1) - 601 / j**2 - sp.Rational(3, 2) / j,
        j,
        GAP_FLOOR_START,
    )
    full_u_floor = shifted_positive_polynomial(
        sp.Rational(753, 125) / (2 * j + 1)
        - 863 / j**2
        - sp.Rational(3, 2) / j,
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
        raise RuntimeError("order-eight transfer polynomial is not positive")
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
            "U(t)=6*B(t)-r(t-1)+2*r(t)-r(t+1); "
            "s(t)=2*r(t)-p(t)+log(1-exp(-U(t)))"
        ),
        "order8_curvature": (
            "s''=2*r''-p''+phi(U)*U''-chi(U)*(U')^2; "
            "phi(z)=1/(exp(z)-1), chi(z)=exp(z)/(exp(z)-1)^2"
        ),
        "discrete_identity": (
            "W_k=s(k-1)-2*s(k)+s(k+1)="
            "integral_[-1,1](1-|v|)*s''(k+v) dv"
        ),
        "moment_error": (
            "a_j=2*((j-1)^(-9)+2*j^(-9)+(j+1)^(-9)); "
            "|B_j-B_j^(1)|<=a_j"
        ),
        "log_defect_error": "L_j=4*j*a_j; |ell_j-ell_j^(1)|<=L_j",
        "first_gap_error": (
            "U1_j=2*a_j+L_(j-1)+2*L_j+L_(j+1); "
            "|J_j-J_j^(1)|<=U1_j"
        ),
        "log_first_gap_error": (
            "V_j=2*L_j+8*j*U1_j; |h_j-h_j^(1)|<=V_j"
        ),
        "second_gap_error": (
            "W1_j=3*a_j+V_(j-1)+2*V_j+V_(j+1); "
            "|R_j-R_j^(1)|<=W1_j"
        ),
        "order5_coordinate_error": (
            "E_j=2*V_j+(5*j/7)*W1_j+L_j; |q_j-q_j^(1)|<=E_j"
        ),
        "third_gap_error": (
            "Z_j=4*a_j+E_(j-1)+2*E_j+E_(j+1); "
            "|S_j-S_j^(1)|<=Z_j"
        ),
        "order6_coordinate_error": (
            "Y_j=2*E_j+V_j+(2*j/3)*Z_j; |p_j-p_j^(1)|<=Y_j"
        ),
        "fourth_gap_error": (
            "O_j=5*a_j+Y_(j-1)+2*Y_j+Y_(j+1); "
            "|T_j-T_j^(1)|<=O_j"
        ),
        "order7_coordinate_error": (
            "N_j=2*Y_j+E_j+(2*j/3)*O_j; |r_j-r_j^(1)|<=N_j"
        ),
        "fifth_gap_error": (
            "P_j=6*a_j+N_(j-1)+2*N_j+N_(j+1); "
            "|U_j-U_j^(1)|<=P_j"
        ),
        "fifth_gap_floor": "min(U_j,U_j^(1))>=3/(2*j), j>=1249",
        "order8_coordinate_error": (
            "C_j=2*N_j+Y_j+(2*j/3)*P_j; |s_j-s_j^(1)|<=C_j"
        ),
        "full_transfer": (
            "|W_k-W_k^(1)|<=C_(k-1)+2*C_k+C_(k+1)<190/k^2, k>=1250"
        ),
        "continuous_target": "s_1''(t)<=4000/t^2 for every real t>=999",
        "tent_transfer": (
            "s_1''(t)<=4000/t^2 => W_k^(1)<=4000*[-log(1-1/k^2)]"
            "<4001/k^2, k>=1250"
        ),
        "conditional_full_ceiling": (
            "s_1''(t)<=4000/t^2 on t>=999 => "
            "W_k<4001/k^2+190/k^2=4191/k^2<4300/k^2, k>=1250"
        ),
        "conditional_endpoint_tail": (
            "s_1''(t)<=4000/t^2 on t>=999 => "
            "Q_(8,n)(-100)>0 for every n>=1243"
        ),
        "first_U_floor_polynomial": first_u_floor,
        "full_U_floor_polynomial": full_u_floor,
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
        "endpoint_scaled_reserve_below_190": str(
            sp.N(FULL_TRANSFER_CONSTANT - endpoint_scaled, 30)
        ),
        "stable_derivative_residuals": [str(first_residual), str(second_residual)],
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_diagnostics()
    rows = [
        BridgeRow(
            "co8fscb_01_coordinate",
            "exact_identity",
            "ready_to_apply",
            "Order eight adds one fifth nested stable logarithm to the completed hierarchy.",
            exact["continuous_coordinates"] + "; " + exact["order8_curvature"],
            "Exact differentiation only.",
        ),
        BridgeRow(
            "co8fscb_02_power9_input",
            "theorem_input",
            "ready_to_apply",
            "The rebalanced first-summand theorem supplies the perturbation power required by the fifth stable logarithm.",
            sources["power9_bound"],
            "Completed lambda=-100 moment-tail theorem.",
        ),
        BridgeRow(
            "co8fscb_03_fifth_floor",
            "exact_analytic_composition",
            "ready_to_apply",
            "First and full fifth-gap coordinates stay uniformly away from zero on the transfer tail.",
            exact["fifth_gap_floor"],
            "Composes completed order-seven first/full curvature and defect theorems.",
            {
                "first": exact["first_U_floor_polynomial"],
                "full": exact["full_U_floor_polynomial"],
            },
        ),
        BridgeRow(
            "co8fscb_04_error_chain",
            "exact_error_propagation",
            "ready_to_apply",
            "The ninth-power moment error survives all five stable logarithms as an inverse-square curvature error.",
            exact["moment_error"] + "; " + exact["order8_coordinate_error"],
            "Exact triangle and mean-value inequalities inside proved positive floors.",
        ),
        BridgeRow(
            "co8fscb_05_full_transfer",
            "analytic_theorem",
            "ready_to_apply",
            "The complete theta kernel changes the order-eight centered curvature by fewer than 190 inverse squares.",
            exact["full_transfer"],
            "Exact degree-133 coefficient-positive rational audit.",
            {
                "endpoint_scaled": exact["endpoint_scaled_transfer_decimal"],
                "reserve": exact["endpoint_scaled_reserve_below_190"],
            },
        ),
        BridgeRow(
            "co8fscb_06_tent",
            "conditional_theorem",
            "ready_to_apply",
            "A continuous first-summand ceiling contributes fewer than 4001 inverse squares discretely.",
            exact["tent_transfer"],
            "Conditional only on the displayed continuous curvature theorem.",
        ),
        BridgeRow(
            "co8fscb_07_full_ceiling",
            "conditional_theorem",
            "ready_to_apply",
            "The continuous theorem and full-kernel transfer fit strictly inside the 4300/k^2 endpoint target.",
            exact["conditional_full_ceiling"],
            "Conditional only on the continuous first-summand ceiling.",
        ),
        BridgeRow(
            "co8fscb_08_endpoint_tail",
            "conditional_theorem",
            "ready_to_apply",
            "The continuous theorem would close the complete missing order-eight endpoint tail.",
            exact["conditional_endpoint_tail"],
            "Conditional only on the continuous first-summand ceiling.",
        ),
        BridgeRow(
            "co8fscb_09_open_continuous",
            "open_handoff",
            "not_ready_to_apply",
            "Prove the continuous fifth-nested first-summand curvature theorem.",
            exact["continuous_target"],
            "This continuous inequality is not proved here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order8_first_summand_curvature_bridge",
        "date": "2026-07-13",
        "status": "exact order-eight first/full curvature transfer with one open continuous theorem",
        "proof_boundary": (
            "This artifact proves the full-kernel transfer and conditional composition. "
            "It does not prove s_1''<=4000/t^2, order-eight entry, PF-infinity, "
            "RH, or Lambda<=0."
        ),
        "sources": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 8,
            "open_rows": 1,
            "fifth_gap_floor_theorems": 1,
            "full_kernel_transfer_theorems": 1,
            "conditional_rows": 3,
            "open_continuous_targets": 1,
            "transfer_polynomial_degree": exact["transfer_polynomial"]["degree"],
            "positive_transfer_coefficients": exact["transfer_polynomial"]["coefficient_count"],
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order8_first_summand_curvature_bridge.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order8_first_summand_curvature_bridge.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Eight First-Summand Curvature Bridge",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact first/full-kernel transfer with one open continuous",
        "first-summand theorem. This is not a proof of order-eight entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order8_first_summand_curvature_bridge.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order8_first_summand_curvature_bridge.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_first_summand_curvature_bridge.py",
        "```",
        "",
        "## Fifth Stable Layer",
        "",
        "```text",
        exact["continuous_coordinates"],
        exact["order8_curvature"],
        exact["discrete_identity"],
        "```",
        "",
        "The completed order-seven first/full bounds prove",
        "",
        "```text",
        exact["fifth_gap_floor"],
        "first shifted numerator: " + exact["first_U_floor_polynomial"]["shifted_polynomial"] + ">0,",
        "full shifted numerator: " + exact["full_U_floor_polynomial"]["shifted_polynomial"] + ">0.",
        "```",
        "",
        "## Full-Kernel Transfer",
        "",
        "The rebalanced 2/k^9 moment-tail theorem propagates as",
        "",
        "```text",
        exact["moment_error"],
        exact["log_defect_error"],
        exact["first_gap_error"],
        exact["log_first_gap_error"],
        exact["second_gap_error"],
        exact["order5_coordinate_error"],
        exact["third_gap_error"],
        exact["order6_coordinate_error"],
        exact["fourth_gap_error"],
        exact["order7_coordinate_error"],
        exact["fifth_gap_error"],
        exact["order8_coordinate_error"],
        exact["full_transfer"],
        "```",
        "",
        "At the splice,",
        "",
        "```text",
        "1250^2*transfer_error=" + exact["endpoint_scaled_transfer_decimal"],
        "reserve below 190=" + exact["endpoint_scaled_reserve_below_190"],
        "```",
        "",
        "The degree-133 shifted numerator has 134 strictly positive",
        "coefficients, so the transfer bound holds on the complete half-line.",
        "",
        "## Remaining Continuous Theorem",
        "",
        "```text",
        exact["continuous_target"],
        exact["tent_transfer"],
        exact["conditional_full_ceiling"],
        exact["conditional_endpoint_tail"],
        "```",
        "",
        "Thus the continuous theorem would contribute fewer than 4001 inverse",
        "squares, the full-kernel transfer fewer than 190, and the total 4191",
        "fits strictly inside the endpoint target 4300. The continuous theorem",
        "itself remains open.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_first_summand_power9_rebalanced_dominance_extension.md",
        "outputs/jensen_window_pf_compound_order8_m100_tail_curvature_reduction.md",
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
        "wrote order-eight first/full curvature bridge: "
        f"{summary['rows']} rows, "
        f"{summary['positive_transfer_coefficients']} positive transfer coefficients, "
        f"{summary['open_continuous_targets']} open continuous target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
