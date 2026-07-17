#!/usr/bin/env python3
"""Transfer the order-seven fourth-nested curvature from first summand to full kernel."""

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
    "jensen_window_pf_compound_order7_first_summand_curvature_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md"
)
TAIL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_m100_tail_curvature_reduction.json"
)
POWER8_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.json"
)
ORDER6_ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_m100_entry_certificate.json"
)
ORDER6_PREFIX_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_m100_prefix_certificate.json"
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
TAIL_FIRST_K = 321
CONTINUOUS_FIRST_T = 320
FIRST_CONTINUOUS_CONSTANT = 600
FULL_TRANSFER_CONSTANT = 262
DISCRETE_FIRST_CONSTANT = 601
FULL_CURVATURE_CONSTANT = 900


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
    power8 = load_json(POWER8_SOURCE)
    order6 = load_json(ORDER6_ENTRY_SOURCE)
    order6_prefix = load_json(ORDER6_PREFIX_SOURCE)
    order4_bridge = load_json(ORDER4_BRIDGE_SOURCE)
    order4_entry = load_json(ORDER4_ENTRY_SOURCE)
    if tail.get("exact", {}).get("sufficient_ceiling") != (
        "R_k<=900/k^2 for every real/integer k>=321"
    ):
        raise RuntimeError("order-seven tail reduction contract changed")
    if power8.get("summary", {}).get("full_tail_power") != 8:
        raise RuntimeError("power-eight dominance source changed")
    if power8.get("summary", {}).get("tail_start_k") != 300:
        raise RuntimeError("power-eight dominance start changed")
    if order6.get("exact", {}).get("global_first_curvature") != (
        "p_1''(t)<=200/t^2 for every real t>=321"
    ):
        raise RuntimeError("completed order-six first curvature changed")
    if order6.get("exact", {}).get("full_ceiling") != (
        "P_k<=P_k^(1)+|P_k-P_k^(1)|<201/k^2+100/k^2=301/k^2<320/k^2, k>=322"
    ):
        raise RuntimeError("completed order-six full ceiling changed")
    if order4_bridge.get("exact", {}).get("integer_B_bounds") != (
        "1/(2*m+1)<=B_1(m)<=3/(2*m-1), integer m>=319"
    ):
        raise RuntimeError("first-summand B floor changed")
    if order4_entry.get("tail_arithmetic", {}).get("defect_buffer") != (
        "-3*log(x_k)>3*d_k>=753/(250*(2*k+1))"
    ):
        raise RuntimeError("full-kernel defect floor changed")
    finite_rows = {
        int(row["n"]): row
        for row in order6_prefix.get("finite", {}).get("rows", [])
    }
    endpoint_lower = {}
    for n in (315, 316):
        row = finite_rows.get(n)
        if row is None:
            raise RuntimeError(f"missing order-six finite endpoint n={n}")
        endpoint_lower[str(n + 5)] = row["relative_H5_margin_lower"]
    return {
        "tail_ceiling": tail["exact"]["sufficient_ceiling"],
        "power8_bound": power8["diagnostics"]["full_tail_relative_bound"],
        "order6_first_curvature": order6["exact"]["global_first_curvature"],
        "order6_full_ceiling": order6["exact"]["full_ceiling"],
        "first_B_floor": order4_bridge["exact"]["integer_B_bounds"],
        "full_B_floor": "B(j)>=d_j>=251/(250*(2*j+1)), j>=320",
        "finite_relative_H5_lower": endpoint_lower,
    }


def rationalized_error_chain() -> dict[str, sp.Expr]:
    j = sp.symbols("j", integer=True, positive=True)

    def shift(expression: sp.Expr, offset: int) -> sp.Expr:
        return expression.subs(j, j + offset)

    def reduce(expression: sp.Expr) -> sp.Expr:
        return sp.cancel(expression)

    a = reduce(2 * ((j - 1) ** -8 + 2 * j**-8 + (j + 1) ** -8))
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
    total_curvature = reduce(
        shift(order7_coordinate, -1)
        + 2 * order7_coordinate
        + shift(order7_coordinate, 1)
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
        "total_curvature": total_curvature,
    }


def exact_diagnostics(sources: dict) -> dict:
    errors = rationalized_error_chain()
    j = errors["symbol"]
    k = sp.symbols("k", integer=True, positive=True)

    first_t_floor = shifted_positive_polynomial(
        5 / (2 * j + 1) - 201 / j**2 - sp.Rational(3, 2) / j,
        j,
        322,
    )
    full_t_floor = shifted_positive_polynomial(
        sp.Rational(251, 50) / (2 * j + 1)
        - 301 / j**2
        - sp.Rational(3, 2) / j,
        j,
        322,
    )

    finite_floor_rows = []
    for endpoint_j_text, lower_text in sorted(
        sources["finite_relative_H5_lower"].items(), key=lambda item: int(item[0])
    ):
        endpoint_j = int(endpoint_j_text)
        relative_lower = sp.Rational(lower_text)
        full_t_lower = sp.cancel(relative_lower / (1 + relative_lower))
        fourth_error = sp.cancel(errors["fourth_gap"].subs(j, endpoint_j))
        first_t_lower = sp.cancel(full_t_lower - fourth_error)
        target = sp.Rational(3, 2 * endpoint_j)
        full_margin = sp.cancel(full_t_lower - target)
        first_margin = sp.cancel(first_t_lower - target)
        if full_margin <= 0 or first_margin <= 0:
            raise RuntimeError(f"finite fourth-gap floor failed at j={endpoint_j}")
        finite_floor_rows.append(
            {
                "j": endpoint_j,
                "relative_H5_lower": lower_text,
                "full_T_lower": str(full_t_lower),
                "full_T_lower_decimal": str(sp.N(full_t_lower, 30)),
                "fourth_gap_error": str(fourth_error),
                "fourth_gap_error_decimal": str(sp.N(fourth_error, 30)),
                "first_T_lower": str(first_t_lower),
                "first_T_lower_decimal": str(sp.N(first_t_lower, 30)),
                "target": str(target),
                "full_margin": str(full_margin),
                "first_margin": str(first_margin),
                "full_margin_decimal": str(sp.N(full_margin, 30)),
                "first_margin_decimal": str(sp.N(first_margin, 30)),
            }
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
        raise RuntimeError("order-seven transfer polynomial is not positive")
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
            "T(t)=5*B(t)-p(t-1)+2*p(t)-p(t+1); "
            "r(t)=2*p(t)-q(t)+log(1-exp(-T(t)))"
        ),
        "order7_curvature": (
            "r''=2*p''-q''+phi(T)*T''-chi(T)*(T')^2; "
            "phi(z)=1/(exp(z)-1), chi(z)=exp(z)/(exp(z)-1)^2"
        ),
        "discrete_identity": (
            "R_k=r(k-1)-2*r(k)+r(k+1)="
            "integral_[-1,1](1-|s|)*r''(k+s) ds"
        ),
        "moment_error": (
            "a_j=2*((j-1)^(-8)+2*j^(-8)+(j+1)^(-8)); "
            "|B_j-B_j^(1)|<=a_j"
        ),
        "log_defect_error": "L_j=4*j*a_j; |ell_j-ell_j^(1)|<=L_j",
        "first_gap_error": (
            "U_j=2*a_j+L_(j-1)+2*L_j+L_(j+1); "
            "|J_j-J_j^(1)|<=U_j"
        ),
        "log_first_gap_error": (
            "V_j=2*L_j+8*j*U_j; |h_j-h_j^(1)|<=V_j"
        ),
        "second_gap_error": (
            "W_j=3*a_j+V_(j-1)+2*V_j+V_(j+1); "
            "|R_j-R_j^(1)|<=W_j"
        ),
        "order5_coordinate_error": (
            "E_j=2*V_j+(5*j/7)*W_j+L_j; |q_j-q_j^(1)|<=E_j"
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
        "fourth_gap_floor": (
            "min(T_j,T_j^(1))>=3/(2*j), j>=320"
        ),
        "order7_coordinate_error": (
            "N_j=2*Y_j+E_j+(2*j/3)*O_j; |r_j-r_j^(1)|<=N_j"
        ),
        "full_transfer": (
            "|R_k-R_k^(1)|<=N_(k-1)+2*N_k+N_(k+1)<262/k^2, k>=321"
        ),
        "continuous_target": "r_1''(t)<=600/t^2 for every real t>=320",
        "tent_transfer": (
            "r_1''(t)<=600/t^2 => R_k^(1)<=600*[-log(1-1/k^2)]"
            "<601/k^2, k>=321"
        ),
        "conditional_full_ceiling": (
            "r_1''(t)<=600/t^2 on t>=320 => R_k<863/k^2<900/k^2 "
            "for k>=321"
        ),
        "conditional_endpoint_tail": (
            "r_1''(t)<=600/t^2 on t>=320 => Q_(7,n)(-100)>0 for every n>=315"
        ),
        "first_T_floor_polynomial": first_t_floor,
        "full_T_floor_polynomial": full_t_floor,
        "finite_T_floor_rows": finite_floor_rows,
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
        "endpoint_scaled_reserve_below_262": str(
            sp.N(FULL_TRANSFER_CONSTANT - endpoint_scaled, 30)
        ),
        "stable_derivative_residuals": [str(first_residual), str(second_residual)],
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_diagnostics(sources)
    rows = [
        BridgeRow(
            "co7fscb_01_coordinate",
            "exact_identity",
            "ready_to_apply",
            "Order seven adds one fourth nested stable logarithm to the completed hierarchy.",
            exact["continuous_coordinates"] + "; " + exact["order7_curvature"],
            "Exact differentiation only.",
        ),
        BridgeRow(
            "co7fscb_02_power8_input",
            "theorem_input",
            "ready_to_apply",
            "The rebalanced first-summand theorem supplies the perturbation power required by the fourth stable logarithm.",
            sources["power8_bound"],
            "Completed lambda=-100 moment-tail theorem.",
        ),
        BridgeRow(
            "co7fscb_03_fourth_floor",
            "exact_interval_composition",
            "ready_to_apply",
            "First and full fourth-gap coordinates stay uniformly away from zero on the complete transfer tail.",
            exact["fourth_gap_floor"],
            "Composes two finite order-six margins with the global order-six curvature and defect theorems.",
            {
                "finite": exact["finite_T_floor_rows"],
                "first_tail": exact["first_T_floor_polynomial"],
                "full_tail": exact["full_T_floor_polynomial"],
            },
        ),
        BridgeRow(
            "co7fscb_04_error_chain",
            "exact_error_propagation",
            "ready_to_apply",
            "The eighth-power moment error survives all four stable logarithms as an inverse-square coordinate error.",
            exact["moment_error"] + "; " + exact["order7_coordinate_error"],
            "Exact triangle and mean-value inequalities inside proved positive floors.",
        ),
        BridgeRow(
            "co7fscb_05_full_transfer",
            "analytic_theorem",
            "ready_to_apply",
            "The complete theta kernel changes the order-seven centered curvature by fewer than 262 inverse squares.",
            exact["full_transfer"],
            "Exact degree-102 coefficient-positive rational audit.",
            {
                "endpoint_scaled": exact["endpoint_scaled_transfer_decimal"],
                "reserve": exact["endpoint_scaled_reserve_below_262"],
            },
        ),
        BridgeRow(
            "co7fscb_06_tent",
            "conditional_theorem",
            "ready_to_apply",
            "A continuous first-summand ceiling contributes fewer than 601 inverse squares discretely.",
            exact["tent_transfer"],
            "Conditional only on the displayed continuous curvature theorem.",
        ),
        BridgeRow(
            "co7fscb_07_full_ceiling",
            "conditional_theorem",
            "ready_to_apply",
            "The continuous theorem and full-kernel transfer fit strictly inside the 900/k^2 endpoint target.",
            exact["conditional_full_ceiling"],
            "Conditional only on the continuous first-summand ceiling.",
        ),
        BridgeRow(
            "co7fscb_08_endpoint_tail",
            "conditional_theorem",
            "ready_to_apply",
            "The continuous theorem would close the complete missing order-seven endpoint tail.",
            exact["conditional_endpoint_tail"],
            "Conditional only on the continuous first-summand ceiling.",
        ),
        BridgeRow(
            "co7fscb_09_open_continuous",
            "open_handoff",
            "not_ready_to_apply",
            "Prove the continuous fourth-nested first-summand curvature theorem.",
            exact["continuous_target"],
            "This continuous inequality is not proved here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order7_first_summand_curvature_bridge",
        "date": "2026-07-13",
        "status": "exact order-seven first/full curvature transfer with one open continuous theorem",
        "proof_boundary": (
            "This artifact proves the full-kernel transfer and conditional composition. "
            "It does not prove r_1''<=600/t^2, order-seven entry, PF-infinity, "
            "RH, or Lambda<=0."
        ),
        "sources": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 8,
            "open_rows": 1,
            "fourth_gap_floor_theorems": 1,
            "full_kernel_transfer_theorems": 1,
            "conditional_rows": 3,
            "open_continuous_targets": 1,
            "transfer_polynomial_degree": exact["transfer_polynomial"]["degree"],
            "positive_transfer_coefficients": exact["transfer_polynomial"]["coefficient_count"],
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order7_first_summand_curvature_bridge.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order7_first_summand_curvature_bridge.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    finite = exact["finite_T_floor_rows"]
    lines = [
        "# Jensen-Window PF Compound Order-Seven First-Summand Curvature Bridge",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact first/full-kernel transfer with one open continuous",
        "first-summand theorem. This is not a proof of order-seven entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order7_first_summand_curvature_bridge.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order7_first_summand_curvature_bridge.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order7_first_summand_curvature_bridge.py",
        "```",
        "",
        "## Fourth Stable Layer",
        "",
        "```text",
        exact["continuous_coordinates"],
        exact["order7_curvature"],
        exact["discrete_identity"],
        "```",
        "",
        "For `j>=322`, the completed first and full order-six curvature bounds",
        "and the first/full defect floors give",
        "",
        "```text",
        exact["fourth_gap_floor"],
        "first shifted numerator: 4*m^2+1769*m+154480>0,",
        "full shifted numerator: 101*m^2+34869*m+740684>0.",
        "```",
        "",
        "The two missing floor indices are rigorous finite compositions:",
        "",
        "```text",
        f"j=320: full margin={finite[0]['full_margin_decimal']}, first margin={finite[0]['first_margin_decimal']},",
        f"j=321: full margin={finite[1]['full_margin_decimal']}, first margin={finite[1]['first_margin_decimal']}.",
        "```",
        "",
        "Here the full gap is `T_j=log(1+K_(j-5))`, bounded below by",
        "`K/(1+K)`, and the first gap is within `O_j` of it.",
        "",
        "## Full-Kernel Transfer",
        "",
        "The rebalanced `2/k^8` moment-tail theorem propagates as",
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
        exact["full_transfer"],
        "```",
        "",
        "At the splice,",
        "",
        "```text",
        "321^2*transfer_error=" + exact["endpoint_scaled_transfer_decimal"],
        "reserve below 262=" + exact["endpoint_scaled_reserve_below_262"],
        "```",
        "",
        "The degree-102 shifted numerator has 103 strictly positive",
        "coefficients, so the transfer bound holds on the entire half-line.",
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
        "Thus the continuous theorem would contribute fewer than 601 inverse",
        "squares, the full-kernel transfer fewer than 262, and the total 863",
        "fits strictly inside the endpoint target 900. The continuous theorem",
        "itself remains open.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.md",
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
        "wrote order-seven first/full curvature bridge: "
        f"{summary['rows']} rows, "
        f"{summary['positive_transfer_coefficients']} positive transfer coefficients, "
        f"{summary['open_continuous_targets']} open continuous target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
