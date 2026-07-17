#!/usr/bin/env python3
"""Transfer the order-six nested curvature from the first summand to the full kernel."""

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
    "jensen_window_pf_compound_order6_first_summand_curvature_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order6_first_summand_curvature_bridge.md"
)
TAIL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_m100_tail_curvature_reduction.json"
)
POWER7_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.json"
)
ORDER5_ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_m100_entry_certificate.json"
)
ORDER5_BRIDGE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_first_summand_curvature_bridge.json"
)
ORDER4_BRIDGE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_first_summand_curvature_bridge.json"
)
TAIL_FIRST_K = 322
CONTINUOUS_FIRST_T = 321
FIRST_CONTINUOUS_CONSTANT = 200
FULL_TRANSFER_CONSTANT = 100
DISCRETE_FIRST_CONSTANT = 201
FULL_CURVATURE_CONSTANT = 320


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
    power7 = load_json(POWER7_SOURCE)
    order5 = load_json(ORDER5_ENTRY_SOURCE)
    order5_bridge = load_json(ORDER5_BRIDGE_SOURCE)
    order4_bridge = load_json(ORDER4_BRIDGE_SOURCE)
    if tail.get("exact", {}).get("sufficient_ceiling") != (
        "P_k<=320/k^2 for every real/integer k>=322"
    ):
        raise RuntimeError("order-six tail reduction contract changed")
    if power7.get("summary", {}).get("full_tail_power") != 7:
        raise RuntimeError("power-seven dominance source changed")
    if power7.get("summary", {}).get("tail_start_k") != 316:
        raise RuntimeError("power-seven dominance start changed")
    if order5.get("exact", {}).get("global_first_curvature") != (
        "q_1''(t)<=60/t^2 for every real t>=320"
    ):
        raise RuntimeError("completed order-five first curvature changed")
    if order5.get("exact", {}).get("full_ceiling") != (
        "C_n<=C_n^(1)+|C_n-C_n^(1)|<63/k^2+37/k^2=100/k^2, k>=321"
    ):
        raise RuntimeError("completed order-five full ceiling changed")
    if order5_bridge.get("exact", {}).get("second_gap_floor") != (
        "min(R_j,R_j^(1))>=7/(5*j), j>=320"
    ):
        raise RuntimeError("order-five second-gap floor changed")
    if order4_bridge.get("exact", {}).get("integer_B_bounds") != (
        "1/(2*m+1)<=B_1(m)<=3/(2*m-1), integer m>=319"
    ):
        raise RuntimeError("first-summand B floor changed")
    return {
        "tail_ceiling": tail["exact"]["sufficient_ceiling"],
        "power7_bound": power7["diagnostics"]["full_tail_relative_bound"],
        "order5_first_curvature": order5["exact"]["global_first_curvature"],
        "order5_full_ceiling": order5["exact"]["full_ceiling"],
        "second_gap_floor": order5_bridge["exact"]["second_gap_floor"],
        "first_B_floor": order4_bridge["exact"]["integer_B_bounds"],
    }


def moment_error(index: sp.Expr) -> sp.Expr:
    return 2 * (
        1 / (index - 1) ** 7
        + 2 / index**7
        + 1 / (index + 1) ** 7
    )


def log_defect_error(index: sp.Expr) -> sp.Expr:
    return 4 * index * moment_error(index)


def first_gap_error(index: sp.Expr) -> sp.Expr:
    return (
        2 * moment_error(index)
        + log_defect_error(index - 1)
        + 2 * log_defect_error(index)
        + log_defect_error(index + 1)
    )


def log_first_gap_error(index: sp.Expr) -> sp.Expr:
    return 2 * log_defect_error(index) + 8 * index * first_gap_error(index)


def second_gap_error(index: sp.Expr) -> sp.Expr:
    return (
        3 * moment_error(index)
        + log_first_gap_error(index - 1)
        + 2 * log_first_gap_error(index)
        + log_first_gap_error(index + 1)
    )


def order5_coordinate_error(index: sp.Expr) -> sp.Expr:
    return (
        2 * log_first_gap_error(index)
        + sp.Rational(5, 7) * index * second_gap_error(index)
        + log_defect_error(index)
    )


def third_gap_error(index: sp.Expr) -> sp.Expr:
    return (
        4 * moment_error(index)
        + order5_coordinate_error(index - 1)
        + 2 * order5_coordinate_error(index)
        + order5_coordinate_error(index + 1)
    )


def order6_coordinate_error(index: sp.Expr) -> sp.Expr:
    return (
        2 * order5_coordinate_error(index)
        + log_first_gap_error(index)
        + sp.Rational(2, 3) * index * third_gap_error(index)
    )


def exact_diagnostics() -> dict:
    j, k = sp.symbols("j k", integer=True, positive=True)
    first_s_floor = shifted_positive_polynomial(
        4 / (2 * j + 1) - 63 / j**2 - sp.Rational(3, 2) / j,
        j,
        321,
    )
    full_s_floor = shifted_positive_polynomial(
        sp.Rational(502, 125) / (2 * j + 1)
        - 100 / j**2
        - sp.Rational(3, 2) / j,
        j,
        321,
    )

    total_error = sp.factor(
        order6_coordinate_error(k - 1)
        + 2 * order6_coordinate_error(k)
        + order6_coordinate_error(k + 1)
    )
    transfer = shifted_positive_polynomial(
        FULL_TRANSFER_CONSTANT / k**2 - total_error,
        k,
        TAIL_FIRST_K,
    )
    endpoint_error = sp.factor(total_error.subs(k, TAIL_FIRST_K))
    endpoint_scaled = sp.factor(TAIL_FIRST_K**2 * endpoint_error)

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
            "S(t)=4*B(t)-q(t-1)+2*q(t)-q(t+1); "
            "p(t)=2*q(t)-h(t)+log(1-exp(-S(t)))"
        ),
        "order6_curvature": (
            "p''=2*q''-h''+phi(S)*S''-chi(S)*(S')^2; "
            "phi(z)=1/(exp(z)-1), chi(z)=exp(z)/(exp(z)-1)^2"
        ),
        "discrete_identity": (
            "P_k=p(k-1)-2*p(k)+p(k+1)="
            "integral_[-1,1](1-|s|)*p''(k+s) ds"
        ),
        "moment_error": (
            "a_j=2*((j-1)^(-7)+2*j^(-7)+(j+1)^(-7)); "
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
        "third_gap_floor": (
            "min(S_j,S_j^(1))>=3/(2*j), j>=321"
        ),
        "order6_coordinate_error": (
            "Y_j=2*E_j+V_j+(2*j/3)*Z_j; |p_j-p_j^(1)|<=Y_j"
        ),
        "full_transfer": (
            "|P_k-P_k^(1)|<=Y_(k-1)+2*Y_k+Y_(k+1)<100/k^2, k>=322"
        ),
        "continuous_target": "p_1''(t)<=200/t^2 for every real t>=321",
        "tent_transfer": (
            "p_1''(t)<=200/t^2 => P_k^(1)<=200*[-log(1-1/k^2)]"
            "<201/k^2, k>=322"
        ),
        "conditional_full_ceiling": (
            "p_1''(t)<=200/t^2 on t>=321 => P_k<301/k^2<320/k^2 "
            "for k>=322"
        ),
        "first_S_floor_polynomial": first_s_floor,
        "full_S_floor_polynomial": full_s_floor,
        "transfer_polynomial": transfer,
        "endpoint_transfer_error": str(endpoint_error),
        "endpoint_scaled_transfer": str(endpoint_scaled),
        "endpoint_scaled_transfer_decimal": str(sp.N(endpoint_scaled, 30)),
        "endpoint_scaled_reserve_below_100": str(
            sp.N(FULL_TRANSFER_CONSTANT - endpoint_scaled, 30)
        ),
        "stable_derivative_residuals": [str(first_residual), str(second_residual)],
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_diagnostics()
    rows = [
        BridgeRow(
            "co6fscb_01_coordinate",
            "exact_identity",
            "ready_to_apply",
            "Order six adds one third nested stable logarithm to the completed hierarchy.",
            exact["continuous_coordinates"] + "; " + exact["order6_curvature"],
            "Exact differentiation only.",
        ),
        BridgeRow(
            "co6fscb_02_power7_input",
            "theorem_input",
            "ready_to_apply",
            "The strengthened first-summand theorem supplies the perturbation power required by one more stable logarithm.",
            sources["power7_bound"],
            "Completed lambda=-100 moment-tail theorem.",
        ),
        BridgeRow(
            "co6fscb_03_third_floor",
            "exact_theorem_composition",
            "ready_to_apply",
            "First and full third-gap coordinates stay uniformly away from zero on the transfer tail.",
            exact["third_gap_floor"],
            "Composes the order-five curvature theorem, defect anchor, and coefficient-positive comparisons.",
            {
                "first": exact["first_S_floor_polynomial"],
                "full": exact["full_S_floor_polynomial"],
            },
        ),
        BridgeRow(
            "co6fscb_04_error_chain",
            "exact_error_propagation",
            "ready_to_apply",
            "The seventh-power moment error survives all three stable logarithms as an inverse-square coordinate error.",
            exact["moment_error"] + "; " + exact["order6_coordinate_error"],
            "Exact triangle and mean-value inequalities inside proved positive floors.",
        ),
        BridgeRow(
            "co6fscb_05_full_transfer",
            "analytic_theorem",
            "ready_to_apply",
            "The complete theta kernel changes the order-six centered curvature by less than one hundred inverse squares.",
            exact["full_transfer"],
            "Exact coefficient-positive rational audit.",
            {
                "endpoint_scaled": exact["endpoint_scaled_transfer_decimal"],
                "reserve": exact["endpoint_scaled_reserve_below_100"],
            },
        ),
        BridgeRow(
            "co6fscb_06_tent",
            "conditional_theorem",
            "ready_to_apply",
            "A continuous first-summand ceiling contributes fewer than 201 inverse squares discretely.",
            exact["tent_transfer"],
            "Conditional only on the displayed continuous curvature theorem.",
        ),
        BridgeRow(
            "co6fscb_07_full_ceiling",
            "conditional_theorem",
            "ready_to_apply",
            "The continuous theorem and full-kernel transfer fit strictly inside the 320/k^2 endpoint target.",
            exact["conditional_full_ceiling"],
            "Conditional only on the continuous first-summand ceiling.",
        ),
        BridgeRow(
            "co6fscb_08_open_continuous",
            "open_handoff",
            "not_ready_to_apply",
            "Prove the continuous third-nested first-summand curvature theorem.",
            exact["continuous_target"],
            "This continuous inequality is not proved here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order6_first_summand_curvature_bridge",
        "date": "2026-07-13",
        "status": "exact order-six first/full curvature transfer with one open continuous theorem",
        "proof_boundary": (
            "This artifact proves the full-kernel transfer and conditional composition. "
            "It does not prove p_1''<=200/t^2, order-six entry, PF-infinity, RH, "
            "or Lambda<=0."
        ),
        "sources": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": 7,
            "open_rows": 1,
            "third_gap_floor_theorems": 1,
            "full_kernel_transfer_theorems": 1,
            "conditional_rows": 2,
            "open_continuous_targets": 1,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order6_first_summand_curvature_bridge.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order6_first_summand_curvature_bridge.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Six First-Summand Curvature Bridge",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact first/full-kernel transfer with one open continuous",
        "first-summand theorem. This is not a proof of order-six entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order6_first_summand_curvature_bridge.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order6_first_summand_curvature_bridge.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_first_summand_curvature_bridge.py",
        "```",
        "",
        "## Third Stable Layer",
        "",
        "```text",
        exact["continuous_coordinates"],
        exact["order6_curvature"],
        exact["discrete_identity"],
        "```",
        "",
        "The completed order-five curvature and defect bounds imply",
        "",
        "```text",
        exact["third_gap_floor"],
        "first floor numerator after j=321+m: 2*m^2+1029*m+124101>0,",
        "full floor numerator after j=321+m: 254*m^2+112693*m+9977039>0.",
        "```",
        "",
        "## Full-Kernel Transfer",
        "",
        "The strengthened `2/k^7` moment-tail theorem propagates as",
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
        exact["full_transfer"],
        "```",
        "",
        "At the splice,",
        "",
        "```text",
        "322^2*transfer_error=" + exact["endpoint_scaled_transfer_decimal"],
        "reserve below 100=" + exact["endpoint_scaled_reserve_below_100"],
        "```",
        "",
        "The degree-75 shifted numerator has 76 strictly positive",
        "coefficients, so the transfer bound holds on the entire half-line.",
        "",
        "## Remaining Continuous Theorem",
        "",
        "```text",
        exact["continuous_target"],
        exact["tent_transfer"],
        exact["conditional_full_ceiling"],
        "```",
        "",
        "Thus the continuous theorem supplies fewer than 201 inverse squares,",
        "the full-kernel transfer supplies fewer than 100, and the total 301",
        "fits strictly inside the endpoint target 320.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.md",
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
        "wrote order-six first/full curvature bridge: "
        f"{summary['rows']} rows, "
        f"{summary['full_kernel_transfer_theorems']} full transfer, "
        f"{summary['open_continuous_targets']} open continuous target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
