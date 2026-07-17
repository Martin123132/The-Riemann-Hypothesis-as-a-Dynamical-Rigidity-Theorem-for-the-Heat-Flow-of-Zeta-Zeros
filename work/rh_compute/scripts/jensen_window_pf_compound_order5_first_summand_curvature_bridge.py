#!/usr/bin/env python3
"""Reduce the order-five endpoint tail to one first-summand curvature bound."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order4_first_summand_curvature_bridge as order4_bridge  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_first_summand_curvature_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md"
)
TAIL_REDUCTION_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_m100_tail_curvature_reduction.json"
)
ORDER4_ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_m100_entry_certificate.json"
)
ORDER4_BRIDGE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_first_summand_curvature_bridge.json"
)
DOMINANCE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_dominance_certificate.json"
)
TAIL_FIRST_K = 321
FIRST_CONTINUOUS_T = TAIL_FIRST_K - 1
FULL_TRANSFER_CONSTANT = 37
FIRST_DISCRETE_CONSTANT = 63
FIRST_CONTINUOUS_CONSTANT = 60
SCOUT_T = (321, 400, 1000, 5000, 100000)
SCOUT_DPS = 50


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
    tail = load_json(TAIL_REDUCTION_SOURCE)
    order4_entry = load_json(ORDER4_ENTRY_SOURCE)
    order4 = load_json(ORDER4_BRIDGE_SOURCE)
    dominance = load_json(DOMINANCE_SOURCE)

    if tail.get("exact", {}).get("sufficient_ceiling") != (
        "C_n<=100/k^2 for every k=n+4>=321"
    ):
        raise RuntimeError("order-five curvature source contract changed")
    if order4_entry.get("exact", {}).get("global_first_summand_curvature") != (
        "K_1(t)<=7/(2*t^2), t>=319"
    ):
        raise RuntimeError("order-four first-summand curvature theorem changed")
    if order4_entry.get("tail_arithmetic", {}).get("total_penalty_bound") != (
        "P_n<=4/k^2"
    ):
        raise RuntimeError("order-four complete penalty theorem changed")
    if order4.get("exact", {}).get("j_floor_lower") != (
        "J_1(t)>=2/(2*t+3)-2/(2*t-1)+(t-3)/(6*t^2)>=1/(7*t), t>=319"
    ):
        raise RuntimeError("order-four first-summand J floor changed")
    dominance_bound = dominance.get("diagnostics", {}).get(
        "full_tail_relative_bound"
    )
    if dominance_bound != (
        "0<=delta_k=(M_k-M_k^(1))/M_k^(1)<=2/k^6 for every integer k>=300"
    ):
        raise RuntimeError("first-summand dominance contract changed")
    return {
        "tail_ceiling": tail["exact"]["sufficient_ceiling"],
        "order4_first_curvature": order4_entry["exact"][
            "global_first_summand_curvature"
        ],
        "order4_first_penalty": order4_entry["tail_arithmetic"][
            "first_summand_bound"
        ],
        "order4_full_penalty": order4_entry["tail_arithmetic"][
            "total_penalty_bound"
        ],
        "first_J_floor": order4["exact"]["j_floor_lower"],
        "integer_B_bounds": order4["exact"]["integer_B_bounds"],
        "moment_dominance": dominance_bound,
    }


def moment_error(index: sp.Expr) -> sp.Expr:
    return 2 * (
        1 / (index - 1) ** 6
        + 2 / index**6
        + 1 / (index + 1) ** 6
    )


def log_defect_error(index: sp.Expr) -> sp.Expr:
    return 4 * index * moment_error(index)


def first_gap_coordinate_error(index: sp.Expr) -> sp.Expr:
    return (
        2 * moment_error(index)
        + log_defect_error(index - 1)
        + 2 * log_defect_error(index)
        + log_defect_error(index + 1)
    )


def log_first_gap_error(index: sp.Expr) -> sp.Expr:
    return (
        2 * log_defect_error(index)
        + 8 * index * first_gap_coordinate_error(index)
    )


def second_gap_coordinate_error(index: sp.Expr) -> sp.Expr:
    return (
        3 * moment_error(index)
        + log_first_gap_error(index - 1)
        + 2 * log_first_gap_error(index)
        + log_first_gap_error(index + 1)
    )


def log_second_gap_over_defect_error(index: sp.Expr) -> sp.Expr:
    return (
        2 * log_first_gap_error(index)
        + sp.Rational(5, 7) * index * second_gap_coordinate_error(index)
        + log_defect_error(index)
    )


def shifted_positive_polynomial(
    expression: sp.Expr, variable: sp.Symbol, start: int
) -> dict:
    numerator, denominator = sp.fraction(sp.factor(expression))
    shift = sp.symbols("m", integer=True, nonnegative=True)
    polynomial = sp.Poly(sp.expand(numerator.subs(variable, start + shift)), shift)
    coefficients = polynomial.all_coeffs()
    if any(coefficient <= 0 for coefficient in coefficients):
        raise RuntimeError(
            f"shifted positivity failed at {variable}>={start}: {polynomial.as_expr()}"
        )
    return {
        "start": start,
        "degree": polynomial.degree(),
        "coefficient_count": len(coefficients),
        "leading_coefficient": str(coefficients[0]),
        "constant_coefficient": str(coefficients[-1]),
        "minimum_coefficient": str(min(coefficients)),
        "denominator_factorization": str(sp.factor(denominator)),
        "shifted_polynomial": str(polynomial.as_expr()),
        "shifted_coefficients": [str(value) for value in coefficients],
    }


def exact_diagnostics() -> dict:
    j, k = sp.symbols("j k", integer=True, positive=True)

    gap_floor = shifted_positive_polynomial(
        1 / (56 * j) - first_gap_coordinate_error(j), j, 319
    )
    first_R_floor = shifted_positive_polynomial(
        3 / (2 * j + 1)
        - sp.Rational(18, 5) / j**2
        - sp.Rational(7, 5) / j,
        j,
        320,
    )
    full_R_floor = shifted_positive_polynomial(
        sp.Rational(753, 250) / (2 * j + 1)
        - 4 / j**2
        - sp.Rational(7, 5) / j,
        j,
        320,
    )

    total_error = sp.factor(
        log_second_gap_over_defect_error(k - 1)
        + 2 * log_second_gap_over_defect_error(k)
        + log_second_gap_over_defect_error(k + 1)
    )
    transfer = shifted_positive_polynomial(
        FULL_TRANSFER_CONSTANT / k**2 - total_error, k, TAIL_FIRST_K
    )
    endpoint_error = sp.factor(total_error.subs(k, TAIL_FIRST_K))
    endpoint_scaled = sp.factor(endpoint_error * TAIL_FIRST_K**2)

    t = sp.symbols("t", positive=True)
    phi = 1 / (sp.exp(t) - 1)
    chi = sp.exp(t) / (sp.exp(t) - 1) ** 2
    derivative_check = sp.simplify(sp.diff(sp.log(1 - sp.exp(-t)), t) - phi)
    second_check = sp.simplify(
        sp.diff(sp.log(1 - sp.exp(-t)), t, 2) + chi
    )
    if derivative_check != 0 or second_check != 0:
        raise RuntimeError("stable logarithm derivative identity failed")

    return {
        "continuous_coordinates": (
            "x=exp(-B), d=1-x, ell=log(d), "
            "g=d^2-x^2*d(t-1)*d(t+1), h=log(g)"
        ),
        "second_gap_coordinate": (
            "R(t)=3*B(t)-h(t-1)+2*h(t)-h(t+1)"
        ),
        "second_gap_factorization": (
            "f(t)=g(t)^2-x(t)^3*g(t-1)*g(t+1)="
            "g(t)^2*(1-exp(-R(t)))"
        ),
        "order5_log_coordinate": (
            "q(t)=log(f(t)/d(t))="
            "2*h(t)-ell(t)+log(1-exp(-R(t)))"
        ),
        "order5_curvature": (
            "q''=2*h''-ell''+phi(R)*R''-chi(R)*(R')^2; "
            "phi(z)=1/(exp(z)-1), chi(z)=exp(z)/(exp(z)-1)^2"
        ),
        "discrete_identity": (
            "C_n=Delta^2 q(k), k=n+4; "
            "C_n^(1)=integral_[-1,1](1-|s|)*q_1''(k+s) ds"
        ),
        "moment_error": (
            "a_j=2*((j-1)^(-6)+2*j^(-6)+(j+1)^(-6)); "
            "|B_j-B_j^(1)|<=a_j"
        ),
        "log_defect_error": "L_j=4*j*a_j; |ell_j-ell_j^(1)|<=L_j",
        "first_gap_error": (
            "U_j=2*a_j+L_(j-1)+2*L_j+L_(j+1); "
            "|J_j-J_j^(1)|<=U_j<=1/(56*j), j>=319"
        ),
        "first_gap_floor": (
            "J_j^(1)>=1/(7*j) and U_j<=1/(56*j) imply "
            "min(J_j,J_j^(1))>=1/(8*j)"
        ),
        "log_first_gap_error": (
            "V_j=2*L_j+8*j*U_j; |h_j-h_j^(1)|<=V_j"
        ),
        "second_gap_error": (
            "W_j=3*a_j+V_(j-1)+2*V_j+V_(j+1); "
            "|R_j-R_j^(1)|<=W_j"
        ),
        "second_gap_floor": (
            "min(R_j,R_j^(1))>=7/(5*j), j>=320"
        ),
        "single_q_error": (
            "E_j=2*V_j+(5*j/7)*W_j+L_j; "
            "|q_j-q_j^(1)|<=E_j"
        ),
        "full_transfer": (
            "|C_n-C_n^(1)|<=E_(k-1)+2*E_k+E_(k+1)"
            "<=37/k^2, k=n+4>=321"
        ),
        "continuous_target": (
            "q_1''(t)<=60/t^2 for every real t>=320"
        ),
        "tent_transfer": (
            "q_1''(t)<=60/t^2 => C_n^(1)<=60*[-log(1-1/k^2)]"
            "<=60/(k^2-1)<63/k^2, k>=321"
        ),
        "conditional_full_ceiling": (
            "q_1''(t)<=60/t^2 on t>=320 => C_n<=100/k^2 on k>=321"
        ),
        "gap_floor_polynomial": gap_floor,
        "first_R_floor_polynomial": first_R_floor,
        "full_R_floor_polynomial": full_R_floor,
        "transfer_polynomial": transfer,
        "endpoint_transfer_error": str(endpoint_error),
        "endpoint_scaled_transfer": str(endpoint_scaled),
        "endpoint_scaled_transfer_decimal": str(sp.N(endpoint_scaled, 30)),
        "endpoint_scaled_reserve_below_37": str(
            sp.N(FULL_TRANSFER_CONSTANT - endpoint_scaled, 30)
        ),
        "stable_derivative_residuals": [str(derivative_check), str(second_check)],
    }


def sci(value: mp.mpf, digits: int = 28) -> str:
    return mp.nstr(value, n=digits, min_fixed=-6, max_fixed=6)


def nested_scout_row(t_value: int) -> dict:
    t = mp.mpf(t_value)
    one = (mp.mpf(1), mp.mpf(0), mp.mpf(0))
    H = {shift: order4_bridge.h_jet(t + shift) for shift in range(-3, 4)}
    B = {
        shift: order4_bridge.jet_add(
            order4_bridge.jet_add(
                H[shift + 1], order4_bridge.jet_scale(H[shift], -2)
            ),
            H[shift - 1],
        )
        for shift in range(-2, 3)
    }
    x = {
        shift: order4_bridge.jet_exp(order4_bridge.jet_scale(B[shift], -1))
        for shift in range(-2, 3)
    }
    d = {
        shift: order4_bridge.jet_sub(one, x[shift]) for shift in range(-2, 3)
    }
    ell = {shift: order4_bridge.jet_log(d[shift]) for shift in range(-2, 3)}
    J: dict[int, tuple] = {}
    h: dict[int, tuple] = {}
    for shift in range(-1, 2):
        centered_ell = order4_bridge.jet_add(
            order4_bridge.jet_add(ell[shift - 1], ell[shift + 1]),
            order4_bridge.jet_scale(ell[shift], -2),
        )
        J[shift] = order4_bridge.jet_sub(
            order4_bridge.jet_scale(B[shift], 2), centered_ell
        )
        stable = order4_bridge.jet_log(
            order4_bridge.jet_sub(
                one,
                order4_bridge.jet_exp(order4_bridge.jet_scale(J[shift], -1)),
            )
        )
        h[shift] = order4_bridge.jet_add(
            order4_bridge.jet_scale(ell[shift], 2), stable
        )

    centered_h = order4_bridge.jet_add(
        order4_bridge.jet_add(h[-1], h[1]), order4_bridge.jet_scale(h[0], -2)
    )
    R = order4_bridge.jet_sub(order4_bridge.jet_scale(B[0], 3), centered_h)
    stable_R = order4_bridge.jet_log(
        order4_bridge.jet_sub(
            one, order4_bridge.jet_exp(order4_bridge.jet_scale(R, -1))
        )
    )
    q = order4_bridge.jet_add(
        order4_bridge.jet_sub(order4_bridge.jet_scale(h[0], 2), ell[0]),
        stable_R,
    )
    scaled = t**2 * q[2]
    return {
        "t": t_value,
        "saddle": sci(order4_bridge.saddle_point(t)),
        "t_times_J": sci(t * J[0][0]),
        "t_times_R": sci(t * R[0]),
        "t2_q_second": sci(scaled),
        "margin_below_60": sci(FIRST_CONTINUOUS_CONSTANT - scaled),
        "below_target": bool(scaled < FIRST_CONTINUOUS_CONSTANT),
        "positive_stable_coordinates": bool(J[0][0] > 0 and R[0] > 0),
        "proof_boundary": (
            "Finite-upper high-precision mpmath saddle quadrature; not an "
            "interval enclosure or a uniform theorem."
        ),
    }


def finite_scout() -> dict:
    mp.mp.dps = SCOUT_DPS
    rows = [nested_scout_row(value) for value in SCOUT_T]
    if not all(
        row["below_target"] and row["positive_stable_coordinates"] for row in rows
    ):
        raise RuntimeError("order-five first-summand curvature scout failed")
    scaled = [mp.mpf(row["t2_q_second"]) for row in rows]
    return {
        "mpmath_dps": SCOUT_DPS,
        "sample_t": list(SCOUT_T),
        "rows": rows,
        "observed_scaled_range": [sci(min(scaled)), sci(max(scaled))],
        "formal_limit_scout": "t^2*q_1''(t) appears to increase to 6 from below",
        "proof_boundary": "Exploratory finite sampling only.",
    }


def build_artifact(*, include_scout: bool = True) -> dict:
    sources = validate_sources()
    exact = exact_diagnostics()
    scout = finite_scout() if include_scout else None
    rows = [
        BridgeRow(
            id="co5fscb_01_stable_second_gap",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The order-four stable margin has a second cancellation-preserving exponential coordinate.",
            formula=exact["second_gap_factorization"] + "; " + exact["second_gap_coordinate"],
            proof_boundary="Exact first-summand and full-kernel algebra.",
        ),
        BridgeRow(
            id="co5fscb_02_continuous_curvature",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The first-summand order-five logarithmic coordinate has an exact stable curvature formula.",
            formula=exact["order5_log_coordinate"] + "; " + exact["order5_curvature"],
            proof_boundary="Exact differentiation; no curvature inequality is inferred.",
        ),
        BridgeRow(
            id="co5fscb_03_gap_floors",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="Existing order-three and order-four theorems keep both stable logarithms uniformly away from zero on the complete tail.",
            formula=exact["first_gap_floor"] + "; " + exact["second_gap_floor"],
            proof_boundary="Composes proved first/full lower-layer bounds with exact rational comparisons.",
        ),
        BridgeRow(
            id="co5fscb_04_nested_error_chain",
            role="exact_perturbation_lemma",
            readiness="ready_to_apply",
            claim="The first-summand moment error propagates through both stable logarithms without subtracting tiny raw determinants.",
            formula=(
                exact["moment_error"]
                + "; "
                + exact["log_defect_error"]
                + "; "
                + exact["log_first_gap_error"]
                + "; "
                + exact["single_q_error"]
            ),
            proof_boundary="Exact Lipschitz chain using only proved positive floors.",
        ),
        BridgeRow(
            id="co5fscb_05_full_transfer",
            role="exact_perturbation_theorem",
            readiness="ready_to_apply",
            claim="The complete theta kernel changes the order-five curvature by at most thirty-seven inverse squares.",
            formula=exact["full_transfer"],
            proof_boundary="Coefficient-positive rational theorem for every k>=321.",
            diagnostics={
                "polynomial_degree": exact["transfer_polynomial"]["degree"],
                "coefficient_count": exact["transfer_polynomial"][
                    "coefficient_count"
                ],
                "endpoint_scaled_transfer": exact["endpoint_scaled_transfer_decimal"],
                "endpoint_reserve": exact["endpoint_scaled_reserve_below_37"],
            },
        ),
        BridgeRow(
            id="co5fscb_06_tent_transfer",
            role="conditional_theorem",
            readiness="conditional_on_open_input",
            claim="A loose continuous first-summand curvature ceiling supplies the remaining sixty-three inverse squares.",
            formula=exact["tent_transfer"],
            proof_boundary="Conditional on the displayed real-parameter curvature ceiling.",
        ),
        BridgeRow(
            id="co5fscb_07_conditional_completion",
            role="conditional_theorem",
            readiness="conditional_on_open_input",
            claim="The continuous first-summand theorem would prove the original complete-kernel scalar ceiling.",
            formula=exact["conditional_full_ceiling"],
            proof_boundary="Conditional; not order-five entry, PF-infinity, RH, or Lambda<=0.",
        ),
        BridgeRow(
            id="co5fscb_08_open_continuous_target",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove one continuous first-summand nested stable-curvature estimate.",
            formula=exact["continuous_target"],
            proof_boundary="Open analytic theorem; finite scout rows are not promoted.",
            diagnostics={
                "scout_scaled_range": None
                if scout is None
                else scout["observed_scaled_range"],
                "formal_limit_scout": None
                if scout is None
                else scout["formal_limit_scout"],
            },
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order5_first_summand_curvature_bridge",
        "date": "2026-07-13",
        "status": (
            "exact order-five first/full curvature transfer with one open "
            "continuous first-summand ceiling"
        ),
        "proof_boundary": (
            "This artifact proves the stable factorization, positive floors, "
            "and complete-kernel perturbation budget. It does not prove the "
            "continuous 60/t^2 ceiling, order-five entry, PF-infinity, RH, or "
            "Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.md",
            "outputs/jensen_window_pf_compound_order4_m100_entry_certificate.md",
            "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order5_first_summand_curvature_bridge.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order5_first_summand_curvature_bridge.py"
        ),
        "parameters": {
            "tail_first_k": TAIL_FIRST_K,
            "continuous_first_t": FIRST_CONTINUOUS_T,
            "full_transfer_constant": FULL_TRANSFER_CONSTANT,
            "first_discrete_constant": FIRST_DISCRETE_CONSTANT,
            "first_continuous_constant": FIRST_CONTINUOUS_CONSTANT,
        },
        "source_contract": sources,
        "exact": exact,
        "scout": scout,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": sum(row.readiness == "ready_to_apply" for row in rows),
            "conditional_rows": sum(
                row.readiness == "conditional_on_open_input" for row in rows
            ),
            "open_rows": sum(
                row.readiness == "not_ready_to_apply" for row in rows
            ),
            "exact_factorizations": 1,
            "positive_floor_theorems": 2,
            "full_kernel_transfer_theorems": 1,
            "scout_rows": 0 if scout is None else len(scout["rows"]),
            "open_continuous_targets": 1,
        },
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    scout = artifact["scout"]
    transfer = exact["transfer_polynomial"]
    lines = [
        "# Jensen-Window PF Compound Order-Five First-Summand Curvature Bridge",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact order-five first/full curvature transfer with one open",
        "continuous first-summand ceiling. This is not a proof of order-five",
        "entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order5_first_summand_curvature_bridge.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order5_first_summand_curvature_bridge.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_first_summand_curvature_bridge.py",
        "```",
        "",
        "## Second Stable Coordinate",
        "",
        "For the continuous first-summand family retain `B,d,ell,g,h=log(g)`",
        "from the completed order-four analysis and put",
        "",
        "```text",
        exact["second_gap_coordinate"],
        exact["second_gap_factorization"],
        exact["order5_log_coordinate"],
        "```",
        "",
        "Then `F_n=f(n+3)` and, with `k=n+4`,",
        "",
        "```text",
        exact["discrete_identity"],
        "```",
        "",
        "Differentiation gives the cancellation-preserving formula",
        "",
        "```text",
        exact["order5_curvature"],
        "```",
        "",
        "The negative square term is retained; no sign is discarded in the",
        "identity itself.",
        "",
        "## Positive Floors",
        "",
        "The proved first-summand `J` floor and the moment perturbation give",
        "",
        "```text",
        exact["first_gap_error"],
        exact["first_gap_floor"],
        "```",
        "",
        "The shifted numerator proving `U_j<=1/(56*j)` has degree",
        f"`{exact['gap_floor_polynomial']['degree']}` and all",
        f"`{exact['gap_floor_polynomial']['coefficient_count']}` coefficients",
        "positive after `j=319+m`.",
        "",
        "At the second stable layer the completed first and full order-four",
        "penalty bounds separately imply",
        "",
        "```text",
        exact["second_gap_floor"],
        "first numerator after j=320+m: m^2+597*m+88622>0",
        "full numerator after j=320+m: 53*m^2+31570*m+4674200>0",
        "```",
        "",
        "Thus neither logarithmic Lipschitz step divides by an uncertified",
        "near-zero raw determinant.",
        "",
        "## Full-Kernel Transfer",
        "",
        "Starting from the proved `0<=delta_j<=2/j^6`, define the explicit",
        "error chain",
        "",
        "```text",
        exact["moment_error"],
        exact["log_defect_error"],
        exact["log_first_gap_error"],
        exact["second_gap_error"],
        exact["single_q_error"],
        "```",
        "",
        "The centered difference therefore satisfies the unconditional theorem",
        "",
        "```text",
        exact["full_transfer"],
        "```",
        "",
        f"After `k=321+m`, its cleared reserve numerator has degree `{transfer['degree']}`,",
        f"all `{transfer['coefficient_count']}` coefficients are positive, its leading",
        f"coefficient is `{transfer['leading_coefficient']}`, and its constant",
        f"coefficient is `{transfer['constant_coefficient']}`.",
        "At the splice, the scaled error and reserve are",
        "",
        "```text",
        f"321^2*error_321={exact['endpoint_scaled_transfer_decimal']}",
        f"37-321^2*error_321={exact['endpoint_scaled_reserve_below_37']}",
        "```",
        "",
        "## Remaining Continuous Target",
        "",
        "The exact tent identity shows that it is now sufficient to prove",
        "",
        "```text",
        exact["continuous_target"],
        "```",
        "",
        "Indeed,",
        "",
        "```text",
        exact["tent_transfer"],
        exact["conditional_full_ceiling"],
        "```",
        "",
        "The constants split the original budget exactly as `37+63=100`.",
        "The displayed continuous theorem is open.",
        "",
        "## High-Precision Scout",
        "",
        "These finite-upper mpmath saddle quadratures are diagnostics, not",
        "interval enclosures or a uniform theorem.",
        "",
        "| t | t J(t) | t R(t) | t^2 q''(t) | margin below 60 |",
        "|---:|---:|---:|---:|---:|",
    ]
    if scout is not None:
        for row in scout["rows"]:
            lines.append(
                f"| `{row['t']}` | `{row['t_times_J']}` | `{row['t_times_R']}` | "
                f"`{row['t2_q_second']}` | `{row['margin_below_60']}` |"
            )
        lines.extend(
            [
                "",
                f"Observed scaled range: `{scout['observed_scaled_range']}`.",
                "The data are consistent with a limit of `6` from below, leaving",
                "a factor-ten reserve under the proposed analytic constant `60`.",
            ]
        )
    lines.extend(
        [
            "",
            "```text",
            "outputs/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.md",
            "outputs/jensen_window_pf_compound_order4_m100_entry_certificate.md",
            "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
            "outputs/signed_hankel_jensen_bridge_target.md",
            "outputs/formal_core.md",
            "```",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument(
        "--skip-scout",
        action="store_true",
        help="Rebuild the exact artifact without the non-rigorous scout rows.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact(include_scout=not args.skip_scout)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-five first-summand curvature bridge: "
        f"{summary['rows']} rows, "
        f"{summary['full_kernel_transfer_theorems']} full-kernel transfer, "
        f"{summary['scout_rows']} scout rows, "
        f"{summary['open_continuous_targets']} open continuous target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
