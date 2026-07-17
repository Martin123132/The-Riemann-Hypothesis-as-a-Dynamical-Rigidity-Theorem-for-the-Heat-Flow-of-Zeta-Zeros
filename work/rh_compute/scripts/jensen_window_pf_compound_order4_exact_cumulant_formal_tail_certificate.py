#!/usr/bin/env python3
"""Certify the two formal Gaussian tails in the epsilon-ten disk contract."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
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
import sympy as sp  # noqa: E402

from jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate import (  # noqa: E402
    sha256,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_upper_text,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate import (  # noqa: E402
    tail_polynomial_ratio,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.md"
)
SOURCE_DISK_CONTRACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.json"
)
SOURCE_SECOND_FINITE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.json"
)
FINITE_Q_FLOOR = 9_000
RAY_Q_FLOOR = 10**35
RAY_START = 20
MAX_EPSILON_ORDER = 10
MAXIMUM_Y_DEGREE = 30
Y_FLOOR = 13
SLOPE_RATIO = Fraction(12, 13)
FORMAL_TAIL_CONSTANT = 2**32


@dataclass(frozen=True)
class TailRow:
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


def formal_exponential_polynomials() -> tuple[
    tuple[sp.Symbol, ...], sp.Symbol, list[sp.Expr]
]:
    y = sp.symbols("y")
    symbols = tuple(sp.symbols("L_3:13"))
    perturbation = [sp.Integer(0) for _ in range(MAX_EPSILON_ORDER + 1)]
    for order, symbol in zip(range(3, 13), symbols):
        perturbation[order - 2] = symbol * y**order / sp.factorial(order)
    exponential = [sp.Integer(0) for _ in range(MAX_EPSILON_ORDER + 1)]
    exponential[0] = sp.Integer(1)
    for degree in range(1, MAX_EPSILON_ORDER + 1):
        exponential[degree] = sp.expand(
            -sum(
                index * perturbation[index] * exponential[degree - index]
                for index in range(1, degree + 1)
            )
            / degree
        )
    if max(sp.Poly(value, y).degree() for value in exponential) != MAXIMUM_Y_DEGREE:
        raise RuntimeError("unexpected epsilon-ten formal weight degree")
    return symbols, y, exponential


def polynomial_norm(
    expression: sp.Expr,
    symbols: tuple[sp.Symbol, ...],
    y: sp.Symbol,
    caps: dict[int, Fraction],
) -> sp.Rational:
    total = sp.Rational(0)
    for powers, coefficient in sp.Poly(expression, *symbols, y).terms():
        term = abs(sp.Rational(coefficient))
        for order, power in zip(range(3, 13), powers[:-1]):
            cap = caps[order]
            term *= sp.Rational(cap.numerator, cap.denominator) ** power
        total += term
    return sp.factor(total)


def formal_weight_majorant() -> dict:
    finite = load_json(SOURCE_SECOND_FINITE)
    caps = {
        int(order): Fraction(value)
        for order, value in finite["normalized_jet_caps"].items()
    }
    symbols, y, exponential = formal_exponential_polynomials()
    norms = [polynomial_norm(value, symbols, y, caps) for value in exponential]
    total = sp.factor(sum(norms))
    if not total < sp.Rational(4, 3):
        raise RuntimeError("formal weight coefficient majorant failed")
    return {
        "normalized_jet_caps": {str(order): str(cap) for order, cap in caps.items()},
        "epsilon_degree_norms": {
            str(degree): str(value) for degree, value in enumerate(norms)
        },
        "total_coefficient_norm": str(total),
        "total_coefficient_cap": "4/3",
        "maximum_y_degree": MAXIMUM_Y_DEGREE,
        "epsilon_powers_discarded_upward": (
            "epsilon^n<=1 is used only for the formal tail upper bound"
        ),
    }


def adaptive_window_geometry() -> dict:
    flint.ctx.prec = 256
    q = flint.arb(FINITE_Q_FLOOR)
    s = (32 * q.log()).sqrt()
    y = 1 + s
    q_quarter = q.sqrt().sqrt()
    endpoint_margin = 2 * q_quarter - y
    derivative_margin = q_quarter / 2 - 16 / s
    if not bool(y > Y_FLOOR and endpoint_margin > 0 and derivative_margin > 0):
        raise RuntimeError("adaptive formal-tail window gate failed")
    ratios = {
        str(order): tail_polynomial_ratio(order, SLOPE_RATIO, Y_FLOOR)
        for order in range(MAXIMUM_Y_DEGREE + 1)
    }
    maximum_order = max(ratios, key=ratios.get)
    maximum_ratio = ratios[maximum_order]
    if maximum_ratio >= Fraction(3, 2):
        raise RuntimeError("formal-tail polynomial ratio gate failed")
    return {
        "window": "Y(q)=1+sqrt(32*log(q))",
        "q_floor": FINITE_Q_FLOOR,
        "Y_at_q_floor_lower": arb_lower_text(y),
        "Y_at_q_floor_upper": arb_upper_text(y),
        "Y_floor": Y_FLOOR,
        "Y_power_bound": "Y(q)<2*q^(1/4), q>=9000",
        "endpoint_margin_lower": arb_lower_text(endpoint_margin),
        "derivative_margin_lower": arb_lower_text(derivative_margin),
        "slope": "Y-1=sqrt(32*log(q))",
        "slope_ratio": str(SLOPE_RATIO),
        "tail_polynomial_ratios": {key: str(value) for key, value in ratios.items()},
        "maximum_tail_polynomial_ratio_order": int(maximum_order),
        "maximum_tail_polynomial_ratio": str(maximum_ratio),
        "completed_square": (
            "exp(-Y^2/2+Y)=exp(1/2)*q^-16"
        ),
    }


def scalar_tail_budgets() -> dict:
    finite_left = 500_000 * FORMAL_TAIL_CONSTANT
    finite_left_fourth = finite_left**4
    finite_right = FINITE_Q_FLOOR**23
    if finite_left_fourth >= finite_right:
        raise RuntimeError("finite formal-tail endpoint comparison failed")

    ray_left = 100_000 * RAY_START * FORMAL_TAIL_CONSTANT
    ray_left_fourth = ray_left**4
    ray_right = RAY_Q_FLOOR**23
    if ray_left_fourth >= ray_right:
        raise RuntimeError("ray formal-tail endpoint comparison failed")
    finite_monotonicity = Fraction(23, 4)
    ray_monotonicity = Fraction(23) - Fraction(4, RAY_START)
    if ray_monotonicity <= 0:
        raise RuntimeError("q^23/u^4 monotonicity failed")
    return {
        "common_one_sided_bound": (
            "formal tail <=2^32*q^(-35/4) for |z|<=1"
        ),
        "finite": {
            "target": "each formal tail <1/(500000*q^3)",
            "endpoint_left_fourth": str(finite_left_fourth),
            "endpoint_right": str(finite_right),
            "q_power_margin": str(finite_monotonicity),
        },
        "ray": {
            "target": "each formal tail <1/(100000*u*q^3)",
            "endpoint_left_fourth": str(ray_left_fourth),
            "endpoint_right": str(ray_right),
            "q23_over_u4_log_derivative_lower": str(ray_monotonicity),
        },
        "allocation": {
            "finite_partition_target": "1/(100000*q^3)",
            "finite_five_equal_parts": "1/(500000*q^3)",
            "ray_partition_target": "1/(20000*u*q^3)",
            "ray_each_formal_tail": "1/(100000*u*q^3)",
            "ray_each_remaining_exact_component": "1/(300000*u*q^3)",
        },
    }


def build_artifact() -> dict:
    disk = load_json(SOURCE_DISK_CONTRACT)
    if disk.get("summary", {}).get("partition_degrees") != 10:
        raise RuntimeError("complex-disk source is not closed")
    weight = formal_weight_majorant()
    window = adaptive_window_geometry()
    budgets = scalar_tail_budgets()
    rows = [
        TailRow(
            id="co4ecftc_01_formal_weight_polynomial",
            role="exact_formal_algebra",
            readiness="ready_to_apply",
            claim="The epsilon-ten formal density factor is a degree-thirty polynomial in y with exact graded coefficients.",
            formula="E^[10](y)=sum_(n=0)^10 epsilon^n E_n(y), deg_y E^[10]=30",
            proof_boundary="Exact formal weight only.",
        ),
        TailRow(
            id="co4ecftc_02_weight_coefficient_majorant",
            role="exact_inequality",
            readiness="ready_to_apply",
            claim="The global finite normalized-jet caps bound the total absolute formal-weight coefficient norm by 4/3.",
            formula="sum_(n,m)|[y^m]E_n|<4/3",
            proof_boundary="Formal coefficient majorant only.",
            diagnostics=weight,
        ),
        TailRow(
            id="co4ecftc_03_adaptive_window",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="The adaptive cutoff has a q^-16 completed-square factor and a controlled outward Gaussian slope.",
            formula=window["window"],
            proof_boundary="Formal Gaussian geometry only.",
            diagnostics=window,
        ),
        TailRow(
            id="co4ecftc_04_polynomial_hazard",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="Exact exponential-tail integration controls every polynomial degree through thirty by less than 3/2 times Y^(m-1).",
            formula="P_m(Y-1,Y)<(3/2)*Y^(m-1), m=0,...,30",
            proof_boundary="One-sided polynomial Gaussian tail calculus only.",
            diagnostics={
                "maximum_order": window["maximum_tail_polynomial_ratio_order"],
                "maximum_ratio": window["maximum_tail_polynomial_ratio"],
            },
        ),
        TailRow(
            id="co4ecftc_05_common_formal_tail_bound",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Both complex unit-disk formal tails satisfy one explicit q-power bound.",
            formula=budgets["common_one_sided_bound"],
            proof_boundary="Formal tails only; not exact-density tails.",
        ),
        TailRow(
            id="co4ecftc_06_finite_formal_tails",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Each formal tail fits one fifth of the finite exact partition residual budget.",
            formula=budgets["finite"]["target"],
            proof_boundary="Formal tails on 2<=u<=20 only.",
            diagnostics=budgets["finite"],
        ),
        TailRow(
            id="co4ecftc_07_ray_formal_tails",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Each formal tail fits its assigned asymptotic exact partition residual budget.",
            formula=budgets["ray"]["target"],
            proof_boundary="Formal tails on u>=20 only.",
            diagnostics=budgets["ray"],
        ),
        TailRow(
            id="co4ecftc_08_exact_residual_handoff",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove only the central exact-minus-formal residual and the two exact-density tails at their remaining allocations.",
            formula="finite each <1/(500000*q^3); ray each <1/(300000*u*q^3)",
            proof_boundary="Open exact-density central and two-tail theorem.",
            diagnostics=budgets["allocation"],
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate",
        "date": "2026-07-13",
        "status": "exact formal two-tail theorem with open exact-density components",
        "proof_boundary": (
            "This artifact proves both epsilon-ten formal Gaussian tails fit the "
            "complex-disk partition budgets. It does not prove the exact central "
            "residual or either exact-density tail, the exact cumulant corridors, "
            "curvature ray, order-four entry, PF-infinity, RH, or Lambda<=0."
        ),
        "weight_majorant": weight,
        "adaptive_window": window,
        "tail_budgets": budgets,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": 7,
            "open_analytic_rows": 1,
            "maximum_y_degree": MAXIMUM_Y_DEGREE,
            "tail_polynomial_orders": MAXIMUM_Y_DEGREE + 1,
            "formal_tails_closed": 2,
            "open_exact_components": 3,
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
            "outputs/formal_core.md",
        ],
        "source_hashes": {
            SOURCE_DISK_CONTRACT.relative_to(REPO_ROOT).as_posix(): sha256(SOURCE_DISK_CONTRACT),
            SOURCE_SECOND_FINITE.relative_to(REPO_ROOT).as_posix(): sha256(SOURCE_SECOND_FINITE),
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.py"
        ),
        "remaining_target": (
            "At Y=1+sqrt(32 log q), prove the central exact-minus-formal partition "
            "residual and both exact-density tails fit their finite and ray allocations."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Jensen-Window PF Order-Four Exact Cumulant Formal-Tail Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact formal two-tail theorem with open exact-density components.",
        "This is not a proof of the exact cumulant corridors, order-four entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.py",
        "```",
        "",
        "## Formal Weight",
        "",
        "The epsilon-ten formal density factor has degree thirty in `y`. The global",
        "normalized-jet caps give the exact coefficient majorant",
        "",
        "```text",
        f"sum_(n,m)|[y^m]E_n|={artifact['weight_majorant']['total_coefficient_norm']}<4/3.",
        "```",
        "",
        "Use the adaptive cutoff",
        "",
        "```text",
        artifact["adaptive_window"]["window"],
        artifact["adaptive_window"]["completed_square"],
        "Y(q)<2*q^(1/4).",
        "```",
        "",
        "Exact exponential-tail integration through polynomial degree thirty has",
        "maximum ratio below `3/2`. Including the complex factors on `|z|<=1` gives",
        "",
        "```text",
        artifact["tail_budgets"]["common_one_sided_bound"],
        artifact["tail_budgets"]["finite"]["target"],
        artifact["tail_budgets"]["ray"]["target"],
        "```",
        "",
        "Both formal tails are therefore closed inside the unit-disk partition",
        "contract.",
        "",
        "## Remaining Boundary",
        "",
        "Three exact-density components remain: the central exact-minus-formal",
        "residual, the exact left tail, and the exact right tail. Their sufficient",
        "allocations are",
        "",
        "```text",
        "finite: each <1/(500000*q^3), 2<=u<=20,",
        "ray:    each <1/(300000*u*q^3), u>=20.",
        "```",
        "",
        "These are open targets. No exact-density tail or cumulant corridor is",
        "promoted by this formal-tail theorem.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "certified order-four exact cumulant formal tails: "
        "8 rows, 7 exact rows, 31 polynomial orders, "
        "2 closed formal tails, 3 open exact components"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
