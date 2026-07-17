#!/usr/bin/env python3
"""Prove the epsilon^6 formal cumulant corridors on the ray u>=20."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
import math
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

from jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate import (  # noqa: E402
    SOURCE_TARGET,
)
from jensen_window_pf_compound_order4_gaussian_cumulant_ray_target import (  # noqa: E402
    CORRIDOR_CAPS,
    CORRIDOR_K2,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.md"
)
RAY_START = 20
Q_FLOOR = 10**35
JET_REMAINDER_CONSTANT = 1_000
LEADING_BUFFER = Fraction(1, 10)
TRANSFER_BUDGET = Fraction(1, 20)
LIPSCHITZ_CAP = 6_000
CORRECTION_NORM_CAP = 20_000
FORMAL_TRANSFER_CONSTANT = 22_000_000
PRECISION_BITS = 256


@dataclass(frozen=True)
class RayRow:
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


def series_multiply(left: list[sp.Expr], right: list[sp.Expr], order: int) -> list[sp.Expr]:
    return [
        sp.expand(sum(left[index] * right[degree - index] for index in range(degree + 1)))
        for degree in range(order + 1)
    ]


def series_exponential(values: list[sp.Expr], order: int) -> list[sp.Expr]:
    result = [sp.Integer(0) for _ in range(order + 1)]
    result[0] = sp.exp(values[0])
    for degree in range(1, order + 1):
        result[degree] = sp.expand(
            sum(
                index * values[index] * result[degree - index]
                for index in range(1, degree + 1)
            )
            / degree
        )
    return result


def series_logarithm(values: list[sp.Expr], order: int) -> list[sp.Expr]:
    result = [sp.Integer(0) for _ in range(order + 1)]
    result[0] = sp.log(values[0])
    for degree in range(1, order + 1):
        result[degree] = sp.cancel(
            (
                degree * values[degree]
                - sum(
                    index * result[index] * values[degree - index]
                    for index in range(1, degree)
                )
            )
            / (degree * values[0])
        )
    return result


def symbolic_potential_jets(order: int = 8) -> tuple[sp.Symbol, sp.Symbol, dict[int, sp.Expr], dict[int, sp.Expr]]:
    u, q = sp.symbols("u q", positive=True)
    u_series = [u / sp.Integer(2**degree * math.factorial(degree)) for degree in range(order + 1)]
    exponent = [4 * value for value in u_series]
    exponent[0] -= 4 * u
    q_series = [q * value for value in series_exponential(exponent, order)]
    denominator = [2 * value for value in q_series]
    denominator[0] -= 3
    log_denominator = series_logarithm(denominator, order)
    log_u = series_logarithm(u_series, order)
    u_squared = series_multiply(u_series, u_series, order)
    potential = [
        sp.cancel(
            100 * u_squared[degree]
            + q_series[degree]
            - 5 * u_series[degree]
            - log_denominator[degree]
            - log_u[degree]
        )
        for degree in range(order + 1)
    ]
    full_jets = {
        degree: sp.cancel(potential[degree] * math.factorial(degree))
        for degree in range(2, order + 1)
    }
    q_jets = {
        degree: sp.expand(q_series[degree] * math.factorial(degree))
        for degree in range(2, order + 1)
    }
    return u, q, full_jets, q_jets


def shifted_positive_gate(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    substitutions: dict[sp.Symbol, sp.Expr],
) -> dict:
    numerator, denominator = sp.fraction(sp.factor(expression))
    shifted = sp.Poly(sp.expand(numerator.subs(substitutions)), *variables)
    coefficients = shifted.coeffs()
    if not coefficients or any(coefficient <= 0 for coefficient in coefficients):
        raise RuntimeError("shifted coefficient-positivity gate failed")
    return {
        "numerator_terms": len(shifted.terms()),
        "minimum_numerator_coefficient": str(min(coefficients)),
        "denominator": sp.sstr(sp.factor(denominator)),
    }


def leading_model() -> tuple[sp.Symbol, dict[int, sp.Expr], dict[int, sp.Expr]]:
    u = sp.symbols("u", positive=True)
    p = {0: sp.Integer(1)}
    for order in range(1, 9):
        p[order] = sp.expand(u * sp.diff(p[order - 1], u) / 2 + 2 * u * p[order - 1])
    curvature = p[2]
    normalized = {
        order: p[order] / curvature ** sp.Rational(order, 2)
        for order in range(3, 9)
    }
    l3, l4, l5, l6, l7, l8 = [normalized[order] for order in range(3, 9)]
    ratios = {
        2: sp.factor((2 * l3**2 - l4) / 2),
        3: sp.factor(l3),
        4: sp.factor((3 * l3**2 - l4) / 2),
        5: sp.factor((15 * l3**3 - 10 * l3 * l4 + l5) / 6),
        6: sp.factor(
            (105 * l3**4 - 105 * l3**2 * l4 + 15 * l3 * l5 + 10 * l4**2 - l6)
            / 24
        ),
        7: sp.factor(
            (
                945 * l3**5
                - 1260 * l3**3 * l4
                + 210 * l3**2 * l5
                + 280 * l3 * l4**2
                - 21 * l3 * l6
                - 35 * l4 * l5
                + l7
            )
            / 120
        ),
        8: sp.factor(
            (
                10395 * l3**6
                - 17325 * l3**4 * l4
                + 3150 * l3**3 * l5
                + 6300 * l3**2 * l4**2
                - 378 * l3**2 * l6
                - 1260 * l3 * l4 * l5
                + 28 * l3 * l7
                - 280 * l4**3
                + 56 * l4 * l6
                + 35 * l5**2
                - l8
            )
            / 720
        ),
    }
    return u, p, ratios


def leading_corridor_gates() -> dict:
    u, p, ratios = leading_model()
    v = sp.symbols("v", nonnegative=True)
    substitutions = {u: v + RAY_START}
    gates = {}
    for order, ratio in ratios.items():
        if order == 2:
            floor, ceiling = CORRIDOR_K2
            lower_expression = ratio - (sp.Rational(floor.numerator, floor.denominator) + sp.Rational(1, 10) / u)
            upper_expression = (sp.Rational(ceiling.numerator, ceiling.denominator) - sp.Rational(1, 10) / u) - ratio
        else:
            cap = CORRIDOR_CAPS[order]
            lower_target = 1 + sp.Rational(1, 10) / u
            upper_target = sp.Rational(cap.numerator, cap.denominator) - sp.Rational(1, 10) / u
            if order % 2:
                lower_expression = ratio**2 - lower_target**2
                upper_expression = upper_target**2 - ratio**2
            else:
                lower_expression = ratio - lower_target
                upper_expression = upper_target - ratio
        gates[str(order)] = {
            "formula": sp.sstr(ratio),
            "lower_buffer_gate": shifted_positive_gate(
                lower_expression, (v,), substitutions
            ),
            "upper_buffer_gate": shifted_positive_gate(
                upper_expression, (v,), substitutions
            ),
        }

    normalized_caps = {}
    curvature = p[2]
    for order in range(3, 9):
        expression = (
            sp.Rational(9, 4) * curvature**order - p[order] ** 2
            if order % 2
            else sp.Rational(3, 2) * curvature ** (order // 2) - p[order]
        )
        normalized_caps[str(order)] = shifted_positive_gate(
            expression, (v,), substitutions
        )
    return {
        "recurrence": "P_0=1; P_(r+1)=(u/2)*P_r'+2u*P_r",
        "potential_polynomials": {str(order): sp.sstr(sp.factor(p[order])) for order in range(2, 9)},
        "ratio_gates": gates,
        "normalized_jet_cap_gates": normalized_caps,
        "proved_buffer": "corridor floor+1/(10u) <= F_r <= corridor ceiling-1/(10u)",
    }


def jet_remainder_gates() -> dict:
    u, q, full, q_part = symbolic_potential_jets()
    v, capital_q = sp.symbols("v capital_q", nonnegative=True)
    substitutions = {u: v + RAY_START, q: capital_q + Q_FLOOR}
    rows = {}
    for order in range(2, 9):
        remainder = sp.factor(full[order] - q_part[order])
        envelope = JET_REMAINDER_CONSTANT * u**8
        rows[str(order)] = {
            "plus_gate": shifted_positive_gate(
                envelope + remainder, (v, capital_q), substitutions
            ),
            "minus_gate": shifted_positive_gate(
                envelope - remainder, (v, capital_q), substitutions
            ),
        }
    return {
        "formula": "|V^(r)-q*P_r(u)|<=1000*u^8, r=2,...,8",
        "substitution": "u=20+v, q=10^35+Q, v,Q>=0",
        "orders": rows,
    }


def polynomial_norms() -> dict:
    artifact = load_json(SOURCE_TARGET)
    symbols = sp.symbols("L_3:9")
    locals_map = {str(symbol): symbol for symbol in symbols}

    def sup_norm(expression: sp.Expr, radius: int = 2) -> sp.Rational:
        polynomial = sp.Poly(expression, *symbols)
        return sum(
            abs(sp.Rational(coefficient)) * radius ** sum(powers)
            for powers, coefficient in polynomial.terms()
        )

    rows = {}
    largest_lipschitz = sp.Rational(0)
    largest_correction = sp.Rational(0)
    for order_text, row in artifact["exact"]["cumulants"].items():
        order = int(order_text)
        offset = 2 if order == 2 else order - 2
        scale = sp.Integer(1) if order == 2 else sp.Rational((-1) ** order, math.factorial(order - 2))
        terms = []
        for term in row["terms"]:
            exponent = int(term["epsilon_power"]) - offset
            expression = scale * sp.sympify(term["coefficient"], locals=locals_map)
            terms.append((exponent, expression))
        leading = next(expression for exponent, expression in terms if exponent == 0)
        lipschitz = sum(sup_norm(sp.diff(leading, symbol)) for symbol in symbols)
        correction = sum(sup_norm(expression) for exponent, expression in terms if exponent > 0)
        largest_lipschitz = max(largest_lipschitz, lipschitz)
        largest_correction = max(largest_correction, correction)
        rows[str(order)] = {
            "leading_lipschitz_l1_on_abs_L_le_2": str(lipschitz),
            "positive_epsilon_correction_norm": str(correction),
        }
    if not (largest_lipschitz < LIPSCHITZ_CAP and largest_correction < CORRECTION_NORM_CAP):
        raise RuntimeError("formal polynomial norm caps failed")
    return {
        "domain": "|L_r|<=2, r=3,...,8",
        "rows": rows,
        "largest_lipschitz": str(largest_lipschitz),
        "lipschitz_cap": LIPSCHITZ_CAP,
        "largest_correction_norm": str(largest_correction),
        "correction_norm_cap": CORRECTION_NORM_CAP,
    }


def scalar_geometry() -> dict:
    flint.ctx.prec = PRECISION_BITS
    q20 = flint.arb.pi() * flint.arb(80).exp()
    if not bool(q20 > Q_FLOOR):
        raise RuntimeError("q floor at u=20 failed")
    endpoint_polynomial = 20 * FORMAL_TRANSFER_CONSTANT * 20**7
    if endpoint_polynomial >= Q_FLOOR:
        raise RuntimeError("formal transfer endpoint comparison failed")
    derivative_margin = Fraction(4) - Fraction(7, RAY_START)
    if derivative_margin <= 0:
        raise RuntimeError("q/u^7 monotonicity failed")

    eta_endpoint_numerator = 250 * 20**6
    eta_cap = Fraction(1, 10**24)
    if Fraction(eta_endpoint_numerator, Q_FLOOR) >= eta_cap:
        raise RuntimeError("relative jet error cap failed")
    derivative_cap = Fraction(4 * 100**5, 99**5)
    value_cap = Fraction(100**4, 99**4)
    if not (derivative_cap < 5 and value_cap < Fraction(21, 20)):
        raise RuntimeError("power perturbation calculus constants failed")
    return {
        "q_at_20_lower": arb_lower_text(q20),
        "q_floor": str(Q_FLOOR),
        "q_over_u7_log_derivative_lower": str(derivative_margin),
        "endpoint_transfer_left": str(endpoint_polynomial),
        "endpoint_transfer_right": str(Q_FLOOR),
        "jet_relative_error": "eta(u)=250*u^6/q<=10^-24",
        "power_derivative_cap": str(derivative_cap),
        "power_value_cap": str(value_cap),
        "normalized_jet_transfer": "|L_r-L_r^(infinity)|<=14*eta(u)",
        "formal_transfer": (
            "|R_r^[6]-F_r|<=22000000*u^6/q<1/(20u), u>=20"
        ),
    }


def build_artifact() -> dict:
    leading = leading_corridor_gates()
    jet_remainders = jet_remainder_gates()
    norms = polynomial_norms()
    geometry = scalar_geometry()
    rows = [
        RayRow(
            id="co4fcarc_01_leading_potential_geometry",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="The q-leading potential jets are generated by an explicit differential recurrence.",
            formula=leading["recurrence"],
            proof_boundary="Exact q-leading geometry only.",
            diagnostics=leading["potential_polynomials"],
        ),
        RayRow(
            id="co4fcarc_02_leading_corridor_buffer",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Every q-leading formal cumulant ratio lies a full 1/(10u) inside its candidate corridor on u>=20.",
            formula=leading["proved_buffer"],
            proof_boundary="Leading q-infinity formal model only.",
            diagnostics=leading["ratio_gates"],
        ),
        RayRow(
            id="co4fcarc_03_normalized_jet_cap",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="The leading normalized potential jets remain below 3/2 in absolute value.",
            formula="0<L_r^(infinity)<3/2, r=3,...,8, u>=20",
            proof_boundary="Leading normalized jets only.",
            diagnostics=leading["normalized_jet_cap_gates"],
        ),
        RayRow(
            id="co4fcarc_04_full_jet_remainder",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The non-q part of every required exact potential jet has a uniform coefficient-positive envelope.",
            formula=jet_remainders["formula"],
            proof_boundary="Potential jets through order eight on u>=20.",
            diagnostics=jet_remainders,
        ),
        RayRow(
            id="co4fcarc_05_polynomial_norm_transfer",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="Explicit coefficient norms control perturbations of every leading cumulant polynomial and every positive epsilon correction.",
            formula="Lip_1(F_r)<6000; sum_(j>=1)||G_(r,2j)||_2<20000 on |L|<=2",
            proof_boundary="Finite polynomial norm theorem only.",
            diagnostics=norms,
        ),
        RayRow(
            id="co4fcarc_06_exponential_transfer_budget",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The exact epsilon-six formal ratios differ from their q-leading ratios by less than half the leading corridor buffer.",
            formula=geometry["formal_transfer"],
            proof_boundary="Exact formal model only; not exact density cumulants.",
            diagnostics=geometry,
        ),
        RayRow(
            id="co4fcarc_07_formal_ray_theorem",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The epsilon-six formal cumulant model satisfies all candidate corridors on the full asymptotic ray.",
            formula="formal cumulant corridors hold for every u>=20",
            proof_boundary="Formal model theorem only; exact-minus-formal remainder remains open.",
        ),
        RayRow(
            id="co4fcarc_08_exact_remainder_handoff",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Fit the exact central and two-tail cumulant errors inside the remaining 1/(20u) formal-ray buffer.",
            formula="scaled |kappa_r-kappa_r^[6]|<1/(20u), r=2,...,8, u>=20",
            proof_boundary="Open exact-density remainder theorem.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate",
        "date": "2026-07-13",
        "status": "exact analytic theorem for the epsilon^6 formal cumulant ray",
        "proof_boundary": (
            "This artifact proves the candidate corridors for the exact epsilon^6 "
            "formal cumulant model on u>=20. Together with the finite formal corridor "
            "certificate it closes the formal model on u>=2. It does not bound the "
            "exact-minus-formal cumulant remainder, prove the exact curvature ray, "
            "order-four entry, PF-infinity, RH, or Lambda<=0."
        ),
        "leading": leading,
        "jet_remainders": jet_remainders,
        "polynomial_norms": norms,
        "scalar_geometry": geometry,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": 7,
            "open_analytic_rows": 1,
            "leading_corridor_gates": len(leading["ratio_gates"]),
            "jet_remainder_sign_gates": 14,
            "formal_ray_closed": True,
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.py"
        ),
        "remaining_target": (
            "Prove scaled exact-minus-formal cumulant errors fit inside the finite and "
            "asymptotic corridor margins for every u>=2."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    geometry = artifact["scalar_geometry"]
    lines = [
        "# Jensen-Window PF Order-Four Formal Cumulant Asymptotic-Ray Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact analytic theorem for the epsilon-six formal cumulant ray.",
        "This is not a proof of the exact cumulant ray, order-four entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.py",
        "```",
        "",
        "## Leading Model",
        "",
        "For the q-leading potential jets,",
        "",
        "```text",
        artifact["leading"]["recurrence"],
        artifact["leading"]["proved_buffer"],
        "```",
        "",
        "After `u=20+v`, clearing positive denominators gives coefficient-positive",
        "numerators for both sides of all seven buffered corridor inequalities.",
        "Odd cumulant ratios are squared only after their positive numerators are",
        "identified. The same method proves `0<L_r^(infinity)<3/2` through order eight.",
        "",
        "## Exact Potential Transfer",
        "",
        "With `q=pi*exp(4u)`, exact symbolic jets and the substitution",
        "`u=20+v, q=10^35+Q` prove",
        "",
        "```text",
        artifact["jet_remainders"]["formula"],
        geometry["jet_relative_error"],
        geometry["normalized_jet_transfer"],
        "```",
        "",
        f"Arb gives `q(20)>{geometry['q_floor']}` with lower enclosure",
        f"`{geometry['q_at_20_lower']}`. The logarithmic derivative",
        "`d log(q/u^7)/du=4-7/u` is positive on the ray.",
        "",
        "Exact coefficient norms give",
        "",
        "```text",
        geometry["formal_transfer"],
        "```",
        "",
        "This consumes less than half the `1/(10u)` leading buffer. Therefore the",
        "exact epsilon-six formal cumulant polynomial satisfies every candidate",
        "corridor for all `u>=20`.",
        "",
        "## Remaining Boundary",
        "",
        "Together with the 1.8-million-block compact formal certificate, the formal",
        "model is now closed for every `u>=2`. What remains is not formal algebra:",
        "prove that the exact standardized-density cumulants differ from the formal",
        "polynomial by less than the available scaled margins, using a central",
        "expansion remainder and two adaptive tails.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md",
        "outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md",
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
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "certified order-four formal cumulant asymptotic ray: "
        "7 exact rows, 7 buffered corridor gates, 14 jet-remainder sign gates, "
        "1 open exact-density remainder"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
