#!/usr/bin/env python3
"""Prove second-next cumulant coefficient bounds on the ray u>=20."""

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

from jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate import (  # noqa: E402
    shifted_positive_gate,
)
from jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate import (  # noqa: E402
    evaluate_polynomial_arb,
)
from jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate import (  # noqa: E402
    sha256,
)
from jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate import (  # noqa: E402
    COEFFICIENT_BOUNDS,
    SOURCE_SECOND_NEXT,
    compile_coefficients,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_upper_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.md"
)
SOURCE_FINITE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.json"
)
SOURCE_FIRST_ASYMPTOTIC = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate.json"
)
RAY_START = 20
Q_FLOOR = 10**35
LEADING_BUFFER = Fraction(1, 100)
TRANSFER_BUDGET = Fraction(1, 200)
JET_REMAINDER_CONSTANT = 100_000
LIPSCHITZ_CAP = 400_000_000
TRANSFER_CONSTANT = 8_000_000_000
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


def leading_model() -> tuple[sp.Symbol, dict[int, sp.Expr], dict[int, sp.Expr]]:
    source = load_json(SOURCE_SECOND_NEXT)
    symbols = sp.symbols("L_3:13")
    locals_map = {str(symbol): symbol for symbol in symbols}
    u = sp.symbols("u", positive=True)
    potential = {0: sp.Integer(1)}
    for order in range(1, 13):
        potential[order] = sp.expand(
            u * sp.diff(potential[order - 1], u) / 2
            + 2 * u * potential[order - 1]
        )
    curvature = potential[2]
    substitutions = {
        symbols[order - 3]: potential[order]
        / curvature ** sp.Rational(order, 2)
        for order in range(3, 13)
    }
    coefficients = {
        int(order_text): sp.factor(
            sp.sympify(row["scaled_coefficient"], locals=locals_map).subs(
                substitutions
            )
        )
        for order_text, row in source["coefficient_rows"].items()
    }
    return u, potential, coefficients


def odd_endpoint_signs(potential: dict[int, sp.Expr]) -> dict:
    flint.ctx.prec = PRECISION_BITS
    u = sp.symbols("u", positive=True)
    curvature = int(potential[2].subs(u, RAY_START))
    normalized = []
    for order in range(3, 13):
        numerator = flint.arb(int(potential[order].subs(u, RAY_START)))
        denominator = flint.arb(curvature) ** (order // 2)
        if order % 2:
            denominator *= flint.arb(curvature).sqrt()
        normalized.append(numerator / denominator)
    compiled = compile_coefficients()
    expected_signs = {3: -1, 5: -1, 7: 1}
    rows = {}
    for order, expected_sign in expected_signs.items():
        value = evaluate_polynomial_arb(compiled[order], normalized)
        if expected_sign < 0 and not bool(value < 0):
            raise RuntimeError(f"negative endpoint sign failed at order {order}")
        if expected_sign > 0 and not bool(value > 0):
            raise RuntimeError(f"positive endpoint sign failed at order {order}")
        rows[str(order)] = {
            "expected_sign": expected_sign,
            "value_lower": arb_lower_text(value),
            "value_upper": arb_upper_text(value),
        }
    return rows


def leading_bound_gates() -> dict:
    u, potential, coefficients = leading_model()
    v = sp.symbols("v", nonnegative=True)
    substitutions = {u: v + RAY_START}
    gates = {}
    for order, coefficient in coefficients.items():
        floor, ceiling = COEFFICIENT_BOUNDS[order]
        floor_value = sp.Rational(floor.numerator, floor.denominator)
        ceiling_value = sp.Rational(ceiling.numerator, ceiling.denominator)
        buffer = sp.Rational(LEADING_BUFFER.numerator, LEADING_BUFFER.denominator) / u
        sign_gate = None
        if order in (3, 5):
            lower_expression = (-(floor_value + buffer)) ** 2 - coefficient**2
            upper_expression = coefficient**2 - (-(ceiling_value - buffer)) ** 2
            sign_gate = shifted_positive_gate(coefficient**2, (v,), substitutions)
        elif order == 7:
            lower_expression = coefficient**2 - (floor_value + buffer) ** 2
            upper_expression = (ceiling_value - buffer) ** 2 - coefficient**2
            sign_gate = shifted_positive_gate(coefficient**2, (v,), substitutions)
        else:
            lower_expression = coefficient - (floor_value + buffer)
            upper_expression = (ceiling_value - buffer) - coefficient
        gates[str(order)] = {
            "formula": sp.sstr(coefficient),
            "nonvanishing_square_gate": sign_gate,
            "lower_buffer_gate": shifted_positive_gate(
                lower_expression, (v,), substitutions
            ),
            "upper_buffer_gate": shifted_positive_gate(
                upper_expression, (v,), substitutions
            ),
        }

    normalized_caps = {}
    curvature = potential[2]
    cap = 2 - sp.Rational(1, 100) / u
    for order in range(3, 13):
        expression = (
            cap**2 * curvature**order - potential[order] ** 2
            if order % 2
            else cap * curvature ** (order // 2) - potential[order]
        )
        normalized_caps[str(order)] = shifted_positive_gate(
            expression, (v,), substitutions
        )
    return {
        "recurrence": "P_0=1; P_(r+1)=(u/2)*P_r'+2u*P_r",
        "potential_polynomials": {
            str(order): sp.sstr(sp.factor(potential[order]))
            for order in range(2, 13)
        },
        "coefficient_bound_gates": gates,
        "odd_endpoint_signs": odd_endpoint_signs(potential),
        "normalized_jet_cap_gates": normalized_caps,
        "proved_buffer": (
            "coefficient floor+1/(100u) < D_r^(infinity) < "
            "coefficient ceiling-1/(100u)"
        ),
        "proved_normalized_cap": (
            "0<L_r^(infinity)<2-1/(100u), r=3,...,12"
        ),
    }


def potential_remainder_jets() -> tuple[sp.Symbol, sp.Symbol, dict[int, sp.Expr]]:
    u, q = sp.symbols("u q", positive=True)
    remainder = 100 * u**2 - 5 * u - sp.log(2 * q - 3) - sp.log(u)
    jets = {}
    for order in range(1, 13):
        remainder = sp.cancel(
            u * (sp.diff(remainder, u) + 4 * q * sp.diff(remainder, q)) / 2
        )
        jets[order] = remainder
    return u, q, jets


def extended_jet_remainder_gates() -> dict:
    source = load_json(SOURCE_FIRST_ASYMPTOTIC)
    summary = source.get("summary", {})
    reused = summary.get("reused_jet_remainder_sign_gates", 0) + summary.get(
        "new_jet_remainder_sign_gates", 0
    )
    if reused != 18:
        raise RuntimeError("the reused order-two-to-ten jet gates are not closed")
    u, q, remainders = potential_remainder_jets()
    v, capital_q = sp.symbols("v capital_q", nonnegative=True)
    substitutions = {u: v + RAY_START, q: capital_q + Q_FLOOR}
    rows = {}
    for order in (11, 12):
        envelope = JET_REMAINDER_CONSTANT * u**12
        rows[str(order)] = {
            "plus_gate": shifted_positive_gate(
                envelope + remainders[order], (v, capital_q), substitutions
            ),
            "minus_gate": shifted_positive_gate(
                envelope - remainders[order], (v, capital_q), substitutions
            ),
        }
    return {
        "reused_source": SOURCE_FIRST_ASYMPTOTIC.relative_to(REPO_ROOT).as_posix(),
        "reused_formula": "|V^(r)-q*P_r(u)|<=10000*u^10, r=2,...,10",
        "reused_sign_gates": reused,
        "dominance": "10000*u^10 < 100000*u^12 for u>=20",
        "new_formula": "|V^(r)-q*P_r(u)|<=100000*u^12, r=11,12",
        "new_sign_gates": 4,
        "new_orders": rows,
        "combined_formula": "|V^(r)-q*P_r(u)|<=100000*u^12, r=2,...,12",
    }


def polynomial_norms() -> dict:
    source = load_json(SOURCE_SECOND_NEXT)
    symbols = sp.symbols("L_3:13")
    locals_map = {str(symbol): symbol for symbol in symbols}

    def sup_norm(expression: sp.Expr, radius: int = 2) -> sp.Rational:
        return sum(
            abs(sp.Rational(coefficient)) * radius ** sum(powers)
            for powers, coefficient in sp.Poly(expression, *symbols).terms()
        )

    rows = {}
    largest = sp.Rational(0)
    for order_text, row in source["coefficient_rows"].items():
        expression = sp.sympify(row["scaled_coefficient"], locals=locals_map)
        lipschitz = sum(sup_norm(sp.diff(expression, symbol)) for symbol in symbols)
        largest = max(largest, lipschitz)
        rows[order_text] = {"lipschitz_l1_on_abs_L_le_2": str(lipschitz)}
    if not largest < LIPSCHITZ_CAP:
        raise RuntimeError("second-next coefficient Lipschitz cap failed")
    return {
        "domain": "|L_r|<=2, r=3,...,12",
        "rows": rows,
        "largest_lipschitz": str(largest),
        "lipschitz_cap": LIPSCHITZ_CAP,
    }


def scalar_geometry() -> dict:
    flint.ctx.prec = PRECISION_BITS
    q20 = flint.arb.pi() * flint.arb(80).exp()
    if not bool(q20 > Q_FLOOR):
        raise RuntimeError("q floor at u=20 failed")
    eta_endpoint = Fraction(25_000 * 20**10, Q_FLOOR)
    if not eta_endpoint < Fraction(1, 10**17):
        raise RuntimeError("second-next relative jet error cap failed")
    monotonicity_margin = Fraction(4) - Fraction(11, RAY_START)
    if monotonicity_margin <= 0:
        raise RuntimeError("q/u^11 monotonicity failed")

    inverse_value_cap = Fraction(100, 99) ** 6
    inverse_derivative_cap = 6 * Fraction(100, 99) ** 7
    if not (inverse_value_cap < Fraction(11, 10) and inverse_derivative_cap < 7):
        raise RuntimeError("relative-power calculus caps failed")
    endpoint_domain_left = 100 * 20 * 20 * 25_000 * 20**10
    if endpoint_domain_left >= Q_FLOOR:
        raise RuntimeError("full normalized domain transfer failed")

    endpoint_transfer_left = (
        200 * TRANSFER_CONSTANT * 25_000 * RAY_START**11
    )
    if endpoint_transfer_left >= Q_FLOOR:
        raise RuntimeError("second-next transfer endpoint comparison failed")
    return {
        "q_at_20_lower": arb_lower_text(q20),
        "q_floor": str(Q_FLOOR),
        "q_over_u11_log_derivative_lower": str(monotonicity_margin),
        "relative_jet_error": "eta(u)=25000*u^10/q<=10^-17",
        "inverse_power_value_cap": str(inverse_value_cap),
        "inverse_power_derivative_cap": str(inverse_derivative_cap),
        "normalized_jet_transfer": (
            "|L_r-L_r^(infinity)|<=20*eta(u), r=3,...,12"
        ),
        "full_normalized_domain": "|L_r|<2, r=3,...,12",
        "coefficient_transfer": (
            "|D_r-D_r^(infinity)|<=8000000000*eta(u)<1/(200u)"
        ),
        "endpoint_domain_left": str(endpoint_domain_left),
        "endpoint_transfer_left": str(endpoint_transfer_left),
        "endpoint_transfer_right": str(Q_FLOOR),
    }


def build_artifact() -> dict:
    leading = leading_bound_gates()
    jet_remainders = extended_jet_remainder_gates()
    norms = polynomial_norms()
    geometry = scalar_geometry()
    finite = load_json(SOURCE_FINITE)
    if finite.get("parameters", {}).get("block_count") != 3600:
        raise RuntimeError("finite second-next source is not closed")
    rows = [
        RayRow(
            id="co4fcsnparc_01_leading_geometry",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="The q-leading potential jets through order twelve obey one explicit differential recurrence.",
            formula=leading["recurrence"],
            proof_boundary="Exact q-leading potential geometry only.",
            diagnostics=leading["potential_polynomials"],
        ),
        RayRow(
            id="co4fcsnparc_02_leading_coefficient_buffer",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Every leading second-next coefficient lies a full 1/(100u) inside its signed finite bound on u>=20.",
            formula=leading["proved_buffer"],
            proof_boundary="Leading q-infinity coefficient model only.",
            diagnostics=leading["coefficient_bound_gates"],
        ),
        RayRow(
            id="co4fcsnparc_03_odd_sign_and_jet_cap",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="Nonvanishing-square gates and rigorous endpoint signs fix all odd coefficient signs, with a quantitative leading normalized-jet cap through order twelve.",
            formula=leading["proved_normalized_cap"],
            proof_boundary="Leading model only.",
            diagnostics={
                "odd_endpoint_signs": leading["odd_endpoint_signs"],
                "normalized_jet_cap_gates": leading["normalized_jet_cap_gates"],
            },
        ),
        RayRow(
            id="co4fcsnparc_04_extended_jet_remainder",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Four new coefficient-positive gates extend the exact potential-jet envelope from order ten through order twelve.",
            formula=jet_remainders["combined_formula"],
            proof_boundary="Exact potential jets through order twelve on u>=20.",
            diagnostics=jet_remainders,
        ),
        RayRow(
            id="co4fcsnparc_05_polynomial_norm_transfer",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="Explicit coefficient norms control every second-next polynomial under the full-jet perturbation.",
            formula="Lip_1(D_r)<400000000 on |L_3|,...,|L_12|<=2",
            proof_boundary="Finite polynomial norm theorem only.",
            diagnostics=norms,
        ),
        RayRow(
            id="co4fcsnparc_06_exponential_transfer",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The full second-next coefficients differ from their leading models by less than half the leading buffer.",
            formula=geometry["coefficient_transfer"],
            proof_boundary="Exact epsilon-nine/ten formal coefficient layer only.",
            diagnostics=geometry,
        ),
        RayRow(
            id="co4fcsnparc_07_global_coefficient_theorem",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Together with the finite Taylor certificate, all seven signed second-next coefficient bounds hold for every u>=2.",
            formula="signed D_r bounds hold globally on u>=2, r=2,...,8",
            proof_boundary="Explicit epsilon-nine/ten coefficient theorem only.",
            diagnostics={
                "finite_source": SOURCE_FINITE.relative_to(REPO_ROOT).as_posix(),
                "finite_blocks": finite["parameters"]["block_count"],
                "asymptotic_interval": "u>=20",
            },
        ),
        RayRow(
            id="co4fcsnparc_08_beyond_epsilon_ten_handoff",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Control the cancellation-preserving exact central remainder beyond epsilon ten and both adaptive density tails.",
            formula="exact cumulant - kappa_r^[10] = central remainder + two tails",
            proof_boundary="Open exact-density theorem; no exact cumulant corridor is asserted.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate",
        "date": "2026-07-13",
        "status": "exact analytic theorem for the global second-next coefficient layer",
        "proof_boundary": (
            "This artifact proves the seven explicit epsilon-nine/ten scaled "
            "coefficient bounds on u>=20 and composes them with the finite Taylor "
            "certificate on 2<=u<=20. It does not bound the exact central density "
            "remainder beyond epsilon ten or either tail, prove the exact cumulant "
            "corridors, curvature ray, order-four entry, PF-infinity, RH, or Lambda<=0."
        ),
        "leading": leading,
        "jet_remainders": jet_remainders,
        "polynomial_norms": norms,
        "scalar_geometry": geometry,
        "bounds": {
            str(order): [str(floor), str(ceiling)]
            for order, (floor, ceiling) in COEFFICIENT_BOUNDS.items()
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": 7,
            "open_analytic_rows": 1,
            "coefficient_bounds": 7,
            "leading_buffer_gates": 14,
            "odd_nonvanishing_gates": 3,
            "normalized_jet_cap_gates": 10,
            "reused_jet_remainder_sign_gates": 18,
            "new_jet_remainder_sign_gates": 4,
            "global_second_next_layer_closed": True,
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
            "outputs/formal_core.md",
        ],
        "source_hashes": {
            SOURCE_SECOND_NEXT.relative_to(REPO_ROOT).as_posix(): sha256(SOURCE_SECOND_NEXT),
            SOURCE_FINITE.relative_to(REPO_ROOT).as_posix(): sha256(SOURCE_FINITE),
            SOURCE_FIRST_ASYMPTOTIC.relative_to(REPO_ROOT).as_posix(): sha256(SOURCE_FIRST_ASYMPTOTIC),
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.py"
        ),
        "remaining_target": (
            "Use the now-global epsilon-ten subtraction to prove the exact central "
            "remainder and both adaptive density tails."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    geometry = artifact["scalar_geometry"]
    lines = [
        "# Jensen-Window PF Order-Four Second-Next Asymptotic-Ray Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact analytic theorem for the global second-next coefficient layer.",
        "This is not a proof of the exact cumulant corridors, order-four entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.py",
        "```",
        "",
        "## Leading Bounds",
        "",
        "```text",
        artifact["leading"]["recurrence"],
        artifact["leading"]["proved_buffer"],
        artifact["leading"]["proved_normalized_cap"],
        "```",
        "",
        "After `u=20+v`, both sides of all seven buffered coefficient bounds",
        "have coefficient-positive numerators. Odd signs are fixed by",
        "nonvanishing-square gates and rigorous Arb endpoint signs.",
        "",
        "## Exact Jet Transfer",
        "",
        "The eighteen sign gates through order ten are reused. Four new gates give",
        "",
        "```text",
        artifact["jet_remainders"]["combined_formula"],
        geometry["relative_jet_error"],
        geometry["normalized_jet_transfer"],
        geometry["full_normalized_domain"],
        geometry["coefficient_transfer"],
        "```",
        "",
        "The transfer consumes less than half the leading `1/(100u)` buffer.",
        "Together with the 3,600-block finite Taylor certificate, every signed",
        "second-next coefficient bound holds for all `u>=2`.",
        "",
        "## Remaining Boundary",
        "",
        "The complete epsilon-ten subtraction is now globally bounded. What remains",
        "is the exact-density theorem itself: intervalize only the residual after",
        "this subtraction, and prove the left and right adaptive tails separately.",
        "No exact cumulant corridor is promoted here.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.md",
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
        "certified order-four second-next asymptotic ray: "
        "8 rows, 7 exact rows, 7 global coefficient bounds, "
        "4 new jet-remainder sign gates, 1 open exact-density remainder"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
