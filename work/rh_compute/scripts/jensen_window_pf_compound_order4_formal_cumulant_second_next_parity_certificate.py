#!/usr/bin/env python3
"""Certify the epsilon-nine/ten layer after the order-four epsilon-eight model."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
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

import sympy as sp  # noqa: E402

from jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate import (  # noqa: E402
    formal_cumulant_expansion as epsilon_eight_expansion,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.md"
)
SOURCE_NEXT_PARITY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.json"
)
MAX_EPSILON_ORDER = 10
MAX_POTENTIAL_ORDER = 12
TARGET_CUMULANT_ORDERS = tuple(range(2, 9))
SECOND_NEXT_POWER = {2: 10, 3: 9, 4: 10, 5: 9, 6: 10, 7: 9, 8: 10}


@dataclass(frozen=True)
class ParityRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def formal_cumulant_expansion() -> dict[int, list[sp.Expr]]:
    z, y = sp.symbols("z y")
    symbols = {
        order: sp.symbols(f"L_{order}")
        for order in range(3, MAX_POTENTIAL_ORDER + 1)
    }
    perturbation = [sp.Integer(0) for _ in range(MAX_EPSILON_ORDER + 1)]
    for order in range(3, MAX_POTENTIAL_ORDER + 1):
        perturbation[order - 2] = symbols[order] * y**order / sp.factorial(order)

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

    maximum_y_degree = max(sp.Poly(value, y).degree() for value in exponential)
    tilted_gaussian_moments = [sp.Integer(1)]
    for _degree in range(maximum_y_degree):
        previous = tilted_gaussian_moments[-1]
        tilted_gaussian_moments.append(sp.expand(z * previous + sp.diff(previous, z)))

    def tilted_expectation(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                coefficient * tilted_gaussian_moments[monomial[0]]
                for monomial, coefficient in sp.Poly(polynomial, y).terms()
            )
        )

    partition = [tilted_expectation(value) for value in exponential]
    logarithm = [sp.Integer(0) for _ in range(MAX_EPSILON_ORDER + 1)]
    for degree in range(1, MAX_EPSILON_ORDER + 1):
        logarithm[degree] = sp.expand(
            partition[degree]
            - sum(
                index * logarithm[index] * partition[degree - index]
                for index in range(1, degree)
            )
            / degree
        )
    return {
        order: [
            sp.factor(sp.diff(logarithm[degree], z, order).subs(z, 0))
            for degree in range(MAX_EPSILON_ORDER + 1)
        ]
        for order in TARGET_CUMULANT_ORDERS
    }


def audit_epsilon_eight(expansion: dict[int, list[sp.Expr]]) -> dict:
    old = epsilon_eight_expansion()
    comparisons = 0
    for order in TARGET_CUMULANT_ORDERS:
        for degree in range(1, 9):
            if sp.expand(expansion[order][degree] - old[order][degree]) != 0:
                raise RuntimeError(
                    f"epsilon-eight audit failed at kappa_{order}, epsilon^{degree}"
                )
            comparisons += 1
    return {
        "coefficient_comparisons": comparisons,
        "result": "all epsilon^1 through epsilon^8 coefficients agree exactly",
    }


def scaled_coefficient(order: int, coefficient: sp.Expr) -> sp.Expr:
    if order == 2:
        return sp.factor(coefficient)
    return sp.factor(
        sp.Rational((-1) ** order, math.factorial(order - 2)) * coefficient
    )


def coefficient_rows(expansion: dict[int, list[sp.Expr]]) -> dict:
    symbols = sp.symbols("L_3:13")
    rows = {}
    for order in TARGET_CUMULANT_ORDERS:
        nonzero = [degree for degree, value in enumerate(expansion[order]) if value != 0]
        if any(degree % 2 != order % 2 for degree in nonzero):
            raise RuntimeError(f"cumulant parity failed at order {order}")
        next_power = SECOND_NEXT_POWER[order]
        coefficient = sp.factor(expansion[order][next_power])
        scaled = scaled_coefficient(order, coefficient)
        offset = 2 if order == 2 else order - 2
        q_power = (next_power - offset) // 2
        rows[str(order)] = {
            "nonzero_epsilon_powers": nonzero,
            "second_next_epsilon_power": next_power,
            "raw_coefficient": sp.sstr(coefficient),
            "raw_coefficient_terms": len(sp.Poly(coefficient, *symbols).terms()),
            "scaled_coefficient": sp.sstr(scaled),
            "scaled_remainder_leading_power": f"q^-{q_power}",
            "scaled_remainder_identity": (
                f"scaled(kappa_{order}^[10]-kappa_{order}^[8])="
                f"q^-{q_power}*D_{order}(L_3,...,L_12)"
            ),
        }
    return rows


def build_artifact() -> dict:
    expansion = formal_cumulant_expansion()
    audit = audit_epsilon_eight(expansion)
    coefficients = coefficient_rows(expansion)
    hierarchy = {"q^-4": [2, 3, 4], "q^-3": [5, 6], "q^-2": [7, 8]}
    rows = [
        ParityRow(
            id="co4fcsnpc_01_exact_epsilon_ten_expansion",
            role="exact_formal_algebra",
            readiness="ready_to_apply",
            claim="The tilted Gaussian recurrence determines every cumulant coefficient through epsilon order ten exactly.",
            formula="log E exp(zY-R_epsilon(Y)) through epsilon^10 with L_3,...,L_12",
            proof_boundary="Exact finite formal algebra only.",
        ),
        ParityRow(
            id="co4fcsnpc_02_epsilon_eight_audit",
            role="exact_reproduction_gate",
            readiness="ready_to_apply",
            claim="The extension reproduces every epsilon-eight cumulant coefficient term for term.",
            formula="kappa_r^[10] mod epsilon^9 = kappa_r^[8]",
            proof_boundary="Exact cross-expansion consistency gate.",
            diagnostics=audit,
        ),
        ParityRow(
            id="co4fcsnpc_03_second_scaled_hierarchy",
            role="exact_formal_theorem",
            readiness="ready_to_apply",
            claim="The next scaled layer after epsilon eight splits into q^-4, q^-3, and q^-2.",
            formula="orders 2-4: q^-4; orders 5-6: q^-3; orders 7-8: q^-2",
            proof_boundary="Explicit epsilon-nine/ten formal terms only.",
            diagnostics={"hierarchy": hierarchy, "coefficients": coefficients},
        ),
        ParityRow(
            id="co4fcsnpc_04_second_coefficient_handoff",
            role="interval_and_analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Bound the seven explicit second-next coefficient functions globally on u>=2.",
            formula="|D_r(L_3(u),...,L_12(u))|<=B_r, u>=2",
            proof_boundary="Open coefficient-bound theorem.",
        ),
        ParityRow(
            id="co4fcsnpc_05_beyond_epsilon_ten_handoff",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Use the epsilon-ten subtraction in the cancellation-preserving central and two-tail theorem.",
            formula="exact cumulant = kappa_r^[10] + central remainder + two tails",
            proof_boundary="Open exact-density theorem.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate",
        "date": "2026-07-13",
        "status": "exact formal epsilon-ten and second-next parity theorem",
        "proof_boundary": (
            "This artifact proves exact formal algebra through epsilon order ten, "
            "audits the epsilon-eight expansion, and identifies the next scaled "
            "coefficient hierarchy. It does not bound those coefficients globally, "
            "control the exact central remainder or tails, prove the exact cumulant "
            "corridors, curvature ray, order-four entry, PF-infinity, RH, or Lambda<=0."
        ),
        "grading": "lambda_r=L_r*epsilon^(r-2), epsilon=q^(-1/2)",
        "epsilon_eight_audit": audit,
        "coefficient_rows": coefficients,
        "scaled_hierarchy": hierarchy,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": 3,
            "open_rows": 2,
            "cumulant_orders": len(coefficients),
            "epsilon_eight_coefficient_comparisons": audit["coefficient_comparisons"],
            "second_next_coefficients": len(coefficients),
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
            "outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.py"
        ),
        "remaining_target": (
            "Bound the explicit D_r functions, then prove the central remainder "
            "beyond epsilon ten and both adaptive exact-density tails."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Jensen-Window PF Order-Four Formal Cumulant Second-Next Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact formal epsilon-ten and second-next parity theorem.",
        "This is not a proof of the exact cumulant corridors, order-four entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.py",
        "```",
        "",
        "## Exact Extension",
        "",
        "The tilted-Gaussian recurrence now includes normalized potential jets",
        "`L_3,...,L_12` through `epsilon^10`. All 56 coefficients at epsilon powers",
        "one through eight agree exactly with the epsilon-eight expansion.",
        "",
        "After corridor scaling, the new layer is",
        "",
        "```text",
        "orders 2,3,4: scaled(kappa_r^[10]-kappa_r^[8])=q^-4*D_r",
        "orders 5,6:   scaled(kappa_r^[10]-kappa_r^[8])=q^-3*D_r",
        "orders 7,8:   scaled(kappa_r^[10]-kappa_r^[8])=q^-2*D_r",
        "```",
        "",
        "The seven exact rational polynomials contain 30 monomials for odd orders",
        "and 42 for even orders. Their complete formulas are stored in the JSON.",
        "",
        "## Proof Boundary",
        "",
        "This is a finite subtraction layer for the central theorem, not a claim that",
        "the asymptotic expansion converges. Global coefficient bounds and the actual",
        "beyond-epsilon-ten central and two-tail estimates remain open.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
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
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "certified order-four formal cumulant second-next parity: "
        "5 rows, 3 exact rows, 56 epsilon-eight audits, "
        "7 second-next coefficients, 2 open analytic rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
