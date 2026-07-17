#!/usr/bin/env python3
"""Certify the next omitted parity terms in the order-four cumulant expansion."""

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


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.md"
)
SOURCE_TARGET = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.json"
)
MAX_EPSILON_ORDER = 8
MAX_POTENTIAL_ORDER = 10
TARGET_CUMULANT_ORDERS = tuple(range(2, 9))
NEXT_PARITY_POWER = {2: 8, 3: 7, 4: 8, 5: 7, 6: 8, 7: 7, 8: 8}


@dataclass(frozen=True)
class ParityRow:
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


def formal_cumulant_expansion() -> dict[int, list[sp.Expr]]:
    """Return exact epsilon coefficients through order eight for kappa_2,...,kappa_8."""
    epsilon, z, y = sp.symbols("epsilon z y")
    del epsilon
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


def audit_epsilon_six(expansion: dict[int, list[sp.Expr]]) -> dict:
    source = load_json(SOURCE_TARGET)
    symbols = sp.symbols("L_3:11")
    locals_map = {str(symbol): symbol for symbol in symbols}
    comparisons = 0
    for order in TARGET_CUMULANT_ORDERS:
        source_terms = {
            int(term["epsilon_power"]): sp.sympify(
                term["coefficient"], locals=locals_map
            )
            for term in source["exact"]["cumulants"][str(order)]["terms"]
        }
        for degree in range(1, 7):
            expected = source_terms.get(degree, sp.Integer(0))
            if sp.expand(expansion[order][degree] - expected) != 0:
                raise RuntimeError(
                    f"epsilon-six audit failed at kappa_{order}, epsilon^{degree}"
                )
            comparisons += 1
    return {
        "source": SOURCE_TARGET.relative_to(REPO_ROOT).as_posix(),
        "coefficient_comparisons": comparisons,
        "result": "all epsilon^1 through epsilon^6 coefficients agree exactly",
    }


def scaled_coefficient(order: int, coefficient: sp.Expr) -> sp.Expr:
    if order == 2:
        return sp.factor(coefficient)
    return sp.factor(sp.Rational((-1) ** order, math.factorial(order - 2)) * coefficient)


def expansion_rows(expansion: dict[int, list[sp.Expr]]) -> dict:
    rows = {}
    for order in TARGET_CUMULANT_ORDERS:
        nonzero = [degree for degree, value in enumerate(expansion[order]) if value != 0]
        parity = order % 2
        if any(degree % 2 != parity for degree in nonzero):
            raise RuntimeError(f"cumulant parity failed at order {order}")
        expected_next = NEXT_PARITY_POWER[order]
        if expected_next not in nonzero:
            raise RuntimeError(f"missing next parity coefficient at order {order}")
        coefficient = sp.factor(expansion[order][expected_next])
        scaled = scaled_coefficient(order, coefficient)
        offset = 2 if order == 2 else order - 2
        q_power = (expected_next - offset) // 2
        rows[str(order)] = {
            "nonzero_epsilon_powers": nonzero,
            "next_epsilon_power": expected_next,
            "raw_coefficient": sp.sstr(coefficient),
            "raw_coefficient_terms": len(sp.Poly(coefficient, *sp.symbols("L_3:11")).terms()),
            "scaled_coefficient": sp.sstr(scaled),
            "scaled_remainder_leading_power": f"q^-{q_power}",
            "scaled_remainder_identity": (
                f"scaled(kappa_{order}^[8]-kappa_{order}^[6])="
                f"q^-{q_power}*C_{order}(L_3,...,L_10)"
            ),
        }
    return rows


def build_artifact() -> dict:
    expansion = formal_cumulant_expansion()
    audit = audit_epsilon_six(expansion)
    coefficients = expansion_rows(expansion)
    hierarchy = {
        "q^-3": [2, 3, 4],
        "q^-2": [5, 6],
        "q^-1": [7, 8],
    }
    rows = [
        ParityRow(
            id="co4fcnpc_01_exact_epsilon_eight_expansion",
            role="exact_formal_algebra",
            readiness="ready_to_apply",
            claim="The tilted Gaussian recurrence determines every cumulant coefficient through epsilon order eight exactly.",
            formula="log E exp(zY-R_epsilon(Y)) through epsilon^8 with L_3,...,L_10",
            proof_boundary="Exact finite formal algebra only; no analytic density remainder.",
            diagnostics={"maximum_epsilon_order": 8, "maximum_potential_jet": 10},
        ),
        ParityRow(
            id="co4fcnpc_02_epsilon_six_audit",
            role="exact_reproduction_gate",
            readiness="ready_to_apply",
            claim="The new expansion reproduces every stored epsilon-six cumulant coefficient term for term.",
            formula="kappa_r^[8] mod epsilon^7 = kappa_r^[6], r=2,...,8",
            proof_boundary="Exact cross-artifact consistency gate.",
            diagnostics=audit,
        ),
        ParityRow(
            id="co4fcnpc_03_parity_theorem",
            role="exact_formal_theorem",
            readiness="ready_to_apply",
            claim="Gaussian reflection parity removes every cumulant coefficient of the wrong epsilon parity through order eight.",
            formula="[epsilon^j] kappa_r=0 when j is not congruent to r modulo 2",
            proof_boundary="Verified exact coefficients through epsilon order eight.",
            diagnostics={"orders": list(TARGET_CUMULANT_ORDERS)},
        ),
        ParityRow(
            id="co4fcnpc_04_scaled_hierarchy",
            role="exact_formal_theorem",
            readiness="ready_to_apply",
            claim="The first omitted scaled terms split into the q^-3, q^-2, and q^-1 hierarchy.",
            formula="orders 2-4: q^-3; orders 5-6: q^-2; orders 7-8: q^-1",
            proof_boundary="First omitted formal term only; the beyond-epsilon-eight remainder remains open.",
            diagnostics={"hierarchy": hierarchy, "coefficients": coefficients},
        ),
        ParityRow(
            id="co4fcnpc_05_coefficient_bound_handoff",
            role="interval_and_analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Bound every explicit next-parity coefficient on the finite interval and asymptotic ray.",
            formula="|C_r(L_3(u),...,L_10(u))| <= B_r(u), u>=2",
            proof_boundary="Open coefficient-bound theorem; no bounds are asserted here.",
        ),
        ParityRow(
            id="co4fcnpc_06_exact_remainder_handoff",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Bound the exact density after subtracting the epsilon-eight central expansion, including both tails.",
            formula="exact cumulant = kappa_r^[8] + central remainder + two tails",
            proof_boundary="Open analytic theorem; this artifact does not prove an exact cumulant corridor.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate",
        "date": "2026-07-13",
        "status": "exact formal epsilon-eight and next-parity coefficient theorem",
        "proof_boundary": (
            "This artifact proves exact formal algebra through epsilon order eight, "
            "audits the prior epsilon-six artifact, and identifies the first omitted "
            "scaled coefficient hierarchy. It does not bound those coefficients on "
            "u>=2, control the exact density remainder or tails, prove the exact "
            "curvature ray, order-four entry, PF-infinity, RH, or Lambda<=0."
        ),
        "grading": "lambda_r=L_r*epsilon^(r-2), epsilon=q^(-1/2)",
        "epsilon_six_audit": audit,
        "coefficient_rows": coefficients,
        "scaled_hierarchy": hierarchy,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": 4,
            "open_rows": 2,
            "cumulant_orders": len(coefficients),
            "epsilon_six_coefficient_comparisons": audit["coefficient_comparisons"],
            "next_parity_coefficients": len(coefficients),
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.py"
        ),
        "remaining_target": (
            "Certify the explicit C_r coefficient functions on u>=2, then prove the "
            "central beyond-epsilon-eight remainder and both exact-density tails."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    lines = [
        "# Jensen-Window PF Order-Four Formal Cumulant Next-Parity Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact formal epsilon-eight and next-parity coefficient theorem.",
        "This is not a proof of the exact cumulant ray, order-four entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.py",
        "```",
        "",
        "## Exact Audit",
        "",
        "The tilted-Gaussian recurrence was extended from `epsilon^6` to",
        "`epsilon^8`, with normalized potential jets `L_3,...,L_10`. All 42",
        "coefficients at epsilon powers one through six agree exactly with the",
        "stored formal-cumulant target.",
        "",
        "Reflection parity gives",
        "",
        "```text",
        "[epsilon^j] kappa_r=0 when j is not congruent to r modulo 2.",
        "```",
        "",
        "## First Omitted Layer",
        "",
        "After applying the corridor normalization, the epsilon-eight extension is",
        "organized as",
        "",
        "```text",
        "orders 2,3,4: scaled(kappa_r^[8]-kappa_r^[6])=q^-3*C_r",
        "orders 5,6:   scaled(kappa_r^[8]-kappa_r^[6])=q^-2*C_r",
        "orders 7,8:   scaled(kappa_r^[8]-kappa_r^[6])=q^-1*C_r",
        "```",
        "",
        "The exact rational polynomials `C_r(L_3,...,L_10)` are stored in the JSON",
        "artifact. In particular, the previously observed order-seven and",
        "order-eight `q^-1` discrepancy is the predicted first omitted parity term,",
        "not a failure of the epsilon-six algebra.",
        "",
        "## Remaining Boundary",
        "",
        "No exact-density remainder has been promoted. The next gates are to certify",
        "the seven explicit coefficient functions on `2<=u<=20` and `u>=20`, then",
        "bound the central remainder beyond epsilon eight and both adaptive tails.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md",
        "outputs/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.md",
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
        "certified order-four formal cumulant next parity: "
        "6 rows, 4 exact rows, 42 epsilon-six audits, "
        "7 next-parity coefficients, 2 open analytic rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
