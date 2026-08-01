#!/usr/bin/env python3
"""Build an exact order-four signed-Hankel nonpromotion gate at quartic contact."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import itertools
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def signed_reshaped_minors(
    coefficients: list[sp.Expr], order: int, base_index: int
) -> list[dict[str, object]]:
    """Enumerate every indexed reshaped minor supported by the finite segment."""
    last = len(coefficients) - 1
    signature = (-1) ** (order * (order - 1) // 2)
    rows: list[dict[str, object]] = []
    for shift in range(len(coefficients)):
        max_column = last - shift - (order - 1)
        if max_column < order - 1:
            continue
        for columns in itertools.combinations(range(max_column + 1), order):
            matrix = sp.Matrix(
                [
                    [coefficients[shift + row + column] for column in columns]
                    for row in range(order)
                ]
            )
            determinant = sp.factor(matrix.det())
            rows.append(
                {
                    "shift": shift + base_index,
                    "columns": list(columns),
                    "determinant": str(determinant),
                    "signed_margin": str(sp.factor(signature * determinant)),
                }
            )
    return rows


def build_exact() -> dict:
    w = sp.symbols("w")
    base_index = 1
    first_contraction_index = 2
    a = sp.Rational(3719, 5000)
    b = sp.Rational(4943, 5000)
    c = sp.factor(4 - 2 * a - b)
    p = sp.factor(b * c)
    if 2 * a + b + c != 4:
        raise RuntimeError("root normalization failed")

    root_A = sp.factor(-3 * a**2 + 8 * a + p)
    root_B = sp.factor(-a**2 + 2 * a + p)
    x = sp.factor(root_A / 6)
    y = sp.factor(18 * a * root_B / root_A**2)
    z = sp.factor(2 * p * root_A / (3 * root_B**2))
    u = sp.Rational(2453, 2500)
    v = sp.Rational(123, 125)
    x7 = sp.Rational(19719, 20000)
    x8 = sp.Rational(987309, 1000000)
    x9 = sp.Rational(98823997, 100000000)
    contractions = [x, y, z, u, v, x7, x8, x9]

    coefficients = [sp.Integer(1), sp.Integer(1)]
    ratio = sp.Integer(1)
    for contraction in contractions:
        ratio = sp.factor(ratio * contraction)
        coefficients.append(sp.factor(coefficients[-1] * ratio))

    quartic = sp.expand(
        sum(sp.binomial(4, j) * coefficients[j] * w**j for j in range(5))
    )
    quartic_factor = sp.expand((1 + a * w) ** 2 * (1 + b * w) * (1 + c * w))
    if sp.expand(quartic - quartic_factor) != 0:
        raise RuntimeError("quartic factorization failed")

    threshold = sp.factor(
        root_B * (3 * a**2 - 5 * a + 5 * p) / (6 * p**2)
    )
    curvature = sp.factor((a - b) * (a - c))
    outer_gap = sp.factor(u - threshold)
    if not (curvature > 0 and outer_gap > 0):
        raise RuntimeError("countermodel is not an outward outer-root contact")

    order_margins = [
        sp.factor(contractions[j + 1] - contractions[j])
        for j in range(len(contractions) - 1)
    ]
    upper_margins = [sp.factor(1 - value) for value in contractions]
    wall_margins = []
    scaled_defects = []
    for offset, value in enumerate(contractions):
        index = first_contraction_index + offset
        wall_margins.append(
            sp.factor(value - sp.Rational(2 * index - 1, 2 * index + 1))
        )
        scaled_defects.append(
            sp.factor(sp.Rational(2 * index + 1, 2) * (1 - value))
        )
    scaled_defect_margins = [
        sp.factor(scaled_defects[j + 1] - scaled_defects[j])
        for j in range(len(scaled_defects) - 1)
    ]
    if not all(value > 0 for value in order_margins + upper_margins + wall_margins):
        raise RuntimeError("ratio walls are not all strict")
    if not all(value > 0 for value in scaled_defect_margins):
        raise RuntimeError("scaled defects are not strictly increasing")

    reciprocal_increment_certificates = []
    for left, right in zip(contractions, contractions[1:]):
        left_defect = sp.factor(1 - left)
        right_defect = sp.factor(1 - right)
        comparison_left = sp.factor(
            1 / right_defect - 1 / left_defect - 1
        )
        squared_margin = sp.factor(4 - left_defect * comparison_left**2)
        if not (comparison_left > 0 and squared_margin > 0):
            raise RuntimeError("reciprocal-defect increment bound failed")
        reciprocal_increment_certificates.append(
            {
                "comparison_left": str(comparison_left),
                "squared_margin": str(squared_margin),
                "conclusion": "1/sqrt(d_next)-1/sqrt(d_current)<1",
            }
        )

    def cubic_frontier(left: sp.Expr, right: sp.Expr) -> sp.Expr:
        return sp.factor(
            left**2 * right**2 - 6 * left * right + 4 * left + 4 * right - 3
        )

    cubic_margins = [
        cubic_frontier(contractions[j], contractions[j + 1])
        for j in range(len(contractions) - 1)
    ]
    if not all(value < 0 for value in cubic_margins):
        raise RuntimeError("an adjacent cubic frontier is not strict")

    minors = {
        str(order): signed_reshaped_minors(coefficients, order, base_index)
        for order in (2, 3, 4)
    }
    expected_counts = {"2": 120, "3": 126, "4": 56}
    for order, entries in minors.items():
        if len(entries) != expected_counts[order]:
            raise RuntimeError(f"unexpected order-{order} minor count")
        if not all(sp.Rational(entry["signed_margin"]) > 0 for entry in entries):
            raise RuntimeError(f"order-{order} signed minor failed")

    quintic = sp.expand(
        sum(sp.binomial(5, j) * coefficients[j] * w**j for j in range(6))
    )
    quintic_at_contact = sp.factor(quintic.subs(w, -1 / a))
    expected_contact = sp.factor(-2 * p**2 * outer_gap / (a**2 * root_B))
    if sp.simplify(quintic_at_contact - expected_contact) != 0:
        raise RuntimeError("quintic contact identity failed")
    quintic_discriminant = sp.factor(sp.discriminant(quintic, w))
    if quintic_discriminant >= 0:
        raise RuntimeError("quintic discriminant is not negative")

    next_x = sp.symbols("x_10")
    next_ratio = sp.factor(ratio * next_x)
    extended_coefficients = [
        *coefficients,
        sp.factor(coefficients[-1] * next_ratio),
    ]

    def contiguous_hankel(order: int, shift: int) -> sp.Expr:
        local_shift = shift - base_index
        return sp.factor(
            sp.det(
                sp.Matrix(
                    [
                        [
                            extended_coefficients[local_shift + row + column]
                            for column in range(order)
                        ]
                        for row in range(order)
                    ]
                )
            )
        )

    next_h3 = contiguous_hankel(3, 7)
    next_h4 = contiguous_hankel(4, 5)
    lower = sp.solve(sp.Eq(next_h3, 0), next_x)[0]
    upper = sp.solve(sp.Eq(next_h4, 0), next_x)[0]
    width = sp.factor(upper - lower)
    continuation_witness = sp.Rational(4944208739, 5000000000)
    if not (
        sp.Poly(next_h3, next_x).LC() < 0
        and sp.Poly(next_h4, next_x).LC() < 0
        and lower < continuation_witness < upper
        and width > 0
    ):
        raise RuntimeError("next signed-Hankel continuation corridor failed")

    previous = contractions[-1]
    next_index = 10
    next_scaled_defect = sp.factor(
        sp.Rational(2 * next_index + 1, 2) * (1 - continuation_witness)
    )
    previous_scaled_defect = scaled_defects[-1]
    next_cubic = cubic_frontier(previous, continuation_witness)
    left_defect = sp.factor(1 - previous)
    right_defect = sp.factor(1 - continuation_witness)
    next_comparison_left = sp.factor(
        1 / right_defect - 1 / left_defect - 1
    )
    next_squared_margin = sp.factor(
        4 - left_defect * next_comparison_left**2
    )
    if not (
        continuation_witness > previous
        and continuation_witness
        > sp.Rational(2 * next_index - 1, 2 * next_index + 1)
        and next_scaled_defect > previous_scaled_defect
        and next_cubic < 0
        and next_comparison_left > 0
        and next_squared_margin > 0
    ):
        raise RuntimeError("next scalar continuation corridor failed")

    return {
        "root_parameters": {
            "base_coefficient_index": base_index,
            "first_contraction_index": first_contraction_index,
            "a": str(a),
            "b": str(b),
            "c": str(c),
            "p": str(p),
            "normalization": "2*a+b+c=4",
            "curvature_C": str(curvature),
        },
        "contractions": {
            f"x_{first_contraction_index + offset}": str(value)
            for offset, value in enumerate(contractions)
        },
        "coefficients_A1_to_A10": [str(value) for value in coefficients],
        "quartic_factorization": (
            "(1+3719*w/5000)^2*(1+4943*w/5000)*(1+7619*w/5000)"
        ),
        "ratio_cone": {
            "successive_increase_margins": [str(value) for value in order_margins],
            "one_minus_coordinate": [str(value) for value in upper_margins],
            "pointwise_wall_margins": [str(value) for value in wall_margins],
            "scaled_defects": [str(value) for value in scaled_defects],
            "scaled_defect_increase_margins": [
                str(value) for value in scaled_defect_margins
            ],
            "reciprocal_defect_increment_certificates": (
                reciprocal_increment_certificates
            ),
        },
        "cubic_frontiers": [str(value) for value in cubic_margins],
        "reshaped_hankel_minors": minors,
        "minor_counts": expected_counts,
        "minimum_signed_minor": {
            order: min(
                (entry for entry in entries),
                key=lambda entry: sp.Rational(entry["signed_margin"]),
            )
            for order, entries in minors.items()
        },
        "outer_threshold": {
            "U": str(threshold),
            "u_minus_U": str(outer_gap),
            "curvature_C": str(curvature),
            "inward_product": str(sp.factor(curvature * outer_gap)),
            "conclusion": "C>0 and u-U>0, so the quartic heat vector is outward",
        },
        "quintic": {
            "P5_at_minus_1_over_a": str(quintic_at_contact),
            "contact_identity": "-2*p^2*(u-U)/(a^2*B)",
            "discriminant": str(quintic_discriminant),
            "discriminant_sign": "negative",
        },
        "next_extension_corridor": {
            "next_contraction": "x_10",
            "order3_lower_root": str(lower),
            "order4_upper_root": str(upper),
            "exact_width": str(width),
            "rational_witness": str(continuation_witness),
            "signed_conditions": "x_10>lower gives H_(3,7)<0; x_10<upper gives H_(4,5)>0",
            "scaled_defect_increase_margin": str(
                sp.factor(next_scaled_defect - previous_scaled_defect)
            ),
            "cubic_frontier": str(next_cubic),
            "reciprocal_increment_squared_margin": str(next_squared_margin),
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    counts = exact["minor_counts"]
    rows = [
        GateRow(
            id="qot4_01_outer_quartic_contact",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="The rational factorization is a normalized hyperbolic quartic boundary point with a double outer root.",
            formula=exact["quartic_factorization"],
            proof_boundary="One abstract quartic contact, not an Xi coefficient window.",
            diagnostics=exact["root_parameters"],
        ),
        GateRow(
            id="qot4_02_contraction_realization",
            role="exact_identity",
            readiness="available_exact",
            claim="The factorization realizes exact contraction coordinates x_2,x_3,x_4; adjoining exact x_5 through x_9 determines a positive segment A_1 through A_10.",
            formula="A_1=A_2=1, r_j=r_(j-1)*x_j, A_(j+1)=A_j*r_j",
            proof_boundary="Finite positive coefficient realization only.",
            diagnostics={
                "contractions": exact["contractions"],
                "coefficients": exact["coefficients_A1_to_A10"],
            },
        ),
        GateRow(
            id="qot4_03_ratio_cone",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="All eight contractions lie strictly between their pointwise Berwald walls and one, increase strictly, have strictly increasing scaled defects, and satisfy every reciprocal-defect increment bound.",
            formula="walls + s_(j+1)>s_j + q_(j+1)-q_j<1, j=2,...,8",
            proof_boundary="Eight local contraction coordinates and seven adjacent scalar-corridor tests only.",
            diagnostics=exact["ratio_cone"],
        ),
        GateRow(
            id="qot4_04_adjacent_cubics",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Every adjacent cubic Jensen frontier supported by the segment is strict.",
            formula="F(x_j,x_(j+1))<0 for j=2,...,8",
            proof_boundary="Seven local cubic tests only.",
            diagnostics={"frontiers": exact["cubic_frontiers"]},
        ),
        GateRow(
            id="qot4_05_order_two_signs",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Every indexed reshaped-Hankel order-two minor supported by A_1 through A_10 has the required negative sign.",
            formula=f"{counts['2']}/{counts['2']} signed order-two margins are positive",
            proof_boundary="Finite-segment minors only.",
            diagnostics=exact["minimum_signed_minor"]["2"],
        ),
        GateRow(
            id="qot4_06_order_three_signs",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Every indexed arbitrary-column reshaped-Hankel order-three minor supported by A_1 through A_10 has the required negative sign.",
            formula=f"{counts['3']}/{counts['3']} signed order-three margins are positive",
            proof_boundary="Finite-segment minors only; not an all-shift infinite theorem.",
            diagnostics=exact["minimum_signed_minor"]["3"],
        ),
        GateRow(
            id="qot4_07_order_four_sign",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Every indexed arbitrary-column reshaped-Hankel order-four minor supported by A_1 through A_10 has the required strict positive sign.",
            formula=f"{counts['4']}/{counts['4']} signed order-four margins are positive",
            proof_boundary="Finite-segment order-four minors only; not an infinite all-shift theorem.",
            diagnostics=exact["minimum_signed_minor"]["4"],
        ),
        GateRow(
            id="qot4_08_outer_threshold_failure",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Despite all preceding strict signs, the live outer-root quartic condition fails.",
            formula="C>0 and u-U(a,p)>0",
            proof_boundary="Local boundary-flow obstruction only.",
            diagnostics=exact["outer_threshold"],
        ),
        GateRow(
            id="qot4_09_quintic_failure",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="The adjacent quintic misses the repeated quartic root and has negative discriminant.",
            formula="P_5(-1/a)=-2*p^2*(u-U)/(a^2*B)<0, Disc(P_5)<0",
            proof_boundary="One degree-five polynomial only.",
            diagnostics=exact["quintic"],
        ),
        GateRow(
            id="qot4_10_forbidden_promotion",
            role="forbidden_promotion",
            readiness="guard_validated",
            claim="The local ratio, scaled-defect, reciprocal-increment, and cubic cones plus every finite signed-Hankel condition through order four visible at contact do not imply u<=U.",
            formula="strengthened local scalar cones + compound order <=4 != outer quartic inward threshold",
            proof_boundary="Does not reject global far-column, infinite-sequence, theta-specific, or degree-five hypotheses.",
        ),
        GateRow(
            id="qot4_11_downstream_handoff",
            role="closed_downstream_handoff",
            readiness="proved_exact",
            claim="The next exact order-three/order-four continuation corridor is nonempty, but the downstream length-13 obstruction proves that no choices through x_12 can preserve both increasing contractions and the next order-four sign at x_13.",
            formula="x_10 in (H3 lower root,H4 upper root), but Delta_13<0 on the complete continuation cube",
            proof_boundary="This closes only the all-length extension of this fixed prefix; a uniform theorem over every outer-contact prefix remains open.",
            diagnostics={
                **exact["next_extension_corridor"],
                "downstream_obstruction": (
                    "outputs/jensen_window_pf_quartic_outer_branch_"
                    "length13_obstruction.md"
                ),
            },
        ),
    ]
    return {
        "kind": "jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate",
        "date": "2026-07-25",
        "status": "exact local compound-order-four nonpromotion gate",
        "proof_boundary": (
            "This finite rational countermodel proves that the local ratio, "
            "scaled-defect, reciprocal-defect-increment, and cubic cones "
            "together with every reshaped-Hankel sign through order four "
            "supported by A_1,...,A_10 do not imply the remaining outer "
            "quartic threshold. It is not an Xi sequence, an infinite "
            "sign-regular countermodel, or a Newman trajectory, and it does not "
            "disprove any global or theta-specific route."
        ),
        "sources": [
            "outputs/jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.md",
            "outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md",
            "outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md",
            "outputs/jensen_window_pf_quartic_quintic_polar_contact_lemma.md",
            "outputs/jensen_window_pf_quartic_outer_branch_length13_obstruction.md",
        ],
        "summary": {
            "rows": 11,
            "exact_quartic_contacts": 1,
            "strict_ratio_coordinates": 8,
            "strict_scaled_defect_steps": 7,
            "reciprocal_defect_increment_bounds": 7,
            "strict_cubic_frontiers": 7,
            "signed_order2_minors": 120,
            "signed_order3_minors": 126,
            "signed_order4_minors": 56,
            "failed_outer_thresholds": 1,
            "negative_quintic_discriminants": 1,
            "forbidden_promotions": 1,
            "closed_downstream_handoffs": 1,
        },
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    threshold = exact["outer_threshold"]
    quintic = exact["quintic"]
    minimum = exact["minimum_signed_minor"]
    extension = exact["next_extension_corridor"]
    return "\n".join(
        [
            "# Quartic Outer-Threshold Order-Four Nonpromotion Gate",
            "",
            "Date: 2026-07-25",
            "",
            "Status: exact finite local scope gate. This is not a proof of",
            "quartic invariance, PF-infinity, `Lambda<=0`, or RH, and it is not",
            "an Xi counterexample or an infinite sign-regular countermodel.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated quartic outer-threshold order-four nonpromotion gate: 11 rows, 0 issues, 1 exact outer quartic contact, 8 strict ratio coordinates, 7 strict scaled-defect steps, 7 reciprocal-defect increment bounds, 7 strict cubic frontiers, 120 signed order-two minors, 126 signed order-three minors, 56 signed order-four minors, 1 failed outer threshold, 1 negative quintic discriminant, 1 forbidden promotion, 1 closed downstream handoff",
            "```",
            "",
            "## Exact Contact",
            "",
            "Take",
            "",
            "```text",
            exact["quartic_factorization"],
            "2*a+b+c=4.",
            "```",
            "",
            "Its contraction coordinates, followed by two exact extensions, are",
            "",
            "```text",
            *(f"{name}={value}" for name, value in exact["contractions"].items()),
            "```",
            "",
            "Regard these as `x_2,...,x_9`. Set `A_1=A_2=1`,",
            "`r_j=r_(j-1)x_j`, and `A_(j+1)=A_j r_j`. This gives a positive",
            "rational segment `A_1,...,A_10`. Every contraction is strictly",
            "below one, strictly above its pointwise wall",
            "`(2j-1)/(2j+1)`, and the eight contractions increase strictly.",
            "In addition, the scaled defects",
            "`s_j=(2j+1)(1-x_j)/2` increase at all seven steps, and exact",
            "one-squaring certificates prove",
            "`(1-x_(j+1))^(-1/2)-(1-x_j)^(-1/2)<1` at all seven steps.",
            "",
            "All seven supported adjacent cubic frontiers are strict:",
            "",
            "```text",
            *(f"F_{index}={value}<0" for index, value in enumerate(exact["cubic_frontiers"], start=1)),
            "```",
            "",
            "## Signed-Hankel Audit",
            "",
            "For every shift and strictly increasing column set whose entries",
            "remain in `A_1,...,A_10`, exact rational enumeration gives",
            "",
            "```text",
            "order 2: 120/120 required signed margins positive",
            "order 3: 126/126 required signed margins positive",
            "order 4:  56/56  required signed margins positive",
            "```",
            "",
            "The closest exact margins are",
            "",
            "```text",
            f"order 2: {minimum['2']['signed_margin']}",
            f"order 3: {minimum['3']['signed_margin']}",
            f"order 4: {minimum['4']['signed_margin']}",
            "```",
            "",
            "Thus the witness includes every local arbitrary-column order-three",
            "test and the available strict contiguous order-four test, not merely",
            "the boundary determinant used to select the outer branch.",
            "",
            "## Failed Outer Threshold",
            "",
            "At the quartic contact,",
            "",
            "```text",
            f"C=(a-b)(a-c)={threshold['curvature_C']}>0",
            f"U(a,p)={threshold['U']}",
            f"u-U={threshold['u_minus_U']}>0.",
            "```",
            "",
            "The repeated root is on the live outer branch, but the required",
            "`u<=U(a,p)` condition fails strictly. The heat vector is outward.",
            "",
            "The adjacent quintic exposes the missing structure exactly:",
            "",
            "```text",
            "P_5(-1/a)=-2*p^2*(u-U)/(a^2*B)",
            f"             ={quintic['P5_at_minus_1_over_a']}<0",
            f"Disc(P_5)={quintic['discriminant']}<0.",
            "```",
            "",
            "## Consequence",
            "",
            "The completed local signed-Hankel information through order four",
            "cannot by itself prove the remaining quartic inward threshold.",
            "This does not challenge the Xi theorems: the witness is finite and",
            "is neither the Xi sequence nor a Newman trajectory. It instead",
            "forces the next proof step to use genuinely nonlocal/far-column",
            "information, a degree-five polar closure, theta arithmetic, or",
            "another Xi-specific heat-compatible invariant.",
            "",
            "The immediate continuation problem is exact rather than merely",
            "numerical. Appending `x_10` preserves the new terminal order-three",
            "and order-four signs whenever",
            "",
            "```text",
            f"{extension['order3_lower_root']} < x_10",
            f"x_10 < {extension['order4_upper_root']}",
            f"corridor width={extension['exact_width']}>0.",
            "```",
            "",
            f"The rational point `x_10={extension['rational_witness']}` lies",
            "strictly inside and also preserves the next pointwise wall,",
            "scaled-defect increase, reciprocal-increment bound, and cubic",
            "frontier. The exact downstream obstruction now settles this",
            "particular all-length question: after arbitrary strict corridor",
            "choices through `x_12`, the best possible compatibility margin",
            "for an increasing `x_13` and the next order-four sign is still",
            "negative. Thus this fixed prefix cannot be an infinite",
            "countermodel. A uniform theorem over all outer-contact prefixes,",
            "or a genuinely global Xi, degree-five, or theta-specific",
            "constraint, remains open.",
            "",
            "```text",
            "outputs/jensen_window_pf_quartic_outer_branch_length13_obstruction.md",
            "```",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(f"wrote quartic outer-threshold order-four gate: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
