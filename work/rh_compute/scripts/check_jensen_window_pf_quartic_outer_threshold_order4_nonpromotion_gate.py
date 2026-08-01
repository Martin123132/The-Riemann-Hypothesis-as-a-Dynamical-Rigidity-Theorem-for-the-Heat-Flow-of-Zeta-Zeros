#!/usr/bin/env python3
"""Check the exact quartic outer-threshold order-four nonpromotion gate."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.md"
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    note = args.note.read_text(encoding="utf-8")
    issues: list[str] = []

    if payload.get("kind") != (
        "jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate"
    ):
        issues.append("artifact kind changed")
    if payload.get("status") != "exact local compound-order-four nonpromotion gate":
        issues.append("artifact status changed")
    rows = payload.get("rows", [])
    if len(rows) != 11:
        issues.append(f"expected 11 rows, found {len(rows)}")
    for index, row in enumerate(rows, start=1):
        if not row.get("id", "").startswith(f"qot4_{index:02d}_"):
            issues.append(f"row order/id changed at {index}")

    expected_summary = {
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
    }
    for key, expected in expected_summary.items():
        if payload.get("summary", {}).get(key) != expected:
            issues.append(f"summary {key} changed")

    w = sp.symbols("w")
    first_contraction_index = 2
    a = sp.Rational(3719, 5000)
    b = sp.Rational(4943, 5000)
    c = sp.Rational(7619, 5000)
    p = b * c
    root_A = -3 * a**2 + 8 * a + p
    root_B = -a**2 + 2 * a + p
    contractions = [
        sp.factor(root_A / 6),
        sp.factor(18 * a * root_B / root_A**2),
        sp.factor(2 * p * root_A / (3 * root_B**2)),
        sp.Rational(2453, 2500),
        sp.Rational(123, 125),
        sp.Rational(19719, 20000),
        sp.Rational(987309, 1000000),
        sp.Rational(98823997, 100000000),
    ]
    coefficients = [sp.Integer(1), sp.Integer(1)]
    ratio = sp.Integer(1)
    for contraction in contractions:
        ratio = sp.factor(ratio * contraction)
        coefficients.append(sp.factor(coefficients[-1] * ratio))

    stored_coefficients = [
        sp.Rational(value)
        for value in payload.get("exact", {}).get("coefficients_A1_to_A10", [])
    ]
    if stored_coefficients != coefficients:
        issues.append("stored coefficient segment changed")

    quartic = sp.expand(
        sum(sp.binomial(4, j) * coefficients[j] * w**j for j in range(5))
    )
    factored = sp.expand((1 + a * w) ** 2 * (1 + b * w) * (1 + c * w))
    if quartic != factored or 2 * a + b + c != 4:
        issues.append("independent quartic contact check failed")

    scaled_defects = []
    for offset, value in enumerate(contractions):
        index = first_contraction_index + offset
        if not value > sp.Rational(2 * index - 1, 2 * index + 1):
            issues.append(f"pointwise wall failed at contraction {index}")
        if not value < 1:
            issues.append(f"upper ratio wall failed at contraction {index}")
        scaled_defects.append(
            sp.factor(sp.Rational(2 * index + 1, 2) * (1 - value))
        )
    if not all(
        contractions[j] < contractions[j + 1]
        for j in range(len(contractions) - 1)
    ):
        issues.append("contractions are not strictly increasing")
    if not all(
        scaled_defects[j] < scaled_defects[j + 1]
        for j in range(len(scaled_defects) - 1)
    ):
        issues.append("scaled defects are not strictly increasing")
    for index, (left, right) in enumerate(
        zip(contractions, contractions[1:]), start=first_contraction_index
    ):
        left_defect = sp.factor(1 - left)
        right_defect = sp.factor(1 - right)
        comparison_left = sp.factor(
            1 / right_defect - 1 / left_defect - 1
        )
        squared_margin = sp.factor(4 - left_defect * comparison_left**2)
        if not (comparison_left > 0 and squared_margin > 0):
            issues.append(
                f"reciprocal-defect increment certificate failed at {index}"
            )

    def cubic_frontier(left: sp.Expr, right: sp.Expr) -> sp.Expr:
        return sp.factor(
            left**2 * right**2 - 6 * left * right + 4 * left + 4 * right - 3
        )

    if not all(
        cubic_frontier(contractions[j], contractions[j + 1]) < 0
        for j in range(len(contractions) - 1)
    ):
        issues.append("an adjacent cubic frontier is not strict")

    counts: dict[int, int] = {}
    for order in (2, 3, 4):
        signature = (-1) ** (order * (order - 1) // 2)
        count = 0
        for shift in range(len(coefficients)):
            max_column = len(coefficients) - 1 - shift - (order - 1)
            if max_column < order - 1:
                continue
            for columns in itertools.combinations(range(max_column + 1), order):
                matrix = sp.Matrix(
                    [
                        [
                            coefficients[shift + row + column]
                            for column in columns
                        ]
                        for row in range(order)
                    ]
                )
                if not signature * sp.factor(matrix.det()) > 0:
                    issues.append(
                        f"wrong order-{order} sign at shift {shift}, columns {columns}"
                    )
                count += 1
        counts[order] = count
    if counts != {2: 120, 3: 126, 4: 56}:
        issues.append(f"minor counts changed: {counts}")

    threshold = sp.factor(
        root_B * (3 * a**2 - 5 * a + 5 * p) / (6 * p**2)
    )
    curvature = sp.factor((a - b) * (a - c))
    outer_gap = sp.factor(contractions[3] - threshold)
    if not (curvature > 0 and outer_gap > 0):
        issues.append("outer-threshold failure was not reproduced")
    stored_gap = sp.Rational(payload["exact"]["outer_threshold"]["u_minus_U"])
    if stored_gap != outer_gap:
        issues.append("stored outer-threshold gap changed")

    quintic = sp.expand(
        sum(sp.binomial(5, j) * coefficients[j] * w**j for j in range(6))
    )
    contact = sp.factor(quintic.subs(w, -1 / a))
    expected_contact = sp.factor(-2 * p**2 * outer_gap / (a**2 * root_B))
    if contact != expected_contact or not contact < 0:
        issues.append("quintic polar-contact identity failed")
    quintic_discriminant = sp.factor(sp.discriminant(quintic, w))
    if not quintic_discriminant < 0:
        issues.append("quintic discriminant is not negative")
    if sp.Rational(payload["exact"]["quintic"]["discriminant"]) != quintic_discriminant:
        issues.append("stored quintic discriminant changed")

    next_x = sp.symbols("x_10")
    extended_coefficients = [
        *coefficients,
        sp.factor(coefficients[-1] * ratio * next_x),
    ]

    def contiguous_hankel(order: int, shift: int) -> sp.Expr:
        local_shift = shift - 1
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
    witness = sp.Rational(4944208739, 5000000000)
    if not (
        sp.Poly(next_h3, next_x).LC() < 0
        and sp.Poly(next_h4, next_x).LC() < 0
        and lower < witness < upper
    ):
        issues.append("next signed-Hankel continuation corridor failed")
    stored_corridor = payload["exact"].get("next_extension_corridor", {})
    if (
        sp.Rational(stored_corridor.get("order3_lower_root", "0")) != lower
        or sp.Rational(stored_corridor.get("order4_upper_root", "0")) != upper
        or sp.Rational(stored_corridor.get("exact_width", "0")) != upper - lower
    ):
        issues.append("stored next-extension corridor changed")

    previous = contractions[-1]
    next_scaled = sp.factor(sp.Rational(21, 2) * (1 - witness))
    previous_scaled = sp.factor(sp.Rational(19, 2) * (1 - previous))
    if not (
        witness > previous
        and witness > sp.Rational(19, 21)
        and next_scaled > previous_scaled
        and cubic_frontier(previous, witness) < 0
    ):
        issues.append("next scalar continuation witness failed")
    left_defect = sp.factor(1 - previous)
    right_defect = sp.factor(1 - witness)
    comparison_left = sp.factor(
        1 / right_defect - 1 / left_defect - 1
    )
    squared_margin = sp.factor(4 - left_defect * comparison_left**2)
    if not (comparison_left > 0 and squared_margin > 0):
        issues.append("next reciprocal-increment witness failed")

    source_markers = {
        "outputs/jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.md": [
            "u<=U(a,p)",
            "D_(3,n)(lambda)<0",
        ],
        "outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md": [
            "H_(4,n)(lambda)>0",
        ],
        "outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md": [
            "R_(4,n)(j_1,j_2,j_3,j_4)>0",
        ],
        "outputs/jensen_window_pf_quartic_quintic_polar_contact_lemma.md": [
            "P_4=P_5-(w/5)*P_5'",
            "u=U(a,p)",
        ],
        "outputs/jensen_window_pf_quartic_outer_branch_length13_obstruction.md": [
            "Delta_13=C_13-M_13",
            "cannot continue",
        ],
    }
    for source, markers in source_markers.items():
        path = REPO_ROOT / source
        if not path.exists():
            issues.append(f"source missing: {source}")
            continue
        text = path.read_text(encoding="utf-8")
        for marker in markers:
            if marker not in text:
                issues.append(f"source marker missing in {source}: {marker}")

    required_note_markers = [
        "120/120 required signed margins positive",
        "126/126 required signed margins positive",
        "56/56  required signed margins positive",
        "scaled defects",
        "one-squaring certificates",
        "u-U=",
        "P_5(-1/a)=-2*p^2*(u-U)/(a^2*B)",
        "Disc(P_5)=",
        "corridor width=",
        "particular all-length question",
        "fixed prefix cannot be an infinite",
        "cannot by itself prove",
        "neither the Xi sequence nor a Newman trajectory",
    ]
    for marker in required_note_markers:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated quartic outer-threshold order-four nonpromotion gate: "
        "11 rows, 0 issues, 1 exact outer quartic contact, 8 strict ratio "
        "coordinates, 7 strict scaled-defect steps, 7 reciprocal-defect "
        "increment bounds, 7 strict cubic frontiers, 120 signed order-two minors, "
        "126 signed order-three minors, 56 signed order-four minors, 1 failed "
        "outer threshold, 1 negative quintic discriminant, 1 forbidden "
        "promotion, 1 closed downstream handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
