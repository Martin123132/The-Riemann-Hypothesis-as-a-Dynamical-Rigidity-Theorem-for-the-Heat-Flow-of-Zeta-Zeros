#!/usr/bin/env python3
"""Build an exact outer-contact tail surviving through x_14."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import sys


if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


PRECISION_BITS = 512
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_contact_length14_survivor_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_quartic_outer_contact_length14_survivor_gate.md"
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


def lower_process_priority() -> None:
    try:
        import psutil

        psutil.Process(os.getpid()).nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass


def arb_exact(value: Fraction) -> flint.arb:
    return flint.arb(value.numerator) / flint.arb(value.denominator)


def strictly_positive(value: flint.arb) -> bool:
    return value.lower() > 0


def fraction_sha256(value: Fraction) -> str:
    return hashlib.sha256(
        f"{value.numerator}/{value.denominator}".encode("ascii")
    ).hexdigest()


def sequence_sha256(values: list[Fraction]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(
            f"{value.numerator}/{value.denominator}\n".encode("ascii")
        )
    return digest.hexdigest()


def decimal_summary(value: Fraction) -> dict[str, object]:
    ball = arb_exact(value)
    return {
        "decimal_ball": ball.str(n=80, more=True),
        "sign": 1 if value > 0 else -1 if value < 0 else 0,
        "numerator_digits": len(str(abs(value.numerator))),
        "denominator_digits": len(str(value.denominator)),
        "sha256": fraction_sha256(value),
    }


def signed_minor_summary(
    coefficients: list[Fraction], order: int
) -> dict[str, object]:
    values = [arb_exact(value) for value in coefficients]
    signature = (-1) ** (order * (order - 1) // 2)
    permutations = list(itertools.permutations(range(order)))
    digest = hashlib.sha256()
    count = 0
    minimum_lower = math.inf
    minimum_location: dict[str, object] | None = None
    for shift in range(len(values)):
        max_column = len(values) - 1 - shift - (order - 1)
        if max_column < order - 1:
            continue
        for columns in itertools.combinations(
            range(max_column + 1), order
        ):
            determinant = flint.arb(0)
            for permutation in permutations:
                inversions = sum(
                    permutation[left] > permutation[right]
                    for left in range(order)
                    for right in range(left + 1, order)
                )
                term = flint.arb(-1 if inversions % 2 else 1)
                for row in range(order):
                    term *= values[
                        shift + row + columns[permutation[row]]
                    ]
                determinant += term
            signed = signature * determinant
            if not strictly_positive(signed):
                raise RuntimeError(
                    f"order-{order} minor is not certified positive at "
                    f"shift {shift + 1}, columns {columns}"
                )
            lower = float(signed.lower())
            if lower < minimum_lower:
                minimum_lower = lower
                minimum_location = {
                    "shift": shift + 1,
                    "columns": list(columns),
                    "lower_ball": signed.lower().str(n=80, more=True),
                }
            digest.update(
                (
                    f"{shift + 1}|"
                    f"{','.join(str(column) for column in columns)}|"
                    f"{signed.str(n=80, more=True)}\n"
                ).encode("ascii")
            )
            count += 1
    return {
        "precision_bits": PRECISION_BITS,
        "count": count,
        "strict_positive": count,
        "minimum_lower_log10": f"{math.log10(minimum_lower):.12f}",
        "minimum_location": minimum_location,
        "arb_sequence_sha256": digest.hexdigest(),
    }


def build_exact() -> dict[str, object]:
    delta = Fraction(13, 200)
    q = Fraction(1, 2)
    epsilon = Fraction(1, 100000)
    corridor = Fraction(99, 100)
    a = 1 - delta
    curvature = 4 * q * delta**2
    p = 4 * a - 3 * a**2 + curvature
    root_A = -3 * a**2 + 8 * a + p
    root_B = -a**2 + 2 * a + p
    threshold = (
        root_B * (3 * a**2 - 5 * a + 5 * p) / (6 * p**2)
    )
    threshold_defect = (1 - threshold) / delta**2
    s = threshold_defect - epsilon
    simple_discriminant = (4 - 2 * a) ** 2 - 4 * p

    x: dict[int, Fraction] = {
        2: root_A / 6,
        3: 18 * a * root_B / root_A**2,
        4: 2 * p * root_A / (3 * root_B**2),
        5: 1 - delta**2 * s,
    }

    def defect(index: int) -> Fraction:
        return 1 - x[index]

    def gap(index: int) -> Fraction:
        return (
            defect(index + 2) ** 2
            - x[index + 2] ** 2
            * defect(index + 1)
            * defect(index + 3)
        )

    def cap(index: int) -> Fraction:
        return (
            gap(index - 4) ** 2
            / (x[index - 2] ** 3 * gap(index - 5))
        )

    def compatibility(index: int) -> Fraction:
        repeated = (
            defect(index - 1) ** 2
            - x[index - 1] ** 2
            * defect(index - 2)
            * defect(index - 1)
        )
        return cap(index) - repeated

    extension_rows: list[dict[str, object]] = []
    for index in range(6, 15):
        current_cap = cap(index)
        target_gap = corridor * current_cap
        current_defect = (
            defect(index - 1) ** 2 - target_gap
        ) / (
            x[index - 1] ** 2 * defect(index - 2)
        )
        x[index] = 1 - current_defect
        if gap(index - 3) != target_gap:
            raise RuntimeError(f"corridor recurrence failed at x_{index}")
        extension_rows.append(
            {
                "index": index,
                "parameter": "99/100",
                "contraction_decimal": arb_exact(x[index]).str(
                    n=40, more=True
                ),
                "contraction_sha256": fraction_sha256(x[index]),
            }
        )

    contractions = [x[index] for index in range(2, 15)]
    increase_margins = [
        right - left
        for left, right in zip(contractions, contractions[1:])
    ]
    upper_margins = [1 - value for value in contractions]
    wall_margins = [
        x[index] - Fraction(2 * index - 1, 2 * index + 1)
        for index in range(2, 15)
    ]
    scaled_margins = [
        (2 * index + 3) * defect(index + 1)
        - (2 * index + 1) * defect(index)
        for index in range(2, 14)
    ]
    if not all(
        value > 0
        for value in (
            increase_margins
            + upper_margins
            + wall_margins
            + scaled_margins
        )
    ):
        raise RuntimeError("a scalar corridor failed")

    cubic_margins: list[Fraction] = []
    reciprocal_rows: list[dict[str, object]] = []
    for index in range(2, 14):
        left = x[index]
        right = x[index + 1]
        cubic = (
            left**2 * right**2
            - 6 * left * right
            + 4 * left
            + 4 * right
            - 3
        )
        if cubic >= 0:
            raise RuntimeError("an adjacent cubic frontier failed")
        cubic_margins.append(-cubic)
        left_defect = defect(index)
        right_defect = defect(index + 1)
        comparison = (
            1 / right_defect - 1 / left_defect - 1
        )
        if comparison <= 0:
            reciprocal_rows.append(
                {
                    "index": index,
                    "branch": "automatic_nonpositive_comparison",
                }
            )
            continue
        squared_margin = 4 - left_defect * comparison**2
        if squared_margin <= 0:
            raise RuntimeError("a reciprocal-defect increment failed")
        reciprocal_rows.append(
            {
                "index": index,
                "branch": "positive_comparison_one_squaring",
                "squared_margin_sha256": fraction_sha256(squared_margin),
            }
        )

    gap_values = [gap(index) for index in range(1, 12)]
    order4_margins = [
        cap(index) - gap(index - 3) for index in range(6, 15)
    ]
    if not all(value > 0 for value in gap_values + order4_margins):
        raise RuntimeError("a contiguous signed layer failed")

    delta14 = compatibility(14)
    delta15 = compatibility(15)
    if not (delta14 > 0 and delta15 > 0):
        raise RuntimeError("terminal compatibility margins changed")

    coefficients = [Fraction(1), Fraction(1)]
    ratio = Fraction(1)
    for contraction in contractions:
        ratio *= contraction
        coefficients.append(coefficients[-1] * ratio)
    if len(coefficients) != 15:
        raise RuntimeError("expected A_1 through A_15")
    minor_summaries = {
        str(order): signed_minor_summary(coefficients, order)
        for order in range(2, 9)
    }
    expected_counts = {
        "2": 455,
        "3": 1001,
        "4": 1287,
        "5": 924,
        "6": 330,
        "7": 45,
        "8": 1,
    }
    if {
        order: summary["count"]
        for order, summary in minor_summaries.items()
    } != expected_counts:
        raise RuntimeError("arbitrary-column minor counts changed")

    outward_gap = x[5] - threshold
    if outward_gap != delta**2 * epsilon:
        raise RuntimeError("outward threshold slack changed")
    if not (
        delta > 0
        and q > 0
        and q < 1
        and simple_discriminant > 0
        and curvature > 0
        and outward_gap > 0
    ):
        raise RuntimeError("outer contact is not strict")

    quintic_coefficients = [
        Fraction(1),
        Fraction(5),
        10 * x[2],
        10 * x[2] ** 2 * x[3],
        5 * x[2] ** 3 * x[3] ** 2 * x[4],
        x[2] ** 4 * x[3] ** 3 * x[4] ** 2 * x[5],
    ]
    double_root = -1 / a
    quintic_value = sum(
        coefficient * double_root**degree
        for degree, coefficient in enumerate(quintic_coefficients)
    )
    expected_quintic_value = (
        -2 * p**2 * outward_gap / (a**2 * root_B)
    )
    if quintic_value != expected_quintic_value or quintic_value >= 0:
        raise RuntimeError("adjacent quintic contact value changed")
    quintic_poly = flint.fmpq_poly(
        [
            flint.fmpq(value.numerator, value.denominator)
            for value in quintic_coefficients
        ]
    )
    quintic_discriminant_fmpq = (
        quintic_poly.resultant(quintic_poly.derivative())
        / quintic_poly[5]
    )
    quintic_discriminant = Fraction(
        int(quintic_discriminant_fmpq.p),
        int(quintic_discriminant_fmpq.q),
    )
    if quintic_discriminant >= 0:
        raise RuntimeError("adjacent quintic discriminant is not negative")

    return {
        "contact": {
            "delta": str(delta),
            "q": str(q),
            "a": str(a),
            "p": str(p),
            "curvature_C": str(curvature),
            "simple_root_discriminant": str(simple_discriminant),
            "threshold_U": str(threshold),
            "threshold_defect_T": str(threshold_defect),
            "s": str(s),
            "epsilon": str(epsilon),
            "x5_minus_U": str(outward_gap),
            "conclusion": "strict outward left-outer-root quartic contact",
        },
        "corridor": {
            "parameters": {
                f"y_{index}": "99/100" for index in range(6, 15)
            },
            "extension_rows": extension_rows,
        },
        "contractions": {
            "count": len(contractions),
            "range": "x_2 through x_14",
            "sequence_sha256": sequence_sha256(contractions),
            "last": decimal_summary(x[14]),
        },
        "scalar_corridors": {
            "strict_increase": len(increase_margins),
            "strict_upper": len(upper_margins),
            "strict_point_walls": len(wall_margins),
            "strict_scaled_defect_steps": len(scaled_margins),
            "strict_cubic_frontiers": len(cubic_margins),
            "reciprocal_defect_steps": len(reciprocal_rows),
            "increase_sha256": sequence_sha256(increase_margins),
            "upper_sha256": sequence_sha256(upper_margins),
            "wall_sha256": sequence_sha256(wall_margins),
            "scaled_sha256": sequence_sha256(scaled_margins),
            "cubic_sha256": sequence_sha256(cubic_margins),
            "reciprocal_rows": reciprocal_rows,
        },
        "signed_layers": {
            "positive_order3_gaps": len(gap_values),
            "positive_order4_margins": len(order4_margins),
            "order3_gap_sha256": sequence_sha256(gap_values),
            "order4_margin_sha256": sequence_sha256(order4_margins),
            "arbitrary_column_minors": minor_summaries,
        },
        "coefficients": {
            "count": len(coefficients),
            "range": "A_1 through A_15",
            "sequence_sha256": sequence_sha256(coefficients),
        },
        "compatibility": {
            "delta_14": decimal_summary(delta14),
            "delta_15": decimal_summary(delta15),
            "length14_survivors": 1,
            "uniform_all_contact_obstructions": 0,
        },
        "adjacent_degree": {
            "normalized_quintic": (
                "1+5*w+10*x_2*w^2+10*x_2^2*x_3*w^3+"
                "5*x_2^3*x_3^2*x_4*w^4+"
                "x_2^4*x_3^3*x_4^2*x_5*w^5"
            ),
            "coefficient_sequence_sha256": sequence_sha256(
                quintic_coefficients
            ),
            "double_root": str(double_root),
            "value_at_double_root_exact": str(quintic_value),
            "value_at_double_root": decimal_summary(quintic_value),
            "value_factorization": (
                "P_5(-1/a)=-2*p^2*(x_5-U)/(a^2*B)"
            ),
            "discriminant_exact": str(quintic_discriminant),
            "discriminant": decimal_summary(quintic_discriminant),
            "distinct_real_roots": 3,
            "nonreal_conjugate_pairs": 1,
            "hyperbolic": False,
            "conclusion": (
                "The exact finite signed-Hankel survivor is excluded by "
                "adjacent degree-five Jensen hyperbolicity."
            ),
        },
    }


def rows(exact: dict[str, object]) -> list[GateRow]:
    return [
        GateRow(
            "qocl14_01_contact",
            "exact_witness",
            "proved_exact",
            "The rational chart parameters define a strict left-outer quartic double-root contact.",
            "delta=13/200, q=1/2",
            "One constructed contact only; not Xi.",
            exact["contact"],
        ),
        GateRow(
            "qocl14_02_outward",
            "exact_witness",
            "proved_exact",
            "The contact violates the inward threshold by an exact positive amount.",
            "x_5-U=delta^2/100000>0",
            "Finite contact witness only; no Xi contact is asserted.",
        ),
        GateRow(
            "qocl14_03_tail",
            "exact_witness",
            "proved_exact",
            "Nine rational corridor parameters construct x_6 through x_14.",
            "y_6=...=y_14=99/100",
            "One selected tail only.",
            exact["corridor"],
        ),
        GateRow(
            "qocl14_04_scalar",
            "exact_certificate",
            "proved_exact",
            "All contractions through x_14 are increasing, below one, above their pointwise walls, and have increasing scaled defects.",
            "12 increase, 13 upper/wall, 12 scaled margins",
            "Finite scalar certificate only.",
            exact["scalar_corridors"],
        ),
        GateRow(
            "qocl14_05_shape",
            "exact_certificate",
            "proved_exact",
            "Every adjacent cubic and reciprocal-defect condition is strict through x_14.",
            "12 cubic and 12 reciprocal steps",
            "Finite local shape conditions only.",
        ),
        GateRow(
            "qocl14_06_contiguous",
            "exact_certificate",
            "proved_exact",
            "Every visible contiguous order-three gap and order-four cap margin is positive.",
            "11 order-three gaps, 9 order-four margins",
            "Finite through A_15 only.",
        ),
        GateRow(
            "qocl14_07_order2",
            "interval_certificate",
            "interval_validated",
            "All supported arbitrary-column order-two signed minors on A_1 through A_15 are positive.",
            "455/455",
            "512-bit Arb enclosures of exact rational inputs.",
            exact["signed_layers"]["arbitrary_column_minors"]["2"],
        ),
        GateRow(
            "qocl14_08_order3",
            "interval_certificate",
            "interval_validated",
            "All supported arbitrary-column order-three signed minors on A_1 through A_15 are positive.",
            "1001/1001",
            "512-bit Arb enclosures of exact rational inputs.",
            exact["signed_layers"]["arbitrary_column_minors"]["3"],
        ),
        GateRow(
            "qocl14_09_order4",
            "interval_certificate",
            "interval_validated",
            "All supported arbitrary-column order-four signed minors on A_1 through A_15 are positive.",
            "1287/1287",
            "512-bit Arb enclosures of exact rational inputs.",
            exact["signed_layers"]["arbitrary_column_minors"]["4"],
        ),
        GateRow(
            "qocl14_10_higher_orders",
            "interval_certificate",
            "interval_validated",
            "All 1300 supported arbitrary-column signed minors of orders five through eight on A_1 through A_15 are positive.",
            "924+330+45+1=1300",
            "512-bit Arb enclosures of exact rational inputs; finite prefix only.",
            {
                order: exact["signed_layers"]["arbitrary_column_minors"][order]
                for order in ("5", "6", "7", "8")
            },
        ),
        GateRow(
            "qocl14_11_delta14",
            "exact_witness",
            "proved_exact",
            "The contact-tail compatibility at length 14 is strictly positive.",
            "Delta_14>0",
            "Rejects a uniform all-contact length-14 obstruction under these finite hypotheses.",
            exact["compatibility"]["delta_14"],
        ),
        GateRow(
            "qocl14_12_delta15",
            "exact_witness",
            "proved_exact",
            "After constructing x_14, the next compatibility Delta_15 also remains positive.",
            "Delta_15>0",
            "No x_15 or infinite extension theorem is claimed.",
            exact["compatibility"]["delta_15"],
        ),
        GateRow(
            "qocl14_13_quintic_contact_value",
            "exact_certificate",
            "proved_exact",
            "The normalized adjacent quintic is nonzero at the quartic double root, with sign fixed by the strict outward slack.",
            "P_5(-1/a)=-2*p^2*(x_5-U)/(a^2*B)<0",
            "Adjacent degree five at this one finite contact only.",
            exact["adjacent_degree"]["value_at_double_root"],
        ),
        GateRow(
            "qocl14_14_quintic_nonhyperbolicity",
            "exact_certificate",
            "proved_exact",
            "The adjacent normalized quintic has negative exact discriminant and hence exactly one nonreal conjugate pair.",
            "Disc(P_5)<0; 3 real roots and 1 nonreal pair",
            "This excludes the survivor by degree-five Jensen hyperbolicity, not by a finite signed-Hankel minor.",
            exact["adjacent_degree"],
        ),
        GateRow(
            "qocl14_15_scope",
            "scope_gate",
            "proved_exact",
            "The one-contact theorem cannot be promoted over every outer contact using only the recorded finite scalar and signed layers, while adjacent-degree hyperbolicity does exclude this witness.",
            "1 exact survivor; 0 uniform all-contact obstructions; 1 nonhyperbolic adjacent quintic",
            "The witness is finite, not Xi, not an all-order sign-regular sequence, and not a Newman trajectory; proving the degree-coupled condition for Xi remains open.",
        ),
    ]


def build_payload() -> dict[str, object]:
    flint.ctx.prec = PRECISION_BITS
    lower_process_priority()
    exact = build_exact()
    gate_rows = rows(exact)
    return {
        "kind": (
            "jensen_window_pf_quartic_outer_contact_"
            "length14_survivor_gate"
        ),
        "date": "2026-07-25",
        "status": "exact finite outer-contact length-14 survivor",
        "proof_boundary": (
            "This exact rational contact and tail survive through x_14 and "
            "all supported arbitrary-column signed minors through the maximum "
            "possible order eight on A_1 through A_15, yet its adjacent "
            "normalized quintic has a negative exact discriminant. It "
            "rejects a uniform all-contact length-14 obstruction from these "
            "finite hypotheses while identifying degree coupling as a "
            "strictly stronger separator. It is not Xi, an infinite "
            "sequence, a Newman trajectory, PF-infinity, Lambda<=0, RH, or "
            "a Clay-prize result."
        ),
        "precision_bits": PRECISION_BITS,
        "exact": exact,
        "rows": [asdict(row) for row in gate_rows],
        "summary": {
            "rows": len(gate_rows),
            "contact_variables": 3,
            "corridor_parameters": 9,
            "contractions": 13,
            "coefficients": 15,
            "order2_minors": 455,
            "order3_minors": 1001,
            "order4_minors": 1287,
            "order5_minors": 924,
            "order6_minors": 330,
            "order7_minors": 45,
            "order8_minors": 1,
            "all_supported_minors": 4043,
            "maximum_supported_order": 8,
            "positive_delta14": 1,
            "positive_delta15": 1,
            "negative_quintic_value_at_double_root": 1,
            "negative_quintic_discriminant": 1,
            "nonreal_quintic_pairs": 1,
            "adjacent_quintic_hyperbolic": 0,
            "uniform_all_contact_obstructions": 0,
        },
    }


def render_note(payload: dict[str, object]) -> str:
    exact = payload["exact"]
    minors = exact["signed_layers"]["arbitrary_column_minors"]
    return "\n".join(
        [
            "# Quartic Outer-Contact Length-14 Survivor Gate",
            "",
            "Date: 2026-07-25",
            "",
            "Status: exact finite outer-contact length-14 survivor; not Xi,",
            "not a uniform theorem, and not a proof of `Lambda<=0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_quartic_outer_contact_length14_survivor_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_quartic_outer_contact_length14_survivor_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_contact_length14_survivor_gate.py",
            "```",
            "",
            "## Contact",
            "",
            "Use the exact normal-form parameters",
            "",
            "```text",
            "delta=13/200, q=1/2,",
            "s=T-1/100000,",
            f"x_5-U={exact['contact']['x5_minus_U']}>0.",
            "```",
            "",
            "Thus this is a strict outward left-outer-root quartic contact.",
            "Set",
            "",
            "```text",
            "y_6=y_7=...=y_14=99/100.",
            "```",
            "",
            "The exact corridor recurrence constructs increasing contractions",
            "through `x_14`. All thirteen contractions are below one and above",
            "their pointwise walls; all twelve scaled-defect, cubic, and",
            "reciprocal-defect steps are strict.",
            "",
            "## Signed Layers",
            "",
            "Exact rational recurrence proves eleven positive contiguous",
            "order-three gaps and nine positive order-four cap margins. At",
            "512-bit Arb precision, direct arbitrary-column enumeration on",
            "`A_1,...,A_15` gives",
            "",
            "```text",
            f"order two:  {minors['2']['strict_positive']}/455 positive",
            f"order three: {minors['3']['strict_positive']}/1001 positive",
            f"order four: {minors['4']['strict_positive']}/1287 positive",
            f"order five: {minors['5']['strict_positive']}/924 positive",
            f"order six:  {minors['6']['strict_positive']}/330 positive",
            f"order seven: {minors['7']['strict_positive']}/45 positive",
            f"order eight: {minors['8']['strict_positive']}/1 positive",
            "```",
            "",
            "Every Arb lower endpoint is strictly positive and is independently",
            "recomputed from the exact rational inputs by the checker.",
            "",
            "## Compatibility",
            "",
            "The terminal exact signs are",
            "",
            "```text",
            f"Delta_14={exact['compatibility']['delta_14']['decimal_ball']}>0,",
            f"Delta_15={exact['compatibility']['delta_15']['decimal_ball']}>0.",
            "```",
            "",
            "Therefore the one-contact interval theorem cannot be promoted to",
            "a uniform all-contact length-14 obstruction under only these",
            "finite scalar and signed-Hankel hypotheses.",
            "",
            "## Adjacent-Degree Separation",
            "",
            "The same exact coefficients define the normalized adjacent",
            "quintic",
            "",
            "```text",
            "P_5(w)=1+5*w+10*x_2*w^2+10*x_2^2*x_3*w^3",
            "       +5*x_2^3*x_3^2*x_4*w^4",
            "       +x_2^4*x_3^3*x_4^2*x_5*w^5.",
            "```",
            "",
            "At the quartic double root, exact arithmetic gives",
            "",
            "```text",
            "P_5(-1/a)=-2*p^2*(x_5-U)/(a^2*B)",
            f"           ={exact['adjacent_degree']['value_at_double_root_exact']}<0,",
            f"Disc(P_5)={exact['adjacent_degree']['discriminant_exact']}<0.",
            "```",
            "",
            "A real quintic with nonzero negative discriminant has exactly",
            "three distinct real roots and one nonreal conjugate pair. Thus",
            "this survivor is not degree-five Jensen hyperbolic. The finite",
            "signed-Hankel conditions tested here are therefore strictly",
            "weaker than the adjacent-degree condition needed at contact.",
            "",
            "This is a finite countermodel gate. It is not the Xi coefficient",
            "sequence, an all-order sign-regular sequence, or a Newman",
            "trajectory. Adjacent-degree closure excludes this witness;",
            "global far-column/all-order structure, theta arithmetic, and",
            "heat compatibility remain candidate ways to establish such",
            "closure for Xi. No",
            "PF-infinity, `Lambda<=0`, RH, or Clay-prize conclusion follows.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote quartic outer-contact length-14 survivor gate: "
        "15 rows, 3 contact variables, 9 corridor parameters, "
        "13 contractions, 15 coefficients, 4043 positive supported "
        "order-2 through order-8 minors, 1 positive Delta_14, "
        "1 positive Delta_15, 1 nonhyperbolic adjacent quintic, "
        "0 uniform all-contact obstructions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
