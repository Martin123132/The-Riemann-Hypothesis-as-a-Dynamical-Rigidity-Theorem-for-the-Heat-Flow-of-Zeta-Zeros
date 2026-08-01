#!/usr/bin/env python3
"""Build an exact alternate length-13 survivor from the outer quartic contact."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys

import sympy as sp


sys.set_int_max_str_digits(0)

REPO_ROOT = Path(__file__).resolve().parents[3]
PARENT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.json"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.md"
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


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational_summary(value: sp.Expr, digits: int = 22) -> dict[str, object]:
    value = sp.cancel(value)
    numerator, denominator = sp.fraction(value)
    encoded = str(value).encode("ascii")
    return {
        "decimal": str(sp.N(value, digits)),
        "numerator_digits": len(str(abs(int(numerator)))),
        "denominator_digits": len(str(abs(int(denominator)))),
        "sha256": hashlib.sha256(encoded).hexdigest(),
    }


def rational_sequence_sha256(values: list[sp.Expr]) -> str:
    encoded = "\n".join(str(sp.cancel(value)) for value in values).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def signed_minor_summary(
    coefficients: list[sp.Expr], order: int, base_index: int = 1
) -> dict[str, object]:
    """Enumerate and hash every supported arbitrary-column reshaped minor."""
    last = len(coefficients) - 1
    signature = (-1) ** (order * (order - 1) // 2)
    digest = hashlib.sha256()
    count = 0
    minimum_log10 = math.inf
    minimum_location: dict[str, object] | None = None
    largest_numerator_digits = 0
    largest_denominator_digits = 0
    for shift in range(len(coefficients)):
        max_column = last - shift - (order - 1)
        if max_column < order - 1:
            continue
        for columns in itertools.combinations(range(max_column + 1), order):
            determinant = sp.cancel(
                sp.det(
                    sp.Matrix(
                        [
                            [
                                coefficients[shift + row + column]
                                for column in columns
                            ]
                            for row in range(order)
                        ]
                    )
                )
            )
            signed = sp.cancel(signature * determinant)
            if signed <= 0:
                raise RuntimeError(
                    f"order-{order} signed minor failed at "
                    f"shift {shift + base_index}, columns {columns}"
                )
            numerator, denominator = sp.fraction(signed)
            numerator_text = str(abs(int(numerator)))
            denominator_text = str(abs(int(denominator)))
            largest_numerator_digits = max(
                largest_numerator_digits, len(numerator_text)
            )
            largest_denominator_digits = max(
                largest_denominator_digits, len(denominator_text)
            )
            leading_numerator = int(numerator_text[:16])
            leading_denominator = int(denominator_text[:16])
            logarithm = (
                len(numerator_text)
                - len(denominator_text)
                + math.log10(leading_numerator / leading_denominator)
            )
            if logarithm < minimum_log10:
                minimum_log10 = logarithm
                minimum_location = {
                    "shift": shift + base_index,
                    "columns": list(columns),
                }
            digest.update(
                (
                    f"{shift + base_index}|"
                    f"{','.join(str(column) for column in columns)}|"
                    f"{signed}\n"
                ).encode("ascii")
            )
            count += 1
    return {
        "count": count,
        "strict_positive": count,
        "sha256": digest.hexdigest(),
        "minimum_log10_approx": f"{minimum_log10:.12f}",
        "minimum_location": minimum_location,
        "largest_numerator_digits": largest_numerator_digits,
        "largest_denominator_digits": largest_denominator_digits,
    }


def build_exact(parent: dict) -> dict[str, object]:
    fixed = {
        2: sp.Rational(24154639, 25000000),
        3: sp.Rational(567331181410000, 583446585220321),
        4: sp.Rational(909681023616163, 930852655574884),
        5: sp.Rational(2453, 2500),
    }
    stored = {
        int(name.split("_")[1]): sp.Rational(value)
        for name, value in parent["exact"]["contractions"].items()
        if int(name.split("_")[1]) <= 5
    }
    if stored != fixed:
        raise RuntimeError("parent outer-contact coordinates changed")

    parameters = {
        6: sp.Rational(37, 100),
        7: sp.Rational(69, 70),
        8: sp.Rational(13, 15),
        9: sp.Rational(7, 10),
        10: sp.Rational(7, 8),
        11: sp.Rational(9, 20),
        12: sp.Rational(1, 9),
        13: sp.Rational(1, 2),
    }
    x: dict[int, sp.Expr] = dict(fixed)

    def defect(index: int) -> sp.Expr:
        return sp.cancel(1 - x[index])

    def gap(index: int) -> sp.Expr:
        return sp.cancel(
            defect(index + 2) ** 2
            - x[index + 2] ** 2
            * defect(index + 1)
            * defect(index + 3)
        )

    def cap(index: int) -> sp.Expr:
        return sp.cancel(
            gap(index - 4) ** 2
            / (x[index - 2] ** 3 * gap(index - 5))
        )

    extension_identities: list[dict[str, object]] = []
    for index, parameter in parameters.items():
        current_cap = cap(index)
        target_gap = sp.cancel(parameter * current_cap)
        current_defect = sp.cancel(
            (defect(index - 1) ** 2 - target_gap)
            / (x[index - 1] ** 2 * defect(index - 2))
        )
        x[index] = sp.cancel(1 - current_defect)
        if sp.cancel(gap(index - 3) - target_gap) != 0:
            raise RuntimeError(f"corridor identity failed at x_{index}")
        extension_identities.append(
            {
                "index": index,
                "parameter": str(parameter),
                "contraction": rational_summary(x[index]),
                "gap_sha256": hashlib.sha256(
                    str(gap(index - 3)).encode("ascii")
                ).hexdigest(),
                "cap_sha256": hashlib.sha256(
                    str(current_cap).encode("ascii")
                ).hexdigest(),
            }
        )

    contractions = [x[index] for index in range(2, 14)]
    increase_margins = [
        sp.cancel(contractions[index + 1] - contractions[index])
        for index in range(len(contractions) - 1)
    ]
    upper_margins = [sp.cancel(1 - value) for value in contractions]
    wall_margins = [
        sp.cancel(x[index] - sp.Rational(2 * index - 1, 2 * index + 1))
        for index in range(2, 14)
    ]
    scaled_defects = [
        sp.cancel(sp.Rational(2 * index + 1, 2) * (1 - x[index]))
        for index in range(2, 14)
    ]
    scaled_margins = [
        sp.cancel(scaled_defects[index + 1] - scaled_defects[index])
        for index in range(len(scaled_defects) - 1)
    ]
    if not all(
        value > 0
        for value in increase_margins
        + upper_margins
        + wall_margins
        + scaled_margins
    ):
        raise RuntimeError("a scalar contraction or scaled-defect wall failed")

    cubic_margins: list[sp.Expr] = []
    reciprocal_certificates: list[dict[str, object]] = []
    for left, right in zip(contractions, contractions[1:]):
        cubic = sp.cancel(
            left**2 * right**2 - 6 * left * right + 4 * left + 4 * right - 3
        )
        if cubic >= 0:
            raise RuntimeError("an adjacent cubic frontier failed")
        cubic_margins.append(cubic)
        left_defect = sp.cancel(1 - left)
        right_defect = sp.cancel(1 - right)
        comparison = sp.cancel(1 / right_defect - 1 / left_defect - 1)
        if comparison <= 0:
            reciprocal_certificates.append(
                {
                    "branch": "automatic_nonpositive_comparison",
                    "comparison_sha256": hashlib.sha256(
                        str(comparison).encode("ascii")
                    ).hexdigest(),
                }
            )
            continue
        squared_margin = sp.cancel(4 - left_defect * comparison**2)
        if squared_margin <= 0:
            raise RuntimeError("a reciprocal-defect increment bound failed")
        reciprocal_certificates.append(
            {
                "branch": "positive_comparison_one_squaring",
                "comparison_sha256": hashlib.sha256(
                    str(comparison).encode("ascii")
                ).hexdigest(),
                "squared_margin_sha256": hashlib.sha256(
                    str(squared_margin).encode("ascii")
                ).hexdigest(),
            }
        )

    gap_values = [gap(index) for index in range(1, 11)]
    h4_margins = [
        sp.cancel(cap(index) - gap(index - 3)) for index in range(6, 14)
    ]
    if not all(value > 0 for value in gap_values + h4_margins):
        raise RuntimeError("a contiguous signed-Hankel layer failed")

    def compatibility(index: int) -> sp.Expr:
        repeated = sp.cancel(
            defect(index - 1) ** 2
            - x[index - 1] ** 2
            * defect(index - 2)
            * defect(index - 1)
        )
        return sp.cancel(cap(index) - repeated)

    delta_13 = compatibility(13)
    delta_14 = compatibility(14)
    if not (delta_13 > 0 and delta_14 < 0):
        raise RuntimeError("terminal compatibility signs changed")

    coefficients = [sp.Integer(1), sp.Integer(1)]
    ratio = sp.Integer(1)
    for contraction in contractions:
        ratio = sp.cancel(ratio * contraction)
        coefficients.append(sp.cancel(coefficients[-1] * ratio))
    if len(coefficients) != 14:
        raise RuntimeError("expected A_1 through A_14")

    minor_summaries = {
        str(order): signed_minor_summary(coefficients, order)
        for order in (2, 3, 4)
    }
    expected_counts = {"2": 364, "3": 715, "4": 792}
    if {
        order: summary["count"] for order, summary in minor_summaries.items()
    } != expected_counts:
        raise RuntimeError("arbitrary-column minor counts changed")

    scout = {
        "status": "numerical_evidence_only",
        "uniform_samples": 800000,
        "uniform_seed": 202607250714,
        "scalar_feasible_through_x13": 6322,
        "best_uniform_delta14": "-1.3662564179373033e-06",
        "optimizer_seeds": [20260725, 20260726, 20260727],
        "best_optimizer_delta14": "-1.2515638495620113e-06",
        "interpretation": (
            "No positive Delta_14 was found, but this bounded search is not "
            "a proof of a uniform length-fourteen obstruction."
        ),
    }
    return {
        "definitions": {
            "defect": "d_j=1-x_j",
            "order3_gap": "G_j=d_(j+2)^2-x_(j+2)^2*d_(j+1)*d_(j+3)",
            "order4_cap": "C_k=G_(k-4)^2/(x_(k-2)^3*G_(k-5))",
            "corridor": "G_(k-3)=y_k*C_k",
            "compatibility": (
                "Delta_k=C_k-[d_(k-1)^2-"
                "x_(k-1)^2*d_(k-2)*d_(k-1)]"
            ),
        },
        "fixed_outer_contact": {
            f"x_{index}": str(value) for index, value in fixed.items()
        },
        "corridor_parameters": {
            f"y_{index}": str(value) for index, value in parameters.items()
        },
        "extension_identities": extension_identities,
        "contractions": {
            f"x_{index}": rational_summary(x[index])
            for index in range(2, 14)
        },
        "contraction_sequence_sha256": rational_sequence_sha256(contractions),
        "coefficient_sequence_sha256": rational_sequence_sha256(coefficients),
        "scalar_corridors": {
            "strict_increase_count": len(increase_margins),
            "strict_upper_count": len(upper_margins),
            "strict_point_wall_count": len(wall_margins),
            "strict_scaled_defect_increase_count": len(scaled_margins),
            "strict_cubic_frontier_count": len(cubic_margins),
            "reciprocal_certificates": reciprocal_certificates,
            "increase_margins_sha256": rational_sequence_sha256(
                increase_margins
            ),
            "wall_margins_sha256": rational_sequence_sha256(wall_margins),
            "scaled_margins_sha256": rational_sequence_sha256(scaled_margins),
            "cubic_margins_sha256": rational_sequence_sha256(cubic_margins),
        },
        "signed_layers": {
            "positive_order3_gaps": len(gap_values),
            "positive_order4_margins": len(h4_margins),
            "order3_gap_sha256": rational_sequence_sha256(gap_values),
            "order4_margin_sha256": rational_sequence_sha256(h4_margins),
            "arbitrary_column_minors": minor_summaries,
        },
        "terminal_compatibility": {
            "delta_13": rational_summary(delta_13),
            "delta_13_sign": "positive",
            "delta_14": rational_summary(delta_14),
            "delta_14_sign": "negative",
            "conclusion": (
                "The exact alternate tail survives through x_13, while no "
                "increasing x_14 can retain its next order-four sign."
            ),
        },
        "bounded_delta14_scout": scout,
    }


def build_payload() -> dict:
    parent = json.loads(PARENT_RESULT.read_text(encoding="utf-8"))
    if parent.get("kind") != (
        "jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate"
    ):
        raise RuntimeError("unexpected parent artifact kind")
    exact = build_exact(parent)
    minors = exact["signed_layers"]["arbitrary_column_minors"]
    rows = [
        GateRow(
            "qoas13_01_parent",
            "provenance",
            "proved_exact",
            "The construction starts at the same exact outer quartic contact as the order-four nonpromotion gate.",
            "x_2,...,x_5 fixed; u>U",
            "It does not retain that parent's later x_6 through x_9 tail.",
        ),
        GateRow(
            "qoas13_02_tail",
            "construction",
            "proved_exact",
            "Eight rational open-corridor parameters construct a distinct exact tail x_6 through x_13.",
            "G_(k-3)=y_k*C_k, k=6,...,13",
            "Finite rational construction only.",
            {"parameters": exact["corridor_parameters"]},
        ),
        GateRow(
            "qoas13_03_ratios",
            "exact_certificate",
            "proved_exact",
            "All twelve contractions increase strictly, stay below one, and clear their pointwise walls.",
            "(2k-1)/(2k+1)<x_k<x_(k+1)<1",
            "Indices 2 through 13 only.",
        ),
        GateRow(
            "qoas13_04_scaled",
            "exact_certificate",
            "proved_exact",
            "All eleven scaled defects increase strictly.",
            "(2k+3)(1-x_(k+1))>(2k+1)(1-x_k)",
            "Adjacent finite scalar corridor only.",
        ),
        GateRow(
            "qoas13_05_reciprocal",
            "exact_certificate",
            "proved_exact",
            "All eleven reciprocal-square-root defect increments are strictly below one.",
            "1/sqrt(d_(k+1))-1/sqrt(d_k)<1",
            "Exact rational one-squaring certificates where needed.",
        ),
        GateRow(
            "qoas13_06_cubic",
            "exact_certificate",
            "proved_exact",
            "All eleven adjacent cubic Jensen frontiers are strict.",
            "x_k^2*x_(k+1)^2-6*x_k*x_(k+1)+4*x_k+4*x_(k+1)-3<0",
            "Adjacent degree-three tests only.",
        ),
        GateRow(
            "qoas13_07_contiguous",
            "exact_certificate",
            "proved_exact",
            "All ten visible order-three gaps and all eight visible order-four cap margins are strictly positive.",
            "G_1,...,G_10>0; C_k-G_(k-3)>0",
            "Contiguous layers through A_14 only.",
        ),
        GateRow(
            "qoas13_08_order2",
            "exact_certificate",
            "proved_exact",
            "Every supported arbitrary-column order-two reshaped minor has the required strict sign.",
            "(-1)^1 det H_2[I,J]>0",
            "Finite A_1 through A_14 segment only.",
            minors["2"],
        ),
        GateRow(
            "qoas13_09_order3",
            "exact_certificate",
            "proved_exact",
            "Every supported arbitrary-column order-three reshaped minor has the required strict sign.",
            "(-1)^3 det H_3[I,J]>0",
            "Finite A_1 through A_14 segment only.",
            minors["3"],
        ),
        GateRow(
            "qoas13_10_order4",
            "exact_certificate",
            "proved_exact",
            "Every supported arbitrary-column order-four reshaped minor has the required strict sign.",
            "(-1)^6 det H_4[I,J]>0",
            "Finite A_1 through A_14 segment only.",
            minors["4"],
        ),
        GateRow(
            "qoas13_11_survival",
            "scope_gate",
            "proved_exact",
            "The same outer quartic contact has an alternate tail that survives all stated gates through x_13.",
            "Delta_13>0",
            "This refutes only a uniform reading of the earlier fixed-prefix obstruction.",
            exact["terminal_compatibility"]["delta_13"],
        ),
        GateRow(
            "qoas13_12_fixed_tail_death",
            "scope_gate",
            "proved_exact",
            "This particular alternate tail cannot continue to an increasing x_14 while retaining the next order-four sign.",
            "Delta_14<0",
            "Fixed-tail result, not a uniform theorem over all corridor choices.",
            exact["terminal_compatibility"]["delta_14"],
        ),
        GateRow(
            "qoas13_13_scout",
            "numerical_diagnostic",
            "numerical_only",
            "A bounded random and optimization scout found no scalar-admissible positive Delta_14.",
            "best sampled Delta_14<0",
            "Search evidence only; it cannot certify a uniform obstruction.",
            exact["bounded_delta14_scout"],
        ),
        GateRow(
            "qoas13_14_boundary",
            "scope_gate",
            "proved_exact",
            "Finite survival through x_13 neither promotes the quartic threshold nor supplies an infinite sign-regular sequence.",
            "finite survivor != PF-infinity",
            "No Xi-specific inequality, Lambda<=0, RH, or Clay-prize conclusion.",
        ),
    ]
    return {
        "kind": (
            "jensen_window_pf_quartic_outer_branch_"
            "alternate_length13_survivor_gate"
        ),
        "date": "2026-07-25",
        "status": "exact alternate length-thirteen survivor and fixed-tail obstruction",
        "proof_boundary": (
            "This artifact constructs one exact tail from the same outer "
            "quartic contact that survives every stated scalar and finite "
            "signed-Hankel gate through x_13, correcting any uniform reading "
            "of the earlier fixed-prefix length-thirteen obstruction. It "
            "also proves that this selected alternate tail dies at x_14. "
            "It does not prove a uniform finite-tail obstruction, an "
            "all-length extension, PF-infinity, an Xi-specific estimate, "
            "Lambda<=0, or RH."
        ),
        "parent": {
            "path": str(PARENT_RESULT.relative_to(REPO_ROOT)),
            "sha256": file_sha256(PARENT_RESULT),
        },
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "fixed_outer_contacts": 1,
            "corridor_parameters": 8,
            "contractions": 12,
            "strict_adjacent_scalar_steps": 11,
            "positive_order3_gaps": 10,
            "positive_order4_margins": 8,
            "signed_order2_minors": minors["2"]["count"],
            "signed_order3_minors": minors["3"]["count"],
            "signed_order4_minors": minors["4"]["count"],
            "positive_length13_compatibilities": 1,
            "negative_fixed_tail_length14_compatibilities": 1,
            "uniform_length14_theorems": 0,
            "repaired_scope_handoffs": 1,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    terminal = exact["terminal_compatibility"]
    minors = exact["signed_layers"]["arbitrary_column_minors"]
    contractions = exact["contractions"]
    return "\n".join(
        [
            "# Quartic Outer-Branch Alternate Length-13 Survivor Gate",
            "",
            "Date: 2026-07-25",
            "",
            "Status: exact alternate length-thirteen survivor and fixed-tail length-fourteen obstruction; not a proof of a uniform finite-tail theorem, `Lambda <= 0`, or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.py",
            "```",
            "",
            "## Exact Alternate Tail",
            "",
            "Keep only the exact `x_2,...,x_5` outer-contact data from the",
            "quartic nonpromotion gate. Define",
            "",
            "```text",
            exact["definitions"]["order3_gap"],
            exact["definitions"]["order4_cap"],
            "G_(k-3)=y_k*C_k.",
            "```",
            "",
            "Use the rational coordinates",
            "",
            "```text",
            *(
                f"{name}={value}"
                for name, value in exact["corridor_parameters"].items()
            ),
            "```",
            "",
            "They give the following exact rationals, displayed only by",
            "certified decimals and hashes because the last fractions contain",
            "thousands of digits:",
            "",
            "```text",
            *(
                f"{name} ~= {row['decimal']}  sha256={row['sha256']}"
                for name, row in contractions.items()
            ),
            "```",
            "",
            "All contractions increase, remain below one, clear the",
            "pointwise walls, have increasing scaled defects, satisfy the",
            "reciprocal-square-root increment bound, and satisfy every",
            "adjacent cubic frontier.",
            "",
            "## Finite Signed Layers",
            "",
            "Exact direct enumeration over `A_1,...,A_14` gives",
            "",
            "```text",
            f"order 2: {minors['2']['count']} strict signed minors",
            f"order 3: {minors['3']['count']} strict signed minors",
            f"order 4: {minors['4']['count']} strict signed minors",
            "```",
            "",
            "The checker rebuilds the rationals, hashes every determinant,",
            "and verifies all signs rather than trusting floating point.",
            "",
            "## Survival And Death",
            "",
            "For the increasing compatibility margin",
            "",
            "```text",
            "Delta_k=C_k-[d_(k-1)^2-x_(k-1)^2*d_(k-2)*d_(k-1)],",
            f"Delta_13 ~= {terminal['delta_13']['decimal']} > 0,",
            f"Delta_14 ~= {terminal['delta_14']['decimal']} < 0.",
            "```",
            "",
            "Thus the same outer quartic contact does have a different tail",
            "surviving through `x_13`; the earlier length-13 obstruction is",
            "genuinely prefix-specific. This selected alternate tail cannot",
            "reach an increasing `x_14` with the next order-four sign.",
            "",
            "A bounded 800,000-point scout plus three deterministic",
            "optimizations also found no positive `Delta_14`, but that",
            "observation is numerical evidence only. It is not promoted to a",
            "uniform theorem.",
            "",
            "## Corrected Handoff",
            "",
            "The live exact target is now a uniform semialgebraic statement",
            "over all admissible earlier corridor choices, or another",
            "explicit tail that survives the length-14 boundary. Neither",
            "finite outcome supplies PF-infinity, an Xi-specific bridge,",
            "`Lambda <= 0`, or RH.",
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
    print(f"wrote alternate length-13 survivor gate: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
