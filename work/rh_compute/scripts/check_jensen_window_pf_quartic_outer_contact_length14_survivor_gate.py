#!/usr/bin/env python3
"""Check the exact outer-contact length-14 survivor gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_contact_length14_survivor_gate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_quartic_outer_contact_length14_survivor_gate.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_quartic_outer_contact_length14_survivor_gate.py"
)


def load_builder():
    spec = importlib.util.spec_from_file_location("contact_survivor14", BUILDER)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load contact survivor builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, default=RESULT)
    parser.add_argument("--note", type=Path, default=NOTE)
    args = parser.parse_args()
    issues: list[str] = []
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    builder = load_builder()
    builder.flint.ctx.prec = builder.PRECISION_BITS
    builder.lower_process_priority()
    rebuilt = builder.build_payload()
    if payload != rebuilt:
        issues.append("stored payload differs from exact reconstruction")
    if payload.get("kind") != (
        "jensen_window_pf_quartic_outer_contact_"
        "length14_survivor_gate"
    ):
        issues.append("artifact kind changed")
    if payload.get("status") != (
        "exact finite outer-contact length-14 survivor"
    ):
        issues.append("artifact status changed")
    if payload.get("precision_bits") != 512:
        issues.append("Arb precision changed")
    if len(payload.get("rows", [])) != 15:
        issues.append("expected fifteen gate rows")

    exact = payload.get("exact", {})
    contact = exact.get("contact", {})
    delta = Fraction(contact.get("delta", "0"))
    q = Fraction(contact.get("q", "0"))
    a = 1 - delta
    curvature = 4 * q * delta**2
    p = 4 * a - 3 * a**2 + curvature
    simple_discriminant = (4 - 2 * a) ** 2 - 4 * p
    root_A = -3 * a**2 + 8 * a + p
    root_B = -a**2 + 2 * a + p
    threshold = (
        root_B * (3 * a**2 - 5 * a + 5 * p) / (6 * p**2)
    )
    outward_gap = Fraction(contact.get("x5_minus_U", "0"))
    if not (
        delta == Fraction(13, 200)
        and q == Fraction(1, 2)
        and curvature > 0
        and simple_discriminant > 0
        and Fraction(contact.get("threshold_U", "0")) == threshold
        and outward_gap == Fraction(169, 4000000000)
    ):
        issues.append("independent contact reconstruction failed")

    x2 = root_A / 6
    x3 = 18 * a * root_B / root_A**2
    x4 = 2 * p * root_A / (3 * root_B**2)
    x5 = threshold + outward_gap
    w = sp.symbols("w")

    def rational(value: Fraction) -> sp.Rational:
        return sp.Rational(value.numerator, value.denominator)

    quintic = (
        1
        + 5 * w
        + 10 * rational(x2) * w**2
        + 10 * rational(x2**2 * x3) * w**3
        + 5 * rational(x2**3 * x3**2 * x4) * w**4
        + rational(x2**4 * x3**3 * x4**2 * x5) * w**5
    )
    quintic_value_sympy = sp.factor(quintic.subs(w, -1 / rational(a)))
    quintic_discriminant_sympy = sp.factor(sp.discriminant(quintic, w))
    quintic_value = Fraction(
        int(sp.numer(quintic_value_sympy)),
        int(sp.denom(quintic_value_sympy)),
    )
    quintic_discriminant = Fraction(
        int(sp.numer(quintic_discriminant_sympy)),
        int(sp.denom(quintic_discriminant_sympy)),
    )
    expected_quintic_value = (
        -2 * p**2 * outward_gap / (a**2 * root_B)
    )
    adjacent_degree = exact.get("adjacent_degree", {})
    if not (
        quintic_value == expected_quintic_value
        and quintic_value < 0
        and quintic_discriminant < 0
        and adjacent_degree.get("value_at_double_root_exact")
        == str(quintic_value)
        and adjacent_degree.get("discriminant_exact")
        == str(quintic_discriminant)
        and adjacent_degree.get("value_at_double_root", {}).get("sign") == -1
        and adjacent_degree.get("discriminant", {}).get("sign") == -1
        and adjacent_degree.get("distinct_real_roots") == 3
        and adjacent_degree.get("nonreal_conjugate_pairs") == 1
        and adjacent_degree.get("hyperbolic") is False
    ):
        issues.append("independent adjacent-quintic separation failed")

    contractions = exact.get("contractions", {})
    coefficients = exact.get("coefficients", {})
    if contractions.get("count") != 13:
        issues.append("expected x_2 through x_14")
    if coefficients.get("count") != 15:
        issues.append("expected A_1 through A_15")
    compatibility = exact.get("compatibility", {})
    if compatibility.get("delta_14", {}).get("sign") != 1:
        issues.append("Delta_14 is not positive")
    if compatibility.get("delta_15", {}).get("sign") != 1:
        issues.append("Delta_15 is not positive")

    minors = exact.get("signed_layers", {}).get(
        "arbitrary_column_minors", {}
    )
    expected_counts = {
        "2": 455,
        "3": 1001,
        "4": 1287,
        "5": 924,
        "6": 330,
        "7": 45,
        "8": 1,
    }
    for order, count in expected_counts.items():
        row = minors.get(order, {})
        if row.get("precision_bits") != 512:
            issues.append(f"order-{order} precision changed")
        if row.get("count") != count or row.get("strict_positive") != count:
            issues.append(f"order-{order} minor count or sign changed")

    if payload.get("summary") != {
        "rows": 15,
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
    }:
        issues.append("summary counts changed")

    note = args.note.read_text(encoding="utf-8")
    for marker in (
        "exact finite outer-contact length-14 survivor",
        "delta=13/200, q=1/2",
        "s=T-1/100000",
        "y_6=y_7=...=y_14=99/100",
        "order two:  455/455 positive",
        "order three: 1001/1001 positive",
        "order four: 1287/1287 positive",
        "order five: 924/924 positive",
        "order six:  330/330 positive",
        "order seven: 45/45 positive",
        "order eight: 1/1 positive",
        "Delta_14=",
        "Delta_15=",
        "P_5(-1/a)=-2*p^2*(x_5-U)/(a^2*B)",
        "Disc(P_5)=",
        "not degree-five Jensen hyperbolic",
        "weaker than the adjacent-degree condition",
        "uniform all-contact length-14 obstruction",
        "not the Xi coefficient",
        "`Lambda<=0`",
        "RH",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated quartic outer-contact length-14 survivor gate: "
        "15 rows, 0 issues, 3 contact variables, 9 corridor parameters, "
        "13 contractions, 15 coefficients, 4043 positive supported "
        "order-two through order-eight minors, "
        "1 positive Delta_14, 1 positive Delta_15, "
        "1 nonhyperbolic adjacent quintic, "
        "0 uniform all-contact obstructions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
