#!/usr/bin/env python3
"""Check the sparse H0-H23 to H0-H16 Taylor propagation core."""

from __future__ import annotations

from fractions import Fraction
import math
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

from jensen_window_pf_compound_order11_sparse_h0_h23_propagation_core import (  # noqa: E402
    ANCHOR_MAXIMUM_H_ORDER,
    OUTPUT_MAXIMUM_H_ORDER,
    propagated_h0_h16_jet,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


DISPLACEMENTS = (
    Fraction(-1),
    Fraction(-1, 2),
    Fraction(0),
    Fraction(1, 2),
    Fraction(1),
)


def support_row(h24: flint.arb) -> dict:
    return {
        "target_t_left": Fraction(-1),
        "target_t_right": Fraction(1),
        "H": {24: h24},
    }


def exact_polynomial_checks() -> list[str]:
    issues = []
    coefficients = [
        arb_rational(Fraction((-1) ** order, order + 1))
        for order in range(ANCHOR_MAXIMUM_H_ORDER + 1)
    ]
    for displacement in DISPLACEMENTS:
        propagated, diagnostics = propagated_h0_h16_jet(
            displacement,
            Fraction(0),
            coefficients,
            [support_row(flint.arb(0))],
            point_diagnostics={"fixture": "degree-23-polynomial"},
        )
        for derivative in range(OUTPUT_MAXIMUM_H_ORDER + 1):
            expected = flint.arb(0)
            for order in range(derivative, ANCHOR_MAXIMUM_H_ORDER + 1):
                expected += (
                    coefficients[order]
                    * math.comb(order, derivative)
                    * arb_rational(displacement) ** (order - derivative)
                )
            if not bool(propagated[derivative].contains(expected)):
                issues.append(
                    f"polynomial H{derivative} misses at displacement {displacement}"
                )
        propagation = diagnostics.get("sparse_h23_propagation", {})
        if (
            propagation.get("anchor_t") != "0"
            or propagation.get("target_t") != str(displacement)
            or propagation.get("displacement") != str(displacement)
        ):
            issues.append(f"diagnostics changed at displacement {displacement}")
    return issues


def exponential_remainder_checks() -> list[str]:
    issues = []
    coefficients = [
        flint.arb(1) / math.factorial(order)
        for order in range(ANCHOR_MAXIMUM_H_ORDER + 1)
    ]
    h24_wall = flint.arb(3)
    for displacement in DISPLACEMENTS:
        propagated, _ = propagated_h0_h16_jet(
            displacement,
            Fraction(0),
            coefficients,
            [support_row(h24_wall)],
            point_diagnostics={"fixture": "exponential"},
        )
        exponential = arb_rational(displacement).exp()
        for derivative in range(OUTPUT_MAXIMUM_H_ORDER + 1):
            expected = exponential / math.factorial(derivative)
            if not bool(propagated[derivative].contains(expected)):
                issues.append(
                    f"exponential H{derivative} misses at displacement {displacement}"
                )
    return issues


def main() -> int:
    flint.ctx.prec = 256
    issues = exact_polynomial_checks() + exponential_remainder_checks()
    if issues:
        print(f"order-eleven sparse H23 propagation: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated order-eleven sparse H23 propagation: "
        "5 polynomial translations + 5 H24 remainder enclosures, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
