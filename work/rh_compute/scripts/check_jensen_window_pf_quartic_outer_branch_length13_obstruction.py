#!/usr/bin/env python3
"""Check the exact fixed-prefix quartic length-13 obstruction."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_branch_length13_obstruction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_quartic_outer_branch_length13_obstruction.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_quartic_outer_branch_length13_obstruction.py"
)


def load_builder():
    spec = importlib.util.spec_from_file_location("quartic_length13_builder", BUILDER)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load obstruction builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def direct_hankel_reconstruction(parent: dict, builder) -> dict[str, object]:
    """Reconstruct the continuation cube from raw terminal Hankel roots."""
    coefficients = [
        sp.Rational(value)
        for value in parent["exact"]["coefficients_A1_to_A10"]
    ]
    variables = sp.symbols("y_10 y_11 y_12")
    contractions: dict[int, sp.Expr] = {}

    def append_contraction(
        current: list[sp.Expr], contraction: sp.Expr
    ) -> list[sp.Expr]:
        return [
            *current,
            sp.cancel(
                current[-1]
                * (current[-1] / current[-2])
                * contraction
            ),
        ]

    def terminal_determinant(
        current: list[sp.Expr], order: int
    ) -> sp.Expr:
        start = len(current) - (2 * order - 1)
        return sp.cancel(
            sp.det(
                sp.Matrix(
                    [
                        [
                            current[start + row + column]
                            for column in range(order)
                        ]
                        for row in range(order)
                    ]
                )
            )
        )

    def terminal_root(
        current: list[sp.Expr], order: int, variable: sp.Symbol
    ) -> sp.Expr:
        extended = append_contraction(current, variable)
        determinant = terminal_determinant(extended, order)
        at_zero = sp.cancel(determinant.subs(variable, 0))
        slope = sp.cancel(determinant.subs(variable, 1) - at_zero)
        if slope == 0:
            raise RuntimeError("terminal determinant lost affine dependence")
        return sp.cancel(-at_zero / slope)

    for index, parameter in zip((10, 11, 12), variables):
        variable = sp.Symbol(f"direct_x_{index}")
        order3_root = terminal_root(coefficients, 3, variable)
        order4_root = terminal_root(coefficients, 4, variable)
        contraction = sp.cancel(
            order3_root + parameter * (order4_root - order3_root)
        )
        contractions[index] = contraction
        coefficients = append_contraction(coefficients, contraction)

    next_variable = sp.Symbol("direct_x_13")
    next_order4_root = terminal_root(coefficients, 4, next_variable)
    direct_gap = sp.cancel(
        (next_order4_root - contractions[12])
        * contractions[12] ** 2
        * (1 - contractions[11])
    )
    derivative_certificates = {
        str(variable): builder.bernstein_certificate(
            sp.fraction(sp.cancel(sp.diff(direct_gap, variable)))[0],
            variables,
        )
        for variable in variables
    }
    return {
        "corner": str(
            sp.cancel(
                direct_gap.subs({variable: 1 for variable in variables})
            )
        ),
        "derivatives": derivative_certificates,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    issues: list[str] = []
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    if payload.get("kind") != (
        "jensen_window_pf_quartic_outer_branch_length13_obstruction"
    ):
        issues.append("artifact kind changed")
    if payload.get("status") != (
        "exact fixed-prefix length-thirteen obstruction"
    ):
        issues.append("artifact status changed")
    rows = payload.get("rows", [])
    if len(rows) != 10:
        issues.append(f"expected 10 rows, found {len(rows)}")

    builder = load_builder()
    parent = json.loads(builder.PARENT_RESULT.read_text(encoding="utf-8"))
    rebuilt_exact = builder.build_symbolic(parent)
    if payload.get("exact") != rebuilt_exact:
        issues.append("stored exact obstruction data did not reproduce")
    direct = direct_hankel_reconstruction(parent, builder)
    if direct["corner"] != rebuilt_exact["compatibility"]["cube_corner"]:
        issues.append("direct Hankel-root cube corner did not reproduce")
    for name, certificate in direct["derivatives"].items():
        stored = rebuilt_exact["compatibility"]["derivatives"][name]
        for key in (
            "degrees",
            "coefficient_count",
            "positive",
            "negative",
            "zero",
            "sha256",
        ):
            if certificate[key] != stored[key]:
                issues.append(
                    f"direct Hankel-root {name} Bernstein {key} changed"
                )
    parent_hash = hashlib.sha256(builder.PARENT_RESULT.read_bytes()).hexdigest()
    if payload.get("parent", {}).get("sha256") != parent_hash:
        issues.append("parent artifact hash changed")

    expected_summary = {
        "rows": 10,
        "extension_parameters": 3,
        "derivative_certificates": 3,
        "derivative_bernstein_coefficients": 1516,
        "denominator_factor_bernstein_coefficients": 8,
        "negative_cube_corner_maxima": 1,
        "forbidden_all_length_promotions": 1,
        "repaired_handoffs": 1,
    }
    if payload.get("summary") != expected_summary:
        issues.append("summary counts changed")

    compatibility = rebuilt_exact["compatibility"]
    if compatibility.get("denominator_sign_on_open_cube") != -1:
        issues.append("compatibility denominator sign changed")
    derivative_expected = {
        "y_10": ([16, 8, 3], 612, 612, 0, 0, 1),
        "y_11": ([16, 7, 3], 544, 0, 480, 64, -1),
        "y_12": ([14, 7, 2], 360, 0, 216, 144, -1),
    }
    for name, expected in derivative_expected.items():
        certificate = compatibility["derivatives"][name]
        observed = (
            certificate["degrees"],
            certificate["coefficient_count"],
            certificate["positive"],
            certificate["negative"],
            certificate["zero"],
            certificate["denominator_sign_on_open_cube"],
        )
        if observed != expected:
            issues.append(f"{name} derivative certificate changed")
        if certificate.get("derivative_sign_on_open_cube") != "nonnegative":
            issues.append(f"{name} derivative conclusion changed")
    if not compatibility["cube_corner"].startswith("-"):
        issues.append("cube-corner maximum is not negative")

    note = args.note.read_text(encoding="utf-8")
    required_note_markers = [
        "This is",
        "not a proof",
        "G_(k-3)=y_k*C_k",
        "Delta_13=C_13-M_13",
        "all three partial derivatives",
        "Delta_13(1,1,1)=",
        "cannot continue",
        "This does not promote the quartic threshold",
        "uniform finite-tail obstruction",
    ]
    for marker in required_note_markers:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated quartic outer-branch length-13 obstruction: "
        "10 rows, 0 issues, 3 exact corridor parameters, "
        "3 derivative certificates, 1516 derivative Bernstein coefficients, "
        "8 denominator-factor Bernstein coefficients, "
        "1 negative cube-corner maximum, 1 forbidden all-length promotion, "
        "1 repaired handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
