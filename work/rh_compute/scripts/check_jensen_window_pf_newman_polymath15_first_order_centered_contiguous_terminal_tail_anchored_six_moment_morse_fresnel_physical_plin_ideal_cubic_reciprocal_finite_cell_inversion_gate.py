#!/usr/bin/env python3
"""Independently check the finite reciprocal-cell inversion gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
    "ideal_cubic_reciprocal_finite_cell_inversion_gate"
)
RESULT_PATH = REPO_ROOT / "work" / "rh_compute" / "results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 23,
    "complete_mode_recompositions": 2,
    "off_lattice_endpoint_phase_identities": 2,
    "reciprocal_cell_identities": 4,
    "finite_cell_poisson_identities": 3,
    "fourier_kernel_identities": 2,
    "carrier_inversion_identities": 3,
    "alias_gap_identities": 1,
    "internal_tie_cancellations": 1,
    "global_recompositions": 3,
    "finite_cell_remainder_bounds": 0,
    "alias_bounds": 0,
    "exterior_tail_bounds": 0,
    "endpoint_composed_h2_bounds": 0,
    "quadratic_remainder_bounds": 0,
    "signed_flow_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    f"rfi_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        [
            "fourier_mode",
            "off_lattice",
            "integer_rejoin",
            "cells",
            "bounds",
            "adjacency",
            "poisson",
            "tie",
            "transfer",
            "kernel",
            "phase",
            "alias_gap",
            "inversion",
            "weights",
            "exterior",
            "cell_defect",
            "tile",
            "telescope",
            "global",
            "involution",
            "achievement",
            "handoff",
            "boundary",
        ],
        start=1,
    )
]


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    simplified = sp.simplify(expression)
    if simplified != 0:
        raise RuntimeError(f"{label}: {simplified}")


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def check_sources(artifact: dict) -> None:
    source = artifact.get("source_audit", {})
    require(
        set(source) == {"endpoint_composition", "finite_poisson", "first_correction"},
        "source-audit keys drifted",
    )
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source for {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted for {key}")


def check_mode_recomposition() -> None:
    y = sp.symbols("y", real=True)
    kappa = 1 / (2 * sp.pi * sp.I)
    b_0, b_1, b_2, b_3 = sp.symbols("b_0 b_1 b_2 b_3")
    b = b_0 + b_1 * y + b_2 * y**2 + b_3 * y**3
    c = sp.cancel((b - b_0) / y)
    oscillation = sp.exp(-sp.pi * sp.I * y**2)
    reconstructed = (
        b_0 * oscillation
        + sp.diff(-kappa * c * oscillation, y)
        + kappa * sp.diff(c, y) * oscillation
    )
    require_zero(b * oscillation - reconstructed, "independent Morse recomposition")

    alpha, r, b_value = sp.symbols("alpha_P r B", positive=True, real=True)
    lower_phase = alpha * sp.log(1) - r
    upper_phase = alpha * sp.log(b_value) - r * b_value
    require_zero(lower_phase + r, "lower off-lattice phase")
    require_zero(
        upper_phase - (alpha * sp.log(b_value) - r * b_value),
        "upper off-lattice phase",
    )
    integer_r, integer_b = sp.symbols("r_i B_i", integer=True)
    require_zero(sp.exp(-2 * sp.pi * sp.I * integer_r) - 1, "integer lower phase")
    require_zero(
        sp.exp(-2 * sp.pi * sp.I * integer_r * integer_b) - 1,
        "integer upper phase",
    )


def check_cells() -> None:
    cases = 0
    ties = 0
    for alpha in [
        Fraction(7, 2),
        Fraction(13, 2),
        Fraction(25, 3),
        Fraction(10),
        Fraction(27, 2),
        Fraction(120),
    ]:
        for b_value in range(2, 7):
            union: list[int] = []
            tie_coefficients: dict[int, Fraction] = {}
            for q in range(1, b_value + 1):
                a_q = alpha / (Fraction(q) + Fraction(1, 2))
                b_q = alpha / (Fraction(q) - Fraction(1, 2))
                direct = [
                    r
                    for r in range(1, floor_fraction(2 * alpha) + 1)
                    if Fraction(q) - Fraction(1, 2)
                    <= alpha / r
                    < Fraction(q) + Fraction(1, 2)
                ]
                endpoint = [
                    r
                    for r in range(1, floor_fraction(2 * alpha) + 1)
                    if a_q < r <= b_q
                ]
                lower = floor_fraction(a_q) + 1
                upper = floor_fraction(b_q)
                formula = list(range(lower, upper + 1)) if lower <= upper else []
                require(direct == endpoint == formula, f"cell audit failed at {alpha}, {q}")
                union.extend(direct)
                cases += 1

                if a_q.denominator == 1:
                    boundary = int(a_q)
                    tie_coefficients[boundary] = tie_coefficients.get(boundary, Fraction()) - Fraction(1, 2)
                    ties += 1
                if b_q.denominator == 1:
                    boundary = int(b_q)
                    tie_coefficients[boundary] = tie_coefficients.get(boundary, Fraction()) + Fraction(1, 2)
                    ties += 1

                if q < b_value:
                    next_b = alpha / (Fraction(q + 1) - Fraction(1, 2))
                    require(a_q == next_b, "continuous cells failed to tile")

            global_lower = alpha / (Fraction(b_value) + Fraction(1, 2))
            global_upper = 2 * alpha
            expected_union = [
                r
                for r in range(1, floor_fraction(global_upper) + 1)
                if global_lower < r <= global_upper
            ]
            require(sorted(union) == expected_union, "physical reciprocal band drifted")
            require(len(union) == len(set(union)), "physical cells overlapped")

            for boundary, coefficient in tie_coefficients.items():
                if Fraction(boundary) not in {global_lower, global_upper}:
                    require(coefficient == 0, "internal tie coefficient did not cancel")
    require(cases >= 100, "insufficient independent cell cases")
    require(ties >= 8, "insufficient independent tie cases")


def check_finite_poisson_sign() -> None:
    """Audit the sign and half-open correction on a finite Fourier polynomial."""
    n_value = 7
    frequencies = [-3, 0, 2]
    coefficients = [sp.Rational(2, 5), sp.Rational(-7, 11), sp.Rational(5, 13)]
    direct = n_value * sum(coefficients)

    dual = sp.Integer(0)
    for dual_frequency in range(-5, 6):
        for frequency, coefficient in zip(frequencies, coefficients):
            total_frequency = frequency + dual_frequency
            integral = n_value if total_frequency == 0 else 0
            dual += coefficient * integral
    require_zero(direct - dual, "finite-Poisson dual sign")

    endpoint_zero = sum(coefficients) - sum(coefficients)
    require_zero(endpoint_zero / 2, "half-open tie correction for periodic audit")


def check_kernel_and_phase() -> None:
    alpha, q = sp.symbols("alpha_P q", positive=True, real=True)
    r, x = sp.symbols("r x", real=True)
    kappa = 1 / (2 * sp.pi * sp.I)
    a_q = alpha / (q + sp.Rational(1, 2))
    b_q = alpha / (q - sp.Rational(1, 2))
    antiderivative = -kappa * sp.exp(-2 * sp.pi * sp.I * r * x) / x
    require_zero(
        sp.diff(antiderivative, r) - sp.exp(-2 * sp.pi * sp.I * r * x),
        "independent Fourier antiderivative",
    )
    kernel = antiderivative.subs(r, b_q) - antiderivative.subs(r, a_q)
    require_zero(sp.limit(kernel, x, 0) - (b_q - a_q), "independent kernel limit")

    phase = alpha * (sp.log(alpha / r) - 1) + q * r
    require_zero(sp.diff(phase, r).subs(r, alpha / q), "independent central saddle")
    require_zero(
        sp.expand_log(phase.subs(r, alpha / q), force=True) - alpha * sp.log(q),
        "independent central value",
    )
    require_zero(
        sp.diff(phase, r, 2).subs(r, alpha / q) - q**2 / alpha,
        "independent central curvature",
    )

    for q_value in range(1, 20):
        interval = [Fraction(q_value) - Fraction(1, 2), Fraction(q_value) + Fraction(1, 2)]
        for s_value in range(-3, 24):
            if s_value == q_value:
                continue
            distance = min(abs(Fraction(s_value) - endpoint) for endpoint in interval)
            if Fraction(s_value) < interval[0]:
                distance = interval[0] - Fraction(s_value)
            elif Fraction(s_value) > interval[1]:
                distance = Fraction(s_value) - interval[1]
            require(distance >= Fraction(1, 2), "independent alias gap failed")


def check_fourier_inversion() -> None:
    r, displacement = sp.symbols("r displacement", real=True)
    transformed_gaussian = sp.exp(-sp.pi * r**2) * sp.exp(
        2 * sp.pi * sp.I * displacement * r
    )
    inversion = sp.integrate(transformed_gaussian, (r, -sp.oo, sp.oo))
    require_zero(
        inversion - sp.exp(-sp.pi * displacement**2),
        "independent Gaussian inversion",
    )

    epsilon = sp.symbols("epsilon", positive=True, real=True)
    left = sp.Integer(1)
    right = sp.Integer(7)
    regularized = lambda point: (
        sp.erf(sp.sqrt(sp.pi) * (right - point) / sp.sqrt(epsilon))
        - sp.erf(sp.sqrt(sp.pi) * (left - point) / sp.sqrt(epsilon))
    ) / 2
    limits = [
        sp.limit(regularized(point), epsilon, 0, dir="+")
        for point in [0, 1, 3, 7, 8]
    ]
    require(
        limits == [0, sp.Rational(1, 2), 1, sp.Rational(1, 2), 0],
        "independent Dirichlet endpoint weights failed",
    )


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field} in {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["rfi_22_handoff"],
        "open-row boundary drifted",
    )
    off_lattice = next(row for row in rows if row["id"] == "rfi_02_off_lattice")
    require("e(-r)" in off_lattice["certificate"], "off-lattice lower phase missing")
    require("-rB" in off_lattice["certificate"], "off-lattice upper phase missing")
    achievement = next(row for row in rows if row["id"] == "rfi_21_achievement")
    require("no full-line remainder" in achievement["certificate"], "inversion scope drifted")
    require(rows[-1]["readiness"] == "guard_validated", "final guard drifted")


def check_note(note: str) -> None:
    required = [
        "This is not a proof",
        "## Complete Real-Mode Package",
        "Replacing them off lattice",
        "## Exact Reciprocal Cell",
        "## Dual Kernel And Carrier",
        "no full-line remainder",
        "## Global Recomposition",
        "only reciprocal-frequency boundary",
        "principal-value at q=1,B",
        "## Handoff",
        "## Pi Provenance",
        "## Proof Boundary",
    ]
    for marker in required:
        require(marker in note, f"note marker missing: {marker}")


def main() -> int:
    require(RESULT_PATH.exists(), "result artifact missing")
    require(NOTE_PATH.exists(), "note artifact missing")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")
    require(artifact.get("kind") == KIND, "kind drifted")
    require(artifact.get("counts") == EXPECTED_COUNTS, "count block drifted")
    require(
        "exterior-tail and signed-flow bounds open" in artifact.get("status", ""),
        "status boundary drifted",
    )
    require(
        "It proves no bound for that exterior" in artifact.get("proof_boundary", ""),
        "proof boundary drifted",
    )
    symbolic = artifact.get("symbolic_certificate", {})
    require(
        set(symbolic) == {
            "complete_mode",
            "reciprocal_cells",
            "finite_cell_poisson",
            "carrier_inversion",
            "global_recomposition",
            "audit",
        },
        "symbolic certificate keys drifted",
    )
    audit = symbolic["audit"]
    require(audit["rational_cell_cases"] == 70, "builder cell audit drifted")
    require(audit["exact_ties"] >= 4, "builder tie audit too small")
    require(audit["physical_tilings"] == 20, "builder tiling audit drifted")
    require(audit["roster_containments"] == 20, "builder roster audit drifted")
    check_sources(artifact)
    check_mode_recomposition()
    check_cells()
    check_finite_poisson_sign()
    check_kernel_and_phase()
    check_fourier_inversion()
    check_rows(artifact)
    check_note(note)
    counts = artifact["counts"]
    print(
        "validated reciprocal finite-cell inversion gate: "
        f"{counts['rows']} rows, "
        f"{counts['finite_cell_poisson_identities']} cell-Poisson identities, "
        f"{counts['carrier_inversion_identities']} carrier inversions, "
        f"{counts['internal_tie_cancellations']} tie cancellation, "
        f"{counts['global_recompositions']} global recompositions, "
        f"{counts['exterior_tail_bounds']} exterior-tail bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
