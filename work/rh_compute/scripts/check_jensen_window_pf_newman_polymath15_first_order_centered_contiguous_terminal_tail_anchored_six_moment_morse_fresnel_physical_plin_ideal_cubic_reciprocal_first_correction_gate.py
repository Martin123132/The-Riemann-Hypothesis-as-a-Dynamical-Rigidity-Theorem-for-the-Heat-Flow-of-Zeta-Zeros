#!/usr/bin/env python3
"""Independently check the reciprocal first-correction cancellation gate."""

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
    "ideal_cubic_reciprocal_first_correction_gate"
)
RESULT_PATH = REPO_ROOT / "work" / "rh_compute" / "results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 18,
    "reciprocal_cell_identities": 4,
    "dual_morse_identities": 5,
    "first_correction_identities": 5,
    "first_correction_cancellations": 1,
    "finite_cell_remainder_bounds": 0,
    "incomplete_fresnel_bounds": 0,
    "endpoint_composed_h2_bounds": 0,
    "q_family_bounds": 0,
    "quadratic_remainder_bounds": 0,
    "signed_flow_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    f"irfc_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        [
            "cells",
            "ties",
            "coverage",
            "phase",
            "jacobian",
            "signature",
            "leading",
            "operator",
            "inner",
            "dual",
            "transport",
            "physical",
            "cancel",
            "achievement",
            "cell_boundary",
            "fresnel_boundary",
            "next",
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
        set(source) == {"endpoint_composition", "physical_plin", "reciprocal"},
        "source-audit keys drifted",
    )
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source for {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted for {key}")


def inverse_morse_series() -> tuple[sp.Symbol, sp.Expr, sp.Expr]:
    """Recover the positive inverse chart from its defining phase equation."""
    z = sp.symbols("z", real=True)
    c_2, c_3, c_4 = sp.symbols("c_2 c_3 c_4", real=True)
    rho_trial = 1 + z + c_2 * z**2 + c_3 * z**3 + c_4 * z**4
    residual = sp.series(
        rho_trial - 1 - sp.log(rho_trial) - z**2 / 2,
        z,
        0,
        6,
    ).removeO().expand()
    solutions = sp.solve(
        [residual.coeff(z, power) for power in (3, 4, 5)],
        (c_2, c_3, c_4),
        dict=True,
    )
    require(len(solutions) == 1, "positive inverse Morse series is not unique")
    rho = sp.expand(rho_trial.subs(solutions[0]))
    jacobian = sp.diff(rho, z)
    require_zero(
        sp.series(rho - 1 - sp.log(rho) - z**2 / 2, z, 0, 6).removeO(),
        "inverse Morse phase",
    )
    require_zero(rho.subs(z, 0) - 1, "inverse Morse saddle")
    require_zero(jacobian.subs(z, 0) - 1, "inverse Morse orientation")
    return z, rho, jacobian


def check_symbolics() -> None:
    z, rho, jacobian = inverse_morse_series()
    alpha, q = sp.symbols("alpha_P q", positive=True, real=True)
    a_0, a_1, a_2 = sp.symbols("A_0 A_1 A_2")

    j_rho_first = (sp.diff(jacobian, z) / jacobian).subs(z, 0)
    j_rho_second = (
        sp.diff(sp.diff(jacobian, z) / jacobian, z) / jacobian
    ).subs(z, 0)
    require_zero(j_rho_first - sp.Rational(2, 3), "independent J'(1)")
    require_zero(j_rho_second + sp.Rational(5, 18), "independent J''(1)")

    u_inner = q * rho
    amplitude_inner = (
        a_0 + a_1 * (u_inner - q) + a_2 * (u_inner - q) ** 2 / 2
    )
    inner_amplitude = q * amplitude_inner * jacobian / sp.sqrt(alpha)
    inner_second = sp.simplify(sp.diff(inner_amplitude, z, 2).subs(z, 0) / alpha)

    u_dual = q / rho
    amplitude_dual = a_0 + a_1 * (u_dual - q) + a_2 * (u_dual - q) ** 2 / 2
    dual_amplitude = jacobian * amplitude_dual / rho
    dual_second = sp.simplify(sp.diff(dual_amplitude, z, 2).subs(z, 0) / alpha)

    differential_symbol = q**2 * a_2 + 2 * q * a_1 + a_0 / 6
    require_zero(
        inner_second - q * differential_symbol / alpha ** sp.Rational(3, 2),
        "independent inner second derivative",
    )
    require_zero(
        dual_second - differential_symbol / alpha,
        "independent dual second derivative",
    )
    transported_inner = sp.sqrt(alpha) * inner_second / q
    require_zero(
        transported_inner - dual_second,
        "independent transported second jet",
    )

    lam, t, sigma = sp.symbols("lambda t sigma", real=True)
    p = sp.Function("P")(lam)
    s = t * lam**2 / 4 - sigma * lam
    weighted = sp.exp(s) * p
    euler_form = sp.diff(weighted, lam, 2) + sp.diff(weighted, lam) + weighted / 6
    g = sp.diff(s, lam)
    physical_form = sp.exp(s) * (
        sp.diff(p, lam, 2)
        + (2 * g + 1) * sp.diff(p, lam)
        + (g**2 + g + t / 2 + sp.Rational(1, 6)) * p
    )
    require_zero(euler_form - physical_form, "independent physical operator")

    a = sp.symbols("a", positive=True, real=True)
    kappa = 1 / (2 * sp.pi * sp.I)
    scaled_fresnel = a ** (-sp.Rational(1, 2))
    negative_moment_ratio = sp.simplify(
        sp.diff(scaled_fresnel, a) / (-sp.pi * sp.I * scaled_fresnel)
    ).subs(a, 1)
    positive_moment_ratio = sp.simplify(
        sp.diff(scaled_fresnel, a) / (sp.pi * sp.I * scaled_fresnel)
    ).subs(a, 1)
    require_zero(negative_moment_ratio - kappa, "negative Gaussian moment")
    require_zero(positive_moment_ratio + kappa, "positive Gaussian moment")
    require_zero(
        positive_moment_ratio * dual_second / 2
        + negative_moment_ratio * transported_inner / 2,
        "independent first-correction cancellation",
    )

    rho_symbol = sp.symbols("rho", positive=True, real=True)
    phase_original = (
        alpha * (sp.log(alpha / (alpha * rho_symbol / q)) - 1)
        + alpha * rho_symbol
    )
    phase_morse = alpha * sp.log(q) + alpha * (
        rho_symbol - 1 - sp.log(rho_symbol)
    )
    require_zero(
        sp.expand_log(phase_original - phase_morse, force=True),
        "independent reciprocal phase",
    )
    require_zero(
        q * a_0 / sp.sqrt(alpha) * sp.sqrt(alpha) / q - a_0,
        "independent leading carrier",
    )
    require_zero(
        sp.exp(-sp.I * sp.pi / 4) * sp.exp(sp.I * sp.pi / 4) - 1,
        "independent signature product",
    )


def in_cell(alpha: Fraction, r: int, q: int) -> bool:
    ratio = alpha / r
    return Fraction(q) - Fraction(1, 2) <= ratio < Fraction(q) + Fraction(1, 2)


def check_reciprocal_cells() -> None:
    alphas = [
        Fraction(1, 2),
        Fraction(1),
        Fraction(3, 2),
        Fraction(5, 3),
        Fraction(7, 2),
        Fraction(13),
        Fraction(120),
    ]
    ties_checked = 0
    for alpha in alphas:
        roster = list(range(1, floor_fraction(2 * alpha) + 1))
        q_max = floor_fraction(alpha + Fraction(1, 2)) + 2
        cell_union: list[int] = []
        for q in range(1, q_max + 1):
            direct = [r for r in roster if in_cell(alpha, r, q)]
            lower = floor_fraction(alpha / (Fraction(q) + Fraction(1, 2))) + 1
            upper = floor_fraction(alpha / (Fraction(q) - Fraction(1, 2)))
            formula = [r for r in roster if lower <= r <= upper]
            require(direct == formula, f"cell-bound audit failed at alpha={alpha}, q={q}")
            cell_union.extend(direct)
        require(sorted(cell_union) == roster, f"cell coverage failed at alpha={alpha}")
        require(len(cell_union) == len(set(cell_union)), f"cell overlap at alpha={alpha}")

        for r in roster:
            ratio = alpha / r
            if ratio.denominator != 2 or ratio < Fraction(3, 2):
                continue
            q = floor_fraction(ratio - Fraction(1, 2))
            require(not in_cell(alpha, r, q), "upper tie remained in lower cell")
            require(in_cell(alpha, r, q + 1), "upper tie did not transfer upward")
            ties_checked += 1
    require(ties_checked >= 4, "insufficient exact half-integer tie coverage")


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field} in {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["irfc_15_cell_boundary", "irfc_16_fresnel_boundary", "irfc_17_next"],
        "open-row boundary drifted",
    )
    cancellation = next(row for row in rows if row["id"] == "irfc_13_cancel")
    require(
        "full-line stationary-symbol" in cancellation["proof_boundary"],
        "cancellation scope drifted",
    )
    final = rows[-1]
    require(final["readiness"] == "guard_validated", "final guard drifted")
    require("No grouped h2" in final["proof_boundary"], "final boundary overpromoted")


def check_note(note: str) -> None:
    required = [
        "This is not a proof of a finite-cell h^2 estimate",
        "## Reciprocal Cells",
        "## Dual Morse Chart",
        "## First Correction",
        "their sum is exactly zero",
        "No O(alpha_P^(-2)) finite-cell remainder is claimed",
        "## Remaining Boundary",
        "## Handoff",
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
        "finite-cell and incomplete-endpoint defects open" in artifact.get("status", ""),
        "status boundary drifted",
    )
    require(
        "It proves no finite reciprocal-cell remainder" in artifact.get("proof_boundary", ""),
        "proof boundary drifted",
    )
    symbolic = artifact.get("symbolic_certificate", {})
    require(
        set(symbolic) == {
            "cell_partition",
            "dual_morse",
            "first_correction",
            "series_audit",
        },
        "symbolic certificate keys drifted",
    )
    check_sources(artifact)
    check_symbolics()
    check_reciprocal_cells()
    check_rows(artifact)
    check_note(note)
    counts = artifact["counts"]
    print(
        "validated reciprocal first-correction gate: "
        f"{counts['rows']} rows, "
        f"{counts['reciprocal_cell_identities']} cell identities, "
        f"{counts['dual_morse_identities']} dual-Morse identities, "
        f"{counts['first_correction_identities']} first-correction identities, "
        f"{counts['first_correction_cancellations']} cancellation, "
        f"{counts['finite_cell_remainder_bounds']} finite-cell bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
