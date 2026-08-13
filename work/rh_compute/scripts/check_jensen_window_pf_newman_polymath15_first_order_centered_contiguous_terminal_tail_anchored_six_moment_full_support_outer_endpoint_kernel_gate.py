#!/usr/bin/env python3
"""Independently check the full-support outer-endpoint kernel gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_outer_endpoint_kernel_gate"
)
RESULT_PATH = REPO_ROOT / "work" / "rh_compute" / "results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 27,
    "outer_band_rosters": 1,
    "pole_free_kernel_families": 3,
    "integer_removable_formulas": 1,
    "twofold_outer_recompositions": 1,
    "stable_endpoint_functionals": 2,
    "absolute_kernel_bounds": 3,
    "full_support_linear_functionals": 1,
    "terminal_plin_evaluations": 1,
    "edge_endpoint_factorizations": 1,
    "edge_terminal_suppression_templates": 1,
    "hermitian_transpose_decompositions": 1,
    "moving_tail_placements": 1,
    "complete_outer_complement_bounds": 0,
    "quadratic_remainder_bounds": 0,
    "signed_full_current_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    "oek_01_band",
    "oek_02_distance",
    "oek_03_tie",
    "oek_04_s1",
    "oek_05_s2",
    "oek_06_s3",
    "oek_07_integer",
    "oek_08_kernel_bounds",
    "oek_09_mode_ibp",
    "oek_10_operator",
    "oek_11_endpoint",
    "oek_12_outer_sum",
    "oek_13_convergence",
    "oek_14_linear",
    "oek_15_terminal_eval",
    "oek_16_edge_split",
    "oek_17_edge_endpoint",
    "oek_18_suppression",
    "oek_19_noncancel",
    "oek_20_quadratic",
    "oek_21_phase",
    "oek_22_move",
    "oek_23_route",
    "oek_24_handoff",
    "oek_25_pi",
    "oek_26_scope",
    "oek_27_boundary",
]


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(expression)
    if value != 0:
        raise RuntimeError(f"{label}: {value}")


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def harmonic(n: int) -> Fraction:
    return sum((Fraction(1, k) for k in range(1, n + 1)), Fraction(0))


def check_sources(artifact: dict) -> None:
    expected = {
        "morse_endpoint",
        "endpoint_composition",
        "observation_image",
        "physical_plin",
        "finite_cell",
        "full_support",
    }
    sources = artifact.get("source_audit", {})
    require(set(sources) == expected, "source-audit keys drifted")
    for key, row in sources.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"source missing: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def check_ibp_algebra() -> None:
    A, A1, A2, q, q1, q2 = sp.symbols("A A1 A2 q q1 q2", nonzero=True)
    d = A1 / q - A * q1 / q**2
    d_over_q = d / q
    differentiated = (
        A2 / q**2
        - 3 * A1 * q1 / q**3
        - A * q2 / q**3
        + 3 * A * q1**2 / q**4
    )
    manual = (
        A2 / q**2
        - 2 * A1 * q1 / q**3
        - A1 * q1 / q**3
        - A * q2 / q**3
        + 3 * A * q1**2 / q**4
    )
    require_zero(d_over_q - A1 / q**2 + A * q1 / q**3, "independent D/q")
    require_zero(differentiated - manual, "independent C operator")


def check_edge_factorization() -> None:
    i = sp.I
    u, F, Fp, g, N, alpha = sp.symbols("u F Fp g N alpha", nonzero=True)
    kappa, S1, S2, S3 = sp.symbols("kappa S1 S2 S3")
    P = -i * u * F
    P_lambda = i * F - i * u * Fp
    A_lambda = P_lambda + g * P
    require_zero(A_lambda - i * (F - u * (Fp + g * F)), "edge derivative")
    direct = (
        kappa * P * S1
        - kappa**2 * A_lambda * S2 / N
        - kappa**2 * alpha * P * S3 / N**2
    )
    factored = i * (
        -kappa * u * F * S1
        - kappa**2 * (F - u * (Fp + g * F)) * S2 / N
        + kappa**2 * alpha * u * F * S3 / N**2
    )
    require_zero(direct - factored, "edge endpoint factorization")
    Fa, Fb = sp.symbols("Fa Fb")
    require_zero(Fa + i * (-u) * Fb - (Fa - i * u * Fb), "terminal P_lin")


def check_quadratic_split() -> None:
    values = sp.symbols("a0:16", real=True)
    z_v = values[0] + sp.I * values[1]
    z_n = values[2] + sp.I * values[3]
    z_a = values[4] + sp.I * values[5]
    z_q = values[6] + sp.I * values[7]
    d_v = values[8] + sp.I * values[9]
    d_n = values[10] + sp.I * values[11]
    d_a = values[12] + sp.I * values[13]
    d_q = values[14] + sp.I * values[15]
    real_form = (
        sp.re(z_v) * sp.re(d_n)
        + sp.re(z_n) * sp.re(d_v)
        - sp.re(z_a) * sp.re(d_q)
        - sp.re(z_q) * sp.re(d_a)
    )
    h_form = (
        z_v * sp.conjugate(d_n)
        + z_n * sp.conjugate(d_v)
        - z_a * sp.conjugate(d_q)
        - z_q * sp.conjugate(d_a)
    )
    t_form = z_v * d_n + z_n * d_v - z_a * d_q - z_q * d_a
    require_zero(
        sp.expand_complex(real_form - sp.re(h_form + t_form) / 2),
        "quadratic H/T split",
    )


def check_band_kernels() -> tuple[int, int, int]:
    mp.mp.dps = 70
    point_count = 0
    numeric_count = 0
    bound_count = 0
    for n_value in range(2, 17):
        for theta in [Fraction(0), Fraction(2, 7), Fraction(5, 6)]:
            a_value = Fraction(n_value) + theta
            for offset in [Fraction(1, 13), Fraction(5, 11), Fraction(12, 13)]:
                alpha = a_value**2 - offset
                a_n = alpha / (Fraction(n_value) + Fraction(1, 2))
                lower = floor_fraction(a_n) + 1
                upper = floor_fraction(2 * alpha)
                require(lower <= upper, "independent band empty")
                for u_value in range(1, n_value + 1):
                    x = alpha / u_value
                    a_x = x - lower + 1
                    b_x = upper + 1 - x
                    require(a_x > Fraction(1, 4), "independent lower distance")
                    require(b_x > alpha, "independent upper distance")
                    if x.denominator == 1:
                        require(lower <= x <= upper, "independent integer containment")
                    a_mp = mp.mpf(a_x.numerator) / a_x.denominator
                    b_mp = mp.mpf(b_x.numerator) / b_x.denominator
                    s2 = mp.zeta(2, a_mp) + mp.zeta(2, b_mp)
                    s3 = mp.zeta(3, a_mp) - mp.zeta(3, b_mp)
                    z4 = mp.zeta(4, a_mp) + mp.zeta(4, b_mp)
                    require(s2 < 21 and abs(s3) < 73 and z4 < 278, "kernel envelope")
                    bound_count += 3

                    if x.denominator != 1 and (u_value in {1, n_value}):
                        x_mp = mp.mpf(x.numerator) / x.denominator
                        sum1 = mp.fsum(1 / (x_mp - r) for r in range(lower, upper + 1))
                        stable = mp.digamma(b_mp) - mp.digamma(a_mp)
                        cot = mp.pi / mp.tan(mp.pi * x_mp) - sum1
                        require(abs(stable - cot) < mp.mpf("1e-55"), "independent S1")
                        numeric_count += 1
                    point_count += 1
    require(point_count > 1100, "insufficient independent kernel points")
    require(numeric_count > 200, "insufficient independent S1 checks")
    return point_count, numeric_count, bound_count


def check_integer_removals() -> int:
    checks = 0
    for lower in range(1, 11):
        for upper in range(lower + 2, lower + 14):
            for x in range(lower, upper + 1):
                stable = harmonic(upper - x) - harmonic(x - lower)
                radius = upper + x + 20
                finite = sum(
                    (
                        Fraction(1, x - r)
                        for r in range(-radius, radius + 1)
                        if r < lower or r > upper
                    ),
                    Fraction(0),
                )
                pv = harmonic(radius + x) - harmonic(radius - x)
                require(finite == stable + pv, "independent removable identity")
                checks += 1
    require(checks > 900, "insufficient removable checks")
    return checks


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty row field: {row.get('id')} {field}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["oek_24_handoff"],
        "open-row boundary drifted",
    )


def check_note(note: str) -> None:
    markers = [
        "This is not a proof of RH",
        "## Reciprocal Band",
        "## Pole-Free Kernels",
        "## Outer Recomposition",
        "## Physical Pullthrough",
        "## Quadratic Channels",
        "P_lin^[N](log N)",
        "u_N<2h",
        "## Pi Provenance",
        "## Proof Boundary",
    ]
    for marker in markers:
        require(marker in note, f"note marker missing: {marker}")


def main() -> int:
    require(RESULT_PATH.exists(), "result artifact missing")
    require(NOTE_PATH.exists(), "note artifact missing")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")
    require(artifact.get("kind") == KIND, "kind drifted")
    require(artifact.get("counts") == EXPECTED_COUNTS, "counts drifted")
    require("complete signed estimate open" in artifact.get("status", ""), "status drifted")
    require(
        "no complete outer-complement bound" in artifact.get("proof_boundary", ""),
        "proof boundary drifted",
    )
    cert = artifact.get("symbolic_certificate", {})
    require(
        set(cert)
        == {
            "outer_band",
            "stable_kernels",
            "outer_identity",
            "physical_pullthrough",
            "quadratic_channels",
            "audit",
        },
        "certificate keys drifted",
    )
    require(
        cert["audit"]
        == {
            "absolute_kernel_checks": 783,
            "band_cases": 36,
            "kernel_points": 261,
            "noninteger_kernel_checks": 780,
            "rational_symmetric_checks": 294,
            "removable_integer_checks": 294,
        },
        "builder audit drifted",
    )
    check_sources(artifact)
    check_ibp_algebra()
    check_edge_factorization()
    check_quadratic_split()
    points, numeric, bounds = check_band_kernels()
    removals = check_integer_removals()
    check_rows(artifact)
    check_note(note)
    print(
        "validated full-support outer endpoint kernel gate: "
        f"{EXPECTED_COUNTS['rows']} rows, {points} kernel points, "
        f"{numeric} S1 checks, {bounds} kernel bounds, "
        f"{removals} removable checks, "
        f"{EXPECTED_COUNTS['complete_outer_complement_bounds']} complete outer bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
