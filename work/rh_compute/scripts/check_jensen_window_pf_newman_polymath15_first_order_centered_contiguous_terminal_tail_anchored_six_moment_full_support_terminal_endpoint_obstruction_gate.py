#!/usr/bin/env python3
"""Independently check the full-support terminal endpoint obstruction gate."""

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
    "terminal_tail_anchored_six_moment_full_support_terminal_endpoint_"
    "obstruction_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 24,
    "terminal_rate_identities": 3,
    "terminal_polynomial_identities": 4,
    "terminal_channel_bounds": 5,
    "row_relative_collapses": 1,
    "lower_endpoint_s1_bounds": 1,
    "conditional_lift_identities": 1,
    "quadratic_conditional_suppression_templates": 1,
    "unique_order_one_obstructions": 1,
    "signed_terminal_kernel_bounds": 0,
    "complete_outer_complement_bounds": 0,
    "quadratic_remainder_bounds": 0,
    "signed_full_current_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    "teo_01_rates",
    "teo_02_rows",
    "teo_03_falpha",
    "teo_04_ideal",
    "teo_05_cbound",
    "teo_06_dbound",
    "teo_07_rowbounds",
    "teo_08_collapse",
    "teo_09_lower_s1",
    "teo_10_lower_geometry",
    "teo_11_lift",
    "teo_12_terminal_lift",
    "teo_13_quadratic",
    "teo_14_unique",
    "teo_15_handoff",
    "teo_16_ties",
    "teo_17_edge",
    "teo_18_phase",
    "teo_19_move",
    "teo_20_correction",
    "teo_21_route",
    "teo_22_pi",
    "teo_23_scope",
    "teo_24_boundary",
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


def check_sources(artifact: dict) -> None:
    expected = {
        "physical_plin",
        "correction_box",
        "rate_box",
        "full_support",
        "outer_endpoint",
    }
    source = artifact.get("source_audit", {})
    require(set(source) == expected, "source-audit keys drifted")
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def check_terminal_algebra() -> None:
    i = sp.I
    u, b, bx, ux, chi, C, D = sp.symbols("u b bx ux chi C D")
    c, cx = sp.symbols("c cx")
    x = -u
    s1 = c + i * b
    s2 = cx + i * bx
    R = -s1 * x - c * u
    Rx = -s2 * x - cx * u + i * b * ux
    require_zero(R - i * b * u, "independent terminal R")
    require_zero(Rx - i * (bx * u + b * ux), "independent terminal Rx")
    delta = chi - i * b * u
    require_zero(R + delta - chi, "independent terminal R+delta")

    A = R * C
    Q = (R + delta) * C + D
    Npoly = R * Q + Rx * C
    np, vp, qp, ap = sp.symbols("np vp qp ap", real=True)
    F = np * C + vp * Npoly - qp * A - ap * Q
    explicit = (
        np * C
        + vp * (i * b * u * (chi * C + D) + i * (bx * u + b * ux) * C)
        - qp * i * b * u * C
        - ap * (chi * C + D)
    )
    require_zero(F - explicit, "independent F_alpha")
    ideal = sp.simplify(
        F.subs({C: 1, D: 0, b: -sp.Rational(1, 2), bx: 0, chi: -i * u / 2})
    )
    expected = np + vp * (-u**2 / 4 - i * ux / 2) + i * u * (qp + ap) / 2
    require_zero(ideal - expected, "independent ideal obstruction")


def check_rational_bounds() -> None:
    for denominator in [72_000_000_00, 72_000_000_000, 10**12, 10**18]:
        h = Fraction(1, denominator)
        c_error = Fraction(4379) + h / 2 + h**2 / 4
        d_size = Fraction(16893) + 5 * h / 8 + h**2 / 4
        require(c_error < 4380, "independent C error")
        require(d_size < 16894, "independent D size")
        require(1 + 4380 * h**2 < Fraction(3, 2), "independent C modulus")

        r = Fraction(3001, 3000)
        chi = r + 4 * h
        require(chi < Fraction(101, 100), "independent chi")
        require(r * Fraction(3, 2) < Fraction(8, 5), "independent A")
        q = Fraction(101, 100) * Fraction(3, 2) + 16894 * h**3
        require(q < Fraction(8, 5), "independent Q")
        rx = h**3 / 8 + Fraction(3001, 144000)
        require(rx < Fraction(1, 20), "independent Rx")
        n = r * Fraction(8, 5) + Fraction(3, 40)
        require(n < Fraction(17, 10), "independent N polynomial")
        require(4380 * h < Fraction(8, 5), "independent row collapse C")
        require(Fraction(17, 10) * h < Fraction(8, 5), "independent row collapse N")


def check_lower_endpoint_s1() -> int:
    mp.mp.dps = 70
    checks = 0
    for n_value in range(2, 41):
        for theta in [Fraction(0), Fraction(1, 7), Fraction(1, 2), Fraction(6, 7), Fraction(1)]:
            a_value = Fraction(n_value) + theta
            for offset in [Fraction(1, 19), Fraction(4, 9), Fraction(18, 19)]:
                alpha = a_value**2 - offset
                a_n = alpha / (Fraction(n_value) + Fraction(1, 2))
                lower = floor_fraction(a_n) + 1
                upper = floor_fraction(2 * alpha)
                a_x = alpha - lower + 1
                b_x = upper + 1 - alpha
                require(a_x > Fraction(3, 5) * alpha, "independent a_1 lower")
                length = b_x - a_x
                require(length > 0, "independent endpoint ordering")
                require(length < 1 + a_n, "independent endpoint length")
                a_mp = mp.mpf(a_x.numerator) / a_x.denominator
                b_mp = mp.mpf(b_x.numerator) / b_x.denominator
                s1 = mp.digamma(b_mp) - mp.digamma(a_mp)
                require(0 < s1 < mp.mpf(5) / n_value, "independent S1 bound")
                checks += 1
    require(checks > 500, "insufficient S1 checks")
    return checks


def check_lift_pairing() -> None:
    i = sp.I
    ell = sp.symbols("ell", real=True)
    phase, amplitude, kernel, kappa = sp.symbols("phase amplitude kernel kappa")
    base = phase * kappa * amplitude * kernel
    lifted = phase * kappa * (i * ell * amplitude) * kernel
    require_zero(lifted - i * ell * base, "independent lifted trace")

    z = sp.symbols("z0:8", real=True)
    p = z[0] + i * z[1]
    q = z[2] + i * z[3]
    p2 = z[4] + i * z[5]
    q2 = z[6] + i * z[7]
    real_product = sp.re(p) * sp.re(i * ell * q)
    split = sp.re(p * sp.conjugate(i * ell * q) + p * i * ell * q) / 2
    require_zero(sp.expand_complex(real_product - split), "independent conditional H/T")
    require_zero(sp.re(p2) * sp.re(q2) - sp.re(p2 * sp.conjugate(q2) + p2 * q2) / 2, "independent real pairing")


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["teo_15_handoff"],
        "open-row boundary drifted",
    )


def check_note(note: str) -> None:
    markers = [
        "This is not a proof of RH",
        "## Terminal Specialization",
        "## Row Bounds",
        "## Opposite Endpoint",
        "## Conditional Quadratic Traces",
        "mathcal N_(p,xi)",
        "(8/5)h",
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
    require("signed kernel open" in artifact.get("status", ""), "status drifted")
    require(
        "no signed terminal kernel" in artifact.get("proof_boundary", ""),
        "proof boundary drifted",
    )
    cert = artifact.get("symbolic_certificate", {})
    require(
        set(cert)
        == {
            "terminal_polynomials",
            "row_bounds",
            "conditional_traces",
            "obstruction",
            "audit",
        },
        "certificate keys drifted",
    )
    require(
        cert["audit"]
        == {"lower_endpoint_s1_cases": 264, "terminal_row_bound_cases": 3},
        "builder audit drifted",
    )
    check_sources(artifact)
    check_terminal_algebra()
    check_rational_bounds()
    s1_checks = check_lower_endpoint_s1()
    check_lift_pairing()
    check_rows(artifact)
    check_note(note)
    print(
        "validated full-support terminal endpoint obstruction gate: "
        f"{EXPECTED_COUNTS['rows']} rows, {s1_checks} S1 checks, "
        f"{EXPECTED_COUNTS['row_relative_collapses']} row collapse, "
        f"{EXPECTED_COUNTS['unique_order_one_obstructions']} order-one obstruction, "
        f"{EXPECTED_COUNTS['signed_terminal_kernel_bounds']} signed kernel bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
