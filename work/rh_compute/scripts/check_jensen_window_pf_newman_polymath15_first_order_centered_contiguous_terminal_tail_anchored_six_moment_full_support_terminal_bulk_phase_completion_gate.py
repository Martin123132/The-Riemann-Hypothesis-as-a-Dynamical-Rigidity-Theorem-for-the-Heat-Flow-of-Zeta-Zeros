#!/usr/bin/env python3
"""Independently check the terminal-bulk phase-completion gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_terminal_bulk_phase_"
    "completion_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 28,
    "terminal_bulk_polarizations": 1,
    "phase_families": 2,
    "anti_conjugate_kernel_relations": 1,
    "ideal_coefficient_pairs": 1,
    "algebraic_noncancellation_witnesses": 1,
    "observation_coefficient_completions": 1,
    "completed_component_pairings": 12,
    "fixed_terminal_linear_splits": 1,
    "pure_terminal_hermitian_nulls": 1,
    "pure_terminal_transpose_collapses": 1,
    "carrier_determinant_collapses": 1,
    "pure_terminal_rational_bounds": 1,
    "signed_completed_terminal_nonterminal_bounds": 0,
    "complete_outer_complement_bounds": 0,
    "quadratic_remainder_bounds": 0,
    "signed_full_current_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    "tbp_01_atom",
    "tbp_02_product",
    "tbp_03_hermitian",
    "tbp_04_transpose",
    "tbp_05_anti",
    "tbp_06_phases",
    "tbp_07_ideal",
    "tbp_08_orientation",
    "tbp_09_witness",
    "tbp_10_current",
    "tbp_11_split",
    "tbp_12_completion",
    "tbp_13_rows",
    "tbp_14_components",
    "tbp_15_leading",
    "tbp_16_affine",
    "tbp_17_amplitudes",
    "tbp_18_lift",
    "tbp_19_hnull",
    "tbp_20_transpose",
    "tbp_21_determinant",
    "tbp_22_collapse",
    "tbp_23_bound",
    "tbp_24_scale",
    "tbp_25_route",
    "tbp_26_pi",
    "tbp_27_join",
    "tbp_28_boundary",
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


def scalar(expression: sp.Matrix) -> sp.Expr:
    return sp.expand(expression[0])


def check_sources(artifact: dict) -> None:
    expected = {
        "terminal_conditional",
        "terminal_obstruction",
        "outer_endpoint",
        "observation_image",
        "physical_plin",
        "full_support",
    }
    source = artifact.get("source_audit", {})
    require(set(source) == expected, "source-audit keys drifted")
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def check_phase_split() -> None:
    imaginary = sp.I
    s, u = sp.symbols("s u", real=True)
    ar, ai, br, bi, cr, ci = sp.symbols("ar ai br bi cr ci", real=True)
    terminal = ar + imaginary * ai
    carrier = br + imaginary * bi
    row = cr + imaginary * ci
    atom = -imaginary * u * row * carrier
    lhs = s * sp.im(terminal) * sp.re(atom) / (2 * sp.pi)
    hermitian_coefficient = s * u * sp.conjugate(row) / (4 * sp.pi)
    transpose_coefficient = -s * u * row / (4 * sp.pi)
    rhs = sp.re(
        hermitian_coefficient * terminal * sp.conjugate(carrier)
        + transpose_coefficient * terminal * carrier
    )
    require_zero(sp.expand_complex(lhs - rhs), "independent H/T split")
    require_zero(
        transpose_coefficient + sp.conjugate(hermitian_coefficient),
        "independent anti-conjugacy",
    )

    ux = sp.symbols("ux", real=True)
    ideal_row = -u**2 / 4 - imaginary * ux / 2
    ideal_h = sp.simplify(hermitian_coefficient.subs({cr: -u**2 / 4, ci: -ux / 2}))
    ideal_t = sp.simplify(transpose_coefficient.subs({cr: -u**2 / 4, ci: -ux / 2}))
    require_zero(
        ideal_h
        - (-s * u**3 / (16 * sp.pi) + imaginary * s * u * ux / (8 * sp.pi)),
        "independent ideal H",
    )
    require_zero(
        ideal_t
        - (s * u**3 / (16 * sp.pi) + imaginary * s * u * ux / (8 * sp.pi)),
        "independent ideal T",
    )
    witness = rhs.subs(
        {
            ar: 0,
            ai: 1,
            br: 0,
            bi: 1,
            cr: sp.re(ideal_row.subs({u: 1, ux: 0})),
            ci: sp.im(ideal_row.subs({u: 1, ux: 0})),
            u: 1,
        }
    )
    require_zero(witness + s / (8 * sp.pi), "independent nonzero witness")


def check_coefficient_completion() -> None:
    j_matrix = sp.Matrix(
        [
            [0, sp.Rational(1, 2), 0, 0],
            [sp.Rational(1, 2), 0, 0, 0],
            [0, 0, 0, -sp.Rational(1, 2)],
            [0, 0, -sp.Rational(1, 2), 0],
        ]
    )
    base = sp.Matrix(sp.symbols("b0:4", real=True))
    dot_base = sp.Matrix(sp.symbols("db0:4", real=True))
    terminal = sp.Matrix(sp.symbols("a0:4", real=True))
    dot_terminal = sp.Matrix(sp.symbols("da0:4", real=True))
    rest = sp.Matrix(sp.symbols("r0:4", real=True))
    dot_rest = sp.Matrix(sp.symbols("dr0:4", real=True))

    def current(error: sp.Matrix, dot_error: sp.Matrix) -> sp.Expr:
        return scalar(
            2 * base.T * j_matrix * dot_error
            + 2 * dot_base.T * j_matrix * error
            + 2 * error.T * j_matrix * dot_error
        )

    complete = base + rest
    dot_complete = dot_base + dot_rest
    terminal_channel = scalar(
        2 * complete.T * j_matrix * dot_terminal
        + 2 * dot_complete.T * j_matrix * terminal
        + 2 * terminal.T * j_matrix * dot_terminal
    )
    require_zero(
        current(terminal + rest, dot_terminal + dot_rest)
        - current(rest, dot_rest)
        - terminal_channel,
        "independent coefficient completion",
    )
    expected_components = (
        complete[0] * dot_terminal[1]
        + complete[1] * dot_terminal[0]
        - complete[2] * dot_terminal[3]
        - complete[3] * dot_terminal[2]
        + dot_complete[0] * terminal[1]
        + dot_complete[1] * terminal[0]
        - dot_complete[2] * terminal[3]
        - dot_complete[3] * terminal[2]
        + terminal[0] * dot_terminal[1]
        + terminal[1] * dot_terminal[0]
        - terminal[2] * dot_terminal[3]
        - terminal[3] * dot_terminal[2]
    )
    require_zero(terminal_channel - expected_components, "independent components")
    mixed = scalar(
        2 * complete.T * j_matrix * dot_terminal
        + 2 * dot_complete.T * j_matrix * terminal
    )
    require_zero(
        sp.diff(mixed, terminal[0]) - dot_complete[1],
        "independent completed coefficient",
    )

    fixed = sp.Matrix(sp.symbols("t0:4", real=True))
    require_zero(
        scalar(fixed.T * (dot_terminal + dot_rest))
        - scalar(fixed.T * dot_terminal)
        - scalar(fixed.T * dot_rest),
        "independent affine split",
    )


def check_pure_terminal() -> None:
    imaginary = sp.I
    u = sp.symbols("u", real=True)
    parts = sp.symbols("vr vi nr ni ar ai qr qi", real=True)
    v = parts[0] + imaginary * parts[1]
    n = parts[2] + imaginary * parts[3]
    a = parts[4] + imaginary * parts[5]
    q = parts[6] + imaginary * parts[7]
    dv, dn, da, dq = [-imaginary * u * z for z in (v, n, a, q)]
    hermitian = (
        v * sp.conjugate(dn)
        + n * sp.conjugate(dv)
        - a * sp.conjugate(dq)
        - q * sp.conjugate(da)
    )
    transpose = v * dn + n * dv - a * dq - q * da
    require_zero(sp.expand_complex(sp.re(hermitian)), "independent Hermitian null")
    require_zero(
        sp.expand(transpose + 2 * imaginary * u * (v * n - a * q)),
        "independent transpose reduction",
    )

    tau, rate, rate_x, carrier, qrow = sp.symbols("tau rate rate_x carrier qrow")
    arow = rate * carrier
    nrow = rate * qrow + rate_x * carrier
    determinant = carrier * nrow - arow * qrow
    require_zero(determinant - rate_x * carrier**2, "independent determinant")
    collapsed = -2 * imaginary * u * tau**2 * determinant
    require_zero(
        collapsed + 2 * imaginary * u * tau**2 * rate_x * carrier**2,
        "independent terminal transpose collapse",
    )


def check_rational_bounds() -> None:
    coefficient = Fraction(2) * Fraction(1, 36) * Fraction(1, 20) * Fraction(9, 4)
    require(coefficient == Fraction(1, 160), "independent 1/160 coefficient")
    for denominator in [72_000_000_000, 10**12, 10**15, 10**18, 10**24]:
        h = Fraction(1, denominator)
        upper = 2 * h * Fraction(1, 36) * h**2 / 20 * Fraction(9, 4)
        require(upper == h**3 / 160, "independent pure-terminal envelope")
        require(upper < h**2, "independent h-cubed separation")


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["tbp_25_route"],
        "open-row boundary drifted",
    )


def check_note(note: str) -> None:
    markers = [
        "not a proof of RH",
        "## Terminal-Bulk Split",
        "K_H(N,q)",
        "K_T(N,q)",
        "## Coefficient Completion",
        "omega_*=omega_p+rho",
        "## Pure Conditional Current",
        "Re(H_aa)=0",
        "C_Nmathcal N_N-A_NQ_N",
        "|w_N|^2/160",
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
    require("signed completed estimate open" in artifact.get("status", ""), "status drifted")
    require(
        "no signed completed terminal/nonterminal estimate"
        in artifact.get("proof_boundary", ""),
        "proof boundary drifted",
    )
    cert = artifact.get("symbolic_certificate", {})
    require(
        set(cert)
        == {
            "terminal_bulk_split",
            "ideal_phase_chart",
            "coefficient_completion",
            "pure_terminal_current",
            "bounds",
            "audit",
            "handoff",
        },
        "certificate keys drifted",
    )
    require(
        cert["audit"]
        == {
            "phase_families": 2,
            "coefficient_completions": 1,
            "pure_hermitian_nulls": 1,
            "pure_transpose_collapses": 1,
        },
        "builder audit drifted",
    )
    check_sources(artifact)
    check_phase_split()
    check_coefficient_completion()
    check_pure_terminal()
    check_rational_bounds()
    check_rows(artifact)
    check_note(note)
    print(
        "validated full-support terminal-bulk phase completion gate: "
        f"{EXPECTED_COUNTS['rows']} rows, "
        f"{EXPECTED_COUNTS['phase_families']} phase families, "
        f"{EXPECTED_COUNTS['observation_coefficient_completions']} coefficient completion, "
        f"{EXPECTED_COUNTS['pure_terminal_hermitian_nulls']} Hermitian null, "
        f"{EXPECTED_COUNTS['pure_terminal_transpose_collapses']} transpose collapse, "
        f"{EXPECTED_COUNTS['signed_completed_terminal_nonterminal_bounds']} signed completed bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
