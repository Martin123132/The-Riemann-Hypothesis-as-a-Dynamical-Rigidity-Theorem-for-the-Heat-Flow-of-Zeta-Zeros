#!/usr/bin/env python3
"""Independently check the terminal full-support reciprocal homotopy gate."""

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
    "ideal_cubic_reciprocal_terminal_full_support_homotopy_gate"
)
RESULT_PATH = REPO_ROOT / "work" / "rh_compute" / "results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 27,
    "exact_current_polarizations": 3,
    "all_carrier_homotopies": 1,
    "physical_anchor_identities": 1,
    "moving_tail_defect_identities": 1,
    "six_moment_full_support_closures": 1,
    "edge_only_matrix_reductions": 1,
    "full_support_linear_functionals": 1,
    "terminal_extension_identities": 3,
    "terminal_reciprocal_strip_tilings": 1,
    "full_support_cell_inversions": 1,
    "outer_complement_gap_identities": 2,
    "terminal_phase_displacement_bounds": 1,
    "outer_complement_bounds": 0,
    "quadratic_remainder_bounds": 0,
    "signed_full_current_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    "rth_01_domain",
    "rth_02_current",
    "rth_03_bulk",
    "rth_04_full",
    "rth_05_anchor",
    "rth_06_derivative",
    "rth_07_move_defect",
    "rth_08_moments",
    "rth_09_matrix",
    "rth_10_edge",
    "rth_11_plin",
    "rth_12_quadratic",
    "rth_13_extension",
    "rth_14_starred",
    "rth_15_mode_extension",
    "rth_16_endpoint_guard",
    "rth_17_cells",
    "rth_18_terminal_strip",
    "rth_19_returned_tail",
    "rth_20_involution",
    "rth_21_low_gap",
    "rth_22_high_gap",
    "rth_23_phase_motion",
    "rth_24_route",
    "rth_25_handoff",
    "rth_26_pi",
    "rth_27_boundary",
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
    expected = {"abel_join", "terminal_prefix", "flow_matrix", "physical_plin", "finite_cell"}
    source = artifact.get("source_audit", {})
    require(set(source) == expected, "source-audit keys drifted")
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def check_current_algebra() -> None:
    J = sp.Rational(1, 2) * sp.Matrix(
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]]
    )
    b = sp.Matrix(sp.symbols("b0:4", real=True))
    c = sp.Matrix(sp.symbols("c0:4", real=True))
    e = sp.Matrix(sp.symbols("e0:4", real=True))
    db = sp.Matrix(sp.symbols("db0:4", real=True))
    dc = sp.Matrix(sp.symbols("dc0:4", real=True))
    c0 = sp.Matrix(sp.symbols("cp0:4", real=True))

    q = lambda w: (w.T * J * w)[0]
    require_zero(q(b) - b[0] * b[1] + b[2] * b[3], "independent current")
    old_physical = q(b) + 2 * ((e + c).T * J * b)[0]
    full_physical = q(e + b + c)
    require_zero(full_physical - old_physical - q(e + c), "independent anchor")

    full_derivative = 2 * ((e + b + c).T * J * (db + dc))[0]
    expanded = (
        (db[0] + dc[0]) * (b[1] + c[1] + e[1])
        + (b[0] + c[0] + e[0]) * (db[1] + dc[1])
        - (db[2] + dc[2]) * (b[3] + c[3] + e[3])
        - (b[2] + c[2] + e[2]) * (db[3] + dc[3])
    )
    require_zero(full_derivative - expanded, "independent full derivative")

    old_derivative = 2 * (b.T * J * db)[0] + 2 * ((e + c0).T * J * db)[0]
    defect = 2 * ((c - c0).T * J * db)[0] + 2 * ((e + b + c).T * J * dc)[0]
    require_zero(full_derivative - old_derivative - defect, "independent moving defect")

    p = sp.Matrix(sp.symbols("p0:4", real=True))
    residual = sp.Matrix(sp.symbols("r0:4", real=True))
    split = q(e + p + residual) - q(e + p)
    gradient = 2 * ((e + p).T * J * residual)[0]
    require_zero(split - gradient - q(residual), "independent split expansion")


def check_matrix_reduction() -> None:
    J = sp.Rational(1, 2) * sp.Matrix(
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]]
    )
    for seed in range(1, 12):
        rows = 5
        u0 = sp.Matrix(rows, 4, lambda i, j: ((seed + 2 * i - 3 * j) % 11) - 5)
        u1 = sp.Matrix(rows, 4, lambda i, j: ((2 * seed - i + 4 * j) % 13) - 6)
        y = sp.Matrix([seed + i - 2 for i in range(rows)])
        edge = sp.Matrix([seed, 1 - seed, seed + 2, 3 - seed])
        matrix = u0 * J * u1.T + u1 * J * u0.T
        linear = u1 * (2 * J * edge)
        omega = u0.T * y
        dot_omega = u1.T * y
        require_zero(
            (y.T * matrix * y)[0] + (linear.T * y)[0]
            - 2 * ((edge + omega).T * J * dot_omega)[0],
            f"matrix reduction seed {seed}",
        )


def check_terminal_extension() -> None:
    for b_value in range(2, 12):
        for n_value in range(b_value + 1, b_value + 8):
            samples = [Fraction((q + 2) * (q + 3), q + 1) for q in range(1, n_value + 1)]
            h_b = sum(samples[:b_value])
            h_n = sum(samples)
            e_b = Fraction(1, 2) * (samples[0] + samples[b_value - 1])
            e_n = Fraction(1, 2) * (samples[0] + samples[n_value - 1])
            m_b = h_b - e_b
            m_n = h_n - e_n
            expected = (
                Fraction(1, 2) * samples[b_value - 1]
                + sum(samples[b_value : n_value - 1])
                + Fraction(1, 2) * samples[n_value - 1]
            )
            require(m_n - m_b == expected, "starred extension failed")


def check_cells_and_gaps() -> None:
    cases = 0
    strips = 0
    gaps = 0
    for n_value in range(2, 18):
        for b_value in range(1, n_value):
            for offset in [Fraction(1, 11), Fraction(2, 5), Fraction(7, 9)]:
                alpha = Fraction(n_value * n_value) - offset
                all_cells: list[int] = []
                terminal_cells: list[int] = []
                for q in range(1, n_value + 1):
                    a_q = alpha / (Fraction(q) + Fraction(1, 2))
                    b_q = alpha / (Fraction(q) - Fraction(1, 2))
                    direct = [
                        r
                        for r in range(1, floor_fraction(2 * alpha) + 1)
                        if a_q < r <= b_q
                    ]
                    lower = floor_fraction(a_q) + 1
                    upper = floor_fraction(b_q)
                    formula = list(range(lower, upper + 1)) if lower <= upper else []
                    require(direct == formula, "independent cell formula failed")
                    all_cells.extend(direct)
                    if q > b_value:
                        terminal_cells.extend(direct)
                    cases += 1

                a_n = alpha / (Fraction(n_value) + Fraction(1, 2))
                a_b = alpha / (Fraction(b_value) + Fraction(1, 2))
                expected_all = [
                    r
                    for r in range(1, floor_fraction(2 * alpha) + 1)
                    if a_n < r <= 2 * alpha
                ]
                expected_terminal = [
                    r
                    for r in range(1, floor_fraction(a_b) + 1)
                    if a_n < r <= a_b
                ]
                require(sorted(all_cells) == expected_all, "independent full tiling failed")
                require(
                    sorted(terminal_cells) == expected_terminal,
                    "independent strip failed",
                )
                require(len(all_cells) == len(set(all_cells)), "independent overlap")
                strips += 1

                low_gap = alpha / n_value - a_n
                require(low_gap > Fraction(1, 4), "independent low gap failed")
                high_gap = (floor_fraction(2 * alpha) + 1) - alpha
                require(high_gap > alpha, "independent high gap failed")
                gaps += 1
    require(cases > 4500, "insufficient independent cell cases")
    require(strips > 300, "insufficient independent terminal strips")
    require(gaps == strips, "gap audit count drifted")


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["rth_25_handoff"],
        "open-row boundary drifted",
    )
    require(rows[-1]["readiness"] == "guard_validated", "final guard drifted")


def check_note(note: str) -> None:
    markers = [
        "This is not a proof of RH",
        "## Current Reassembly",
        "## Full-Support Observation Image",
        "## Terminal Interval",
        "## Reciprocal Terminal Strip",
        "## Outer Complement",
        "only the genuine endpoint edge",
        "tilde(Psi)(T_0)",
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
    require(artifact.get("counts") == EXPECTED_COUNTS, "count block drifted")
    require("signed-current bounds open" in artifact.get("status", ""), "status drifted")
    require("no bound for the edge-composed outer complement" in artifact.get("proof_boundary", ""), "proof boundary drifted")
    cert = artifact.get("symbolic_certificate", {})
    require(
        set(cert)
        == {
            "current_homotopy",
            "full_support_observations",
            "terminal_extension",
            "reciprocal_full_support",
            "outer_complement",
            "audit",
        },
        "certificate keys drifted",
    )
    audit = cert["audit"]
    require(audit["cell_cases"] == 828, "builder cell count drifted")
    require(audit["terminal_strip_tilings"] == 90, "builder strip count drifted")
    require(audit["full_support_tilings"] == 90, "builder tiling count drifted")
    require(audit["roster_containments"] == 90, "builder roster count drifted")
    require(audit["low_gap_cases"] == 90, "builder gap count drifted")
    check_sources(artifact)
    check_current_algebra()
    check_matrix_reduction()
    check_terminal_extension()
    check_cells_and_gaps()
    check_rows(artifact)
    check_note(note)
    counts = artifact["counts"]
    print(
        "validated reciprocal terminal full-support homotopy gate: "
        f"{counts['rows']} rows, "
        f"{counts['all_carrier_homotopies']} all-carrier homotopy, "
        f"{counts['terminal_extension_identities']} terminal extensions, "
        f"{counts['terminal_reciprocal_strip_tilings']} terminal strip, "
        f"{counts['outer_complement_bounds']} outer-complement bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
