#!/usr/bin/env python3
"""Independently check the ideal-cubic symmetric outer-pairing gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_ideal_cubic_"
    "symmetric_outer_pairing_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 34,
    "recombination_identities": 5,
    "finite_roster_checks": 88,
    "character_split_checks": 7,
    "near_kernel_mean_checks": 8,
    "remote_pair_formula_checks": 5,
    "remote_tail_bound_checks": 10,
    "ideal_endpoint_identities": 6,
    "one_sided_divergence_witnesses": 1,
    "reflection_identities": 2,
    "conjugation_checks": 4,
    "fourier_mode_checks": 14,
    "fourier_mode_witnesses": 2,
    "signed_complete_ideal_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [f"sop_{index:02d}_{suffix}" for index, suffix in enumerate([
    "sources", "poisson", "exterior", "hermitian_join", "transpose_join",
    "roster", "finite_split", "near", "remote", "near_kernel", "kernel_mean",
    "pair_cosine", "pair_ibp", "pair_bound", "scaling_guard", "endpoint_ibp",
    "split_guard", "h_terminal", "h_lower", "h_delta", "h_positive", "h_negative",
    "h_pair", "linear_witness", "cubic_reflection", "coordinate_reflection",
    "phase_reflection", "conjugation", "abel", "fourier_inside", "fourier_outside",
    "near_open", "terminal_open", "boundary",
], start=1)]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(sp.factor(sp.together(expression)))
    if value != 0 and value.equals(0) is not True:
        raise RuntimeError(f"{label}: {value}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(artifact: dict) -> None:
    expected = {
        "finite_cell_inversion",
        "outer_endpoint_kernel",
        "joined_recombination",
        "ridge_carrier",
    }
    source = artifact.get("source_audit", {})
    require(set(source) == expected, "source keys drifted")
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def check_recombination() -> int:
    m, f, o, d, terminal = sp.symbols("m f o d terminal")
    o_value = m - f
    d_value = f - m
    expressions = [
        f + o_value - m,
        o_value + d_value,
        m - d_value - (m + o_value),
        m + o_value + terminal - (2 * m - f + terminal),
        (m + o_value) - (2 * m - f),
    ]
    for index, expression in enumerate(expressions):
        require_zero(expression, f"independent recombination {index}")
    return len(expressions)


def check_cutoff() -> tuple[int, int]:
    roster_checks = 0
    for m_value in range(1, 6):
        for n_value in range(m_value, 10):
            for radius in [n_value, n_value + 2, n_value + 5]:
                direct = {
                    r
                    for r in range(-radius, radius + 1)
                    if not m_value <= r <= n_value
                }
                near = set(range(1, m_value)) | {0}
                near |= {-k for k in range(1, n_value + 1)}
                remote = {
                    sign * k
                    for k in range(n_value + 1, radius + 1)
                    for sign in (-1, 1)
                }
                require(not near & remote, "independent split overlap")
                require(direct == near | remote, "independent cutoff split")
                roster_checks += 1

    mp.mp.dps = 70
    mean_checks = 0
    for m_value, n_value in [(1, 3), (2, 6), (4, 9)]:
        for cell in [-2, 0, 3]:
            mean = mp.quad(
                lambda u: mp.fsum(
                    mp.e ** (2j * mp.pi * j * u)
                    for j in range(-(m_value - 1), n_value + 1)
                ),
                [cell, cell + 1],
            )
            require(abs(mean - 1) < mp.mpf("1e-58"), "independent kernel mean")
            mean_checks += 1
    return roster_checks, mean_checks


def check_remote_pairs() -> int:
    mp.mp.dps = 72
    checks = 0
    terminal = mp.mpf(5)
    for k in [6, 8, 11, 14, 19, 23]:
        direct = 2 * mp.quad(
            lambda u: u * u * mp.cos(2 * mp.pi * k * u),
            [1, 2, 3, 4, 5],
        )
        expected = (terminal - 1) / (mp.pi**2 * k * k)
        require(abs(direct - expected) < mp.mpf("1e-58"), "independent remote pair")
        checks += 1
    return checks


def check_endpoints() -> int:
    x, a_log, n_log, ux = sp.symbols("x a_log n_log ux", real=True)
    u_n = a_log - n_log
    p_h = -sp.I * (x + u_n) * (x - u_n) ** 2 / 4
    p_t = -sp.I * (x - u_n) * (x + u_n) ** 2 / 4 + ux * (x - u_n)
    c_h = n_log * (2 * a_log - n_log) ** 2 / 4
    delta_h = -sp.I * c_h
    kappa = 1 / (2 * sp.pi * sp.I)
    expressions = [
        p_h.subs(x, -u_n),
        p_h.subs(x, -a_log) - sp.I * c_h,
        -p_h.subs(x, -a_log) - delta_h,
        -kappa * delta_h - c_h / (2 * sp.pi),
        kappa * delta_h + c_h / (2 * sp.pi),
        p_t + p_h.subs(x, -x) - ux * (x - u_n),
    ]
    for index, expression in enumerate(expressions):
        require_zero(expression, f"independent endpoint {index}")
    return len(expressions)


def check_reflected_phase() -> int:
    v, alpha, r, a_value, candidate = sp.symbols(
        "v alpha r a candidate", nonzero=True
    )
    reflected_derivative = -alpha / v + r * a_value**2 / v**2
    ordinary_derivative = alpha / v - candidate
    polynomial = sp.Poly(
        sp.expand(v**2 * (reflected_derivative - ordinary_derivative)), v
    )
    require(polynomial.coeff_monomial(v) == -2 * alpha, "reflected linear coefficient")
    require(polynomial.coeff_monomial(v**2) == candidate, "reflected quadratic coefficient")
    return 2


def check_conjugation() -> int:
    mp.mp.dps = 65
    alpha = mp.mpf("3.17")

    def s_value(log_u):
        return -mp.mpf("0.07") * log_u**2

    def p_value(log_u):
        return (1 - mp.mpf("0.3") * 1j) + (mp.mpf("0.8") + mp.mpf("0.2") * 1j) * log_u

    checks = 0
    for r_value in [-4, 1, 6]:
        left = mp.conj(mp.quad(
            lambda u: mp.e ** (2j * mp.pi * alpha * mp.log(u))
            * mp.e ** s_value(mp.log(u)) * p_value(mp.log(u))
            * mp.e ** (-2j * mp.pi * r_value * u),
            [1, 2, 3],
        ))
        right = mp.quad(
            lambda u: mp.e ** (-2j * mp.pi * alpha * mp.log(u))
            * mp.e ** s_value(mp.log(u)) * mp.conj(p_value(mp.log(u)))
            * mp.e ** (2j * mp.pi * r_value * u),
            [1, 2, 3],
        )
        require(abs(left - right) < mp.mpf("1e-50"), "independent conjugation")
        checks += 1
    return checks


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["sop_32_near_open", "sop_33_terminal_open"],
        "open rows drifted",
    )


def check_note(note: str) -> None:
    markers = [
        "not a proof of RH",
        "## Exterior Pullback",
        "## Canonical Symmetric Split",
        "K_near,N",
        "## Absolutely Convergent Remote Pairs",
        "2pi k",
        "## One-Sided Obstruction",
        "Delta_H=-i c_H",
        "## Symmetry Audit",
        "a^2/u",
        "## Coefficient-Blind Guard",
        "## Route Decision",
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
    require("one-sided-divergence obstruction" in artifact.get("status", ""), "status drifted")
    require("no signed complete ideal join" in artifact.get("proof_boundary", ""), "boundary drifted")
    require(
        set(artifact.get("symbolic_certificate", {}))
        == {"recombination", "cutoff", "paired_tail", "endpoint", "symmetry", "fourier_guard", "handoff"},
        "certificate keys drifted",
    )
    check_sources(artifact)
    recombinations = check_recombination()
    roster_checks, mean_checks = check_cutoff()
    pair_checks = check_remote_pairs()
    endpoint_checks = check_endpoints()
    reflection_checks = check_reflected_phase()
    conjugation_checks = check_conjugation()
    check_rows(artifact)
    check_note(note)
    print(
        "validated ideal-cubic symmetric outer-pairing gate: "
        f"{EXPECTED_COUNTS['rows']} rows, {recombinations} recombinations, "
        f"{roster_checks} independent roster checks, {mean_checks} independent means, "
        f"{pair_checks} independent remote pairs, {endpoint_checks} endpoint identities, "
        f"{reflection_checks} reflected-phase checks, {conjugation_checks} conjugation checks, "
        f"{EXPECTED_COUNTS['signed_complete_ideal_bounds']} signed complete ideal bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
