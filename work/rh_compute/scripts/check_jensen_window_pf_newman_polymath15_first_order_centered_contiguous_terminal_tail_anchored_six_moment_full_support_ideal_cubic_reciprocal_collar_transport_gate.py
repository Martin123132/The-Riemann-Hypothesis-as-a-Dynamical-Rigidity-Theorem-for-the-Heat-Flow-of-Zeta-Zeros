#!/usr/bin/env python3
"""Independently check the ideal-cubic reciprocal collar transport gate."""

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
    "terminal_tail_anchored_six_moment_full_support_ideal_cubic_reciprocal_"
    "collar_transport_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 31,
    "collar_decompositions": 107,
    "far_interval_compressions": 1475,
    "internal_tie_checks": 2150,
    "radius_one_tie_checks": 432,
    "noninvariance_witnesses": 1400,
    "fixed_cell_transports": 2,
    "phase_identities": 3,
    "corner_orientation_checks": 752,
    "far_ibp_identities": 1,
    "far_aggregate_bounds": 0,
    "finite_moment_identities": 2,
    "full_line_cancellations": 1,
    "finite_correction_bounds": 0,
    "ideal_cubic_identities": 5,
    "joined_recombinations": 2,
    "signed_collar_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    "rct_01_atom",
    "rct_02_collar",
    "rct_03_radius_one",
    "rct_04_far_compression",
    "rct_05_internal_tie",
    "rct_06_radius_one_tie",
    "rct_07_far_tie",
    "rct_08_witness",
    "rct_09_outer_tie",
    "rct_10_fixed_cell",
    "rct_11_diagonal",
    "rct_12_right_corner",
    "rct_13_left_corner",
    "rct_14_corner_sign",
    "rct_15_far_gap",
    "rct_16_ibp",
    "rct_17_telescope",
    "rct_18_far_bound",
    "rct_19_full_line",
    "rct_20_negative_moment",
    "rct_21_positive_moment",
    "rct_22_boundary_guard",
    "rct_23_finite_correction",
    "rct_24_cubics",
    "rct_25_atoms",
    "rct_26_reflection",
    "rct_27_terminal",
    "rct_28_joined",
    "rct_29_decision",
    "rct_30_obligation",
    "rct_31_boundary",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(sp.expand(expression))
    if value != 0 and value.equals(0) is not True:
        raise RuntimeError(f"{label}: {value}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(artifact: dict) -> None:
    source = artifact.get("source_audit", {})
    require(
        set(source) == {"finite_band", "first_correction", "terminal_quadrature"},
        "source keys drifted",
    )
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def check_collar_combinatorics() -> tuple[int, int]:
    decomposition_checks = 0
    compression_checks = 0
    for n_value in range(2, 22):
        for radius in range(0, n_value):
            atoms = {
                (q, p): Fraction(31 * q * q + 17 * p + 3 * q * p, 127)
                for q in range(1, n_value + 1)
                for p in range(1, n_value + 1)
            }
            full = sum(atoms.values(), Fraction(0))
            near = sum(
                value
                for (q, p), value in atoms.items()
                if abs(q - p) <= radius
            )
            far = sum(
                value
                for (q, p), value in atoms.items()
                if abs(q - p) > radius
            )
            require(full == near + far, "independent collar split")
            decomposition_checks += 1

            for p in range(1, n_value + 1):
                direct = [q for q in range(1, n_value + 1) if abs(q - p) > radius]
                left = list(range(1, max(1, p - radius)))
                right = list(range(min(n_value + 1, p + radius + 1), n_value + 1))
                require(direct == left + right, "independent far compression")
                compression_checks += 1
    return decomposition_checks, compression_checks


def check_tie_transport() -> tuple[int, int]:
    tie_checks = 0
    witness_checks = 0
    for n_value in range(5, 28):
        for radius in range(0, min(6, n_value - 2) + 1):
            for p in range(1, n_value):
                old = {q for q in range(1, n_value + 1) if abs(q - p) <= radius}
                new = {
                    q for q in range(1, n_value + 1) if abs(q - (p + 1)) <= radius
                }
                entering = ({p + radius + 1} if p + radius + 1 <= n_value else set())
                leaving = ({p - radius} if p - radius >= 1 else set())
                require(new - old == entering, "independent entering rail")
                require(old - new == leaving, "independent leaving rail")

                old_far = set(range(1, n_value + 1)) - old
                new_far = set(range(1, n_value + 1)) - new
                require(new_far - old_far == leaving, "independent far entering rail")
                require(old_far - new_far == entering, "independent far leaving rail")
                tie_checks += 1

                if leaving and entering:
                    q_in = next(iter(entering))
                    q_out = next(iter(leaving))
                    # For f(u)=u^2 e(ru), the mode atom is integral_Jq u^2 du.
                    def second_moment(q: int) -> Fraction:
                        left = Fraction(1) if q == 1 else Fraction(q) - Fraction(1, 2)
                        right = (
                            Fraction(n_value)
                            if q == n_value
                            else Fraction(q) + Fraction(1, 2)
                        )
                        return (right**3 - left**3) / 3

                    require(
                        second_moment(q_in) != second_moment(q_out),
                        "independent collar nonclosure witness",
                    )
                    witness_checks += 1
    return tie_checks, witness_checks


def check_phase_charts() -> int:
    mp.mp.dps = 70
    checks = 0
    alphas = [mp.mpf("17.25"), mp.mpf("43.75")]
    offsets = [mp.mpf("-0.4"), mp.mpf("-0.13"), mp.mpf("0.19"), mp.mpf("0.41")]
    for alpha in alphas:
        for q in range(2, 11):
            for z in offsets:
                v = q + z
                r = alpha / v
                for w in offsets:
                    u = q + w
                    y = u / v - 1
                    direct = alpha * mp.log(u) - r * u
                    reduced = alpha * (mp.log(v) - 1) + alpha * (mp.log(1 + y) - y)
                    require(abs(direct - reduced) < mp.mpf("1e-62"), "diagonal phase chart")
                    checks += 1

        corner_offsets = [mp.mpf("0"), mp.mpf("0.07"), mp.mpf("0.21"), mp.mpf("0.46")]
        for q in range(2, 11):
            right_c = mp.mpf(q) + mp.mpf("0.5")
            left_c = mp.mpf(q) - mp.mpf("0.5")
            for a in corner_offsets:
                for b in corner_offsets:
                    right_u = right_c - a
                    right_v = right_c + b
                    right_r = alpha / right_v
                    d_right = (a + b) / right_v
                    direct_right = alpha * mp.log(right_u) - right_r * right_u
                    reduced_right = alpha * (mp.log(right_v) - 1) + alpha * (
                        mp.log(1 - d_right) + d_right
                    )
                    require(
                        abs(direct_right - reduced_right) < mp.mpf("1e-62"),
                        "right corner chart",
                    )

                    left_u = left_c + a
                    left_v = left_c - b
                    left_r = alpha / left_v
                    d_left = (a + b) / left_v
                    direct_left = alpha * mp.log(left_u) - left_r * left_u
                    reduced_left = alpha * (mp.log(left_v) - 1) + alpha * (
                        mp.log(1 + d_left) - d_left
                    )
                    require(
                        abs(direct_left - reduced_left) < mp.mpf("1e-62"),
                        "left corner chart",
                    )
                    checks += 2
    return checks


def check_far_integration_by_parts() -> int:
    mp.mp.dps = 60
    cases = [
        (mp.mpf("12.5"), 3, mp.mpf("1.1"), mp.mpf("2.2")),
        (mp.mpf("12.5"), 8, mp.mpf("2.4"), mp.mpf("3.7")),
        (mp.mpf("30.25"), 4, mp.mpf("1.0"), mp.mpf("3.0")),
        (mp.mpf("30.25"), 20, mp.mpf("2.7"), mp.mpf("4.6")),
    ]

    def amplitude(u: mp.mpf) -> mp.mpc:
        return (1 + u / 7 + 1j * u**2 / 19) * mp.e ** (-u / 23)

    def amplitude_prime(u: mp.mpf) -> mp.mpc:
        polynomial = 1 + u / 7 + 1j * u**2 / 19
        polynomial_prime = mp.mpf(1) / 7 + 2j * u / 19
        return mp.e ** (-u / 23) * (polynomial_prime - polynomial / 23)

    checks = 0
    kappa = 1 / (2 * mp.pi * 1j)
    for alpha, r, lower, upper in cases:
        saddle = alpha / r
        require(not (lower <= saddle <= upper), "test interval contains saddle")

        def phase(u: mp.mpf) -> mp.mpc:
            return mp.e ** (2j * mp.pi * (alpha * mp.log(u) - r * u))

        def gradient(u: mp.mpf) -> mp.mpf:
            return alpha / u - r

        def gradient_prime(u: mp.mpf) -> mp.mpf:
            return -alpha / u**2

        direct = mp.quad(lambda u: amplitude(u) * phase(u), [lower, upper])
        boundary = kappa * (
            amplitude(upper) * phase(upper) / gradient(upper)
            - amplitude(lower) * phase(lower) / gradient(lower)
        )
        remainder = kappa * mp.quad(
            lambda u: phase(u)
            * (
                amplitude_prime(u) / gradient(u)
                - amplitude(u) * gradient_prime(u) / gradient(u) ** 2
            ),
            [lower, upper],
        )
        require(abs(direct - (boundary - remainder)) < mp.mpf("2e-48"), "numeric far IBP")
        checks += 1
    return checks


def check_finite_moments() -> int:
    mp.mp.dps = 60
    intervals = [
        (mp.mpf("-1.3"), mp.mpf("0.7")),
        (mp.mpf("-0.5"), mp.mpf("1.2")),
        (mp.mpf("0.1"), mp.mpf("1.6")),
    ]
    kappa = 1 / (2 * mp.pi * 1j)
    checks = 0
    for lower, upper in intervals:
        minus = lambda y: mp.e ** (-mp.pi * 1j * y**2)
        plus = lambda y: mp.e ** (mp.pi * 1j * y**2)
        minus_second = mp.quad(lambda y: y**2 * minus(y), [lower, upper])
        minus_right = kappa * mp.quad(minus, [lower, upper]) - kappa * (
            upper * minus(upper) - lower * minus(lower)
        )
        require(abs(minus_second - minus_right) < mp.mpf("2e-49"), "negative moment")

        plus_second = mp.quad(lambda y: y**2 * plus(y), [lower, upper])
        plus_right = -kappa * mp.quad(plus, [lower, upper]) + kappa * (
            upper * plus(upper) - lower * plus(lower)
        )
        require(abs(plus_second - plus_right) < mp.mpf("2e-49"), "positive moment")
        checks += 2
    return checks


def check_ideal_cubics() -> None:
    x, u_n, u_x, u_q = sp.symbols("x u_N u_x u_q")
    p_h = -sp.I * (x + u_n) * (x - u_n) ** 2 / 4
    p_t = -sp.I * (x - u_n) * (x + u_n) ** 2 / 4 + u_x * (x - u_n)
    require_zero(p_t + p_h.subs(x, -x) - u_x * (x - u_n), "independent reflection")
    require_zero(p_h.subs(x, -u_n), "independent H terminal")
    require_zero(p_t.subs(x, -u_n) + 2 * u_n * u_x, "independent T terminal")
    require_zero(
        p_h.subs(x, -u_q) - sp.I * (u_q - u_n) * (u_q + u_n) ** 2 / 4,
        "independent H atom",
    )
    carrier, collar, far, terminal = sp.symbols("M C E Z")
    defect = collar + far - carrier
    require_zero(carrier - defect - (2 * carrier - collar - far), "independent H join")
    require_zero(
        carrier - defect + terminal - (2 * carrier - collar - far + terminal),
        "independent T join",
    )


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["rct_18_far_bound", "rct_23_finite_correction", "rct_30_obligation"],
        "open-row boundary drifted",
    )


def check_note(note: str) -> None:
    markers = [
        "not a proof of RH",
        "## Collar Decomposition",
        "mathcal C_(N,1)",
        "## Rail Transport",
        "I_(p+2,r)",
        "2L+1",
        "## Exact Local Phases",
        "log(1-d)+d",
        "log(1+d)-d",
        "## Far Rails",
        "## Finite First Correction",
        "[y e(-y^2/2)]",
        "## Ideal Cubics",
        "mathcal J_H^0=2mathscr M_N",
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
    require("signed tie-complete ideal estimate open" in artifact.get("status", ""), "status drifted")
    require("no signed tie-complete collar" in artifact.get("proof_boundary", ""), "boundary drifted")
    require(
        set(artifact.get("symbolic_certificate", {}))
        == {"collar", "ties", "phases", "far", "moments", "ideal", "handoff"},
        "certificate keys drifted",
    )
    check_sources(artifact)
    collar_checks, compression_checks = check_collar_combinatorics()
    tie_checks, witness_checks = check_tie_transport()
    phase_checks = check_phase_charts()
    ibp_checks = check_far_integration_by_parts()
    moment_checks = check_finite_moments()
    check_ideal_cubics()
    check_rows(artifact)
    check_note(note)
    print(
        "validated ideal-cubic reciprocal collar transport gate: "
        f"{EXPECTED_COUNTS['rows']} rows, {collar_checks} independent collar checks, "
        f"{compression_checks} independent far compressions, "
        f"{tie_checks} independent tie checks, {witness_checks} independent nonclosure witnesses, "
        f"{phase_checks} independent phase checks, {ibp_checks} numerical IBP checks, "
        f"{moment_checks} numerical moment checks, "
        f"{EXPECTED_COUNTS['signed_collar_bounds']} signed collar bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
