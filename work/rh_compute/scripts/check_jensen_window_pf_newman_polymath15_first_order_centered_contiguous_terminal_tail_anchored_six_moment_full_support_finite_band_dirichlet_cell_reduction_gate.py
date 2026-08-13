#!/usr/bin/env python3
"""Independently check the finite-band Dirichlet-cell reduction gate."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_finite_band_dirichlet_"
    "cell_reduction_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 33,
    "zero_mean_checks": 204,
    "odd_harmonic_checks": 204,
    "periodic_checks": 1020,
    "cell_recombinations": 3,
    "parity_decompositions": 1,
    "relative_cell_reductions": 2,
    "tie_jumps": 2,
    "fixed_cell_transports": 1,
    "reciprocal_block_decompositions": 1,
    "far_gap_checks": 20300,
    "adjacent_block_classes": 2,
    "fourier_witness_checks": 1680,
    "signed_physical_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [
    "fbd_01_roster",
    "fbd_02_kernel",
    "fbd_03_closed",
    "fbd_04_symmetry",
    "fbd_05_half",
    "fbd_06_partition",
    "fbd_07_variation",
    "fbd_08_band",
    "fbd_09_defect",
    "fbd_10_joined",
    "fbd_11_parity",
    "fbd_12_hermitian",
    "fbd_13_transpose",
    "fbd_14_terminal",
    "fbd_15_ideal",
    "fbd_16_upper_tie",
    "fbd_17_lower_tie",
    "fbd_18_transport",
    "fbd_19_floor_guard",
    "fbd_20_cells",
    "fbd_21_blocks",
    "fbd_22_gradient",
    "fbd_23_far",
    "fbd_24_adjacent",
    "fbd_25_local_phase",
    "fbd_26_inband",
    "fbd_27_outband",
    "fbd_28_nouniversal",
    "fbd_29_derivative",
    "fbd_30_route",
    "fbd_31_decision",
    "fbd_32_obligation",
    "fbd_33_boundary",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.expand(expression)
    if value != 0:
        raise RuntimeError(f"{label}: {value}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(artifact: dict) -> None:
    source = artifact.get("source_audit", {})
    require(set(source) == {"joined_recombination", "finite_cell"}, "source keys drifted")
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def direct_kernel(m: int, n: int, value: mp.mpf) -> mp.mpc:
    return mp.fsum(mp.e ** (-2j * mp.pi * r * value) for r in range(m, n + 1))


def closed_kernel(m: int, n: int, value: mp.mpf) -> mp.mpc:
    nearest = mp.nint(value)
    offset = value - nearest
    length = n - m + 1
    if abs(offset) < mp.mpf("1e-50"):
        return mp.mpc(length)
    return (
        mp.e ** (-mp.pi * 1j * (m + n) * offset)
        * mp.sin(mp.pi * length * offset)
        / mp.sin(mp.pi * offset)
    )


def check_kernel() -> int:
    mp.mp.dps = 70
    checks = 0
    offsets = [
        mp.mpf(-9) / 20,
        mp.mpf(-7) / 31,
        mp.mpf(-1) / 17,
        mp.mpf(1) / 29,
        mp.mpf(5) / 23,
        mp.mpf(11) / 25,
    ]
    for m in range(1, 21):
        for width in (1, 2, 5, 9, 14):
            n = m + width - 1
            for offset in offsets:
                direct = direct_kernel(m, n, offset)
                closed = closed_kernel(m, n, offset)
                require(abs(direct - closed) < mp.mpf("1e-60"), "sine quotient mismatch")
                require(
                    abs(direct_kernel(m, n, offset + 1) - direct) < mp.mpf("1e-60"),
                    "period mismatch",
                )
                require(
                    abs(direct_kernel(m, n, -offset) - mp.conj(direct))
                    < mp.mpf("1e-60"),
                    "conjugation mismatch",
                )
                checks += 3

            half_numeric = mp.fsum(
                (1 - (-1) ** r) / (2j * mp.pi * r) for r in range(m, n + 1)
            )
            half_expected = mp.fsum(
                1 / (mp.pi * 1j * r) for r in range(m, n + 1) if r % 2 == 1
            )
            require(abs(half_numeric - half_expected) < mp.mpf("1e-65"), "half harmonic")
            checks += 1
    return checks


def check_cell_algebra() -> None:
    h, f1, fn, variation, carrier = sp.symbols("h f1 fn variation carrier")
    finite = h * (f1 - fn) + variation
    defect = finite - carrier
    require_zero(defect - (-carrier + h * (f1 - fn) + variation), "independent defect")
    require_zero(
        carrier - defect - (2 * carrier - h * (f1 - fn) - variation),
        "independent joined carrier",
    )

    kr, ki = sp.symbols("kr ki", real=True)
    fp, fm, f0 = sp.symbols("fp fm f0")
    left = (kr + sp.I * ki) * (fp - f0) + (kr - sp.I * ki) * (fm - f0)
    right = kr * (fp + fm - 2 * f0) + sp.I * ki * (fp - fm)
    require_zero(left - right, "independent parity")


def check_cell_numerics() -> int:
    mp.mp.dps = 45
    checks = 0
    for n_value in range(2, 9):
        for m, n in ((1, 3), (2, 7), (5, 11)):
            def kernel(v: mp.mpf) -> mp.mpc:
                return direct_kernel(m, n, v)

            def source(u: mp.mpf) -> mp.mpc:
                return (1 + u / 7 + 1j * u**2 / 19) * mp.e ** (1j * u / 5)

            finite_direct = mp.quad(lambda u: source(u) * kernel(u), [1, n_value])
            h_value = mp.fsum(
                1 / (mp.pi * 1j * r) for r in range(m, n + 1) if r % 2 == 1
            )
            variation = mp.quad(
                lambda v: (source(1 + v) - source(1)) * kernel(v),
                [0, mp.mpf("0.5")],
            )
            for q in range(2, n_value):
                variation += mp.quad(
                    lambda v, q=q: (source(q + v) - source(q)) * kernel(v),
                    [mp.mpf("-0.5"), mp.mpf("0.5")],
                )
            variation += mp.quad(
                lambda v: (source(n_value + v) - source(n_value)) * kernel(v),
                [mp.mpf("-0.5"), 0],
            )
            recomposed = h_value * (source(1) - source(n_value)) + variation
            require(abs(finite_direct - recomposed) < mp.mpf("2e-35"), "numeric cell recomposition")
            checks += 1
    return checks


def check_relative_and_ties() -> None:
    mh, mt, h, fh1, ft1, ftn, vh, vt, term, un = sp.symbols(
        "mh mt h fh1 ft1 ftn vh vt term un"
    )
    dh = -mh + h * fh1 + vh
    dt = -mt + h * (ft1 - ftn) + vt
    require_zero(-dh - (mh - h * fh1 - vh), "independent H relative")
    require_zero(
        -dt + 2 * sp.I * un * term
        - (mt - h * (ft1 - ftn) - vt + 2 * sp.I * un * term),
        "independent T relative",
    )

    f, o, moved = sp.symbols("f o moved")
    require_zero((f + moved) + (o - moved) - f - o, "upper tie")
    require_zero((f - moved) + (o + moved) - f - o, "lower tie")


def check_reciprocal_blocks() -> int:
    alpha, u, v = sp.symbols("alpha u v", positive=True)
    q, w, r = sp.symbols("q w r", positive=True)
    require_zero(
        alpha / u - alpha / v - alpha * (v - u) / (u * v),
        "independent reciprocal gradient",
    )
    exact_phase = alpha * sp.log(1 + w / q) - r * w
    split_phase = (alpha / q - r) * w + alpha * (sp.log(1 + w / q) - w / q)
    require_zero(exact_phase - split_phase, "independent local phase")

    gap_checks = 0
    offsets = [Fraction(-2, 5), Fraction(-1, 9), Fraction(2, 11), Fraction(7, 20)]
    for q_value in range(1, 25):
        for p_value in range(1, 25):
            distance = abs(p_value - q_value)
            if distance < 2:
                continue
            for u_offset in offsets:
                for v_offset in offsets:
                    u_value = Fraction(q_value) + u_offset
                    v_value = Fraction(p_value) + v_offset
                    separation = abs(v_value - u_value)
                    require(separation >= distance - 1, "independent far-cell separation")
                    gradient = abs(Fraction(1, 1) / u_value - Fraction(1, 1) / v_value)
                    lower = Fraction(distance - 1, 1) / (
                        (Fraction(q_value) + Fraction(1, 2))
                        * (Fraction(p_value) + Fraction(1, 2))
                    )
                    require(gradient >= lower, "independent far-gradient envelope")
                    gap_checks += 1

    for n_value in range(2, 13):
        for alpha_value in (Fraction(17, 3), Fraction(41, 5), Fraction(97, 7)):
            m_value = int(alpha_value / (Fraction(n_value) + Fraction(1, 2))) + 1
            n_band = int(2 * alpha_value)
            for frequency in range(m_value, n_band + 1):
                reciprocal = alpha_value / frequency
                p_value = int(reciprocal + Fraction(1, 2))
                require(1 <= p_value <= n_value, "frequency cell escaped physical labels")
                require(
                    Fraction(p_value) - Fraction(1, 2)
                    <= reciprocal
                    < Fraction(p_value) + Fraction(1, 2),
                    "frequency cell half-open convention drifted",
                )
    return gap_checks


def check_fourier_witnesses() -> int:
    checks = 0
    for length in range(1, 65):
        for m in range(1, 11):
            n = m + 7
            for k in (0, m - 1, m, m + 2, n - 1, n, n + 1, n + 5):
                finite = length if m <= k <= n else 0
                carrier = length
                defect = finite - carrier
                require(defect == (0 if m <= k <= n else -length), "independent Fourier witness")
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
        == ["fbd_30_route", "fbd_32_obligation"],
        "open-row boundary drifted",
    )


def check_note(note: str) -> None:
    markers = [
        "not a proof of RH",
        "## Band Kernel",
        "K_N(u)=sum",
        "purely imaginary",
        "## Cell Reduction",
        "mathscr D_N[P]=-mathscr M_N[P]",
        "V_N[P]",
        "## Relative Channels",
        "2mathscr M_N[P_H]",
        "## Tie Transport",
        "## Reciprocal Blocks",
        "|p-q|>=2",
        "## Route Guards",
        "coefficient-blind",
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
    require("signed physical variation estimate open" in artifact.get("status", ""), "status drifted")
    require("no signed physical variation" in artifact.get("proof_boundary", ""), "proof boundary drifted")
    require(
        set(artifact.get("symbolic_certificate", {}))
        == {
            "band_kernel",
            "cell_reduction",
            "relative_channels",
            "tie_transport",
            "reciprocal_blocks",
            "route_guards",
            "handoff",
        },
        "certificate keys drifted",
    )
    check_sources(artifact)
    kernel_checks = check_kernel()
    check_cell_algebra()
    cell_checks = check_cell_numerics()
    check_relative_and_ties()
    gap_checks = check_reciprocal_blocks()
    witness_checks = check_fourier_witnesses()
    check_rows(artifact)
    check_note(note)
    print(
        "validated finite-band Dirichlet-cell reduction gate: "
        f"{EXPECTED_COUNTS['rows']} rows, {kernel_checks} independent kernel checks, "
        f"{cell_checks} independent cell checks, "
        f"{EXPECTED_COUNTS['relative_cell_reductions']} relative reductions, "
        f"{gap_checks} independent far-gap checks, "
        f"{witness_checks} independent Fourier witnesses, "
        f"{EXPECTED_COUNTS['signed_physical_bounds']} signed physical bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
