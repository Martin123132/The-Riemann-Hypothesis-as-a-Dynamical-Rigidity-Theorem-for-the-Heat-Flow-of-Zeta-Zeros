#!/usr/bin/env python3
"""Independently replay the pole-safe P_W/unowned-cell contour join."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp
from sympy import Rational, euler


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_PW_unowned_cell_pole_safe_two_boundary_contour_join_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_direct(z: mp.mpc, alpha_min: int, count: int) -> mp.mpc:
    return mp.fsum(mp.exp(1j * mp.pi * (alpha_min + 2 * j) * z) for j in range(count))


def source_paired(z: mp.mpc, n_minus: int, n_plus: int) -> mp.mpc:
    return (
        0.5j
        * (mp.exp(2j * mp.pi * n_minus * z) - mp.exp(2j * mp.pi * n_plus * z))
        / mp.sin(mp.pi * z)
    )


def regularized_packet_discrepancy(z: mp.mpc, s: mp.mpc, terms: int) -> mp.mpc:
    common = mp.power(z, -s) * mp.exp(-1j * mp.pi * z * z) / (2j * mp.sin(mp.pi * z))
    csc_form = common * (mp.exp(2j * mp.pi * terms * z) - 1)
    finite_form = (
        mp.power(z, -s)
        * mp.exp(-1j * mp.pi * z * z)
        * mp.fsum(mp.exp(1j * mp.pi * (2 * j + 1) * z) for j in range(terms))
    )
    return csc_form - finite_form


def euler_alternating_moment(alpha_min: int, count: int, order: int) -> int:
    x = Rational(alpha_min, 2)
    value = Rational(2) ** (order - 1) * (euler(order, x) - ((-1) ** count) * euler(order, x + count))
    require(value.q == 1, "fresh Euler moment is not integral")
    return int(value)


def closed_label(t: mp.mpf, alpha: int) -> tuple[mp.mpc, mp.mpc]:
    s_a = mp.mpf("0.25") - 0.5j * t
    s_b = mp.mpf("0.75") - 0.5j * t
    xi = 1j * mp.pi * alpha * alpha / 4
    c = mp.pi * alpha * (1 + 1j) / mp.sqrt(2)
    d = mp.exp(-mp.pi * t / 4 - 1j * mp.pi / 8)
    a_term = mp.gamma(s_a) * mp.pi ** (-s_a) * mp.hyp1f1(s_a, mp.mpf("0.5"), xi) / 2
    b_term = c * mp.gamma(s_b) * mp.pi ** (-s_b) * mp.hyp1f1(s_b, mp.mpf("1.5"), xi) / 2
    return 1j * mp.exp(mp.pi * t) * d * (a_term - b_term), d * (a_term + b_term)


def contour_prefix(t: mp.mpf, alpha_min: int, count: int) -> mp.mpc:
    s = mp.mpf("0.5") + 1j * t
    direction = (1 - 1j) / mp.sqrt(2)

    def integrand(q: mp.mpf) -> mp.mpc:
        z = mp.mpf("0.5") + q * direction
        return -mp.power(z, -s) * mp.exp(-1j * mp.pi * z * z) * source_direct(z, alpha_min, count) * direction

    return mp.quad(integrand, [-mp.inf, -7, -2, 0, 2, 5, 8, 12, 18, 25, mp.inf])


def integrated_lower_series(
    t: mp.mpf,
    alpha_min: int,
    count: int,
    u_cut: mp.mpf,
    terms: int,
) -> mp.mpc:
    s = mp.mpf("0.5") + 1j * t
    x_cut = mp.exp(u_cut)
    lam = mp.pi * (1 + 1j) / mp.sqrt(2)
    total = mp.mpc(0)
    for alpha in range(alpha_min, alpha_min + 2 * count, 2):
        linear = -lam * alpha
        quadratic = -mp.pi
        coefficients = [mp.mpc(1)]
        for order in range(1, terms):
            previous = coefficients[order - 1]
            previous_two = coefficients[order - 2] if order > 1 else mp.mpc(0)
            coefficients.append((linear * previous + 2 * quadratic * previous_two) / order)
        total += mp.fsum(
            coefficient * mp.power(x_cut, order + 1 - s) / (order + 1 - s)
            for order, coefficient in enumerate(coefficients)
        )
    return total


def rotated_negative(t: mp.mpf, alpha_min: int, count: int) -> tuple[mp.mpc, mp.mpf, mp.mpf]:
    s = mp.mpf("0.5") + 1j * t
    lam = mp.pi * (1 + 1j) / mp.sqrt(2)
    prefactor = mp.exp(3 * mp.pi * t / 4 + 3j * mp.pi / 8)
    u_lower = mp.mpf("-9.5")
    u_upper = mp.mpf("2.25")

    def integrand(u: mp.mpf) -> mp.mpc:
        x = mp.exp(u)
        labels = mp.fsum(mp.exp(-lam * (alpha_min + 2 * j) * x) for j in range(count))
        return mp.exp((1 - s) * u) * mp.exp(-mp.pi * x * x) * labels

    lower_38 = integrated_lower_series(t, alpha_min, count, u_lower, 38)
    lower_46 = integrated_lower_series(t, alpha_min, count, u_lower, 46)
    lower_delta = abs(prefactor * (lower_46 - lower_38))
    middle = mp.quad(integrand, [u_lower, -6.25, -3.75, -1.5, 0.25, u_upper])

    x_upper = mp.exp(u_upper)
    decay_rate = mp.pi / mp.sqrt(2)
    exponent = mp.pi * x_upper * x_upper + decay_rate * alpha_min * x_upper
    derivative_floor = 2 * mp.pi * x_upper + decay_rate * alpha_min
    upper_bound = (
        abs(prefactor)
        * count
        * mp.power(x_upper, mp.mpf("-0.5"))
        * mp.exp(-exponent)
        / derivative_floor
    )
    return prefactor * (lower_46 + middle), lower_delta, upper_bound


def finite_interior(t: mp.mpf, alpha_min: int, count: int, lower: mp.mpf, upper: mp.mpf) -> mp.mpc:
    s = mp.mpf("0.5") + 1j * t

    def integrand(y: mp.mpf) -> mp.mpc:
        return mp.power(y, -s) * mp.exp(-1j * mp.pi * y * y) * source_direct(y, alpha_min, count)

    points = [lower]
    cursor = lower + mp.mpf("0.2")
    while cursor < upper:
        points.append(cursor)
        cursor += mp.mpf("0.2")
    points.append(upper)
    return mp.quad(integrand, points)


def main() -> int:
    priority = gate.set_low_priority()
    require(priority == "below_normal_one_cpu", f"checker resource cap unavailable: {priority}")
    require(gate.RESULT.is_file(), "missing two-boundary result")
    require(gate.NOTE.is_file(), "missing two-boundary note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("kind") == gate.STEM, "kind drift")
    require(
        artifact.get("status") == "exact_PW_unowned_cells_pole_safe_two_boundary_contour_join_certified",
        "status drift",
    )
    require(artifact.get("passed") is True, "gate did not pass")

    for name, path in gate.DEPENDENCIES.items():
        require(file_hash(path) == artifact["dependencies"][name]["sha256"], f"dependency hash drift: {name}")

    roster = artifact["actual_roster"]
    require(roster["alpha_min"] == 2 * roster["n_minus"] + 1, "lower index map drift")
    require(roster["alpha_max"] == 2 * roster["n_plus"] - 1, "upper index map drift")
    require(roster["n_plus"] - roster["n_minus"] == roster["label_count"], "count drift")

    mp.mp.dps = 82
    for z in (mp.mpc("0.19", "0.33"), mp.mpc("1.61", "-0.11"), mp.mpc("3.5", "0.07")):
        discrepancy = abs(source_direct(z, 5, 5) - source_paired(z, 2, 7))
        require(discrepancy < mp.mpf("1e-71"), "fresh paired-source identity failed")
        for terms in (1, 3, 6, 10):
            require(
                abs(regularized_packet_discrepancy(z, mp.mpf("0.5") + mp.mpf("2.75") * 1j, terms))
                < mp.mpf("1e-70"),
                "fresh common-G0 finite-prefix regularization failed",
            )

    s = mp.mpf("0.5") + mp.mpf("5.75") * 1j
    for k in (1, 2, 5):
        expected_residue = mp.power(k, -s) / (2j * mp.pi)
        residues = []
        for n in (4, 13):
            eps = mp.mpf("1e-25")
            z = mp.mpf(k) + eps
            g_n = mp.power(z, -s) * mp.exp(-1j * mp.pi * z * z + 2j * mp.pi * n * z) / (2j * mp.sin(mp.pi * z))
            residues.append(eps * g_n)
        require(abs(residues[0] - expected_residue) < mp.mpf("1e-23"), "first residue limit failed")
        require(abs(residues[1] - expected_residue) < mp.mpf("1e-23"), "second residue limit failed")
        require(abs(residues[0] - residues[1]) < mp.mpf("1e-23"), "common residue cancellation failed")

    for k in (-3, -1, 0, 2, 6):
        value = source_direct(mp.mpf(k) + mp.mpf("0.5"), 9, 9)
        target = 1j * (1 if k % 2 == 0 else -1)
        require(abs(value - target) < mp.mpf("1e-68"), "fresh half-integer endpoint phase failed")

    for order in range(9):
        direct_moment = sum(((-1) ** j) * (5 + 2 * j) ** order for j in range(5))
        require(direct_moment == euler_alternating_moment(5, 5, order), "fresh Euler endpoint jet failed")
        recorded = artifact["half_integer_endpoints"]["production_alternating_moments_orders_0_through_8"][order]
        require(recorded["order"] == order, "recorded endpoint-jet order drift")
        require(
            int(recorded["alternating_power_sum"])
            == euler_alternating_moment(roster["alpha_min"], roster["label_count"], order),
            "recorded production Euler moment drift",
        )

    t = mp.mpf("4.5")
    alpha_min = 5
    alpha_max = 13
    count = 5
    lower = mp.mpf("1.5")
    upper = mp.mpf("3.5")
    source_closed = mp.mpc(0)
    correction_closed = mp.mpc(0)
    for alpha in range(alpha_min, alpha_max + 1, 2):
        source, correction = closed_label(t, alpha)
        source_closed += source
        correction_closed += correction

    source_quad, lower_delta, upper_bound = rotated_negative(t, alpha_min, count)
    prefix_quad = contour_prefix(t, alpha_min, count)
    interior_quad = finite_interior(t, alpha_min, count, lower, upper)
    joined_source = source_closed - interior_quad
    joined_contour = prefix_quad + correction_closed - interior_quad
    require(abs(source_closed - source_quad) < mp.mpf("1e-54"), "fresh rotated-negative replay failed")
    require(lower_delta < mp.mpf("1e-66"), "fresh lower-series stabilization failed")
    require(upper_bound < mp.mpf("1e-90"), "fresh upper-tail bound failed")
    require(abs(source_closed - prefix_quad - correction_closed) < mp.mpf("1e-54"), "fresh branch identity failed")
    require(abs(joined_source - joined_contour) < mp.mpf("1e-54"), "fresh two-boundary join failed")

    regularized_j = {}
    for terms in (2, 7):
        prefix_source = mp.fsum(closed_label(t, alpha)[0] for alpha in range(1, 2 * terms, 2))
        prefix_interior = finite_interior(t, 1, terms, lower, upper)
        regularized_j[terms] = prefix_interior - prefix_source
    require(
        abs(joined_source - (regularized_j[2] - regularized_j[7])) < mp.mpf("1e-54"),
        "fresh regularized J_N join failed",
    )

    decision = artifact["decision"]
    require(decision["P_W_joined_exactly_to_lower_and_upper_unowned_cells"] is True, "join decision lost")
    require(decision["paired_integer_pole_cancellation_proved"] is True, "pole decision lost")
    require(decision["separate_csc_tail_packets_rejected"] is True, "unsafe split guard lost")
    require(decision["common_G0_regularization_proved"] is True, "common regularization lost")
    require(decision["regularized_finite_prefix_J_N_join_proved"] is True, "regularized J_N join lost")
    require(decision["complete_three_packet_fresnel_cancellation_proved"] is True, "packet cancellation lost")
    require(decision["quantitative_endpoint_complete_Mordell_compression_proved"] is False, "Mordell overclaim")
    require(decision["joined_unowned_packet_enclosed_at_actual_height"] is False, "actual packet overclaim")
    require(decision["D_K_enclosed"] is False and decision["rh_implication"] is False, "conclusion overclaim")
    require(
        artifact["half_integer_endpoints"]["constant_work_full_integrand_endpoint_jets"] is True,
        "full endpoint-jet recurrence lost",
    )

    note = gate.NOTE.read_text(encoding="utf-8")
    for phrase in (
        "canonical common subtraction",
        "quantitative Mordell compression open",
        "must be assembled before a",
        "No fitted circle constant",
        "RH, or a prize-level conclusion",
    ):
        require(phrase in note, f"missing note guard: {phrase}")

    print("independently checked pole-safe P_W/unowned-cell two-boundary contour join", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
