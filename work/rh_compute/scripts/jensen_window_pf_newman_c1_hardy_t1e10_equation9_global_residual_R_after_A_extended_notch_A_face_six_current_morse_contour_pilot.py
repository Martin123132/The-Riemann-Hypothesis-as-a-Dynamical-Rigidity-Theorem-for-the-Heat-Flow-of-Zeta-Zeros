#!/usr/bin/env python3
"""Scout the weighted Morse-contour integral of the six-current A defect."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import time


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_morse_contour_pilot"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
DEPENDENCY = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_midpoint_defect_reduction_gate.json"

A = 159_577
T = 10_000_000_000
FIRST_MODE = 39_895
LAST_MODE = 39_936
LEFT = mp.mpf(79_789) / 2
RIGHT = mp.mpf(79_873) / 2
TERMS = 6
Q_CUT = 16
Q_SPLICE_TEXT = "0.0083"
CONTOUR_HEIGHT_TEXT = "0.001"
PRECISION = 80


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def odd_double_factorial(index: int) -> int:
    result = 1
    for value in range(1, index + 1, 2):
        result *= value
    return result


def g_six(y: mp.mpf | mp.mpc, x: mp.mpf | mp.mpc) -> mp.mpc:
    delta = y - mp.mpf(A) * x / 2
    value = mp.mpf(A) * x / (2j * mp.pi * delta)
    for index in range(1, TERMS):
        value += (
            y
            * odd_double_factorial(2 * index - 1)
            * (x / 2) ** index
            / ((1j * mp.pi) ** (index + 1) * delta ** (2 * index + 1))
        )
    return value


def integral_g_six(x: mp.mpf | mp.mpc) -> mp.mpc:
    a = mp.mpf(A) * x / 2
    dl = LEFT - a
    dr = RIGHT - a
    value = mp.mpf(A) * x * mp.log(dr / dl) / (2j * mp.pi)
    for index in range(1, TERMS):
        coefficient = odd_double_factorial(2 * index - 1) * (x / 2) ** index / (1j * mp.pi) ** (index + 1)
        primitive = (
            (dr ** (1 - 2 * index) - dl ** (1 - 2 * index)) / (1 - 2 * index)
            - a * (dr ** (-2 * index) - dl ** (-2 * index)) / (2 * index)
        )
        value += coefficient * primitive
    return value


def defect_six(x: mp.mpf | mp.mpc) -> mp.mpc:
    return mp.fsum(g_six(mp.mpf(mode), x) for mode in range(FIRST_MODE, LAST_MODE + 1)) - integral_g_six(x)


def logistic_x(q: mp.mpf | mp.mpc) -> mp.mpf | mp.mpc:
    return 1 / (1 + mp.exp(q))


def phase(q: mp.mpf | mp.mpc) -> mp.mpf | mp.mpc:
    return mp.pi * A**2 * logistic_x(q) / 4 + mp.mpf(T) * q / 2


def relative_integrand(q: mp.mpf | mp.mpc, phase_star: mp.mpf) -> mp.mpc:
    x = logistic_x(q)
    # E_42=-D_42.  This pilot replaces D_42 by D_6 only on the tail-safe core.
    amplitude = -defect_six(x) * x ** (-mp.mpf(1) / 4) * (1 - x) ** (mp.mpf(3) / 4)
    return amplitude * mp.exp(1j * (phase(q) - phase_star))


def gauss_integral(function, start, stop, panels: int, nodes: int) -> mp.mpc:
    abscissae, weights = mp.gauss_quadrature(nodes, "legendre")
    width = (stop - start) / panels
    total = mp.mpc(0)
    for panel in range(panels):
        left = start + panel * width
        right = left + width
        midpoint = (left + right) / 2
        half_width = width / 2
        total += half_width * mp.fsum(
            weights[index] * function(midpoint + half_width * abscissae[index])
            for index in range(nodes)
        )
    return total


def complex_record(value: mp.mpc, digits: int = 50) -> dict[str, str]:
    return {
        "real": mp.nstr(mp.re(value), digits),
        "imag": mp.nstr(mp.im(value), digits),
        "modulus": mp.nstr(abs(value), digits),
    }


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    mp.mp.dps = PRECISION
    dependency = json.loads(DEPENDENCY.read_text(encoding="utf-8"))
    require(dependency.get("passed") is True, "six-current dependency is not passed")
    require(dependency["decision"]["uniform_signed_defect_error_below_1_point_025e_minus_8_proved"] is True, "six-current remainder drift")
    q_splice = mp.mpf(Q_SPLICE_TEXT)
    contour_height = mp.mpf(CONTOUR_HEIGHT_TEXT)

    sqrt_x_cut = (mp.sqrt(2 * Q_CUT**2 + 8 * A * LEFT) - Q_CUT * mp.sqrt(2)) / (2 * A)
    x_cut = sqrt_x_cut**2
    q_cut = mp.log((1 - x_cut) / x_cut)
    x_star = (1 - mp.sqrt(1 - 8 * T / (mp.pi * A**2))) / 2
    q_star = mp.log((1 - x_star) / x_star)
    phase_star = phase(q_star)

    compact_coarse = gauss_integral(lambda q: relative_integrand(q, phase_star), q_cut, q_splice, 48, 16)
    compact_fine = gauss_integral(lambda q: relative_integrand(q, phase_star), q_cut, q_splice, 96, 24)

    vertical_coarse = 1j * gauss_integral(
        lambda eta: relative_integrand(q_splice + 1j * eta, phase_star),
        mp.mpf(0),
        contour_height,
        16,
        16,
    )
    vertical_fine = 1j * gauss_integral(
        lambda eta: relative_integrand(q_splice + 1j * eta, phase_star),
        mp.mpf(0),
        contour_height,
        32,
        24,
    )

    ray_stop = q_splice + mp.mpf("0.02")
    shifted_ray = gauss_integral(
        lambda real_q: relative_integrand(real_q + 1j * contour_height, phase_star),
        q_splice,
        ray_stop,
        40,
        16,
    )
    full_relative = compact_fine + vertical_fine + shifted_ray
    full_complex = mp.exp(1j * phase_star) * full_relative
    normalization = 2 * (mp.pi / (32 * T)) ** (mp.mpf(1) / 4)
    physical = normalization * mp.re(mp.exp(-1j * mp.pi / 8) * full_complex)

    # Integrable absolute majorant for replacing D_42 by D_6 on the core.
    tail_constant = 2 * odd_double_factorial(11) * LEFT / (mp.pi**7 * 2**6)

    def defect_error_weight(x: mp.mpf) -> mp.mpf:
        delta = LEFT - mp.mpf(A) * x / 2
        defect_error = 84 * tail_constant * x**6 / delta**13
        return defect_error * x ** (-mp.mpf(5) / 4) * (1 - x) ** (-mp.mpf(1) / 4)

    replacement_absolute = normalization * mp.quad(defect_error_weight, [0, x_cut])

    compact_difference = abs(compact_fine - compact_coarse)
    vertical_difference = abs(vertical_fine - vertical_coarse)
    ray_terminal_modulus = abs(relative_integrand(ray_stop + 1j * contour_height, phase_star))
    print(
        "route differences",
        mp.nstr(compact_difference, 12),
        mp.nstr(vertical_difference, 12),
        mp.nstr(ray_terminal_modulus, 12),
        flush=True,
    )
    require(compact_difference < mp.mpf("2e-25"), "compact coarse/fine disagreement")
    require(vertical_difference < mp.mpf("1e-25"), "vertical coarse/fine disagreement")
    require(ray_terminal_modulus < mp.mpf("1e-40"), "shifted ray did not damp sufficiently for route selection")

    artifact = {
        "kind": STEM,
        "status": "exploratory_high_precision_six_current_A_face_Morse_contour_route_selection_not_interval_certified",
        "passed": True,
        "scope": {
            "height": T,
            "A": A,
            "x_region": ["0", mp.nstr(x_cut, 50)],
            "q_region": [mp.nstr(q_cut, 50), "infinity"],
            "q_splice": Q_SPLICE_TEXT,
            "contour_height": CONTOUR_HEIGHT_TEXT,
            "retained_currents": TERMS,
        },
        "formula": {
            "oriented_core": "I_core=integral_0^x16 W_t(x)[-D_6(x)/x]exp(i*pi*A^2*x/4)dx",
            "logistic_amplitude": "-D_6(x(q))*x(q)^(-1/4)*(1-x(q))^(3/4)",
            "relative_phase": "Phi_A(q)-Phi_A(q_star)",
            "physical_projection": "2*(pi/(32t))^(1/4)*Re[exp(-i*pi/8)*I_core]",
        },
        "quadrature": {
            "compact_coarse": complex_record(compact_coarse),
            "compact_fine": complex_record(compact_fine),
            "compact_difference": mp.nstr(compact_difference, 30),
            "vertical_coarse": complex_record(vertical_coarse),
            "vertical_fine": complex_record(vertical_fine),
            "vertical_difference": mp.nstr(vertical_difference, 30),
            "shifted_ray_to_q_splice_plus_0_point_02": complex_record(shifted_ray),
            "shifted_ray_terminal_modulus": mp.nstr(ray_terminal_modulus, 30),
            "full_relative": complex_record(full_relative),
            "full_complex_with_carrier": complex_record(full_complex),
            "candidate_physical_core": mp.nstr(physical, 50),
            "six_current_replacement_physical_absolute_majorant": mp.nstr(replacement_absolute, 40),
        },
        "decision": {
            "floating_values_used_as_proof": False,
            "interval_certified": False,
            "coarse_fine_route_agreement": True,
            "six_current_core_worth_interval_certifying": True,
            "Morse_core_integral_bounded": False,
            "endpoint_layer_included": False,
            "full_signed_A_face_bound_proved": False,
            "R_after_A_bound_proved": False,
            "rh_implication": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "precision_decimal_digits": PRECISION,
            "process_priority": priority,
        },
        "dependency": {"path": str(DEPENDENCY.relative_to(REPO_ROOT)).replace("\\", "/"), "sha256": file_hash(DEPENDENCY)},
        "source": {"path": str(Path(__file__).resolve().relative_to(REPO_ROOT)).replace("\\", "/"), "sha256": file_hash(Path(__file__).resolve())},
        "proof_boundary": "High-precision floating route selection for the six-current tail-safe core only. No interval enclosure, endpoint-layer contribution, full A-face channel, R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"scouted six-current Morse core; candidate physical value {mp.nstr(physical, 16)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
