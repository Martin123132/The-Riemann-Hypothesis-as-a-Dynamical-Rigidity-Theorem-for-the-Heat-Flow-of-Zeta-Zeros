#!/usr/bin/env python3
"""Scout the exact A-face endpoint layer after the six-current handoff."""

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
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import erfcx


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_layer_pilot"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
DEPENDENCY = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_midpoint_Peano_partition_gate.json"

A = 159_577.0
T = 10_000_000_000.0
LEFT = 39_894.5
RIGHT = 39_936.5
MODES = np.arange(39_895.0, 39_937.0)
Q_CUT = 16.0


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


def phase_stripped_tail(y: np.ndarray | float, x: float) -> np.ndarray | complex:
    values = np.asarray(y, dtype=np.float64)
    delta = values - A * x / 2.0
    q = -delta * math.sqrt(2.0 / x)
    require(bool(np.all(q < 0.0)), "endpoint scout left the lower Fresnel sheet")
    w = np.exp(-1j * math.pi / 4.0) * math.sqrt(math.pi / 2.0) * (-q)
    h = np.exp(1j * math.pi / 4.0) * erfcx(w) / math.sqrt(2.0)
    result = -1.0 / (1j * math.pi) - values * math.sqrt(2.0 / x) * h
    if np.ndim(y) == 0:
        return complex(result)
    return result


def mapped_legendre(start: float, stop: float, nodes: int) -> tuple[np.ndarray, np.ndarray]:
    abscissae, weights = leggauss(nodes)
    midpoint = (start + stop) / 2.0
    half_width = (stop - start) / 2.0
    return midpoint + half_width * abscissae, half_width * weights


def panelled_legendre(start: float, stop: float, panels: int, nodes_per_panel: int) -> tuple[np.ndarray, np.ndarray]:
    abscissae, weights = leggauss(nodes_per_panel)
    panel_width = (stop - start) / panels
    node_parts = []
    weight_parts = []
    for panel in range(panels):
        panel_left = start + panel * panel_width
        midpoint = panel_left + panel_width / 2.0
        node_parts.append(midpoint + (panel_width / 2.0) * abscissae)
        weight_parts.append((panel_width / 2.0) * weights)
    return np.concatenate(node_parts), np.concatenate(weight_parts)


def make_defect(inner_nodes_per_cell: int):
    cell_abscissae, cell_weights = leggauss(inner_nodes_per_cell)
    strip_y_parts = []
    strip_weight_parts = []
    for cell in range(42):
        cell_left = LEFT + cell
        cell_midpoint = cell_left + 0.5
        strip_y_parts.append(cell_midpoint + 0.5 * cell_abscissae)
        strip_weight_parts.append(0.5 * cell_weights)
    strip_y = np.concatenate(strip_y_parts)
    strip_weights = np.concatenate(strip_weight_parts)

    def defect(x: float) -> complex:
        point_sum = complex(np.sum(phase_stripped_tail(MODES, x)))
        strip_integral = complex(np.dot(strip_weights, phase_stripped_tail(strip_y, x)))
        return point_sum - strip_integral

    return defect


def six_current_defect(x: float) -> complex:
    def g6(y: np.ndarray | float) -> np.ndarray | complex:
        values = np.asarray(y, dtype=np.float64)
        delta = values - A * x / 2.0
        result = A * x / (2j * math.pi * delta)
        double_factorial = 1
        for index in range(1, 6):
            double_factorial *= 2 * index - 1
            result += values * double_factorial * (x / 2.0) ** index / ((1j * math.pi) ** (index + 1) * delta ** (2 * index + 1))
        if np.ndim(y) == 0:
            return complex(result)
        return result

    a = A * x / 2.0
    dl = LEFT - a
    dr = RIGHT - a
    strip = A * x * math.log(dr / dl) / (2j * math.pi)
    double_factorial = 1
    for index in range(1, 6):
        double_factorial *= 2 * index - 1
        coefficient = double_factorial * (x / 2.0) ** index / (1j * math.pi) ** (index + 1)
        strip += coefficient * (
            (dr ** (1 - 2 * index) - dl ** (1 - 2 * index)) / (1 - 2 * index)
            - a * (dr ** (-2 * index) - dl ** (-2 * index)) / (2 * index)
        )
    return complex(np.sum(g6(MODES))) - strip


def stable_phase_relative(q: float, q_star: float, x_star: float) -> float:
    u = q - q_star
    exponential_minus_one = math.expm1(u)
    term = u * u / 2.0
    exponential_minus_one_minus_u = term
    for order in range(3, 14):
        term *= u / order
        exponential_minus_one_minus_u += term
    one_minus_x_star = 1.0 - x_star
    denominator = 1.0 + one_minus_x_star * exponential_minus_one
    numerator = -exponential_minus_one_minus_u + one_minus_x_star * u * exponential_minus_one
    return (math.pi * A * A / 4.0) * x_star * one_minus_x_star * numerator / denominator


def endpoint_integral(inner_nodes_per_cell: int, outer_panels: int, outer_nodes_per_panel: int, q_cut: float, x_star: float, q_star: float) -> tuple[complex, complex, complex]:
    defect = make_defect(inner_nodes_per_cell)
    q_nodes, q_weights = panelled_legendre(0.0, q_cut, outer_panels, outer_nodes_per_panel)
    values = []
    defects = []
    for q in q_nodes:
        x = 1.0 / (1.0 + math.exp(q))
        d_value = defect(float(x))
        phase_relative = stable_phase_relative(float(q), q_star, x_star)
        amplitude = -d_value * x ** (-0.25) * (1.0 - x) ** 0.75
        values.append(amplitude * np.exp(1j * phase_relative))
        defects.append(d_value)
    integral = complex(np.dot(q_weights, np.asarray(values, dtype=np.complex128)))
    return integral, defects[0], defects[-1]


def complex_record(value: complex) -> dict[str, float]:
    return {"real": value.real, "imag": value.imag, "modulus": abs(value)}


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependency = json.loads(DEPENDENCY.read_text(encoding="utf-8"))
    require(dependency.get("passed") is True, "endpoint partition dependency is not passed")
    require(dependency["decision"]["endpoint_layer_phase_derivative_above_21000_proved"] is True, "endpoint nonstationarity drift")

    sqrt_x_cut = (-Q_CUT * math.sqrt(2.0) + math.sqrt(2.0 * Q_CUT**2 + 8.0 * A * LEFT)) / (2.0 * A)
    x_cut = sqrt_x_cut**2
    x_star = (1.0 - math.sqrt(1.0 - 8.0 * T / (math.pi * A * A))) / 2.0
    q_cut = math.log((1.0 - x_cut) / x_cut)
    q_star = math.log((1.0 - x_star) / x_star)

    coarse, coarse_first, coarse_last = endpoint_integral(8, 32, 16, q_cut, x_star, q_star)
    inner_refined, _, _ = endpoint_integral(16, 32, 16, q_cut, x_star, q_star)
    fine, fine_first, fine_last = endpoint_integral(16, 64, 24, q_cut, x_star, q_star)
    difference = abs(fine - coarse)
    inner_difference = abs(inner_refined - coarse)
    outer_difference = abs(fine - inner_refined)

    exact_defect = make_defect(24)
    handoff_exact = exact_defect(x_cut)
    handoff_six = six_current_defect(x_cut)
    handoff_difference = abs(handoff_exact - handoff_six)
    endpoint_defect = exact_defect(0.5)

    mp.mp.dps = 80
    aa = mp.mpf(159_577)
    tt = mp.mpf(10_000_000_000)
    xs = (1 - mp.sqrt(1 - 8 * tt / (mp.pi * aa**2))) / 2
    phi_star = mp.pi * aa**2 * xs / 4 + tt * mp.log((1 - xs) / xs) / 2
    full_complex = mp.exp(1j * phi_star) * mp.mpc(fine.real, fine.imag)
    normalization = 2 * (mp.pi / (32 * tt)) ** (mp.mpf(1) / 4)
    physical = normalization * mp.re(mp.exp(-1j * mp.pi / 8) * full_complex)

    print(
        "endpoint route differences",
        f"{difference:.12e}",
        f"{inner_difference:.12e}",
        f"{outer_difference:.12e}",
        f"{handoff_difference:.12e}",
        flush=True,
    )
    print("endpoint route values", repr(coarse), repr(inner_refined), repr(fine), flush=True)
    require(difference < 2.0e-10, "endpoint coarse/fine disagreement")
    require(handoff_difference < 2.0e-8, "six-current/exact handoff mismatch")

    artifact = {
        "kind": STEM,
        "status": "exploratory_double_precision_exact_Fresnel_endpoint_layer_route_selection_not_interval_certified",
        "passed": True,
        "scope": {
            "height": int(T),
            "A": int(A),
            "x_region": [x_cut, 0.5],
            "phase_derivative_floor_from_exact_dependency": ">21000",
            "coarse_nodes": {"inner_per_unit_cell": 8, "inner_total": 336, "outer_panels": 32, "outer_per_panel": 16, "outer_total": 512},
            "fine_nodes": {"inner_per_unit_cell": 16, "inner_total": 672, "outer_panels": 64, "outer_per_panel": 24, "outer_total": 1536},
        },
        "formula": {
            "oriented_endpoint_layer": "I_end=integral_x16^(1/2) W_t(x)[-D_42(x)/x]exp(i*pi*A^2*x/4)dx",
            "D_42": "sum_(m=39895)^39936 G_A(m,x)-integral_(39894.5)^39936.5 G_A(y,x)dy",
            "G_A": "-1/(i*pi)-y*sqrt(2/x)*exp(i*pi/4)erfcx[e^(-i*pi/4)sqrt(pi/2)(-q_A)]/sqrt(2)",
        },
        "quadrature": {
            "coarse_relative_integral": complex_record(coarse),
            "fine_relative_integral": complex_record(fine),
            "coarse_fine_difference": difference,
            "inner_rule_difference_on_coarse_outer_panels": inner_difference,
            "outer_panel_rule_difference_at_16_inner_nodes_per_cell": outer_difference,
            "candidate_physical_endpoint_layer": mp.nstr(physical, 50),
            "handoff_exact_defect": complex_record(handoff_exact),
            "handoff_six_current_defect": complex_record(handoff_six),
            "handoff_difference": handoff_difference,
            "endpoint_defect": complex_record(endpoint_defect),
            "first_outer_node_defect_coarse": complex_record(coarse_first),
            "last_outer_node_defect_coarse": complex_record(coarse_last),
            "first_outer_node_defect_fine": complex_record(fine_first),
            "last_outer_node_defect_fine": complex_record(fine_last),
        },
        "decision": {
            "floating_values_used_as_proof": False,
            "interval_certified": False,
            "coarse_fine_route_agreement": True,
            "endpoint_layer_worth_interval_certifying": True,
            "endpoint_layer_IBP_remainder_bounded": False,
            "full_signed_A_face_bound_proved": False,
            "R_after_A_bound_proved": False,
            "rh_implication": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "blas_threads": 1,
            "process_priority": priority,
        },
        "dependency": {"path": str(DEPENDENCY.relative_to(REPO_ROOT)).replace("\\", "/"), "sha256": file_hash(DEPENDENCY)},
        "source": {"path": str(Path(__file__).resolve().relative_to(REPO_ROOT)).replace("\\", "/"), "sha256": file_hash(Path(__file__).resolve())},
        "proof_boundary": "Double-precision exact-Fresnel route selection for the short endpoint layer only. No interval enclosure, integration-by-parts remainder, full A-face channel, R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"scouted exact A-face endpoint layer; candidate physical value {mp.nstr(physical, 16)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
