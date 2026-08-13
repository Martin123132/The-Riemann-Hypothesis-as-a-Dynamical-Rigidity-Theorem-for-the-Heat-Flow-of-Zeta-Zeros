#!/usr/bin/env python3
"""Build the fixed-cell low-N determinant evaluator contract gate."""

from __future__ import annotations

import ctypes
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Callable

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad_vec


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_fixed_cell_low_n_determinant_evaluator_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "determinant_chart": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_actual_source_finite_determinant_reciprocal_block_symmetry_gate.json",
    "finite_band": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_finite_band_"
        "dirichlet_cell_reduction_gate.json"
    ),
    "outer_pairing": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_ideal_cubic_"
        "symmetric_outer_pairing_gate.json"
    ),
    "calibrated_roots": REPO_ROOT
    / "work/rh_compute/results"
    / "jensen_window_pf_newman_c1_calibrated_roster_interior_phase_root_scout.json",
}

MAX_DEGREE = 6
GAUSS_ORDER = 512
QUAD_EPS = 2.0e-11
ALPHA_STEP = 1.0e-5
Z0 = complex(3.0 / 5.0, 4.0 / 5.0)
EDGE = np.array([2.0 / 7.0, -3.0 / 8.0, 5.0 / 11.0, -7.0 / 13.0])

FIXTURES = (
    {"id": "n4", "N": 4, "a": 4.25, "alpha": 17.8, "t": 0.08, "sigma": 0.55},
    {"id": "n5", "N": 5, "a": 5.4, "alpha": 28.9, "t": 0.05, "sigma": 0.52},
    {"id": "n7", "N": 7, "a": 7.35, "alpha": 53.6, "t": 0.03, "sigma": 0.51},
)

# Rows are V, mathcal N, A, Q. Columns are powers lambda^0,...,lambda^5.
COEFFICIENT_ROWS = np.array(
    [
        [1 / 3 + 1j / 5, -2 / 7 + 1j / 11, 1 / 13 - 1j / 17, -1 / 19, 1j / 23, 1 / 29],
        [-2 / 5 + 1j / 7, 3 / 11 - 1j / 13, -1 / 17 + 2j / 19, 1 / 23, -1j / 29, -1 / 31],
        [3 / 8 - 1j / 9, -1 / 10 + 2j / 15, 1 / 21 + 1j / 22, -1 / 26, 1 / 33, 1j / 35],
        [-1 / 6 + 2j / 9, 1 / 14 - 1j / 16, -2 / 25 + 1j / 27, 1j / 30, -1 / 34, 1 / 39],
    ],
    dtype=np.complex128,
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
    proof_boundary: str


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def set_below_normal_priority() -> None:
    if os.name != "nt":
        return
    try:
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        kernel32.SetPriorityClass.argtypes = (ctypes.c_void_p, ctypes.c_uint32)
        kernel32.SetPriorityClass.restype = ctypes.c_int
        kernel32.SetPriorityClass(kernel32.GetCurrentProcess(), 0x00004000)
    except (AttributeError, OSError):
        pass


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(content, encoding="utf-8")
    os.replace(temporary, path)


def load_and_audit_sources() -> dict[str, dict[str, str]]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.is_file(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))

    require(
        payloads["determinant_chart"]["exact"]["reciprocal_block_classes"] == 3,
        "determinant partition drifted",
    )
    finite = payloads["finite_band"]["symbolic_certificate"]
    require(
        "mathscr F_N[P]=sum_(q=1)^N sum_(p=1)^N"
        in finite["reciprocal_blocks"]["blocks"],
        "finite block identity drifted",
    )
    pairing = payloads["outer_pairing"]["symbolic_certificate"]
    require(
        "mathscr C_N[P]=sum_(k=n_N+1)^infinity"
        in pairing["cutoff"]["remote"],
        "outer paired-remote definition drifted",
    )
    require(
        "mathscr F_N[P]+mathscr O_N[P]=mathscr M_N[P]"
        in pairing["recombination"]["poisson"],
        "outer recombination identity drifted",
    )
    require(
        payloads["calibrated_roots"]["summary"]["calibrated_roots"] == 3,
        "calibrated-root count drifted",
    )
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def roster(n_value: int, alpha: float) -> tuple[int, int, np.ndarray]:
    lower = math.floor(alpha / (n_value + 0.5)) + 1
    upper = math.floor(2.0 * alpha)
    require(lower <= upper, "empty fixture roster")
    return lower, upper, np.arange(lower, upper + 1, dtype=np.int64)


def reciprocal_labels(alpha: float, modes: np.ndarray, n_value: int) -> np.ndarray:
    labels = np.floor(alpha / modes + 0.5).astype(np.int64)
    require(bool(np.all((1 <= labels) & (labels <= n_value))), "reciprocal label escaped")
    return labels


def tie_margin(n_value: int, alpha: float, modes: np.ndarray) -> float:
    upper = abs(2.0 * alpha - round(2.0 * alpha)) / 2.0
    lower_coordinate = alpha / (n_value + 0.5)
    lower = abs(lower_coordinate - round(lower_coordinate)) * (n_value + 0.5)
    internal = []
    for mode in modes:
        coordinate = alpha / float(mode)
        nearest_half = round(coordinate - 0.5) + 0.5
        internal.append(abs(alpha - float(mode) * nearest_half))
    return min([upper, lower, *internal])


def cell_bounds(q_value: int, n_value: int) -> tuple[float, float]:
    return max(1.0, q_value - 0.5), min(float(n_value), q_value + 0.5)


def powers(lam: np.ndarray | float) -> np.ndarray:
    values = np.asarray(lam)
    return np.stack([values**degree for degree in range(MAX_DEGREE + 1)], axis=-1)


def source_basis(
    u: np.ndarray | float, alpha: float, t_value: float, sigma: float
) -> np.ndarray:
    u_array = np.asarray(u)
    lam = np.log(u_array)
    exponent = t_value * lam**2 / 4.0 - sigma * lam + 2j * np.pi * alpha * lam
    return np.exp(exponent)[..., None] * powers(lam)


def source_prime_basis(
    u: float, alpha: float, t_value: float, sigma: float
) -> np.ndarray:
    lam = math.log(u)
    poly = powers(lam)
    first = np.zeros(MAX_DEGREE + 1, dtype=np.complex128)
    for degree in range(1, MAX_DEGREE + 1):
        first[degree] = degree * lam ** (degree - 1)
    g_value = t_value * lam / 2.0 - sigma + 2j * np.pi * alpha
    scalar = np.exp(t_value * lam**2 / 4.0 - sigma * lam + 2j * np.pi * alpha * lam)
    return scalar * (first + g_value * poly) / u


def source_second_basis(
    u: float, alpha: float, t_value: float, sigma: float
) -> np.ndarray:
    lam = math.log(u)
    poly = powers(lam)
    first = np.zeros(MAX_DEGREE + 1, dtype=np.complex128)
    second = np.zeros(MAX_DEGREE + 1, dtype=np.complex128)
    for degree in range(1, MAX_DEGREE + 1):
        first[degree] = degree * lam ** (degree - 1)
    for degree in range(2, MAX_DEGREE + 1):
        second[degree] = degree * (degree - 1) * lam ** (degree - 2)
    g_value = t_value * lam / 2.0 - sigma + 2j * np.pi * alpha
    scalar = np.exp(t_value * lam**2 / 4.0 - sigma * lam + 2j * np.pi * alpha * lam)
    bracket = second + (2.0 * g_value - 1.0) * first
    bracket += (g_value**2 - g_value + t_value / 2.0) * poly
    return scalar * bracket / u**2


def integrate_cells(
    function: Callable[[float], np.ndarray], n_value: int, eps: float = QUAD_EPS
) -> np.ndarray:
    total = np.zeros(MAX_DEGREE + 1, dtype=np.complex128)
    for q_value in range(1, n_value + 1):
        left, right = cell_bounds(q_value, n_value)
        value, _ = quad_vec(function, left, right, epsabs=eps, epsrel=eps, limit=500)
        total += value
    return total


def carrier_basis(
    n_value: int, alpha: float, t_value: float, sigma: float
) -> np.ndarray:
    total = np.zeros(MAX_DEGREE + 1, dtype=np.complex128)
    for q_value in range(1, n_value + 1):
        weight = 0.5 if q_value in (1, n_value) else 1.0
        total += weight * source_basis(float(q_value), alpha, t_value, sigma)
    return total


def block_basis(
    n_value: int,
    alpha: float,
    t_value: float,
    sigma: float,
    modes: np.ndarray,
    labels: np.ndarray,
    order: int = GAUSS_ORDER,
) -> np.ndarray:
    nodes, weights = leggauss(order)
    blocks = np.zeros(
        (n_value, n_value, MAX_DEGREE + 1), dtype=np.complex128
    )
    modes_float = modes.astype(np.float64)
    for q_value in range(1, n_value + 1):
        left, right = cell_bounds(q_value, n_value)
        u = (right - left) * nodes / 2.0 + (right + left) / 2.0
        scaled_weights = (right - left) * weights / 2.0
        amplitude = source_basis(u, alpha, t_value, sigma)
        oscillation = np.exp(-2j * np.pi * modes_float[:, None] * u[None, :])
        atoms = oscillation @ (scaled_weights[:, None] * amplitude)
        for mode_index, p_value in enumerate(labels):
            blocks[q_value - 1, p_value - 1] += atoms[mode_index]
    return blocks


def direct_band_basis(
    n_value: int,
    alpha: float,
    t_value: float,
    sigma: float,
    modes: np.ndarray,
) -> np.ndarray:
    modes_float = modes.astype(np.float64)

    def integrand(u: float) -> np.ndarray:
        kernel = np.exp(-2j * np.pi * modes_float * u).sum()
        return source_basis(u, alpha, t_value, sigma) * kernel

    return integrate_cells(integrand, n_value)


def variation_band_basis(
    n_value: int,
    alpha: float,
    t_value: float,
    sigma: float,
    modes: np.ndarray,
) -> tuple[np.ndarray, complex, np.ndarray]:
    modes_float = modes.astype(np.float64)
    odd = modes[modes % 2 == 1].astype(np.float64)
    harmonic = (1.0 / (np.pi * 1j)) * np.sum(1.0 / odd) if len(odd) else 0j

    def kernel(v: float) -> complex:
        return np.exp(-2j * np.pi * modes_float * v).sum()

    lower_value = source_basis(1.0, alpha, t_value, sigma)
    upper_value = source_basis(float(n_value), alpha, t_value, sigma)
    variation = np.zeros(MAX_DEGREE + 1, dtype=np.complex128)

    value, _ = quad_vec(
        lambda v: (source_basis(1.0 + v, alpha, t_value, sigma) - lower_value)
        * kernel(v),
        0.0,
        0.5,
        epsabs=QUAD_EPS,
        epsrel=QUAD_EPS,
        limit=500,
    )
    variation += value
    for q_value in range(2, n_value):
        centre = source_basis(float(q_value), alpha, t_value, sigma)
        value, _ = quad_vec(
            lambda v, q=q_value, c=centre: (
                source_basis(q + v, alpha, t_value, sigma) - c
            )
            * kernel(v),
            -0.5,
            0.5,
            epsabs=QUAD_EPS,
            epsrel=QUAD_EPS,
            limit=500,
        )
        variation += value
    value, _ = quad_vec(
        lambda v: (
            source_basis(float(n_value) + v, alpha, t_value, sigma) - upper_value
        )
        * kernel(v),
        -0.5,
        0.0,
        epsabs=QUAD_EPS,
        epsrel=QUAD_EPS,
        limit=500,
    )
    variation += value
    return harmonic * (lower_value - upper_value) + variation, harmonic, variation


def near_basis(
    n_value: int,
    alpha: float,
    t_value: float,
    sigma: float,
    lower: int,
    upper: int,
) -> np.ndarray:
    frequencies = np.arange(-(lower - 1), upper + 1, dtype=np.float64)

    def integrand(u: float) -> np.ndarray:
        kernel = np.exp(2j * np.pi * frequencies * u).sum()
        return source_basis(u, alpha, t_value, sigma) * kernel

    return integrate_cells(integrand, n_value)


def remote_basis(
    n_value: int,
    alpha: float,
    t_value: float,
    sigma: float,
    upper: int,
) -> np.ndarray:
    positive = np.arange(1, upper + 1, dtype=np.float64)
    inverse_squares = 1.0 / positive**2
    zeta_tail = np.pi**2 / 6.0 - inverse_squares.sum()
    endpoint = source_prime_basis(float(n_value), alpha, t_value, sigma)
    endpoint -= source_prime_basis(1.0, alpha, t_value, sigma)

    def tail_kernel(u: float) -> float:
        fractional = u - math.floor(u)
        bernoulli = fractional**2 - fractional + 1.0 / 6.0
        finite = np.sum(np.cos(2.0 * np.pi * positive * u) * inverse_squares)
        return np.pi**2 * bernoulli - finite

    integral = integrate_cells(
        lambda u: source_second_basis(u, alpha, t_value, sigma) * tail_kernel(u),
        n_value,
    )
    return (endpoint * zeta_tail - integral) / (2.0 * np.pi**2)


def split_blocks(blocks: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    diagonal = np.zeros(MAX_DEGREE + 1, dtype=np.complex128)
    adjacent = np.zeros(MAX_DEGREE + 1, dtype=np.complex128)
    far = np.zeros(MAX_DEGREE + 1, dtype=np.complex128)
    n_value = blocks.shape[0]
    for q_index in range(n_value):
        for p_index in range(n_value):
            distance = abs(q_index - p_index)
            if distance == 0:
                diagonal += blocks[q_index, p_index]
            elif distance == 1:
                adjacent += blocks[q_index, p_index]
            else:
                far += blocks[q_index, p_index]
    return diagonal, adjacent, far


def evaluate_functionals(config: dict, order: int = GAUSS_ORDER) -> dict:
    n_value = int(config["N"])
    alpha = float(config["alpha"])
    t_value = float(config["t"])
    sigma = float(config["sigma"])
    lower, upper, modes = roster(n_value, alpha)
    labels = reciprocal_labels(alpha, modes, n_value)
    blocks = block_basis(n_value, alpha, t_value, sigma, modes, labels, order)
    diagonal, adjacent, far = split_blocks(blocks)
    carrier = carrier_basis(n_value, alpha, t_value, sigma)
    direct = direct_band_basis(n_value, alpha, t_value, sigma, modes)
    variation, harmonic, variation_only = variation_band_basis(
        n_value, alpha, t_value, sigma, modes
    )
    near = near_basis(n_value, alpha, t_value, sigma, lower, upper)
    remote = remote_basis(n_value, alpha, t_value, sigma, upper)
    return {
        "lower": lower,
        "upper": upper,
        "modes": modes,
        "labels": labels,
        "tie_margin": tie_margin(n_value, alpha, modes),
        "carrier": carrier,
        "diagonal": diagonal,
        "adjacent": adjacent,
        "far": far,
        "band_blocks": diagonal + adjacent + far,
        "band_direct": direct,
        "band_variation": variation,
        "harmonic": harmonic,
        "variation": variation_only,
        "near": near,
        "remote": remote,
    }


def lift_coefficients(coefficients: np.ndarray, log_a: float) -> np.ndarray:
    lifted = np.zeros(MAX_DEGREE + 1, dtype=np.complex128)
    for degree, coefficient in enumerate(coefficients):
        lifted[degree] += -1j * log_a * coefficient
        lifted[degree + 1] += 1j * coefficient
    return lifted


def project(functional: np.ndarray, coefficients: np.ndarray, z_value: complex) -> float:
    return float(np.real(z_value * np.dot(coefficients, functional[: len(coefficients)])))


def observation(functional: np.ndarray, z_value: complex, log_a: float, *, tangent: bool) -> np.ndarray:
    result = np.zeros(4, dtype=np.float64)
    for row_index, coefficients in enumerate(COEFFICIENT_ROWS):
        if tangent:
            lifted = lift_coefficients(coefficients, log_a)
            result[row_index] = project(functional, lifted, z_value)
        else:
            result[row_index] = project(functional, coefficients, z_value)
    return result


def quadratic(vector: np.ndarray) -> float:
    return float(vector[0] * vector[1] - vector[2] * vector[3])


def bilinear(left: np.ndarray, right: np.ndarray) -> float:
    return float(
        (
            left[0] * right[1]
            + left[1] * right[0]
            - left[2] * right[3]
            - left[3] * right[2]
        )
        / 2.0
    )


def primitive_increment(base: np.ndarray, defect: np.ndarray) -> float:
    return quadratic(base - defect) - quadratic(base)


def current_increment(
    base: np.ndarray,
    defect: np.ndarray,
    base_dot: np.ndarray,
    defect_dot: np.ndarray,
) -> float:
    return float(
        -2.0 * bilinear(base, defect_dot)
        - 2.0 * bilinear(base_dot, defect)
        + 2.0 * bilinear(defect, defect_dot)
    )


def projected_state(functionals: dict, config: dict, z_value: complex) -> dict:
    log_a = math.log(float(config["a"]))
    names = ("carrier", "diagonal", "adjacent", "far")
    state = {}
    for name in names:
        state[name] = observation(functionals[name], z_value, log_a, tangent=False)
        state[f"{name}_dot"] = observation(
            functionals[name], z_value, log_a, tangent=True
        )
    p_value = state["carrier"]
    p_dot = state["carrier_dot"]
    d0 = state["diagonal"] - p_value
    d1 = state["adjacent"]
    d2 = state["far"]
    dd0 = state["diagonal_dot"] - p_dot
    dd1 = state["adjacent_dot"]
    dd2 = state["far_dot"]
    base = EDGE + p_value
    base_dot = p_dot
    defect = d0 + d1 + d2
    defect_dot = dd0 + dd1 + dd2
    primitive = primitive_increment(base, defect)
    current = current_increment(base, defect, base_dot, defect_dot)
    primitive_parts = (
        primitive_increment(base, d0),
        primitive_increment(base - d0, d1),
        primitive_increment(base - d0 - d1, d2),
    )
    current_parts = (
        current_increment(base, d0, base_dot, dd0),
        current_increment(base - d0, d1, base_dot - dd0, dd1),
        current_increment(base - d0 - d1, d2, base_dot - dd0 - dd1, dd2),
    )
    return {
        **state,
        "base": base,
        "base_dot": base_dot,
        "d0": d0,
        "d1": d1,
        "d2": d2,
        "dd0": dd0,
        "dd1": dd1,
        "dd2": dd2,
        "defect": defect,
        "defect_dot": defect_dot,
        "primitive": primitive,
        "current": current,
        "primitive_parts": primitive_parts,
        "current_parts": current_parts,
    }


def finite_difference_certificate(config: dict, centre: dict, centre_state: dict) -> dict:
    alpha = float(config["alpha"])
    log_a = math.log(float(config["a"]))
    delta = ALPHA_STEP
    states = []
    labels = []
    for direction in (-1.0, 1.0):
        shifted = dict(config)
        shifted["alpha"] = alpha + direction * delta
        functionals = evaluate_functionals(shifted)
        z_value = Z0 * np.exp(-2j * np.pi * direction * delta * log_a)
        states.append(projected_state(functionals, shifted, z_value))
        labels.append(tuple(int(value) for value in functionals["labels"]))
    require(labels[0] == labels[1], "finite difference crossed an internal tie")
    require(
        labels[0] == tuple(int(value) for value in centre["labels"]),
        "finite difference changed the reciprocal assignment",
    )
    xi_denominator = 4.0 * np.pi * delta
    primitive_fd = (states[1]["primitive"] - states[0]["primitive"]) / xi_denominator
    sector_errors = []
    for name in ("carrier", "diagonal", "adjacent", "far"):
        numerical = (states[1][name] - states[0][name]) / xi_denominator
        sector_errors.append(float(np.max(np.abs(numerical - centre_state[f"{name}_dot"]))))
    return {
        "alpha_step": delta,
        "primitive_derivative": primitive_fd,
        "analytic_current": centre_state["current"],
        "primitive_current_error": abs(primitive_fd - centre_state["current"]),
        "maximum_sector_tangent_error": max(sector_errors),
        "assignments_unchanged": True,
    }


def ideal_coefficients(config: dict) -> tuple[np.ndarray, np.ndarray]:
    log_a = math.log(float(config["a"]))
    u_n = log_a - math.log(float(config["N"]))
    u_n_x = 1.0 / (8.0 * np.pi * float(config["a"]) ** 2)
    x = np.array([-log_a, 1.0], dtype=np.complex128)
    plus = x.copy()
    plus[0] += u_n
    minus = x.copy()
    minus[0] -= u_n
    hermitian = -0.25j * np.polynomial.polynomial.polymul(
        plus, np.polynomial.polynomial.polymul(minus, minus)
    )
    transpose = -0.25j * np.polynomial.polynomial.polymul(
        minus, np.polynomial.polynomial.polymul(plus, plus)
    )
    transpose[: len(minus)] += u_n_x * minus
    return hermitian, transpose


def complex_payload(values: np.ndarray | complex) -> dict | list[dict]:
    array = np.asarray(values)
    if array.ndim == 0:
        scalar = complex(array)
        return {"re": f"{scalar.real:.17e}", "im": f"{scalar.imag:.17e}"}
    return [complex_payload(value) for value in array]


def real_payload(values: np.ndarray) -> list[str]:
    return [f"{float(value):.17e}" for value in values]


def fixture_certificate(config: dict) -> dict:
    functionals = evaluate_functionals(config)
    state = projected_state(functionals, config, Z0)
    finite_difference = finite_difference_certificate(config, functionals, state)
    band = functionals["band_blocks"]
    outer = functionals["carrier"] - band
    remote_from_partition = outer - functionals["near"]
    hermitian, transpose = ideal_coefficients(config)

    errors = {
        "block_vs_direct": float(np.max(np.abs(band - functionals["band_direct"]))),
        "block_vs_variation": float(np.max(np.abs(band - functionals["band_variation"]))),
        "band_partition": float(
            np.max(
                np.abs(
                    band
                    - functionals["diagonal"]
                    - functionals["adjacent"]
                    - functionals["far"]
                )
            )
        ),
        "outer_near_remote": float(
            np.max(np.abs(remote_from_partition - functionals["remote"]))
        ),
        "primitive_telescope": abs(state["primitive"] - sum(state["primitive_parts"])),
        "current_telescope": abs(state["current"] - sum(state["current_parts"])),
        "primitive_current_finite_difference": finite_difference["primitive_current_error"],
        "sector_tangent_finite_difference": finite_difference[
            "maximum_sector_tangent_error"
        ],
    }
    require(functionals["tie_margin"] > 1000.0 * ALPHA_STEP, "fixture too near a tie")
    require(errors["block_vs_direct"] < 2.0e-9, "direct band mismatch")
    require(errors["block_vs_variation"] < 2.0e-9, "variation band mismatch")
    require(errors["outer_near_remote"] < 3.0e-9, "near/remote mismatch")
    require(errors["primitive_telescope"] < 2.0e-10, "primitive telescope mismatch")
    require(errors["current_telescope"] < 2.0e-9, "current telescope mismatch")
    require(
        errors["primitive_current_finite_difference"] < 2.0e-5,
        "primitive current finite difference mismatch",
    )
    require(
        errors["sector_tangent_finite_difference"] < 2.0e-5,
        "sector tangent finite difference mismatch",
    )

    basis_names = (
        "carrier",
        "diagonal",
        "adjacent",
        "far",
        "band_blocks",
        "near",
        "remote",
    )
    ideal = {}
    for name, coefficients in (("P_H0", hermitian), ("P_T0", transpose)):
        ideal[name] = {
            key: complex_payload(np.dot(coefficients, functionals[key][: len(coefficients)]))
            for key in basis_names
        }
        ideal[name]["coefficients"] = complex_payload(coefficients)

    return {
        "id": config["id"],
        "parameters": {key: config[key] for key in ("N", "a", "alpha", "t", "sigma")},
        "roster": {
            "m_N": functionals["lower"],
            "n_N": functionals["upper"],
            "mode_count": len(functionals["modes"]),
            "minimum_alpha_tie_margin": functionals["tie_margin"],
            "labels": [int(value) for value in functionals["labels"]],
        },
        "basis": {key: complex_payload(functionals[key]) for key in basis_names},
        "endpoint_harmonic": complex_payload(functionals["harmonic"]),
        "ideal_relative_channels": ideal,
        "determinant": {
            "edge": real_payload(EDGE),
            "base": real_payload(state["base"]),
            "defect": real_payload(state["defect"]),
            "primitive": f"{state['primitive']:.17e}",
            "primitive_parts": [f"{value:.17e}" for value in state["primitive_parts"]],
            "current": f"{state['current']:.17e}",
            "current_parts": [f"{value:.17e}" for value in state["current_parts"]],
        },
        "finite_difference": finite_difference,
        "errors": errors,
    }


def exact_contract() -> dict:
    return {
        "basis": "For P(lambda)=sum_(k=0)^5 c_k lambda^k, tabulate mathscr M_N, the p=q, |p-q|=1, and |p-q|>=2 pieces of mathscr F_N, mathscr N_N, and mathscr C_N on lambda^0,...,lambda^6. Degree six is retained because the fixed-cell tangent lift i(lambda-log a)P raises degree by one.",
        "blocks": "J_q=[max(1,q-1/2),min(N,q+1/2)] and C_p={r>0:p-1/2<=alpha_P/r<p+1/2}. Every roster mode is assigned once. The block sum is checked against both the direct finite Dirichlet kernel and H_N{f_P(1)-f_P(N)}+V_N[P].",
        "remote": "The common-cutoff paired remote functional is evaluated without truncating an infinite series: sum_(k>n)1/k^2 and sum_(k>n)cos(2pi k u)/k^2 are replaced by zeta(2) minus finite sums and the periodic Bernoulli B_2 kernel after two exact integrations by parts.",
        "lift": "With S and P fixed and z=nu c_xi carrying exp(-i xi log a), partial_xi Re{z L[P]}=Re{z L[i(lambda-log a)P]} on a fixed roster-and-reciprocal-cell chart.",
        "determinant": "For b=omega_E+p, d_0=f^(0)-p, d_1=f^(1), d_2=f^(2), the evaluator computes Q(b-d)-Q(b) and its diagonal-first primitive/current telescope. A centered xi finite difference independently checks the lifted current while every mode label remains fixed.",
        "fixture_scope": "The fixtures retain the physical scale relation a^2-1<alpha_P<=a^2, the exact ideal cubics, the completed Fourier normalization, and all endpoint halves, but use small N and synthetic four-row polynomial coefficients. They validate the evaluator contract; they are not low-height Xi data and carry no physical sign conclusion.",
        "next_action": "Extract the four actual physical coefficient rows P_X and the genuine edge vector at each calibrated N about 7.2e10. Replace only the direct block engine by a tie-complete compressed diagonal/adjacent plus rail-compressed far evaluator, while preserving this basis, lift, near/remote, and determinant API as regression oracles.",
        "pi_provenance": "Every pi is forced by e(x)=exp(2pi i x), alpha_P=xi/(2pi), kappa=1/(2pi i), and the exact Fourier/Bernoulli identities. No geometric or fitted pi is introduced.",
        "proof_boundary": "This gate proves an executable low-N contract and cross-representation numerical identities to recorded tolerances. It is floating-point quadrature, not interval arithmetic. It proves no physical-N block value or sign, no retained determinant bound, no complete-current inequality, no all-q transport, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }


def build_rows(contract: dict, fixtures: list[dict]) -> list[GateRow]:
    maximum = lambda key: max(row["errors"][key] for row in fixtures)
    rows = [
        GateRow("lnd_01_domain", "fixed-cell domain", "validated", "All fixtures remain inside one reciprocal chart under the xi finite difference.", f"minimum alpha tie margin={min(row['roster']['minimum_alpha_tie_margin'] for row in fixtures):.6g}; step={ALPHA_STEP}", "No tie derivative is taken."),
        GateRow("lnd_02_basis", "functional basis", "validated", "The evaluator stores a reusable degree-six monomial basis.", contract["basis"], "Actual physical coefficients are not inserted yet."),
        GateRow("lnd_03_carrier", "starred carrier", "validated", "Both endpoint half weights are retained in mathscr M_N.", "The q=1 and q=N atoms have weight 1/2 in every fixture.", "This is finite Fourier endpoint ownership."),
        GateRow("lnd_04_blocks", "reciprocal blocks", "validated", "Every band mode belongs to exactly one C_p and every physical cell to one J_q.", contract["blocks"], "Boundary measure-zero choices do not alter integrals."),
        GateRow("lnd_05_direct", "direct kernel audit", "validated", "The block matrix equals the direct finite Dirichlet-kernel integral.", f"maximum discrepancy={maximum('block_vs_direct'):.3e}", "This is numerical quadrature."),
        GateRow("lnd_06_variation", "cell variation audit", "validated", "The block matrix equals the endpoint-harmonic plus variation formula.", f"maximum discrepancy={maximum('block_vs_variation'):.3e}", "Cancellation is retained before moduli."),
        GateRow("lnd_07_partition", "three sectors", "validated", "Diagonal, adjacent, and far blocks partition the finite band.", f"maximum discrepancy={maximum('band_partition'):.3e}", "No sector sign is inferred."),
        GateRow("lnd_08_near", "finite near", "validated", "The asymmetric mean-one near roster is evaluated as one contiguous kernel.", "frequencies -(m_N-1),...,n_N are summed before integration", "The kernel is complex, not positive."),
        GateRow("lnd_09_remote", "paired remote", "validated", "The common-cutoff remote tail is finite-evaluated through the Bernoulli kernel.", contract["remote"], "The one-sided conditional tails are never separated."),
        GateRow("lnd_10_outer", "outer recomposition", "validated", "mathscr M_N-mathscr F_N equals near plus paired remote.", f"maximum discrepancy={maximum('outer_near_remote'):.3e}", "This is a low-N numerical audit of the exact identity."),
        GateRow("lnd_11_lift", "tangent lift", "validated", "The degree-raising lift supplies every observation tangent.", contract["lift"], "S and P are held fixed on the auxiliary chart."),
        GateRow("lnd_12_tangent_fd", "tangent finite difference", "validated", "Sector tangents agree with a centered xi difference.", f"maximum discrepancy={maximum('sector_tangent_finite_difference'):.3e}", "The common phase z rotates with xi."),
        GateRow("lnd_13_primitive", "determinant primitive", "validated", "The projected four-row finite defect enters the division-free determinant chart.", contract["determinant"], "The coefficient rows are synthetic regression data."),
        GateRow("lnd_14_primitive_tel", "primitive telescope", "validated", "The diagonal-first primitive owns every cross exactly once.", f"maximum discrepancy={maximum('primitive_telescope'):.3e}", "Reusing the original base is not allowed."),
        GateRow("lnd_15_current_tel", "current telescope", "validated", "The same ownership survives the xi derivative.", f"maximum discrepancy={maximum('current_telescope'):.3e}", "The normalized physical factor is not used in the fixture."),
        GateRow("lnd_16_current_fd", "current finite difference", "validated", "The analytic determinant current differentiates the computed primitive.", f"maximum discrepancy={maximum('primitive_current_finite_difference'):.3e}", "This is a centered floating-point test."),
        GateRow("lnd_17_ideal", "ideal relative channels", "validated", "The exact P_H^0 and P_T^0 cubics pass through the same stored basis.", "Both ideal coefficient vectors and all seven functional projections are recorded per fixture.", "They are relative channels, not four physical observation rows."),
        GateRow("lnd_18_scope", "fixture scope", "guard_validated", "The low-N tests validate code paths but are not physical Xi samples.", contract["fixture_scope"], "No sign extrapolation to N about 7.2e10 is admissible."),
        GateRow("lnd_19_next", "physical lift", "open", "Lift the actual physical rows through a compressed evaluator at the calibrated roots.", contract["next_action"], "Direct O(N) enumeration is forbidden at physical N."),
        GateRow("lnd_20_boundary", "proof boundary", "guard_validated", "The evaluator contract is not a proof of RH.", contract["proof_boundary"], "Every prize-level conclusion remains open."),
    ]
    require(len(rows) == 20, "gate row count")
    return rows


def render_note(payload: dict) -> str:
    contract = payload["exact"]
    summary = payload["summary"]
    lines = [
        "# Newman C1 Fixed-Cell Low-N Determinant Evaluator Gate",
        "",
        "Date: 2026-08-05",
        "",
        "Status: finite evaluator contract numerically cross-validated; physical-scale compression and signed determinant open; not a proof of RH.",
        "",
        "## Contract",
        "",
        contract["basis"],
        "",
        contract["blocks"],
        "",
        contract["remote"],
        "",
        contract["lift"],
        "",
        contract["determinant"],
        "",
        "## Audit Summary",
        "",
        f"- Fixtures: `{summary['fixtures']}`.",
        f"- Total finite modes: `{summary['total_modes']}`.",
        f"- Largest block/direct discrepancy: `{summary['maximum_block_direct_error']:.3e}`.",
        f"- Largest block/variation discrepancy: `{summary['maximum_block_variation_error']:.3e}`.",
        f"- Largest outer near/remote discrepancy: `{summary['maximum_outer_error']:.3e}`.",
        f"- Largest determinant-current finite-difference discrepancy: `{summary['maximum_current_fd_error']:.3e}`.",
        "",
        "The fixtures use `N=4,5,7`, retain `a^2-1<alpha_P<=a^2`, and stay at least `0.1` in alpha from every relevant half-open transfer. The stored basis includes the exact ideal relative cubics and a synthetic four-row determinant projection.",
        "",
        "## Scope",
        "",
        contract["fixture_scope"],
        "",
        "## Next Action",
        "",
        contract["next_action"],
        "",
        "## Pi Provenance",
        "",
        contract["pi_provenance"],
        "",
        "## Proof Boundary",
        "",
        contract["proof_boundary"],
        "",
        payload["success"],
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    set_below_normal_priority()
    sources = load_and_audit_sources()
    fixtures = [fixture_certificate(config) for config in FIXTURES]
    contract = exact_contract()
    rows = build_rows(contract, fixtures)
    summary = {
        "fixtures": len(fixtures),
        "total_modes": sum(row["roster"]["mode_count"] for row in fixtures),
        "basis_degree": MAX_DEGREE,
        "determinant_rows": 4,
        "maximum_block_direct_error": max(row["errors"]["block_vs_direct"] for row in fixtures),
        "maximum_block_variation_error": max(row["errors"]["block_vs_variation"] for row in fixtures),
        "maximum_outer_error": max(row["errors"]["outer_near_remote"] for row in fixtures),
        "maximum_current_fd_error": max(row["errors"]["primitive_current_finite_difference"] for row in fixtures),
        "physical_signed_bounds": 0,
    }
    success = (
        "built fixed-cell low-N determinant evaluator gate: "
        f"{len(rows)} rows, {len(fixtures)} fixtures, {summary['total_modes']} finite modes, "
        "3 independent finite-band representations, exact finite paired-remote kernel, "
        "primitive/current telescopes and xi finite differences validated, 0 physical signed bounds"
    )
    payload = {
        "kind": KIND,
        "date": "2026-08-05",
        "status": "fixed_cell_low_n_determinant_evaluator_contract_validated",
        "sources": sources,
        "configuration": {
            "max_degree": MAX_DEGREE,
            "gauss_order": GAUSS_ORDER,
            "quad_eps": QUAD_EPS,
            "alpha_step": ALPHA_STEP,
            "nu_cxi_at_centre": complex_payload(Z0),
            "coefficient_rows": [complex_payload(row) for row in COEFFICIENT_ROWS],
        },
        "exact": contract,
        "fixtures": fixtures,
        "rows": [asdict(row) for row in rows],
        "summary": summary,
        "proof_boundary": contract["proof_boundary"],
        "success": success,
    }
    atomic_write(RESULT_PATH, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(payload))
    print(success)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
