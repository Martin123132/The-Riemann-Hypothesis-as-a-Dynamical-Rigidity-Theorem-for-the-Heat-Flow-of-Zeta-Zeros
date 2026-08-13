#!/usr/bin/env python3
"""Independently check the fixed-cell low-N determinant evaluator gate."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import quad_vec


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_fixed_cell_low_n_determinant_evaluator_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
MAX_DEGREE = 6
CHECK_EPS = 8.0e-11
BASIS_TOLERANCE = 2.0e-8


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decode_complex(value: dict) -> complex:
    return complex(float(value["re"]), float(value["im"]))


def decode_vector(values: list[dict]) -> np.ndarray:
    return np.array([decode_complex(value) for value in values], dtype=np.complex128)


def powers(lam: np.ndarray | float) -> np.ndarray:
    values = np.asarray(lam)
    return np.stack([values**degree for degree in range(MAX_DEGREE + 1)], axis=-1)


def source_basis(u: float, alpha: float, t_value: float, sigma: float) -> np.ndarray:
    lam = math.log(u)
    exponent = t_value * lam**2 / 4.0 - sigma * lam + 2j * np.pi * alpha * lam
    return np.exp(exponent) * powers(lam)


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


def cell_bounds(q_value: int, n_value: int) -> tuple[float, float]:
    return max(1.0, q_value - 0.5), min(float(n_value), q_value + 0.5)


def integrate_cells(function: Callable[[float], np.ndarray], n_value: int) -> np.ndarray:
    total = np.zeros(MAX_DEGREE + 1, dtype=np.complex128)
    for q_value in range(1, n_value + 1):
        left, right = cell_bounds(q_value, n_value)
        value, _ = quad_vec(
            function, left, right, epsabs=CHECK_EPS, epsrel=CHECK_EPS, limit=700
        )
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


def adaptive_block_sectors(
    n_value: int,
    alpha: float,
    t_value: float,
    sigma: float,
    lower: int,
    upper: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[int, ...]]:
    modes = np.arange(lower, upper + 1, dtype=np.int64)
    labels = np.floor(alpha / modes + 0.5).astype(np.int64)
    sectors = [
        np.zeros(MAX_DEGREE + 1, dtype=np.complex128),
        np.zeros(MAX_DEGREE + 1, dtype=np.complex128),
        np.zeros(MAX_DEGREE + 1, dtype=np.complex128),
    ]
    for q_value in range(1, n_value + 1):
        left, right = cell_bounds(q_value, n_value)
        for p_value in range(1, n_value + 1):
            cell_modes = modes[labels == p_value].astype(np.float64)
            if not len(cell_modes):
                continue

            def integrand(u: float, selected=cell_modes) -> np.ndarray:
                kernel = np.exp(-2j * np.pi * selected * u).sum()
                return source_basis(u, alpha, t_value, sigma) * kernel

            value, _ = quad_vec(
                integrand,
                left,
                right,
                epsabs=CHECK_EPS,
                epsrel=CHECK_EPS,
                limit=700,
            )
            distance = abs(q_value - p_value)
            sectors[0 if distance == 0 else 1 if distance == 1 else 2] += value
    return sectors[0], sectors[1], sectors[2], tuple(int(value) for value in labels)


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


def lift_coefficients(coefficients: np.ndarray, log_a: float) -> np.ndarray:
    lifted = np.zeros(MAX_DEGREE + 1, dtype=np.complex128)
    for degree, coefficient in enumerate(coefficients):
        lifted[degree] -= 1j * log_a * coefficient
        lifted[degree + 1] += 1j * coefficient
    return lifted


def observation(
    functional: np.ndarray,
    coefficients: np.ndarray,
    z_value: complex,
    log_a: float,
    *,
    tangent: bool,
) -> np.ndarray:
    result = np.zeros(4, dtype=np.float64)
    for row_index, row in enumerate(coefficients):
        selected = lift_coefficients(row, log_a) if tangent else row
        result[row_index] = np.real(z_value * np.dot(selected, functional[: len(selected)]))
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


def check_fixture(payload: dict, coefficients: np.ndarray, z_value: complex) -> dict:
    parameters = payload["parameters"]
    n_value = int(parameters["N"])
    a_value = float(parameters["a"])
    alpha = float(parameters["alpha"])
    t_value = float(parameters["t"])
    sigma = float(parameters["sigma"])
    lower = math.floor(alpha / (n_value + 0.5)) + 1
    upper = math.floor(2.0 * alpha)
    diagonal, adjacent, far, labels = adaptive_block_sectors(
        n_value, alpha, t_value, sigma, lower, upper
    )
    carrier = carrier_basis(n_value, alpha, t_value, sigma)
    near = near_basis(n_value, alpha, t_value, sigma, lower, upper)
    remote = remote_basis(n_value, alpha, t_value, sigma, upper)
    computed = {
        "carrier": carrier,
        "diagonal": diagonal,
        "adjacent": adjacent,
        "far": far,
        "band_blocks": diagonal + adjacent + far,
        "near": near,
        "remote": remote,
    }
    discrepancies = {}
    for name, value in computed.items():
        stored = decode_vector(payload["basis"][name])
        discrepancies[name] = float(np.max(np.abs(value - stored)))
        require(discrepancies[name] < BASIS_TOLERANCE, f"{payload['id']} {name} drifted")
    require(labels == tuple(payload["roster"]["labels"]), f"{payload['id']} labels drifted")
    outer_discrepancy = float(
        np.max(np.abs(carrier - (diagonal + adjacent + far) - near - remote))
    )
    require(outer_discrepancy < BASIS_TOLERANCE, f"{payload['id']} outer identity drifted")

    log_a = math.log(a_value)
    values = {
        name: observation(computed[name], coefficients, z_value, log_a, tangent=False)
        for name in ("carrier", "diagonal", "adjacent", "far")
    }
    tangents = {
        name: observation(computed[name], coefficients, z_value, log_a, tangent=True)
        for name in ("carrier", "diagonal", "adjacent", "far")
    }
    edge = np.array([float(value) for value in payload["determinant"]["edge"]])
    p_value = values["carrier"]
    p_dot = tangents["carrier"]
    d0 = values["diagonal"] - p_value
    d1 = values["adjacent"]
    d2 = values["far"]
    dd0 = tangents["diagonal"] - p_dot
    dd1 = tangents["adjacent"]
    dd2 = tangents["far"]
    base = edge + p_value
    defect = d0 + d1 + d2
    defect_dot = dd0 + dd1 + dd2
    primitive = primitive_increment(base, defect)
    current = current_increment(base, defect, p_dot, defect_dot)
    primitive_parts = (
        primitive_increment(base, d0),
        primitive_increment(base - d0, d1),
        primitive_increment(base - d0 - d1, d2),
    )
    current_parts = (
        current_increment(base, d0, p_dot, dd0),
        current_increment(base - d0, d1, p_dot - dd0, dd1),
        current_increment(base - d0 - d1, d2, p_dot - dd0 - dd1, dd2),
    )
    stored_primitive = float(payload["determinant"]["primitive"])
    stored_current = float(payload["determinant"]["current"])
    require(abs(primitive - stored_primitive) < 5.0e-8, f"{payload['id']} primitive drifted")
    require(abs(current - stored_current) < 2.0e-7, f"{payload['id']} current drifted")
    require(abs(primitive - sum(primitive_parts)) < 5.0e-12, "primitive telescope failed")
    require(abs(current - sum(current_parts)) < 5.0e-11, "current telescope failed")
    require(
        payload["finite_difference"]["assignments_unchanged"],
        f"{payload['id']} finite-difference assignment guard missing",
    )
    require(
        payload["errors"]["primitive_current_finite_difference"] < 2.0e-8,
        f"{payload['id']} stored current finite difference too large",
    )
    return {
        "basis_error": max(discrepancies.values()),
        "outer_error": outer_discrepancy,
        "primitive_error": abs(primitive - stored_primitive),
        "current_error": abs(current - stored_current),
    }


def main() -> int:
    require(RESULT_PATH.is_file(), f"missing result: {RESULT_PATH}")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind drifted")
    require(payload["status"] == "fixed_cell_low_n_determinant_evaluator_contract_validated", "status drifted")
    require(len(payload["rows"]) == 20, "row count drifted")
    require(len(payload["fixtures"]) == 3, "fixture count drifted")
    require(payload["summary"]["total_modes"] == 184, "mode count drifted")
    require(payload["summary"]["physical_signed_bounds"] == 0, "signed-bound count drifted")
    require("not interval arithmetic" in payload["proof_boundary"], "proof boundary drifted")
    require("No geometric or fitted pi" in payload["exact"]["pi_provenance"], "pi provenance drifted")
    for source in payload["sources"].values():
        path = REPO_ROOT / source["path"]
        require(path.is_file(), f"missing audited source: {path}")
        require(file_hash(path) == source["sha256"], f"source hash drifted: {path}")

    coefficients = np.array(
        [decode_vector(row) for row in payload["configuration"]["coefficient_rows"]],
        dtype=np.complex128,
    )
    require(coefficients.shape == (4, 6), "coefficient shape drifted")
    z_value = decode_complex(payload["configuration"]["nu_cxi_at_centre"])
    require(abs(abs(z_value) - 1.0) < 1.0e-15, "centre phase is not unit")
    audits = [check_fixture(fixture, coefficients, z_value) for fixture in payload["fixtures"]]
    success = (
        "validated fixed-cell low-N determinant evaluator gate: "
        f"{len(payload['rows'])} rows, {len(audits)} fixtures, "
        f"{payload['summary']['total_modes']} finite modes, "
        f"maximum alternate-basis error {max(row['basis_error'] for row in audits):.3e}, "
        f"maximum alternate-current error {max(row['current_error'] for row in audits):.3e}, "
        "0 physical signed bounds"
    )
    print(success)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
