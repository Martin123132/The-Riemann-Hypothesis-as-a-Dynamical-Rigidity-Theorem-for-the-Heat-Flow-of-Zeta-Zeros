#!/usr/bin/env python3
"""Non-rigorous scout of the grouped selected beta^-4 finite-integral remainder."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import time

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import numpy as np
from scipy.integrate import quad_vec
from scipy.special import airy


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/exploratory_selected_beta4_finite_integral_remainder.json"
C = 159_577
Q = (C - 1) // 4
MODE_LO = 39_696
MODE_HI = 40_094


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


def weights(lam: float, y: float, beta: float) -> tuple[float, float]:
    x = lam + y
    epsilon = beta**-2
    u1 = -(13 * x - 10 * y) / 60
    v1 = (8 * x * x - 20 * x * y + 15 * y * y) / 60
    u2 = -(
        448 * x**5 - 2240 * x**4 * y + 4480 * x**3 * y**2
        - 4200 * x**2 * y**3 + 4565 * x**2 + 1575 * x * y**4
        - 9100 * x * y + 4130 * y**2
    ) / 50400
    v2 = (40 * x**3 - 140 * x**2 * y + 161 * x * y**2 - 70 * y**3 + 27) / 1680
    return 1 + epsilon * u1 + epsilon**2 * u2, epsilon * v1 + epsilon**2 * v2


def branch_value(lam: float, y: float, beta: float, sigma: int) -> complex:
    x = lam + y
    ai, aip, bi, bip = airy(-x)
    w = ai - 1j * sigma * bi
    wx = -aip + 1j * sigma * bip
    u, v = weights(lam, y, beta)
    return (u * w + v * wx) / 2


def calculate(mu_fraction: float) -> dict[str, object]:
    pi = math.pi
    c = float(C)
    tstar = pi * c * c / 8
    beta = tstar ** (1 / 3)
    mu_radius = pi / (16 * tstar)
    mu = mu_fraction * mu_radius
    lam = beta * beta * mu
    t = tstar * (1 - mu)
    y_max = pi * c / (2 * beta)
    modes = np.array([mode for mode in range(MODE_LO, MODE_HI + 1) if mode != Q], dtype=float)
    s = 4 * modes / c
    labels = 4 * modes - c
    sigma = np.sign(labels).astype(int)
    d = beta * (s - 1)
    root = np.sqrt(2 * s + mu - 1)
    saddle_y = 2 * beta * beta * (s - root)
    saddle_x = lam + saddle_y
    phase_second = 1 / (2 * beta) + sigma / (2 * np.sqrt(saddle_x))
    scale = np.sqrt(root / modes)

    def grouped_integrand(y: float) -> np.ndarray:
        common = np.exp(1j * y * y / (4 * beta))
        minus_value = branch_value(lam, y, beta, -1)
        plus_value = branch_value(lam, y, beta, 1)
        branches = np.where(sigma < 0, minus_value, plus_value)
        return common * scale * branches * np.exp(-1j * d * y)

    finite_vector, error = quad_vec(grouped_integrand, 0.0, y_max, epsabs=2e-11, epsrel=2e-11, limit=1200)
    finite = np.sum(finite_vector)
    stationary = 0j
    stationary_rows = []
    individual_modes = {39_696, 39_750, 39_794, 39_893, 39_895, 39_994, 40_050, 40_092, 40_093, 40_094}
    for index, mode in enumerate(modes.astype(int)):
        y = saddle_y[index]
        x = saddle_x[index]
        sig = int(sigma[index])
        xi = 2 * x**1.5 / 3
        carrier_phase = y * y / (4 * beta) - d[index] * y + sig * (xi - pi / 4)
        scaled_amplitude = np.exp(-1j * sig * (xi - pi / 4)) * branch_value(lam, y, beta, sig)
        lead = scaled_amplitude * np.exp(1j * carrier_phase)
        lead *= math.sqrt(2 * pi / abs(phase_second[index])) * np.exp(1j * np.sign(phase_second[index]) * pi / 4)
        lead *= scale[index]
        stationary += lead
        if mode in individual_modes:
            def one_integrand(y_value: float) -> complex:
                return (
                    scale[index]
                    * np.exp(1j * (y_value * y_value / (4 * beta) - d[index] * y_value))
                    * branch_value(lam, y_value, beta, sig)
                )

            individual_finite, individual_error = quad_vec(
                one_integrand,
                0.0,
                y_max,
                epsabs=2e-12,
                epsrel=2e-12,
                limit=1200,
            )
            stationary_rows.append({
                "mode": int(mode),
                "saddle_y": float(y),
                "scaled_finite_value": [float(individual_finite.real), float(individual_finite.imag)],
                "scaled_finite_quadrature_error_estimate": float(individual_error),
                "scaled_finite_stationary_lead": [float(lead.real), float(lead.imag)],
                "finite_minus_stationary_absolute": float(abs(individual_finite - lead)),
                "finite_plus_stationary_absolute": float(abs(individual_finite + lead)),
            })
    remainder = finite - stationary
    return {
        "mu_fraction": mu_fraction,
        "height": t,
        "lambda": lam,
        "selected_mode_count": int(modes.size),
        "finite_grouped_value": [float(finite.real), float(finite.imag)],
        "quadrature_error_estimate": float(error),
        "stationary_grouped_value": [float(stationary.real), float(stationary.imag)],
        "finite_minus_stationary": [float(remainder.real), float(remainder.imag)],
        "finite_minus_stationary_absolute": float(abs(remainder)),
        "stationary_rows": stationary_rows,
    }


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    rows = [calculate(fraction) for fraction in (0.0, 0.25, 0.5, 0.75, 1.0)]
    artifact = {
        "kind": "exploratory_selected_beta4_finite_integral_remainder",
        "status": "nonrigorous_diagnostic_only",
        "proof": False,
        "rows": rows,
        "maximum_finite_minus_stationary_absolute": max(row["finite_minus_stationary_absolute"] for row in rows),
        "warning": (
            "This scout tests a natural full-saddle normalization. It is not yet identified with the exact endpoint-retaining "
            "ordinary target and is not an interval enclosure. A large value rejects only this naive finite-integral comparison."
        ),
        "runtime": {"workers": 1, "process_priority": priority, "elapsed_seconds": round(time.perf_counter() - started, 3)},
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"maximum": artifact["maximum_finite_minus_stationary_absolute"], "rows": rows}, indent=2), flush=True)


if __name__ == "__main__":
    main()
