#!/usr/bin/env python3
"""Independently check the grouped selected-branch leading-defect gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad_vec
from scipy.special import airy


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_selected_branch_grouped_integral_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
C = 159_577
Q = 39_894
L = 39_696
U = 40_094


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ball_upper(text: str) -> mp.mpf:
    body = text.strip("[] ")
    if "+/-" not in body:
        return mp.mpf(body)
    midpoint, radius = body.split("+/-", 1)
    return mp.mpf(midpoint) + mp.mpf(radius)


def ball_midpoint(text: str) -> float:
    return float(text.strip("[] ").split("+/-", 1)[0])


def beta4_weights(lam: float, y: float, beta: float) -> tuple[float, float]:
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
    ai, ai_prime, bi, bi_prime = airy(-(lam + y))
    w = ai - 1j * sigma * bi
    w_x = -ai_prime + 1j * sigma * bi_prime
    u, v = beta4_weights(lam, y, beta)
    return (u * w + v * w_x) / 2


def action_defect(s: np.ndarray, mu: float) -> np.ndarray:
    root = np.sqrt(2 * s + mu - 1)
    return (
        3 * s * s - 8 * s * root + 12 * s - 4 * mu * root
        + 6 * mu * np.log(s) - 3 * mu * np.log(1 - mu) + 9 * mu
        + 4 * root - 6 * np.log(s) + 3 * np.log(1 - mu) - 11
    ) / 6


def independent_grouped_value(fraction: float) -> complex:
    pi = np.pi
    tstar = pi * C * C / 8
    beta = tstar ** (1 / 3)
    mu = fraction * pi / (16 * tstar)
    lam = beta * beta * mu
    y_max = pi * C / (2 * beta)
    modes = np.array([mode for mode in range(L, U + 1) if mode != Q], dtype=float)
    s = 4 * modes / C
    sigma = np.sign(4 * modes - C).astype(int)
    detuning = beta * (s - 1)
    ratio = ((1 - mu) * (2 * s + mu - 1)) ** 0.25 / np.sqrt(s)
    coefficients = (ratio * np.exp(1j * beta**3 * action_defect(s, mu)) - 1) / np.sqrt(modes)

    def integrand(y: float) -> complex:
        minus = branch_value(lam, y, beta, -1)
        plus = minus.conjugate()
        branches = np.where(sigma < 0, minus, plus)
        return np.exp(1j * y * y / (4 * beta)) * np.sum(
            coefficients * branches * np.exp(-1j * detuning * y)
        )

    value, error = quad_vec(integrand, 0.0, y_max, epsabs=2e-12, epsrel=2e-12, limit=1200)
    require(float(error) < 3e-11, "independent grouped quadrature failed")
    return complex(value)


def main() -> None:
    require(RESULT.is_file(), "missing result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "stored gate is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")
    for source in artifact["sources"].values():
        path = REPO_ROOT / source["path"]
        require(path.is_file(), f"missing source: {path}")
        require(file_hash(path) == source["sha256"], f"source hash drift: {path}")
    mp.mp.dps = 60
    c = artifact["certificate"]
    require(c["y_panel_count"] == 65536, "panel count drift")
    require(c["height_node_count"] == 5, "height node count drift")
    maximum = max(ball_upper(row["absolute_ball"]) for row in c["node_rows"])
    require(abs(maximum - ball_upper(c["maximum_node_absolute_ball"])) < mp.mpf("1e-55"), "maximum node mismatch")
    require(maximum < mp.mpf("3.3e-5"), "node maximum failed")
    require(ball_upper(c["y_midpoint_error_ball"]) < mp.mpf("5e-6"), "midpoint error failed")
    require(ball_upper(c["uniform_selected_weighted_bound_ball"]) < mp.mpf("0.001"), "uniform bound failed")
    midpoint_error = float(ball_upper(c["y_midpoint_error_ball"]))
    for row in c["node_rows"]:
        diagnostic = independent_grouped_value(row["height_fraction_index_over_4"] / 4)
        stored = complex(
            ball_midpoint(row["value_ball"]["real_ball"]),
            ball_midpoint(row["value_ball"]["imag_ball"]),
        )
        require(abs(diagnostic - stored) < midpoint_error + 3e-10, "independent grouped value misses certified midpoint enclosure")
    require(artifact["decision"]["local_normalized_headroom_closed"] is False, "proof boundary drift")
    print("independently checked grouped selected-branch leading-defect integral", flush=True)


if __name__ == "__main__":
    main()
