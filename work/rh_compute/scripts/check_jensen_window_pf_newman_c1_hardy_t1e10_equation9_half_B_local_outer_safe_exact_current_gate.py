#!/usr/bin/env python3
"""Independent replay of the exact local B-window signed current."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp
import numpy as np
from scipy.special import erfc, roots_legendre


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_outer_safe_exact_current_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

T = 10_000_000_000
B = 5_122_421
SPLITS = {621: 0.00024238523659, 622: 0.000242916437095}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ball_bounds(text: str) -> tuple[mp.mpf, mp.mpf]:
    match = re.fullmatch(r"\[([^\]]+) \+/- ([^\]]+)\]", text)
    require(match is not None, f"unparsed Arb ball: {text}")
    midpoint, radius = (mp.mpf(value) for value in match.groups())
    return midpoint - radius, midpoint + radius


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    cache_path = REPO_ROOT / artifact["sources"]["row_cache"]["path"]
    require(cache_path.is_file() and file_hash(cache_path) == artifact["sources"]["row_cache"]["sha256"], "row-cache hash drift")

    pi, imaginary_unit = np.pi, 1j
    t, endpoint = float(T), float(B)
    x0 = (1 - np.sqrt(1 - 8 * t / (pi * endpoint**2))) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_hessian = np.sqrt(hessian)
    full = (1 + imaginary_unit) / 2

    mp.mp.dps = 80
    t_mp, endpoint_mp = mp.mpf(T), mp.mpf(B)
    x0_mp = (1 - mp.sqrt(1 - 8 * t_mp / (mp.pi * endpoint_mp**2))) / 2
    phase0 = mp.pi * endpoint_mp**2 * x0_mp / 4 + t_mp * mp.log((1 - x0_mp) / x0_mp) / 2
    rotation = complex(mp.exp(1j * (phase0 - mp.pi / 8)))
    normalizer = (pi / (32 * t)) ** 0.25

    def phase_relative(x: np.ndarray) -> np.ndarray:
        delta = x - x0
        return pi * endpoint**2 * delta / 4 + t * (
            np.log1p(-delta / (1 - x0)) - np.log1p(delta / x0)
        ) / 2

    def direct_terms(mode: int, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        c = endpoint * x / 2
        positive = endpoint / (imaginary_unit * pi * (endpoint * x - 2 * mode)) - mode / (2 * pi**2 * (c - mode) ** 3) + 3 * imaginary_unit * mode * x / (4 * pi**3 * (c - mode) ** 5)
        negative = endpoint / (imaginary_unit * pi * (endpoint * x + 2 * mode)) + mode / (2 * pi**2 * (c + mode) ** 3) - 3 * imaginary_unit * mode * x / (4 * pi**3 * (c + mode) ** 5)
        return positive, negative

    def exact_positive(mode: int, x: np.ndarray, side: int) -> np.ndarray:
        q = np.sqrt(x / 2) * (endpoint - 2 * mode / x)
        exponential = np.exp(imaginary_unit * pi * q**2 / 2)
        argument = np.exp(-imaginary_unit * pi / 4) * np.sqrt(pi / 2) * np.abs(q)
        tail = full * erfc(argument)
        amplitude = mode * np.sqrt(2 / x)
        sign_adapted = exponential / (imaginary_unit * pi) - side * amplitude * tail
        step = ((1 if mode <= 621 else 0) - (1 if side < 0 else 0)) * amplitude * (1 + imaginary_unit)
        return (sign_adapted + step) / (x * exponential)

    def integrand(xi: np.ndarray, mode: int, regime: str, side: int) -> np.ndarray:
        x = x0 + xi / root_hessian
        positive, negative = direct_terms(mode, x)
        coefficient = positive + negative if regime == "normal" else exact_positive(mode, x, side) + negative
        weight = (x * (1 - x)) ** -0.25
        return weight * np.exp(imaginary_unit * phase_relative(x)) * coefficient / root_hessian

    nodes, weights = roots_legendre(96)

    def integrate_interval(left: float, right: float, mode: int, regime: str, side: int) -> complex:
        panel_count = max(1, int(np.ceil((right - left) / 0.125)))
        edges = np.linspace(left, right, panel_count + 1)
        total = 0j
        for a, b in zip(edges[:-1], edges[1:]):
            points = (a + b) / 2 + (b - a) * nodes / 2
            total += (b - a) * np.dot(weights, integrand(points, mode, regime, side)) / 2
        return total

    reduced = 0j
    for mode in (621, 622):
        crossing = (2 * mode / endpoint - x0) * root_hessian
        split = (SPLITS[mode] - x0) * root_hessian
        if mode == 621:
            intervals = [(-70, crossing, "outer", -1), (crossing, 70, "outer", 1)]
        else:
            intervals = [(-70, crossing, "outer", -1), (crossing, split, "outer", 1), (split, 70, "normal", 1)]
        reduced += sum(integrate_interval(a, b, mode, regime, side) for a, b, regime, side in intervals)

    physical = 2 * normalizer * (rotation * reduced).real
    require(-1.39e-4 < physical < -1.37e-4, "independent Gauss replay left the signed bracket")
    stored = artifact["interval_certificate"]
    stored_low, stored_high = ball_bounds(stored["combined_principal_physical_signed_projection_ball"])
    require(stored_low < physical < stored_high or abs(physical - float((stored_low + stored_high) / 2)) < 5e-10, "independent replay disagrees with stored interval")
    require(ball_bounds(stored["exact_local_B_window_physical_upper_ball"])[1] < mp.mpf("-1.369e-4"), "stored exact upper bound is not negative enough")

    decision = artifact["decision"]
    require(decision["target_step_and_Fresnel_tail_kept_joined"] is True, "crossing guard drift")
    require(decision["outside_window_grouped_trace_tails_bounded"] is False, "outside-window overclaim")
    require(decision["complete_B_face_estimate_proved"] is False, "complete-B overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print(f"independent local B current replay: physical={physical:.15e}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
