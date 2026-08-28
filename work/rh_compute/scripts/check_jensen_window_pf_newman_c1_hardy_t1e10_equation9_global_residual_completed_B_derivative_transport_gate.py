#!/usr/bin/env python3
"""Numerically and algebraically replay the completed B derivative transport."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_derivative_transport_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    import mpmath as mp

    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    mp.mp.dps = 75
    physical_endpoint = mp.mpf(5_122_421)
    rho = 1 / (1j * mp.pi)

    def fresnel(q: mp.mpf) -> mp.mpc:
        return mp.fresnelc(q) + 1j * mp.fresnels(q)

    def y_value(x: mp.mpf, endpoint: mp.mpf, mode: int, sign: int, tau: int) -> mp.mpc:
        a = mp.mpf(mode) * mp.sqrt(2 / x)
        b = endpoint * mp.sqrt(x / 2)
        a_sign = sign * a
        q = b - sign * a
        phase = mp.e ** (1j * mp.pi * q**2 / 2)
        p_b = phase / (1j * mp.pi) - a_sign * ((1 + 1j) / 2 - fresnel(q))
        if sign == 1:
            p_b += (1 - tau) * a * (1 + 1j)
        return p_b / phase

    def transport(
        x: mp.mpf, endpoint: mp.mpf, mode: int, sign: int, tau: int
    ) -> tuple[mp.mpc, mp.mpc, mp.mpc]:
        a = mp.mpf(mode) * mp.sqrt(2 / x)
        b = endpoint * mp.sqrt(x / 2)
        a_sign = sign * a
        q = b - sign * a
        q_first = (b + sign * a) / (2 * x)
        q_second = (-b - 3 * sign * a) / (4 * x**2)
        kappa = 1 / (2 * x) + 1j * mp.pi * q * q_first
        kappa_first = -1 / (2 * x**2) + 1j * mp.pi * (q_first**2 + q * q_second)
        y = y_value(x, endpoint, mode, sign, tau)
        y_first = a_sign * q_first + (rho - y) * kappa
        y_second = a_sign * (q_second - q_first / (2 * x)) - y_first * kappa + (rho - y) * kappa_first
        return y, y_first, y_second

    test_rows = (
        (physical_endpoint, 620, 0, mp.mpf(1240) / physical_endpoint),
        (physical_endpoint, 623, 1, mp.mpf(1246) / physical_endpoint),
        (mp.mpf(101), 7, 0, mp.mpf("0.2")),
    )
    tolerance = mp.mpf("1e-48")
    for endpoint, mode, tau, x in test_rows:
        values: dict[int, tuple[mp.mpc, mp.mpc, mp.mpc]] = {}
        for sign in (1, -1):
            y, y_first, y_second = transport(x, endpoint, mode, sign, tau)
            numeric_first = mp.diff(lambda xx: y_value(xx, endpoint, mode, sign, tau), x, 1)
            numeric_second = mp.diff(lambda xx: y_value(xx, endpoint, mode, sign, tau), x, 2)
            first_scale = max(mp.mpf(1), abs(y_first), abs(numeric_first))
            second_scale = max(mp.mpf(1), abs(y_second), abs(numeric_second))
            require(abs(y_first - numeric_first) / first_scale < tolerance, f"Y first derivative mismatch for m={mode}, s={sign}")
            require(abs(y_second - numeric_second) / second_scale < tolerance, f"Y second derivative mismatch for m={mode}, s={sign}")
            values[sign] = (y, y_first, y_second)

        y_sum = values[1][0] + values[-1][0]
        y_first_sum = values[1][1] + values[-1][1]
        y_second_sum = values[1][2] + values[-1][2]
        pair = y_sum / x
        pair_first = y_first_sum / x - y_sum / x**2
        pair_second = y_second_sum / x - 2 * y_first_sum / x**2 + 2 * y_sum / x**3

        def pair_value(xx: mp.mpf) -> mp.mpc:
            return (
                y_value(xx, endpoint, mode, 1, tau)
                + y_value(xx, endpoint, mode, -1, tau)
            ) / xx

        numeric_pair_first = mp.diff(pair_value, x, 1)
        numeric_pair_second = mp.diff(pair_value, x, 2)
        require(abs(pair_first - numeric_pair_first) / max(1, abs(pair_first)) < tolerance, f"pair first derivative mismatch for m={mode}")
        require(abs(pair_second - numeric_pair_second) / max(1, abs(pair_second)) < tolerance, f"pair second derivative mismatch for m={mode}")
        require(abs(pair - pair_value(x)) < tolerance * max(1, abs(pair)), f"pair value mismatch for m={mode}")

    decision = artifact["decision"]
    require(decision["phase_stripped_completed_mode_transport_proved"] is True, "mode transport missing")
    require(decision["bulk_completion_is_homogeneous_transport_solution"] is True, "homogeneous completion missing")
    require(decision["exact_first_second_pair_derivatives_proved"] is True, "pair derivatives missing")
    require(decision["uniform_M_epsilon_derivative_bounds_proved"] is False, "uniform derivative bound overclaim")
    require(decision["joined_remainder_below_1p4058e_minus_4_proved"] is False, "joined bound overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently checked completed B first/second derivative transport at crossing and exterior witnesses", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
