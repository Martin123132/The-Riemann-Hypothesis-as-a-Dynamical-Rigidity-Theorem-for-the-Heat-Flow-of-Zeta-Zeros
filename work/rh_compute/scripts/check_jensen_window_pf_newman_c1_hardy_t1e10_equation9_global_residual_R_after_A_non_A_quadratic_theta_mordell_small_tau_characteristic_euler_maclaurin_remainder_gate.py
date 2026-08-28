#!/usr/bin/env python3
"""Independently replay the small-t characteristic Euler--Maclaurin gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
import flint
from flint import acb, arb, fmpq


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_small_tau_characteristic_euler_maclaurin_remainder_gate"
)
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
CHECKER = Path(__file__).resolve()
GROUPED_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_small_tau_grouped_first_wall_bridge_gate"
)
GROUPED_RESULT = REPO_ROOT / f"work/rh_compute/results/{GROUPED_STEM}.json"
GROUPED_BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{GROUPED_STEM}.py"

N = 496_283
T_MAX = Fraction(1, 2 * N)
RIGHT = (Fraction(95_747, 6), Fraction(-1, 3))
LEFT = (Fraction(47_873, 3), Fraction(-2, 3))
ORDER = 48
PRECISION_BITS = 448
BLOCK_SIZE = 13
SAMPLES = (Fraction(0), Fraction(1, 16), Fraction(1, 4), Fraction(1, 2), Fraction(1))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rat(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def cis(value: Fraction) -> acb:
    return (2 * acb.pi() * acb(0, 1) * acb(rat(value % 1))).exp()


def pair(side: str) -> tuple[Fraction, Fraction]:
    require(side in ("right", "left"), "invalid side")
    return RIGHT if side == "right" else LEFT


def envelope_reversed(order: int, d: Fraction) -> arb:
    total = arb(0)
    two_pi = 2 * arb.pi()
    factorial = math.factorial(order)
    for j in reversed(range(order // 2 + 1)):
        coefficient = factorial // (math.factorial(j) * math.factorial(order - 2 * j))
        total += coefficient * (two_pi * rat(d)) ** (order - 2 * j) * (two_pi * rat(T_MAX)) ** j
    return total


def independent_remainder(side: str) -> arb:
    p, b = pair(side)
    endpoint_derivatives = (b, b + 2 * T_MAX * (p + N))
    require(endpoint_derivatives[0] < 0 < endpoint_derivatives[1], "stationary inventory drift")
    d = max(abs(endpoint_derivatives[0]), abs(endpoint_derivatives[1]))
    require(d < Fraction(7, 10), "frequency strip drift")
    m = 2 * ORDER
    h_m = envelope_reversed(m, d)
    h_previous = envelope_reversed(m - 1, d)
    integral_bound = rat(p * N + Fraction(N * N, 2)) * h_m + m * N * h_previous
    return 2 * arb(m).zeta() * integral_bound / (2 * arb.pi()) ** m


def e_factor_reversed(order: int, q_prime: Fraction, t: Fraction) -> acb:
    total = acb(0)
    i_two_pi = 2 * acb.pi() * acb(0, 1)
    factorial = math.factorial(order)
    for j in reversed(range(order // 2 + 1)):
        coefficient = factorial // (math.factorial(j) * math.factorial(order - 2 * j))
        total += (
            coefficient
            * (i_two_pi * acb(rat(q_prime))) ** (order - 2 * j)
            * (i_two_pi * acb(rat(t))) ** j
        )
    return total


def endpoint(order: int, p: Fraction, b: Fraction, t: Fraction, x: int) -> acb:
    y = p + x
    q_prime = b + 2 * t * y
    phase = t * x * x + (b + 2 * t * p) * x
    factor = acb(rat(y)) * e_factor_reversed(order, q_prime, t)
    if order:
        factor += order * e_factor_reversed(order - 1, q_prime, t)
    return cis(phase) * factor


def integral(p: Fraction, b: Fraction, t: Fraction) -> acb:
    if t == 0:
        lam = 2 * acb.pi() * acb(0, 1) * acb(rat(b))
        return cis(b * N) * (acb(rat(p + N)) / lam - 1 / lam**2) - (
            acb(rat(p)) / lam - 1 / lam**2
        )

    sqrt_t = rat(t).sqrt()
    z0 = p + b / (2 * t)
    z1 = p + N + b / (2 * t)
    u0 = acb(2 * sqrt_t * rat(z0))
    u1 = acb(2 * sqrt_t * rat(z1))
    f0 = u0.fresnel_c() + acb(0, 1) * u0.fresnel_s()
    f1 = u1.fresnel_c() + acb(0, 1) * u1.fresnel_s()
    j0 = (f1 - f0) / (2 * sqrt_t)
    j1 = (cis(t * z1 * z1) - cis(t * z0 * z0)) / (
        4 * acb.pi() * acb(0, 1) * acb(rat(t))
    ) - acb(rat(b / (2 * t))) * j0
    return cis(-t * p * p - b * p - b * b / (4 * t)) * j1


def em_reversed(side: str, t: Fraction) -> acb:
    p, b = pair(side)
    corrections = acb(0)
    for r in reversed(range(1, ORDER + 1)):
        derivative_order = 2 * r - 1
        coefficient = acb(arb(fmpq.bernoulli(2 * r)) / math.factorial(2 * r))
        corrections += coefficient * (
            endpoint(derivative_order, p, b, t, N)
            - endpoint(derivative_order, p, b, t, 0)
        )
    return integral(p, b, t) + (endpoint(0, p, b, t, 0) + endpoint(0, p, b, t, N)) / 2 + corrections


def direct_reversed_blocks(side: str, t: Fraction) -> acb:
    p, b = pair(side)
    a = b + 2 * t * p
    step = cis(2 * t)
    total = acb(0)
    starts = list(range(0, N + 1, BLOCK_SIZE))
    for start in reversed(starts):
        stop = min(N + 1, start + BLOCK_SIZE)
        term = cis(a * start + t * start * start)
        ratio = cis(a + t * (2 * start + 1))
        subtotal = acb(0)
        for k in range(start, stop):
            subtotal += acb(rat(p + k)) * term
            term *= ratio
            ratio *= step
        total += subtotal
    return total


def main() -> int:
    require(RESULT.is_file() and BUILDER.is_file(), "missing production artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact did not pass")
    require(artifact["scope"]["euler_maclaurin_order"] == ORDER, "EM order drift")
    require(artifact["decision"]["finite_parameter_cover_built"] is False, "proof boundary overpromoted")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(CHECKER), "checker hash drift")
    require(artifact["dependencies"]["grouped_wall_result"]["sha256"] == file_hash(GROUPED_RESULT), "grouped result hash drift")
    require(artifact["dependencies"]["grouped_wall_builder"]["sha256"] == file_hash(GROUPED_BUILDER), "grouped builder hash drift")

    flint.ctx.prec = PRECISION_BITS
    independent = {side: independent_remainder(side) for side in ("right", "left")}
    for side, remainder in independent.items():
        require(remainder < arb("0.001"), f"{side} independent remainder misses target")
        reported = arb(artifact["uniform_remainder_contract"][side]["remainder_radius_upper"])
        ratio = reported / remainder
        require(arb("0.9999999999") < ratio < arb("1.0000000001"), f"{side} remainder replay drift")

    rows = {(row["side"], row["t_over_first_wall"]): row for row in artifact["pointwise_feasibility_matrix"]}
    require(len(rows) == 10, "production sample matrix size drift")
    for side in ("right", "left"):
        for sample in SAMPLES:
            key = (side, str(sample))
            require(key in rows, f"missing production sample {key}")
            t = sample * T_MAX
            direct = direct_reversed_blocks(side, t)
            approximation = em_reversed(side, t)
            radius = independent[side].abs_upper()
            enclosure = approximation + acb(arb(0, radius), arb(0, radius))
            require(enclosure.overlaps(direct), f"independent EM enclosure misses {key}")

    print("independently replayed characteristic Euler-Maclaurin remainder and ten direct overlaps", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
