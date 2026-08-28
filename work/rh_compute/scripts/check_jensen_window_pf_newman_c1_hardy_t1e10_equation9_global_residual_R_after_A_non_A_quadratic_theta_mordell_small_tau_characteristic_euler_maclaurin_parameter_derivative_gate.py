#!/usr/bin/env python3
"""Altered-order replay of the characteristic parameter-derivative gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
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
    "non_A_quadratic_theta_mordell_small_tau_characteristic_euler_maclaurin_parameter_derivative_gate"
)
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
CHECKER = Path(__file__).resolve()
BASE_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_small_tau_characteristic_euler_maclaurin_remainder_gate"
)
BASE_RESULT = REPO_ROOT / f"work/rh_compute/results/{BASE_STEM}.json"
BASE_BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{BASE_STEM}.py"
BASE_CHECKER = REPO_ROOT / f"work/rh_compute/scripts/check_{BASE_STEM}.py"

N = 496_283
T_MAX = Fraction(1, 2 * N)
RIGHT = (Fraction(95_747, 6), Fraction(-1, 3), 1)
LEFT = (Fraction(47_873, 3), Fraction(-2, 3), -1)
ORDER = 80
PRECISION_BITS = 512
BLOCK_SIZE = 13
SAMPLES = (Fraction(0), Fraction(1, 16), Fraction(1, 4), Fraction(1, 2), Fraction(1))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_builder():
    spec = importlib.util.spec_from_file_location("characteristic_parameter_derivative_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load production builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def rat(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def cis(value: Fraction) -> acb:
    return (2 * acb.pi() * acb(0, 1) * acb(rat(value % 1))).exp()


def pair(side: str) -> tuple[Fraction, Fraction, int]:
    return RIGHT if side == "right" else LEFT


def envelope_reversed(order: int, d: Fraction) -> arb:
    total = arb(0)
    factorial = math.factorial(order)
    two_pi = 2 * arb.pi()
    for j in reversed(range(order // 2 + 1)):
        coefficient = factorial // (math.factorial(j) * math.factorial(order - 2 * j))
        total += coefficient * (two_pi * rat(d)) ** (order - 2 * j) * (two_pi * rat(T_MAX)) ** j
    return total


def masses(kind: str, p: Fraction) -> list[Fraction]:
    y0 = p
    y1 = p + N
    if kind == "one":
        return [Fraction(N)]
    if kind == "quadratic":
        return [
            (y1**3 - y0**3) / 3 - p * (y1**2 - y0**2) / 2,
            (y1**2 - p * y1) - (y0**2 - p * y0),
            Fraction(2 * N),
        ]
    if kind == "cubic":
        return [
            (y1**4 - y0**4) / 4 - p * p * (y1**2 - y0**2) / 2,
            (y1**3 - p * p * y1) - (y0**3 - p * p * y0),
            3 * (y1**2 - y0**2),
            Fraction(6 * N),
        ]
    raise RuntimeError("unknown mass kind")


def polynomial_remainder(side: str, kind: str) -> arb:
    p, b, _ = pair(side)
    d = max(abs(b), abs(b + 2 * T_MAX * (p + N)))
    require(d < Fraction(7, 10), "derivative strip drift")
    m = 2 * ORDER
    integral = arb(0)
    for r, mass in reversed(list(enumerate(masses(kind, p)))):
        integral += math.comb(m, r) * rat(mass) * envelope_reversed(m - r, d)
    return 2 * arb(m).zeta() * integral / (2 * arb.pi()) ** m


def derivative_remainders(side: str) -> tuple[arb, arb]:
    rem_one = polynomial_remainder(side, "one")
    rem_quadratic = polynomial_remainder(side, "quadratic")
    rem_cubic = polynomial_remainder(side, "cubic")
    return (
        2 * arb.pi() * rem_cubic,
        rem_one + 4 * arb.pi() * (1 + rat(T_MAX)) * rem_quadratic,
    )


def em_polynomial_reversed(module, p, b, t, coefficients) -> acb:
    corrections = acb(0)
    for r in reversed(range(1, ORDER + 1)):
        derivative_order = 2 * r - 1
        coefficient = acb(arb(fmpq.bernoulli(2 * r)) / math.factorial(2 * r))
        corrections += coefficient * (
            module.endpoint_derivative(derivative_order, p, b, t, coefficients, N)
            - module.endpoint_derivative(derivative_order, p, b, t, coefficients, 0)
        )
    return (
        module.polynomial_integral(p, b, t, coefficients)
        + (
            module.endpoint_derivative(0, p, b, t, coefficients, 0)
            + module.endpoint_derivative(0, p, b, t, coefficients, N)
        )
        / 2
        + corrections
    )


def direct_reversed(side: str, t: Fraction) -> tuple[acb, acb, acb]:
    p, b, _ = pair(side)
    a = b + 2 * t * p
    step = cis(2 * t)
    one = acb(0)
    quadratic = acb(0)
    cubic = acb(0)
    starts = list(range(0, N + 1, BLOCK_SIZE))
    for start in reversed(starts):
        stop = min(N + 1, start + BLOCK_SIZE)
        term = cis(a * start + t * start * start)
        ratio = cis(a + t * (2 * start + 1))
        block_one = acb(0)
        block_quadratic = acb(0)
        block_cubic = acb(0)
        for k in range(start, stop):
            y = p + k
            block_one += term
            block_quadratic += acb(rat(y * y - p * y)) * term
            block_cubic += acb(rat(y**3 - p * p * y)) * term
            term *= ratio
            ratio *= step
        one += block_one
        quadratic += block_quadratic
        cubic += block_cubic
    return one, quadratic, cubic


def main() -> int:
    require(RESULT.is_file() and BUILDER.is_file(), "missing production derivative artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production derivative artifact did not pass")
    require(artifact["scope"]["euler_maclaurin_order"] == ORDER, "EM order drift")
    require(artifact["decision"]["current_level_transport_promoted_to_source_cover"] is False, "diagnostic overpromoted")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(CHECKER), "checker hash drift")
    require(artifact["dependencies"]["base_result"]["sha256"] == file_hash(BASE_RESULT), "base result hash drift")
    require(artifact["dependencies"]["base_builder"]["sha256"] == file_hash(BASE_BUILDER), "base builder hash drift")
    require(artifact["dependencies"]["base_checker"]["sha256"] == file_hash(BASE_CHECKER), "base checker hash drift")

    flint.ctx.prec = PRECISION_BITS
    module = load_builder()
    remainders = {side: derivative_remainders(side) for side in ("right", "left")}
    for side, (rem_t, rem_s) in remainders.items():
        require(rem_t < arb("0.05"), f"{side} altered t remainder misses target")
        require(rem_s < arb("0.000001"), f"{side} altered s remainder misses target")
        reported = artifact["uniform_derivative_remainder_contract"][side]
        ratio_t = arb(reported["partial_t_W_remainder_upper"]) / rem_t
        ratio_s = arb(reported["partial_s_W_remainder_upper"]) / rem_s
        require(arb("0.9999999999") < ratio_t < arb("1.0000000001"), f"{side} t remainder drift")
        require(arb("0.9999999999") < ratio_s < arb("1.0000000001"), f"{side} s remainder drift")

    rows = {(row["side"], row["t_over_first_wall"]): row for row in artifact["derivative_overlap_matrix"]}
    require(len(rows) == 10, "production derivative matrix size drift")
    for side in ("right", "left"):
        p, b, epsilon = pair(side)
        one_coefficients = (Fraction(1),)
        quadratic_coefficients = (Fraction(0), -p, Fraction(1))
        cubic_coefficients = (Fraction(0), -p * p, Fraction(0), Fraction(1))
        for sample in SAMPLES:
            require((side, str(sample)) in rows, f"missing production derivative row {side} {sample}")
            t = sample * T_MAX
            em_one = em_polynomial_reversed(module, p, b, t, one_coefficients)
            em_quadratic = em_polynomial_reversed(module, p, b, t, quadratic_coefficients)
            em_cubic = em_polynomial_reversed(module, p, b, t, cubic_coefficients)
            em_t = 2 * acb.pi() * acb(0, 1) * em_cubic
            em_s = epsilon * (
                em_one + 4 * acb.pi() * acb(0, 1) * (1 + acb(rat(t))) * em_quadratic
            )
            direct_one, direct_quadratic, direct_cubic = direct_reversed(side, t)
            direct_t = 2 * acb.pi() * acb(0, 1) * direct_cubic
            direct_s = epsilon * (
                direct_one
                + 4 * acb.pi() * acb(0, 1) * (1 + acb(rat(t))) * direct_quadratic
            )
            rem_t, rem_s = remainders[side]
            require((em_t + acb(arb(0, rem_t.abs_upper()), arb(0, rem_t.abs_upper()))).overlaps(direct_t), f"altered t overlap failed {side} {sample}")
            require((em_s + acb(arb(0, rem_s.abs_upper()), arb(0, rem_s.abs_upper()))).overlaps(direct_s), f"altered s overlap failed {side} {sample}")

    print("independently replayed characteristic t and s derivatives with ten direct overlaps", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
