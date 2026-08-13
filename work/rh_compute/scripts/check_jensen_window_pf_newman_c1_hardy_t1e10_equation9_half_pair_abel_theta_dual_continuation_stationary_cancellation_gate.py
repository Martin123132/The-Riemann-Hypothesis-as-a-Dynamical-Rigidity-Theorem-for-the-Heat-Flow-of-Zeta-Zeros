#!/usr/bin/env python3
"""Independent replay of the dual-continuation stationary cancellation gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_continuation_stationary_cancellation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

T = 10_000_000_000
A = 159_577
B = 5_122_421
L = (B - A) // 2
K = (A - 1) // 2


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def arctan_partial_inverse(q: int, terms: int) -> Fraction:
    return sum(
        (Fraction(1, (2 * j + 1) * q ** (2 * j + 1)) if j % 2 == 0 else -Fraction(1, (2 * j + 1) * q ** (2 * j + 1)))
        for j in range(terms)
    )


def machin_pi_bounds() -> tuple[Fraction, Fraction]:
    # Even partial sums are lower bounds and odd partial sums are upper bounds
    # for arctan(x), 0<x<1.
    a_low = arctan_partial_inverse(5, 30)
    a_high = arctan_partial_inverse(5, 31)
    b_low = arctan_partial_inverse(239, 10)
    b_high = arctan_partial_inverse(239, 11)
    lower = 16 * a_low - 4 * b_high
    upper = 16 * a_high - 4 * b_low
    require(lower < upper, "Machin enclosure ordering failed")
    return lower, upper


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")

    # Replay in face-distance y=x-z coordinates, unlike the production
    # z-coordinate calculation.
    x, alpha, n = sp.symbols("x alpha n", positive=True, real=True)
    endpoint = alpha + 2 * n
    y = sp.symbols("y", positive=True, real=True)
    z = x - y
    r = y / (x * z)
    phase = sp.pi * endpoint**2 * z / 4 - sp.pi * n**2 / r
    y_star = 2 * n * x / endpoint
    require(sp.simplify(sp.diff(phase, y).subs(y, y_star)) == 0, "face-distance saddle failed")

    hessian = sp.factor(sp.simplify(sp.diff(phase, y, 2).subs(y, y_star)))
    expected_hessian = -sp.pi * endpoint**3 / (4 * n * x)
    require(sp.simplify(hessian - expected_hessian) == 0, "face-distance Hessian failed")

    z_star = sp.simplify(z.subs(y, y_star))
    r_star = sp.simplify(r.subs(y, y_star))
    multiplier = (
        r_star ** (-sp.Rational(1, 2))
        * z_star ** (-sp.Rational(1, 2))
        * (2 + sp.I * sp.pi * endpoint**2 * z_star)
        * sp.sqrt(2 * sp.pi / (-hessian))
    )
    expected = 4 * sp.sqrt(x) / endpoint + 2 * sp.I * sp.pi * alpha * x ** sp.Rational(3, 2)
    replay = sp.powsimp(sp.powdenest(multiplier, force=True), force=True)
    require(sp.simplify(replay - expected) == 0, "independent stationary multiplier failed")

    endpoint_1, endpoint_2 = sp.symbols("E_1 E_2", positive=True, real=True)
    paired = (
        4 * sp.sqrt(x) / endpoint_1
        + 2 * sp.I * sp.pi * alpha * x ** sp.Rational(3, 2)
        - 4 * sp.sqrt(x) / endpoint_2
        - 2 * sp.I * sp.pi * alpha * x ** sp.Rational(3, 2)
    )
    require(sp.simplify(paired - 4 * sp.sqrt(x) * (1 / endpoint_1 - 1 / endpoint_2)) == 0, "paired cancellation replay failed")

    # Certify A-2 < sqrt(8t/pi) < A without Arb or floating-point pi.
    pi_low, pi_high = machin_pi_bounds()
    eight_t = 8 * T
    require(pi_high * (A - 2) ** 2 < eight_t, "lower threshold inequality failed")
    require(eight_t < pi_low * A**2, "upper threshold inequality failed")
    derivative_margin_lower = Fraction(2 * T, 1) - pi_high * Fraction((A - 2) ** 2, 4)
    require(derivative_margin_lower > 479_304, "exact rational derivative margin failed")

    require(B - 2 * L == A, "continuation offset failed")
    require(A - 2 * K == 1, "positive-label endpoint failed")
    require(A - 2 * (K + 1) == -1, "negative-label start failed")
    for k in (1, 2, 97, 10_000, K):
        require(B - 2 * (L + k) == A - 2 * k, "matched continuation identity failed")
        require(A - 2 * k > 0, "positive continuation label failed")

    decision = artifact["decision"]
    require(decision["dual_continuation_saddle_roster_partition_complete"] is True, "partition decision drift")
    require(decision["matched_positive_continuation_pairs"] == K, "continuation count drift")
    require(decision["leading_universal_z_stationary_channel_cancels_for_k_1_through_K"] is True, "leading cancellation drift")
    require(decision["scalar_endpoint_channel_survives"] is True, "scalar-channel guard drift")
    require(decision["all_matched_positive_continuation_reduced_phases_are_x_nonstationary"] is True, "x-nonstationarity drift")
    require(decision["exceptional_k_0_face_fold_term_cancelled"] is False, "k=0 overclaim")
    require(decision["uniform_z_stationary_expansion_remainder_bounded"] is False, "remainder overclaim")
    require(decision["complete_modular_continuation_bounded"] is False, "continuation overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked dual-continuation cancellation in face-distance coordinates with rational pi bounds", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
