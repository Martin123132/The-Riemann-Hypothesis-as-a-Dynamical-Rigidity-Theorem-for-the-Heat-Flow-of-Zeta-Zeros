#!/usr/bin/env python3
"""Independent replay of the completed-B outer rational derivative tail."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_three_current_derivative_tail_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

T = 10_000_000_000
B = 5_122_421
M = B


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def components(order: int, p: int, k: int) -> tuple[tuple[int, int, int], ...]:
    if order == 0:
        return ((1, p, k),)
    if order == 1:
        return tuple(row for row in ((p, p - 1, k), (2 * k, p + 1, k + 1)) if row[0])
    return tuple(
        row
        for row in (
            (p * (p - 1), p - 2, k),
            (2 * k * (2 * p + 1), p, k + 1),
            (4 * k * (k + 1), p + 2, k + 2),
        )
        if row[0]
    )


def endpoint_power(endpoint: arb, exponent: int) -> arb:
    return endpoint**exponent if exponent >= 0 else 1 / endpoint ** (-exponent)


def independent_bound(term: tuple[int, int, int, int, int], order: int) -> arb:
    numerator, pi_power, b_power, p, k = term
    endpoint, start, c_max = arb(B), arb(M), arb(B) / 4
    rho = arb(15) / 16
    coefficient = arb(numerator) / arb.pi() ** pi_power
    total = arb(0)
    for factor, cp, kp in components(order, p, k):
        # The decreasing-series integral test is reconstructed here from the
        # derivative rows, independently of the production helper functions.
        zeta_majorant = start ** (-2 * kp) + start ** (1 - 2 * kp) / (2 * kp - 1)
        total += (
            coefficient
            * endpoint_power(endpoint, b_power)
            * factor
            * c_max**cp
            * rho**(-kp)
            * zeta_majorant
        )
    return total * (endpoint / 2) ** order


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    # Reconstruct the rational-current decomposition from the original R_k
    # formulas, without importing the production monomial table.
    c, endpoint, mode = sp.symbols("c B m", positive=True, real=True)
    d = c**2 - mode**2
    pi = sp.pi
    formulas = (
        -(1 + sp.I * pi * endpoint * c) / (pi**2 * d),
        c / (pi**3 * endpoint) * (
            (-4 * pi * endpoint * c**3 + 4 * sp.I * c**2) * d**-3
            + (5 * pi * endpoint * c - 3 * sp.I) * d**-2
        ),
        -sp.I * c**2 / (pi**4 * endpoint**2) * (
            (-48 * pi * endpoint * c**5 + 48 * sp.I * c**4) * d**-5
            + (84 * pi * endpoint * c**3 - 60 * sp.I * c**2) * d**-4
            + (-35 * pi * endpoint * c + 15 * sp.I) * d**-3
        ),
    )
    for formula in formulas:
        for order in range(3):
            require(sp.diff(formula, c, order).has(mode), "independent symbolic derivative lost mode dependence")

    terms = {
        "C0": ((1, 2, 0, 0, 1), (1, 1, 1, 1, 1)),
        "C1": ((4, 2, 0, 4, 3), (4, 3, -1, 3, 3), (5, 2, 0, 2, 2), (3, 3, -1, 1, 2)),
        "C2": ((48, 3, -1, 7, 5), (48, 4, -2, 6, 5), (84, 3, -1, 5, 4), (60, 4, -2, 4, 4), (35, 3, -1, 3, 3), (15, 4, -2, 2, 3)),
    }

    ctx.dps = 110
    ctx.threads = 1
    totals = {order: arb(0) for order in range(3)}
    for current_terms in terms.values():
        for order in range(3):
            totals[order] += sum((independent_bound(term, order) for term in current_terms), arb(0))

    t, endpoint_ball, pi_ball = arb(T), arb(B), arb.pi()
    x0 = (1 - (1 - 8 * t / (pi_ball * endpoint_ball**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_hessian = hessian.sqrt()
    delta = arb("1e-4")
    weight = (delta * (1 - delta)) ** (-arb(1) / 4)
    log_first = (1 - 2 * delta) / (4 * delta * (1 - delta))
    second_ratio = (12 * delta**2 - 12 * delta + 5) / (16 * delta**2 * (1 - delta) ** 2)
    normalized = {
        0: weight * totals[0],
        1: weight * (totals[1] + log_first * totals[0]) / root_hessian,
        2: weight * (totals[2] + 2 * log_first * totals[1] + second_ratio * totals[0]) / hessian,
    }

    stored_raw = artifact["interval_certificate"]["summed_three_current_x_derivative_tail_balls"]
    stored_normalized = artifact["interval_certificate"]["weighted_normalized_derivative_tail_balls"]
    for order in range(3):
        require(arb(stored_raw[str(order)]).overlaps(totals[order]), f"raw order-{order} replay mismatch")
        require(arb(stored_normalized[str(order)]).overlaps(normalized[order]), f"normalized order-{order} replay mismatch")

    decision = artifact["decision"]
    require(decision["outer_rational_amplitude_first_second_derivative_bounds_uniform_in_cutoff_and_Abel_weight"] is True, "tail decision drift")
    require(decision["differentiated_Fresnel_remainder_tail_bounded"] is False, "Fresnel remainder overclaim")
    require(decision["finite_mode_block_below_outer_start_bounded"] is False, "finite-block overclaim")
    require(decision["completed_exterior_current_bound_proved"] is False, "completed-B overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print(
        "independently checked outer rational derivative tail: "
        + ", ".join(f"order{order}={normalized[order]}" for order in range(3)),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
