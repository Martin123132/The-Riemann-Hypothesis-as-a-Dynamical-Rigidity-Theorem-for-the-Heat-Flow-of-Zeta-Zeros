#!/usr/bin/env python3
"""Independently check the gamma-balanced Mellin-Hurwitz packet gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
from typing import Callable


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_gamma_balanced_mellin_hurwitz_finite_difference_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
ENTIRE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_riemann_auxiliary_entire_kernel_finite_difference_gate.json"
HALF_LINE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_QK_geometric_half_line_and_saddle_route_gate.json"
PI_DIGITS = "3141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067"
PI_DENOMINATOR = 10 ** (len(PI_DIGITS) - 1)
PI_LOWER = Fraction(int(PI_DIGITS), PI_DENOMINATOR)
PI_UPPER = Fraction(int(PI_DIGITS) + 1, PI_DENOMINATOR)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            process.cpu_affinity([process.cpu_affinity()[0]])
            return "below_normal_one_cpu"
        process.nice(10)
        process.cpu_affinity([process.cpu_affinity()[0]])
        return "nice_10_one_cpu"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def delta(u: mp.mpf | mp.mpc, q: mp.mpf, count: int) -> mp.mpc:
    if u == 0:
        return mp.mpc(count)
    return mp.exp(-q * u) * mp.expm1(-count * u) / mp.expm1(-u)


def integrate(
    s: mp.mpc,
    rest: Callable[[mp.mpf], mp.mpc],
    low_cuts: tuple[mp.mpf, ...],
    high_cuts: tuple[mp.mpf, ...],
) -> mp.mpc:
    def low(v: mp.mpf) -> mp.mpc:
        if not mp.isfinite(v):
            return mp.mpc(0)
        u = mp.exp(-v)
        return mp.exp(-(1 - s) * v) * rest(u)

    require(low_cuts[-1] == mp.inf, "logarithmic cuts must end at infinity")
    require(high_cuts[-1] != mp.inf, "upper cutoff must be finite and explicit")
    low_value = mp.quad(low, low_cuts[:-1])
    low_value += mp.quadosc(low, [low_cuts[-2], mp.inf], omega=abs(mp.im(s)))
    return low_value + mp.quad(lambda u: mp.power(u, -s) * rest(u), high_cuts)


def prefactor(s: mp.mpc) -> mp.mpc:
    return mp.exp(mp.pi * mp.im(s) / 2 + 0.25j * mp.pi) * mp.power(2 * mp.pi, s - 1)


def main() -> int:
    priority = set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")
    for path in (RESULT, NOTE, BUILDER, ENTIRE, HALF_LINE):
        require(path.is_file(), f"missing gate artifact: {path.name}")

    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "saved gate is not passed")
    require(artifact["source_hashes"]["builder"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["source_hashes"]["checker"] == file_hash(Path(__file__)), "checker hash drift")
    require(
        artifact["dependencies"]["entire_kernel"]["sha256"] == file_hash(ENTIRE),
        "entire-kernel dependency drift",
    )
    require(
        artifact["dependencies"]["geometric_half_line"]["sha256"] == file_hash(HALF_LINE),
        "half-line dependency drift",
    )

    mp.mp.dps = 105
    s = mp.mpc("0.5", "4.75")
    first_label = 7
    count = 4
    q = mp.mpf(first_label) / 2
    q_plus = q + count
    low_cuts = (mp.mpf(0), mp.mpf("0.75"), mp.mpf(2), mp.mpf(5), mp.mpf(11), mp.mpf(20), mp.inf)
    high_cuts = (
        mp.mpf(1),
        mp.mpf("1.5"),
        mp.mpf(3),
        mp.mpf(6),
        mp.mpf(10),
        mp.mpf(20),
        mp.mpf(40),
        mp.mpf(64),
    )
    omega = mp.exp(0.25j * mp.pi)
    k0 = mp.exp(3 * mp.pi * mp.im(s) / 4 + 3j * mp.pi / 8)

    direct_integral = integrate(
        s,
        lambda x: mp.exp(-mp.pi * x * x)
        * mp.fsum(
            mp.exp(-mp.pi * omega * (first_label + 2 * index) * x)
            for index in range(count)
        ),
        low_cuts,
        high_cuts,
    )
    transformed_integral = integrate(
        s,
        lambda u: mp.exp(1j * u * u / (4 * mp.pi)) * delta(u, q, count),
        low_cuts,
        high_cuts,
    )
    require(
        abs(k0 * direct_integral - prefactor(s) * transformed_integral) < mp.mpf("1e-66"),
        "changed-row scaled Mellin identity failed",
    )

    anchor = integrate(s, lambda u: delta(u, q, count), low_cuts, high_cuts)
    anchor_exact = mp.gamma(1 - s) * (mp.zeta(1 - s, q) - mp.zeta(1 - s, q_plus))
    require(abs(anchor - anchor_exact) < mp.mpf("1e-66"), "changed-row Hurwitz anchor failed")

    order = 4
    expansion = mp.fsum(
        mp.power(1j / (4 * mp.pi), index)
        / mp.factorial(index)
        * mp.gamma(1 - s + 2 * index)
        * (mp.zeta(1 - s + 2 * index, q) - mp.zeta(1 - s + 2 * index, q_plus))
        for index in range(order)
    )

    def remainder(u: mp.mpf) -> mp.mpc:
        x = u * u / (4 * mp.pi)
        polynomial = mp.fsum(mp.power(1j * x, index) / mp.factorial(index) for index in range(order))
        return delta(u, q, count) * (mp.exp(1j * x) - polynomial)

    remainder_value = integrate(s, remainder, low_cuts, high_cuts)
    require(
        abs(transformed_integral - expansion - remainder_value) < mp.mpf("1e-65"),
        "changed-row finite expansion failed",
    )
    explicit_bound = (
        abs(prefactor(s))
        * mp.gamma(2 * order + mp.mpf("0.5"))
        * (mp.zeta(2 * order + mp.mpf("0.5"), q) - mp.zeta(2 * order + mp.mpf("0.5"), q_plus))
        / (mp.power(4 * mp.pi, order) * mp.factorial(order))
    )
    require(abs(prefactor(s) * remainder_value) <= explicit_bound, "changed-row remainder bound failed")

    balance = abs(prefactor(s) * mp.gamma(1 - s)) ** 2
    require(
        abs(balance - 1 / (1 + mp.exp(-2 * mp.pi * mp.im(s)))) < mp.mpf("1e-88"),
        "changed-row gamma balance failed",
    )

    threshold = artifact["actual_height_scale_audit"][
        "last_order_with_certified_ratio_upper_bound_below_one"
    ]
    require(threshold == 20_000_022_018, "saved threshold drift")
    q_actual = Fraction(159577, 2)
    coordinate = lambda item: Fraction(item, 1) + Fraction(3, 16 * (item + 1))
    require(coordinate(threshold) < PI_LOWER * q_actual * q_actual, "threshold lower proof failed")
    require(
        coordinate(threshold + 1) > PI_UPPER * q_actual * q_actual,
        "threshold upper proof failed",
    )
    require(len(artifact["actual_height_scale_audit"]["bound_rows"]) == 5, "scale rows drift")
    require(artifact["decision"]["absolute_taylor_remainder_route_rejected"] is True, "route decision drift")
    require(artifact["decision"]["J_Z_enclosed"] is False, "J_Z overclaim")
    require(artifact["decision"]["D_K_enclosed"] is False, "D_K overclaim")
    print("independently checked gamma-balanced Mellin-Hurwitz finite difference and absolute Taylor barrier")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
