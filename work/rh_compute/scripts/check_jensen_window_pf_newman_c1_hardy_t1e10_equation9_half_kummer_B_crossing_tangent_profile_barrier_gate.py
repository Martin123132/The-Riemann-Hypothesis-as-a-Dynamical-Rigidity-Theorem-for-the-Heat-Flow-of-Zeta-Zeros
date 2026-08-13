#!/usr/bin/env python3
"""Independent reverse-order replay of the B tangent crossing profile."""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_B_crossing_tangent_profile_barrier_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
CHECKER = Path(__file__).resolve()

T = 10_000_000_000
B = 5_122_421
TARGET_END = 39_894
PRECISION = 140


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def direct_tail(u: arb, imaginary_unit: acb, pi: arb) -> acb:
    argument = (-imaginary_unit * pi / 4).exp() * (pi / 2).sqrt() * u
    return (1 + imaginary_unit) / 2 * argument.erfc()


def replay_profile(mode: int, pi: arb, imaginary_unit: acb) -> tuple[acb, acb, arb, arb]:
    m = arb(mode)
    t = arb(T)
    endpoint = arb(B)
    x = 2 * pi * m**2 / (t + 2 * pi * m**2)
    q = m * (pi / (t + 2 * pi * m**2)).sqrt() * (endpoint - 2 * m - t / (pi * m))
    kappa = -t.sqrt() * (pi * endpoint * m + 2 * pi * m**2 + t) / (
        pi.sqrt() * (2 * pi * m**2 + t) ** (arb(3) / 2)
    )
    delta = 1 - pi * kappa**2 / 2
    require(not delta.contains(0), f"Delta contains zero at mode {mode}")
    abs_delta = abs(delta)
    u = abs(q) / abs_delta.sqrt()
    tail_plus = direct_tail(u, imaginary_unit, pi)
    if delta > 0:
        rotation = acb(1)
        tail = tail_plus
    else:
        rotation = imaginary_unit
        tail = tail_plus.conjugate()
    sign_q = -1 if q < 0 else 1
    scalar = -sign_q * rotation * tail / (1 + imaginary_unit)
    coefficient = x.sqrt() / (m * arb(2).sqrt())
    q_density = coefficient * rotation * (imaginary_unit * pi * q**2 / (2 * delta)).exp() / (
        abs_delta.sqrt() * imaginary_unit * pi * (1 + imaginary_unit)
    )
    return scalar, q_density, q, delta


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "B tangent artifact is not passed")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(CHECKER), "checker hash drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(dependency["sha256"] == file_hash(path), f"dependency hash drift: {path}")

    y, rho, k = sp.symbols("y rho k", real=True)
    delta = 1 - sp.pi * k**2 / 2
    original = -y**2 + sp.pi * (rho + k * y) ** 2 / 2
    completed = -delta * (y - sp.pi * k * rho / (2 * delta)) ** 2 + sp.pi * rho**2 / (2 * delta)
    require(sp.simplify(original - completed) == 0, "independent quadratic completion failed")
    require(sp.simplify(rho + k * sp.pi * k * rho / (2 * delta) - rho / delta) == 0, "independent q moment failed")

    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    imaginary_unit = acb(0, 1)
    t = arb(T)
    theta_exact = acb(arb("0.25"), t / 2).lgamma().imag - t * pi.log() / 2
    scalar_sum = acb(0)
    q_sum = acb(0)
    triangle = arb(0)
    minimum_delta = None
    minimum_mode = None
    witnesses = {row["mode"]: row for row in artifact["certificate"]["witnesses"]}

    for mode in range(TARGET_END, 0, -1):
        scalar, q_density, q, mode_delta = replay_profile(mode, pi, imaginary_unit)
        carrier = (imaginary_unit * (theta_exact - t * arb(mode).log())).exp() / arb(mode).sqrt()
        scalar_sum += carrier * scalar
        q_sum += carrier * q_density
        triangle += abs(carrier * (scalar + q_density))
        if minimum_delta is None or abs(mode_delta) < minimum_delta:
            minimum_delta = abs(mode_delta)
            minimum_mode = mode
        if mode in witnesses:
            require(q.overlaps(arb(witnesses[mode]["q_B_star_ball"])), f"q witness mismatch at mode {mode}")
            require(mode_delta.overlaps(arb(witnesses[mode]["Delta_B_ball"])), f"Delta witness mismatch at mode {mode}")
            require((scalar + q_density).overlaps(parse_complex(witnesses[mode]["grouped_profile_ball"])), f"profile witness mismatch at mode {mode}")

    total = scalar_sum + q_sum
    certificate = artifact["certificate"]
    require(total.overlaps(parse_complex(certificate["grouped_exact_theta_aggregate_ball"])), "grouped aggregate mismatch")
    require(scalar_sum.overlaps(parse_complex(certificate["scalar_only_exact_theta_aggregate_ball"])), "scalar aggregate mismatch")
    require(q_sum.overlaps(parse_complex(certificate["q_density_exact_theta_aggregate_ball"])), "q-density aggregate mismatch")
    require(triangle.overlaps(arb(certificate["termwise_complex_triangle_ball"])), "triangle mismatch")
    require(minimum_mode == certificate["minimum_integer_abs_Delta_B_mode"] == 257, "minimum Delta mode mismatch")
    require(minimum_delta.overlaps(arb(certificate["minimum_integer_abs_Delta_B_ball"])), "minimum Delta mismatch")

    decision = artifact["decision"]
    require(decision["scalar_only_counterfactual_below_normalized_target"] is True, "scalar diagnostic drift")
    require(decision["grouped_tangent_profile_below_normalized_target"] is False, "grouped target overclaim")
    require(decision["dropping_q_density_channel_admissible"] is False, "q-density omission overclaim")
    require(decision["tangent_profile_is_exact_half_Kummer_residual"] is False, "model promotion overclaim")
    require(decision["complete_T_upper_proved"] is False, "T_upper overclaim")
    print("checked B tangent profile in reverse order at 140 digits using direct erfc tails", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
