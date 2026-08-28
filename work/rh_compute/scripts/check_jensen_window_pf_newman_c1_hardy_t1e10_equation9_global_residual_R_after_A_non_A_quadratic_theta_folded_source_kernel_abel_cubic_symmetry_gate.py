#!/usr/bin/env python3
"""Independently replay the folded-source derivative and cubic symmetry gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
import flint
from flint import acb, arb


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_folded_source_kernel_abel_cubic_symmetry_gate"
)
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
CHECKER = Path(__file__).resolve()
WALL_RESULT = REPO_ROOT / (
    "work/rh_compute/results/"
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_small_tau_grouped_first_wall_bridge_gate.json"
)

A = 159_577
B = 5_122_421
L = 2_481_422
FIRST = 622
LAST = 39_936
Q = LAST - FIRST + 1
R = Q // 3
N = 496_283
T_WALL = Fraction(1, 2 * N)
S0 = Fraction(1, 3)
PRECISION_BITS = 416
BLOCK_SIZE = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rat(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def cis_pi(value: Fraction) -> acb:
    return (acb.pi() * acb(0, 1) * acb(rat(value % 2))).exp()


def parse_ball(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def x_of_t(side: str, t: Fraction) -> Fraction:
    coefficient = 4 if side == "right" else 6
    return 2 * (1 + t) / (5 + coefficient * t)


def dx_dt(side: str, t: Fraction) -> Fraction:
    coefficient = 4 if side == "right" else 6
    sign = 1 if side == "right" else -1
    return sign * 2 / (5 + coefficient * t) ** 2


def reverse_block_source_moments(x: Fraction, s: Fraction) -> tuple[acb, acb, acb, acb]:
    totals = [acb(0), acb(0), acb(0), acb(0)]
    starts = list(range(0, L, BLOCK_SIZE))
    step = cis_pi(2 * x)
    for start in reversed(starts):
        stop = min(L, start + BLOCK_SIZE)
        alpha = Fraction(A) + 2 * s + 2 * start
        term = cis_pi(x * alpha * alpha / 4)
        ratio = cis_pi(x * (alpha + 1))
        alpha_ball = rat(alpha)
        local = [acb(0), acb(0), acb(0), acb(0)]
        for _ in range(start, stop):
            alpha2 = alpha_ball * alpha_ball
            local[0] += term
            local[1] += acb(alpha_ball) * term
            local[2] += acb(alpha2) * term
            local[3] += acb(alpha2 * alpha_ball) * term
            term *= ratio
            ratio *= step
            alpha_ball += 2
        for degree in reversed(range(4)):
            totals[degree] += local[degree]
    return tuple(totals)  # type: ignore[return-value]


def direct_root_sums() -> tuple[acb, acb]:
    root_sum = acb(0)
    weighted_sum = acb(0)
    omega = cis_pi(Fraction(-2, 3))
    term = omega ** FIRST
    for m in range(FIRST, LAST + 1):
        root_sum += term
        weighted_sum += m * term
        term *= omega
    return root_sum, weighted_sum


def verify_dependency_hashes(artifact: dict[str, Any]) -> None:
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency {dependency['path']}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {dependency['path']}")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(CHECKER), "checker hash drift")


def main() -> int:
    require(RESULT.is_file() and BUILDER.is_file() and WALL_RESULT.is_file(), "missing replay input")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    wall = json.loads(WALL_RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact did not pass")
    verify_dependency_hashes(artifact)

    require(B - A == 2 * L, "source roster identity failed")
    require(Q == 39_315 and Q == 3 * R and R == 13_105, "cubic roster identity failed")
    require(FIRST % 3 == 1 and LAST % 3 == 0, "cubic endpoint residues failed")
    require(A % 2 == 1, "Fourier unfolding parity failed")

    flint.ctx.prec = PRECISION_BITS
    sqrt_three = arb(3).sqrt()
    exact_weighted = arb(R) * acb(arb(3) / 2, sqrt_three / 2)
    c_prime = 2 * acb.pi() * acb(0, 1) * exact_weighted
    recorded_c_prime = parse_ball(artifact["kernel_symmetry_certificate"]["first_derivative_ball"])
    require(c_prime.overlaps(recorded_c_prime), "kernel derivative ball mismatch")
    c_prime_abs = 2 * arb.pi() * R * sqrt_three
    require(
        (c_prime.real**2 + c_prime.imag**2).overlaps(c_prime_abs**2),
        "kernel derivative magnitude mismatch",
    )

    root_sum, weighted_sum = direct_root_sums()
    require(root_sum.contains(acb(0)), "direct cubic root sum misses zero")
    require(weighted_sum.overlaps(exact_weighted), "direct weighted root sum mismatch")

    # Independently check the coefficient reduction: y=2(phi'_m(y)+m)/x.
    x_symbol = Fraction(7, 13)
    y_symbol = Fraction(17, 5)
    m_symbol = 622
    phi_prime = x_symbol * y_symbol / 2 - m_symbol
    require(y_symbol == 2 * (phi_prime + m_symbol) / x_symbol, "Fresnel coefficient algebra failed")

    require(len(artifact["wall_source_derivative_rows"]) == 2, "wall row count drift")
    for row in artifact["wall_source_derivative_rows"]:
        side = row["side"]
        x = x_of_t(side, T_WALL)
        require(str(x) == row["x"], f"{side} x map mismatch")
        m0, m1, m2, m3 = reverse_block_source_moments(x, S0)
        source = m1
        source_x = acb.pi() * acb(0, 1) * m3 / 4
        source_s = 2 * m0 + acb.pi() * acb(0, 1) * acb(rat(x)) * m2
        source_t = acb(rat(dx_dt(side, T_WALL))) * source_x
        product_s = source * c_prime
        mixed_ts = source_t * c_prime

        require(source.overlaps(parse_ball(row["direct_source_F"])), f"{side} source replay mismatch")
        require(source_x.overlaps(parse_ball(row["direct_partial_x_F"])), f"{side} x derivative mismatch")
        require(source_t.overlaps(parse_ball(row["direct_partial_t_F"])), f"{side} t derivative mismatch")
        require(source_s.overlaps(parse_ball(row["direct_partial_s_F"])), f"{side} s derivative mismatch")
        require(
            product_s.overlaps(parse_ball(row["interior_limit_partial_s_product_at_s0"])),
            f"{side} product s derivative mismatch",
        )
        require(
            mixed_ts.overlaps(parse_ball(row["interior_limit_mixed_partial_t_partial_s_product_at_s0"])),
            f"{side} mixed product derivative mismatch",
        )
        grouped = parse_ball(wall[f"{side}_first_wall_source_box"]["complete_source_box"])
        require(source.overlaps(grouped), f"{side} source misses grouped dependency")
        require(parse_ball(row["interior_limit_product_F_times_C_at_s0"]).contains(acb(0)), "product zero drift")
        require(
            parse_ball(row["interior_limit_partial_t_product_at_s0"]).contains(acb(0)),
            "product t zero drift",
        )

    decision = artifact["decision"]
    require(decision["pointwise_interior_limit_interchange_across_integer_boundary_forbidden"] is True, "boundary guard lost")
    require(decision["complete_cell_Abel_limit_reduced_to_endpoint_half_sum_and_finite_modes"] is True, "Abel handoff lost")
    require(decision["near_constant_source_product_cover_selected"] is False, "value-box route reintroduced")
    require(decision["physical_quadrature_completed"] is False, "physical quadrature overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently replayed folded-source wall derivatives and cubic kernel symmetry", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
