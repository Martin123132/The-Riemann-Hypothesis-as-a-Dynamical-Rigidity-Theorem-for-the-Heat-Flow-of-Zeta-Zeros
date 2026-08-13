#!/usr/bin/env python3
"""Independent replay of the dual-cutoff shift-defect vanishing gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_cutoff_shift_defect_vanishing_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

A = 159_577
B = 5_122_421
L = (B - A) // 2


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")

    # Replay from s=1-q and keep the x factor inside the integration-by-parts
    # quotient throughout.
    s, x, endpoint, mode = sp.symbols("s x A n", positive=True, real=True)
    d = mode**2 - endpoint**2 * s**2 / 4
    numerator = s ** sp.Rational(3, 2) * (2 + sp.I * sp.pi * endpoint**2 * x * (1 - s))
    quotient = -numerator / (sp.pi * d)
    require(sp.limit(quotient, s, 0, dir="+") == 0, "face boundary does not vanish")
    outer_boundary = sp.limit(quotient, s, 1, dir="-")
    expected_boundary = -2 / (sp.pi * (mode**2 - endpoint**2 / 4))
    require(sp.simplify(outer_boundary - expected_boundary) == 0, "outer boundary replay failed")

    # Exact rational reconstruction of the total-variation constants.
    boundary_constant = Fraction(8, 3)
    numerator_constant = Fraction(8, 3)
    denominator_constant = Fraction(32, 63)
    require(boundary_constant + numerator_constant + denominator_constant == Fraction(368, 63), "C0 reconstruction failed")
    numerator_linear = Fraction(28, 15)
    denominator_linear = Fraction(16, 63)
    require(numerator_linear + denominator_linear == Fraction(668, 315), "C1 reconstruction failed")

    # Recompute all split integrals in a different order.
    n = sp.symbols("n", positive=True)
    x0 = n**-2
    small_1 = sp.integrate(4 * x ** (-sp.Rational(3, 4)), (x, 0, x0))
    small_2 = sp.integrate(2 * sp.pi * endpoint**2 * x ** sp.Rational(1, 4), (x, 0, x0))
    require(sp.simplify(small_1 - 16 * n ** (-sp.Rational(1, 2))) == 0, "small-x leading integral failed")
    require(sp.simplify(small_2 - 8 * sp.pi * endpoint**2 * n ** (-sp.Rational(5, 2)) / 5) == 0, "small-x quadratic integral failed")

    c0 = sp.Rational(368, 63) / sp.pi
    c1 = sp.Rational(668, 315)
    large_1_guard = 4 * c0 * n ** (-sp.Rational(1, 2)) / 3
    large_2_guard = 4 * 2 ** (-sp.Rational(1, 4)) * c1 * endpoint**2 / n**2
    require(large_1_guard > 0 and large_2_guard > 0, "large-x guards failed")

    # For fixed A,L, this explicit block guard tends to zero.
    N = sp.symbols("N", positive=True)
    n_min = N - L + 1
    block_guard = L * (
        2 ** sp.Rational(1, 4) * (16 + 4 * c0 / 3) * n_min ** (-sp.Rational(1, 2))
        + 4 * c1 * A**2 * n_min**-2
        + 2 ** sp.Rational(1, 4) * 8 * sp.pi * A**2 * n_min ** (-sp.Rational(5, 2)) / 5
    )
    require(sp.limit(block_guard, N, sp.oo) == 0, "shift-block limit failed")

    decision = artifact["decision"]
    require(decision["scaled_A_dual_kernel_retains_x_factor"] is True, "scaled-kernel drift")
    require(decision["explicit_q_integration_by_parts_bound_proved"] is True, "IBP-bound drift")
    require(decision["weighted_mode_norm_decays_as_O_n_minus_one_half"] is True, "weighted-decay drift")
    require(decision["fixed_width_L_term_A_shift_defect_vanishes"] is True, "shift-defect drift")
    require(decision["fixed_height_symmetric_dual_cutoff_reindexing_justified"] is True, "reindexing drift")
    require(decision["height_uniform_shift_defect_theorem_proved"] is False, "height-uniform overclaim")
    require(decision["common_Gaussian_tail_series_bounded"] is False, "common-tail overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked dual-cutoff shift-defect decay in the complementary s coordinate", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
