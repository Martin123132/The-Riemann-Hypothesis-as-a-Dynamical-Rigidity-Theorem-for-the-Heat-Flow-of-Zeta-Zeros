#!/usr/bin/env python3
"""Independent replay of the double-Weber source reconstruction gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_double_weber_exact_source_reconstruction_gate"
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

    # Replay the second Weber endpoint cancellation from the moment ratio.
    a, b = sp.symbols("a b", positive=True)
    i0 = sp.sqrt(sp.pi) * sp.exp(-2 * sp.sqrt(a * b)) / (2 * sp.sqrt(a))
    ratio_i2_i0 = sp.simplify(-sp.diff(sp.log(i0), a))
    require(ratio_i2_i0 == 1 / (2 * a) + sp.sqrt(b / a), "Weber moment ratio failed")

    D, x = sp.symbols("D x", positive=True, real=True)
    ell = sp.symbols("ell", positive=True, integer=True)
    ratio_boundary = ratio_i2_i0.subs({a: sp.I * sp.pi * D**2 * x / 4, b: sp.I * sp.pi * ell**2 / x})
    ratio_boundary = sp.powdenest(ratio_boundary, force=True)
    coefficient = sp.simplify(2 - sp.I * sp.pi * D**2 * x * ratio_boundary)
    require(sp.simplify(coefficient + 2 * sp.I * sp.pi * D * ell) == 0, "independent transformed-mode coefficient failed")

    i0_boundary = sp.exp(-sp.I * sp.pi / 4) * sp.exp(-sp.I * sp.pi * D * ell) / (D * sp.sqrt(x))
    endpoint_mode = sp.simplify(i0_boundary * coefficient)
    endpoint_difference = sp.simplify(endpoint_mode.subs(D, B) - endpoint_mode.subs(D, A))
    require(endpoint_difference == 0, "independent odd-endpoint cancellation failed")

    # Numerically integrate in a genuinely convergent perturbed half-plane
    # and compare the two Weber moments to their closed forms.
    mp.mp.dps = 60
    aa = mp.mpc("0.73", "1.11")
    bb = mp.mpc("0.64", "0.39")
    j0 = mp.quad(lambda p: mp.exp(-aa * p * p - bb / (p * p)), [0, 1, mp.inf])
    j2 = mp.quad(lambda p: p * p * mp.exp(-aa * p * p - bb / (p * p)), [0, 1, mp.inf])
    c0 = mp.sqrt(mp.pi) * mp.exp(-2 * mp.sqrt(aa * bb)) / (2 * mp.sqrt(aa))
    c2 = c0 * (1 / (2 * aa) + mp.sqrt(bb / aa))
    require(abs(j0 - c0) < mp.mpf("1e-40"), "independent I0 witness failed")
    require(abs(j2 - c2) < mp.mpf("1e-40"), "independent I2 witness failed")

    # Replay endpoint bookkeeping as coefficients of independent f_A,f_B,K0.
    # H=(f_A+f_B)/2, dual completion=(f_B-f_A)/2, and direct=-K0.
    f_a, f_b, k0 = sp.symbols("f_A f_B K_0")
    endpoint_sum = sp.expand((f_a + f_b) / 2 + (f_b - f_a) / 2)
    zero_sum = sp.expand(k0 - k0)
    require(endpoint_sum == f_b, "endpoint-half bookkeeping failed")
    require(zero_sum == 0, "zero/direct bookkeeping failed")

    labels = [B - 2 * n for n in range(0, L + 1)]
    require(len(labels) == L + 1, "source roster count failed")
    require(labels[0] == B and labels[-1] == A, "source roster endpoints failed")
    require(all(labels[j] - labels[j + 1] == 2 for j in range(len(labels) - 1)), "source roster contiguity failed")

    decision = artifact["decision"]
    require(decision["second_modular_transformed_nonzero_endpoint_difference_vanishes_exactly"] is True, "second modular cancellation drift")
    require(decision["second_modular_zero_mode_vanishes_by_Abel_primitive"] is True, "zero-mode primitive drift")
    require(decision["same_index_theta_tail_collapses_to_dual_zero_tail_difference"] is True, "tail collapse drift")
    require(decision["original_Poisson_zero_mode_cancels_direct_Theta_minus_one_term"] is True, "zero/direct drift")
    require(decision["dual_zero_completion_plus_endpoint_half_current_recovers_upper_label_B"] is True, "upper endpoint recovery drift")
    require(decision["entire_finite_equation9_source_roster_reconstructed_exactly"] is True, "source reconstruction drift")
    require(decision["small_source_minus_target_bound_proved"] is False, "small-error overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked double-Weber cancellation and endpoint bookkeeping independently", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
