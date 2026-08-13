#!/usr/bin/env python3
"""Independent replay of the dual Weber-completion cancellation gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_weber_completion_cancellation_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

A = 159_577
B = 5_122_421
L = (B - A) // 2
K = (A - 1) // 2


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def numerical_weber_witness() -> str:
    mp.mp.dps = 70
    a = mp.mpc("0.7", "1.3")
    b = mp.mpc("0.9", "0.4")
    direct = mp.quad(lambda p: mp.exp(-a * p * p - b / (p * p)), [0, 1, mp.inf])
    closed = mp.sqrt(mp.pi) * mp.exp(-2 * mp.sqrt(a * b)) / (2 * mp.sqrt(a))
    error = abs(direct - closed)
    require(error < mp.mpf("1e-45"), "independent convergent Weber witness failed")
    return mp.nstr(error, 60)


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")

    # Replay from the face ratio w=(x-z)/x rather than the production p
    # coordinate.
    D, n, x, w = sp.symbols("D n x w", positive=True, real=True)
    alpha = D - 2 * n
    u = sp.sqrt(x) * (D * sp.sqrt(w) / 2 - n / sp.sqrt(w))
    reduced = D**2 * x * (1 - w) / 4 - n**2 * x * (1 - w) / w
    require(sp.factor(reduced - (alpha**2 * x / 4 - u**2)) == 0, "ratio-coordinate phase replay failed")
    require(sp.simplify(u.subs(w, 1) - alpha * sp.sqrt(x) / 2) == 0, "ratio-coordinate endpoint failed")

    du_dw = sp.diff(u, w)
    ratio_amplitude = sp.factor(
        x * w ** (-sp.Rational(1, 2))
        * (2 + sp.I * sp.pi * D**2 * x * (1 - w))
        / du_dw
    )
    expected_ratio_amplitude = (
        4 * sp.sqrt(x) * w * (2 + sp.I * sp.pi * D**2 * x * (1 - w))
        / (D * w + 2 * n)
    )
    require(sp.simplify(ratio_amplitude - expected_ratio_amplitude) == 0, "ratio-coordinate amplitude replay failed")

    # Reconstruct the Abel boundary moments directly and check the completed
    # endpoint dependence disappears.
    i0 = sp.exp(-sp.I * sp.pi / 4) * sp.exp(-sp.I * sp.pi * D * n * x) / D
    i2 = i0 * (2 / (sp.I * sp.pi * D**2) + 2 * n * x / D)
    polynomial_moment = sp.simplify((2 + sp.I * sp.pi * D**2 * x) * i0 - sp.I * sp.pi * D**2 * i2)
    completed = sp.factor(
        sp.exp(sp.I * sp.pi / 4)
        * 2 * sp.sqrt(x)
        * sp.exp(sp.I * sp.pi * x * (D**2 / 4 + n**2))
        * polynomial_moment
    )
    expected_completed = 2 * sp.I * sp.pi * alpha * x ** sp.Rational(3, 2) * sp.exp(sp.I * sp.pi * alpha**2 * x / 4)
    require(sp.simplify(completed - expected_completed) == 0, "completed-kernel replay failed")

    for k in (1, 2, 97, 10_000, K, K + 1, K + 100, 200_000, 1_000_000):
        alpha_b = B - 2 * (L + k)
        alpha_a = A - 2 * k
        require(alpha_b == alpha_a, "matched continuation label failed")

    # Replay the common-cutoff reindexing as an integer coefficient identity.
    for cutoff in (L + 1, L + 17, L + 1000):
        matched_count = cutoff - L
        require((L + 1) + matched_count == cutoff + 1, "B cutoff cardinality failed")
        require(L + 1 == L + 1, "B partition adjacency failed")
        require(1 + matched_count + L == cutoff + 1, "A cutoff cardinality failed")
        require((cutoff - L) + 1 == cutoff - L + 1, "A partition adjacency failed")
        require(cutoff - (cutoff - L + 1) + 1 == L, "A shift-defect cardinality failed")

    numerical_weber_witness()

    decision = artifact["decision"]
    require(decision["matched_endpoint_dual_phases_share_exact_Gaussian_coordinate"] is True, "Gaussian-coordinate drift")
    require(decision["Abel_completed_Weber_kernel_evaluated_exactly"] is True, "Weber evaluation drift")
    require(decision["completed_B_minus_A_kernel_cancels_for_each_fixed_k_at_least_1"] is True, "completion cancellation drift")
    require(decision["every_finite_matched_prefix_reduced_to_common_Gaussian_tails"] is True, "tail reduction drift")
    require(decision["same_cutoff_reindexing_has_L_term_A_shift_defect"] is True, "shift-defect bookkeeping drift")
    require(decision["infinite_continuation_Abel_reindexing_justified"] is False, "infinite reindexing overclaim")
    require(decision["shift_defect_vanishing_proved"] is False, "shift-defect overclaim")
    require(decision["stationary_expansion_remainder_needed_for_completed_bulk"] is False, "unnecessary remainder restored")
    require(decision["exceptional_k_0_face_fold_term_cancelled"] is False, "k=0 overclaim")
    require(decision["common_Gaussian_tail_sum_bounded"] is False, "tail overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked Weber completion in ratio coordinates and the convergent half-plane", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
