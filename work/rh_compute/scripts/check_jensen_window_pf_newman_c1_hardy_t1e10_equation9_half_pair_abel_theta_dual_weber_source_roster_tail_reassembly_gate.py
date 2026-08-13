#!/usr/bin/env python3
"""Independent replay of Weber source-roster and theta-tail reassembly."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_weber_source_roster_tail_reassembly_gate"
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

    # Check the reindexing on a finite symbolic model with distinct atoms.
    small_l = 4
    cutoff = 13
    b_atoms = [f"B{n}" for n in range(1, cutoff + 1)]
    a_atoms = [f"A{n}" for n in range(1, cutoff - small_l + 1)]
    roster = [f"B{n}" for n in range(1, small_l + 1)]
    continuation_b = [f"B{small_l + k}" for k in range(1, cutoff - small_l + 1)]
    require(roster + continuation_b == b_atoms, "finite B reindexing replay failed")
    require([f"A{k}" for k in range(1, cutoff - small_l + 1)] == a_atoms, "finite A reindexing replay failed")

    # Replay source normalization starting from W_t rather than a_t.
    x, alpha = sp.symbols("x alpha", positive=True, real=True)
    source_weight = (x * (1 - x)) ** (-sp.Rational(1, 4))
    completed = 2 * sp.I * sp.pi * alpha * x ** sp.Rational(3, 2)
    pair_prefactor = sp.Rational(1, 2) / (sp.I * sp.pi)
    triangle_to_source = sp.simplify(pair_prefactor * x ** (-sp.Rational(7, 4)) * (1 - x) ** (-sp.Rational(1, 4)) * completed)
    require(sp.simplify(triangle_to_source - alpha * source_weight) == 0, "independent source normalization failed")

    labels = range(B - 2, A - 1, -2)
    require(len(labels) == L, "completed roster count failed")
    require(labels[0] == B - 2 and labels[-1] == A, "completed roster endpoints failed")
    require(B not in labels, "upper label incorrectly included")

    # Derive the tail in v=x/p^2=1/s coordinates, unlike production.
    D, n, v = sp.symbols("D n v", positive=True, real=True)
    p2 = x / v
    tail_phase = D**2 * x / 4 + n**2 * x - D**2 * p2 / 4 - n**2 * x**2 / p2
    expected = x * (D**2 * (1 - 1 / v) / 4 + n**2 * (1 - v))
    require(sp.simplify(tail_phase - expected) == 0, "reciprocal-tail phase failed")

    g_a = sp.exp(sp.I * sp.pi * A**2 * x * (1 - 1 / v) / 4) * (2 + sp.I * sp.pi * A**2 * x * (1 - 1 / v))
    g_b = sp.exp(sp.I * sp.pi * B**2 * x * (1 - 1 / v) / 4) * (2 + sp.I * sp.pi * B**2 * x * (1 - 1 / v))
    difference = sp.simplify(g_a - g_b)
    require(sp.simplify(difference.subs(v, 1)) == 0, "tail face cancellation replay failed")
    derivative_s = -sp.diff(difference, v).subs(v, 1)
    require(sp.simplify(derivative_s - 3 * sp.I * sp.pi * (B**2 - A**2) * x / 2) == 0, "tail face slope replay failed")

    decision = artifact["decision"]
    require(decision["positive_dual_sector_reindexed_with_cutoff_guard"] is True, "reindexing decision drift")
    require(decision["completed_finite_B_block_has_exact_equation9_source_normalization"] is True, "source normalization drift")
    require(decision["completed_B_block_recovers_labels_B_minus_2_through_A"] is True, "roster recovery drift")
    require(decision["upper_source_label_B_recovered_by_this_gate"] is False, "upper-label overclaim")
    require(decision["all_remaining_positive_dual_tails_reassembled_at_same_index"] is True, "same-index tail drift")
    require(decision["summed_tail_driver_vanishes_at_s_equal_1"] is True, "face cancellation drift")
    require(decision["summed_theta_tail_quantitatively_bounded"] is False, "tail overclaim")
    require(decision["zero_mode_direct_exceptional_sector_closed"] is False, "exceptional-sector overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked source normalization and theta-tail reassembly in reciprocal-tail coordinates", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
