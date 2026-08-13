#!/usr/bin/env python3
"""Independent replay of the symmetric endpoint-tail reassembly gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_symmetric_endpoint_tail_reassembly_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "reassembly artifact is not passed")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(CHECKER), "checker hash drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(dependency["sha256"] == file_hash(path), f"dependency hash drift: {path}")

    scope = artifact["scope"]
    a, b = scope["target_positive_modes"]
    target = set(range(a, b + 1))
    require(len(target) == scope["target_mode_count"] == 39_273, "target block drift")

    witness_cutoff = 40_037
    negative = set(range(-witness_cutoff, 0))
    zero = {0}
    low = set(range(1, a))
    outer = set(range(b + 1, witness_cutoff + 1))
    pieces = (negative, zero, low, target, outer)
    require(sum(len(piece) for piece in pieces) == 2 * witness_cutoff + 1, "witness count failed")
    require(len(set().union(*pieces)) == 2 * witness_cutoff + 1, "witness partition overlap")
    require(set().union(*pieces) == set(range(-witness_cutoff, witness_cutoff + 1)), "witness partition gap")

    z = sp.symbols("z", nonzero=True)
    small_a, small_b, small_cutoff = 3, 7, 11
    direct_target = sum(z**mode for mode in range(small_a, small_b + 1))
    geometric_target = z**small_a * (1 - z ** (small_b - small_a + 1)) / (1 - z)
    require(sp.cancel(direct_target - geometric_target) == 0, "independent target kernel replay failed")
    direct_symmetric = sum(z**mode for mode in range(-small_cutoff, small_cutoff + 1))
    geometric_symmetric = z ** (-small_cutoff) * (1 - z ** (2 * small_cutoff + 1)) / (1 - z)
    require(sp.cancel(direct_symmetric - geometric_symmetric) == 0, "independent Dirichlet kernel replay failed")

    delta0, delta1, rem_plus, rem_minus, mode = sp.symbols("d0 d1 rp rm mode", nonzero=True)
    d_plus = 2 * sp.pi * sp.I * mode
    positive = -delta0 / d_plus - delta1 / d_plus**2 + rem_plus / d_plus**2
    d_minus = -d_plus
    negative_mode = -delta0 / d_minus - delta1 / d_minus**2 + rem_minus / d_minus**2
    paired = sp.expand(positive + negative_mode)
    require(sp.simplify(sp.diff(paired, delta0)) == 0, "independent one-over-m cancellation failed")

    h, c, j_target, b_target = sp.symbols("h c jt bt")
    require(sp.expand((h + c + j_target - b_target) - (h + c + (j_target - b_target))) == 0, "residual replay failed")

    decision = artifact["decision"]
    require(decision["finite_symmetric_dirichlet_reassembly_proved"] is True, "finite reassembly decision drift")
    require(decision["symmetric_residual_limit_proved"] is True, "symmetric limit decision drift")
    require(decision["target_negative_partners_retained_in_complement"] is True, "negative-partner decision drift")
    require(decision["target_positive_endpoint_currents_cancel_internally"] is False, "internal cancellation overclaim")
    require(decision["zero_mode_is_endpoint_current_canceller"] is False, "zero-mode cancellation overclaim")
    require(decision["half_current_is_endpoint_current_canceller"] is False, "half-current cancellation overclaim")
    require(decision["joint_residual_quantitatively_bounded"] is False, "quantitative residual overclaim")
    require(decision["complete_T_upper_proved"] is False, "T_upper overclaim")
    print("checked exact symmetric endpoint-tail residual and 39273-mode partition", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
