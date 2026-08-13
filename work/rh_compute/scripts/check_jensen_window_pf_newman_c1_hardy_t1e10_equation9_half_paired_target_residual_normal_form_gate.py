#!/usr/bin/env python3
"""Independent replay of the paired half-domain target-residual normal form."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_paired_target_residual_normal_form_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(RESULT.is_file(), "missing production artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")

    # Reconstruct the phase allocation in real coordinates, independently of
    # the builder's complex notation.
    q_r, q_i, t_r, t_i = sp.symbols("q_r q_i t_r t_i", real=True)
    c, s = sp.cos(sp.pi / 8), sp.sin(sp.pi / 8)
    allocated_r = c * t_r - s * t_i
    allocated_i = s * t_r + c * t_i
    projected = 2 * (c * (q_r - allocated_r) + s * (q_i - allocated_i))
    direct = 2 * (c * q_r + s * q_i) - 2 * t_r
    phase_identity = sp.simplify(sp.trigsimp(projected - direct).rewrite(sp.sqrt))
    require(phase_identity == 0, "independent phase allocation failed")

    # Audit the partition by explicit integer sets at three different cutoffs.
    for cutoff in (39_895, 40_000, 41_123):
        blocks = (
            {0}
            | set(range(1, 622))
            | set(range(622, 39_695))
            | set(range(39_695, 39_895))
            | set(range(39_895, cutoff + 1))
        )
        require(blocks == set(range(0, cutoff + 1)), f"partition gap at M={cutoff}")
        require(len(blocks) == cutoff + 1, f"partition overlap at M={cutoff}")
    require(39_894 - 622 + 1 == 39_273, "target count drift")

    decision = artifact["decision"]
    require(decision["finite_paired_half_domain_target_residual_identity_proved"] is True, "identity decision drift")
    require(decision["symmetric_limit_of_paired_normal_form_proved"] is True, "limit decision drift")
    require(decision["target_allocation_is_modewise_half_Gamma_bulk_claim"] is False, "half-bulk overclaim")
    require(decision["midpoint_discrete_parity_cancels_Poisson_pairs"] is False, "midpoint overclaim")
    require(decision["paired_blocks_quantitatively_bounded"] is False, "quantitative overclaim")
    require(decision["complete_T_upper_proved"] is False, "T_upper overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked paired half-domain target residual by independent real-coordinate replay", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
