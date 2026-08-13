#!/usr/bin/env python3
"""Independent replay of the discrete B-face joint-gradient gap."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_discrete_joint_gradient_gap_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

T = 10_000_000_000
B = 5_122_421


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")

    mp.mp.dps = 100
    t = mp.mpf(T)
    endpoint = mp.mpf(B)
    x0 = (1 - mp.sqrt(1 - 8 * t / (mp.pi * endpoint**2))) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_h = mp.sqrt(hessian)

    def g(x: mp.mpf) -> mp.mpf:
        return (mp.pi * endpoint**2 / 4 - t / (2 * x * (1 - x))) / root_h

    def q(mode: int, x: mp.mpf) -> mp.mpf:
        return (endpoint * x - 2 * mode) / mp.sqrt(2 * x)

    left = mp.findroot(lambda x: -g(x) - q(621, x), (2 * mp.mpf(621) / endpoint, x0))
    right = mp.findroot(lambda x: g(x) + q(622, x), (x0, 2 * mp.mpf(622) / endpoint))
    left_margin = -g(left)
    right_margin = g(right)
    require(left_margin > mp.mpf("28.08"), "mode-621 minimax replay failed")
    require(right_margin > mp.mpf("22.40"), "mode-622 minimax replay failed")
    require(right_margin < left_margin, "adjacent margin ordering replay failed")

    # Directly minimize a dense set of neighboring mode profiles.  The
    # analytic monotonic argument in the production gate promotes this local
    # check to every positive integer.
    for mode in range(615, 629):
        crossing = 2 * mp.mpf(mode) / endpoint
        if mode <= 621:
            balance = mp.findroot(lambda x: -g(x) - q(mode, x), (crossing, x0))
            margin = -g(balance)
        else:
            balance = mp.findroot(lambda x: g(x) + q(mode, x), (x0, crossing))
            margin = g(balance)
        require(margin >= right_margin, f"neighbor mode {mode} beats mode 622")

    require(abs(root_h * (left - x0) + mp.mpf("28.0699707266557")) < mp.mpf("1e-12"), "left xi replay drift")
    require(abs(root_h * (right - x0) - mp.mpf("22.4197118656154")) < mp.mpf("1e-12"), "right xi replay drift")

    decision = artifact["decision"]
    require(decision["all_positive_integer_face_minimax_margin_above_22p40"] is True, "global gap decision drift")
    require(decision["smooth_two_direction_partition_denominator_above_501"] is True, "partition decision drift")
    require(decision["full_triangle_has_no_interior_saddle"] is False, "interior-saddle guard lost")
    require(decision["two_direction_IBP_constant_completed"] is False, "IBP overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked discrete B-face joint-gradient gap by minimax replay", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
