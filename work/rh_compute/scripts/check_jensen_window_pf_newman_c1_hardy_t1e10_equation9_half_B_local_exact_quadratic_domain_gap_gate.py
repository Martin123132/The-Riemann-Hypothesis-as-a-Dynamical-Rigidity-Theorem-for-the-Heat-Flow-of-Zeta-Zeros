#!/usr/bin/env python3
"""Independent replay of the exact local B quadratic-domain gap."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_exact_quadratic_domain_gap_gate"
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

    def coordinates(mode: int, x: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
        r = t / (2 * mp.pi * mode**2)
        v = (1 - x) / (r * x)
        s = mp.sign(v - 1) * mp.sqrt(2 * (v - 1 - mp.log(v)))
        return mp.sqrt(t) * s / 2, (endpoint * x - 2 * mode) / mp.sqrt(2 * x)

    def solve(mode: int, side: str) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
        x_m = 2 * mp.pi * mode**2 / (t + 2 * mp.pi * mode**2)
        x_cross = 2 * mp.mpf(mode) / endpoint
        if side == "left":
            f = lambda x: -2 * coordinates(mode, x)[0] + mp.pi * coordinates(mode, x)[1]
        else:
            f = lambda x: 2 * coordinates(mode, x)[0] - mp.pi * coordinates(mode, x)[1]
        root = mp.findroot(f, (x_m, x_cross))
        y, q = coordinates(mode, root)
        return root, y, q

    x621, y621, q621 = solve(621, "left")
    x622, y622, q622 = solve(622, "right")
    require(2 * abs(y621) > mp.mpf("57.17"), "mode-621 exact margin replay failed")
    require(2 * abs(y622) > mp.mpf("45.63"), "mode-622 exact margin replay failed")
    require(abs(2 * abs(y621) - mp.pi * abs(q621)) < mp.mpf("1e-80"), "mode-621 balance replay failed")
    require(abs(2 * abs(y622) - mp.pi * abs(q622)) < mp.mpf("1e-80"), "mode-622 balance replay failed")

    decision = artifact["decision"]
    require(decision["outer_and_normal_phase_is_exactly_quadratic"] is True, "quadratic-phase decision drift")
    require(decision["local_exact_max_gradient_component_above_45p63"] is True, "local-gap decision drift")
    require(decision["tangent_half_plane_degeneracy_is_exact_domain_degeneracy"] is False, "tangent guard lost")
    require(decision["local_nonstationary_integral_bound_completed"] is False, "integral overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked exact local B quadratic-domain gaps by Lambert-free replay", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
