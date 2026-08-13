#!/usr/bin/env python3
"""Independent replay of the pole-subtracted B cotangent current gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_pole_subtracted_cotangent_current_gate"
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


def exact_sum(c: mp.mpf) -> mp.mpf:
    value = mp.pi / mp.tan(mp.pi * c) / (8 * c) - 1 / (8 * c**2)
    return value - 1 / (4 * (c**2 - 621**2)) - 1 / (4 * (c**2 - 622**2))


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")

    mp.mp.dps = 100
    t = mp.mpf(T)
    endpoint = mp.mpf(B)
    x0 = (1 - mp.sqrt(1 - 8 * t / (mp.pi * endpoint**2))) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_h = mp.sqrt(hessian)
    c0 = endpoint * x0 / 2
    c_low = endpoint * (x0 - 70 / root_h) / 2
    c_high = endpoint * (x0 + 70 / root_h) / 2
    require(exact_sum(c_low) > 0 > exact_sum(c0) > exact_sum(c_high), "monotone endpoint replay failed")
    require(max(abs(exact_sum(c_low)), abs(exact_sum(c_high))) < mp.mpf("0.000288"), "uniform endpoint bound replay failed")

    root = mp.findroot(exact_sum, (mp.mpf("621.49"), mp.mpf("621.51")))
    require(mp.mpf("621.49827876663") < root < mp.mpf("621.49827876665"), "root bracket replay failed")
    xi_root = 2 * root_h / endpoint * (root - c0)
    require(mp.mpf("-6.568") < xi_root < mp.mpf("-6.566"), "standardized root replay failed")

    # Nonmatching direct series plus an explicit tail enclosure at the trace.
    cutoff = 200_000
    direct = mp.fsum(
        1 / (4 * (c0**2 - mode**2))
        for mode in range(1, cutoff + 1)
        if mode not in (621, 622)
    )
    tail_bound = 1 / (4 * cutoff * (1 - (c0 / cutoff) ** 2))
    require(abs(direct - exact_sum(c0)) < tail_bound, "direct-series cotangent replay failed")

    limit_621 = -mp.mpf(3) / (16 * 621**2) - 1 / (4 * (621**2 - 622**2))
    limit_622 = -mp.mpf(3) / (16 * 622**2) - 1 / (4 * (622**2 - 621**2))
    require(limit_621 > 0 > limit_622, "removable-limit sign replay failed")

    current_abs = 2 * mp.sqrt(4 + (mp.pi * endpoint**2 * x0) ** 2) * abs(exact_sum(c0)) / mp.pi**2
    require(current_abs > 88_000 and current_abs / root_h > mp.mpf("0.000302"), "analytic-current size replay failed")

    decision = artifact["decision"]
    require(decision["local_621_622_poles_removed_exactly"] is True, "pole-removal decision drift")
    require(decision["pole_subtracted_background_strictly_decreasing"] is True, "monotonicity decision drift")
    require(decision["raw_absolute_value_bound_is_cancellation_compatible"] is False, "triangle guard lost")
    require(decision["complete_B_face_estimate_proved"] is False, "B-face overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked pole-subtracted B cotangent current by direct-series replay", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
