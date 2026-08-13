#!/usr/bin/env python3
"""Independent replay of the bare full-saddle interface barrier."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from flint import acb, arb, ctx


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_full_saddle_interface_barrier_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fresnel(z: arb, pi: arb, imaginary_unit: acb) -> acb:
    return (imaginary_unit * pi / 4).exp() / arb(2).sqrt() * (((-imaginary_unit * pi / 4).exp() * (pi / 2).sqrt() * z).erf())


def replay(mode: int) -> tuple[arb, arb, acb]:
    pi = arb.pi()
    imaginary_unit = acb(0, 1)
    m = arb(mode)
    t = arb(10_000_000_000)
    x = 2 * pi * m**2 / (t + 2 * pi * m**2)
    q_a = (x / 2).sqrt() * (arb(159_577) - 2 * m / x)
    q_b = (x / 2).sqrt() * (arb(5_122_421) - 2 * m / x)
    i0 = fresnel(q_b, pi, imaginary_unit) - fresnel(q_a, pi, imaginary_unit)
    i1 = ((imaginary_unit * pi * q_b**2 / 2).exp() - (imaginary_unit * pi * q_a**2 / 2).exp()) / (imaginary_unit * pi)
    ratio = i0 / (1 + imaginary_unit) + x.sqrt() * i1 / (m * arb(2).sqrt() * (1 + imaginary_unit))
    normalized = abs(ratio - 1) / m.sqrt()
    physical = arb(2).sqrt() / pi * normalized
    return normalized, physical, ratio


def main() -> int:
    ctx.dps = 125
    ctx.threads = 1
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "barrier artifact is not passed")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(CHECKER), "checker hash drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(dependency["sha256"] == file_hash(path), f"dependency hash drift: {path}")

    rows = {row["mode"]: row for row in artifact["certificate"]["rows"]}
    require(set(rows) == {39694, 39695, 39852, 39853}, "witness roster drift")
    for mode, row in rows.items():
        normalized, physical, ratio = replay(mode)
        require(normalized.overlaps(arb(row["carrier_scaled_normalized_defect_ball"])), f"normalized replay mismatch at {mode}")
        require(physical.overlaps(arb(row["carrier_scaled_physical_defect_ball"])), f"physical replay mismatch at {mode}")
        stored_ratio = acb(arb(row["incomplete_grouped_current_ratio_ball"]["real_ball"]), arb(row["incomplete_grouped_current_ratio_ball"]["imag_ball"]))
        require(ratio.real.overlaps(stored_ratio.real) and ratio.imag.overlaps(stored_ratio.imag), f"ratio replay mismatch at {mode}")

    normalized_39695, _, _ = replay(39695)
    normalized_39852, physical_39852, _ = replay(39852)
    require(normalized_39695 > 55 * arb("0.000019"), "atlas-edge barrier replay failed")
    require(normalized_39852 > 131 * arb("0.000019"), "structural normalized barrier replay failed")
    require(physical_39852 > 130 * arb("0.0000086"), "structural physical barrier replay failed")
    require(artifact["decision"]["complete_x_integrated_mode_error_lower_bound_proved"] is False, "full-mode lower-bound overclaim")
    require(artifact["decision"]["complete_T_upper_proved"] is False, "T_upper overclaim")

    print("checked bare full-saddle interface barrier at four modes with 125-digit replay", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
