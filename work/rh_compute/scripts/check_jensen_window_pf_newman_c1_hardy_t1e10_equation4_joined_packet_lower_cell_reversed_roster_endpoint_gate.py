#!/usr/bin/env python3
"""Independently replay the reversed-roster lower-cell endpoint gate."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_lower_cell_reversed_roster_endpoint_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = ROOT / "work" / "rh_compute" / "scripts" / f"{STEM}.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            process.cpu_affinity([process.cpu_affinity()[0]])
            return "below_normal_one_cpu"
        process.nice(10)
        process.cpu_affinity([process.cpu_affinity()[0]])
        return "nice_10_one_cpu"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def load_builder():
    spec = importlib.util.spec_from_file_location("lower_cell_endpoint_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load builder module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def altered_direct_integral() -> tuple[mp.mpc, mp.mpf]:
    """Directly integrate a small changed finite roster in the logarithmic chart."""

    mp.mp.dps = 90
    t = mp.mpf(50)
    boundary = mp.mpf("1.5")
    q_plus = mp.mpf("5.5")
    count = 3
    labels = [q_plus - r for r in range(1, count + 1)]

    def integrand(v: mp.mpf) -> mp.mpc:
        if not mp.isfinite(v):
            return mp.mpc(0)
        y = boundary * mp.exp(-v)
        roster = mp.fsum(
            mp.exp(-1j * mp.pi * y**2 + 2j * mp.pi * label * y)
            for label in labels
        )
        return y ** (mp.mpf("0.5") - 1j * t) * roster

    first = mp.quad(integrand, [0, 1, 2, 4, 8, 16, 32, 64, mp.inf])
    mp.mp.dps = 110
    second = mp.quad(integrand, [0, mp.mpf("0.5"), 1, 3, 6, 12, 24, 48, 96, mp.inf])
    return second, abs(second - first)


def main() -> int:
    resource_mode = set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True and NOTE.is_file(), "saved gate or note missing")
    require(
        artifact["decision"]["lower_cell_complex_ball_certified"] is True,
        "lower-cell decision drift",
    )
    for dependency in artifact["dependencies"].values():
        path = ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], "dependency hash drift")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")

    b = load_builder()
    b.ctx.prec = 448
    b.ctx.threads = 1
    replay = b.lower_cell_certificate(
        b.arb(b.HEIGHT),
        b.arb_from_fraction(b.LOWER_BOUNDARY),
        b.arb_from_fraction(b.Q_PLUS),
        b.SOURCE_COUNT,
        20,
        b.arb(280),
        8192,
        resource_mode,
    )
    saved = b.complex_from_record(artifact["certificate"]["lower_cell_integral_ball"])
    fresh = b.complex_from_record(replay["lower_cell_integral_ball"])
    require(fresh.overlaps(saved), "altered-order production replay misses saved lower-cell ball")
    require(
        b.arb(replay["total_remainder_bound_ball"]["ball"]).upper()
        < b.arb(artifact["certificate"]["total_remainder_bound_ball"]["ball"]).upper(),
        "higher-order replay did not tighten the lower-cell remainder",
    )

    # The changed finite roster catches the half-integer parity and recurrence signs.
    altered = b.lower_cell_certificate(
        b.arb(50),
        b.arb("1.5"),
        b.arb("5.5"),
        3,
        4,
        b.arb("0.5"),
        2048,
        resource_mode,
    )
    direct, direct_drift = altered_direct_integral()
    expanded = b.complex_from_record(altered["lower_cell_integral_ball"])
    midpoint = mp.mpc(float(expanded.real.mid()), float(expanded.imag.mid()))
    allowed = mp.mpf(str(float(b.arb(altered["total_remainder_bound_ball"]["ball"]).upper())))
    require(abs(direct - midpoint) < allowed, "altered direct roster misses endpoint enclosure")
    require(direct_drift < mp.mpf("1e-6"), "altered direct quadrature is unstable")

    print(
        "independently checked reversed-roster lower cell; "
        f"production remainder {replay['total_remainder_bound_ball']['ball']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
