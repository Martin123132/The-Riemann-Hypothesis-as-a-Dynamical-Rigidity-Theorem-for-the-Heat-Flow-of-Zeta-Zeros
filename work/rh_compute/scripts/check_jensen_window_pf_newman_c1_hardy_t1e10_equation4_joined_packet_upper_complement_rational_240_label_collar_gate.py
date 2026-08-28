#!/usr/bin/env python3
"""Independently check the rationally closed 240-label upper collar."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx

import check_jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_endpoint_saddle_core_gate as base_check


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_rational_240_label_collar_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()
BASE_RESULT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_endpoint_saddle_core_gate.json"
BASE_BUILDER = ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_endpoint_saddle_core_gate.py"
BASE_CHECKER = ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_endpoint_saddle_core_gate.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def altered_dual_endpoint_witness() -> tuple[acb, acb]:
    pi = arb.pi()
    p = arb("4.5725")
    t = 2 * pi * p
    lower = arb("1.5")
    split = arb("1.625")
    q0 = arb("4.5")
    count = 8
    width = split - lower
    grouped = base_check.integrate(
        base_check.sine_integrand_factory(t=t, lower=lower, q0=q0, count=count),
        width,
        24,
        arb("1e-55"),
    )
    direct = acb(0)
    for index in range(count):
        direct += base_check.integrate(
            base_check.direct_label_integrand_factory(t=t, lower=lower, q=q0 + index),
            width,
            10,
            arb("1e-55"),
        )
    return grouped, direct


def main() -> int:
    resource_mode = base_check.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = 448
    ctx.threads = 1
    for path in (RESULT, NOTE, BUILDER, CHECKER, BASE_RESULT, BASE_BUILDER, BASE_CHECKER):
        require(path.is_file(), f"missing artifact: {path}")

    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production gate not passed")
    for path in (BASE_RESULT, BASE_BUILDER, BASE_CHECKER, BUILDER, CHECKER):
        require(artifact["sources"][relative(path)] == file_hash(path), f"source hash drift: {path}")

    c = artifact["certificate"]
    require(c["stationary_label_count"] == 230, "stationary count drift")
    require(c["collar_label_count"] == 240, "collar count drift")
    require(c["root_of_unity_order_at_split"] == 16, "root order drift")
    require(c["least_closed_collar_count"] == 240, "least collar drift")
    require(artifact["decision"]["upper_remainders_enclosed"] is False, "remainders overpromoted")

    replay = base_check.integrate(
        base_check.sine_integrand_factory(
            t=arb("10000000000"),
            lower=arb("621.5"),
            q0=arb("2561211.5"),
            count=240,
        ),
        arb("0.0625"),
        48,
        arb("1e-55"),
    )
    saved = base_check.saved_complex(c["grouped_240_label_collar_ball"])
    require(replay.overlaps(saved), "changed-formula production replay does not overlap")

    grouped, direct = altered_dual_endpoint_witness()
    require(grouped.overlaps(direct), "altered grouped/direct collar witness does not overlap")
    discrepancy = abs(grouped - direct)
    require(discrepancy.contains(0), "altered collar discrepancy excludes zero")

    text = " ".join(NOTE.read_text(encoding="utf-8").split())
    for fragment in (
        "240=15*16",
        "T_upper=P_240+R_B+R_L",
        "No quantitative enclosure of `R_B`, `R_L`",
        "no fitted geometric constant",
    ):
        require(fragment in text, f"missing note fragment: {fragment}")

    print(
        "independently checked rational 240-label collar; "
        f"altered discrepancy {discrepancy.str(10, more=True)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
