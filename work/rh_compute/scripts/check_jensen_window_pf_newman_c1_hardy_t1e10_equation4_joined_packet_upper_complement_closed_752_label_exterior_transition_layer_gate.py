#!/usr/bin/env python3
"""Independently check the closed 752-label exterior-transition gate."""

from __future__ import annotations

import hashlib
import json
import math
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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_label_exterior_transition_layer_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def saved_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def saved_real(record: dict[str, object]) -> arb:
    return arb(str(record["ball"]))


def q_value(y: arb, p: arb) -> arb:
    return y + p / y


def q_prime_absolute(y: arb, p: arb) -> arb:
    derivative = 1 - p / (y * y)
    require(derivative < 0, "checker endpoint escaped the lower saddle branch")
    return -derivative


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
        28,
        arb("1e-55"),
    )
    direct = acb(0)
    for index in range(count):
        direct += base_check.integrate(
            base_check.direct_label_integrand_factory(t=t, lower=lower, q=q0 + index),
            width,
            12,
            arb("1e-55"),
        )
    return grouped, direct


def main() -> int:
    resource_mode = base_check.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = 448
    ctx.threads = 1

    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production gate not passed")
    for relpath, expected in artifact["sources"].items():
        path = ROOT / relpath
        require(path.is_file(), f"missing source: {path}")
        require(file_hash(path) == expected, f"source hash drift: {path}")

    c = artifact["certificate"]
    require(c["stationary_label_count"] == 230, "stationary count drift")
    require(c["transition_label_count"] == 752, "transition count drift")
    require(c["previous_multiple_of_16_count"] == 736, "previous packet count drift")
    require(c["root_of_unity_order_at_split"] == 16, "root order drift")
    require(c["packet_factorization_exact"] == "752=47*16", "factorization drift")
    require(math.gcd(11, 16) == 1, "right endpoint is not primitive order 16")
    require(752 % 16 == 0 and 752 % 2 == 0, "dual endpoint closure drift")
    require(artifact["decision"]["finite_continuation_enclosed"] is False, "F_C overpromoted")
    require(artifact["decision"]["remote_tail_enclosed"] is False, "I_L overpromoted")

    t = arb("10000000000")
    pi = arb.pi()
    p = t / (2 * pi)
    lower = arb("621.5")
    split = arb("621.6875")
    previous_split = arb("621.625")
    upper = arb("39852.5")
    q0 = arb("2561211.5")
    q_lower = q_value(lower, p)
    require(upper * upper < p, "physical interval crosses the minimum of Q")

    finite_detuning = q0 - q_value(split, p)
    previous_finite_detuning = q0 - q_value(previous_split, p)
    finite_margin = finite_detuning**2 - 64 * q_prime_absolute(split, p)
    previous_finite_margin = (
        previous_finite_detuning**2 - 64 * q_prime_absolute(previous_split, p)
    )
    remote_margin = (q0 + 752 - q_lower) ** 2 - 64 * q_prime_absolute(lower, p)
    previous_remote_margin = (
        (q0 + 736 - q_lower) ** 2 - 64 * q_prime_absolute(lower, p)
    )
    require(previous_finite_detuning > 0, "previous grid detuning is not positive")
    require(q0 + 736 - q_lower > 0, "previous remote detuning is not positive")
    require(finite_margin > 0, "changed-precision finite margin is not positive")
    require(previous_finite_margin < 0, "changed-precision previous grid margin is not negative")
    require(remote_margin > 0, "changed-precision remote margin is not positive")
    require(previous_remote_margin < 0, "changed-precision previous remote margin is not negative")
    require(
        finite_margin.overlaps(saved_real(c["finite_continuation_eight_width_margin_ball"])),
        "finite margin does not overlap production",
    )
    require(
        previous_finite_margin.overlaps(saved_real(c["previous_grid_eight_width_margin_ball"])),
        "previous grid margin does not overlap production",
    )
    require(
        remote_margin.overlaps(saved_real(c["remote_tail_eight_width_margin_ball"])),
        "remote margin does not overlap production",
    )
    require(
        previous_remote_margin.overlaps(saved_real(c["previous_multiple_eight_width_margin_ball"])),
        "previous remote margin does not overlap production",
    )

    replay = base_check.integrate(
        base_check.sine_integrand_factory(t=t, lower=lower, q0=q0, count=752),
        split - lower,
        256,
        arb("1e-55"),
    )
    saved = saved_complex(c["grouped_752_label_transition_ball"])
    require(replay.overlaps(saved), "changed-formula production replay does not overlap")

    grouped, direct = altered_dual_endpoint_witness()
    require(grouped.overlaps(direct), "altered grouped/direct transition witness does not overlap")
    discrepancy = abs(grouped - direct)
    require(discrepancy.contains(0), "altered grouped/direct discrepancy excludes zero")

    text = " ".join(NOTE.read_text(encoding="utf-8").split())
    for fragment in (
        "752=47*16",
        "gcd(11,16)=1",
        "T_upper=G_752+F_C+I_L",
        "delta^2-64*abs(Q'(y))>0",
        "No quantitative enclosure of `F_C`",
        "No circle, polygon, fitted constant, or visual pattern supplies",
    ):
        require(fragment in text, f"missing note fragment: {fragment}")

    print(
        "independently checked closed 752-label transition layer; "
        f"altered discrepancy {discrepancy.str(10, more=True)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
