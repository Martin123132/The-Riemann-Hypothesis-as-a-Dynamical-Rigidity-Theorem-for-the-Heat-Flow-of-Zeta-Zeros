#!/usr/bin/env python3
"""Independently check the grouped upper-complement endpoint-saddle core."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Callable


ROOT = Path(__file__).resolve().parents[3]
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_endpoint_saddle_core_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()
COMPLEMENT = ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_ordinary_complementary_half_lattice_gate.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def saved_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


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


def sine_integrand_factory(
    *, t: arb, lower: arb, q0: arb, count: int
) -> Callable[[acb, bool], acb]:
    pi = arb.pi()
    imaginary = acb(0, 1)
    endpoint_phase = -t * lower.log() - pi * lower * lower + 2 * pi * q0 * lower
    endpoint_carrier = lower ** arb("-0.5") * (imaginary * endpoint_phase).exp()
    linear = 2 * pi * (q0 - lower) - t / lower

    def integrand(x: acb, analytic: bool) -> acb:
        w = x / lower
        root = (1 + w).sqrt(analytic=analytic)
        relative_phase = linear * x - t * ((1 + w).log() - w) - pi * x * x
        y = lower + x
        geometric = (
            (imaginary * pi * (count - 1) * y).exp()
            * (pi * count * y).sin()
            / (pi * y).sin()
        )
        return endpoint_carrier / root * (imaginary * relative_phase).exp() * geometric

    return integrand


def direct_label_integrand_factory(
    *, t: arb, lower: arb, q: arb
) -> Callable[[acb, bool], acb]:
    pi = arb.pi()
    imaginary = acb(0, 1)
    endpoint_phase = -t * lower.log() - pi * lower * lower + 2 * pi * q * lower
    endpoint_carrier = lower ** arb("-0.5") * (imaginary * endpoint_phase).exp()
    linear = 2 * pi * (q - lower) - t / lower

    def integrand(x: acb, analytic: bool) -> acb:
        w = x / lower
        root = (1 + w).sqrt(analytic=analytic)
        relative_phase = linear * x - t * ((1 + w).log() - w) - pi * x * x
        return endpoint_carrier / root * (imaginary * relative_phase).exp()

    return integrand


def integrate(
    integrand: Callable[[acb, bool], acb], width: arb, panels: int, tolerance: arb
) -> acb:
    total = acb(0)
    for index in range(panels):
        left = width * index / panels
        right = width * (index + 1) / panels
        total += acb.integral(
            integrand,
            left,
            right,
            abs_tol=tolerance,
            rel_tol=tolerance,
            eval_limit=350_000,
            depth_limit=42,
        )
    require(total.is_finite(), "independent integral is not finite")
    return total


def altered_witness() -> tuple[acb, acb]:
    pi = arb.pi()
    p = arb("4.5725")
    t = 2 * pi * p
    lower = arb("1.5")
    split = arb("1.625")
    q0 = arb("4.5")
    count = 4
    width = split - lower
    grouped = integrate(
        sine_integrand_factory(t=t, lower=lower, q0=q0, count=count),
        width,
        24,
        arb("1e-55"),
    )
    direct = acb(0)
    for index in range(count):
        direct += integrate(
            direct_label_integrand_factory(t=t, lower=lower, q=q0 + index),
            width,
            12,
            arb("1e-55"),
        )
    return grouped, direct


def main() -> int:
    resource_mode = set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = 448
    ctx.threads = 1
    for path in (RESULT, NOTE, BUILDER, CHECKER, COMPLEMENT):
        require(path.is_file(), f"missing artifact: {path}")

    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production gate not passed")
    for path in (COMPLEMENT, BUILDER, CHECKER):
        expected = artifact["sources"][relative(path)]
        require(expected == file_hash(path), f"source hash drift: {path}")

    c = artifact["certificate"]
    require(c["label_count"] == 230, "label count drift")
    require(c["rational_split"].startswith("621.5625"), "rational split drift")
    require(c["endpoint_dirichlet_value_exact"] == "sum_(j=0)^229 (-1)^j=0", "endpoint parity drift")
    require(artifact["decision"]["remaining_upper_piece_enclosed"] is False, "upper remainder overpromoted")

    t = arb("10000000000")
    lower = arb("621.5")
    split = arb("621.5625")
    q0 = arb("2561211.5")
    replay = integrate(
        sine_integrand_factory(t=t, lower=lower, q0=q0, count=230),
        split - lower,
        48,
        arb("1e-55"),
    )
    saved = saved_complex(c["grouped_endpoint_saddle_core_ball"])
    require(replay.overlaps(saved), "changed-formula production replay does not overlap")

    grouped, direct = altered_witness()
    require(grouped.overlaps(direct), "altered grouped/direct roster witness does not overlap")
    discrepancy = abs(grouped - direct)
    require(discrepancy.contains(0), "altered grouped/direct discrepancy excludes zero")

    required_note_fragments = (
        "sum_(j=0)^229 (-1)^j=0",
        "T_upper=P_core+R_upper",
        "No quantitative enclosure of `R_upper`",
        "no fitted geometric constant",
    )
    text = NOTE.read_text(encoding="utf-8")
    for fragment in required_note_fragments:
        require(fragment in text, f"missing note fragment: {fragment}")

    print(
        "independently checked grouped 230-label endpoint-saddle core; "
        f"altered discrepancy {discrepancy.str(10, more=True)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
