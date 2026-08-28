#!/usr/bin/env python3
"""Independently check the quarter-disk vertical/arc/tail gate."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
from typing import Callable


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_quarter_disk_vertical_arc_tail_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
ARB_HELPER = Path(__file__).with_name(STEM.removesuffix("_gate") + "_arb.py")


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


def finite_difference(u: mp.mpf | mp.mpc, q: mp.mpf, count: int) -> mp.mpc:
    if u == 0:
        return mp.mpc(count)
    return mp.exp(-q * u) * mp.expm1(-count * u) / mp.expm1(-u)


def direct_roster(y: mp.mpf, q: mp.mpf, count: int) -> mp.mpc:
    return mp.fsum(mp.exp(2j * mp.pi * (q + index) * y) for index in range(count))


def prefactor(s: mp.mpc) -> mp.mpc:
    return mp.exp(mp.pi * mp.im(s) / 2 + 0.25j * mp.pi) * mp.power(2 * mp.pi, s - 1)


def zero_to_bound(
    s: mp.mpc,
    rest: Callable[[mp.mpf], mp.mpc],
    bound: mp.mpf,
    direct_points: tuple[mp.mpf, ...],
) -> mp.mpc:
    split = min(bound, mp.mpf(1))
    start = -mp.log(split)

    def logarithmic(v: mp.mpf) -> mp.mpc:
        if not mp.isfinite(v):
            return mp.mpc(0)
        y = mp.exp(-v)
        return mp.exp(-(1 - s) * v) * rest(y)

    cuts = (start, start + 1, start + 2.5, start + 6, start + 13, start + 26)
    value = mp.quad(logarithmic, cuts)
    value += mp.quadosc(logarithmic, [cuts[-1], mp.inf], omega=abs(mp.im(s)))
    if bound > 1:
        require(direct_points[0] == 1 and direct_points[-1] == bound, "bad direct partition")
        value += mp.quad(lambda y: mp.power(y, -s) * rest(y), direct_points)
    return value


def independent_contour_replay() -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    mp.mp.dps = 78
    s = mp.mpc("0.5", "4.75")
    q = mp.mpf("3.5")
    count = 5
    y_end = mp.mpf("4.5")
    radius = 2 * mp.pi * y_end
    e_s = prefactor(s)
    positive_points = (
        mp.mpf(1),
        mp.mpf("1.75"),
        mp.mpf(3),
        mp.mpf(6),
        mp.mpf(10),
        mp.mpf(16),
        mp.mpf(22),
        radius,
    )
    finite = zero_to_bound(
        s,
        lambda u: mp.exp(1j * u * u / (4 * mp.pi)) * finite_difference(u, q, count),
        radius,
        positive_points,
    )
    tail = mp.quad(
        lambda u: mp.power(u, -s)
        * mp.exp(1j * u * u / (4 * mp.pi))
        * finite_difference(u, q, count),
        (radius, mp.mpf(34), mp.mpf(42), mp.mpf(55), mp.mpf(72), mp.mpf(96)),
    )
    tail_cut = mp.mpf(96)
    omitted_tail_bound = (
        abs(e_s)
        * tail_cut ** mp.mpf("-0.5")
        * mp.exp(-q * tail_cut)
        / (q * (1 - mp.exp(-tail_cut)))
    )
    vertical = zero_to_bound(
        s,
        lambda y: mp.exp(-1j * mp.pi * y * y) * direct_roster(y, q, count),
        y_end,
        (mp.mpf(1), mp.mpf("1.75"), mp.mpf(3), y_end),
    )

    u0 = -1j * radius
    endpoint = (
        e_s
        * mp.power(u0, -s)
        * mp.exp(1j * u0 * u0 / (4 * mp.pi))
        * finite_difference(u0, q, count)
        * radius
    )

    def normalized(delta: mp.mpf) -> mp.mpc:
        w_exponent = 1j * radius * mp.expm1(1j * delta)
        w = mp.exp(w_exponent)
        exponent = (
            mp.im(s) * delta
            + 0.5j * delta
            - 1j * mp.pi * y_end * y_end * mp.expm1(2j * delta)
            + 1j * q * radius * mp.expm1(1j * delta)
        )
        return mp.exp(exponent) * (1 + w**count) / (1 + w)

    arc = endpoint * mp.quad(
        normalized,
        (
            mp.mpf(0),
            mp.mpf("1e-9"),
            mp.mpf("1e-7"),
            mp.mpf("1e-5"),
            mp.mpf("0.001"),
            mp.mpf("0.02"),
            mp.mpf("0.2"),
            mp.mpf("0.7"),
            mp.pi / 2,
        ),
    )
    finite_error = abs(e_s * finite - vertical - arc)
    full_error = abs(e_s * (finite + tail) - vertical - arc - e_s * tail)
    scale_error = abs(abs(endpoint) - mp.sqrt(y_end))
    require(finite_error < mp.mpf("1e-38"), "independent finite contour replay failed")
    require(full_error < mp.mpf("1e-38"), "independent full contour replay failed")
    require(scale_error < mp.mpf("1e-62"), "independent endpoint scale failed")
    require(omitted_tail_bound < mp.mpf("1e-100"), "independent omitted tail is too large")
    return finite_error, full_error, scale_error, omitted_tail_bound


def main() -> int:
    priority = set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")
    for path in (RESULT, NOTE, BUILDER, ARB_HELPER):
        require(path.is_file(), f"missing artifact: {path.name}")

    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "saved gate did not pass")
    require(
        artifact.get("status")
        == "exact_quarter_disk_identity_and_production_upper_arc_arb_enclosure_certified_lower_cell_and_complete_DK_open",
        "status drift",
    )
    require(artifact["source_hashes"]["builder"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["source_hashes"]["checker"] == file_hash(Path(__file__)), "checker hash drift")
    require(artifact["source_hashes"]["arb_helper"] == file_hash(ARB_HELPER), "Arb helper hash drift")
    for row in artifact["dependencies"].values():
        path = REPO_ROOT / row["path"]
        require(path.is_file() and file_hash(path) == row["sha256"], "dependency hash drift")

    finite_error, full_error, scale_error, omitted_tail = independent_contour_replay()
    command = [
        "py",
        "-3.13",
        str(ARB_HELPER),
        "--variant",
        "independent",
        "--reference",
        str(RESULT),
    ]
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        timeout=1800,
        check=False,
    )
    require(completed.returncode == 0, f"independent Arb helper failed: {completed.stderr[-2000:]}")
    independent = json.loads(completed.stdout)
    require(independent.get("passed") is True, "independent Arb replay did not pass")
    require(independent.get("overlaps_reference") is True, "independent Arb balls do not overlap")
    require(independent["degree"] != artifact["arb_certificate"]["degree"], "degree was not changed")
    require(
        independent["precision_bits"] != artifact["arb_certificate"]["precision_bits"],
        "precision was not changed",
    )
    require(
        independent["delta_min"] != artifact["arb_certificate"]["delta_min"],
        "initial split was not changed",
    )
    require(
        independent["rotated_arc_absolute"]["upper_float"] < 0.057665,
        "independent arc scale guard failed",
    )
    require("No lower-cell enclosure" in artifact["scope"]["not_certified"], "proof boundary weakened")
    require("arc is not D_K" in artifact["route_decision"]["comparison_only"], "D_K guard missing")

    print(
        "independently checked quarter-disk vertical arc tail representation; "
        f"altered errors {mp.nstr(finite_error, 6)}, {mp.nstr(full_error, 6)}, "
        f"{mp.nstr(scale_error, 6)}; omitted tail {mp.nstr(omitted_tail, 6)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
