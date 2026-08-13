#!/usr/bin/env python3
"""Independently check the event-ordered 398-term selector pairing."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
C = 159_577
Q = 39_894


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def original_action_defect(s: arb, mu: arb) -> arb:
    ratio = (1 - mu) / (s * s)
    q = ratio.log() / 2
    tanh_q = (ratio - 1) / (ratio + 1)
    x = (1 - tanh_q) / 2
    y_exact = ((s - 1) + tanh_q) / x
    exact = (1 - mu) * q - tanh_q - x * y_exact * y_exact / 2
    w = 1 - (2 * s + mu - 1).sqrt()
    y_canonical = w * w - mu
    canonical = w**3 / 3 - mu * w - y_canonical * y_canonical / 4
    return exact - canonical


def coefficient(mode: int, c: arb, beta: arb, mu: arb) -> acb:
    s = 4 * arb(mode) / c
    ratio = ((1 - mu) * (2 * s + mu - 1)) ** (arb(1) / 4) / s.sqrt()
    phase = beta**3 * original_action_defect(s, mu)
    return (ratio * acb(0, phase).exp() - 1) / arb(mode).sqrt()


def height_derivative(s: arb, mu: arb) -> arb:
    root = (2 * s + mu - 1).sqrt()
    numerator = 4 * s + 2 * mu - 2 * root * s.log() + root * (1 - mu).log() - 2 * root - 2
    return -numerator / (2 * root)


def panel_coefficient(mode: int, c: arb, beta: arb, center: arb, interval: arb, half_width: arb) -> acb:
    s = 4 * arb(mode) / c
    ratio = ((1 - interval) * (2 * s + interval - 1)) ** (arb(1) / 4) / s.sqrt()
    phase_center = beta**3 * original_action_defect(s, center)
    phase = arb(phase_center, beta**3 * half_width * abs(height_derivative(s, interval)).upper())
    return (ratio * acb(0, phase).exp() - 1) / arb(mode).sqrt()


def main() -> None:
    require(RESULT.is_file(), "missing pairing artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "pairing artifact did not pass")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file() and file_hash(path) == record["sha256"], f"hash drift: {path}")

    selected = {(mode, 4 * mode - C) for mode in range(39_696, 40_095)}
    remaining = selected - {(Q, -1)}
    pairs = {
        item
        for j in range(1, 199)
        for item in ((Q - j, -(4 * j + 1)), (Q + j, 4 * j - 1))
    }
    edge = {(40_093, 795), (40_094, 799)}
    require(remaining == pairs | edge and pairs.isdisjoint(edge), "independent mode pairing failed")
    require(len(pairs) == 396 and len(edge) == 2, "independent pairing count failed")

    ctx.dps = 110
    c = arb(C)
    pi = arb.pi()
    tstar = pi * c * c / 8
    beta = tstar ** (arb(1) / 3)
    radius = pi / (16 * tstar)
    maximum = arb(0)
    pair_defect = arb(0)
    sections = 512
    for section in range(sections):
        lower = radius * section / sections
        upper = radius * (section + 1) / sections
        mu = arb((lower + upper) / 2, (upper - lower) / 2)
        values = {
            mode: coefficient(mode, c, beta, mu)
            for mode, _ in remaining
        }
        maximum = max(maximum, *(abs(value).upper() for value in values.values()))
        pair_defect = max(
            pair_defect,
            *(abs(values[Q + j] - values[Q - j].conjugate()).upper() for j in range(1, 199)),
        )
    # The original-action formula has severe interval cancellation near s=1.
    # This deliberately looser cap cross-checks scale and mode coverage; the
    # sharper 7.9e-6 enclosure is certified by the stable derivative chart.
    require(maximum < arb("8.8e-6"), f"independent leading-defect cap failed: {maximum}")
    require(pair_defect < arb("1.6e-5"), f"independent conjugacy cap failed: {pair_defect}")
    aggregate_maximum = arb(0)
    aggregate_sections = 512
    for section in range(aggregate_sections):
        lower = radius * section / aggregate_sections
        upper = radius * (section + 1) / aggregate_sections
        center = (lower + upper) / 2
        half_width = (upper - lower) / 2
        interval = arb(center, half_width)
        t_interval = tstar * (1 - interval)
        aggregate = acb(0)
        for mode, _ in sorted(remaining):
            a = panel_coefficient(mode, c, beta, center, interval, half_width)
            carrier = acb(0, -t_interval * (arb(mode) / Q).log()).exp()
            aggregate += a * carrier
        aggregate_maximum = max(aggregate_maximum, abs(aggregate).upper())
    require(aggregate_maximum < arb("1.64e-5"), f"independent grouped aggregate cap failed: {aggregate_maximum}")
    saved_aggregate = arb(artifact["leading_multiplier_certificate"]["continuous_relative_carrier_aggregate_bound_ball"])
    require(saved_aggregate < arb("1.64e-5"), "saved grouped aggregate cap failed")
    require(artifact["decision"]["finite_integral_kernel_bound_proved"] is False, "finite-integral claim overpromoted")
    print("independently checked event-ordered 398-term pairing and leading multipliers", flush=True)


if __name__ == "__main__":
    main()
