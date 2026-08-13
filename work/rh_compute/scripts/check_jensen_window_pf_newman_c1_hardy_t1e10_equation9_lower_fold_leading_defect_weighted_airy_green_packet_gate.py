#!/usr/bin/env python3
"""Independently check the leading-defect weighted Airy packet gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_leading_defect_weighted_airy_green_packet_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
C = 159_577
Q = (C - 1) // 4


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ball_upper(text: str) -> mp.mpf:
    if "+/-" not in text:
        return mp.mpf(text.strip("[] "))
    midpoint, radius = text.strip("[] ").split("+/-")
    return mp.mpf(midpoint.strip()) + mp.mpf(radius.strip())


def action_defect(s: mp.mpf, mu: mp.mpf) -> mp.mpf:
    root = mp.sqrt(mu + 2 * s - 1)
    return (
        -2 * mu * root / 3
        + mu * mp.log(s)
        - mu * mp.log(1 - mu) / 2
        + 3 * mu / 2
        + s**2 / 2
        - 4 * s * root / 3
        + 2 * s
        + 2 * root / 3
        - mp.log(s)
        + mp.log(1 - mu) / 2
        - mp.mpf(11) / 6
    )


def defect(mode: int, mu: mp.mpf, beta_cubed: mp.mpf) -> mp.mpc:
    s = 4 * mp.mpf(mode) / C
    ratio = ((1 - mu) * (2 * s + mu - 1)) ** mp.mpf("0.25") / mp.sqrt(s)
    return (ratio * mp.exp(1j * beta_cubed * action_defect(s, mu)) - 1) / mp.sqrt(mode)


def effective_variation(modes: list[int], values: dict[int, mp.mpc]) -> mp.mpf:
    best = mp.mpf("0")
    for start in range(len(modes)):
        block_maximum = mp.mpf("0")
        total_variation = mp.mpf("0")
        for stop in range(start, len(modes)):
            if stop > start:
                total_variation += abs(values[modes[stop]] - values[modes[stop - 1]])
            block_maximum = max(block_maximum, abs(values[modes[stop]]))
            coefficient_variation = min(abs(values[modes[start]]), abs(values[modes[stop]])) + total_variation
            best = max(best, coefficient_variation + block_maximum)
    return best


def main() -> None:
    require(RESULT.is_file(), "missing result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "stored gate is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")

    mp.mp.dps = 90
    beta_cubed = mp.pi * C**2 / 8
    mu_max = mp.pi / (16 * beta_cubed)
    minus_modes = list(range(Q - 198, Q))
    plus_modes = list(range(Q + 1, Q + 201))
    certificate = artifact["certificate"]
    stored_minus = ball_upper(certificate["minus_branch"]["maximum_any_contiguous_modulus_adjusted_variation_ball"])
    stored_plus = ball_upper(certificate["plus_branch"]["maximum_any_contiguous_modulus_adjusted_variation_ball"])
    for mu in (mp.mpf("0"), mu_max / 2, mu_max):
        values = {mode: defect(mode, mu, beta_cubed) for mode in minus_modes + plus_modes}
        require(effective_variation(minus_modes, values) < stored_minus, "minus sampled variation escaped interval bound")
        require(effective_variation(plus_modes, values) < stored_plus, "plus sampled variation escaped interval bound")

    weighted_k = ball_upper(certificate["weighted_two_branch_K_packet_bound_ball"])
    weighted_derivative = ball_upper(certificate["weighted_two_branch_partial_s_K_packet_bound_ball"])
    require(weighted_k < mp.mpf("7.58e-7"), "weighted K threshold failed")
    require(weighted_derivative < mp.mpf("0.001634"), "weighted derivative threshold failed")
    require(artifact["decision"]["finite_integral_amplitude_remainder_proved"] is False, "proof-boundary drift")
    print("independently checked leading-defect weighted Airy Green packets", flush=True)


if __name__ == "__main__":
    main()
