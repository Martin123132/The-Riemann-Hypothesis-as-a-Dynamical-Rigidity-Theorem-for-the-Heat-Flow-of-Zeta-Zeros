#!/usr/bin/env python3
"""Independently check the selected/logistic stationary-action bridge."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_logistic_stationary_action_bridge_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
PRECISION = 120


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
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def f_center(s: arb) -> arb:
    root = (2 * s - 1).sqrt()
    return (3 * s * s - 8 * s * root + 12 * s + 4 * root - 6 * s.log() - 11) / 6


def f_mu(s: arb, mu: arb) -> arb:
    root = (2 * s + mu - 1).sqrt()
    return -(4 * s + 2 * mu - 2 * root * s.log() + root * (1 - mu).log() - 2 * root - 2) / (2 * root)


def main() -> None:
    priority = set_low_priority()
    require(RESULT.is_file(), "missing stationary-action bridge artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "stationary-action bridge did not pass")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file() and file_hash(path) == record["sha256"], f"hash drift: {path}")

    s, mu = sp.symbols("s mu", positive=True)
    root = sp.sqrt(2 * s + mu - 1)
    formula = (3 * s**2 - 8 * s * root + 12 * s - 4 * mu * root + 6 * mu * sp.log(s) - 3 * mu * sp.log(1 - mu) + 9 * mu + 4 * root - 6 * sp.log(s) + 3 * sp.log(1 - mu) - 11) / 6
    center = sp.simplify(formula.subs(mu, 0))
    require([sp.simplify(sp.diff(center, s, k).subs(s, 1)) for k in range(6)] == [0, 0, 0, 0, 0, 6], "independent fifth-order tangency failed")
    height = sp.simplify(sp.diff(formula, mu).subs(mu, 0))
    require([sp.simplify(sp.diff(height, s, k).subs(s, 1)) for k in range(4)] == [0, 0, 0, -1], "independent height tangency failed")
    amplitude = ((1 - mu) * (2 * s + mu - 1)) ** sp.Rational(1, 4) / sp.sqrt(s)
    amplitude_center = sp.simplify(amplitude.subs(mu, 0))
    require(
        [sp.simplify(sp.diff(amplitude_center, s, k).subs(s, 1)) for k in range(4)]
        == [1, 0, -sp.Rational(1, 2), 3],
        "independent amplitude tangency failed",
    )

    ctx.dps = PRECISION
    pi = arb.pi()
    c = arb(159_577)
    tstar = pi * c * c / 8
    beta = tstar ** (arb(1) / 3)
    mu_radius = pi / (16 * tstar)
    mu_ball = arb(mu_radius / 2, mu_radius / 2)
    maximum = arb(0)
    maximum_mode = None
    maximum_carrier_scaled = arb(0)
    maximum_carrier_mode = None
    for mode in range(39_696, 40_095):
        scale = 4 * arb(mode) / c
        center_action = beta**3 * f_center(scale)
        variation = beta**3 * mu_radius * abs(f_mu(scale, mu_ball))
        bound = abs(center_action) + variation
        amplitude_ratio = ((1 - mu_ball) * (2 * scale + mu_ball - 1)) ** (arb(1) / 4) / scale.sqrt()
        combined = abs(amplitude_ratio - 1) + amplitude_ratio.upper() * bound
        carrier_scaled = combined / arb(mode).sqrt()
        if bound.upper() > maximum:
            maximum = bound.upper()
            maximum_mode = mode
        if carrier_scaled.upper() > maximum_carrier_scaled:
            maximum_carrier_scaled = carrier_scaled.upper()
            maximum_carrier_mode = mode
        saved = next(row for row in artifact["certificate"]["rows"] if row["mode"] == mode)
        require(center_action.overlaps(arb(saved["center_action_defect_ball"])), f"mode {mode} center action misses saved ball")
        delta = scale - 1
        exact_y = beta * beta * (delta * delta - mu_ball) / scale
        canonical_y = 2 * beta * beta * (scale - (2 * scale + mu_ball - 1).sqrt())
        require(exact_y.lower() > 0 and canonical_y.lower() > 0, f"mode {mode} saddle positivity failed")
    require(maximum_mode == 40_094 and maximum < arb("0.001555"), "independent maximum action bound failed")
    require(maximum_carrier_mode == 40_094 and maximum_carrier_scaled < arb("7.8e-6"), "independent carrier-scaled bound failed")
    saved_carrier = arb(artifact["certificate"]["maximum_carrier_scaled_combined_leading_defect_bound_ball"])
    require(maximum_carrier_scaled <= saved_carrier.upper(), "independent carrier bound exceeds saved enclosure")
    require(arb(artifact["certificate"]["remaining_local_normalized_headroom_lower_ball"]) > arb("9e-6"), "local headroom guard failed")
    require(artifact["decision"]["selected_branch_amplitude_match_proved"] is False, "amplitude match overpromoted")
    require(artifact["decision"]["phase_multiplier_must_remain_inside_grouped_projection"] is True, "grouping guard lost")
    print(f"independently validated stationary-action bridge on 399 modes; priority={priority}", flush=True)


if __name__ == "__main__":
    main()
