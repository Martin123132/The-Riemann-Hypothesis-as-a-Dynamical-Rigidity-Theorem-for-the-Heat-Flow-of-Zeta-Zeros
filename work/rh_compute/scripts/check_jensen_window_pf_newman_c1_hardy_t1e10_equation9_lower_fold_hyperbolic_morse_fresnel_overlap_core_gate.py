#!/usr/bin/env python3
"""Independently validate the hyperbolic Morse-Fresnel overlap core."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, REPO_ROOT / "work/rh_compute/scripts"):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_overlap_core_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_overlap_core_gate.md"
PRECISION = 110
C = 159_577
MODE = 39_894
DELTA_RADIUS = arb("0.00284")
PANELS = 98_304
SERIES_DEGREE = 20


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


def ball(left: arb, right: arb) -> arb:
    return arb((left + right) / 2, (right - left) / 2)


def uabs(value: arb) -> arb:
    return abs(value).upper()


def k_value(delta: arb) -> arb:
    value = arb(0)
    for power in range(SERIES_DEGREE, -1, -1):
        value = value * delta + arb(1 if power % 2 == 0 else -1) / (power + 2)
    tail = DELTA_RADIUS ** (SERIES_DEGREE + 1) / ((SERIES_DEGREE + 3) * (1 - DELTA_RADIUS))
    return value + arb(0, tail)


def row(t: arb, delta: arb, pi: arb, tau: arb, beta: arb, ymax: arb) -> dict[str, arb]:
    m = arb(MODE)
    r = t / (2 * pi * m**2)
    d0 = 1 + r
    den = d0 + r * delta
    x = 1 / den
    k = k_value(delta)
    require(k.lower() > 0, "independent k cover lost positivity")
    u = delta * (t * k).sqrt()
    a0 = (tau - t) / (pi * m)
    b = t / (pi * m)
    qc = (a0 - b * delta) * (pi / (2 * den)).sqrt()
    inside = -pi * d0 * a0 * b - pi * a0**2 * r / 2 + t * d0 * delta * (r - k * den)
    chi = delta * inside / (2 * d0 * den)
    ell = 2 * (t - tau) / (C * pi.sqrt())
    rb = chi - ell * u - u**2 / (2 * C) - u**3 / (3 * C * pi.sqrt())
    rho = (2 * x).sqrt()
    cross = (
        pi.sqrt() * a0 / den
        + delta * t * (k - 2 * r / den**2)
        / ((t * k).sqrt() + b * pi.sqrt() / den)
    )
    quad = x - arb(1) / 2
    g = (1 + delta) * (2 * k).sqrt()
    base = r ** (arb(1) / 4) * (1 + delta) ** (-arb(1) / 4) * g * rho
    eps = (r * x / t).sqrt()
    return {
        "u": u,
        "rb": rb,
        "cross": cross,
        "quad": quad,
        "amp0": base * (1 + qc * eps),
        "amp1": base * (1 + (qc + rho * ymax) * eps),
    }


def phase_bound(r: dict[str, arb], ymax: arb) -> arb:
    endpoints = [uabs(r["rb"]), uabs(r["rb"] + r["cross"] * ymax + r["quad"] * ymax**2)]
    result = max(endpoints)
    if not r["quad"].contains(0):
        yc = -r["cross"] / (2 * r["quad"])
        if yc.lower() >= 0 and yc.upper() <= ymax.lower():
            result = max(result, uabs(r["rb"] - r["cross"] ** 2 / (4 * r["quad"])))
        elif not (yc.upper() < 0 or yc.lower() > ymax.upper()):
            result = max(result, uabs(r["rb"]) + ymax * uabs(r["cross"]) + ymax**2 * uabs(r["quad"]))
    elif r["cross"].contains(0):
        result = max(result, uabs(r["rb"]) + ymax * uabs(r["cross"]) + ymax**2 * uabs(r["quad"]))
    return result


def scan(t: arb, pi: arb, tau: arb, beta: arb, ymax: arb) -> dict[str, arb]:
    maxima = {name: arb(0) for name in ("boundary_phase_remainder", "linear_inner_remainder", "quadratic_inner_remainder", "complete_phase_remainder", "normalized_amplitude_remainder")}
    for index in range(PANELS):
        left = -DELTA_RADIUS + 2 * DELTA_RADIUS * index / PANELS
        right = -DELTA_RADIUS + 2 * DELTA_RADIUS * (index + 1) / PANELS
        current = row(t, ball(left, right), pi, tau, beta, ymax)
        values = {
            "boundary_phase_remainder": uabs(current["rb"]),
            "linear_inner_remainder": uabs(current["cross"]),
            "quadratic_inner_remainder": uabs(current["quad"]),
            "complete_phase_remainder": phase_bound(current, ymax),
            "normalized_amplitude_remainder": max(uabs(current["amp0"] - 1), uabs(current["amp1"] - 1)),
        }
        for key, value in values.items():
            maxima[key] = max(maxima[key], value)
    require(row(t, -DELTA_RADIUS, pi, tau, beta, ymax)["u"].upper() < -200, "negative u coverage failed")
    require(row(t, DELTA_RADIUS, pi, tau, beta, ymax)["u"].lower() > 200, "positive u coverage failed")
    return maxima


def main() -> None:
    priority = set_low_priority()
    require(RESULT.is_file() and NOTE.is_file(), "missing overlap-core artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "overlap-core artifact is not passed")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file() and file_hash(path) == record["sha256"], f"hash drift: {path}")

    x, m, endpoint = sp.symbols("x m D", positive=True)
    q = sp.sqrt(x / 2) * (endpoint - 2 * m / x)
    require(sp.simplify(-sp.pi * m**2 / x + sp.pi * q**2 / 2 - (sp.pi * x * endpoint**2 / 4 - sp.pi * m * endpoint)) == 0, "independent completion failed")
    t, v = sp.symbols("t v", positive=True)
    rr = t / (2 * sp.pi * m**2)
    xv = 1 / (1 + rr * v)
    psi = -sp.pi * m**2 / xv + t * sp.log(rr * v) / 2
    require(sp.simplify(sp.expand_log(psi - psi.subs(v, 1) + t * (v - 1 - sp.log(v)) / 2, force=True)) == 0, "independent hyperbolic phase failed")

    ctx.dps = PRECISION
    pi = arb.pi()
    tstar = pi * C**2 / 8
    beta = tstar ** (arb(1) / 3)
    tau = pi * MODE * (C - 2 * MODE)
    ymax = arb(64) / (2 * beta).sqrt()
    heights = [tau - pi / 16, tau, tau + pi / 16]
    fresh = [scan(height, pi, tau, beta, ymax) for height in heights]
    saved = artifact["compact_overlap_core"]["uniform_bounds"]
    for key in fresh[0]:
        bound = max(record[key] for record in fresh)
        require(bound <= arb(saved[key]), f"independent {key} exceeds saved bound")
    require(max(record["complete_phase_remainder"] for record in fresh) < arb("0.082"), "independent phase target failed")
    require(max(record["normalized_amplitude_remainder"] for record in fresh) < arb("0.001"), "independent amplitude target failed")

    decisions = artifact["decision"]
    require(decisions["exact_hyperbolic_Morse_Fresnel_phase_proved"] is True, "exact phase decision missing")
    require(decisions["same_physical_alpha_slice_used"] is True, "alpha-slice decision missing")
    require(decisions["three_height_compact_overlap_core_proved"] is True, "compact-core decision missing")
    require(decisions["core_integral_remainder_proved"] is False, "integral remainder overpromoted")
    require(decisions["outer_u_tail_proved"] is False, "outer tail overpromoted")
    require(decisions["propagation_across_399_events_proved"] is False, "atlas propagation overpromoted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")

    note = NOTE.read_text(encoding="utf-8")
    for token in (
        "(Q^2-Q_0^2-u^2)/2",
        "same physical alpha",
        "does not integrate the exact-minus-canonical core",
        "not a proof of the complete Airy/logistic integral join",
    ):
        require(token in note, f"note token missing: {token}")
    print(
        "validated hyperbolic Morse-Fresnel overlap core independently: "
        "faces=3, |u|<=200, phase<0.082, amplitude<0.001; "
        f"priority={priority}"
    )


if __name__ == "__main__":
    main()
