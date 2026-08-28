#!/usr/bin/env python3
"""Independently replay the translated A-face endpoint interval theorem."""

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
import mpmath as mp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_peano_interval_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
A = 159_577
T = 10_000_000_000
X_LOWER_TEXT = "0.4999028779462474"
X_SLABS = 320
HALF_CELL_SUBDIVISIONS = 20
TERMS = 6
PRECISION = 120


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def set_low_priority() -> None:
    try:
        import psutil

        process = psutil.Process()
        process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS if os.name == "nt" else 10)
    except Exception:
        pass


def odd_df(index: int) -> int:
    value = 1
    for factor in range(1, index + 1, 2):
        value *= factor
    return value


def box(left: mp.mpf, right: mp.mpf) -> arb:
    middle = (left + right) / 2
    radius = (right - left) / 2 + mp.mpf("1e-125")
    return arb(mp.nstr(middle, 135), mp.nstr(radius, 135))


def finite(value: arb | acb) -> bool:
    return "nan" not in value.str(25, more=True).lower()


def symmetric_ball(radius_source: arb) -> arb:
    midpoint, radius, exponent = radius_source.upper().mid_rad_10exp()
    return arb(0, f"{midpoint + radius}e{exponent}")


def independent_gyy_bound(y: arb, x: arb) -> arb:
    pi = arb.pi()
    c = (2 / x).sqrt()
    delta = y - arb(A) * x / 2
    q = -delta * c
    q_abs = -q
    require(q_abs.lower() > 0, "independent q cover crossed zero")

    if q_abs.lower() >= 2:
        i_pi = acb(0, pi)
        direct_gyy = acb(arb(A) * x) / (i_pi * delta**3)
        for n in range(1, TERMS):
            power = 2 * n + 1
            coefficient_n = arb(odd_df(2 * n - 1)) * (x / 2) ** n / i_pi ** (n + 1)
            direct_gyy += coefficient_n * power * (
                arb(power + 1) * y / delta ** (power + 2)
                - 2 / delta ** (power + 1)
            )

        # Independent terminal-current recurrence, rather than the builder's
        # factored delta-polynomial correction.
        a_five = acb(arb(odd_df(9))) / i_pi**6
        first_ode_residual = -11 * a_five / q**12
        residual_derivative = 132 * a_five / q**13
        correction = (
            -acb(2 * c**2) * first_ode_residual
            - acb(y * c**3) * (acb(0, pi * q) * first_ode_residual - residual_derivative)
        )
        coefficient = -acb(0, 2 * pi * c**2 * q) + acb(y * c**3) * (acb(0, pi) + acb(pi**2 * q**2))
        remainder = 2 * arb(odd_df(11)) / (pi**7 * q_abs**13)
        approximation = direct_gyy + correction
        if finite(approximation) and finite(coefficient) and finite(remainder):
            result = abs(approximation) + abs(coefficient) * remainder
            if finite(result):
                return result

    q_center = q.mid()
    radius = q.rad()
    z_scale = (-q_center) * pi.sqrt() / 2
    z = acb(z_scale, -z_scale)
    h_center = acb(arb(1) / 2, arb(1) / 2) * (z**2).exp() * z.erfc()
    require(finite(h_center), "independent exact H center failed")
    phase_radius = pi * radius * (abs(q_center) + radius / 2)
    variation = radius + abs(h_center) * phase_radius
    error = symmetric_ball(variation)
    h = acb(h_center.real + error, h_center.imag + error)
    h_prime = 1 - acb(0, pi * q) * h
    h_second = -acb(0, pi * q) - (acb(0, pi) + acb(pi**2 * q**2)) * h
    result = abs(acb(2 * c**2) * h_prime - acb(y * c**3) * h_second)
    require(finite(result), "independent exact G_yy box failed")
    return result


def g_six_independent(y: arb, x: arb) -> acb:
    pi_i = acb(0, arb.pi())
    shift = y - arb(A) * x / 2
    answer = acb(arb(A) * x) / (2 * pi_i * shift)
    for n in range(1, TERMS):
        answer += acb(y) * arb(odd_df(2 * n - 1)) * (x / 2) ** n / (pi_i ** (n + 1) * shift ** (2 * n + 1))
    return answer


def primitive_difference_independent(left: arb, right: arb, x: arb) -> acb:
    pi_i = acb(0, arb.pi())
    shift = arb(A) * x / 2
    u = left - shift
    v = right - shift
    answer = acb(arb(A) * x) * (v / u).log() / (2 * pi_i)
    for n in range(1, TERMS):
        coefficient = arb(odd_df(2 * n - 1)) * (x / 2) ** n / pi_i ** (n + 1)
        primitive_v = v ** (1 - 2 * n) / (1 - 2 * n) - shift * v ** (-2 * n) / (2 * n)
        primitive_u = u ** (1 - 2 * n) / (1 - 2 * n) - shift * u ** (-2 * n) / (2 * n)
        answer += coefficient * (primitive_v - primitive_u)
    return answer


def tail_defect_independent(x: arb) -> acb:
    left = arb(79_805) / 2
    right = arb(79_873) / 2
    points = sum((g_six_independent(arb(mode), x) for mode in range(39_903, 39_937)), acb(0))
    return points - primitive_difference_independent(left, right, x)


def replay(localization: dict) -> arb:
    ctx.dps = PRECISION
    ctx.threads = 1
    mp.mp.dps = 150
    x0 = mp.mpf(X_LOWER_TEXT)
    x1 = mp.mpf("0.5")
    dx = (x1 - x0) / X_SLABS
    dy = mp.mpf("0.5") / HALF_CELL_SUBDIVISIONS
    normalization = 2 * (arb.pi() / (32 * arb(T))) ** (arb(1) / 4)
    eight_integral = arb(0)
    tail_integral = arb(0)

    for slab in range(X_SLABS):
        xb = box(x0 + slab * dx, x0 + (slab + 1) * dx)
        peano = arb(0)
        for mode in range(39_895, 39_903):
            left = mp.mpf(mode) - mp.mpf("0.5")
            for side in range(2):
                start = left + mp.mpf(side) / 2
                for part in range(HALF_CELL_SUBDIVISIONS):
                    yb = box(start + part * dy, start + (part + 1) * dy)
                    if side == 0:
                        kernel_upper = ((part + 1) * dy) ** 2 / 2
                    else:
                        kernel_upper = ((HALF_CELL_SUBDIVISIONS - part) * dy) ** 2 / 2
                    peano += arb(mp.nstr(dy * kernel_upper, 135)) * independent_gyy_bound(yb, xb).upper()
        weight = 1 / (xb * (xb * (1 - xb)) ** (arb(1) / 4))
        scale = normalization * arb(mp.nstr(dx, 135)) * weight.upper()
        eight_integral += scale * peano
        tail_bound = abs(tail_defect_independent(xb))
        require(finite(tail_bound), "independent tail defect failed")
        tail_integral += scale * tail_bound

    replacement = arb(localization["numerical_certificate"]["physical_tail_replacement_error_ball"])
    total = eight_integral + tail_integral + replacement
    require(eight_integral.upper() < arb("0.0028"), "independent eight-cell bound failed")
    require(tail_integral.upper() < arb("0.000030"), "independent rational-tail bound failed")
    require(total.upper() < arb("0.0030"), "independent complete endpoint bound failed")
    return total


def main() -> int:
    set_low_priority()
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decision = artifact["decision"]
    for key in (
        "eight_cell_transition_interval_bound_proved",
        "tail_six_current_interval_bound_included",
        "tail_replacement_error_included",
        "complete_A_face_endpoint_layer_absolute_bound_below_0_point_0037_proved",
    ):
        require(decision.get(key) is True, f"missing proved decision: {key}")
    for key in (
        "pre_endpoint_Morse_integral_proved",
        "full_signed_A_face_bound_proved",
        "R_after_A_bound_proved",
        "R_Dir_bound_proved",
        "QK_minus_T_bound_proved",
        "rh_implication",
    ):
        require(decision.get(key) is False, f"proof-boundary drift: {key}")

    certificate = artifact["interval_certificate"]
    builder_total = arb(certificate["physical_complete_endpoint_bound"])
    require(builder_total.upper() < arb("0.0037"), "builder endpoint threshold drift")
    require(certificate["partition"]["x_slabs"] == 256, "builder x partition drift")
    require(certificate["partition"]["half_cell_subdivisions"] == 16, "builder y partition drift")

    dependency_record = artifact["dependency"]
    localization_path = REPO_ROOT / dependency_record["path"]
    require(file_hash(localization_path) == dependency_record["sha256"], "localization dependency hash drift")
    localization = load_json(localization_path)
    independent_total = replay(localization)

    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("The 34-cell tail is not discarded" in note, "tail-accounting statement missing")
    require("< 0.0037" in note, "endpoint theorem threshold missing")
    require("does not prove the separate" in note and "RH" in note, "proof boundary missing")
    print(
        "independently checked complete translated A-face endpoint interval bound; "
        f"replay {independent_total.str(12, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
