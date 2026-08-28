#!/usr/bin/env python3
"""Independently check the six-current A-face midpoint-defect gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_midpoint_defect_reduction_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
A = 159_577
T = 10_000_000_000
LEFT_NUMERATOR = 79_789
RIGHT_NUMERATOR = 79_873
FIRST_MODE = 39_895
LAST_MODE = 39_936
TERMS = 6


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def odd_double_factorial(index: int) -> int:
    result = 1
    for value in range(1, index + 1, 2):
        result *= value
    return result


def g_six(y: arb, x: arb) -> acb:
    i_pi = acb(0, arb.pi())
    delta = y - arb(A) * x / 2
    value = acb(arb(A) * x) / (2 * i_pi * delta)
    for index in range(1, TERMS):
        value += (
            acb(y)
            * arb(odd_double_factorial(2 * index - 1))
            * (x / 2) ** index
            / (i_pi ** (index + 1) * delta ** (2 * index + 1))
        )
    return value


def integral_g_six(x: arb) -> acb:
    i_pi = acb(0, arb.pi())
    left = arb(LEFT_NUMERATOR) / 2
    right = arb(RIGHT_NUMERATOR) / 2
    a = arb(A) * x / 2
    dl = left - a
    dr = right - a
    value = acb(arb(A) * x * (dr / dl).log()) / (2 * i_pi)
    for index in range(1, TERMS):
        b = arb(odd_double_factorial(2 * index - 1)) * (x / 2) ** index / i_pi ** (index + 1)
        primitive = (
            (dr ** (1 - 2 * index) - dl ** (1 - 2 * index)) / (1 - 2 * index)
            - a * (dr ** (-2 * index) - dl ** (-2 * index)) / (2 * index)
        )
        value += b * primitive
    return value


def defect_six(x: arb) -> acb:
    value = acb(0)
    for mode in range(FIRST_MODE, LAST_MODE + 1):
        value += g_six(arb(mode), x)
    return value - integral_g_six(x)


def main() -> int:
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decision = artifact["decision"]
    for key in (
        "six_current_Fresnel_tail_with_finite_remainder_proved",
        "phase_stripped_G6_rational_formula_proved",
        "G6_strip_integral_reduced_to_rational_log_endpoints",
        "uniform_G_minus_G6_below_1_point_22e_minus_10_proved",
        "uniform_signed_defect_error_below_1_point_025e_minus_8_proved",
        "D6_over_x_regular_at_zero_proved",
    ):
        require(decision.get(key) is True, f"missing decision: {key}")
    for key in (
        "Morse_core_integral_bounded",
        "endpoint_layer_IBP_remainder_bounded",
        "full_signed_A_face_bound_proved",
        "R_after_A_bound_proved",
        "R_Dir_bound_proved",
        "QK_minus_T_bound_proved",
        "rh_implication",
    ):
        require(decision.get(key) is False, f"proof-boundary drift: {key}")

    x, y, delta = sp.symbols("x y delta", positive=True, real=True)
    c = sp.sqrt(2 / x)
    q = -c * delta
    h_six = sum(
        sp.factorial2(2 * index - 1) / ((sp.I * sp.pi) ** (index + 1) * q ** (2 * index + 1))
        for index in range(TERMS)
    )
    direct = -1 / (sp.I * sp.pi) - y * c * h_six
    reduced = A * x / (2 * sp.I * sp.pi * delta)
    for index in range(1, TERMS):
        reduced += y * sp.factorial2(2 * index - 1) * (x / 2) ** index / ((sp.I * sp.pi) ** (index + 1) * delta ** (2 * index + 1))
    require(sp.simplify(direct.subs(y, delta + sp.Rational(A, 2) * x) - reduced.subs(y, delta + sp.Rational(A, 2) * x)) == 0, "independent G6 reduction drift")

    a = sp.symbols("a", positive=True, real=True)
    for index in range(1, TERMS):
        dy = y - a
        primitive = dy ** (1 - 2 * index) / (1 - 2 * index) - a * dy ** (-2 * index) / (2 * index)
        require(sp.simplify(sp.diff(primitive, y) - y / dy ** (2 * index + 1)) == 0, f"independent strip primitive drift at {index}")

    shape = y * x**TERMS / (y - sp.Rational(A, 2) * x) ** (2 * TERMS + 1)
    denominator = (y - sp.Rational(A, 2) * x) ** (2 * TERMS + 2)
    dx_positive_form = x ** (TERMS - 1) * y * (TERMS * y + sp.Rational(A * (TERMS + 1), 2) * x) / denominator
    dy_negative_form = -x**TERMS * (2 * TERMS * y + sp.Rational(A, 2) * x) / denominator
    require(sp.simplify(sp.diff(shape, x) - dx_positive_form) == 0, "independent x monotonicity factorization failed")
    require(sp.simplify(sp.diff(shape, y) - dy_negative_form) == 0, "independent y monotonicity factorization failed")

    ctx.dps = 130
    ctx.threads = 1
    pi = arb.pi()
    aa = arb(A)
    tt = arb(T)
    left = arb(LEFT_NUMERATOR) / 2
    right = arb(RIGHT_NUMERATOR) / 2
    qq = arb(16)
    sx = ((2 * qq**2 + 8 * aa * left).sqrt() - qq * arb(2).sqrt()) / (2 * aa)
    x_cut = sx**2
    x_star = (1 - (1 - 8 * tt / (pi * aa**2)).sqrt()) / 2
    delta_cut = left - aa * x_cut / 2
    point_error = 2 * arb(odd_double_factorial(11)) * left * (x_cut / 2) ** TERMS / (pi**7 * delta_cut**13)
    defect_error = 84 * point_error
    numerical = artifact["numerical_certificate"]
    require(arb(numerical["uniform_pointwise_G_minus_G6_ball"]).overlaps(point_error), "pointwise error ball drift")
    require(arb(numerical["uniform_D42_minus_D6_ball"]).overlaps(defect_error), "defect error ball drift")
    require(point_error.upper() < arb("1.22e-10"), "independent pointwise threshold failed")
    require(defect_error.upper() < arb("1.025e-8"), "independent defect threshold failed")

    saddle = defect_six(x_star)
    handoff = defect_six(x_cut)
    saddle_record = numerical["six_current_saddle_defect"]
    handoff_record = numerical["six_current_x16_defect"]
    require(arb(saddle_record["real_ball"]).overlaps(saddle.real), "saddle real drift")
    require(arb(saddle_record["imag_ball"]).overlaps(saddle.imag), "saddle imag drift")
    require(arb(saddle_record["modulus_ball"]).overlaps(abs(saddle)), "saddle modulus drift")
    require(arb(handoff_record["real_ball"]).overlaps(handoff.real), "handoff real drift")
    require(arb(handoff_record["imag_ball"]).overlaps(handoff.imag), "handoff imag drift")
    require(arb(handoff_record["modulus_ball"]).overlaps(abs(handoff)), "handoff modulus drift")

    sum_1 = sum((arb(1) / mode for mode in range(FIRST_MODE, LAST_MODE + 1)), arb(0))
    sum_2 = sum((arb(1) / arb(mode) ** 2 for mode in range(FIRST_MODE, LAST_MODE + 1)), arb(0))
    i_pi = acb(0, pi)
    endpoint = acb(aa) * (sum_1 - (right / left).log()) / (2 * i_pi) + (sum_2 - (1 / left - 1 / right)) / (2 * i_pi**2)
    endpoint_record = numerical["regularized_x0_D6_over_x_limit"]
    require(arb(endpoint_record["real_ball"]).overlaps(endpoint.real), "endpoint real drift")
    require(arb(endpoint_record["imag_ball"]).overlaps(endpoint.imag), "endpoint imag drift")
    require(arb(endpoint_record["modulus_ball"]).overlaps(abs(endpoint)), "endpoint modulus drift")
    require(abs(endpoint).upper() < arb("1.4e-9"), "independent endpoint regularity threshold failed")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")

    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("finite rational-log" in note, "elementary reduction missing")
    require("The factor `84`" in note, "defect remainder accounting missing")
    require("is harmless at the" in note, "lower-endpoint regularity missing")
    require("No Morse-core integral" in note and "RH" in note, "proof boundary missing")
    print("independently checked six-current A-face midpoint defect", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
