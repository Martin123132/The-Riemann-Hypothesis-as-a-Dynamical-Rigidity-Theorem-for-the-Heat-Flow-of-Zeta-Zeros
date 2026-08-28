#!/usr/bin/env python3
"""Independently check the eight-cell A-face endpoint localization."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_transition_localization_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
A = 159_577
T = 10_000_000_000
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
        value += acb(y) * arb(odd_double_factorial(2 * index - 1)) * (x / 2) ** index / (i_pi ** (index + 1) * delta ** (2 * index + 1))
    return value


def integral_g_six(left: arb, right: arb, x: arb) -> acb:
    i_pi = acb(0, arb.pi())
    a = arb(A) * x / 2
    dl = left - a
    dr = right - a
    value = acb(arb(A) * x * (dr / dl).log()) / (2 * i_pi)
    for index in range(1, TERMS):
        coefficient = arb(odd_double_factorial(2 * index - 1)) * (x / 2) ** index / i_pi ** (index + 1)
        value += coefficient * (
            (dr ** (1 - 2 * index) - dl ** (1 - 2 * index)) / (1 - 2 * index)
            - a * (dr ** (-2 * index) - dl ** (-2 * index)) / (2 * index)
        )
    return value


def tail_defect_six(x: arb) -> acb:
    left = arb(79_805) / 2
    right = arb(79_873) / 2
    value = acb(0)
    for mode in range(39_903, 39_937):
        value += g_six(arb(mode), x)
    return value - integral_g_six(left, right, x)


def main() -> int:
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decision = artifact["decision"]
    for key in (
        "endpoint_defect_split_into_exact_eight_cell_transition_and_34_cell_tail",
        "tail_uniform_abs_q_at_least_16_point_5_proved",
        "tail_six_current_defect_error_below_5_point_56e_minus_9_proved",
        "physical_tail_replacement_error_below_6e_minus_15_proved",
        "special_function_obligation_reduced_to_eight_cells",
    ):
        require(decision.get(key) is True, f"missing decision: {key}")
    for key in (
        "eight_cell_transition_interval_value_proved",
        "complete_endpoint_layer_bound_proved",
        "full_signed_A_face_bound_proved",
        "R_after_A_bound_proved",
        "R_Dir_bound_proved",
        "QK_minus_T_bound_proved",
        "rh_implication",
    ):
        require(decision.get(key) is False, f"proof-boundary drift: {key}")

    full_left = sp.Rational(79_789, 2)
    tail_left = sp.Rational(79_805, 2)
    right = sp.Rational(79_873, 2)
    require(tail_left - full_left == 8 and right - tail_left == 34, "independent cell split drift")
    require(tuple(range(39_895, 39_903))[-1] == 39_902, "independent transition roster drift")
    require(tuple(range(39_903, 39_937))[-1] == 39_936, "independent tail roster drift")
    require(2 * (tail_left - sp.Rational(A, 4)) == sp.Rational(33, 2), "independent tail q floor drift")

    x, y = sp.symbols("x y", positive=True, real=True)
    delta = y - sp.Rational(A, 2) * x
    shape = y * x**TERMS / delta ** (2 * TERMS + 1)
    dx_form = x ** (TERMS - 1) * y * (TERMS * y + sp.Rational(A * (TERMS + 1), 2) * x) / delta ** (2 * TERMS + 2)
    dy_form = -x**TERMS * (2 * TERMS * y + sp.Rational(A, 2) * x) / delta ** (2 * TERMS + 2)
    require(sp.simplify(sp.diff(shape, x) - dx_form) == 0, "independent error x monotonicity drift")
    require(sp.simplify(sp.diff(shape, y) - dy_form) == 0, "independent error y monotonicity drift")
    weight = x ** (-sp.Rational(5, 4)) * (1 - x) ** (-sp.Rational(1, 4))
    require(sp.simplify(sp.diff(weight, x) - weight * (6 * x - 5) / (4 * x * (1 - x))) == 0, "independent physical weight drift")

    ctx.dps = 130
    ctx.threads = 1
    pi = arb.pi()
    aa = arb(A)
    tt = arb(T)
    ell = arb(79_789) / 2
    tail = arb(79_805) / 2
    qq = arb(16)
    sx = ((2 * qq**2 + 8 * aa * ell).sqrt() - qq * arb(2).sqrt()) / (2 * aa)
    x_cut = sx**2
    endpoint = arb(1) / 2
    delta_tail = tail - aa / 4
    point_error = 2 * arb(odd_double_factorial(11)) * tail * (endpoint / 2) ** TERMS / (pi**7 * delta_tail**13)
    defect_error = 68 * point_error
    normalization = 2 * (pi / (32 * tt)) ** (arb(1) / 4)
    weight_cut = 1 / (x_cut * (x_cut * (1 - x_cut)) ** (arb(1) / 4))
    physical_error = normalization * (endpoint - x_cut) * weight_cut * defect_error
    numerical = artifact["numerical_certificate"]
    require(arb(numerical["uniform_tail_pointwise_G_minus_G6_ball"]).overlaps(point_error), "tail point error drift")
    require(arb(numerical["uniform_tail_defect_error_ball"]).overlaps(defect_error), "tail defect error drift")
    require(arb(numerical["physical_tail_replacement_error_ball"]).overlaps(physical_error), "tail physical error drift")
    require(defect_error.upper() < arb("5.56e-9"), "independent tail defect threshold failed")
    require(physical_error.upper() < arb("6e-15"), "independent physical replacement threshold failed")

    at_cut = tail_defect_six(x_cut)
    at_endpoint = tail_defect_six(endpoint)
    for label, value in (("tail_D6_at_x16", at_cut), ("tail_D6_at_endpoint", at_endpoint)):
        record = numerical[label]
        require(arb(record["real_ball"]).overlaps(value.real), f"{label} real drift")
        require(arb(record["imag_ball"]).overlaps(value.imag), f"{label} imag drift")
        require(arb(record["modulus_ball"]).overlaps(abs(value)), f"{label} modulus drift")

    telemetry = artifact["floating_route_telemetry"]
    require(telemetry["used_as_proof"] is False, "floating telemetry proof-boundary drift")
    require(file_hash(REPO_ROOT / telemetry["path"]) == telemetry["sha256"], "floating telemetry hash drift")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")

    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("precisely the midpoint modes `39895,...,39902`" in note, "eight-cell roster missing")
    require("33/2=16.5" in note, "tail floor missing")
    require("Only" in note and "D_(tr,8)" in note, "remaining obligation missing")
    require("No interval value for the eight-cell transition" in note and "RH" in note, "proof boundary missing")
    print("independently checked eight-cell A-face endpoint localization", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
