#!/usr/bin/env python3
"""Independently check the translated A-face two-current expansion gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_two_current_expansion_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
A = 159_577
D_NUMERATOR = 79_873


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decision = artifact["decision"]
    for key in (
        "continuous_edge_A_tail_defined_by_exact_finite_endpoint_current",
        "first_two_normal_currents_derived_exactly",
        "uniform_abs_q_at_least_84_point_5_certified",
        "two_current_pointwise_remainder_below_4e_minus_6",
        "half_integer_phase_equals_minus_i",
        "translated_edge_adjacent_to_exact_42_mode_block",
    ):
        require(decision.get(key) is True, f"missing decision: {key}")
    for key in (
        "leading_current_small_enough_for_separate_norm",
        "signed_edge_42_mode_cancellation_bound_proved",
        "tangential_Morse_integral_bounded",
        "R_after_A_bound_proved",
        "R_Dir_bound_proved",
        "QK_minus_T_bound_proved",
        "rh_implication",
    ):
        require(decision.get(key) is False, f"proof-boundary drift: {key}")

    x, d, delta = sp.symbols("x d delta", positive=True, real=True)
    q = -delta * sp.sqrt(2 / x)
    a = d * sp.sqrt(2 / x)
    first = -1 / (sp.I * sp.pi) - a / (sp.I * sp.pi * q)
    second = -a / ((sp.I * sp.pi) ** 2 * q**3)
    require(sp.simplify(first.subs(delta, d - A * x / 2) - A * x / (2 * sp.I * sp.pi * (d - A * x / 2))) == 0, "independent first-current reduction failed")
    require(sp.simplify(second + d * x / (2 * sp.pi**2 * delta**3)) == 0, "independent second-current reduction failed")

    d_exact = sp.Rational(D_NUMERATOR, 2)
    delta_exact = d_exact - sp.Rational(A, 2) * x
    q_exact = -delta_exact * sp.sqrt(2 / x)
    require(sp.simplify(q_exact.subs(x, sp.Rational(1, 2)) + sp.Rational(169, 2)) == 0, "independent q endpoint failed")
    require((A * D_NUMERATOR) % 4 == 1, "independent half-integer phase failed")

    ctx.dps = 120
    ctx.threads = 1
    pi = arb.pi()
    aa = arb(A)
    dd = arb(D_NUMERATOR) / 2
    xx = arb(1) / 2
    gap = dd - aa * xx / 2
    q_abs = gap * (2 / xx).sqrt()
    c0 = aa * xx / (2 * pi * gap)
    c1 = dd * xx / (2 * pi**2 * gap**3)
    remainder = 3 * dd * xx**2 / (2 * pi**3 * gap**5)
    numerical = artifact["numerical_certificate"]
    require(arb(numerical["corner_delta_ball"]).overlaps(gap), "corner gap drift")
    require(arb(numerical["minimum_abs_q_ball"]).overlaps(q_abs), "q bound drift")
    require(arb(numerical["maximum_first_current_modulus_ball"]).overlaps(c0), "first current bound drift")
    require(arb(numerical["maximum_second_current_modulus_ball"]).overlaps(c1), "second current bound drift")
    require(arb(numerical["uniform_two_current_remainder_ball"]).overlaps(remainder), "remainder bound drift")
    require(remainder.upper() < arb("4e-6"), "remainder threshold failed")
    require(c0.lower() > arb(300), "leading-current scale guard failed")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")

    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("exp(-i*pi*d_A*A)=-i" in note, "common phase orientation missing")
    require("starts one normal-coordinate unit" in note, "roster adjacency missing")
    require("does not make the leading current small" in note, "scale guard missing")
    require("No signed edge/42-mode cancellation bound" in note, "proof boundary missing")
    print("independently checked translated A-face two-current expansion", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
