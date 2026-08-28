#!/usr/bin/env python3
"""Independently check the A-face midpoint-Peano partition gate."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_midpoint_Peano_partition_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
A = 159_577
T = 10_000_000_000
LEFT_NUMERATOR = 79_789
RIGHT_NUMERATOR = 79_873
Q_CUT = 16


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
        "fixed_regulator_42_mode_and_translated_strip_channel_oriented_exactly",
        "midpoint_Peano_identity_proved",
        "zero_regulator_compact_strip_limit_proved",
        "phase_stripped_normal_derivative_ODE_proved",
        "continuous_y_phase_cancellation_proved",
        "uniform_abs_q_at_least_16_before_x_16_proved",
        "endpoint_layer_strictly_after_tangential_saddle_proved",
        "endpoint_layer_phase_derivative_above_21000_proved",
    ):
        require(decision.get(key) is True, f"missing decision: {key}")
    for key in (
        "uniform_quantitative_G_yy_envelope_proved",
        "endpoint_layer_IBP_remainder_bounded",
        "tangential_Morse_integral_bounded",
        "signed_42_mode_cancellation_bound_proved",
        "R_after_A_bound_proved",
        "R_Dir_bound_proved",
        "QK_minus_T_bound_proved",
        "rh_implication",
    ):
        require(decision.get(key) is False, f"proof-boundary drift: {key}")

    left = sp.Rational(LEFT_NUMERATOR, 2)
    right = sp.Rational(RIGHT_NUMERATOR, 2)
    modes = tuple(range(39_895, 39_937))
    require(right - left == len(modes) == 42, "independent strip cardinality drift")
    for index, mode in enumerate(modes):
        require(sp.Rational(mode) == left + sp.Rational(2 * index + 1, 2), "independent midpoint roster drift")

    z = sp.symbols("z", real=True)
    half = sp.Rational(1, 2)
    k_minus = (half + z) ** 2 / 2
    k_plus = (half - z) ** 2 / 2
    mass = sp.integrate(k_minus, (z, -half, 0)) + sp.integrate(k_plus, (z, 0, half))
    require(mass == sp.Rational(1, 24) and 42 * mass == sp.Rational(7, 4), "independent Peano mass drift")
    for degree in range(7):
        f = z**degree
        direct = f.subs(z, 0) - sp.integrate(f, (z, -half, half))
        peano = -sp.integrate(k_minus * sp.diff(f, z, 2), (z, -half, 0)) - sp.integrate(k_plus * sp.diff(f, z, 2), (z, 0, half))
        require(sp.simplify(direct - peano) == 0, f"independent Peano polynomial check failed at degree {degree}")

    q = sp.symbols("q", real=True)
    h = sp.Function("H")(q)
    h_q = 1 - sp.I * sp.pi * q * h
    h_qq = sp.diff(h_q, q).subs(sp.diff(h, q), h_q)
    require(sp.simplify(h_qq + sp.I * sp.pi * q + (sp.I * sp.pi + sp.pi**2 * q**2) * h) == 0, "independent Fresnel ODE drift")

    y, c = sp.symbols("y c", positive=True, real=True)

    def d_y(expression: sp.Expr) -> sp.Expr:
        return sp.diff(expression, y) - c * sp.diff(expression, q)

    g = -1 / (sp.I * sp.pi) - y * c * h
    g_yy = sp.simplify(d_y(d_y(g)))
    expected_g_yy = 2 * c**2 * sp.diff(h, q) - y * c**3 * sp.diff(h, q, 2)
    require(sp.simplify(g_yy - expected_g_yy) == 0, "independent G_yy reduction drift")

    epsilon = sp.symbols("epsilon", nonnegative=True, real=True)
    gg = sp.Function("G")(y)
    weighted = sp.exp(-sp.pi * epsilon * y**2) * gg
    weighted_expected = sp.exp(-sp.pi * epsilon * y**2) * (
        sp.diff(gg, y, 2)
        - 4 * sp.pi * epsilon * y * sp.diff(gg, y)
        + (4 * sp.pi**2 * epsilon**2 * y**2 - 2 * sp.pi * epsilon) * gg
    )
    require(sp.simplify(sp.diff(weighted, y, 2) - weighted_expected) == 0, "independent weighted derivative drift")

    x = sp.symbols("x", positive=True, real=True)
    phase_completion = sp.expand((y - sp.Rational(A, 2) * x) ** 2 / x - y**2 / x + A * y)
    require(sp.simplify(phase_completion - sp.Rational(A**2, 4) * x) == 0, "independent continuous phase cancellation drift")

    normal_size = (y - sp.Rational(A, 2) * x) * sp.sqrt(2 / x)
    require(sp.diff(normal_size, y).is_positive is True, "independent y monotonicity drift")
    require(sp.simplify(sp.diff(normal_size, x) + sp.sqrt(2) * y / (2 * x ** sp.Rational(3, 2)) + A / (2 * sp.sqrt(2 * x))) == 0, "independent x monotonicity drift")
    require(-normal_size.subs({y: left, x: sp.Rational(1, 2)}) == -sp.Rational(1, 2), "independent left endpoint q drift")
    require(-normal_size.subs({y: right, x: sp.Rational(1, 2)}) == -sp.Rational(169, 2), "independent right endpoint q drift")

    ctx.dps = 120
    ctx.threads = 1
    pi = arb.pi()
    aa = arb(A)
    tt = arb(T)
    ell = arb(LEFT_NUMERATOR) / 2
    qq = arb(Q_CUT)
    sqrt_x_cut = ((2 * qq**2 + 8 * aa * ell).sqrt() - qq * arb(2).sqrt()) / (2 * aa)
    x_cut = sqrt_x_cut**2
    x_star = (1 - (1 - 8 * tt / (pi * aa**2)).sqrt()) / 2
    phase_prime = pi * aa**2 / 4 - tt / (2 * x_cut * (1 - x_cut))
    width = arb(1) / 2 - x_cut
    saddle_q = (ell - aa * x_star / 2) * (2 / x_star).sqrt()
    numerical = artifact["numerical_certificate"]
    require(arb(numerical["x_16_ball"]).overlaps(x_cut), "x_16 ball drift")
    require(arb(numerical["x_star_ball"]).overlaps(x_star), "x_star ball drift")
    require(arb(numerical["endpoint_layer_width_ball"]).overlaps(width), "endpoint width ball drift")
    require(arb(numerical["phase_prime_at_x_16_ball"]).overlaps(phase_prime), "phase derivative ball drift")
    require(arb(numerical["minimum_abs_q_on_strip_at_saddle_ball"]).overlaps(saddle_q), "saddle q ball drift")
    require(x_cut.lower() > x_star.upper(), "independent saddle partition failed")
    require(phase_prime.lower() > arb(21_000), "independent phase floor failed")
    require(width.upper() < arb("0.0001"), "independent endpoint width failed")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    telemetry = artifact["floating_route_telemetry"]
    telemetry_path = REPO_ROOT / telemetry["path"]
    require(telemetry["used_as_proof"] is False, "floating telemetry proof-boundary drift")
    require(file_hash(telemetry_path) == telemetry["sha256"], "floating telemetry hash drift")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")

    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("total kernel mass exactly `42/24 = 7/4`" not in note, "unexpected prose drift")
    require("sum_(42 cells) integral_cell K_m(y)dy=7/4" in note, "Peano mass missing")
    require("Phi_A'(x)>=Phi_A'(x_16)" in note, "endpoint phase floor missing")
    require("numbers select" in note and "neither is used in the proof" in note, "pilot boundary missing")
    require("No certified" in note and "RH" in note, "proof boundary missing")
    print("independently checked A-face midpoint-Peano reduction and partition", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
