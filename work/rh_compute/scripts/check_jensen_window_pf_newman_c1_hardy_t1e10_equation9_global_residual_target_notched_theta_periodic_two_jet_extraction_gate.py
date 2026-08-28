#!/usr/bin/env python3
"""Independently check the periodic two-jet target-notch extraction gate."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_periodic_two_jet_extraction_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
A = 159_577
B = 5_122_421
L = (B - A) // 2
T_LO = 622
T_HI = 39_894


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
    decisions = artifact["decision"]
    for key in (
        "theta_value_jump_extracted_exactly",
        "theta_first_derivative_jump_extracted_exactly",
        "two_jet_remainder_periodic_C1",
        "zero_x_extension_regular_and_residual_vanishes",
        "target_asymmetry_reduced_to_closed_harmonic_current",
        "differentiated_leakage_absolutely_convergent",
        "finite_regulator_and_limit_order_preserved",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in (
        "floating_data_used_as_proof",
        "A_B_extraction_cancellation_completed",
        "quantitative_R_after_A_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    s, x, alpha = sp.symbols("s x alpha", real=True)
    phase = sp.exp(sp.I * sp.pi * x * (alpha + 2 * s) ** 2 / 4)
    require(sp.simplify(sp.diff(phase, s) / (sp.I * sp.pi * x) - (alpha + 2 * s) * phase) == 0, "independent current derivative failed")

    source_sum = sp.Integer(L) * sp.Integer(A + L - 1)
    psi_zero = sp.factor((4 * s * source_sum + 4 * L * s**2 - s * (B * B - A * A)) / 4)
    require(psi_zero == L * s * (s - 1), "independent zero-x extension failed")
    p = s**2 - s
    require(sp.expand(psi_zero - sp.Rational(B - A, 2) * p) == 0, "independent zero-x residual failed")
    require(p.subs(s, 0) == p.subs(s, 1) == 0, "independent bridge endpoints failed")
    require(sp.diff(p, s).subs(s, 1) - sp.diff(p, s).subs(s, 0) == 2, "independent bridge jump failed")

    a, b = sp.symbols("a b", real=True)
    midpoint = sp.pi * x * (a**2 + b**2) / 8
    delta = sp.pi * x * (b**2 - a**2) / 8
    raw = (sp.exp(sp.I * sp.pi * x * b**2 / 4) - sp.exp(sp.I * sp.pi * x * a**2 / 4)) / (sp.I * sp.pi * x)
    stable = (b**2 - a**2) * sp.exp(sp.I * midpoint) * sp.sin(delta) / (4 * delta)
    require(sp.simplify(sp.expand_complex(raw - stable)) == 0, "independent endpoint divided difference failed")

    m = sp.symbols("m", integer=True, nonzero=True)
    z = 2 * sp.pi * sp.I * m
    poly_fourier = sp.integrate(p * sp.exp(-z * s), (s, 0, 1))
    require(sp.simplify(poly_fourier - 1 / (2 * sp.pi**2 * m**2)) == 0, "independent polynomial Fourier coefficient failed")

    ctx.dps = 95
    ctx.threads = 1
    pi = arb.pi()
    rows = artifact["harmonic_certificate"]["rows"]
    require([row["epsilon"] for row in rows] == ["0", "1e-6", "1e-8", "1e-10"], "harmonic row roster drift")
    for row in rows:
        epsilon_text = row["epsilon"]
        epsilon = arb(epsilon_text)
        total = arb(0)
        for k in range(T_LO, T_HI + 1):
            weight = arb(1) if epsilon_text == "0" else (-pi * epsilon * k * k).exp()
            total += weight / k
        require(arb(row["H_T_ball"]).overlaps(total), f"harmonic current drift at epsilon={epsilon_text}")
        require(arb(row["H_T_over_pi_ball"]).overlaps(total / pi), f"scaled harmonic current drift at epsilon={epsilon_text}")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("Psi_0(s)=L(s^2-s)" in note and "Phi_0(s)=0" in note, "zero-x extension missing")
    require("i*K_x*H_T(epsilon)/(2*pi)" in note, "closed harmonic current missing")
    require("No quantitative finite-current evaluator" in note, "proof boundary missing")
    print("independently checked periodic two-jet target-notch extraction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
