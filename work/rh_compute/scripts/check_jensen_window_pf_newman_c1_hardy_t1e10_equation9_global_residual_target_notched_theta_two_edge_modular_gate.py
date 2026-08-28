#!/usr/bin/env python3
"""Independently check the two-edge target-notch modular gate."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_two_edge_modular_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
PILOT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_dual_pilot.json"


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
        "target_complement_realized_as_half_integer_continuous_notch",
        "Poisson_dual_two_edge_transform_exact_at_fixed_positive_epsilon",
        "first_two_boundary_currents_summed_in_closed_form",
        "remaining_dual_series_absolutely_convergent_with_cubic_tail",
        "hybrid_direct_near_integer_and_dual_away_evaluator_licensed",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in (
        "floating_pilot_used_as_proof",
        "Kummer_source_integration_completed",
        "R_after_A_evaluated",
        "R_after_A_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    c = sp.Rational(1243, 2)
    d = sp.Rational(79789, 2)
    require(d - c == 39_273, "independent edge-width check failed")
    k = sp.symbols("k", integer=True)
    require(sp.simplify(sp.exp(-2 * sp.pi * sp.I * k * c) - (-1) ** k) == 0, "independent left phase check failed")
    require(sp.simplify(sp.exp(-2 * sp.pi * sp.I * k * d) - (-1) ** k) == 0, "independent right phase check failed")

    y, epsilon = sp.symbols("y epsilon", positive=True, real=True)
    q = sp.exp(-sp.pi * epsilon * y**2)
    q3 = sp.simplify(sp.diff(q, y, 3) / q)
    require(sp.simplify(q3 - (12 * sp.pi**2 * epsilon**2 * y - 8 * sp.pi**3 * epsilon**3 * y**3)) == 0, "independent q''' check failed")
    moment1 = 2 * sp.integrate(y * q, (y, 0, sp.oo))
    moment3 = 2 * sp.integrate(y**3 * q, (y, 0, sp.oo))
    majorant = sp.simplify(12 * sp.pi**2 * epsilon**2 * moment1 + 8 * sp.pi**3 * epsilon**3 * moment3)
    require(sp.simplify(majorant - 20 * sp.pi * epsilon) == 0, "independent q''' L1 majorant failed")

    transform = artifact["exact_transform_certificate"]
    require("erfc" in transform["Fourier_transform"], "Fourier-transform formula missing")
    require(transform["boundary_sum_0"] == "B_0(u)=E_0(u)/(2i*sin(pi*u))", "first boundary sum drift")
    require(transform["boundary_sum_1"] == "B_1(u)=-E_1(u)cos(pi*u)/(4sin(pi*u)^2)", "second boundary sum drift")
    require("/(2pi|xi|)^3" in transform["pointwise_remainder_bound"], "cubic tail drift")

    ctx.dps = 120
    ctx.threads = 1
    rows = artifact["tail_certificate"]["rows"]
    require([row["epsilon"] for row in rows] == ["1e-6", "1e-8", "1e-10"], "tail epsilon roster drift")
    pi = arb.pi()
    c_ball = arb(1243) / 2
    d_ball = arb(79789) / 2
    cutoff = arb(64)
    delta = arb(1) / 10
    reciprocal = arb(1) / (2 * (cutoff + delta) ** 2) + arb(1) / (2 * (cutoff - (1 - delta)) ** 2)
    for row in rows:
        eps = arb(row["epsilon"])

        def q2(edge: arb) -> arb:
            value = (-pi * eps * edge**2).exp()
            return (4 * pi**2 * eps**2 * edge**2 - 2 * pi * eps) * value

        k3 = abs(q2(c_ball)) + abs(q2(d_ball)) + 20 * pi * eps
        tail = k3 * reciprocal / (2 * pi) ** 3
        require(arb(row["K3_ball"]).overlaps(k3), f"K3 interval drift at epsilon={row['epsilon']}")
        require(arb(row["uniform_dual_tail_ball"]).overlaps(tail), f"tail interval drift at epsilon={row['epsilon']}")
    require(arb(rows[-1]["uniform_dual_tail_ball"]).upper() < arb("7e-15"), "small-epsilon dual tail guard failed")

    pilot = load_json(PILOT)
    require(pilot["summary"]["row_count"] == 9, "pilot row count drift")
    require(pilot["summary"]["largest_direct_cutoff"] == 328_014, "pilot direct cutoff drift")
    require(pilot["summary"]["fixed_dual_cutoff"] == 64, "pilot dual cutoff drift")
    require(pilot["summary"]["maximum_absolute_discrepancy"] < 4.19e-11, "pilot discrepancy drift")
    require(pilot["decision"]["interval_certified"] is False and pilot["decision"]["R_after_A_evaluated"] is False, "pilot was promoted improperly")
    require(artifact["exploratory_context"]["used_as_proof"] is False, "exploratory-context guard lost")
    require(file_hash(PILOT) == artifact["exploratory_context"]["sha256"], "pilot hash drift")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("c=621.5" in note and "d=39894.5" in note, "half-integer edges missing from note")
    require("20pi*epsilon" in note, "third-derivative majorant missing from note")
    require("No integration" in note, "proof boundary missing from note")
    print("independently checked two-edge target-notch modular transform", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
