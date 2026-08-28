#!/usr/bin/env python3
"""Independently check the punctured-cell target-notch fold gate."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_fold_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
PILOT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_pilot.json"
A = 159_577
B = 5_122_421
L = 2_481_422


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
        "k0_residual_recombined_before_evaluation",
        "apparent_integer_lattice_poles_removed_exactly",
        "single_dual_formula_continuous_on_closed_centred_cell",
        "uniform_cubic_tail_on_closed_centred_cell_certified",
        "long_u_integral_folded_exactly_to_one_periodic_cell",
        "folded_source_is_exact_finite_quadratic_theta_current",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in (
        "floating_pilot_used_as_proof",
        "finite_Gauss_current_error_theorem_proved",
        "Kummer_x_integration_completed",
        "R_after_A_evaluated",
        "R_after_A_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    z = sp.symbols("z", real=True)
    g0_core = sp.csc(z) - 1 / z
    g1 = (1 / z**2 - sp.cos(z) / sp.sin(z) ** 2) / 4
    require(sp.limit(g0_core, z, 0) == 0, "independent g0 limit failed")
    require(sp.limit(g1, z, 0) == sp.Rational(1, 24), "independent g1 limit failed")
    require(sp.series(g0_core, z, 0, 6).coeff(z, 1) == sp.Rational(1, 6), "independent g0 series failed")
    require(sp.series(4 * g1, z, 0, 6).coeff(z, 2) == sp.Rational(7, 120), "independent g1 series failed")

    Hhat, B0, B1, b0, b1 = sp.symbols("Hhat B0 B1 b0 b1")
    require(sp.expand(B0 + B1 + (Hhat - b0 - b1) - (Hhat + B0 - b0 + B1 - b1)) == 0, "independent recombination failed")

    n, s, x = sp.symbols("n s x", real=True)
    alpha = A + 2 * n + 2 * s
    phase = sp.exp(sp.I * sp.pi * x * alpha**2 / 4)
    require(sp.simplify(sp.diff(phase, s) / (sp.I * sp.pi * x) - alpha * phase) == 0, "independent theta-current derivative failed")
    require((B - A) // 2 == L and A + 2 * (L - 1) == B - 2, "independent fold roster failed")

    fold = artifact["fold_certificate"]
    require(fold["roster_length"] == L, "fold length drift")
    require(fold["unshifted_fold_roster"] == [A, B - 2], "unshifted roster drift")
    require(fold["shifted_fold_roster_at_s_1"] == [A + 2, B], "shifted roster drift")
    require("integral_0^1" in fold["one_cell_fold"], "one-cell identity missing")

    ctx.dps = 120
    ctx.threads = 1
    pi = arb.pi()
    c = arb(1243) / 2
    d = arb(79789) / 2
    cutoff = arb(64)
    rows = artifact["tail_certificate"]["rows"]
    for row in rows:
        eps = arb(row["epsilon"])

        def q2(edge: arb) -> arb:
            q = (-pi * eps * edge**2).exp()
            return (4 * pi**2 * eps**2 * edge**2 - 2 * pi * eps) * q

        k3 = abs(q2(c)) + abs(q2(d)) + 20 * pi * eps
        tail = k3 / (2 * pi) ** 3 / (cutoff - arb(1) / 2) ** 2
        require(arb(row["K3_ball"]).overlaps(k3), f"K3 drift at epsilon={row['epsilon']}")
        require(arb(row["uniform_punctured_tail_ball"]).overlaps(tail), f"punctured tail drift at epsilon={row['epsilon']}")
    require(arb(rows[-1]["uniform_punctured_tail_ball"]).upper() < arb("7e-15"), "punctured tail threshold failed")

    pilot = load_json(PILOT)
    require(pilot["summary"]["row_count"] == 21 and pilot["summary"]["integer_rows"] == 3, "pilot row roster drift")
    require(pilot["summary"]["near_integer_rows"] == 9, "pilot near-integer roster drift")
    require(pilot["summary"]["maximum_relative_discrepancy"] < 3.91e-11, "pilot discrepancy drift")
    require(pilot["decision"]["interval_certified"] is False and pilot["decision"]["R_after_A_evaluated"] is False, "pilot promoted improperly")
    require(artifact["exploratory_context"]["used_as_proof"] is False, "exploratory guard lost")
    require(file_hash(PILOT) == artifact["exploratory_context"]["sha256"], "pilot hash drift")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("g_0(0)=0" in note and "g_1(0)=1/24" in note, "removable limits missing")
    require("integral_0^1 F_x(s)C_epsilon(s)ds" in note, "one-cell fold missing")
    require("No finite-Gauss-current error theorem" in note, "proof boundary missing")
    print("independently checked pole-free target-notch cell fold", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
