#!/usr/bin/env python3
"""Independently check the post-A normal/tangential split gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_tangential_split_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
A = 159_577
T = 10_000_000_000


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
        "edge_translation_equals_explicit_42_mode_block",
        "zero_regulator_translation_is_finite_trigonometric_polynomial",
        "new_upper_edge_uniformly_nonstationary_in_A_face_normal_direction",
        "normal_denominator_gain_factor_169_certified",
        "tangential_A_Morse_stationary_point_persists",
        "normal_then_tangential_analysis_licensed",
    ):
        require(decision.get(key) is True, f"missing decision: {key}")
    for key in (
        "blanket_x_integration_by_parts_licensed",
        "signed_42_mode_cancellation_bound_proved",
        "R_after_A_bound_proved",
        "R_Dir_bound_proved",
        "QK_minus_T_bound_proved",
        "rh_implication",
    ):
        require(decision.get(key) is False, f"proof-boundary drift: {key}")

    translation = artifact["translation_certificate"]
    require(translation["shifted_block"] == [39895, 39936], "shifted block drift")
    require(translation["shifted_mode_count"] == 42, "shifted count drift")
    require(translation["edge_shift"] == "42", "edge shift drift")

    old_edge = Fraction(79789, 2)
    new_edge = Fraction(79873, 2)
    require(old_edge - Fraction(A, 4) == Fraction(1, 4), "independent old gap failed")
    require(new_edge - Fraction(A, 4) == Fraction(169, 4), "independent new gap failed")
    require((new_edge - Fraction(A, 4)) / (old_edge - Fraction(A, 4)) == 169, "independent gain failed")

    alpha, x, m = sp.symbols("alpha x m", positive=True)
    phase = sp.pi * alpha**2 * x / 4 - sp.pi * m * alpha + sp.Rational(T, 2) * sp.log((1 - x) / x)
    require(sp.simplify(sp.diff(phase, alpha) - sp.pi * (alpha * x / 2 - m)) == 0, "independent normal derivative failed")
    require(sp.simplify(sp.diff(phase.subs(alpha, A), x) - (sp.pi * A**2 / 4 - sp.Rational(T, 2) / (x * (1 - x)))) == 0, "independent tangential derivative failed")

    ctx.dps = 120
    ctx.threads = 1
    pi = arb.pi()
    discriminant = 1 - arb(8 * T) / (pi * arb(A) ** 2)
    x_star = (1 - discriminant.sqrt()) / 2
    old_face_gap = arb(79789) / 2 - arb(A) * x_star / 2
    new_face_gap = arb(79873) / 2 - arb(A) * x_star / 2
    curvature = arb(T) * (1 - 2 * x_star) / (2 * (x_star * (1 - x_star)) ** 2)
    geometry = artifact["phase_geometry_certificate"]
    require(arb(geometry["x_star_ball"]).overlaps(x_star), "x-star drift")
    require(arb(geometry["old_edge_gap_at_x_star_ball"]).overlaps(old_face_gap), "old face gap drift")
    require(arb(geometry["new_edge_gap_at_x_star_ball"]).overlaps(new_face_gap), "new face gap drift")
    require(arb(geometry["tangential_curvature_ball"]).overlaps(curvature), "curvature drift")
    require(old_face_gap.lower() > arb(42) and new_face_gap.lower() > arb(84), "face gap thresholds failed")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")

    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("explicit 42-mode" in note, "finite-block ownership guard missing")
    require("exact factor `169`" in note, "normal gain missing")
    require("Integrating by\nparts through `x_*` is invalid" in note, "tangential route guard missing")
    require("No endpoint-expansion remainder" in note, "proof boundary missing")
    print("independently checked post-A normal/tangential split", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
