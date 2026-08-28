#!/usr/bin/env python3
"""Independently check the null-coordinate A-wedge unfolding gate."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_null_coordinate_unfolding_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
HEIGHT = 10_000_000_000
A = 159_577


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def check_symbolic() -> None:
    P, S, R, a, delta = sp.symbols("P S R a delta", real=True)
    phase = sp.expand(((P**2 - S**2) / 2).subs(P, R - S))
    require(sp.simplify(phase - (R**2 / 2 - R * S)) == 0, "fresh null phase failed")
    rho = delta - 1
    require(sp.simplify((P - a - rho * S).subs(P, R - S) - (R - a - delta * S)) == 0, "fresh face failed")

    t, endpoint = sp.symbols("t D", positive=True, real=True)
    dynamic = sp.sqrt(8 * t / sp.pi)
    eta = sp.sqrt(sp.pi) * (endpoint - dynamic) / 2
    alt = 4 * (sp.pi * endpoint**2 / 8 - t) / (sp.sqrt(sp.pi) * (endpoint + dynamic))
    require(sp.simplify(eta - alt) == 0, "fresh unfolding identity failed")
    require(sp.simplify(eta.subs(t, sp.pi * endpoint**2 / 8)) == 0, "fresh null height failed")

    m = endpoint / 4
    event_rho = -sp.sqrt(2) * sp.sqrt(endpoint - 2 * m) / sp.sqrt(endpoint)
    event_curvature = -2 * (-endpoint + 6 * m) / (3 * sp.sqrt(sp.pi) * endpoint ** sp.Rational(3, 2) * sp.sqrt(m))
    require(sp.simplify(event_rho + 1) == 0, "fresh null rho failed")
    require(sp.simplify(event_rho * event_curvature / 2 - 1 / (3 * endpoint * sp.sqrt(sp.pi))) == 0, "fresh cubic failed")


def main() -> int:
    ctx.dps = 130
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decisions = artifact["decision"]
    for key in (
        "near_null_wedge_reduced_to_one_dimensional_Abel_integrals",
        "reduction_avoids_division_by_characteristic_defect",
        "half_jump_retained_explicitly",
        "A_corner_unfolding_parameter_identified_exactly",
        "Airy_cubic_recovered_from_continuous_null_triangle",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in ("old_fold_atlas_embedded_in_global_projector", "canonical_A_wedge_numerically_bounded", "R_Dir_bound_proved", "rh_implication"):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")
    check_symbolic()

    geometry = artifact["geometry_certificate"]
    pi, t, endpoint = arb.pi(), arb(HEIGHT), arb(A)
    dynamic = (8 * t / pi).sqrt()
    eta = pi.sqrt() * (endpoint - dynamic) / 2
    require(eta.overlaps(arb(geometry["corner_unfolding_eta_A_ball"])), "saved eta misses fresh value")
    require(eta.lower() > arb("0.07784") and eta.upper() < arb("0.07785"), "fresh eta drift")
    rows = {row["mode"]: row for row in geometry["rows"]}
    require(set(rows) == {39_894, 39_895}, "mode roster drift")
    require(arb(rows[39_894]["delta_one_plus_rho_ball"]).upper() < 0, "fresh 39894 sign drift")
    require(arb(rows[39_895]["delta_one_plus_rho_ball"]).lower() > 0, "fresh 39895 sign drift")
    require(abs(arb(rows[39_894]["tangent_corner_minus_eta_ball"])).upper() < arb("1e-10"), "39894 tangent mismatch")
    require(abs(arb(rows[39_895]["tangent_corner_minus_eta_ball"])).upper() < arb("1e-9"), "39895 tangent mismatch")
    require(abs(arb(rows[39_894]["exact_curved_corner_minus_eta_ball"])).upper() < arb("1e-6"), "39894 curved mismatch")
    require(abs(arb(rows[39_895]["exact_curved_corner_minus_eta_ball"])).upper() < arb("1e-5"), "39895 curved mismatch")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("eta_D=sqrt(pi)(D-D_t)/2" in note, "note unfolding missing")
    require("No\ncanonical wedge value" in note, "note proof boundary missing")
    print("independently checked null-coordinate A-wedge unfolding", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
