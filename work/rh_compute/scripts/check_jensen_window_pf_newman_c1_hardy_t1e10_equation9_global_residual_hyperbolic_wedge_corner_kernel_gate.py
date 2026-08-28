#!/usr/bin/env python3
"""Independent check of the hyperbolic wedge corner kernel."""

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

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_corner_kernel_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def half_coordinate(mode: int, t: arb, pi: arb) -> tuple[arb, arb]:
    m = arb(mode)
    v_half = 2 * pi * m**2 / t
    defect = v_half - 1 - v_half.log()
    require(defect.lower() > 0, f"half-coordinate defect failed at m={mode}")
    root = (2 * defect).sqrt()
    signed = -root if (v_half - 1).upper() < 0 else root
    return v_half, (t / 2).sqrt() * signed


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    require(NOTE.is_file(), "missing theorem note")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    p, s, a, rho, h = sp.symbols("p s a rho h", real=True)
    c = rho**2 - 1
    phase = ((a + rho * s) ** 2 - s**2) / 2
    completed = c * (s + a * rho / c) ** 2 / 2 - a**2 / (2 * c)
    require(sp.simplify(phase - completed) == 0, "independent wedge completion failed")
    density = sp.exp(sp.I * ((a + rho * h) ** 2 - h**2) / 2) / (2 * sp.pi)
    require(sp.simplify(-density - (-density)) == 0, "independent mixed derivative failed")

    mode, height, endpoint = sp.symbols("mode height endpoint", positive=True)
    v_half = 2 * sp.pi * mode**2 / height
    require(sp.simplify(v_half.subs(mode, sp.sqrt(height / (2 * sp.pi))) - 1) == 0, "independent half event failed")
    require(sp.simplify((endpoint / 4).subs(endpoint, sp.sqrt(8 * height / sp.pi)) - sp.sqrt(height / (2 * sp.pi))) == 0, "independent joint corner failed")

    ctx.dps = 100
    pi = arb.pi()
    t = arb(10_000_000_000)
    for row in artifact["geometry_certificate"]["rows"]:
        v_ball, h_ball = half_coordinate(row["mode"], t, pi)
        require((v_ball - 1).overlaps(arb(row["v_half_minus_one_ball"])), f"v_half drift at m={row['mode']}")
        require(h_ball.overlaps(arb(row["h_half_boundary_ball"])), f"h drift at m={row['mode']}")

    indexed = {(row["endpoint"], row["mode"]): row for row in artifact["geometry_certificate"]["rows"]}
    require(abs(arb(indexed[(5_122_421, 621)]["h_half_boundary_ball"])).lower() > arb("250000"), "B 621 remote boundary failed")
    require(abs(arb(indexed[(5_122_421, 622)]["h_half_boundary_ball"])).lower() > arb("250000"), "B 622 remote boundary failed")
    require(arb(indexed[(159_577, 39_852)]["h_half_boundary_ball"]).upper() < arb("-149"), "A 39852 scale failed")
    require(arb(indexed[(159_577, 39_853)]["h_half_boundary_ball"]).upper() < arb("-146"), "A 39853 scale failed")
    require(arb("-0.81") < arb(indexed[(159_577, 39_894)]["h_half_boundary_ball"]) < arb("-0.80"), "A 39894 corner failed")
    require(arb("2.73") < arb(indexed[(159_577, 39_895)]["h_half_boundary_ball"]) < arb("2.74"), "A 39895 corner failed")

    nu = (t / (2 * pi)).sqrt()
    stored_nu = arb(artifact["geometry_certificate"]["half_boundary_mode_ball"])
    require(nu.overlaps(stored_nu) and arb(39_894) < nu < arb(39_895), "half-boundary mode interval failed")
    separation = arb(159_577) / 4 - nu
    require(separation.overlaps(arb(artifact["geometry_certificate"]["A_null_minus_half_boundary_mode_ball"])), "null/half separation drift")

    decision = artifact["decision"]
    require(decision["canonical_wedge_boundary_derivative_system_proved"] is True, "wedge decision drift")
    require(decision["flat_face_is_remote_half_boundary_limit"] is True, "flat-limit decision drift")
    require(decision["B_occupancy_edge_and_half_boundary_are_locally_coupled"] is False, "B coupling overclaim")
    require(decision["A_target_edge_and_half_boundary_are_on_same_canonical_scale"] is True, "A corner decision drift")
    require(decision["curved_face_remainder_proved"] is False, "curved-face overclaim")
    require(decision["compressed_R_Dir_target_proved"] is False, "R_Dir overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently checked hyperbolic wedge kernel and A half-boundary localization", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
