#!/usr/bin/env python3
"""Independent check of the hyperbolic flat-face projector kernel."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_flat_face_projector_kernel_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interval_sign(value: arb) -> int:
    if value.lower() > 0:
        return 1
    if value.upper() < 0:
        return -1
    return 0


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    require(NOTE.is_file(), "missing theorem note")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    # Independent canonical algebra.
    p, s, a, rho = sp.symbols("p s a rho", real=True)
    c = rho**2 - 1
    boundary_phase = ((a + rho * s) ** 2 - s**2) / 2
    completed = c * (s + a * rho / c) ** 2 / 2 - a**2 / (2 * c)
    require(sp.simplify(boundary_phase - completed) == 0, "independent square completion failed")

    q, y, q0, kappa = sp.symbols("q y q0 kappa", real=True)
    require(sp.simplify(((sp.sqrt(sp.pi) * q) ** 2 - (sp.sqrt(2) * y) ** 2) / 2 - (-y**2 + sp.pi * q**2 / 2)) == 0, "independent phase dictionary failed")
    require(sp.simplify((kappa * sp.sqrt(sp.pi / 2)) ** 2 - 1 + (1 - sp.pi * kappa**2 / 2)) == 0, "independent defect dictionary failed")
    require(sp.simplify((sp.sqrt(sp.pi) * q - sp.sqrt(sp.pi) * q0 - kappa * sp.sqrt(sp.pi / 2) * sp.sqrt(2) * y) / sp.sqrt(sp.pi) - (q - q0 - kappa * y)) == 0, "independent boundary dictionary failed")

    # Check both branch constants in the exact U0 dictionary without asking
    # a CAS to choose an oscillatory square-root branch.
    require(sp.simplify((1 - sp.I) / 2 - 1 / (1 + sp.I)) == 0, "Delta-positive U0 constant failed")
    require(sp.simplify((1 + sp.I) / 2 - sp.I / (1 + sp.I)) == 0, "Delta-negative U0 constant failed")

    mode, endpoint, height = sp.symbols("mode endpoint height", positive=True)
    alpha = 2 * mode + height / (sp.pi * mode)
    rho2 = height * (endpoint + alpha) ** 2 / (2 * sp.pi * mode * alpha**3)
    event = sp.pi * mode * (endpoint - 2 * mode)
    require(sp.simplify((rho2 - 1).subs(height, event) - (1 - 4 * mode / endpoint)) == 0, "independent event defect failed")
    require(sp.simplify(event.subs(mode, endpoint / 4) - sp.pi * endpoint**2 / 8) == 0, "independent joint-fold height failed")

    # Recompute every saved-height scale row using Arb.
    ctx.dps = 100
    pi = arb.pi()
    t = arb(10_000_000_000)
    for row in artifact["geometry_certificate"]["rows"]:
        m = arb(row["mode"])
        endpoint_ball = arb(row["endpoint"])
        alpha_ball = 2 * m + t / (pi * m)
        p0 = (endpoint_ball - alpha_ball) / (endpoint_ball * alpha_ball).sqrt()
        a_ball = (pi * m * endpoint_ball).sqrt() * p0
        rho2_ball = t * (endpoint_ball + alpha_ball) ** 2 / (2 * pi * m * alpha_ball**3)
        c_ball = rho2_ball - 1
        scaled = a_ball / abs(c_ball).sqrt()
        require(a_ball.overlaps(arb(row["a_face_detuning_ball"])), f"a interval drift at m={row['mode']}")
        require(c_ball.overlaps(arb(row["c_characteristic_defect_ball"])), f"c interval drift at m={row['mode']}")
        require(scaled.overlaps(arb(row["a_over_sqrt_abs_c_ball"])), f"scaled interval drift at m={row['mode']}")
        require(interval_sign(c_ball) == row["c_sign"] != 0, f"c sign drift at m={row['mode']}")

    by_key = {(row["endpoint"], row["mode"]): row for row in artifact["geometry_certificate"]["rows"]}
    require(abs(arb(by_key[(5_122_421, 257)]["a_over_sqrt_abs_c_ball"])).lower() > arb("1.8e6"), "remote B characteristic audit failed")
    require(abs(arb(by_key[(159_577, 39_852)]["a_over_sqrt_abs_c_ball"])).upper() < arb("0.05"), "A lower face audit failed")
    require(abs(arb(by_key[(159_577, 39_853)]["a_over_sqrt_abs_c_ball"])).upper() < arb("0.08"), "A upper face audit failed")
    require(by_key[(159_577, 39_894)]["c_sign"] == 1 and by_key[(159_577, 39_895)]["c_sign"] == -1, "A null bracket failed")

    b_audit = artifact["inherited_B_audit"]
    require(b_audit["minimum_integer_abs_Delta_B_mode"] == 257, "B inherited minimum mode drift")
    require(arb(b_audit["minimum_integer_abs_Delta_B_ball"]).lower() > arb("0.0009"), "B inherited defect separation failed")
    require(arb(b_audit["mode_257_B_tail_argument_ball"]).lower() > arb("1e6"), "B inherited tail scale failed")

    decision = artifact["decision"]
    require(decision["flat_face_projector_profile_evaluated_exactly"] is True, "flat profile decision drift")
    require(decision["B_scalar_tangent_U0_is_same_canonical_kernel"] is True, "B dictionary decision drift")
    require(decision["B_q_density_U1_may_be_discarded"] is False, "q-density overclaim")
    require(decision["remote_B_slope_near_degeneracy_is_joint_fold"] is False, "remote B fold overclaim")
    require(decision["joint_fold_requires_face_incidence_and_null_slope"] is True, "joint-fold criterion drift")
    require(decision["curved_face_remainder_proved"] is False, "curved-face overclaim")
    require(decision["compressed_R_Dir_target_proved"] is False, "R_Dir overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently checked hyperbolic flat-face projector kernel and dictionary", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
