#!/usr/bin/env python3
"""Independently check the post-A extended-notch pole-free kernel gate."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_pole_free_kernel_gate"
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
        "extended_upper_edge_instantiated_exactly",
        "shared_half_integer_dual_phase_certified",
        "pole_free_closed_centred_cell_formula_certified",
        "uniform_cubic_dual_tail_certified",
        "A_half_boundary_projector_collar_removed",
        "surviving_common_A_phase_has_interior_stationary_point",
    ):
        require(decision.get(key) is True, f"missing decision: {key}")
    for key in (
        "blanket_nonstationary_IBP_allowed",
        "Kummer_weighted_cell_integral_evaluated",
        "R_after_A_bound_proved",
        "R_Dir_bound_proved",
        "QK_minus_T_bound_proved",
        "rh_implication",
    ):
        require(decision.get(key) is False, f"proof-boundary drift: {key}")

    kernel = artifact["exact_kernel_certificate"]
    require(kernel["extended_target"] == [622, 39936], "extended target drift")
    require(kernel["target_count"] == 39315, "extended target count drift")
    require(kernel["half_integer_edges"] == ["1243/2", "79873/2"], "edge drift")
    require(kernel["edge_width"] == "39315", "edge width drift")

    z = sp.symbols("z", real=True)
    g0_core = sp.csc(z) - 1 / z
    g1 = (1 / z**2 - sp.cos(z) / sp.sin(z) ** 2) / 4
    require(sp.limit(g0_core, z, 0) == 0, "independent g0 limit failed")
    require(sp.limit(g1, z, 0) == sp.Rational(1, 24), "independent g1 limit failed")

    ctx.dps = 120
    ctx.threads = 1
    pi = arb.pi()
    c = arb(1243) / 2
    d_a = arb(79873) / 2
    cutoff = arb(64)
    for row in artifact["tail_certificate"]["rows"]:
        epsilon = arb(row["epsilon"])

        def q2(edge: arb) -> arb:
            q = (-pi * epsilon * edge**2).exp()
            return (4 * pi**2 * epsilon**2 * edge**2 - 2 * pi * epsilon) * q

        q2_c = q2(c)
        q2_d = q2(d_a)
        k3 = abs(q2_c) + abs(q2_d) + 20 * pi * epsilon
        tail = k3 / (2 * pi) ** 3 / (cutoff - arb(1) / 2) ** 2
        require(arb(row["q2_lower_edge_ball"]).overlaps(q2_c), f"lower q2 drift at {row['epsilon']}")
        require(arb(row["q2_upper_edge_ball"]).overlaps(q2_d), f"upper q2 drift at {row['epsilon']}")
        require(arb(row["K3_ball"]).overlaps(k3), f"K3 drift at {row['epsilon']}")
        require(arb(row["uniform_punctured_tail_ball"]).overlaps(tail), f"tail drift at {row['epsilon']}")
    require(arb(artifact["tail_certificate"]["rows"][-1]["uniform_punctured_tail_ball"]).upper() < arb("7.1e-15"), "tail threshold drift")

    discriminant = 1 - arb(8 * T) / (pi * arb(A) ** 2)
    x_star = (1 - discriminant.sqrt()) / 2
    curvature = arb(T) * (1 - 2 * x_star) / (2 * (x_star * (1 - x_star)) ** 2)
    guard = artifact["stationary_guard_certificate"]
    require(arb(guard["x_star_ball"]).overlaps(x_star), "stationary point drift")
    require(arb(guard["curvature_ball"]).overlaps(curvature), "stationary curvature drift")
    require(x_star.lower() > arb(0) and x_star.upper() < arb(1) / 2, "stationary point domain drift")
    require(curvature.lower() > arb(0), "stationary curvature sign drift")
    require("not on the notch edge d_A" in guard["edge_independence"], "edge-independence guard missing")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")

    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("d_A=39936.5" in note and "d_A-c=39315=|U|" in note, "extended-edge statement missing")
    require("g_0(0)=0" in note and "g_1(0)=1/24" in note, "pole-free limits missing")
    require("blanket\nnonstationary integration by parts is invalid" in note, "stationary route guard missing")
    require("No Kummer-weighted cell integral" in note, "proof boundary missing")
    print("independently checked post-A extended-notch pole-free kernel", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
