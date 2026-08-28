#!/usr/bin/env python3
"""Independent finer-partition check of the outer phase-coupled B bound."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_phase_coupled_tangential_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_builder():
    spec = importlib.util.spec_from_file_location("outer_phase_gate", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load production formulas")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def subdivide(points: list[arb]) -> list[arb]:
    refined = [points[0]]
    for left, right in zip(points, points[1:]):
        refined.extend([(left + right) / 2, right])
    return refined


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    gate = load_builder()
    ctx.dps = 110
    ctx.threads = 1
    _, _, x_low, transition, left, right = gate.partition_geometry()
    transition, left, right = subdivide(transition), subdivide(left), subdivide(right)
    integral = sum((gate.slab_bound(a, b, True, True) for a, b in zip(transition, transition[1:])), arb(0))
    integral += sum((gate.slab_bound(a, b, True, False) for a, b in zip(left, left[1:])), arb(0))
    integral += sum((gate.slab_bound(a, b, False, False) for a, b in zip(right, right[1:])), arb(0))

    boundary0 = arb(0)
    boundary1 = arb(0)
    for point in (x_low, right[0], arb(1) / 2):
        first, second = gate.boundary_bounds(point)
        boundary0 += first
        boundary1 += second
    physical_factor = 2 * (arb.pi() / (32 * arb(gate.T))) ** (arb(1) / 4)
    physical = physical_factor * (boundary0 + boundary1 + integral)
    require(physical < arb("8e-10"), "independent refined outer bound exceeds 8e-10")
    require(arb(artifact["interval_certificate"]["physical_outer_bound_ball"]) < arb("8e-10"), "stored outer bound exceeds target")

    decision = artifact["decision"]
    require(decision["complete_analytic_m_ge_B_exterior_physical_contribution_below_8e_minus_10"] is True, "outer decision drift")
    require(decision["outer_block_separate_Abel_limit_justified_by_summable_derivative_majorants"] is True, "outer limit drift")
    require(decision["finite_mode_block_below_B_bounded"] is False, "finite-block overclaim")
    require(decision["completed_exterior_current_bound_proved"] is False, "completed-B overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print(f"independently checked refined phase-coupled outer B bound: {physical}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
