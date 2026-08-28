#!/usr/bin/env python3
"""Independently replay the grouped m3=0|1 first-wall bridge."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
import flint
from flint import acb, arb


BUILDER = Path(__file__).with_name(Path(__file__).name.removeprefix("check_"))
STEM = BUILDER.stem
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
ALTERED_BLOCK_SIZE = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_acb(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def independent_phase_moments(
    base_a: Fraction,
    base_tau: Fraction,
    n: int,
    block_size: int,
    prior: object,
) -> list[acb]:
    """Use independently ordered exact-reset blocks and Horner-style powers."""

    values = [acb(0) for _ in range(4)]
    starts = list(range(0, n + 1, block_size))
    for start in reversed(starts):
        stop = min(n + 1, start + block_size)
        term = prior.acb_cis_pi(2 * (base_a * start + base_tau * start * start))
        ratio = prior.acb_cis_pi(2 * (base_a + base_tau * (2 * start + 1)))
        curvature_ratio = prior.acb_cis_pi(4 * base_tau)
        local = [acb(0) for _ in range(4)]
        for k in range(start, stop):
            local[0] += term
            local[1] += k * term
            local[2] += k * k * term
            local[3] += k * k * k * term
            term *= ratio
            ratio *= curvature_ratio
        for degree in range(4):
            values[degree] += local[degree]
    require(all(value.is_finite() for value in values), "independent phase moments are nonfinite")
    return values


def independent_jump(side: str, builder: object) -> dict[str, Fraction]:
    n = builder.CHILD_N
    t = Fraction(1, 2 * n)
    if side == "right":
        p = Fraction(95_747, 6)
        z = Fraction(-896_819, 6 * n)
        expected_child_p = Fraction(1, 3)
    else:
        p = Fraction(47_873, 3)
        z = Fraction(544_156, 3 * n)
        expected_child_p = Fraction(-1, 3)
    child_p = 2 * t * p - z
    require(child_p == expected_child_p, f"{side} independent child p mismatch")
    u0 = z + (2 * n + 1) * t - Fraction(1, 2)
    u1 = z + (2 * n + 1) * t - 1 - Fraction(1, 2)
    require(u0 - u1 == 1, f"{side} unit recurrence mismatch")
    endpoint_weight = 2 * p + 2 * n + 1 - (u1 + Fraction(1, 2)) / t
    require(endpoint_weight == (child_p + 1) / t, f"{side} endpoint coefficient mismatch")
    w = z / (2 * t)
    sigma = -Fraction(1, 4 * t)
    upper_phase = 2 * (Fraction(n) + Fraction(1, 2)) * (
        z + t * (Fraction(n) + Fraction(1, 2))
    )
    endpoint_phase = Fraction(1, 4) + upper_phase - (u1 + Fraction(1, 2)) ** 2 / (2 * t)
    child_phase = Fraction(5, 4) - z * z / (2 * t) + 2 * (w + sigma)
    difference = endpoint_phase - child_phase
    require(difference == 2 * n and difference % 2 == 0, f"{side} jump phase mismatch")
    return {
        "child_p": child_p,
        "u0": u0,
        "u1": u1,
        "endpoint_weight": endpoint_weight,
        "phase_difference": difference,
    }


def main() -> int:
    require(BUILDER.is_file() and RESULT.is_file() and NOTE.is_file(), "missing production artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact not passed")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(Path(__file__)), "checker hash drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], "dependency hash drift")

    builder = load_module("check_grouped_wall_builder", BUILDER)
    builder.set_low_priority()
    adaptive = load_module("check_grouped_wall_adaptive", builder.ADAPTIVE_BUILDER)
    rational = load_module("check_grouped_wall_rational", adaptive.RATIONAL_BUILDER)
    two_sided = load_module("check_grouped_wall_two_sided", rational.TWO_SIDED_BUILDER)
    right_dependency = load_module("check_grouped_wall_right", two_sided.RIGHT_BUILDER)
    small = load_module("check_grouped_wall_small", right_dependency.SMALL_BUILDER)
    box = load_module("check_grouped_wall_box", right_dependency.BOX_BUILDER)
    prior = load_module("check_grouped_wall_prior", right_dependency.PRIOR_BUILDER)
    altered = box.Settings(precision_bits=384, cutoff=11, abs_tol="1e-65")
    flint.ctx.prec = altered.precision_bits

    for side in ("right", "left"):
        audit = independent_jump(side, builder)
        production_jump = artifact["exact_grouped_floor_jump"][side]
        require(str(audit["child_p"]) == production_jump["child_p_at_wall"], "child p audit drift")
        require(str(audit["u0"]) == production_jump["upper_argument_m0"], "m0 argument audit drift")
        require(str(audit["u1"]) == production_jump["upper_argument_m1"], "m1 argument audit drift")
        require(
            str(audit["phase_difference"]) == production_jump["endpoint_minus_child_phase_difference"],
            "phase audit drift",
        )

        base_a = builder.RIGHT_A if side == "right" else builder.LEFT_A
        independent = independent_phase_moments(
            base_a, builder.TAU0, builder.CHILD_N, ALTERED_BLOCK_SIZE, prior
        )
        production_source = artifact[f"{side}_first_wall_source_box"]
        production_moments = production_source["parent"]["common_wall_taylor_child"]["F_balls"]
        for degree, value in enumerate(independent):
            require(
                value.overlaps(parse_acb(production_moments[f"F{degree}"])),
                f"{side} F{degree} altered block recurrence mismatch",
            )

        replay = builder.source_wall_box(
            side,
            altered,
            box,
            small,
            prior,
            rational,
            adaptive,
            ALTERED_BLOCK_SIZE,
        )
        replay_source = parse_acb(replay["complete_source_box"])
        production_ball = parse_acb(production_source["complete_source_box"])
        require(replay_source.overlaps(production_ball), f"{side} altered source replay mismatch")
        require(
            arb(replay["complete_source_radius_upper"]) < arb(builder.SOURCE_RADIUS_TARGET),
            f"{side} altered source replay misses target",
        )
        floor_box = replay["parent"]["third_floor_argument_box"]
        require(Fraction(floor_box["lo"]) < 1 < Fraction(floor_box["hi"]), "replay does not cross floor wall")

    note = NOTE.read_text(encoding="utf-8")
    require("Delta upper endpoint = - Delta third main" in note, "note omits grouped jump")
    require("No cover between the normalization wall" in note, "note omits proof boundary")
    require(artifact["decision"]["rh_implication"] is False, "RH implication flag drift")
    print("independently replayed grouped m3=0|1 first-wall bridge", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
