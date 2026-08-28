#!/usr/bin/env python3
"""Independently replay the mirrored adaptive-width crossover gate."""

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


def same_box(record: dict[str, str], value: object) -> bool:
    return Fraction(record["lo"]) == value.lo and Fraction(record["hi"]) == value.hi


def parse_acb(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def independent_left_child_boxes(units: int, builder: object, box: object):
    x_lo = builder.X0 - units * builder.WIDTH_UNIT
    x_hi = builder.X0 - builder.WIDTH_UNIT
    s_lo = builder.S0 - builder.WIDTH_UNIT
    s_hi = builder.S0 + builder.WIDTH_UNIT
    c_lo = Fraction(3, 2) + Fraction(builder.A, 2) + s_lo
    c_hi = Fraction(3, 2) + Fraction(builder.A, 2) + s_hi
    p_box = box.RationalBox(
        3 * builder.FIRST_P - Fraction(3, 2) - Fraction(builder.A, 2) - s_hi,
        3 * builder.FIRST_P - Fraction(3, 2) - Fraction(builder.A, 2) - s_lo,
    )
    raw_a_lo = (c_lo * x_lo - builder.FIRST_P) / (3 * x_lo - 1)
    raw_a_hi = (c_hi * x_hi - builder.FIRST_P) / (3 * x_hi - 1)
    raw_tau_lo = -x_lo / (2 * (3 * x_lo - 1))
    raw_tau_hi = -x_hi / (2 * (3 * x_hi - 1))
    return (
        p_box,
        box.RationalBox(1 - raw_a_hi, 1 - raw_a_lo),
        box.RationalBox(-raw_tau_hi - 1, -raw_tau_lo - 1),
    )


def independent_wall_period(a0: Fraction, n: int) -> tuple[int, int]:
    tau = Fraction(1, 2 * n)
    for q in range(1, 13):
        period = n * q
        require((2 * tau * period).denominator == 1, "period is not a tau multiple")
        if (a0 * period + tau * period * period).denominator == 1:
            return period, q
    raise RuntimeError("independent wall period search failed")


def main() -> int:
    require(BUILDER.is_file() and RESULT.is_file() and NOTE.is_file(), "missing production artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact not passed")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(Path(__file__)), "checker hash drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], "dependency hash drift")

    builder = load_module("check_adaptive_width_builder", BUILDER)
    rational = load_module("check_adaptive_width_rational", builder.RATIONAL_BUILDER)
    two_sided = load_module("check_adaptive_width_two_sided", rational.TWO_SIDED_BUILDER)
    right_dependency = load_module("check_adaptive_width_right", two_sided.RIGHT_BUILDER)
    small = load_module("check_adaptive_width_small", right_dependency.SMALL_BUILDER)
    affine = load_module("check_adaptive_width_affine", right_dependency.AFFINE_BUILDER)
    box = load_module("check_adaptive_width_box", right_dependency.BOX_BUILDER)
    prior = load_module("check_adaptive_width_prior", right_dependency.PRIOR_BUILDER)
    altered = box.Settings(precision_bits=384, cutoff=11, abs_tol="1e-65")
    flint.ctx.prec = altered.precision_bits
    previous_rational = json.loads(builder.RATIONAL_RESULT.read_text(encoding="utf-8"))
    previous_two_sided = json.loads(rational.TWO_SIDED_RESULT.read_text(encoding="utf-8"))

    jacobian = artifact["local_jacobian_and_wall_anchor"]["wall_jacobian_at_x_2_over_5_s_1_over_3"]
    right_c = Fraction(builder.A, 2) + builder.S0 + 1
    right_raw_derivative = (-right_c + 2 * builder.FIRST_P) / (1 - 2 * builder.X0) ** 2
    require(Fraction(jacobian["right"]["d_a3_d_x"]) == -right_raw_derivative, "right a Jacobian drift")
    require(Fraction(jacobian["right"]["d_t3_d_x"]) == Fraction(25, 2), "right t Jacobian drift")
    left_c = Fraction(3, 2) + Fraction(builder.A, 2) + builder.S0
    left_raw_derivative = (3 * builder.FIRST_P - left_c) / (3 * builder.X0 - 1) ** 2
    require(Fraction(jacobian["left"]["d_a3_d_x"]) == -left_raw_derivative, "left a Jacobian drift")
    require(Fraction(jacobian["left"]["d_t3_d_x"]) == Fraction(-25, 2), "left t Jacobian drift")

    anchor = artifact["local_jacobian_and_wall_anchor"]["first_positive_floor_wall_anchor"]
    n = builder.OPEN_CHILD_N
    right_period, right_multiple = independent_wall_period(Fraction(anchor["right_a"]), n)
    left_period, left_multiple = independent_wall_period(Fraction(anchor["left_a"]), n)
    require(right_period == anchor["right_minimal_phase_period"], "right wall period drift")
    require(left_period == anchor["left_minimal_phase_period"], "left wall period drift")
    require(right_multiple == anchor["right_period_multiple_of_child_n"], "right wall multiple drift")
    require(left_multiple == anchor["left_period_multiple_of_child_n"], "left wall multiple drift")
    require(right_period > n and left_period > n, "wall phase unexpectedly compresses a full cycle")

    for side, search_name in (("right", "right_width_search"), ("left", "left_width_search")):
        search = artifact[search_name]
        pass_units = search["largest_passing_integer_units"]
        fail_units = search["first_failing_integer_units"]
        require(fail_units == pass_units + 1, f"{side} search boundary is not consecutive")
        production_rows = {row["width_units"]: row for row in search["rows"]}
        previous_radius = None
        for row in search["rows"]:
            radius = arb(row["complete_source_radius_upper"])
            if previous_radius is not None:
                require(radius > previous_radius, f"{side} production radii not increasing")
            previous_radius = radius
        for units, expected_pass in ((pass_units, True), (fail_units, False)):
            replay = builder.first_source_assembly(
                side,
                units,
                altered,
                box,
                small,
                affine,
                prior,
                rational,
                previous_rational,
                previous_two_sided,
            )
            require(replay["meets_source_radius_target"] is expected_pass, f"{side} altered target classification drift")
            require(
                parse_acb(replay["complete_source_box"]).overlaps(
                    parse_acb(production_rows[units]["complete_source_box"])
                ),
                f"{side} altered source box misses production",
            )
            if side == "left":
                p_box, a_box, tau_box = independent_left_child_boxes(units, builder, box)
                parent = production_rows[units]["parent"]
                require(same_box(parent["normalized_child_p_box"], p_box), "left child p box drift")
                require(same_box(parent["normalized_child_a_box"], a_box), "left child a box drift")
                require(same_box(parent["normalized_child_tau_box"], tau_box), "left child tau box drift")

    require(artifact["decision"]["finite_recursive_cover_built"] is False, "finite-cover overclaim")
    require(artifact["decision"]["rh_implication"] is False, "RH overclaim")
    require(artifact["proof_boundary"] in NOTE.read_text(encoding="utf-8"), "note proof boundary drift")
    print("independently replayed mirrored adaptive-width crossovers and wall-anchor periods", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
