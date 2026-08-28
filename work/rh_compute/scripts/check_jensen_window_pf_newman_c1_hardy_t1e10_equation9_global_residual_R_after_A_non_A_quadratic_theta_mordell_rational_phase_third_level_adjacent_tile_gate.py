#!/usr/bin/env python3
"""Independently replay the rational-phase adjacent-tile certificate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
from typing import Any


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


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def fractional_part(value: Fraction) -> Fraction:
    return value - floor_fraction(value)


def independent_period(a: Fraction, tau: Fraction) -> int:
    for period in range(1, 131_073):
        linear_increment = 2 * tau * period
        constant_increment = a * period + tau * period * period
        if linear_increment.denominator == 1 and constant_increment.denominator == 1:
            for smaller in range(1, period):
                smaller_linear = 2 * tau * smaller
                smaller_constant = a * smaller + tau * smaller * smaller
                require(
                    smaller_linear.denominator != 1 or smaller_constant.denominator != 1,
                    "reported phase period is not minimal",
                )
            return period
    raise RuntimeError("independent period search failed")


def direct_coefficients(
    n: int,
    degree: int,
    a0: Fraction,
    tau0: Fraction,
) -> dict[Fraction, Fraction]:
    period = independent_period(a0, tau0)
    result: dict[Fraction, Fraction] = {}
    for residue in range(period):
        coefficient = sum(Fraction(k**degree) for k in range(residue, n + 1, period))
        phase = fractional_part(a0 * residue + tau0 * residue * residue)
        result[phase] = result.get(phase, Fraction(0)) + coefficient
    return {phase: value for phase, value in result.items() if value}


def add_coefficients(
    left: dict[Fraction, Fraction],
    right: dict[Fraction, Fraction],
) -> dict[Fraction, Fraction]:
    result = dict(left)
    for phase, coefficient in right.items():
        result[phase] = result.get(phase, Fraction(0)) + coefficient
    return {phase: value for phase, value in result.items() if value}


def scale_coefficients(
    values: dict[Fraction, Fraction], scalar: Fraction
) -> dict[Fraction, Fraction]:
    return {phase: scalar * coefficient for phase, coefficient in values.items() if scalar * coefficient}


def coefficient_ball(values: dict[Fraction, Fraction], prior: Any) -> acb:
    result = acb(0)
    for phase, coefficient in sorted(values.items()):
        result += acb(prior.arb_rational(coefficient)) * prior.acb_cis_pi(2 * phase)
    return result


def parse_coefficients(rows: list[dict[str, str]]) -> dict[Fraction, Fraction]:
    return {
        Fraction(row["phase_turns"]): Fraction(row["coefficient"])
        for row in rows
    }


def independent_terminal(
    p_box: Any,
    a_box: Any,
    tau_box: Any,
    n: int,
    base_p: Fraction,
    base_a: Fraction,
    base_tau: Fraction,
    box: Any,
    small: Any,
    prior: Any,
) -> tuple[acb, dict[str, Any]]:
    f_coefficients = [
        direct_coefficients(n, degree, base_a, base_tau) for degree in range(4)
    ]
    m_coefficients = [
        add_coefficients(
            scale_coefficients(f_coefficients[degree], base_p),
            f_coefficients[degree + 1],
        )
        for degree in range(3)
    ]
    f_balls = [coefficient_ball(values, prior) for values in f_coefficients]
    m_balls = [coefficient_ball(values, prior) for values in m_coefficients]
    rho = box.arb_box(p_box) - prior.arb_rational(base_p)
    alpha = box.arb_box(a_box) - prior.arb_rational(base_a)
    beta = box.arb_box(tau_box) - prior.arb_rational(base_tau)
    b0 = m_balls[0] + acb(rho) * f_balls[0]
    b1 = m_balls[1] + acb(rho) * f_balls[1]
    b2 = m_balls[2] + acb(rho) * f_balls[2]
    first_order = b0 + 2 * acb.pi() * acb(0, 1) * (
        acb(alpha) * b1 + acb(beta) * b2
    )
    delta_a = max(abs(a_box.lo - base_a), abs(a_box.hi - base_a))
    delta_tau = max(abs(tau_box.lo - base_tau), abs(tau_box.hi - base_tau))
    masses = {degree: small.weighted_mass(p_box.hi, n, degree) for degree in range(2, 5)}
    remainder_mass = (
        delta_a * delta_a * masses[2]
        + 2 * delta_a * delta_tau * masses[3]
        + delta_tau * delta_tau * masses[4]
    )
    remainder_radius = 2 * arb.pi() ** 2 * prior.arb_rational(remainder_mass)
    enclosure = first_order + box.symmetric_complex_error(remainder_radius)
    return enclosure, {
        "F": f_coefficients,
        "M": m_coefficients,
        "remainder_radius": remainder_radius,
    }


def independently_derive_child_boxes(width: Fraction, box: Any, builder: Any) -> tuple[Any, Any, Any]:
    x_lo = builder.X0 + width
    x_hi = builder.X0 + 2 * width
    s_lo = builder.S0 - width
    s_hi = builder.S0 + width
    c_lo = Fraction(builder.A, 2) + s_lo + 1
    c_hi = Fraction(builder.A, 2) + s_hi + 1
    p_box = box.RationalBox(c_lo - 2 * builder.FIRST_P, c_hi - 2 * builder.FIRST_P)
    raw_a_lo = (builder.FIRST_P - c_hi * x_hi) / (1 - 2 * x_hi)
    raw_a_hi = (builder.FIRST_P - c_lo * x_lo) / (1 - 2 * x_lo)
    raw_tau_lo = -x_hi / (2 * (1 - 2 * x_hi))
    raw_tau_hi = -x_lo / (2 * (1 - 2 * x_lo))
    a_box = box.RationalBox(-raw_a_hi, -raw_a_lo)
    tau_box = box.RationalBox(-raw_tau_hi - 1, -raw_tau_lo - 1)
    return p_box, a_box, tau_box


def same_box(record: dict[str, str], value: Any) -> bool:
    return Fraction(record["lo"]) == value.lo and Fraction(record["hi"]) == value.hi


def main() -> int:
    require(BUILDER.is_file() and RESULT.is_file() and NOTE.is_file(), "missing production artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact not passed")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(Path(__file__)), "checker hash drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], "dependency hash drift")

    builder = load_module("check_rational_phase_builder", BUILDER)
    two_sided = load_module("check_rational_phase_two_sided", builder.TWO_SIDED_BUILDER)
    right = load_module("check_rational_phase_right", two_sided.RIGHT_BUILDER)
    small = load_module("check_rational_phase_small", right.SMALL_BUILDER)
    affine = load_module("check_rational_phase_affine", right.AFFINE_BUILDER)
    box = load_module("check_rational_phase_box", right.BOX_BUILDER)
    prior = load_module("check_rational_phase_prior", right.PRIOR_BUILDER)
    altered = box.Settings(precision_bits=384, cutoff=11, abs_tol="1e-65")
    flint.ctx.prec = altered.precision_bits

    reproduction_cases = (
        (Fraction(39_894), Fraction(-1, 2), Fraction(0), 1_240_710, 2),
        (Fraction(95_747, 6), Fraction(1, 3), Fraction(-1), 496_283, 3),
        (Fraction(47_873, 3), Fraction(2, 3), Fraction(-1), 496_283, 3),
    )
    for p0, a0, tau0, n, expected_period in reproduction_cases:
        require(independent_period(a0, tau0) == expected_period, "independent phase period drift")
        production_box, production_record = builder.rational_phase_taylor_box(
            box.RationalBox(p0, p0),
            box.RationalBox(a0, a0),
            box.RationalBox(tau0, tau0),
            n,
            p0,
            a0,
            tau0,
            box,
            small,
            affine,
            prior,
        )
        independent_box, independent_record = independent_terminal(
            box.RationalBox(p0, p0),
            box.RationalBox(a0, a0),
            box.RationalBox(tau0, tau0),
            n,
            p0,
            a0,
            tau0,
            box,
            small,
            prior,
        )
        for degree in range(4):
            require(
                parse_coefficients(production_record["F_phase_coefficients"][f"F{degree}"])
                == independent_record["F"][degree],
                "exact F residue coefficients drift",
            )
        for degree in range(3):
            require(
                parse_coefficients(production_record["M_phase_coefficients"][f"M{degree}"])
                == independent_record["M"][degree],
                "exact M residue coefficients drift",
            )
        require(production_box.overlaps(independent_box), "generic point replay mismatch")

    width = builder.ADJACENT_WIDTH
    p_child, a_child, tau_child = independently_derive_child_boxes(width, box, builder)
    production_child = artifact["certified_adjacent_third_level_tile"]["complete_normalized_parent"]
    require(same_box(production_child["normalized_child_p_box"], p_child), "adjacent p box drift")
    require(same_box(production_child["normalized_child_a_box"], a_child), "adjacent a box drift")
    require(same_box(production_child["normalized_child_tau_box"], tau_child), "adjacent tau box drift")
    independent_child, independent_record = independent_terminal(
        p_child,
        a_child,
        tau_child,
        builder.OPEN_CHILD_N,
        builder.RIGHT_BASE_P,
        builder.RIGHT_BASE_A,
        Fraction(0),
        box,
        small,
        prior,
    )
    production_terminal = production_child["rational_phase_terminal"]
    for degree in range(4):
        require(
            parse_coefficients(production_terminal["F_phase_coefficients"][f"F{degree}"])
            == independent_record["F"][degree],
            "adjacent F residue coefficients drift",
        )
    require(
        independent_child.overlaps(acb(arb(production_terminal["rational_phase_enclosure"]["real_ball"]), arb(production_terminal["rational_phase_enclosure"]["imag_ball"]))),
        "independent adjacent terminal misses production",
    )

    n = builder.OPEN_CHILD_N
    right_wall = Fraction(2 * n + 1, 5 * n + 2)
    left_wall = Fraction(2 * n + 1, 5 * n + 3)
    inventory = artifact["recursive_branch_wall_inventory"]
    require(Fraction(inventory["right_first_positive_floor_wall"]) == right_wall, "right wall formula drift")
    require(Fraction(inventory["left_first_positive_floor_wall"]) == left_wall, "left wall formula drift")
    require(builder.X0 + 2 * width < right_wall, "adjacent tile crossed first right floor wall")
    require(
        Fraction(inventory["adjacent_width_units_to_first_right_floor_wall"])
        == (right_wall - builder.X0) / width,
        "uniform tile-count scale drift",
    )

    previous = json.loads(builder.TWO_SIDED_RESULT.read_text(encoding="utf-8"))
    replay = builder.adjacent_source_tile(
        width, altered, box, small, affine, prior, previous
    )
    replay_source = acb(
        arb(replay["complete_source_box"]["real_ball"]),
        arb(replay["complete_source_box"]["imag_ball"]),
    )
    production_source = acb(
        arb(artifact["certified_adjacent_third_level_tile"]["complete_source_box"]["real_ball"]),
        arb(artifact["certified_adjacent_third_level_tile"]["complete_source_box"]["imag_ball"]),
    )
    require(replay_source.overlaps(production_source), "altered-settings source replay mismatch")
    require(arb(replay["complete_radius_absolute_upper"]) < arb("0.02"), "replay radius target failed")
    require(artifact["decision"]["finite_recursive_cover_built"] is False, "proof boundary overclaim")
    require(artifact["decision"]["rh_implication"] is False, "RH overclaim")
    require(artifact["proof_boundary"] in NOTE.read_text(encoding="utf-8"), "note proof boundary drift")
    print("independently replayed the rational-phase terminal and adjacent third-level source tile", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
