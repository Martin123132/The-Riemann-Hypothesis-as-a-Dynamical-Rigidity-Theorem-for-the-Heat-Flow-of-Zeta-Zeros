#!/usr/bin/env python3
"""Independently replay the two-sided affine Mordell wall box."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
from flint import acb, arb
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_two_sided_wall_period_three_recursive_box_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
REPLAY_BITS = 384
REPLAY_CUTOFF = 11
REPLAY_ABS_TOL = "1e-65"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_builder():
    spec = importlib.util.spec_from_file_location("mordell_two_sided_wall_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "builder import failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_acb(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file() and BUILDER.is_file(), "two-sided artifact missing")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("kind") == STEM and artifact.get("passed") is True, "artifact identity drift")
    builder = load_builder()
    right = builder.load_module("mordell_two_sided_right_replay_dependency", builder.RIGHT_BUILDER)
    small = builder.load_module("mordell_two_sided_small_replay_dependency", right.SMALL_BUILDER)
    affine = builder.load_module("mordell_two_sided_affine_replay_dependency", right.AFFINE_BUILDER)
    box = builder.load_module("mordell_two_sided_box_replay_dependency", right.BOX_BUILDER)
    prior = builder.load_module("mordell_two_sided_point_replay_dependency", right.PRIOR_BUILDER)
    settings = box.Settings(REPLAY_BITS, REPLAY_CUTOFF, REPLAY_ABS_TOL)

    omega = prior.acb_cis_pi(Fraction(2, 3))
    omega2 = omega * omega
    for n in range(0, 31):
        for degree in range(4):
            right_pair = right.period_three_power_pair(n, degree, small)
            left_pair = builder.left_period_three_power_pair(n, degree, right, small)
            require(left_pair == builder.conjugate_pair(right_pair), "conjugate pair drift")
            encoded = right.pair_ball(left_pair, prior)
            direct = sum((acb(k**degree) * omega2**k for k in range(n + 1)), acb(0))
            require(encoded.overlaps(direct), "left period-three pair evaluation drift")
        for degree in range(3):
            pair = builder.left_period_three_weighted_pair(
                builder.LEFT_BASE_P, n, degree, right, small
            )
            encoded = right.pair_ball(pair, prior)
            direct = sum(
                (
                    acb(prior.arb_rational(builder.LEFT_BASE_P + k))
                    * acb(k**degree)
                    * omega2**k
                    for k in range(n + 1)
                ),
                acb(0),
            )
            require(encoded.overlaps(direct), "left weighted pair evaluation drift")

    a_r, tau_r = sp.symbols("a_r tau_r", real=True)
    a_left = sp.Rational(1, 2) - a_r
    tau_left = sp.Rational(1, 2) - tau_r
    require(a_left.subs(a_r, sp.Rational(1, 6)) == sp.Rational(1, 3), "left a wall map drift")
    require(tau_left.subs(tau_r, sp.Rational(1, 4)) == sp.Rational(1, 4), "left tau wall map drift")
    x = builder.X0 - builder.TWO_SIDED_WIDTH
    tau_open = Fraction(3, 2) - Fraction(1, 2) / x
    require(496_283 < 2 * builder.FIRST_N * tau_open < 496_284, "left open floor drift")
    require(2 * builder.FIRST_N * Fraction(1, 4) == 496_284, "wall floor drift")

    replay = builder.two_sided_source_box(
        builder.TWO_SIDED_WIDTH, settings, box, right, small, affine, prior
    )
    stored = artifact["certified_two_sided_wall_box"]
    require(stored["left_open_child_floor"] == replay["left_open_child_floor"] == 496_283, "left floor drift")
    require(stored["wall_child_floor"] == replay["wall_child_floor"] == 496_284, "wall floor drift")
    require(
        stored["right_dependency_projection"]["open_branch_child_floor"]
        == replay["right_dependency_projection"]["open_branch_child_floor"]
        == 496_283,
        "right projection floor drift",
    )

    stored_complete = parse_acb(stored["complete_source_box"])
    replay_complete = parse_acb(replay["complete_source_box"])
    stored_source = parse_acb(stored["exact_center_source"])
    replay_source = parse_acb(replay["exact_center_source"])
    require(stored_complete.contains(stored_source), "stored complete box misses source")
    require(replay_complete.contains(replay_source), "replay complete box misses source")
    require(stored_complete.overlaps(replay_complete), "altered-precision complete box mismatch")
    require(stored_source.overlaps(replay_source), "altered-precision source mismatch")

    for branch_name in ("left_open_branch", "left_wall_branch"):
        stored_branch = stored[branch_name]
        replay_branch = replay[branch_name]
        stored_parent = parse_acb(stored_branch["complete_original_parent"])
        replay_parent = parse_acb(replay_branch["complete_original_parent"])
        stored_child = parse_acb(stored_branch["child_Taylor"]["child_enclosure"])
        replay_child = parse_acb(replay_branch["child_Taylor"]["child_enclosure"])
        stored_base = parse_acb(stored_branch["child_Taylor"]["exact_base_current"])
        replay_base = parse_acb(replay_branch["child_Taylor"]["exact_base_current"])
        require(stored_parent.overlaps(replay_parent), f"{branch_name} parent replay mismatch")
        require(stored_child.contains(stored_base), f"{branch_name} stored child misses base")
        require(replay_child.contains(replay_base), f"{branch_name} replay child misses base")
        require(stored_child.overlaps(replay_child), f"{branch_name} child replay mismatch")
        require(stored_branch["child_Taylor"]["base_period"] == 3, f"{branch_name} period drift")

    left_open = parse_acb(stored["left_open_branch"]["complete_original_parent"])
    left_wall = parse_acb(stored["left_wall_branch"]["complete_original_parent"])
    right_wall = parse_acb(stored["right_dependency_projection"]["right_wall_original_parent"])
    right_open = parse_acb(stored["right_dependency_projection"]["right_open_original_parent"])
    require(left_open.overlaps(left_wall), "stored left closure mismatch")
    require(left_wall.overlaps(right_wall), "stored wall orientation mismatch")
    require(right_wall.overlaps(right_open), "stored right closure mismatch")
    transformed = parse_acb(stored["transformed_two_sided_hull"])
    center_transformed = parse_acb(stored["exact_center_transformed_current"])
    require(transformed.contains(center_transformed), "stored transformed hull misses center")

    exact_left_parent, left_parent_period = affine.weighted_period_ball(
        builder.FIRST_P, Fraction(1, 3), Fraction(1, 4), builder.FIRST_N, prior
    )
    exact_right_parent, right_parent_period = affine.weighted_period_ball(
        builder.FIRST_P, Fraction(1, 6), Fraction(1, 4), builder.FIRST_N, prior
    )
    require(left_parent_period > 0 and right_parent_period > 0, "wall parent period missing")
    require(left_wall.contains(exact_left_parent), "left wall misses exact left parent")
    require(right_wall.contains(exact_right_parent.conjugate()), "right wall misses conjugated exact parent")
    require(exact_left_parent.overlaps(exact_right_parent.conjugate()), "exact wall orientation identity failed")

    replay_audit = builder.left_scale_audit(settings, box, right, small, affine, prior)
    require(len(replay_audit) == len(artifact["left_scale_audit"]) == 7, "scale-audit count drift")
    for stored_row, replay_row in zip(artifact["left_scale_audit"], replay_audit):
        require(stored_row["half_width"] == replay_row["half_width"], "scale width drift")
        for key in ("delta_p_max", "delta_a_max", "delta_tau_max"):
            require(stored_row[key] == replay_row[key], f"exact scale drift: {key}")
        stored_bound = arb(stored_row["two_level_child_variation_upper"])
        replay_bound = arb(replay_row["two_level_child_variation_upper"])
        require(bool(stored_bound < 2 * replay_bound and replay_bound < 2 * stored_bound), "scale replay drift")

    k, alpha, beta = sp.symbols("k alpha beta", real=True)
    phase_increment = 2 * sp.pi * (alpha * k + beta * k**2)
    require(
        sp.expand(phase_increment**2 / 2 - 2 * sp.pi**2 * (alpha * k + beta * k**2) ** 2) == 0,
        "Taylor remainder coefficient drift",
    )

    decisions = artifact["decision"]
    for key in (
        "left_parity_half_shift_derived",
        "double_conjugation_orientation_proved",
        "left_period_three_child_Taylor_boxes_built",
        "left_endpoint_complete_recursive_box_built",
        "left_and_right_wall_orientations_joined",
        "two_sided_interior_source_box_certified",
    ):
        require(decisions.get(key) is True, f"missing positive decision: {key}")
    for key in (
        "recursive_tiling_beyond_local_box_built",
        "physical_quadrature_completed",
        "non_A_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")

    note = NOTE.read_text(encoding="utf-8")
    for token in ("(TW1)", "(TW2)", "(TW3)", "(TW4)", "(TW5)", "e(u)=exp(2*pi*i*u)"):
        require(token in note, f"note token missing: {token}")
    normalized = " ".join(note.split()).lower()
    require("two-sided local interval lemma" in normalized, "two-sided status missing")
    require("no recursive tiling beyond this local box" in normalized, "proof boundary missing")
    print(
        "independently replayed the two-sided Mordell wall and period-three recursive source box",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
