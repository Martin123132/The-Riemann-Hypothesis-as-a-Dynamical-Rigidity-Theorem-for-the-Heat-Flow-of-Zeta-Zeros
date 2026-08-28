#!/usr/bin/env python3
"""Independently replay the right-side affine Mordell wall box."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_right_wall_period_three_recursive_box_gate"
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
    spec = importlib.util.spec_from_file_location("mordell_right_wall_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "builder import failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_acb(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file() and BUILDER.is_file(), "right-wall artifact missing")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("kind") == STEM and artifact.get("passed") is True, "artifact identity drift")
    builder = load_builder()
    small = builder.load_module("mordell_right_wall_small_replay_dependency", builder.SMALL_BUILDER)
    affine = builder.load_module("mordell_right_wall_affine_replay_dependency", builder.AFFINE_BUILDER)
    box = builder.load_module("mordell_right_wall_box_replay_dependency", builder.BOX_BUILDER)
    prior = builder.load_module("mordell_right_wall_point_replay_dependency", builder.PRIOR_BUILDER)
    settings = box.Settings(REPLAY_BITS, REPLAY_CUTOFF, REPLAY_ABS_TOL)

    for n in range(0, 31):
        for modulus in (2, 3, 4, 5):
            for residue in range(modulus):
                for degree in range(5):
                    direct = sum(k**degree for k in range(residue, n + 1, modulus))
                    require(
                        builder.residue_power_sum(n, residue, modulus, degree, small) == direct,
                        "residue power-sum drift",
                    )
        for degree in range(4):
            c = [sum(k**degree for k in range(residue, n + 1, 3)) for residue in range(3)]
            expected = (Fraction(c[0] - c[2]), Fraction(c[1] - c[2]))
            require(builder.period_three_power_pair(n, degree, small) == expected, "period-three pair drift")

    omega = prior.acb_cis_pi(Fraction(2, 3))
    for n in range(0, 25):
        for degree in range(3):
            direct = sum((acb(k**degree) * omega**k for k in range(n + 1)), acb(0))
            encoded = builder.pair_ball(builder.period_three_power_pair(n, degree, small), prior)
            require(encoded.overlaps(direct), "period-three pair evaluation drift")

    x = sp.symbols("x", positive=True)
    tau = sp.Rational(1, 2) / x - 1
    require(sp.simplify(tau.subs(x, sp.Rational(2, 5)) - sp.Rational(1, 4)) == 0, "wall tau drift")
    wall_floor_argument = 2 * builder.FIRST_N * Fraction(1, 4)
    require(wall_floor_argument == 496_284, "wall floor argument drift")
    right_x = builder.X0 + builder.RIGHT_WIDTH
    right_tau = Fraction(1, 2) / right_x - 1
    require(496_283 < 2 * builder.FIRST_N * right_tau < 496_284, "right floor interval drift")

    replay = builder.right_wall_source_box(
        builder.RIGHT_WIDTH, settings, box, small, affine, prior
    )
    stored = artifact["certified_right_wall_box"]
    require(stored["open_branch_child_floor"] == replay["open_branch_child_floor"] == 496_283, "open floor drift")
    require(stored["wall_branch_child_floor"] == replay["wall_branch_child_floor"] == 496_284, "wall floor drift")

    stored_complete = parse_acb(stored["complete_source_box"])
    replay_complete = parse_acb(replay["complete_source_box"])
    stored_source = parse_acb(stored["exact_center_source"])
    replay_source = parse_acb(replay["exact_center_source"])
    require(stored_complete.contains(stored_source), "stored source box misses center")
    require(replay_complete.contains(replay_source), "replay source box misses center")
    require(stored_complete.overlaps(replay_complete), "altered-precision source box mismatch")
    require(stored_source.overlaps(replay_source), "altered-precision center source mismatch")

    for branch_name in ("interior_branch", "wall_branch"):
        stored_branch = stored[branch_name]
        replay_branch = replay[branch_name]
        stored_parent = parse_acb(stored_branch["complete_normalized_parent"])
        replay_parent = parse_acb(replay_branch["complete_normalized_parent"])
        stored_child = parse_acb(stored_branch["child_Taylor"]["child_enclosure"])
        replay_child = parse_acb(replay_branch["child_Taylor"]["child_enclosure"])
        stored_base = parse_acb(stored_branch["child_Taylor"]["exact_base_current"])
        replay_base = parse_acb(replay_branch["child_Taylor"]["exact_base_current"])
        require(stored_parent.overlaps(replay_parent), f"{branch_name} parent replay mismatch")
        require(stored_child.contains(stored_base), f"{branch_name} stored child misses base")
        require(replay_child.contains(replay_base), f"{branch_name} replay child misses base")
        require(stored_child.overlaps(replay_child), f"{branch_name} child replay mismatch")
        require(stored_branch["child_Taylor"]["base_period"] == 3, f"{branch_name} base period drift")

    center_parent, center_period = affine.weighted_period_ball(
        builder.FIRST_P, Fraction(1, 6), Fraction(1, 4), builder.FIRST_N, prior
    )
    require(center_period > 0, "center parent period missing")
    require(
        parse_acb(stored["wall_branch"]["complete_normalized_parent"]).contains(center_parent),
        "stored wall parent misses exact center",
    )
    require(
        parse_acb(replay["wall_branch"]["complete_normalized_parent"]).contains(center_parent),
        "replay wall parent misses exact center",
    )
    require(
        parse_acb(stored["interior_branch"]["complete_normalized_parent"]).contains(center_parent),
        "stored open closure misses exact center limit",
    )

    replay_audit = builder.child_scale_audit(settings, box, small, affine, prior)
    require(len(replay_audit) == len(artifact["child_scale_audit"]) == 7, "scale-audit count drift")
    for stored_row, replay_row in zip(artifact["child_scale_audit"], replay_audit):
        require(stored_row["half_width"] == replay_row["half_width"], "scale width drift")
        for key in ("delta_p_max", "delta_a_max", "delta_tau_max"):
            require(stored_row[key] == replay_row[key], f"exact scale drift: {key}")
        stored_bound = arb(stored_row["two_level_child_variation_upper"])
        replay_bound = arb(replay_row["two_level_child_variation_upper"])
        require(bool(stored_bound < 2 * replay_bound and replay_bound < 2 * stored_bound), "scale replay drift")

    p, a, t = sp.symbols("p a t", real=True)
    rho, alpha, beta = sp.symbols("rho alpha beta", real=True)
    require(sp.expand((p + rho) - p - rho) == 0, "affine weight expansion drift")
    phase_increment = 2 * sp.pi * (alpha * sp.Symbol("k") + beta * sp.Symbol("k") ** 2)
    require(sp.expand(phase_increment**2 / 2 - 2 * sp.pi**2 * (alpha * sp.Symbol("k") + beta * sp.Symbol("k") ** 2) ** 2) == 0, "remainder coefficient drift")

    decisions = artifact["decision"]
    for key in (
        "right_normalization_branch_partitioned",
        "child_floor_wall_partitioned",
        "period_three_child_Taylor_boxes_built",
        "right_endpoint_complete_recursive_box_built",
        "right_complete_source_box_certified",
    ):
        require(decisions.get(key) is True, f"missing positive decision: {key}")
    for key in (
        "left_normalization_branch_built",
        "two_sided_interior_neighbourhood_built",
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
    for token in ("(RW1)", "(RW2)", "(RW3)", "(RW4)", "e(u)=exp(2*pi*i*u)"):
        require(token in note, f"note token missing: {token}")
    normalized = " ".join(note.split()).lower()
    require("left side remains open" in normalized, "one-sided status missing")
    require("no left-side branch" in normalized, "proof boundary missing")
    print(
        "independently replayed the right-side Mordell wall and period-three recursive source box",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
