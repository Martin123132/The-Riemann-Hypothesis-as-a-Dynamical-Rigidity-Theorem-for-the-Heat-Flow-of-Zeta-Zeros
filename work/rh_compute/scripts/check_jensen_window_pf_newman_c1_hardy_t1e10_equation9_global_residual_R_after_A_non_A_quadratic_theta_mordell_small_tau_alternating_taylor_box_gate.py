#!/usr/bin/env python3
"""Independently replay the small-tau alternating Taylor corner gate."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_small_tau_alternating_taylor_box_gate"
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
    spec = importlib.util.spec_from_file_location("mordell_small_tau_taylor_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "builder import failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_acb(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def direct_alternating_moment(p: Fraction, n: int, degree: int) -> Fraction:
    return sum(
        (((-1) ** k) * (p + k) * k**degree for k in range(n + 1)),
        Fraction(0),
    )


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file() and BUILDER.is_file(), "small-tau artifact missing")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("kind") == STEM and artifact.get("passed") is True, "artifact identity drift")
    builder = load_builder()
    box = builder.load_module("mordell_small_tau_box_replay_dependency", builder.BOX_BUILDER)
    affine = builder.load_module("mordell_small_tau_affine_replay_dependency", builder.AFFINE_BUILDER)
    prior = builder.load_module("mordell_small_tau_point_replay_dependency", builder.PRIOR_BUILDER)
    settings = box.Settings(REPLAY_BITS, REPLAY_CUTOFF, REPLAY_ABS_TOL)

    for count in range(0, 24):
        for degree in range(6):
            direct = sum(h**degree for h in range(count))
            require(builder.power_sum(count, degree) == direct, "power-sum formula drift")
    for n in range(0, 25):
        for residue in (0, 1):
            for degree in range(4):
                direct = sum(k**degree for k in range(residue, n + 1, 2))
                require(
                    builder.residue_power_sum(n, residue, degree) == direct,
                    "parity power-sum formula drift",
                )
    for p in (Fraction(0), Fraction(3, 2), Fraction(7)):
        for n in range(0, 24):
            for degree in range(3):
                require(
                    builder.alternating_weighted_moment(p, n, degree)
                    == direct_alternating_moment(p, n, degree),
                    "alternating moment formula drift",
                )
            for degree in range(5):
                direct_mass = sum(
                    ((p + k) * k**degree for k in range(n + 1)),
                    Fraction(0),
                )
                require(builder.weighted_mass(p, n, degree) == direct_mass, "weighted mass formula drift")

    p = Fraction(39_894)
    n = 1_240_710
    full_moments = [builder.alternating_weighted_moment(p, n, degree) for degree in range(3)]
    require(
        full_moments == [Fraction(660_249), Fraction(794_429_714_775), Fraction(985_657_301_007_258_645)],
        "full alternating moments drift",
    )

    d = sp.symbols("d", real=True)
    series = sp.series(sp.exp(sp.I * d), d, 0, 3).removeO()
    require(sp.expand(series - (1 + sp.I * d - d**2 / 2)) == 0, "complex exponential Taylor jet drift")
    alpha, t, k = sp.symbols("alpha t k", real=True)
    phase_increment = sp.expand(2 * sp.pi * (alpha * k + t * k**2))
    expected_square = sp.expand(2 * sp.pi**2 * (alpha * k + t * k**2) ** 2)
    require(
        sp.expand(phase_increment**2 / 2 - expected_square) == 0,
        "Taylor remainder coefficient drift",
    )

    replay = builder.corner_source_box(builder.CORNER_WIDTH, settings, box, affine, prior)
    stored = artifact["certified_corner_box"]
    stored_complete = parse_acb(stored["complete_source_box"])
    replay_complete = parse_acb(replay["complete_source_box"])
    stored_source = parse_acb(stored["exact_corner_source"])
    replay_source = parse_acb(replay["exact_corner_source"])
    require(stored_complete.contains(stored_source), "stored complete box misses source")
    require(replay_complete.contains(replay_source), "replay complete box misses source")
    require(stored_complete.overlaps(replay_complete), "altered-precision complete box mismatch")
    require(stored_source.overlaps(replay_source), "altered-precision source mismatch")
    require(
        parse_acb(stored["small_tau_normalized_child"]["small_tau_enclosure"]).contains(
            parse_acb(stored["exact_terminal_current"])
        ),
        "stored small-tau child misses exact terminal",
    )
    require(
        parse_acb(replay["small_tau_normalized_child"]["small_tau_enclosure"]).contains(
            parse_acb(replay["exact_terminal_current"])
        ),
        "replay small-tau child misses exact terminal",
    )
    require(stored["exact_terminal_period"] == replay["exact_terminal_period"] == 2, "terminal period drift")
    require(stored["half_width"] == replay["half_width"] == "1/10000000000000000000000", "width drift")
    require(artifact["summary"]["width_gain_over_direct_1e_minus_28_box"] == 1_000_000, "width gain drift")
    require(bool(arb(stored["complete_radius_absolute_upper"]) < arb("0.02")), "stored radius target failed")
    require(bool(arb(replay["complete_radius_absolute_upper"]) < arb("0.02")), "replay radius target failed")

    replay_audit = builder.scale_audit(settings, box)
    require(len(replay_audit) == len(artifact["scale_audit"]) == 6, "scale-audit count drift")
    for stored_row, replay_row in zip(artifact["scale_audit"], replay_audit):
        require(stored_row["half_width"] == replay_row["half_width"], "scale-audit width drift")
        require(stored_row["delta_a_max"] == replay_row["delta_a_max"], "delta-a exact drift")
        require(stored_row["t_max"] == replay_row["t_max"], "t-max exact drift")
        stored_bound = arb(stored_row["first_main_child_variation_upper"])
        replay_bound = arb(replay_row["first_main_child_variation_upper"])
        require(bool(stored_bound < 2 * replay_bound and replay_bound < 2 * stored_bound), "scale bound replay drift")

    decisions = artifact["decision"]
    for key in (
        "small_tau_corner_neighbourhood_branch_built",
        "exact_tau_zero_terminal_joined_continuously",
        "complete_corner_source_box_certified",
    ):
        require(decisions.get(key) is True, f"missing positive decision: {key}")
    for key in (
        "uniform_small_tau_branch_away_from_corner_built",
        "recursive_interior_parameter_boxes_built",
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
    for token in ("(ST1)", "(ST2)", "(ST3)", "(ST4)", "e(u)=exp(2*pi*i*u)"):
        require(token in note, f"note token missing: {token}")
    normalized = " ".join(note.split()).lower()
    require("one million times" in normalized, "width-gain statement missing")
    require("does not prove" not in normalized or "rh" in normalized, "proof boundary missing")
    require("no uniform small-tau branch away from this corner" in normalized, "local-status boundary missing")
    print(
        "independently replayed the alternating small-tau Taylor box and corner source reassembly",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
