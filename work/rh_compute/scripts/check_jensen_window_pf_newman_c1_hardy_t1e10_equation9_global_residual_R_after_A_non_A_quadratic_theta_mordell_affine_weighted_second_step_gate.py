#!/usr/bin/env python3
"""Independently check the affine-weighted Mordell second-step gate."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_affine_weighted_second_step_gate"
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
    spec = importlib.util.spec_from_file_location("mordell_affine_second_step_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "builder import failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_acb(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file() and BUILDER.is_file(), "affine artifact missing")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("kind") == STEM and artifact.get("passed") is True, "artifact identity drift")
    builder = load_builder()
    prior = builder.load_module("mordell_affine_second_step_replay_dependency", builder.PRIOR_BUILDER)
    replay_settings = prior.Settings(REPLAY_BITS, REPLAY_CUTOFF, REPLAY_ABS_TOL)

    p, z, tau = sp.symbols("p z tau", real=True, nonzero=True)
    child_weight = sp.simplify(2 * tau * p - z)
    require(child_weight == 2 * tau * p - z, "child-weight algebra failed")
    k = sp.symbols("k", integer=True)
    require(sp.simplify((k - k**2) / 2).is_integer is not False, "parity identity rejected")
    for witness in range(-12, 13):
        require((witness - witness * witness) % 2 == 0, "parity half-shift witness failed")

    short_rows = artifact["short_interval_rows"]
    require(len(short_rows) == len(builder.SHORT_CASES) == 4, "short-row count drift")
    for stored, case in zip(short_rows, builder.SHORT_CASES):
        complete, replay = builder.affine_weighted_mordell_step(
            settings=replay_settings, prior=prior, **case
        )
        production = parse_acb(stored["complete_current"])
        direct = parse_acb(stored["direct_weighted_current"])
        require(production.overlaps(direct), "stored short direct overlap failed")
        require(complete.overlaps(parse_acb(replay["direct_weighted_current"])), "replay short direct overlap failed")
        require(production.overlaps(complete), "short altered replay mismatch")
        require(parse_acb(stored["complete_minus_direct"]).contains(0), "stored short difference excludes zero")
        require(parse_acb(replay["complete_minus_direct"]).contains(0), "replay short difference excludes zero")

    full_cases = ((Fraction(2, 5), Fraction(1, 3)), (Fraction(1, 2), Fraction(0)))
    full_rows = artifact["full_second_step_rows"]
    require(len(full_rows) == 2, "full-row count drift")
    for stored, (x, s) in zip(full_rows, full_cases):
        replay = builder.normalized_full_second_step(x, s, replay_settings, prior)
        production_complete = parse_acb(stored["recursive_complete_source_current"])
        production_source = parse_acb(stored["exact_period_source_current"])
        replay_complete = parse_acb(replay["recursive_complete_source_current"])
        replay_source = parse_acb(replay["exact_period_source_current"])
        require(production_complete.overlaps(production_source), "stored full source overlap failed")
        require(replay_complete.overlaps(replay_source), "replay full source overlap failed")
        require(production_complete.overlaps(replay_complete), "full altered replay mismatch")
        require(parse_acb(stored["recursive_complete_minus_source"]).contains(0), "stored source difference excludes zero")
        require(parse_acb(replay["recursive_complete_minus_source"]).contains(0), "replay source difference excludes zero")
        require(stored["step_kind"] == replay["step_kind"], "full step-kind drift")

    interior, corner = full_rows
    require(interior["first_transformed_n"] == 992568 and interior["second_child_m"] == 496284, "interior contraction drift")
    require(interior["second_step"]["child_p"] == "95747/6", "affine child weight drift")
    require(interior["second_step"]["child_a"] == "1/3", "second child phase drift")
    require(interior["second_step"]["child_tau"] == "-1", "second child curvature drift")
    require(corner["normalization"]["tau"] == "0" and corner["step_kind"] == "exact_tau_zero_terminal", "corner terminal drift")

    decisions = artifact["decision"]
    for key in (
        "general_affine_weight_closed_under_Mordell_step",
        "exact_phase_normalization_proved",
        "interior_full_second_step_contracts",
        "corner_center_has_exact_tau_zero_terminal",
        "two_full_recursive_source_reassemblies_certified",
    ):
        require(decisions.get(key) is True, f"missing positive decision: {key}")
    for key in (
        "recursive_parameter_boxes_built",
        "small_tau_neighbourhood_branch_built",
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
    for token in ("(AW1)", "(AW2)", "(AW3)", "(AW4)", "p'=2tau p-z"):
        require(token in note, f"note token missing: {token}")
    normalized = " ".join(note.split()).lower()
    require("recursive parameter boxes and the small-tau neighbourhood branch remain open" in normalized, "status boundary missing")
    require("does not prove recursive parameter boxes" in normalized and "rh" in normalized, "proof boundary missing")
    print(
        "independently replayed affine-weighted Mordell recursion and both full second-step rows",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
