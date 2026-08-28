#!/usr/bin/env python3
"""Independently replay the certified Mordell interval one-step gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
import flint
from flint import acb, arb
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_interval_one_step_current_gate"
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
    spec = importlib.util.spec_from_file_location("Mordell_interval_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "builder import failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_acb(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file() and BUILDER.is_file(), "interval artifact missing")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("kind") == STEM and artifact.get("passed") is True, "artifact identity drift")
    builder = load_builder()

    y, tau = sp.symbols("Y tau", positive=True, real=True)
    h_integral = sp.erfc(sp.sqrt(sp.pi * tau) * y) / (2 * sp.sqrt(tau))
    hp_integral = sp.exp(-sp.pi * tau * y**2) / (2 * sp.pi * tau)
    require(sp.simplify(sp.diff(h_integral, y) + sp.exp(-sp.pi * tau * y**2)) == 0, "h tail antiderivative failed")
    require(sp.simplify(sp.diff(hp_integral, y) + y * sp.exp(-sp.pi * tau * y**2)) == 0, "h-prime tail antiderivative failed")

    decisions = artifact["decision"]
    for key in (
        "uniform_central_strip_tail_bounds_proved",
        "pointwise_rational_Mordell_h_and_h_prime_encloser_built",
        "exact_unit_recurrence_propagated_with_balls",
        "transformed_current_enclosed_by_exact_rational_period",
        "source_current_enclosed_by_exact_rational_period",
        "two_full_one_step_current_enclosures_certified",
    ):
        require(decisions.get(key) is True, f"missing certified decision: {key}")
    require(decisions.get("nonrigorous_Gauss_Laguerre_route_used") is False, "quadrature route drift")
    for key in (
        "uniform_parameter_box_Mordell_evaluator_built",
        "small_tau_Euler_Maclaurin_branch_built",
        "recursive_interval_theta_evaluator_built",
        "physical_quadrature_completed",
        "non_A_bound_proved",
        "R_after_A_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    production_rows = artifact["full_interval_rows"]
    require(len(production_rows) == 2, "production row count drift")
    replay_settings = builder.Settings(
        precision_bits=REPLAY_BITS,
        cutoff=REPLAY_CUTOFF,
        abs_tol=REPLAY_ABS_TOL,
    )
    for stored in production_rows:
        x = Fraction(stored["x"])
        s = Fraction(stored["s"])
        replay = builder.compute_full_row(x, s, replay_settings)
        flint.ctx.prec = REPLAY_BITS
        production_complete = parse_acb(stored["complete_current"])
        production_source = parse_acb(stored["exact_period_source_current"])
        replay_complete = parse_acb(replay["complete_current"])
        replay_source = parse_acb(replay["exact_period_source_current"])
        require(production_complete.overlaps(production_source), "stored production source mismatch")
        require(replay_complete.overlaps(replay_source), "altered replay source mismatch")
        require(production_complete.overlaps(replay_complete), "production/replay current mismatch")
        require(production_source.overlaps(replay_source), "production/replay source mismatch")
        require(parse_acb(stored["complete_minus_source"]).contains(0), "stored difference excludes zero")
        require(parse_acb(replay["complete_minus_source"]).contains(0), "replay difference excludes zero")
        require(stored["transformed_period"] == replay["transformed_period"], "transformed period drift")
        require(stored["source_period"] == replay["source_period"], "source period drift")

    source = artifact["primary_source"]
    require(source["arxiv"] == "1306.4081v2", "primary-source drift")
    require(source["practical_Gauss_Laguerre_used"] is False, "source boundary drift")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")

    note = NOTE.read_text(encoding="utf-8")
    for token in ("(MI1)", "(MI2)", "(MI3)", "(MI4)", "coth(pi Y/sqrt(2))"):
        require(token in note, f"note token missing: {token}")
    normalized = " ".join(note.split()).lower()
    require("parameter-box recursion open" in normalized and "rh" in normalized, "proof boundary missing")
    print(
        "independently replayed certified Mordell tails and two full interval one-step currents",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
