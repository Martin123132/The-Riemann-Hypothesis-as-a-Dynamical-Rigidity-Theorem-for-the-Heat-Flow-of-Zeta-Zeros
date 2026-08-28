#!/usr/bin/env python3
"""Independently replay the Mordell parameter-box and phase-scale gate."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_parameter_box_phase_scale_gate"
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
    spec = importlib.util.spec_from_file_location("mordell_parameter_box_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "builder import failed")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_acb(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def parse_box(builder, record: dict[str, str]):
    return builder.RationalBox(Fraction(record["lo"]), Fraction(record["hi"]))


def exact_masses(r: int, m: int) -> tuple[int, int, int]:
    direct0 = sum(r + k for k in range(m + 1))
    direct1 = sum((r + k) * k for k in range(m + 1))
    direct2 = sum((r + k) * k * k for k in range(m + 1))
    return direct0, direct1, direct2


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file() and BUILDER.is_file(), "parameter-box artifact missing")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("kind") == STEM and artifact.get("passed") is True, "artifact identity drift")
    builder = load_builder()
    prior = builder.load_module("mordell_parameter_box_replay_dependency", builder.PRIOR_BUILDER)
    replay_settings = builder.Settings(REPLAY_BITS, REPLAY_CUTOFF, REPLAY_ABS_TOL)

    y, tau = sp.symbols("Y tau", positive=True, real=True)
    h_tail = sp.erfc(sp.sqrt(sp.pi * tau) * y) / (2 * sp.sqrt(tau))
    hz_tail = sp.exp(-sp.pi * tau * y**2) / (2 * sp.pi * tau)
    require(sp.simplify(sp.diff(h_tail, y) + sp.exp(-sp.pi * tau * y**2)) == 0, "h tail primitive failed")
    require(sp.simplify(sp.diff(hz_tail, y) + y * sp.exp(-sp.pi * tau * y**2)) == 0, "h_z tail primitive failed")

    target_rows = artifact["target_parameter_boxes"]
    require(len(target_rows) == 3, "target box count drift")
    for stored in target_rows:
        enclosure = stored["enclosure"]
        target = parse_box(builder, enclosure["target_box"])
        tau_box = parse_box(builder, enclosure["tau_box"])
        replay_h, replay_hz, replay_record = builder.mordell_target_parameter_box(
            target, tau_box, enclosure["signed_tau"], replay_settings, prior
        )
        production_h = parse_acb(enclosure["h"])
        production_hz = parse_acb(enclosure["h_z"])
        require(production_h.overlaps(replay_h), f"altered h box replay mismatch: {stored['name']}")
        require(production_hz.overlaps(replay_hz), f"altered h_z box replay mismatch: {stored['name']}")
        require(replay_record["branch_count"] == enclosure["branch_count"], "branch-count replay drift")

        z_halves = (
            builder.RationalBox(target.lo, target.mid),
            builder.RationalBox(target.mid, target.hi),
        )
        tau_halves = (
            builder.RationalBox(tau_box.lo, tau_box.mid),
            builder.RationalBox(tau_box.mid, tau_box.hi),
        )
        for z_sub in z_halves:
            for tau_sub in tau_halves:
                sub_h, sub_hz, _ = builder.mordell_target_parameter_box(
                    z_sub, tau_sub, enclosure["signed_tau"], replay_settings, prior
                )
                require(production_h.overlaps(sub_h), "subdivision h replay mismatch")
                require(production_hz.overlaps(sub_hz), "subdivision h_z replay mismatch")

    decisions = artifact["decision"]
    for key in (
        "uniform_central_Mordell_parameter_boxes_built",
        "half_integer_recurrence_wall_split_exactly",
        "physical_r_and_m_branch_guards_enforced",
        "joined_endpoint_parameter_boxes_built",
        "two_complete_current_microboxes_certified",
    ):
        require(decisions.get(key) is True, f"missing positive decision: {key}")
    for key in (
        "direct_source_length_first_difference_boxes_scale_to_physical_quadrature",
        "small_tau_branch_built",
        "recursive_interval_theta_evaluator_built",
        "physical_quadrature_completed",
        "non_A_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    for stored, case in zip(artifact["complete_current_microboxes"], builder.COMPLETE_BOXES):
        replay = builder.complete_current_parameter_box(case, replay_settings, prior)
        production_complete = parse_acb(stored["complete_current_box"])
        replay_complete = parse_acb(replay["complete_current_box"])
        production_source = parse_acb(stored["center_source_current"])
        require(production_complete.contains(production_source), "production box misses center source")
        require(replay_complete.contains(production_source), "replay box misses center source")
        require(production_complete.overlaps(replay_complete), "complete microbox replay mismatch")
        r = stored["endpoint"]["integer_shift_r"]
        m = stored["endpoint"]["transformed_m"]
        closed = builder.polynomial_masses(r, m)
        # Small m witnesses independently verify the closed polynomial identities.
        for witness_m in (0, 1, 7, 19):
            require(builder.polynomial_masses(r, witness_m) == exact_masses(r, witness_m), "moment formula failed")
        require(tuple(int(stored["main_variation"][f"S{j}"]) for j in range(3)) == closed, "stored moment masses drift")

    audit = artifact["direct_phase_scale_audit"]
    require(len(audit) == 10, "scale-audit row count drift")
    for center in (("1/2", "0"), ("2/5", "1/3")):
        rows = [row for row in audit if (row["center_x"], row["center_s"]) == center]
        require([row["half_width_decimal"] for row in rows] == ["1e-12", "1e-16", "1e-20", "1e-24", "1e-28"], "scale widths drift")
        errors = [arb(row["main_variation_bound"]) for row in rows]
        require(all(bool(errors[i] > errors[i + 1]) for i in range(4)), "scale audit is not strictly descending")
        require(bool(errors[3] > arb(1)) and bool(errors[4] < arb("0.01")), "recorded scale barrier drift")

    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")

    note = NOTE.read_text(encoding="utf-8")
    for token in ("(PB1)", "(PB2)", "(PB4)", "(PB5)", "first-difference box strategy"):
        require(token in note, f"note token missing: {token}")
    normalized = " ".join(note.split()).lower()
    require("scalable recursion and physical quadrature remain open" in normalized, "status boundary missing")
    require("not an impossibility theorem" in normalized and "rh" in normalized, "proof boundary missing")
    print(
        "independently replayed Mordell parameter boxes, recurrence walls, and direct phase-scale barrier",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
