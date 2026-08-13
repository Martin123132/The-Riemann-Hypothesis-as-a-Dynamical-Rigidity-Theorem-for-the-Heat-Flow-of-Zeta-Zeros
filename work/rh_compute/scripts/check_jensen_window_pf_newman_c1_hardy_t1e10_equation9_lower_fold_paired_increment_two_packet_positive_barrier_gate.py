#!/usr/bin/env python3
"""Check the rigorous positive two-packet selector-increment barrier."""

from __future__ import annotations

import hashlib
import importlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = REPO_ROOT / "work/rh_compute/scripts"
sys.path.insert(0, str(SCRIPT_ROOT))
gate = importlib.import_module(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_increment_two_packet_positive_barrier_gate"
)
RESULT = REPO_ROOT / f"work/rh_compute/results/{gate.STEM}.json"
CHECK_PANELS = 32_768


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    gate.ctx.dps = gate.PRECISION
    require(RESULT.is_file(), "missing result artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "saved gate did not pass")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")

    c = gate.arb(gate.C)
    tstar = gate.arb.pi() * c * c / 8
    beta = tstar ** (gate.arb(1) / 3)
    lambda_max = gate.arb.pi() / (16 * beta)
    lam = gate.arb(lambda_max / 2, lambda_max / 2)
    y_max = gate.arb.pi() * c / (2 * beta)
    derivatives = {
        "g": gate.arb(artifact["certificate"]["profile_modulus_bound_ball"]).upper(),
        "g_s": gate.arb(artifact["certificate"]["profile_first_derivative_bound_ball"]).upper(),
        "g_ss": gate.arb(artifact["certificate"]["profile_second_derivative_bound_ball"]).upper(),
    }
    low = gate.low_packet_value(lam, beta, y_max, panels=CHECK_PANELS)
    low_error = gate.low_packet_midpoint_error(derivatives, panels=CHECK_PANELS)
    low_lower = low.real.lower() - low_error
    low_upper = low.real.upper() + low_error
    require(low_lower > 4, "alternate-grid low packet lost positive margin 4")

    high = gate.high_packet_by_parts(lam, beta, y_max, derivatives["g_ss"])
    saved_high_lower = gate.arb(artifact["certificate"]["high_beta4_real_lower_ball"]).lower()
    saved_high_upper = gate.arb(artifact["certificate"]["high_beta4_real_upper_ball"]).upper()
    require(high["center"].real.lower() - high["remainder"] >= saved_high_lower, "remote packet replay drift")
    require(high["center"].real.upper() + high["remainder"] <= saved_high_upper, "remote packet upper replay drift")

    exact_error = gate.arb(artifact["certificate"]["exact_profile_two_packet_error_ball"]).upper()
    physical_lower = 4 * gate.arb(2).sqrt() * (low_lower + saved_high_lower - exact_error)
    physical_upper = 4 * gate.arb(2).sqrt() * (low_upper + saved_high_upper + exact_error)
    require(physical_lower.lower() > 20, "alternate-grid physical barrier lost margin 20")
    require(physical_upper.upper() < 40, "alternate-grid physical upper barrier lost margin 40")

    decision = artifact["decision"]
    require(decision["exact_two_packet_fold_increment_identity_proved"] is True, "packet identity drift")
    require(decision["selector_fold_increment_uniformly_positive_above_20"] is True, "positive barrier drift")
    require(decision["selector_fold_increment_uniformly_below_35"] is True, "upper barrier drift")
    require(decision["selector_increment_is_small_residual_route"] is False, "invalid small-increment route promoted")
    require(decision["fixed_state_fold_residual_bounded"] is False, "fixed-state residual overpromoted")
    require(decision["complete_T_upper_proved"] is False, "T_upper overpromoted")
    print("PASS: alternate-grid two-packet positive selector-increment barrier; fixed-state residual remains open")


if __name__ == "__main__":
    main()
