#!/usr/bin/env python3
"""Independently replay the selector completion-sector nonlocality bound."""

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
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_selector_completion_sector_nonlocality_gate"
)
packet = gate.packet_gate
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

    dependencies = {
        name: json.loads((REPO_ROOT / record["path"]).read_text(encoding="utf-8"))
        for name, record in artifact["dependencies"].items()
    }
    c_value = gate.arb(packet.C)
    beta = (gate.arb.pi() * c_value * c_value / 8) ** (gate.arb(1) / 3)
    lambda_max = gate.arb.pi() / (16 * beta)
    lam = gate.arb(lambda_max / 2, lambda_max / 2)
    y_max = gate.arb.pi() * c_value / (2 * beta)

    # Rebuild the endpoint formula directly rather than calling profile_derivatives.
    u, v, _, _, _, _ = packet.beta4_weights(lam, gate.arb(0), beta)
    ai, ai_prime, _, _ = packet.acb(-lam).airy()
    g4_origin = y_max * (u * ai - v * ai_prime)
    exact_error = gate.arb(
        dependencies["exact_profile"]["interval_certificate"]["uniform_canonical_exact_minus_beta4_profile_bound"]
    ).upper()
    physical_source_lower = 4 * gate.arb(2).sqrt() * (g4_origin.real.lower() - exact_error)
    require(physical_source_lower.lower() > 230, "independent completed source lost lower barrier 230")

    packet_artifact = dependencies["two_packet_barrier"]
    derivatives = {
        "g": gate.arb(packet_artifact["certificate"]["profile_modulus_bound_ball"]).upper(),
        "g_s": gate.arb(packet_artifact["certificate"]["profile_first_derivative_bound_ball"]).upper(),
        "g_ss": gate.arb(packet_artifact["certificate"]["profile_second_derivative_bound_ball"]).upper(),
    }
    low = packet.low_packet_value(lam, beta, y_max, panels=CHECK_PANELS)
    low_error = packet.low_packet_midpoint_error(derivatives, panels=CHECK_PANELS)
    high = packet.high_packet_by_parts(lam, beta, y_max, derivatives["g_ss"])
    exact_packet_error = gate.arb(packet_artifact["certificate"]["exact_profile_two_packet_error_ball"]).upper()
    alternate_fold_upper = 4 * gate.arb(2).sqrt() * (
        low.real.upper() + low_error + high["center"].real.upper() + high["remainder"] + exact_packet_error
    )
    require(alternate_fold_upper.upper() < 40, "alternate-grid fold upper barrier lost margin 40")
    require((physical_source_lower - alternate_fold_upper).lower() > 190, "alternate outside completion lost margin 190")

    saved = artifact["certificate"]
    require(gate.arb(saved["physical_outside_completion_lower_ball"]).lower() > 200, "saved outside barrier drift")
    require(gate.arb(saved["fold_to_source_fraction_upper_ball"]).upper() < gate.arb("0.14"), "saved fraction drift")
    decision = artifact["decision"]
    require(decision["fold_and_outside_completion_balance_proved"] is True, "completion identity drift")
    require(decision["fold_selector_packet_is_cancellation_closed"] is False, "invalid local closure promoted")
    require(decision["fixed_state_endpoint_residual_bounded"] is False, "fixed-state residual overpromoted")
    require(decision["complete_T_upper_proved"] is False, "T_upper overpromoted")
    print("PASS: independent selector completion-sector barrier; global fixed-state endpoint residual remains open")


if __name__ == "__main__":
    main()
