#!/usr/bin/env python3
"""Independently check the offset-20 grouped-752 phase-transport gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_offset20_grouped_752_phase_transport_monotonicity_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    require(path.is_file(), f"missing artifact: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    ctx.prec = 512
    ctx.threads = 1
    artifact = load_json(gate.RESULT)
    require(artifact.get("passed") is True, "phase-transport gate did not pass")
    require(
        artifact.get("status")
        == "offset20_grouped_752_phase_transport_derivative_and_sign_extension_certified",
        "status drift",
    )
    for dependency in artifact["dependencies"].values():
        path = ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")
    for relative, expected in artifact["sources"].items():
        path = ROOT / relative
        require(path.is_file(), f"missing source: {path}")
        require(file_hash(path) == expected, f"source hash drift: {path}")

    panel_cache = gate.transport.load_jsonl(ROOT / artifact["panel_cache_path"])
    component_cache = gate.transport.load_jsonl(ROOT / artifact["component_cache_path"])
    for selected in artifact["selected_component_rows"].values():
        key = selected["cache_key"]
        require(key in component_cache, f"selected component row disappeared: {key}")
        require(
            gate.complete.canonical_hash(component_cache[key])
            == selected["canonical_row_sha256"],
            f"selected component row changed: {key}",
        )

    production_panels = gate.select_panel_rows(
        panel_cache,
        panels=gate.PRODUCTION_PANELS,
        precision_bits=gate.PRODUCTION_PANEL_PRECISION,
        variant="production",
    )
    independent_panels = gate.select_panel_rows(
        panel_cache,
        panels=gate.INDEPENDENT_PANELS,
        precision_bits=gate.INDEPENDENT_PANEL_PRECISION,
        variant="independent",
    )
    require(
        gate.panel_manifest(production_panels)
        == artifact["certificate"]["production"]["panel_manifest"],
        "production panel manifest drift",
    )
    require(
        gate.panel_manifest(independent_panels)
        == artifact["certificate"]["independent"]["panel_manifest"],
        "independent panel manifest drift",
    )

    replay = gate.assemble_certificate(panel_cache, component_cache)
    require(replay == artifact["certificate"], "structural replay drift")
    production = gate.real(replay["production"]["Q_K_minus_T_derivative"])
    independent = gate.real(replay["independent"]["Q_K_minus_T_derivative"])
    require(production.lower() > 0, "production derivative lost positivity")
    require(independent.lower() > 0, "independent derivative lost positivity")
    require(production.overlaps(independent), "production/independent derivatives miss")

    identity = replay["panel_transport_identity"]
    require(identity["base_panel_error"] == "r*epsilon_j*B_j", "base error identity drift")
    require(
        identity["l1_majorant"] == "B_j=752*(y_right-y_left)/sqrt(y_left)",
        "L1 majorant drift",
    )
    require(
        identity["stable_phase"] == "theta(t_c)+h*(theta'(D)-s0)",
        "stable phase identity drift",
    )

    extension = replay["monotonicity_sign_extension"]
    require(extension["closed_interval"][0] == "10000000019.93", "extended left endpoint drift")
    require(extension["closed_interval"][1] == "10000000020.00012", "extended right endpoint drift")
    require(extension["remaining_open_gap"] == "19.9299", "remaining gap drift")
    require(extension["Q_K_minus_T_strictly_negative_on_closed_interval"] is True, "sign extension lost")
    loss = gate.real(replay["radius_diagnostic"]["production_loss_derivative"])
    require(loss.lower() < 0 < loss.upper(), "loss diagnostic no longer straddles zero")

    decision = artifact["decision"]
    for key in (
        "grouped_752_common_phase_transport_proved",
        "offset20_radius_0p07_complete_derivative_positive",
        "production_independent_reproduction_passed",
        "offset20_leftward_sign_extension_to_19p93_proved",
    ):
        require(decision[key] is True, f"proved decision lost at {key}")
    for key in (
        "connected_interval_from_first_subcell_proved",
        "maximal_event_cell_sign_proved",
        "event_wall_handoff_proved",
        "all_height_transport_theorem_proved",
        "rh_implication",
    ):
        require(decision[key] is False, f"overclaim at {key}")

    note = " ".join(gate.NOTE.read_text(encoding="utf-8").split()).lower()
    for token in (
        "common rotation",
        "rigorous representation change, not sampled interpolation",
        "real mean-value theorem",
        "exact remaining gap",
        "enclosure-width diagnostic only",
        "ordinary circumference-to-diameter constant",
        "rh, or a clay-prize conclusion",
    ):
        require(token in note, f"note token missing: {token}")

    print(
        "independently checked offset-20 grouped-752 phase transport; "
        f"production_lower={production.lower().str(20, more=True)}, "
        f"independent_lower={independent.lower().str(20, more=True)}, "
        "extended_left=10000000019.93",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
