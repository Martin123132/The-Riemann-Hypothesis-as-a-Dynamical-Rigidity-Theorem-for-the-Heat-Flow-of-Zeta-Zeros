#!/usr/bin/env python3
"""Independently check the offset-20 complete K_T monotonicity bridge."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_offset20_complete_K_T_monotonicity_gate as gate


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
    require(artifact.get("passed") is True, "monotonicity gate did not pass")
    require(
        artifact.get("status")
        == "offset20_complete_K_T_positive_derivative_and_left_sign_extension_certified",
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

    rows = gate.load_cache(ROOT / artifact["cache_path"])
    for selected in artifact["selected_append_only_component_rows"].values():
        key = selected["cache_key"]
        require(key in rows, f"selected component row disappeared: {key}")
        require(
            gate.scout.canonical_hash(rows[key]) == selected["canonical_row_sha256"],
            f"selected component row changed: {key}",
        )

    replay = gate.assemble_certificate(rows)
    require(replay == artifact["certificate"], "structural replay drift")
    c = replay
    for label in (
        "production_derivative_direct",
        "production_derivative_component_sum",
        "independent_derivative_direct",
        "independent_derivative_component_sum",
    ):
        require(gate.real(c[label]).lower() > 0, f"positive derivative lost at {label}")

    anchor = c["value_anchor"]
    require(gate.real(anchor["production_Q_K_minus_T"]).upper() < 0, "production anchor lost sign")
    require(gate.real(anchor["independent_Q_K_minus_T"]).upper() < 0, "independent anchor lost sign")
    extension = c["monotonicity_sign_extension"]
    left = arb(extension["closed_interval"][0])
    right = arb(extension["closed_interval"][1])
    require(extension["closed_interval"][0] == "10000000019.994", "extended left endpoint drift")
    require(extension["closed_interval"][1] == "10000000020.00012", "extended right endpoint drift")
    require(extension["Q_K_minus_T_strictly_negative_on_closed_interval"] is True, "sign extension lost")
    require(extension["remaining_open_gap"] == "19.9939", "remaining gap drift")

    loss = c["radius_diagnostic"]
    require(
        gate.real(loss["production_loss_direct"]).lower() < 0
        < gate.real(loss["production_loss_direct"]).upper(),
        "direct loss diagnostic no longer straddles zero",
    )
    require(
        gate.real(loss["production_loss_component_sum"]).lower() < 0
        < gate.real(loss["production_loss_component_sum"]).upper(),
        "component loss diagnostic no longer straddles zero",
    )

    decision = artifact["decision"]
    for key in (
        "offset20_complete_K_T_derivative_positive",
        "production_independent_reproduction_passed",
        "offset20_leftward_sign_extension_proved",
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
        "fundamental theorem of calculus",
        "strictly increasing throughout",
        "extends the offset-20 sign island leftward",
        "interval-enclosure width obstruction only",
        "where pi comes from",
        "ordinary circle constant",
        "remaining open gap",
        "rh, or a clay-prize conclusion",
    ):
        require(token in note, f"note token missing: {token}")

    print(
        "independently checked offset-20 complete K_T monotonicity bridge; "
        f"production_lower={gate.real(c['production_derivative_component_sum']).lower().str(20, more=True)}, "
        f"independent_lower={gate.real(c['independent_derivative_component_sum']).lower().str(20, more=True)}, "
        f"extended_left={left.str(20, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
