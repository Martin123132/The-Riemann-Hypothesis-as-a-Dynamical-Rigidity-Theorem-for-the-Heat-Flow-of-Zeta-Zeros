#!/usr/bin/env python3
"""Independently check the offset-20 local sign-box promotion."""

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

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_offset20_local_box_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    require(path.is_file(), f"missing artifact: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def complex_from(record: dict) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> int:
    ctx.prec = 640
    ctx.threads = 1
    artifact = load_json(gate.RESULT)
    require(artifact.get("passed") is True, "local-box gate did not pass")
    require(
        artifact["status"] == "offset20_closed_local_Q_K_minus_T_negative_interval_certified",
        "status drift",
    )

    for row in artifact["dependencies"].values():
        path = ROOT / row["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == row["sha256"], f"dependency hash drift: {path}")
    for relative, expected in artifact["sources"].items():
        path = ROOT / relative
        require(path.is_file(), f"missing source: {path}")
        require(file_hash(path) == expected, f"source hash drift: {path}")

    rows = gate.load_cache(ROOT / artifact["cache_path"])
    for selected in artifact["selected_append_only_cache_rows"].values():
        key = selected["cache_key"]
        require(key in rows, f"selected cache row disappeared: {key}")
        require(
            gate.row_hash(rows[key]) == selected["canonical_row_sha256"],
            f"selected cache row changed: {key}",
        )

    replay = gate.assemble_certificate(rows)
    require(replay == artifact["certificate"], "higher-precision structural replay drift")
    c = replay
    q_production = gate.real_from(c["production_Q_K_minus_T"])
    q_independent = gate.real_from(c["independent_Q_K_minus_T"])
    require(q_production.upper() < 0 and q_independent.upper() < 0, "replay lost negativity")
    require(q_production.overlaps(q_independent), "replay Q_K-T boxes miss")
    require(
        complex_from(c["production_complete_packet"]).overlaps(
            complex_from(c["independent_complete_packet"])
        ),
        "replay complete packets miss",
    )
    require(
        complex_from(c["production_joined_transition"]).overlaps(
            complex_from(c["independent_joined_transition"])
        ),
        "replay transition packets miss",
    )
    loss = gate.real_from(c["direct_radius_bracket"]["loss_Q_K_minus_T"])
    require(loss.lower() < 0 < loss.upper(), "radius-loss box no longer straddles zero")

    decision = artifact["decision"]
    for key in (
        "offset20_local_closed_box_sign_proved",
        "production_independent_reproduction_passed",
        "direct_radius_width_bracket_measured",
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
        "q_k(t)-t(t)<0 for every t in i_20",
        "interval-dependency/width obstruction only",
        "where pi comes from",
        "ordinary circle constant",
        "uncovered gap of almost 20 height units",
        "no derivative sign on that gap",
        "rh, or clay-prize conclusion",
    ):
        require(token in note, f"note token missing: {token}")

    print(
        "independently checked offset-20 local sign box; "
        f"production_upper={q_production.upper().str(20, more=True)}, "
        f"independent_upper={q_independent.upper().str(20, more=True)}, "
        f"first_loss_upper={loss.upper().str(20, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
