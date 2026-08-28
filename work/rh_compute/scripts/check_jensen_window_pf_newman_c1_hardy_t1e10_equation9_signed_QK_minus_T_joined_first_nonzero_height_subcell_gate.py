#!/usr/bin/env python3
"""Independently check the first nonzero joined Q_K-T height subcell."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_first_nonzero_height_subcell_gate as gate
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_first_subcell_interval_scout as scout


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    require(path.is_file(), f"missing artifact: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def validate_hashes(artifact: dict) -> None:
    for record in artifact["dependencies"].values():
        path = ROOT / record["path"]
        require(path.is_file(), f"missing dependency source: {path}")
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for relative, expected in artifact["sources"].items():
        path = ROOT / relative
        require(path.is_file(), f"missing gate source: {path}")
        require(file_hash(path) == expected, f"source hash drift: {path}")


def main() -> int:
    artifact = load_json(gate.RESULT)
    require(artifact.get("passed") is True, "production artifact did not pass")
    require(
        artifact["status"]
        == "first_nonzero_fixed_roster_joined_QK_minus_T_negative_height_subcell_certified",
        "status drift",
    )
    validate_hashes(artifact)
    decision = artifact["decision"]
    require(decision["first_nonzero_radius_Q_K_minus_T_sign_theorem_proved"] is True, "sign theorem missing")
    require(decision["A_free_join_retained_before_projection"] is True, "A-free join drift")
    require(decision["maximal_event_cell_sign_proved"] is False, "event-cell overclaim")
    require(decision["all_height_theorem"] is False, "all-height overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    ctx.prec = 448
    ctx.threads = 1
    left = scout.build(
        "0.00005",
        448,
        center_text="9999999999.99995",
        variant="independent",
    )
    right = scout.build(
        "0.00005",
        448,
        center_text="10000000000.00005",
        variant="independent",
    )
    production = arb(
        artifact["certificate"]["Q_K_minus_T_height_box"]["ball"]
    )
    left_qkt = arb(left["projection"]["Q_K_minus_T_height_box"]["ball"])
    right_qkt = arb(right["projection"]["Q_K_minus_T_height_box"]["ball"])
    require(left_qkt.upper() < 0, "left independent half-box lost negativity")
    require(right_qkt.upper() < 0, "right independent half-box lost negativity")
    require(left_qkt.overlaps(production), "left half-box misses production enclosure")
    require(right_qkt.overlaps(production), "right half-box misses production enclosure")

    height = arb("10000000000")
    radius = arb("0.0001")
    left_box = arb(left["height_ball"])
    right_box = arb(right["height_ball"])
    require(left_box.lower() <= (height - radius).lower(), "left half-box misses lower endpoint")
    require(left_box.upper() >= height.lower(), "left half-box misses center")
    require(right_box.lower() <= height.upper(), "right half-box misses center")
    require(right_box.upper() >= (height + radius).upper(), "right half-box misses upper endpoint")

    for name, replay in (("left", left), ("right", right)):
        config = replay["configuration"]
        require(config["arc_variant"] == "independent", f"{name} arc variant drift")
        require(config["ordinary_rounds"] == 4, f"{name} ordinary recurrence drift")
        require(config["ordinary_slabs"] == 20000, f"{name} ordinary slabs drift")
        require(config["transition_jet_order"] == 14, f"{name} transition jet drift")
        require(config["upper_transition_panels"] == 256, f"{name} grouped panels drift")
        action = replay["components"]["upper_arc_action_guard"]
        require(action["identity"] == "d A_max(t)/dt=delta_*(t)", f"{name} action identity drift")
        require(action["all_donor_error_bounds_retained"] is True, f"{name} arc error guard drift")
        H = arb(replay["projection"]["H_ball"]["ball"])
        require(H.lower() > 0, f"{name} H box is not positive")
        require(H.rad() < arb("1e-25"), f"{name} stable H box widened")

    note = " ".join(gate.NOTE.read_text(encoding="utf-8").split()).lower()
    for token in (
        "q_k(t)-t(t)=hardy_t[y_t(t)]/h(t)",
        "d a_max(t)/dt=delta_*(t)",
        "upper=",
        "does not extend the sign",
        "no fitted geometric constant",
    ):
        require(token in note, f"note token missing: {token}")

    print(
        "independently checked first nonzero joined Q_K-T height subcell; "
        f"left_upper={left_qkt.upper().str(16, more=True)}, "
        f"right_upper={right_qkt.upper().str(16, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
