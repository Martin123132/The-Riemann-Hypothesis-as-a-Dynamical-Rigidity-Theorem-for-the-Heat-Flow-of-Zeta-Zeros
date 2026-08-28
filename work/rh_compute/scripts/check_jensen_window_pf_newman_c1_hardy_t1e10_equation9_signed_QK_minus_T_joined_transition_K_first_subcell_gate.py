#!/usr/bin/env python3
"""Independently check the joined-transition K first-subcell gate."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_transition_K_first_subcell_gate as gate
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_transition_height_derivative_scout as scout


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
    artifact = load_json(gate.RESULT)
    require(artifact.get("passed") is True, "production artifact did not pass")
    require(
        artifact["status"] == "joined_transition_K_first_nonzero_height_subcell_interval_certified",
        "status drift",
    )
    for record in artifact["dependencies"].values():
        path = ROOT / record["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for relative, expected in artifact["sources"].items():
        path = ROOT / relative
        require(path.is_file(), f"missing source: {path}")
        require(file_hash(path) == expected, f"source hash drift: {path}")

    decision = artifact["decision"]
    require(decision["transition_K_nonzero_radius_interval_proved"] is True, "transition claim lost")
    require(decision["transition_Hardy_derivative_contribution_strictly_positive"] is True, "sign lost")
    require(decision["complete_K_T_interval_proved"] is False, "complete K_T overclaim")
    require(decision["wider_Q_K_minus_T_sign_interval_proved"] is False, "sign extension overclaim")
    require(decision["maximal_event_cell_sign_proved"] is False, "event-cell overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    ctx.prec = 448
    ctx.threads = 1
    replay = scout.build("0.0001", 448, "independent")
    production = artifact["certificate"]
    production_K = complex_from(production["transition_K_direct"])
    replay_K = complex_from(replay["transition_K_direct"])
    production_projection = arb(production["transition_Q_derivative_contribution"]["ball"])
    replay_projection = arb(replay["transition_Q_derivative_contribution"]["ball"])
    require(production_K.overlaps(replay_K), "independent transition K misses production")
    require(production_projection.overlaps(replay_projection), "independent projection misses production")
    require(production_projection.lower() > 0 and replay_projection.lower() > 0, "positive sign lost")
    require(replay["variant"] == "independent", "independent variant drift")
    config = replay["configuration"]
    require(config["jet_order"] == 16, "independent jet order drift")
    require(config["derivative_slabs"] == 384, "independent slab count drift")
    require(config["integral_panels"] == 64, "independent panel count drift")
    require(config["series_terms"] == 40, "independent series order drift")
    require(replay["altered_five_label_fixture"]["overlap"] is True, "altered replay failed")
    require(replay["direct_and_assembled_Gamma_K_overlap"] is True, "Gamma overlap lost")
    require(replay["direct_and_assembled_transition_K_overlap"] is True, "transition overlap lost")

    note = " ".join(gate.NOTE.read_text(encoding="utf-8").split()).lower()
    for token in (
        "i(theta'-log(y))-h'/h".replace("h'/h", "h'/h"),
        "no separately widened",
        "complete `k_t` is not yet enclosed",
        "do not prove a complete",
        "no fitted geometric constant",
    ):
        require(token in note, f"note token missing: {token}")

    print(
        "independently checked joined-transition K first height subcell; "
        f"production={production_projection.str(16, more=True)}, "
        f"replay={replay_projection.str(16, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
