#!/usr/bin/env python3
"""Independently check the transition--upper-arc K pair on I_1."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_transition_height_derivative_scout as transition_scout
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_transition_upper_arc_K_pair_first_subcell_gate as gate
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_upper_arc_K_first_subcell_scout as arc_scout


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    require(path.is_file(), f"missing artifact: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    artifact = load_json(gate.RESULT)
    require(artifact.get("passed") is True, "production pair artifact did not pass")
    require(
        artifact["status"]
        == "joined_transition_upper_arc_K_pair_first_height_subcell_interval_certified",
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
    require(decision["transition_upper_arc_K_pair_interval_proved"] is True, "pair claim lost")
    require(
        decision["paired_Hardy_derivative_absolute_upper_below_1_3e_minus_6"] is True,
        "pair cap claim lost",
    )
    require(decision["complete_K_T_interval_proved"] is False, "complete K_T overclaim")
    require(decision["wider_Q_K_minus_T_sign_interval_proved"] is False, "sign overclaim")
    require(decision["maximal_event_cell_sign_proved"] is False, "event-cell overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    ctx.prec = 448
    ctx.threads = 1
    transition = transition_scout.build("0.0001", 448, "independent")
    arc = arc_scout.build("0.0001", 448, "independent")
    transition_projection = arb(transition["transition_Q_derivative_contribution"]["ball"])
    arc_projection = arb(arc["upper_arc_Q_derivative_contribution"]["ball"])
    pair = transition_projection + arc_projection
    production = artifact["certificate"]
    production_pair = arb(production["paired_Hardy_derivative_contribution"]["ball"])
    require(pair.overlaps(production_pair), "independent pair misses production")
    require(transition_projection.lower() > 0, "independent transition sign lost")
    require(arc_projection.upper() < 0, "independent arc sign lost")
    require(abs(pair).upper() < arb("1.3e-6"), "independent pair misses cap")
    require(transition["configuration"]["jet_order"] == 16, "transition jet change lost")
    require(transition["configuration"]["derivative_slabs"] == 384, "transition slab change lost")
    require(arc["configuration"]["degree"] == 26, "arc degree change lost")
    require(arc["configuration"]["delta_cut"] == "[0.003200000000000000000000000000000000000000000000000000000 +/- 2.36e-138]" or str(arc["configuration"]["delta_cut"]).startswith("[0.0032"), "arc cutoff change lost")
    require(transition["direct_and_assembled_transition_K_overlap"] is True, "transition overlap lost")
    require(arc["direct_and_assembled_K_overlap"] is True, "arc overlap lost")

    note = " ".join(gate.NOTE.read_text(encoding="utf-8").split()).lower()
    for token in (
        "sup_i1 |hardy_t[k_tr+k_u]/h| < 1.3e-6",
        "theta'-log(y)",
        "does not prove a complete",
        "no circle fit",
    ):
        require(token in note, f"note token missing: {token}")

    print(
        "independently checked transition-upper-arc K pair first subcell; "
        f"production={production_pair.str(16, more=True)}, replay={pair.str(16, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
