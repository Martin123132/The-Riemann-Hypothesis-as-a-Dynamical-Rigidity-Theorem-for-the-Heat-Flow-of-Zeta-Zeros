#!/usr/bin/env python3
"""Independently check the complete K_T first-subcell assembly."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_complete_K_T_first_subcell_gate as gate
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_height_derivative_identity_gate as joined
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_lower_finite_cell_K_first_subcell_scout as lower_scout


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


def altered_linearity_fixture() -> tuple[arb, arb]:
    theta = arb("3.125", "0.0002")
    H = arb("1.125", "0.0001")
    packets = (
        acb(arb("0.002", "1e-7"), arb("-0.003", "2e-7")),
        acb(arb("-0.0008", "3e-8"), arb("0.0011", "4e-8")),
        acb(arb("2e-7", "1e-9"), arb("-3e-7", "1e-9")),
        acb(arb(0, arb("2e-12")), arb(0, arb("2e-12"))),
    )
    direct = joined.hardy_projection(theta, sum(packets, acb(0))) / H
    separate = sum(
        (joined.hardy_projection(theta, packet) / H for packet in packets),
        arb(0),
    )
    require(direct.overlaps(separate), "altered Hardy linearity fixture failed")
    return direct, separate


def main() -> int:
    priority = joined.set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")
    ctx.prec = 512
    ctx.threads = 1

    artifact = load_json(gate.RESULT)
    require(artifact.get("passed") is True, "production complete K_T gate did not pass")
    require(
        artifact["status"]
        == "complete_K_T_and_positive_Q_K_minus_T_derivative_first_height_subcell_certified",
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

    decision = artifact["decision"]
    for key in (
        "complete_K_T_interval_proved",
        "complete_K_T_Hardy_projection_interval_proved",
        "Q_K_minus_T_derivative_strictly_positive_on_I1",
        "Q_K_minus_T_strictly_increasing_on_I1",
        "first_nonzero_radius_Q_K_minus_T_sign_theorem_preserved",
    ):
        require(decision[key] is True, f"proved decision lost at {key}")
    for key in (
        "wider_Q_K_minus_T_sign_interval_proved",
        "maximal_event_cell_sign_proved",
        "event_wall_handoff_proved",
        "all_height_transport_theorem_proved",
        "rh_implication",
    ):
        require(decision[key] is False, f"overclaim at {key}")

    dependencies = {
        name: load_json(ROOT / row["path"])
        for name, row in artifact["dependencies"].items()
    }
    pair = dependencies["transition_upper_arc"]["certificate"]
    lower = dependencies["lower_ordinary"]["certificate"]
    tail = dependencies["positive_real_tail"]["certificate"]
    correction = dependencies["tiny_target_correction"]["certificate"]

    transition_K = complex_from(pair["transition_K_direct"])
    upper_arc_K = complex_from(pair["upper_arc_K_direct"])
    lower_K = complex_from(lower["joined_lower_ordinary_K"])
    correction_K = complex_from(correction["compact_weighted_mode_sum_K_corr"])
    tail_modulus = arb(tail["modulus_bound"]["upper"]).upper()
    tail_K = acb(arb(0, tail_modulus), arb(0, tail_modulus))

    # Deliberately change the production addition order and independently
    # reconstruct theta and H at higher precision.
    complete_K = correction_K + upper_arc_K + tail_K + lower_K + transition_K
    t = arb(arb(10_000_000_000), arb("0.0001"))
    theta, _ = joined.theta_and_derivative(t)
    H, _ = lower_scout.stable_H_box(t)
    require(H.lower() > 0, "independent H transport lost positivity")
    direct_projection = joined.hardy_projection(theta, complete_K) / H

    lower_projection = arb(
        lower["joined_lower_ordinary_Q_derivative_contribution"]["ball"]
    )
    pair_projection = arb(pair["paired_Hardy_derivative_contribution"]["ball"])
    correction_projection = arb(correction["Hardy_projection_over_H"]["ball"])
    tail_projection_bound = arb(
        tail["Hardy_projection_over_H_absolute_bound"]["upper"]
    ).upper()
    replay_projection = (
        correction_projection
        + arb(0, tail_projection_bound)
        + pair_projection
        + lower_projection
    )
    envelope_projection = (
        correction_projection
        + lower_projection
        + arb(0, arb("1.3e-6") + tail_projection_bound)
    )

    production = artifact["certificate"]
    production_K = complex_from(production["complete_K_T"])
    production_direct = arb(
        production["complete_Hardy_projection_over_H_direct"]["ball"]
    )
    production_sum = arb(
        production["complete_Hardy_projection_over_H_component_sum"]["ball"]
    )
    production_envelope = arb(
        production["complete_Hardy_projection_over_H_pair_envelope"]["ball"]
    )
    require(complete_K.overlaps(production_K), "changed-order K_T misses production")
    require(direct_projection.overlaps(production_direct), "direct projection misses production")
    require(replay_projection.overlaps(production_sum), "projection sum misses production")
    require(envelope_projection.overlaps(production_envelope), "pair envelope misses production")
    require(direct_projection.overlaps(replay_projection), "direct and summed projections miss")
    require(direct_projection.lower() > 0, "independent direct projection lost positivity")
    require(replay_projection.lower() > 0, "independent projection sum lost positivity")
    require(envelope_projection.lower() > 0, "independent pair envelope lost positivity")

    fixture_direct, fixture_sum = altered_linearity_fixture()
    note = " ".join(gate.NOTE.read_text(encoding="utf-8").split()).lower()
    for token in (
        "k_t=k_(l+o)+(k_tr+k_u)+k_tail+k_corr",
        "(q_k-t)'=hardy_t[k_t]/h",
        "only its certified absolute envelope `1.3e-6`",
        "no circle, polygon, visual symmetry",
        "no full event-cell theorem",
        "rh, or prize-level conclusion",
    ):
        require(token in note, f"note token missing: {token}")

    print(
        "independently checked complete K_T and positive joined derivative first subcell; "
        f"production={production_sum.str(18, more=True)}, "
        f"replay={replay_projection.str(18, more=True)}, "
        f"envelope={envelope_projection.str(18, more=True)}, "
        f"fixture_overlap={fixture_direct.overlaps(fixture_sum)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
