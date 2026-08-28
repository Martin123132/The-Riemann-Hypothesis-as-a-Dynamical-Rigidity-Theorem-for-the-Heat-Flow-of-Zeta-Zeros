#!/usr/bin/env python3
"""Independently check the joined lower-plus-ordinary K certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

from flint import acb, arb, ctx
import mpmath as mp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_lower_ordinary_K_first_subcell_gate as gate
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_lower_ordinary_K_first_subcell_scout as scout


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing artifact: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def complex_from(record: dict[str, Any]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def endpoint_record(labels: list[arb], p: arb, y: arb, maximum_power: int) -> dict[str, Any]:
    """Build endpoint denominator sums directly, without a Hurwitz formula."""

    qy = scout.ibp.q_curve(y, p)
    rows: list[dict[str, Any]] = []
    for power in range(1, maximum_power + 1):
        value = sum(
            (
                (acb(0, 2 * arb.pi() * q * y)).exp() / (q - qy) ** power
                for q in labels
            ),
            acb(0),
        )
        rows.append(
            {
                "power": power,
                "phase_weighted_sum_ball": scout.lower_k.complex_record(value),
            }
        )
    return {"endpoint_y": y.str(80, more=True), "rows": rows}


def direct_weighted_integral(labels: list[str], dps: int, cuts: list[str]) -> mp.mpc:
    """Fresh direct quadrature for (C-i log y)y^-1/2 times each phase."""

    with mp.workdps(dps):
        t = mp.mpf("100")
        hardy_shift = mp.mpc("0.3", "0.7")
        q_values = [mp.mpf(label) for label in labels]
        points = [mp.mpf(cut) for cut in cuts]

        def one_label(q: mp.mpf) -> mp.mpc:
            def integrand(y: mp.mpf) -> mp.mpc:
                phase = -t * mp.log(y) - mp.pi * y * y + 2 * mp.pi * q * y
                return (
                    (hardy_shift - 1j * mp.log(y))
                    / mp.sqrt(y)
                    * mp.exp(1j * phase)
                )

            return mp.quad(integrand, points)

        return +sum((one_label(q) for q in q_values), mp.mpc(0))


def altered_orientation_fixture(orientation: str) -> dict[str, Any]:
    """Check one sign of D=q-y-p/y against a changed direct quadrature."""

    require(orientation in ("upper", "lower"), "unknown fixture orientation")
    order = 7
    p = arb(100) / (2 * arb.pi())
    t = arb(100)
    left = arb("1.5")
    right = arb("2.0")
    hardy_shift = acb(arb("0.3"), arb("0.7"))
    label_text = (
        ["20.5", "21.5", "22.5", "23.5"]
        if orientation == "upper"
        else ["0.5", "-0.5", "-1.5", "-2.5"]
    )
    labels = [arb(label) for label in label_text]
    q0 = labels[0]
    base_table = scout.ibp.coefficient_table(order)
    log_table = scout.ordinary_logarithmic_table(order)
    maximum_power = max(
        denominator_power + 1
        for level in base_table[:order]
        for denominator_power in level
    )
    endpoints = {
        "left": endpoint_record(labels, p, left, maximum_power),
        "right": endpoint_record(labels, p, right, maximum_power),
    }
    recurrence = scout.weighted_family(
        name=f"altered_{orientation}_four_label_weighted_fixture",
        base_table=base_table,
        log_table=log_table,
        p=p,
        t=t,
        endpoints=endpoints,
        left_key="left",
        right_key="right",
        left=left,
        right=right,
        q0=q0,
        orientation=orientation,
        finite_count=len(labels),
        slabs=2048,
        cluster_power=4,
        order=order,
        hardy_shift=hardy_shift,
    )
    enclosure = complex_from(recurrence["weighted_tail_ball"])
    direct_90 = direct_weighted_integral(
        label_text, 90, ["1.5", "1.625", "1.75", "1.875", "2.0"]
    )
    direct_120 = direct_weighted_integral(
        label_text, 120, ["1.5", "1.55", "1.7", "1.85", "1.95", "2.0"]
    )
    stability = abs(direct_120 - direct_90)
    direct_ball = acb(
        arb(mp.nstr(direct_120.real, 105)),
        arb(mp.nstr(direct_120.imag, 105)),
    )
    midpoint = mp.mpc(
        mp.mpf(repr(float(enclosure.real.mid()))),
        mp.mpf(repr(float(enclosure.imag.mid()))),
    )
    discrepancy = abs(direct_120 - midpoint)
    require(stability < mp.mpf("1e-35"), f"{orientation} direct quadrature drift")
    require(discrepancy < mp.mpf("1e-8"), f"{orientation} recurrence midpoint drift")
    require(enclosure.overlaps(direct_ball), f"{orientation} recurrence misses direct quadrature")
    return {
        "orientation": orientation,
        "labels": label_text,
        "direct_stability": mp.nstr(stability, 20),
        "direct_to_recurrence_discrepancy": mp.nstr(discrepancy, 20),
        "recurrence_remainder": recurrence["remainder_census"]
        ["scaled_remainder_bound_ball"]["ball"],
    }


def main() -> int:
    artifact = load_json(gate.RESULT)
    require(artifact.get("passed") is True, "production artifact did not pass")
    require(
        artifact["status"]
        == "joined_lower_plus_ordinary_K_first_nonzero_height_subcell_interval_certified",
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
    require(
        decision["joined_lower_plus_ordinary_K_nonzero_radius_interval_proved"] is True,
        "joined K claim lost",
    )
    require(
        decision["complete_lower_plus_ordinary_K_interval_proved"] is True,
        "complete lower-plus-ordinary claim lost",
    )
    require(
        decision["joined_lower_plus_ordinary_Hardy_contribution_strictly_positive"]
        is True,
        "joined contribution sign lost",
    )
    for key in (
        "complete_K_T_interval_proved",
        "wider_Q_K_minus_T_sign_interval_proved",
        "maximal_event_cell_sign_proved",
        "rh_implication",
    ):
        require(decision[key] is False, f"overclaim at {key}")

    ctx.prec = 448
    ctx.threads = 1
    replay = scout.build("0.0001", 448, "independent")
    production = artifact["certificate"]
    production_K = complex_from(production["joined_lower_ordinary_K"])
    replay_K = complex_from(replay["joined_lower_ordinary_K"])
    production_projection = arb(
        production["joined_lower_ordinary_Q_derivative_contribution"]["ball"]
    )
    replay_projection = arb(
        replay["joined_lower_ordinary_Q_derivative_contribution"]["ball"]
    )
    production_error = arb(production["joined_lower_ordinary_total_error"]["ball"])
    replay_error = arb(replay["joined_lower_ordinary_total_error"]["ball"])
    require(production_K.overlaps(replay_K), "independent joined K misses production")
    require(
        production_projection.overlaps(replay_projection),
        "independent joined projection misses production",
    )
    require(
        production_projection.lower() > 0 and replay_projection.lower() > 0,
        "strict positive contribution lost",
    )
    require(production_error.upper() < arb("1e-5"), "production remainder widened")
    require(replay_error.upper() < arb("1e-5"), "independent remainder widened")
    require(replay_K.real.rad() < arb("3e-5"), "independent real radius widened")
    require(replay_K.imag.rad() < arb("3e-5"), "independent imag radius widened")
    require(replay["variant"] == "independent", "independent variant drift")
    require(replay["component_sum_overlap"] is True, "independent component overlap lost")
    config = replay["configuration"]
    expected = {
        "lower_order": 20,
        "lower_split": "280",
        "lower_slabs": 6144,
        "ordinary_lower_order": 4,
        "ordinary_lower_slabs": 20000,
        "upper_order": 9,
        "upper_slabs": 10000,
        "cluster_power": 5,
        "grouped_panels": 256,
        "endpoint_max_power": 20,
    }
    require(config == expected, "independent configuration drift")
    require(
        len(replay["complete_lattice_L_endpoint_currents"]) == 4,
        "independent L-current census drift",
    )

    upper_fixture = altered_orientation_fixture("upper")
    lower_fixture = altered_orientation_fixture("lower")

    note = " ".join(gate.NOTE.read_text(encoding="utf-8").split()).lower()
    for token in (
        "complete signed complex coefficient",
        "four-label direct quadratures",
        "positive-real tail derivative and tiny target correction",
        "no fitted geometric constant",
        "no complete `k_t` bound",
    ):
        require(token in note, f"note token missing: {token}")

    print(
        "independently checked joined lower-plus-ordinary K first height subcell; "
        f"production={production_projection.str(16, more=True)}, "
        f"replay={replay_projection.str(16, more=True)}, "
        f"upper_fixture={upper_fixture['direct_to_recurrence_discrepancy']}, "
        f"lower_fixture={lower_fixture['direct_to_recurrence_discrepancy']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
