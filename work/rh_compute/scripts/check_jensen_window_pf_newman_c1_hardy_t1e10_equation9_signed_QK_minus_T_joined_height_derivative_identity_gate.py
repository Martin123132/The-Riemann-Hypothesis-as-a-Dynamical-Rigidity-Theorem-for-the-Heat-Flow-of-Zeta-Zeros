#!/usr/bin/env python3
"""Independently check the A-free joined Q_K-T height derivative gate."""

from __future__ import annotations

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
import mpmath as mp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_height_derivative_identity_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def changed_derivative_fixture() -> None:
    mp.mp.dps = 92
    s = mp.mpc("0.5", "2.875")
    first_label = 7
    count = 4
    low_cuts = (
        mp.mpf(0),
        mp.mpf("0.75"),
        mp.mpf("2.5"),
        mp.mpf("6.5"),
        mp.mpf("14.5"),
        mp.inf,
    )
    high_cuts = tuple(mp.mpf(item) for item in (1, 3, 6, 12, 24, 48, 72))
    direct = gate.direct_source_derivative(s, first_label, count, low_cuts, high_cuts)
    mellin = gate.mellin_source_derivative(s, first_label, count, low_cuts, high_cuts)
    require(abs(direct - mellin) < mp.mpf("1e-55"), "changed derivative fixture failed")


def main() -> int:
    priority = gate.set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")
    require(gate.RESULT.is_file() and gate.NOTE.is_file(), "missing result or note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "artifact is not passed")
    require(
        artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER),
        "builder hash mismatch",
    )
    require(
        artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER),
        "checker hash mismatch",
    )
    for name, path in gate.DEPENDENCIES.items():
        require(
            artifact["dependencies"][name]["sha256"] == gate.file_hash(path),
            f"dependency hash mismatch: {name}",
        )

    deps = {name: gate.load_json(path) for name, path in gate.DEPENDENCIES.items()}
    ctx.prec = 512
    ctx.threads = 1
    t = arb(gate.HEIGHT)
    theta, theta_prime = gate.theta_and_derivative(t)
    log_h_direct = gate.direct_log_H(t)
    log_h_duplicate = gate.duplication_log_H(t)
    h_log_prime_direct = gate.direct_log_H_derivative(t)
    h_log_prime_duplicate = gate.duplication_log_H_derivative(t)
    require(log_h_direct.overlaps(log_h_duplicate), "512-bit H formulas miss")
    require(
        h_log_prime_direct.overlaps(h_log_prime_duplicate),
        "512-bit H derivative formulas miss",
    )
    H = log_h_duplicate.exp()
    saved_phase = artifact["certificate"]["phase_and_prefactor"]
    require(theta.overlaps(arb(saved_phase["theta_ball"]["ball"])), "theta drift")
    require(
        theta_prime.overlaps(arb(saved_phase["theta_prime_ball"]["ball"])),
        "theta-prime drift",
    )
    require(H.overlaps(arb(saved_phase["H_ball"]["ball"])), "H drift")
    require(
        h_log_prime_duplicate.overlaps(
            arb(saved_phase["H_logarithmic_derivative_direct_ball"]["ball"])
        ),
        "H logarithmic derivative drift",
    )

    D_target, D_target_prime = gate.finite_dirichlet(
        t, gate.TARGET_START, gate.TARGET_END, reverse=True
    )
    saved_target = artifact["certificate"]["classical_target"]
    require(D_target.overlaps(parse_complex(saved_target["D_T_ball"])), "D_T reverse sum drift")
    require(
        D_target_prime.overlaps(parse_complex(saved_target["D_T_prime_ball"])),
        "D_T-prime reverse sum drift",
    )
    T = gate.hardy_projection(theta, D_target)
    T_prime = gate.hardy_projection(
        theta, D_target_prime + acb(0, theta_prime) * D_target
    )
    require(T.overlaps(arb(saved_target["T_ball"]["ball"])), "T reverse projection drift")
    require(
        T_prime.overlaps(arb(saved_target["T_prime_ball"]["ball"])),
        "T-prime reverse projection drift",
    )

    transition = deps["transition_packet"]
    C_G = (1 + (-2 * arb.pi() * t).exp()) ** (-arb(1) / 2)
    _, _ = gate.finite_dirichlet(t, gate.LOWER_TRANSITION_START, gate.TARGET_END, reverse=True)
    D_extra, _ = gate.finite_dirichlet(t, gate.EXTRA_START, gate.EXTENDED_END, reverse=True)
    signed = deps["signed_JZ"]["certificate"]
    J_argument = gate.interval_base.saved_complex(signed["complete_preprojection_argument_ball"])
    natural_A = gate.interval_base.saved_complex(
        transition["joined_packet"]["grouped_natural_A_lift_ball"]
    )
    alternate_target = J_argument + natural_A + C_G * D_extra + (C_G - H) * D_target
    alternate_QKT = gate.hardy_projection(theta, alternate_target) / H
    saved_join = artifact["certificate"]["A_free_saved_height_assembly"]
    require(
        alternate_target.overlaps(
            parse_complex(saved_join["joined_classical_target_preprojection_ball"])
        ),
        "changed-order A-free preprojection misses",
    )
    require(
        alternate_QKT.overlaps(arb(saved_join["direct_Q_K_minus_T_ball"]["ball"])),
        "changed-order A-free Q_K-T misses",
    )
    require(alternate_QKT.upper() < 0, "changed-order Q_K-T lost negativity")

    changed_derivative_fixture()
    decision = artifact["decision"]
    require(decision["A_endpoint_cancels_exactly_from_QK_minus_T"] is True, "A cancellation missing")
    require(decision["exact_joined_first_height_derivative_proved"] is True, "derivative identity missing")
    require(
        decision["nonzero_radius_joined_derivative_bound_proved"] is False,
        "nonzero-radius bound overclaimed",
    )
    require(decision["off_height_QK_minus_T_sign_proved"] is False, "off-height sign overclaimed")
    require(decision["rh_implication"] is False, "RH boundary corrupted")
    budget = artifact["certificate"]["first_subcell_budget"]
    require(budget["radius_selected"] is False, "radius selected without a joined bound")
    require(len(budget["candidate_radius_rows"]) == 4, "transport ladder drift")
    print(
        "independently checked A-free joined Q_K-T derivative identity, "
        "reverse target sum, and changed derivative fixture",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
