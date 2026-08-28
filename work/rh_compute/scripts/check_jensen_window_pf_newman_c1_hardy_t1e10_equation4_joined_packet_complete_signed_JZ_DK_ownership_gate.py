#!/usr/bin/env python3
"""Independently check the complete signed J_Z and D_K ownership gate."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx

import check_jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate as interval_check


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_complete_signed_JZ_DK_ownership_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()
HEIGHT = arb("10000000000")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing artifact: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def saved_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def add_error(value: acb, radius: arb) -> acb:
    return acb(arb(value.real, radius), arb(value.imag, radius))


def hardy(theta: arb, value: acb) -> arb:
    return 2 * (acb(0, theta).exp() * value).real


def radius_from_log10(log_record: str) -> arb:
    log_ball = arb(log_record)
    radius = (log_ball.upper() * arb(10).log()).exp().upper()
    require(radius > 0, "logarithmic radius vanished")
    return radius


def central_theta(calibration: dict[str, Any]) -> arb:
    rows = [row for row in calibration["output_rows"] if row["output_label"] == "t+0.00"]
    require(len(rows) == 1, "central calibration row drift")
    return arb(rows[0]["high_precision"]["theta_ball"])


def disjoint(first: arb, second: arb) -> bool:
    return first.upper() < second.lower() or second.upper() < first.lower()


def main() -> int:
    resource_mode = interval_check.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = 512
    ctx.threads = 1
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing artifact: {path}")

    artifact = load_json(RESULT)
    require(artifact.get("passed") is True, "production gate not passed")
    require(
        artifact["decision"]["complete_production_height_J_Z_enclosed"] is True
        and artifact["decision"]["complete_production_height_D_K_enclosed"] is True,
        "production scalar enclosure decision drift",
    )
    require(
        artifact["decision"]["D_K_sufficient_corridor_proved"] is False,
        "failed corridor was overpromoted",
    )
    require(
        artifact["decision"]["natural_A_lift_counted_exactly_once_inside_transition"] is True,
        "natural A ownership drift",
    )

    for source in artifact["sources"].values():
        path = ROOT / source["path"]
        require(path.is_file() and file_hash(path) == source["sha256"], f"source hash drift: {path}")
    deps: dict[str, dict[str, Any]] = {}
    for name, source in artifact["dependencies"].items():
        path = ROOT / source["path"]
        require(path.is_file() and file_hash(path) == source["sha256"], f"dependency hash drift: {path}")
        deps[name] = load_json(path)
    require(all(dep.get("passed") is True for dep in deps.values()), "dependency pass drift")

    lower = saved_complex(deps["lower_cell"]["certificate"]["lower_cell_integral_ball"])
    quarter = deps["quarter_arc"]
    arc = saved_complex(quarter["arb_certificate"]["actual_arc"])
    tail_radius = radius_from_log10(
        quarter["arb_certificate"]["error_budgets"]["positive_real_tail_log10_upper"]["ball"]
    )
    tail = acb(arb(0, tail_radius), arb(0, tail_radius))

    ordinary_artifact = deps["ordinary_packet"]
    lower_tail = saved_complex(ordinary_artifact["certificate"]["lower_complementary_tail_ball"])
    upper_tail = saved_complex(
        deps["upper_complement"]["certificate"]["complete_upper_complement_ball"]
    )
    gamma_radius = radius_from_log10(
        deps["ordinary_half_lattice"]["certificate"][
            "ordinary_Gamma_defect_packet_log10_upper_ball"
        ]["ball"]
    )
    ordinary = add_error(-lower_tail - upper_tail, gamma_radius)
    saved_ordinary = saved_complex(ordinary_artifact["certificate"]["complete_ordinary_packet_ball"])
    require(ordinary.overlaps(saved_ordinary), "independent ordinary reconstruction misses saved ball")

    joined = deps["transition_packet"]["joined_packet"]
    transition_source = saved_complex(
        deps["transition_packet"]["source_certificate"][
            "complete_transition_source_integral_ball"
        ]
    )
    transition_gamma = saved_complex(joined["transition_Gamma_lift_ball"])
    natural_A = saved_complex(joined["grouped_natural_A_lift_ball"])
    transition = transition_source - transition_gamma - natural_A
    saved_transition = saved_complex(joined["joined_transition_packet_ball"])
    require(transition.overlaps(saved_transition), "independent transition reconstruction misses saved ball")

    theta = acb(arb("0.25"), HEIGHT / 2).lgamma().imag - HEIGHT * arb.pi().log() / 2
    saved_theta = central_theta(deps["theta_calibration"])
    production_theta = arb(artifact["certificate"]["theta_ball"]["ball"])
    require(theta.overlaps(saved_theta) and theta.overlaps(production_theta), "changed-precision theta drift")
    require(
        hardy(theta, transition).overlaps(arb(joined["joined_transition_Hardy_projection_ball"])),
        "independent transition phase/sign check failed",
    )

    # Deliberately use a different addition order from the production builder.
    argument = transition + arc + ordinary + lower + tail
    H = arb(deps["exact_H"]["exact_H_certificate"]["high_precision_H_ball"])
    J_Z = hardy(theta, argument) / H
    saved_J_Z = arb(artifact["certificate"]["J_Z_ball"]["ball"])
    require(J_Z.overlaps(saved_J_Z), "changed-precision J_Z misses production")

    projection_sum = (
        hardy(theta, transition)
        + hardy(theta, arc)
        + hardy(theta, ordinary)
        + hardy(theta, lower)
        + hardy(theta, tail)
    )
    require(projection_sum.overlaps(hardy(theta, argument)), "Hardy linearity reconstruction drift")
    saved_theta_J_Z = hardy(saved_theta, argument) / H
    require(saved_theta_J_Z.overlaps(saved_J_Z), "saved-theta independent J_Z misses production")

    A_transition = arb(
        deps["A_transition"]["certificate"]["components"][
            "projector_completed_exact_A_transition_ball"
        ]
    )
    R_KGamma = J_Z + A_transition
    gamma_bound = abs(
        arb(
            deps["Gamma_correction"]["correction_certificate"][
                "Gamma_to_classical_physical_absolute_bound_ball"
            ]
        )
    ).upper()
    G_minus_T = arb(0, gamma_bound)
    bridge = deps["QK_bridge"]["bridge_target_certificate"]
    classical_main = arb(bridge["classical_upper_main_ball"])
    exact_upper = arb(bridge["exact_complementary_upper_ball"])
    rho = exact_upper - classical_main
    QK_minus_T = R_KGamma + G_minus_T
    Delta_KU = QK_minus_T - rho
    D_K = (1 - H) * exact_upper - H * Delta_KU
    D_K_direct = exact_upper - H * (classical_main + QK_minus_T)
    require(D_K.overlaps(D_K_direct), "independent D_K identities do not overlap")

    for name, value in (
        ("R_KGamma_ball", R_KGamma),
        ("Q_K_minus_T_ball", QK_minus_T),
        ("Delta_KU_ball", Delta_KU),
        ("D_K_ball", D_K),
    ):
        require(
            value.overlaps(arb(artifact["certificate"][name]["ball"])),
            f"changed-precision {name} misses production",
        )

    require(
        J_Z.lower() > arb(2_939_975) / 100_000_000,
        "absolute J_Z miss is no longer rigorous",
    )
    require(
        R_KGamma.lower() > -arb(727_437) / 100_000_000,
        "R_KGamma miss is no longer rigorous",
    )
    require(
        D_K.upper() < arb(72_743_687) / 10_000_000_000,
        "D_K miss is no longer rigorous",
    )

    for name, value in (
        ("V_L", lower),
        ("C_U", arc),
        ("O_join", ordinary),
        ("T_A_join", transition),
    ):
        flipped = hardy(theta, argument - 2 * value) / H
        require(disjoint(J_Z, flipped), f"changed sign guard failed for {name}")
    require(
        disjoint(J_Z, hardy(theta, argument + natural_A) / H)
        and disjoint(J_Z, hardy(theta, argument - natural_A) / H),
        "changed natural A duplication guard failed",
    )

    ledger = artifact["certificate"]["ownership_ledger"]
    require(
        [row["packet"] for row in ledger] == ["V_L", "C_U", "E_s_T_U", "O_join", "T_A_join"],
        "top-level ledger order or membership drift",
    )
    require(
        artifact["certificate"]["nested_ownership_guards"][
            "grouped_natural_A_lift_occurrences"
        ]
        == 1,
        "natural A occurrence count drift",
    )

    note = " ".join(NOTE.read_text(encoding="utf-8").split())
    for fragment in (
        "H(t)J_Z=Hardy_t[V_L+C_U+E_s*T_U+O_join+T_A^join]",
        "subtracted exactly once inside",
        "The stronger sufficient route does not close",
        "failure due to enclosure width",
        "does not disprove the original one-sided route",
        "No circle, polygon, fitted constant, or visual pattern supplies",
    ):
        require(fragment.lower() in note.lower(), f"missing note fragment: {fragment}")

    print(
        "independently checked complete signed J_Z and D_K ownership; "
        f"J_Z={J_Z.str(12, more=True)}, D_K={D_K.str(12, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
