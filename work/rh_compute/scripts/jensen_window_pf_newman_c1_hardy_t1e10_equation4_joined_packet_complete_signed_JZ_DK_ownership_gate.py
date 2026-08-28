#!/usr/bin/env python3
"""Assemble the complete signed production-height J_Z and D_K packets."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate as interval_base


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_complete_signed_JZ_DK_ownership_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)


def result_path(stem: str) -> Path:
    return ROOT / "work" / "rh_compute" / "results" / f"{stem}.json"


DEPENDENCIES = {
    "normalization_repair": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_domain_Kummer_normalization_repair_gate"
    ),
    "carrier_partition": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_cell_common_carrier_subtraction_reduction_gate"
    ),
    "pole_safe_join": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_PW_unowned_cell_pole_safe_two_boundary_contour_join_gate"
    ),
    "lower_cell": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_lower_cell_reversed_roster_endpoint_gate"
    ),
    "quarter_arc": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_quarter_disk_vertical_arc_tail_gate"
    ),
    "ordinary_half_lattice": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_ordinary_complementary_half_lattice_gate"
    ),
    "upper_complement": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_far_remainder_eight_round_gate"
    ),
    "ordinary_packet": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_lower_complement_signed_coefficient_interval_gate"
    ),
    "transition_packet": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_transition_alternating_boundary_jet_packet_gate"
    ),
    "A_transition": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_complete_endpoint_projector_assembly_gate"
    ),
    "exact_H": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_exact_H_continuation_defect_target_gate"
    ),
    "QK_bridge": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_QK_exact_hardy_truncation_bridge_target_gate"
    ),
    "Gamma_correction": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_uniform_gamma_target_insertion_gate"
    ),
    "theta_calibration": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_exact_hardy_arb_calibration_gate"
    ),
}

PRECISION_BITS = 384
HEIGHT = arb("10000000000")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def central_theta(calibration: dict[str, Any]) -> arb:
    rows = [row for row in calibration["output_rows"] if row["output_label"] == "t+0.00"]
    require(len(rows) == 1, "central theta row drift")
    require(rows[0]["target_t"] == "10000000000.00000000000000000", "theta height drift")
    return arb(rows[0]["high_precision"]["theta_ball"])


def hardy_projection(theta: arb, value: acb) -> arb:
    return 2 * (acb(0, theta).exp() * value).real


def gamma_defect_radius(complement: dict[str, Any]) -> tuple[arb, arb]:
    log_ball = arb(
        complement["certificate"]["ordinary_Gamma_defect_packet_log10_upper_ball"]["ball"]
    )
    radius = (log_ball.upper() * arb(10).log()).exp().upper()
    require(radius > 0, "ordinary Gamma defect radius vanished")
    return log_ball, radius


def upper_tail_radius(quarter: dict[str, Any]) -> tuple[arb, arb]:
    log_ball = arb(
        quarter["arb_certificate"]["error_budgets"]["positive_real_tail_log10_upper"]["ball"]
    )
    radius = (log_ball.upper() * arb(10).log()).exp().upper()
    require(radius > 0, "upper positive-real tail radius vanished")
    return log_ball, radius


def disjoint_separation(first: arb, second: arb) -> arb:
    if first.upper() < second.lower():
        return second.lower() - first.upper()
    if second.upper() < first.lower():
        return first.lower() - second.upper()
    return arb(0)


def build_certificate(deps: dict[str, dict[str, Any]]) -> dict[str, Any]:
    jz_absolute_target = arb(2_939_975) / 100_000_000
    r_kgamma_upper_target = -arb(727_437) / 100_000_000
    dk_lower_target = arb(72_743_687) / 10_000_000_000
    dk_upper_target = arb(663_407_986) / 10_000_000_000
    normalization = deps["normalization_repair"]
    carrier = deps["carrier_partition"]
    join = deps["pole_safe_join"]
    require(
        normalization["decision"]["joined_J_Z_identity_retained_in_corrected_normalization"] is True,
        "corrected J_Z normalization drift",
    )
    require(
        carrier["decision"]["finite_Fresnel_cell_partition_exact"] is True,
        "finite cell partition drift",
    )
    require(
        join["decision"]["complete_three_packet_fresnel_cancellation_proved"] is True,
        "three-packet contour join drift",
    )

    lower_cell = interval_base.saved_complex(
        deps["lower_cell"]["certificate"]["lower_cell_integral_ball"]
    )
    quarter = deps["quarter_arc"]
    upper_arc = interval_base.saved_complex(quarter["arb_certificate"]["actual_arc"])
    tail_log10, tail_radius = upper_tail_radius(quarter)
    upper_tail = acb(arb(0, tail_radius), arb(0, tail_radius))
    unowned = lower_cell + upper_arc + upper_tail

    ordinary_artifact = deps["ordinary_packet"]
    ordinary = interval_base.saved_complex(
        ordinary_artifact["certificate"]["complete_ordinary_packet_ball"]
    )
    lower_tail = interval_base.saved_complex(
        ordinary_artifact["certificate"]["lower_complementary_tail_ball"]
    )
    upper_complement = interval_base.saved_complex(
        deps["upper_complement"]["certificate"]["complete_upper_complement_ball"]
    )
    ordinary_gamma_log10, ordinary_gamma_radius = gamma_defect_radius(
        deps["ordinary_half_lattice"]
    )
    reconstructed_ordinary = interval_base.add_complex_error(
        -lower_tail - upper_complement, ordinary_gamma_radius
    )
    require(reconstructed_ordinary.overlaps(ordinary), "ordinary packet sign reconstruction drift")

    transition_artifact = deps["transition_packet"]
    joined = transition_artifact["joined_packet"]
    transition = interval_base.saved_complex(joined["joined_transition_packet_ball"])
    transition_source = interval_base.saved_complex(
        transition_artifact["source_certificate"]["complete_transition_source_integral_ball"]
    )
    transition_gamma = interval_base.saved_complex(joined["transition_Gamma_lift_ball"])
    natural_A = interval_base.saved_complex(joined["grouped_natural_A_lift_ball"])
    reconstructed_transition = transition_source - transition_gamma - natural_A
    require(reconstructed_transition.overlaps(transition), "transition packet sign reconstruction drift")

    theta = acb(arb("0.25"), HEIGHT / 2).lgamma().imag - HEIGHT * arb.pi().log() / 2
    saved_theta = central_theta(deps["theta_calibration"])
    require(theta.overlaps(saved_theta), "production theta misses independent saved calibration")
    transition_projection = hardy_projection(theta, transition)
    saved_transition_projection = arb(joined["joined_transition_Hardy_projection_ball"])
    require(
        transition_projection.overlaps(saved_transition_projection),
        "transition Hardy projection sign or phase drift",
    )

    packet_rows = [
        ("V_L", lower_cell, "+1", "lower finite cell inside U_unowned"),
        ("C_U", upper_arc, "+1", "compact upper arc inside U_unowned"),
        ("E_s_T_U", upper_tail, "+1", "positive-real upper-tail error disk"),
        ("O_join", ordinary, "+1", "complete Gamma-subtracted ordinary packet"),
        (
            "T_A_join",
            transition,
            "+1",
            "transition source minus transition Gamma lift minus grouped natural A lift",
        ),
    ]
    require(len({row[0] for row in packet_rows}) == len(packet_rows), "duplicate top-level packet")
    complete_argument = sum((row[1] for row in packet_rows), acb(0))
    complete_HJZ = hardy_projection(theta, complete_argument)

    H = arb(deps["exact_H"]["exact_H_certificate"]["high_precision_H_ball"])
    require(H.lower() > 0, "H positivity lost")
    J_Z = complete_HJZ / H

    A_transition = arb(
        deps["A_transition"]["certificate"]["components"][
            "projector_completed_exact_A_transition_ball"
        ]
    )
    R_KGamma = J_Z + A_transition
    gamma_correction_bound = abs(
        arb(
            deps["Gamma_correction"]["correction_certificate"][
                "Gamma_to_classical_physical_absolute_bound_ball"
            ]
        )
    ).upper()
    G_minus_T = arb(0, gamma_correction_bound)
    bridge = deps["QK_bridge"]["bridge_target_certificate"]
    T = arb(bridge["classical_upper_main_ball"])
    U = arb(bridge["exact_complementary_upper_ball"])
    rho = arb(bridge["rho_RS_ball"])
    require((U - T).overlaps(rho), "rho_RS source identity drift")
    QK_minus_T = R_KGamma + G_minus_T
    Delta_KU = QK_minus_T - rho
    D_K = (1 - H) * U - H * Delta_KU
    direct_D_K = U - H * (T + QK_minus_T)
    require(D_K.overlaps(direct_D_K), "two exact D_K assemblies do not overlap")

    require(J_Z.lower() > jz_absolute_target, "J_Z no longer rigorously misses absolute target")
    require(
        R_KGamma.lower() > r_kgamma_upper_target,
        "R_KGamma no longer rigorously misses stronger-route upper wall",
    )
    require(D_K.upper() < dk_lower_target, "D_K no longer rigorously misses corridor lower wall")
    require(D_K.lower() < dk_upper_target, "D_K interval unexpectedly moved above whole corridor")

    sign_flip_rows: list[dict[str, Any]] = []
    for name, value, _, _ in packet_rows:
        if name == "E_s_T_U":
            continue
        flipped = hardy_projection(theta, complete_argument - 2 * value) / H
        separation = disjoint_separation(J_Z, flipped)
        require(separation > 0, f"single sign flip not detected for {name}")
        sign_flip_rows.append(
            {
                "packet": name,
                "flipped_J_Z_ball": interval_base.arb_record(flipped),
                "separation_from_certified_J_Z_ball": interval_base.arb_record(separation),
            }
        )

    duplicated_A = hardy_projection(theta, complete_argument + natural_A) / H
    doubled_subtraction_A = hardy_projection(theta, complete_argument - natural_A) / H
    require(
        disjoint_separation(J_Z, duplicated_A) > 0
        and disjoint_separation(J_Z, doubled_subtraction_A) > 0,
        "natural A ownership duplication guard failed",
    )

    projection_rows = []
    for name, value, coefficient, role in packet_rows:
        projection = hardy_projection(theta, value)
        projection_rows.append(
            {
                "packet": name,
                "top_level_coefficient": coefficient,
                "ownership_role": role,
                "complex_ball": interval_base.acb_record(value),
                "Hardy_projection_ball": interval_base.arb_record(projection),
                "Hardy_projection_radius_ball": interval_base.arb_record(projection.rad()),
            }
        )
    dominant_uncertainty = max(
        projection_rows,
        key=lambda row: arb(row["Hardy_projection_radius_ball"]["ball"]).upper(),
    )["packet"]

    return {
        "passed": True,
        "precision_bits": PRECISION_BITS,
        "height": HEIGHT.str(30, more=True),
        "theta_ball": interval_base.arb_record(theta),
        "saved_theta_calibration_ball": interval_base.arb_record(saved_theta),
        "H_ball": interval_base.arb_record(H),
        "exact_identity_chain": {
            "cell_partition": "H(t)J_Z=Hardy_t[U_unowned+O_join+T_A^join]",
            "unowned_reduction": "U_unowned=V_L+C_U+E_s*T_U",
            "ordinary_packet": "O_join=-T_lower-T_upper+Gamma_defect",
            "transition_packet": (
                "T_A^join=transition_source_integral-transition_Gamma_lift-"
                "grouped_natural_A_lift"
            ),
            "physical_residual": "R_KGamma=J_Z+A_transition=Q_K-G",
            "classical_bridge": "Q_K-T=R_KGamma+(G-T)",
            "finite_roster_defect": "Delta_KU=(Q_K-T)-rho_RS",
            "continuation_defect": "D_K=(1-H)U-H*Delta_KU=U-HQ_K",
        },
        "ownership_ledger": projection_rows,
        "nested_ownership_guards": {
            "grouped_natural_A_lift_occurrences": 1,
            "grouped_natural_A_lift_location": "subtracted inside T_A_join only",
            "standalone_top_level_natural_A_lift_added": False,
            "ordinary_reconstruction_ball": interval_base.acb_record(reconstructed_ordinary),
            "transition_reconstruction_ball": interval_base.acb_record(reconstructed_transition),
            "transition_projection_reconstruction_ball": interval_base.arb_record(
                transition_projection
            ),
            "duplicated_natural_A_J_Z_ball": interval_base.arb_record(duplicated_A),
            "doubled_subtraction_natural_A_J_Z_ball": interval_base.arb_record(
                doubled_subtraction_A
            ),
            "single_top_level_sign_flip_guards": sign_flip_rows,
        },
        "upper_tail_log10_upper_ball": interval_base.arb_record(tail_log10),
        "upper_tail_component_radius_ball": interval_base.arb_record(tail_radius),
        "ordinary_Gamma_defect_log10_upper_ball": interval_base.arb_record(
            ordinary_gamma_log10
        ),
        "ordinary_Gamma_defect_component_radius_ball": interval_base.arb_record(
            ordinary_gamma_radius
        ),
        "complete_preprojection_argument_ball": interval_base.acb_record(complete_argument),
        "complete_H_times_J_Z_ball": interval_base.arb_record(complete_HJZ),
        "J_Z_ball": interval_base.arb_record(J_Z),
        "A_transition_ball": interval_base.arb_record(A_transition),
        "R_KGamma_ball": interval_base.arb_record(R_KGamma),
        "G_minus_T_absolute_bound_ball": interval_base.arb_record(gamma_correction_bound),
        "Q_K_minus_T_ball": interval_base.arb_record(QK_minus_T),
        "rho_RS_ball": interval_base.arb_record(rho),
        "Delta_KU_ball": interval_base.arb_record(Delta_KU),
        "D_K_ball": interval_base.arb_record(D_K),
        "direct_D_K_ball": interval_base.arb_record(direct_D_K),
        "target_audit": {
            "absolute_J_Z_target": "0.02939975",
            "absolute_J_Z_target_holds": False,
            "J_Z_upper_wall_miss_margin_ball": interval_base.arb_record(
                J_Z.lower() - jz_absolute_target
            ),
            "R_KGamma_sufficient_corridor": "-0.0663408<R_KGamma<-0.00727437",
            "R_KGamma_sufficient_corridor_holds": False,
            "R_KGamma_upper_wall_miss_margin_ball": interval_base.arb_record(
                R_KGamma.lower() - r_kgamma_upper_target
            ),
            "D_K_sufficient_corridor": "0.0072743687<D_K<0.0663407986",
            "D_K_sufficient_corridor_holds": False,
            "D_K_lower_wall_miss_margin_ball": interval_base.arb_record(
                dk_lower_target - D_K.upper()
            ),
            "failure_due_to_enclosure_width": False,
            "dominant_projection_uncertainty_packet": dominant_uncertainty,
        },
    }


def note_text(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    target = c["target_audit"]
    rows = c["ownership_ledger"]
    row_text = "\n".join(
        f"{row['packet']:8s}  {row['Hardy_projection_ball']['ball']}"
        for row in rows
    )
    return f"""# Complete signed production-height J_Z and D_K ownership gate

Date: 2026-08-28

Status: the complete saved-height packet, `J_Z`, and `D_K` are rigorously
enclosed.  The previously selected sufficient corridor is rigorously missed;
this is a route result, not a proof of RH.

The exact preprojection ownership identity is

```text
H(t)J_Z=Hardy_t[V_L+C_U+E_s*T_U+O_join+T_A^join],
Hardy_t[X]=2 Re[exp(i*theta(t))X].                    (FJ1)
```

Here `theta(t)=Im log Gamma(1/4+it/2)-t log(pi)/2`, evaluated directly by
Arb and overlapped with the independent saved Hardy calibration.  The natural
A lift is not a sixth top-level packet: it is subtracted exactly once inside

```text
T_A^join=transition source-transition Gamma lift-mathcal A_A^nat. (FJ2)
```

The ordinary packet is independently reconstructed as
`-T_lower-T_upper+Gamma_defect`, and the transition packet is independently
reconstructed from all three terms in (FJ2).  Every non-negligible top-level
single-sign flip is disjoint from the certified result; adding the natural A
lift again or subtracting it twice is also rejected.

The retained Hardy projections are

```text
{row_text}
```

Their signed sum gives

```text
J_arg = {c['complete_preprojection_argument_ball']['real_ball']}
       +i{c['complete_preprojection_argument_ball']['imag_ball']},

H J_Z = {c['complete_H_times_J_Z_ball']['ball']},
J_Z   = {c['J_Z_ball']['ball']}.                    (FJ3)
```

Using the independently certified physical `A_transition`, exact `H`,
`rho_RS=U-T`, and the rigorous absolute Gamma-to-classical correction gives

```text
R_KGamma = {c['R_KGamma_ball']['ball']},
Q_K-T    = {c['Q_K_minus_T_ball']['ball']},
Delta_KU = {c['Delta_KU_ball']['ball']},
D_K      = {c['D_K_ball']['ball']}.                 (FJ4)
```

The two algebraically exact evaluations
`D_K=(1-H)U-H*Delta_KU` and `D_K=U-H[T+(Q_K-T)]` overlap.

The stronger sufficient route does not close:

```text
|J_Z|<0.02939975                         false,
-0.0663408<R_KGamma<-0.00727437          false,
0.0072743687<D_K<0.0663407986            false.     (FJ5)
```

The rigorous misses of the relevant walls are respectively
`{target['J_Z_upper_wall_miss_margin_ball']['ball']}`,
`{target['R_KGamma_upper_wall_miss_margin_ball']['ball']}`, and
`{target['D_K_lower_wall_miss_margin_ball']['ball']}`.  The dominant interval
uncertainty is the `{target['dominant_projection_uncertainty_packet']}` row,
but its width is far smaller than the miss.  This is not a failure due to
enclosure width.  Sharpening the same numerical
enclosures cannot reverse (FJ5); the next route must use additional signed
structure or a different sufficient theorem.

Pi provenance: every `pi` comes from the inherited Fresnel/Fourier phase,
the standard Riemann--Siegel theta normalization, the Mellin quarter turn, or
`p=t/(2*pi)`.  No circle, polygon, fitted constant, or visual pattern supplies
`pi`.

Proof boundary: rigorous final signed assembly and scalar enclosures at the
single height `t=10^10`.  This rejects only the displayed sufficient corridor.
It does not disprove the original one-sided route, prove the non-A bound,
certify equation (4) as a global infinite-series identity, give an all-height
theorem, prove `Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    resource_mode = interval_base.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = PRECISION_BITS
    ctx.threads = 1
    for path in (*DEPENDENCIES.values(), BUILDER, CHECKER):
        require(path.is_file(), f"missing source or dependency: {path}")
    deps = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(dep.get("passed") is True for dep in deps.values()), "dependency failure")

    certificate = build_certificate(deps)
    artifact = {
        "kind": "rh_c1_hardy_equation4_complete_signed_JZ_DK_ownership_gate",
        "date": "2026-08-28",
        "status": "complete_signed_production_height_JZ_and_DK_enclosed_but_sufficient_corridor_rigorously_missed",
        "passed": True,
        "resource_mode": resource_mode,
        "certificate": certificate,
        "decision": {
            "complete_top_level_packet_ownership_certified": True,
            "natural_A_lift_counted_exactly_once_inside_transition": True,
            "complete_production_height_J_Z_enclosed": True,
            "complete_production_height_D_K_enclosed": True,
            "absolute_J_Z_sufficient_target_proved": False,
            "R_KGamma_sufficient_corridor_proved": False,
            "D_K_sufficient_corridor_proved": False,
            "same_enclosure_sharpening_can_close_displayed_corridor": False,
            "non_A_bound_proved": False,
            "all_height_theorem_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Treat the failed saved-height sufficient corridor as a route diagnostic, not an RH disproof. "
            "Return to the exact signed formula R_nonA=J_Z-E_Btr,win-E_outer-I_(A,42), retain the signs "
            "of the B trace, outer term, and A integral, and determine whether the original one-sided "
            "criterion can still close. In parallel, audit which all-height theorem would transport the "
            "finite split-contour identity without importing the heuristic A21 continuation."
        ),
        "proof_boundary": (
            "Rigorous ownership, complex packet, J_Z, R_KGamma, Q_K-T, Delta_KU, and D_K enclosures "
            "at t=10^10 only. The displayed sufficient corridor is rigorously missed. This neither "
            "disproves the original one-sided route nor proves the non-A bound, equation-(4) global "
            "infinite-series exactness, an all-height theorem, Lambda<=0, PF-infinity, RH, or a "
            "prize-level conclusion."
        ),
    }
    interval_base.atomic_write(RESULT, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    interval_base.atomic_write(NOTE, note_text(artifact))
    print("certified complete signed J_Z and D_K ownership; sufficient corridor rigorously missed", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
