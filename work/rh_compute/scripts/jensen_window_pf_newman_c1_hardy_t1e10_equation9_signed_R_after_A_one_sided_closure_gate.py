#!/usr/bin/env python3
"""Close the original signed one-sided residual at the saved height."""

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

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate as interval_base


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_one_sided_closure_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)


def result_path(stem: str) -> Path:
    return ROOT / "work" / "rh_compute" / "results" / f"{stem}.json"


DEPENDENCIES = {
    "signed_JZ": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_complete_signed_JZ_DK_ownership_gate"
    ),
    "physical_ownership": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_physical_transform_ownership_ledger_gate"
    ),
    "A_face_embedding": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_fixed_regulator_embedding_gate"
    ),
    "R_Dir_ownership": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_Dir_projector_ownership_ledger_gate"
    ),
    "Gamma_insertion": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_uniform_gamma_target_insertion_gate"
    ),
    "A_transition": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_complete_endpoint_projector_assembly_gate"
    ),
    "B_local": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_outer_safe_exact_current_gate"
    ),
    "B_window": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_complete_gaussian_window_signed_gate"
    ),
    "B_outer": result_path(
        "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_phase_coupled_tangential_gate"
    ),
}

PRECISION_BITS = 384
HEIGHT = "10000000000"
R_AFTER_WORKING_TARGET = arb("0.0368147039947")
R_AFTER_NEGATIVE_TARGET = arb("0.0368061039947")
R_DIR_WORKING_TARGET = arb("0.00014057919999999995")
R_DIR_NEGATIVE_TARGET = arb("0.00013197919999999995")
QKT_WORKING_TARGET = arb("0.0000086")


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


def signed_b_trace(deps: dict[str, dict[str, Any]]) -> tuple[arb, dict[str, Any]]:
    local = deps["B_local"]["interval_certificate"]
    window = deps["B_window"]["interval_certificate"]
    principal = arb(local["combined_principal_physical_signed_projection_ball"])
    signed_first = arb(window["signed_nonlocal_first_current_ball"])
    allowances = {
        "local_positive_normal_remainder": arb(
            local["positive_normal_remainder_absolute_allowance_ball"]
        ).upper(),
        "local_negative_partner_remainder": arb(
            local["negative_partner_remainder_absolute_allowance_ball"]
        ).upper(),
        "nonlocal_higher_currents": arb(
            window["nonlocal_higher_current_absolute_allowance_ball"]
        ).upper(),
        "nonlocal_post_third_current": arb(
            window["nonlocal_post_third_current_absolute_allowance_ball"]
        ).upper(),
        "Fresnel_boundary_dictionary": arb(
            window["Fresnel_boundary_dictionary_absolute_allowance_ball"]
        ).upper(),
    }
    allowance = sum(allowances.values(), arb(0))
    result = principal + signed_first + arb(0, allowance)
    require(result.lower() > arb("-0.000134948"), "B trace lower guard failed")
    require(result.upper() < arb("-0.00013198"), "B trace upper guard failed")
    return result, {
        "construction": (
            "signed local principal plus signed nonlocal first current, with every "
            "remaining certified B allowance inserted as one symmetric radius"
        ),
        "local_principal_ball": interval_base.arb_record(principal),
        "nonlocal_signed_first_ball": interval_base.arb_record(signed_first),
        "absolute_allowances": {
            name: interval_base.arb_record(value) for name, value in allowances.items()
        },
        "absolute_allowance_sum_ball": interval_base.arb_record(allowance),
        "complete_B_trace_ball": interval_base.arb_record(result),
        "safe_enclosure": "-0.000134948<E_Btr,win<-0.00013198",
    }


def coefficient_certificate() -> dict[str, Any]:
    r_non_a = {"J_Z": 1, "E_Btr,win": -1, "E_outer": -1, "I_(A,42)": -1}
    r_after_a = dict(r_non_a)
    r_after_a["I_(A,42)"] += 1
    require(
        r_after_a == {"J_Z": 1, "E_Btr,win": -1, "E_outer": -1, "I_(A,42)": 0},
        "A-face cancellation coefficient drift",
    )
    route_two_left = {"J_Z": 1, "A_transition": 1}
    route_two_right_fixed = {"E_Btr,win": 1, "E_outer": 1, "A_transition": 1}
    solved = {
        "J_Z": route_two_left["J_Z"],
        "E_Btr,win": -route_two_right_fixed["E_Btr,win"],
        "E_outer": -route_two_right_fixed["E_outer"],
        "A_transition": route_two_left["A_transition"]
        - route_two_right_fixed["A_transition"],
    }
    require(
        solved == {"J_Z": 1, "E_Btr,win": -1, "E_outer": -1, "A_transition": 0},
        "independent R_KGamma cancellation drift",
    )
    return {
        "route_one_inputs": [
            "R_after_A=I_(A,42)+R_nonA",
            "R_nonA=J_Z-E_Btr,win-E_outer-I_(A,42)",
        ],
        "route_one_collected_coefficients": r_after_a,
        "route_two_inputs": [
            "R_KGamma=J_Z+A_transition",
            "R_KGamma=E_Btr,win+E_outer+A_transition+R_after_A",
        ],
        "route_two_solved_coefficients": solved,
        "common_exact_result": "R_after_A=J_Z-E_Btr,win-E_outer",
        "A_face_coefficient_after_collection": 0,
        "A_transition_coefficient_after_collection": 0,
    }


def build_certificate(deps: dict[str, dict[str, Any]]) -> dict[str, Any]:
    signed = deps["signed_JZ"]["certificate"]
    J_Z = arb(signed["J_Z_ball"]["ball"])
    saved_R_KGamma = arb(signed["R_KGamma_ball"]["ball"])
    saved_QKT = arb(signed["Q_K_minus_T_ball"]["ball"])
    B_trace, B_trace_record = signed_b_trace(deps)

    outer_cert = deps["B_outer"]["interval_certificate"]
    outer_bound = arb(outer_cert["physical_outer_bound_ball"]).upper()
    require(outer_bound < arb("8e-10"), "outer bound lost 8e-10 guard")
    E_outer = arb(0, outer_bound)

    A_transition = arb(
        deps["A_transition"]["certificate"]["components"][
            "projector_completed_exact_A_transition_ball"
        ]
    )
    require(
        A_transition.overlaps(arb(signed["A_transition_ball"]["ball"])),
        "A-transition dependency drift",
    )
    R_KGamma = J_Z + A_transition
    require(R_KGamma.overlaps(saved_R_KGamma), "J_Z+A_transition misses saved R_KGamma")

    R_after_A = J_Z - B_trace - E_outer
    R_Dir = A_transition + R_after_A
    R_Dir_alternate = saved_R_KGamma - B_trace - E_outer
    require(R_Dir.overlaps(R_Dir_alternate), "two R_Dir constructions do not overlap")

    gamma_bound = abs(
        arb(
            deps["Gamma_insertion"]["correction_certificate"][
                "Gamma_to_classical_physical_absolute_bound_ball"
            ]
        )
    ).upper()
    require(gamma_bound < arb("5e-20"), "Gamma-to-classical guard failed")
    G_minus_T = arb(0, gamma_bound)
    QKT = saved_R_KGamma + G_minus_T
    QKT_rejoined = B_trace + E_outer + R_Dir + G_minus_T
    require(QKT.overlaps(saved_QKT), "reconstructed Q_K-T misses saved enclosure")
    require(QKT_rejoined.overlaps(saved_QKT), "rejoined Q_K-T misses saved enclosure")

    require(R_after_A.upper() < R_AFTER_WORKING_TARGET, "working R_after_A target failed")
    require(
        R_after_A.upper() < R_AFTER_NEGATIVE_TARGET,
        "stricter R_after_A negativity target failed",
    )
    require(R_Dir.upper() < R_DIR_WORKING_TARGET, "working R_Dir target failed")
    require(R_Dir.upper() < R_DIR_NEGATIVE_TARGET, "negative-route R_Dir target failed")
    require(R_Dir.upper() < 0, "R_Dir itself is not certified negative")
    require(QKT.upper() < QKT_WORKING_TARGET, "Q_K-T working target failed")
    require(QKT.upper() < 0, "Q_K-T negativity failed")
    require(QKT_rejoined.upper() < 0, "rejoined Q_K-T negativity failed")

    return {
        "passed": True,
        "precision_bits": PRECISION_BITS,
        "height": HEIGHT,
        "coefficient_cancellation": coefficient_certificate(),
        "B_trace_certificate": B_trace_record,
        "J_Z_ball": interval_base.arb_record(J_Z),
        "B_trace_ball": interval_base.arb_record(B_trace),
        "E_outer_absolute_bound_ball": interval_base.arb_record(outer_bound),
        "E_outer_symmetric_ball": interval_base.arb_record(E_outer),
        "A_transition_ball": interval_base.arb_record(A_transition),
        "R_after_A_ball": interval_base.arb_record(R_after_A),
        "R_Dir_ball": interval_base.arb_record(R_Dir),
        "R_Dir_alternate_ball": interval_base.arb_record(R_Dir_alternate),
        "R_KGamma_ball": interval_base.arb_record(saved_R_KGamma),
        "Gamma_to_classical_absolute_bound_ball": interval_base.arb_record(gamma_bound),
        "Q_K_minus_T_ball": interval_base.arb_record(QKT),
        "Q_K_minus_T_rejoined_ball": interval_base.arb_record(QKT_rejoined),
        "target_audit": {
            "R_after_A_working_target": R_AFTER_WORKING_TARGET.str(30, more=True),
            "R_after_A_working_margin_ball": interval_base.arb_record(
                R_AFTER_WORKING_TARGET - R_after_A.upper()
            ),
            "R_after_A_negative_route_target": R_AFTER_NEGATIVE_TARGET.str(30, more=True),
            "R_after_A_negative_route_margin_ball": interval_base.arb_record(
                R_AFTER_NEGATIVE_TARGET - R_after_A.upper()
            ),
            "R_Dir_working_target": R_DIR_WORKING_TARGET.str(30, more=True),
            "R_Dir_working_margin_ball": interval_base.arb_record(
                R_DIR_WORKING_TARGET - R_Dir.upper()
            ),
            "R_Dir_negative_route_target": R_DIR_NEGATIVE_TARGET.str(30, more=True),
            "R_Dir_negative_route_margin_ball": interval_base.arb_record(
                R_DIR_NEGATIVE_TARGET - R_Dir.upper()
            ),
            "Q_K_minus_T_working_target": QKT_WORKING_TARGET.str(20, more=True),
            "Q_K_minus_T_working_margin_ball": interval_base.arb_record(
                QKT_WORKING_TARGET - QKT.upper()
            ),
            "Q_K_minus_T_negativity_margin_ball": interval_base.arb_record(-QKT.upper()),
            "rejoined_Q_K_minus_T_negativity_margin_ball": interval_base.arb_record(
                -QKT_rejoined.upper()
            ),
        },
    }


def note_text(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    t = c["target_audit"]
    return f"""# Signed saved-height one-sided residual closure gate

Date: 2026-08-28

Status: interval certificate for the original saved-height one-sided residual;
not an all-height proof, not a proof of `Lambda<=0`, and not an RH proof.

Two already certified exact identities give

```text
R_after_A=I_(A,42)+R_nonA,
R_nonA=J_Z-E_Btr,win-E_outer-I_(A,42).
```

Collecting coefficients before any interval norm cancels the translated
A-face channel exactly:

```text
R_after_A=J_Z-E_Btr,win-E_outer.                    (SC1)
```

The same identity follows independently from

```text
R_KGamma=J_Z+A_transition,
R_KGamma=E_Btr,win+E_outer+A_transition+R_after_A. (SC2)
```

Thus both `I_(A,42)` and `A_transition` have coefficient zero in the final
saved-height `R_after_A` formula.  No sign or quadrature estimate for the
translated A-face is needed for this closure.

The signed B trace is rebuilt from the raw local principal, signed nonlocal
first current, and every remaining certified absolute allowance:

```text
E_Btr,win={c['B_trace_ball']['ball']},
-0.000134948<E_Btr,win<-0.00013198.                (SC3)
```

Together with

```text
J_Z={c['J_Z_ball']['ball']},
|E_outer|< {c['E_outer_absolute_bound_ball']['ball']},
```

outward-rounded Arb arithmetic gives

```text
R_after_A={c['R_after_A_ball']['ball']}.            (SC4)
```

Its upper endpoint is below both `0.0368147039947` and
`0.0368061039947`.  The corresponding rigorous margins are
`{t['R_after_A_working_margin_ball']['ball']}` and
`{t['R_after_A_negative_route_margin_ball']['ball']}`.

Restoring the independently certified projector-completed A transition gives

```text
R_Dir={c['R_Dir_ball']['ball']}.                    (SC5)
```

This is below zero, hence below both inherited `R_Dir` sufficient targets.
The alternate construction `R_KGamma-E_Btr,win-E_outer` overlaps (SC5).

Finally, the sub-`5e-20` Gamma-to-classical correction gives

```text
Q_K-T={c['Q_K_minus_T_ball']['ball']}.              (SC6)
```

Its upper endpoint is negative with margin
`{t['Q_K_minus_T_negativity_margin_ball']['ball']}` and therefore is also
below the working `8.6e-6` threshold.  Rejoining the B trace, outer ball,
`R_Dir`, and Gamma correction with dependency loss still remains negative;
its margin is `{t['rejoined_Q_K_minus_T_negativity_margin_ball']['ball']}`.

Pi provenance: this gate introduces no new occurrence of `pi`.  All inherited
occurrences belong to the already certified Fresnel, Fourier, Kummer, Gamma,
or Riemann--Siegel normalizations.

Proof boundary: this certifies the original signed one-sided residual and
`Q_K-T<0` only at the single saved height `t=10^10`, for the exact finite
roster and common-regulator identities pinned here.  It does not prove the
stronger absolute non-A bound, an all-height finite split-contour transport
theorem, equation-(4) as a global infinite-series identity, `Lambda<=0`,
PF-infinity, RH, or a prize-level conclusion.
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
    require(
        deps["A_face_embedding"]["decision"]["physical_Abel_zero_limit_split_proved"] is True,
        "A-face physical split drift",
    )
    require(
        deps["R_Dir_ownership"]["decision"]["post_A_one_sided_targets_derived"] is True,
        "R_Dir target provenance drift",
    )
    require(
        deps["Gamma_insertion"]["decision"]["Gamma_residual_identified_as_QK_minus_G"] is True,
        "Gamma residual identity drift",
    )
    require(
        deps["signed_JZ"]["decision"]["complete_production_height_J_Z_enclosed"] is True,
        "saved J_Z enclosure drift",
    )

    certificate = build_certificate(deps)
    artifact = {
        "kind": "rh_c1_hardy_equation9_signed_saved_height_one_sided_residual_closure_gate",
        "date": "2026-08-28",
        "status": "saved_height_original_signed_one_sided_residual_and_QK_minus_T_negativity_certified",
        "passed": True,
        "resource_mode": resource_mode,
        "certificate": certificate,
        "decision": {
            "translated_A_face_cancels_exactly_from_R_after_A": True,
            "translated_A_face_sign_or_new_quadrature_needed_for_closure": False,
            "R_after_A_working_upper_target_proved_at_saved_height": True,
            "R_after_A_stricter_upper_target_proved_at_saved_height": True,
            "R_Dir_working_target_proved_at_saved_height": True,
            "R_Dir_negative_at_saved_height": True,
            "Q_K_minus_T_below_8_point_6e_minus_6_at_saved_height": True,
            "Q_K_minus_T_negative_at_saved_height": True,
            "original_saved_height_one_sided_route_closed": True,
            "stronger_absolute_non_A_route_proved": False,
            "all_height_transport_theorem_proved": False,
            "equation4_global_infinite_series_identity_proved": False,
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
            "Promote the exact finite split-contour identities from the isolated height t=10^10 "
            "to a nontrivial height interval. Define the event-cell partition and prove uniform "
            "crossing, endpoint, and remainder transport without importing the quarantined A21 "
            "infinite interchange. Keep the stronger absolute non-A route retired unless a later "
            "theorem genuinely requires it."
        ),
        "proof_boundary": (
            "Rigorous original one-sided residual closure and Q_K-T negativity at t=10^10 only. "
            "No stronger absolute non-A estimate, all-height transport theorem, equation-(4) global "
            "infinite-series identity, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
        ),
    }
    interval_base.atomic_write(RESULT, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    interval_base.atomic_write(NOTE, note_text(artifact))
    print(
        "certified saved-height signed one-sided closure and Q_K-T negativity; "
        f"R_after_A={arb(certificate['R_after_A_ball']['ball']).str(12, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
