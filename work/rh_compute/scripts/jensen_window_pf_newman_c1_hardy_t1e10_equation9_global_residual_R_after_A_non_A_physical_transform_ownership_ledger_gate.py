#!/usr/bin/env python3
"""Certify the physical-transform ownership ledger for the non-A residual."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import arb, ctx
import sympy as sp


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_physical_transform_ownership_ledger_gate"
)
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "source_interchange": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json",
    "source_nonidentification": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_endpoint_residual_source_hybrid_nonidentification_guard_gate.json",
    "finite_regulator": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_finite_regulator_equivalence_gate.json",
    "Gamma_insertion": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_uniform_gamma_target_insertion_gate.json",
    "A_transition": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_complete_endpoint_projector_assembly_gate.json",
    "B_window": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_complete_gaussian_window_signed_gate.json",
    "B_local": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_outer_safe_exact_current_gate.json",
    "B_outer": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_phase_coupled_tangential_gate.json",
    "A_face": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_complete_absolute_assembly_gate.json",
    "A_face_embedding": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_fixed_regulator_embedding_gate.json",
    "source_reassembly": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_folded_abel_source_roster_reassembly_gate.json",
    "sparse_pilot": REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_source_owned_carrier_sparse_pilot_gate.json",
}

HEIGHT = 10_000_000_000
A = 159_577
B = 5_122_421
L = 2_481_422
TARGET_FIRST = 622
TARGET_LAST = 39_894
EXTENDED_LAST = 39_936
A_PAIR_FIRST = 39_853
NON_A_TARGET_TEXT = "0.0331747039947"
B_OUTER_BOUND_TEXT = "8e-10"
A_FACE_BOUND_TEXT = "0.00364"
GAMMA_CORRECTION_BOUND_TEXT = "5e-20"
PUBLISHED_JZ_LOWER_TEXT = "-0.02966668"
PUBLISHED_JZ_UPPER_TEXT = "0.02939975"
PUBLISHED_JZ_ABSOLUTE_TEXT = "0.02939975"
PUBLISHED_RKGAMMA_LOWER_TEXT = "-0.0663408"
PUBLISHED_RKGAMMA_UPPER_TEXT = "-0.00727437"
PRECISION = 110


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {relative(path)}")
    return json.loads(path.read_text(encoding="utf-8"))


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def source_transform_certificate() -> dict[str, str]:
    t = sp.symbols("t", positive=True, real=True)
    a = sp.Rational(3, 4) - sp.I * t / 2
    b = sp.Rational(3, 4) + sp.I * t / 2
    require(sp.simplify(a + b - sp.Rational(3, 2)) == 0, "Kummer denominator drift")
    require(
        sp.simplify((a - 1) - (-sp.Rational(1, 4) - sp.I * t / 2)) == 0,
        "lower endpoint exponent drift",
    )
    require(
        sp.simplify((b - 1) - (-sp.Rational(1, 4) + sp.I * t / 2)) == 0,
        "upper endpoint exponent drift",
    )
    return {
        "physical_functional": "P_t[F]=2(pi/(32t))^(1/4) Re[e^(-i*pi/8) integral_0^1 W_t(x)F(x)dx]",
        "weight": "W_t(x)=x^(-1/4-it/2)(1-x)^(-1/4+it/2)",
        "parameters": "a=3/4-it/2, b=3/4+it/2, a+b=3/2",
        "Euler_Kummer_identity": "integral_0^1 x^(a-1)(1-x)^(b-1)e^(zx)dx=B(a,b) 1F1(a;3/2;z)",
        "one_source_label": "K_t(alpha)=2(pi/(32t))^(1/4) Re[e^(-i*pi/8) alpha B(a,b) 1F1(a;3/2;i*pi*alpha^2/4)]",
        "complete_source_transform": "Q_K=sum_(alpha=A,A+2,...,B) K_t(alpha)=P_t[sum_(n=0)^L f_x(n)]",
        "roster": "A=159577, B=5122421, L=2481422, labels=L+1=2481423",
    }


def identity_certificate() -> dict[str, str]:
    qk, g, a_transition, b_trace, b_outer, a_face = sp.symbols(
        "Q_K G A_transition E_Btr E_outer I_A42"
    )
    j_z = qk - g - a_transition
    r_non_a = j_z - b_trace - b_outer - a_face
    expanded = qk - g - a_transition - b_trace - b_outer - a_face
    require(sp.expand(r_non_a - expanded) == 0, "physical ownership identity failed")
    return {
        "target_Gamma_residual": "R_KGamma=Q_K-G",
        "extended_P_block": "P_t[sum_(m=622)^39936 P_m]=G+G_extra, where G_extra uses 39895..39936",
        "paired_A_block": "P_t[sum_(m=39853)^39936(A_m+A_-m)]=A_endpoint",
        "overlap_reassembly": "A_transition=A_endpoint+G_extra",
        "joined_source_owned_transform": "J_Z=P_t[Z]=Q_K-G-A_transition=R_KGamma-A_transition",
        "non_A_identity": "R_nonA=J_Z-E_Btr,win-E_outer-I_(A,42)",
        "expanded_non_A_identity": "R_nonA=Q_K-G-A_transition-E_Btr,win-E_outer-I_(A,42)",
        "sign_provenance": "The three final minus signs are inherited from -chi_W Btr-O-mathfrak E_(A,42) in the common non-A kernel.",
    }


def b_trace_interval(dependencies: dict[str, dict[str, Any]]) -> dict[str, Any]:
    local = dependencies["B_local"]["interval_certificate"]
    window = dependencies["B_window"]["interval_certificate"]
    principal = arb(local["combined_principal_physical_signed_projection_ball"])
    signed_first = arb(window["signed_nonlocal_first_current_ball"])
    absolute_allowances = {
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
    radius = sum(absolute_allowances.values(), arb(0))
    interval = principal + signed_first + arb(0, radius)
    require(interval.upper() < arb("-0.00013198"), "B trace lost certified negative sign")
    require(interval.lower() > arb("-0.000134948"), "B trace lower scale guard failed")
    require(abs(interval).upper() < arb("0.000134948"), "B trace absolute envelope failed")
    return {
        "construction": "local signed principal + nonlocal signed first current +/- every remaining certified absolute allowance",
        "local_principal_ball": principal.str(90, more=True),
        "nonlocal_signed_first_ball": signed_first.str(90, more=True),
        "absolute_allowances": {name: value.str(70, more=True) for name, value in absolute_allowances.items()},
        "absolute_allowance_sum_ball": radius.str(80, more=True),
        "derived_two_sided_ball": interval.str(100, more=True),
        "derived_lower_ball": interval.lower().str(90, more=True),
        "derived_upper_ball": interval.upper().str(90, more=True),
        "derived_absolute_upper_ball": abs(interval).upper().str(90, more=True),
        "published_safe_enclosure": "-0.000134948<E_Btr,win<-0.00013198",
        "stored_one_sided_upper_ball": window[
            "complete_B_trace_package_window_physical_upper_ball"
        ],
        "interval": interval,
    }


def physical_ownership(
    dependencies: dict[str, dict[str, Any]], b_trace: arb
) -> tuple[list[dict[str, Any]], arb]:
    a_components = dependencies["A_transition"]["certificate"]["components"]
    a_endpoint = arb(a_components["complete_exact_A_endpoint_carrier_ball"])
    gamma_extra = arb(a_components["outer_positive_full_line_Gamma_bulk_ball"])
    a_transition = arb(a_components["projector_completed_exact_A_transition_ball"])
    assembled = a_endpoint + gamma_extra
    require(
        assembled.lower() <= a_transition.upper() and a_transition.lower() <= assembled.upper(),
        "A endpoint/Gamma-extra reassembly drift",
    )
    gamma = dependencies["Gamma_insertion"]["correction_certificate"]
    return (
        [
            {
                "channel": "complete_source",
                "pointwise_owner": "sum_(n=0)^L f_x(n)",
                "physical_owner": "Q_K exact finite Kummer roster",
                "level": "exact physical transform",
                "certified_value_or_bound": None,
                "open": True,
            },
            {
                "channel": "target_P_Gamma_block",
                "pointwise_owner": "sum_(m=622)^39894 P_m",
                "physical_owner": "G",
                "level": "exact physical transform plus certified classical correction",
                "certified_value_or_bound": gamma[
                    "Gamma_to_classical_physical_absolute_bound_ball"
                ],
                "bound_interpretation": "|G-T|<5e-20",
                "open": False,
            },
            {
                "channel": "extra_P_Gamma_block",
                "pointwise_owner": "sum_(m=39895)^39936 P_m",
                "physical_owner": "G_extra",
                "level": "signed physical interval",
                "certified_value_or_bound": gamma_extra.str(90, more=True),
                "open": False,
            },
            {
                "channel": "paired_A_endpoint_block",
                "pointwise_owner": "sum_(m=39853)^39936(A_m+A_-m)",
                "physical_owner": "A_endpoint",
                "level": "signed physical interval",
                "certified_value_or_bound": a_endpoint.str(90, more=True),
                "open": False,
            },
            {
                "channel": "projector_completed_A_transition",
                "pointwise_owner": "paired A endpoint block plus extra P/Gamma block",
                "physical_owner": "A_transition=A_endpoint+G_extra",
                "level": "signed physical interval",
                "certified_value_or_bound": a_transition.str(90, more=True),
                "open": False,
            },
            {
                "channel": "B_trace_window",
                "pointwise_owner": "chi_W Btr in the common regulator",
                "physical_owner": "E_Btr,win",
                "level": "derived signed physical interval",
                "certified_value_or_bound": b_trace.str(90, more=True),
                "open": False,
            },
            {
                "channel": "completed_B_outer",
                "pointwise_owner": "O_(M,epsilon)",
                "physical_owner": "E_outer",
                "level": "absolute physical bound",
                "certified_value_or_bound": "|E_outer|<8e-10",
                "open": False,
            },
            {
                "channel": "translated_A_face",
                "pointwise_owner": "mathfrak E_(A,42),epsilon",
                "physical_owner": "I_(A,42)",
                "level": "absolute physical bound",
                "certified_value_or_bound": "|I_(A,42)|<0.00364",
                "open": False,
            },
            {
                "channel": "joined_source_owned_transform",
                "pointwise_owner": "Z=complete source-P_(622..39936)-paired A_(39853..39936)",
                "physical_owner": "J_Z=Q_K-G-A_transition",
                "level": "exact identity, numerical interval open",
                "certified_value_or_bound": None,
                "open": True,
            },
        ],
        a_transition,
    )


def corridor_certificate(b_trace: arb, a_transition: arb) -> dict[str, str]:
    non_a_target = arb(NON_A_TARGET_TEXT)
    b_outer_bound = arb(B_OUTER_BOUND_TEXT)
    a_face_bound = arb(A_FACE_BOUND_TEXT)
    gamma_correction_bound = arb(GAMMA_CORRECTION_BOUND_TEXT)
    published_jz_lower = arb(PUBLISHED_JZ_LOWER_TEXT)
    published_jz_upper = arb(PUBLISHED_JZ_UPPER_TEXT)
    published_jz_absolute = arb(PUBLISHED_JZ_ABSOLUTE_TEXT)
    published_rkgamma_lower = arb(PUBLISHED_RKGAMMA_LOWER_TEXT)
    published_rkgamma_upper = arb(PUBLISHED_RKGAMMA_UPPER_TEXT)

    known_after_jz = b_trace + arb(0, b_outer_bound + a_face_bound)
    jz_lower = known_after_jz.upper() - non_a_target
    jz_upper = known_after_jz.lower() + non_a_target
    require(jz_lower < jz_upper, "J_Z corridor is empty")
    require(published_jz_lower > jz_lower, "published J_Z lower endpoint is not conservative")
    require(published_jz_upper < jz_upper, "published J_Z upper endpoint is not conservative")
    require(published_jz_absolute < jz_upper, "published symmetric J_Z target is not conservative")
    require(-published_jz_absolute > jz_lower, "published symmetric lower J_Z target failed")

    extracted = a_transition + known_after_jz
    rkg_lower = extracted.upper() - non_a_target
    rkg_upper = extracted.lower() + non_a_target
    require(published_rkgamma_lower > rkg_lower, "published R_KGamma lower endpoint failed")
    require(published_rkgamma_upper < rkg_upper, "published R_KGamma upper endpoint failed")

    qkt_lower = rkg_lower + gamma_correction_bound
    qkt_upper = rkg_upper - gamma_correction_bound
    require(published_rkgamma_lower > qkt_lower, "published Q_K-T lower endpoint failed")
    require(published_rkgamma_upper < qkt_upper, "published Q_K-T upper endpoint failed")
    return {
        "non_A_absolute_target": non_a_target.str(60, more=True),
        "known_post_JZ_channel_ball": known_after_jz.str(100, more=True),
        "exact_JZ_safe_lower_ball": jz_lower.str(90, more=True),
        "exact_JZ_safe_upper_ball": jz_upper.str(90, more=True),
        "published_asymmetric_JZ_sufficient_corridor": "-0.02966668<J_Z<0.02939975",
        "published_symmetric_JZ_sufficient_target": "|J_Z|<0.02939975",
        "combined_extracted_channel_ball": extracted.str(100, more=True),
        "exact_R_KGamma_safe_lower_ball": rkg_lower.str(90, more=True),
        "exact_R_KGamma_safe_upper_ball": rkg_upper.str(90, more=True),
        "published_R_KGamma_sufficient_corridor": "-0.0663408<R_KGamma<-0.00727437",
        "Gamma_to_classical_absolute_bound": "5e-20",
        "exact_Q_K_minus_T_safe_lower_ball": qkt_lower.str(90, more=True),
        "exact_Q_K_minus_T_safe_upper_ball": qkt_upper.str(90, more=True),
        "published_Q_K_minus_T_sufficient_corridor": "-0.0663408<Q_K-T<-0.00727437",
        "scope_guard": "These corridors are sufficient for the stronger absolute non-A route, not necessary conditions for the one-sided Q_K-T theorem.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    bcert = artifact["B_trace_two_sided_certificate"]
    corridor = artifact["corridor_certificate"]
    return f"""# Non-A physical-transform ownership ledger

Date: 2026-08-26

Status: exact physical ownership certified; one joined Kummer transform remains open

Use the inherited physical functional

```text
P_t[F]=2(pi/(32t))^(1/4)
 Re[e^(-i*pi/8) integral_0^1 W_t(x)F(x)dx].          (PT1)
```

For `a=3/4-it/2`, `b=3/4+it/2`, Euler's Kummer integral gives, label by
label,

```text
K_t(alpha)=2(pi/(32t))^(1/4)
 Re[e^(-i*pi/8) alpha B(a,b)
                    1F1(a;3/2;i*pi*alpha^2/4)],

Q_K=sum_(alpha=A,A+2,...,B)K_t(alpha)
   =P_t[sum_(n=0)^L f_x(n)].                          (PT2)
```

Thus the complete source has an exact physical transform: it is the finite
equation-(9) Kummer roster `Q_K`.  It does not have a certified numerical
interval here.  The implemented or published hybrid telemetry cannot be
substituted for it.

The finite carrier ownership now closes without overlap.  The target
`P` block transforms to `G`; modes `39895..39936` transform to `G_extra`;
and the paired `A` block transforms to `A_endpoint`.  The previously
certified identity `A_transition=A_endpoint+G_extra` therefore gives

```text
J_Z=P_t[Z]=Q_K-G-A_transition=R_KGamma-A_transition. (PT3)
```

The remaining native physical subtractions retain their original signs:

```text
R_nonA=J_Z-E_Btr,win-E_outer-I_(A,42).                (PT4)
```

Reassembling the signed local and nonlocal first B currents with every
remaining absolute allowance gives the new two-sided enclosure

```text
E_Btr,win={bcert['derived_two_sided_ball']},
-0.000134948<E_Btr,win<-0.00013198.                  (PT5)
```

Together with `|E_outer|<8e-10` and `|I_(A,42)|<0.00364`, exact interval
arithmetic gives the sufficient transformed target

```text
{corridor['published_asymmetric_JZ_sufficient_corridor']},

or, more simply,

{corridor['published_symmetric_JZ_sufficient_target']}.              (PT6)
```

Either statement implies `|R_nonA|<0.0331747039947`.  In terms of the
unextracted Gamma residual, the same sufficient corridor is

```text
{corridor['published_R_KGamma_sufficient_corridor']}. (PT7)
```

The sub-`5e-20` Gamma-to-classical correction leaves the same displayed
rounded corridor for `Q_K-T`.  These are sufficient conditions for the
stronger absolute non-A route, not necessary conditions for the original
one-sided upper theorem.

The selected next object is therefore the single joined transform `J_Z`,
not six separate pointwise integrals and not the source-hybrid telemetry.
An admissible next method must evaluate or enclose `J_Z` in its physical
normalization while preserving the source/Gamma/A cancellation.  Direct
event-aware `x` panels remain only a fallback.

Pi provenance: every `pi` in (PT1)--(PT7) is inherited from the equation-(9)
physical normalization, Kummer quadratic phase, Gaussian Abel regulator, or
the already certified Gamma/A/B identities.  No geometric or fitted
occurrence is introduced.

Proof boundary: exact saved-height physical-transform ownership, a derived
two-sided B-trace interval, and sufficient scalar target arithmetic only.
No numerical enclosure of `J_Z`, non-A bound, joined `R_after_A`, `R_Dir`,
`Q_K-T`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    ctx.dps = PRECISION
    ctx.threads = 1
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "dependency failure")
    require(
        dependencies["source_interchange"]["decision"]["symmetric_poisson_interchange_proved"] is True,
        "source interchange drift",
    )
    require(
        dependencies["source_nonidentification"]["decision"]["exact_endpoint_object_is_finite_Kummer_roster"]
        is True,
        "Q_K ownership drift",
    )
    require(
        dependencies["finite_regulator"]["decision"]["finite_R_after_A_defined_on_one_common_regulator"]
        is True,
        "finite regulator drift",
    )
    require(
        dependencies["Gamma_insertion"]["decision"]["Gamma_residual_identified_as_QK_minus_G"] is True,
        "Gamma residual drift",
    )
    require(
        dependencies["A_transition"]["decision"]["projector_completed_exact_A_transition_certified"]
        is True,
        "A transition drift",
    )
    require(
        dependencies["B_window"]["decision"]["all_nonlocal_B_trace_current_window_channels_covered"]
        is True,
        "B trace coverage drift",
    )
    require(
        dependencies["B_local"]["decision"]["target_step_and_Fresnel_tail_kept_joined"] is True,
        "local B grouping drift",
    )
    require(
        dependencies["B_outer"]["decision"][
            "complete_analytic_m_ge_B_exterior_physical_contribution_below_8e_minus_10"
        ]
        is True,
        "B outer drift",
    )
    require(
        dependencies["A_face"]["decision"][
            "complete_translated_42_mode_A_face_absolute_bound_below_0_point_00364_proved"
        ]
        is True,
        "A-face bound drift",
    )
    require(
        dependencies["A_face_embedding"]["decision"]["physical_Abel_zero_limit_split_proved"] is True,
        "A-face embedding drift",
    )
    require(
        dependencies["source_reassembly"]["decision"]["zero_regulator_source_block_matches_post_A_ownership"]
        is True,
        "source reassembly drift",
    )
    require(
        dependencies["sparse_pilot"]["decision"]["physical_transform_ownership_ledger_selected_next"]
        is True,
        "sparse-pilot route drift",
    )

    b_trace_record = b_trace_interval(dependencies)
    b_trace = b_trace_record.pop("interval")
    ownership, a_transition = physical_ownership(dependencies, b_trace)
    artifact = {
        "kind": STEM,
        "status": "exact_non_A_physical_transform_ownership_ledger_certified_single_joined_J_Z_transform_open",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "source_labels": [A, B, 2],
            "source_label_count": L + 1,
            "target_P_modes": [TARGET_FIRST, TARGET_LAST],
            "extended_P_modes": [TARGET_FIRST, EXTENDED_LAST],
            "paired_A_modes": [A_PAIR_FIRST, EXTENDED_LAST],
            "physical_normalization": "2(pi/(32t))^(1/4) Re[e^(-i*pi/8) integral_0^1 W_t(x)(.)dx]",
        },
        "source_transform_certificate": source_transform_certificate(),
        "identity_certificate": identity_certificate(),
        "physical_transform_ownership": ownership,
        "B_trace_two_sided_certificate": b_trace_record,
        "corridor_certificate": corridor_certificate(b_trace, a_transition),
        "decision": {
            "complete_source_transform_identified_exactly_as_Q_K": True,
            "target_and_extra_Gamma_blocks_owned_without_overlap": True,
            "paired_A_and_extra_Gamma_blocks_reassembled_as_A_transition": True,
            "B_trace_two_sided_physical_interval_derived": True,
            "B_outer_and_translated_A_face_reused_in_native_physical_normalization": True,
            "single_joined_J_Z_transform_isolated": True,
            "source_hybrid_telemetry_transfer_permitted": False,
            "naive_pointwise_panel_route_selected": False,
            "J_Z_interval_enclosed": False,
            "non_A_bound_proved": False,
            "complete_Q_K_minus_T_proved": False,
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
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "arb_threads": 1,
            "process_priority": priority,
            "precision_decimal_digits": PRECISION,
        },
        "next_obligation": "Derive a cancellation-preserving physical representation or rigorous enclosure for J_Z=P_t[complete source-P_(622..39936)-paired A_(39853..39936)]. Prefer an exact Kummer/contiguous-recurrence or common-contour transform that keeps Q_K-G-A_transition joined. Target |J_Z|<0.02939975. Do not substitute published or implemented source-hybrid telemetry for Q_K; use event-aware x panels only if the transformed route leaves a genuinely unevaluated term.",
        "proof_boundary": "Exact saved-height physical-transform ownership, a conservative two-sided B-trace interval, and sufficient scalar target arithmetic only. No numerical J_Z enclosure, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified non-A physical-transform ownership ledger and isolated J_Z", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
