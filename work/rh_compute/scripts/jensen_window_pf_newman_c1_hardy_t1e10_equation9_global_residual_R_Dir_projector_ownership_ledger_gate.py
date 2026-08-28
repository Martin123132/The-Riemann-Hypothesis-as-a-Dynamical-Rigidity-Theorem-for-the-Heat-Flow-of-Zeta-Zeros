#!/usr/bin/env python3
"""Certify the global R_Dir ownership ledger after extracting A."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_Dir_projector_ownership_ledger_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "symmetric_reassembly": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_symmetric_endpoint_tail_reassembly_gate.json",
    "double_weber_reconstruction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_double_weber_exact_source_reconstruction_gate.json",
    "B_window_extraction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_B_window_extraction_target_gate.json",
    "B_window_bound": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_complete_gaussian_window_signed_gate.json",
    "B_outer_bound": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_phase_coupled_tangential_gate.json",
    "finite_dirichlet_rejoin": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_finite_block_dirichlet_rejoin_gate.json",
    "Gamma_insertion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_uniform_gamma_target_insertion_gate.json",
    "projector_defect": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_phase_space_projector_defect_gate.json",
    "A_transition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_complete_endpoint_projector_assembly_gate.json",
}

HEIGHT = 10_000_000_000
TARGET_START = 622
TARGET_END = 39_894
A_START = 39_853
A_END = 39_936
B_OUTER_START = 5_122_421
WORKING_TARGET = "0.00014057919999999995"
NEGATIVE_TARGET = "0.00013197919999999995"
PRECISION = 100
ATOM_ORDER = ("P_plus", "P_minus", "A_plus", "A_minus", "B_plus", "B_minus")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
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


def vector(*values: int) -> dict[str, int]:
    require(len(values) == len(ATOM_ORDER), "coefficient-vector length drift")
    return dict(zip(ATOM_ORDER, values, strict=True))


def add_vectors(left: dict[str, int], right: dict[str, int]) -> dict[str, int]:
    return {atom: left[atom] + right[atom] for atom in ATOM_ORDER}


def mode_ownership_rows() -> list[dict[str, Any]]:
    zero = vector(0, 0, 0, 0, 0, 0)
    rows = [
        {
            "sector_id": "zero_mode",
            "mode_range": [0, 0],
            "count": "1",
            "target_indicator": 0,
            "pre_A_coefficients": {"I_zero": 1},
            "A_transition_coefficients": {"I_zero": 0},
            "post_A_coefficients": {"I_zero": 1},
            "post_A_owner": "joined_R_after_A",
            "separate_norm_admissible": False,
        },
        {
            "sector_id": "low_positive_pairs",
            "mode_range": [1, 621],
            "count": "621",
            "target_indicator": 0,
            "pre_A_coefficients": vector(1, 1, 1, 1, 1, 1),
            "A_transition_coefficients": zero,
            "post_A_coefficients": vector(1, 1, 1, 1, 1, 1),
            "post_A_owner": "joined_R_after_A",
            "separate_norm_admissible": False,
        },
        {
            "sector_id": "ordinary_target_pairs_before_A_window",
            "mode_range": [622, 39_852],
            "count": "39231",
            "target_indicator": 1,
            "pre_A_coefficients": vector(0, 1, 1, 1, 1, 1),
            "A_transition_coefficients": zero,
            "post_A_coefficients": vector(0, 1, 1, 1, 1, 1),
            "post_A_owner": "joined_R_after_A",
            "separate_norm_admissible": False,
        },
        {
            "sector_id": "A_window_target_pairs",
            "mode_range": [39_853, 39_894],
            "count": "42",
            "target_indicator": 1,
            "pre_A_coefficients": vector(0, 1, 1, 1, 1, 1),
            "A_transition_coefficients": vector(0, 0, 1, 1, 0, 0),
            "post_A_coefficients": vector(0, 1, 0, 0, 1, 1),
            "post_A_owner": "joined_R_after_A_except_extracted_A_endpoint",
            "separate_norm_admissible": False,
        },
        {
            "sector_id": "A_window_outer_pairs",
            "mode_range": [39_895, 39_936],
            "count": "42",
            "target_indicator": 0,
            "pre_A_coefficients": vector(1, 1, 1, 1, 1, 1),
            "A_transition_coefficients": vector(1, 0, 1, 1, 0, 0),
            "post_A_coefficients": vector(0, 1, 0, 0, 1, 1),
            "post_A_owner": "joined_R_after_A_except_extracted_A_endpoint_and_positive_bulk",
            "separate_norm_admissible": False,
        },
        {
            "sector_id": "remote_positive_pairs",
            "mode_range": [39_937, "M"],
            "count": "M-39936",
            "target_indicator": 0,
            "pre_A_coefficients": vector(1, 1, 1, 1, 1, 1),
            "A_transition_coefficients": zero,
            "post_A_coefficients": vector(1, 1, 1, 1, 1, 1),
            "post_A_owner": "joined_R_after_A",
            "separate_norm_admissible": False,
        },
    ]
    for row in rows[1:]:
        require(
            add_vectors(row["A_transition_coefficients"], row["post_A_coefficients"])
            == row["pre_A_coefficients"],
            f"coefficient split failed: {row['sector_id']}",
        )
    return rows


def symbolic_certificate() -> dict[str, str]:
    Pp, Pm, Ap, Am, Bp, Bm, chi = sp.symbols("Pp Pm Ap Am Bp Bm chi")
    Ip = Pp + Ap + Bp
    Im = Pm + Am + Bm
    paired = sp.expand(Ip + Im - chi * Pp)
    projector = sp.expand((1 - chi) * Pp + Pm + Ap + Am + Bp + Bm)
    require(sp.expand(paired - projector) == 0, "paired projector identity failed")
    EB, EO, AD, RJ, RKG, RDir = sp.symbols("EB EO AD RJ RKG RDir")
    require(sp.expand((EB + EO + RDir).subs(RDir, AD + RJ) - (EB + EO + AD + RJ)) == 0, "global join identity failed")
    return {
        "paired_projector": "I_m+I_-m-chi_T(m)P_m=(1-chi_T(m))P_m+P_-m+A_m+A_-m+B_m+B_-m",
        "finite_defect": "Delta_(M,epsilon)=H_x+I_0+sum_(m=1)^M w_m[I_m+I_-m-chi_T(m)P_m]",
        "A_extraction": "A_(M,epsilon)=sum_(m=39853)^39936 w_m[A_m+A_-m+(1-chi_T(m))P_m]",
        "post_A_finite_join": "J_A,(M,epsilon)=Delta_(M,epsilon)-A_(M,epsilon)-Btr_win,(M,epsilon)-B_outer,(M,epsilon)",
        "physical_join": "R_KGamma=E_Btr,win+E_outer+A_transition+R_after_A",
        "R_Dir_split": "R_Dir=A_transition+R_after_A",
        "atom_dictionary": "P_+ and P_- are full-line bulk atoms; A_+/- and B_+/- are the two endpoint atoms in Q_+/-",
    }


def partition_certificate() -> dict[str, Any]:
    M = sp.symbols("M", integer=True, positive=True)
    positive_count = 621 + 39_231 + 42 + 42 + (M - 39_936)
    require(sp.expand(positive_count - M) == 0, "positive-mode partition count failed")
    samples = []
    for cutoff in (39_936, 5_122_421, 5_122_438):
        counts = [621, 39_231, 42, 42, cutoff - 39_936]
        require(sum(counts) == cutoff, f"sample partition failed at M={cutoff}")
        samples.append({"cutoff": cutoff, "class_counts": counts, "positive_total": sum(counts), "full_symmetric_total": 2 * cutoff + 1})
    return {
        "cutoff_condition": "integer M>=39936",
        "positive_count_identity": "621+39231+42+42+(M-39936)=M",
        "full_symmetric_count_identity": "1+2M=2M+1",
        "samples": samples,
    }


def extraction_certificate(dependencies: dict[str, dict[str, Any]]) -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    A_ball = arb(
        dependencies["A_transition"]["certificate"]["components"]["projector_completed_exact_A_transition_ball"]
    )
    working_target = arb(WORKING_TARGET)
    negative_target = arb(NEGATIVE_TARGET)
    working_room = working_target - A_ball
    negative_room = negative_target - A_ball
    # Publish rounded-down scalars that remain conservative under a
    # higher-precision reconstruction of the dependency interval.
    working_safe_text = "0.0368147039947"
    negative_safe_text = "0.0368061039947"
    working_safe = arb(working_safe_text)
    negative_safe = arb(negative_safe_text)
    require(working_safe < arb(working_room.lower()), "published working threshold is not conservative")
    require(negative_safe < arb(negative_room.lower()), "published negative threshold is not conservative")
    require(working_safe > arb("0.0368"), "post-A working room scale guard lost")
    require(negative_safe > arb("0.03679"), "post-A negative room scale guard lost")
    B_window = dependencies["B_window_bound"]["interval_certificate"]
    B_outer = dependencies["B_outer_bound"]["interval_certificate"]
    Gamma = dependencies["Gamma_insertion"]["correction_certificate"]
    return {
        "extractions": [
            {
                "sector_id": "B_trace_window",
                "owner": "outside_R_Dir",
                "support_or_grouping": "chi_W Btr; all nonlocal trace currents plus exact local modes 621 and 622",
                "certified_value_or_bound": B_window["complete_B_trace_package_window_physical_upper_ball"],
                "bound_type": "signed_upper",
                "overlap_guard": "A extraction uses A endpoint atoms and positive bulk only; Btr uses B trace atoms and local bulk only at modes 621 and 622",
            },
            {
                "sector_id": "analytic_B_outer",
                "owner": "outside_R_Dir",
                "support_or_grouping": "completed B current on the released exterior for modes m>=5122421",
                "certified_value_or_bound": B_outer["physical_outer_bound_ball"],
                "bound_type": "absolute_upper",
                "overlap_guard": "mode range starts above every A-transition mode",
            },
            {
                "sector_id": "projector_completed_A_transition",
                "owner": "inside_R_Dir_extracted_signed",
                "support_or_grouping": "A endpoint atoms on 39853..39936 plus positive Gamma bulk on 39895..39936",
                "certified_value_or_bound": A_ball.str(80, more=True),
                "bound_type": "signed_interval",
                "overlap_guard": "target Gamma modes 622..39894 remain subtracted; only the outside-target positive bulk coefficient is restored",
            },
            {
                "sector_id": "joined_R_after_A",
                "owner": "inside_R_Dir_open",
                "support_or_grouping": "half-current, zero mode, all unextracted negative/endpoint/remote-positive atoms, and the finite Dirichlet/Abel completion after B extractions",
                "certified_value_or_bound": None,
                "bound_type": "open_one_sided_join",
                "overlap_guard": "no independent norm is permitted on zero, half-current, negative, or remote-positive subsectors",
            },
        ],
        "A_transition_ball": A_ball.str(80, more=True),
        "working_R_Dir_target": WORKING_TARGET,
        "negative_R_Dir_target": NEGATIVE_TARGET,
        "working_R_after_A_threshold_ball": working_room.str(80, more=True),
        "negative_R_after_A_threshold_ball": negative_room.str(80, more=True),
        "working_R_after_A_safe_upper": working_safe_text,
        "negative_R_after_A_safe_upper": negative_safe_text,
        "working_threshold_expansion_factor_ball": (working_safe / working_target).str(50, more=True),
        "Gamma_to_classical_absolute_bound_ball": Gamma["Gamma_to_classical_physical_absolute_bound_ball"],
        "target_implications": [
            "R_after_A < inf(working_target-A_transition) implies R_Dir < 0.00014057919999999995",
            "R_after_A < inf(negative_target-A_transition) implies R_Dir < 0.00013197919999999995",
        ],
    }


def render_note(artifact: dict[str, Any]) -> str:
    cert = artifact["extraction_certificate"]
    return f"""# Global R_Dir projector ownership ledger

Date: 2026-08-23

Status: exact saved-height ownership and post-A target certified; joined
remainder still open

For a common finite cutoff `M>=39936` and the inherited Abel weight, the
positive/negative pair identity is

```text
I_m+I_-m-chi_T(m)P_m
 =(1-chi_T(m))P_m+P_-m+A_m+A_-m+B_m+B_-m.          (OL1)
```

The six positive-label classes are `1..621`, `622..39852`,
`39853..39894`, `39895..39936`, and `39937..M`, together with the zero
mode.  Their counts satisfy

```text
621+39231+42+42+(M-39936)=M.                        (OL2)
```

The certified A extraction removes `A_m+A_-m` on all 84 A-window modes and
also removes `P_m` on the 42 outside-target modes `39895..39936`.  It removes
no `P_-m` or B endpoint atom.  Thus

```text
R_KGamma=E_Btr,win+E_outer+A_transition+R_after_A,
R_Dir=A_transition+R_after_A,                        (OL3)

A_transition={cert['A_transition_ball']}.            (OL4)
```

The two certified one-sided continuation targets are

```text
R_after_A < {cert['working_R_after_A_safe_upper']}
  ==> R_Dir < {cert['working_R_Dir_target']},

R_after_A < {cert['negative_R_after_A_safe_upper']}
  ==> R_Dir < {cert['negative_R_Dir_target']}.        (OL5)
```

The first threshold is about
`{cert['working_threshold_expansion_factor_ball']}` times the original
`R_Dir` working scale.  This is only target arithmetic.  It does not say the
remaining joined current is small.

The half-current, zero mode, negative full-line bulk, remaining endpoint
atoms, and remote-positive modes stay inside one `R_after_A`.  The exact
double-Weber reconstruction shows that the zero, direct, modular, and
endpoint-half sectors are one copy of the finite source, while symmetric
pairing cancels the apparent one-over-`m` current.  Taking separate norms on
those labels would destroy the proved cancellation.

The B trace window is already outside `R_Dir`.  The completed analytic B
outer block is also outside it.  The finite completed B block is not a third
numerical contribution: it is the common Dirichlet/Abel representation from
which the A positive-bulk modes `39895..39936` are now reassigned exactly.

Pi provenance: every `pi` in (OL1)--(OL5) is inherited from the equation-(9)
Fourier kernel, the Gaussian Abel regulator, and the previously certified
Gamma/A endpoint identities.  No geometric or fitted occurrence is added.

Proof boundary: exact mode/component ownership, no-overlap guards, and
post-A sufficient target arithmetic at `t=10^10` only.  No bound for
`R_after_A`, complete `R_Dir`, `Q_K-T`, all-height theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["symmetric_reassembly"]["decision"]["target_negative_partners_retained_in_complement"] is True, "negative-partner ownership drift")
    require(dependencies["double_weber_reconstruction"]["decision"]["entire_finite_equation9_source_roster_reconstructed_exactly"] is True, "double-Weber reconstruction drift")
    require(dependencies["B_window_extraction"]["decision"]["remaining_object_is_one_grouped_global_remainder"] is True, "B extraction grouping drift")
    require(dependencies["B_window_bound"]["decision"]["complete_saved_height_B_trace_package_window_upper_bound_proved"] is True, "B window bound drift")
    require(dependencies["B_outer_bound"]["decision"]["complete_analytic_m_ge_B_exterior_physical_contribution_below_8e_minus_10"] is True, "B outer bound drift")
    require(dependencies["finite_dirichlet_rejoin"]["decision"]["finite_completed_B_block_rejoined_before_norms"] is True, "finite rejoin drift")
    require(dependencies["Gamma_insertion"]["decision"]["Dirichlet_target_kernel_is_exact_Gamma_bulk_insertion"] is True, "Gamma insertion drift")
    require(dependencies["projector_defect"]["decision"]["target_full_line_carrier_subtracted_before_norms"] is True, "projector drift")
    require(dependencies["A_transition"]["decision"]["projector_completed_exact_A_transition_certified"] is True, "A transition drift")

    artifact = {
        "kind": STEM,
        "status": "exact_global_R_Dir_projector_ownership_and_post_A_join_target_certified_quantitative_join_open",
        "passed": True,
        "scope": {
            "height": HEIGHT,
            "finite_cutoff": "M>=39936",
            "target_positive_modes": [TARGET_START, TARGET_END],
            "A_transition_modes": [A_START, A_END],
            "analytic_B_outer_start": B_OUTER_START,
            "common_regulator": "w_m=exp(-pi*epsilon*m^2), symmetric M limit before epsilon down to zero",
        },
        "symbolic_certificate": symbolic_certificate(),
        "partition_certificate": partition_certificate(),
        "mode_ownership_rows": mode_ownership_rows(),
        "extraction_certificate": extraction_certificate(dependencies),
        "decision": {
            "all_finite_mode_classes_partitioned_without_gap_or_overlap": True,
            "A_endpoint_and_positive_bulk_coefficients_transferred_exactly": True,
            "B_window_A_transition_and_analytic_B_outer_no_overlap_guards_certified": True,
            "zero_half_negative_and_remote_positive_sectors_remain_joined": True,
            "post_A_one_sided_targets_derived": True,
            "joined_R_after_A_bound_proved": False,
            "complete_R_Dir_proved": False,
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
        "next_obligation": "Construct a cancellation-preserving evaluator or analytic enclosure for R_after_A directly from the common Dirichlet/Abel representation with the 84 A endpoint atoms and 42 outside-target positive bulk atoms deleted algebraically. Keep the half-current, zero, negative, and remote-positive sectors joined. Target the certified one-sided upper bound recorded here before attempting any all-height promotion.",
        "proof_boundary": "Exact saved-height mode/component ownership, no-overlap guards, and post-A sufficient target arithmetic only. No quantitative R_after_A or R_Dir bound, complete Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified global R_Dir ownership ledger and post-A targets", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
