#!/usr/bin/env python3
"""Run the core proof-programme reproducibility gates.

This is a compact smoke/ledger runner for the RH dynamical-rigidity corpus.
It does not prove PF-infinity, Laguerre-Polya membership, RH, or Lambda <= 0.
It checks that the promoted finite-evidence manifests and executable
countermodel guards still match the advertised proof-programme status.
"""

from __future__ import annotations

import argparse
import ctypes
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from time import perf_counter


REPO_ROOT = Path(__file__).resolve().parents[3]


def set_below_normal_priority() -> None:
    """Keep a long umbrella run responsive to other desktop work on Windows."""
    if os.name != "nt":
        return
    try:
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        kernel32.SetPriorityClass.argtypes = (ctypes.c_void_p, ctypes.c_uint32)
        kernel32.SetPriorityClass.restype = ctypes.c_int
        kernel32.SetPriorityClass(kernel32.GetCurrentProcess(), 0x00004000)
    except (AttributeError, OSError):
        pass


@dataclass(frozen=True)
class GateSpec:
    name: str
    command: tuple[str, ...]
    expected: tuple[str, ...]
    category: str
    slow: bool = False


GATES: tuple[GateSpec, ...] = (
    GateSpec(
        name="countermodel proof-safety gates",
        command=("work/rh_compute/scripts/countermodel_gate_examples.py",),
        expected=("validated 11 countermodel gate examples",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="result-language boundary scan",
        command=("work/rh_compute/scripts/check_result_language_boundaries.py",),
        expected=("validated result-language boundaries: scanned", "0 overclaims"),
        category="non-promotion guards",
    ),
    GateSpec(
        name="output artifact status manifest",
        command=("work/rh_compute/scripts/check_output_status_manifest.py",),
        expected=("validated output artifact statuses: scanned", "0 status issues"),
        category="non-promotion guards",
    ),
    GateSpec(
        name="output reference integrity",
        command=("work/rh_compute/scripts/check_output_reference_integrity.py",),
        expected=("validated output references: scanned", "0 missing required paths"),
        category="non-promotion guards",
    ),
    GateSpec(
        name="proof-claim ledger",
        command=("work/rh_compute/scripts/check_proof_claim_ledger.py",),
        expected=("validated proof-claim ledger:", "0 issues, 10 open theorem targets"),
        category="non-promotion guards",
    ),
    GateSpec(
        name="signed-Hankel/Jensen dependency graph",
        command=("work/rh_compute/scripts/check_signed_hankel_jensen_dependency_graph.py",),
        expected=("validated signed-Hankel/Jensen dependency graph with 0 issues",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="signed-Hankel finite certificates",
        command=("work/rh_compute/scripts/check_hankel_certificate_manifest.py",),
        expected=("validated 2500 signed-Hankel finite certificates",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Toeplitz/PF finite certificates",
        command=("work/rh_compute/scripts/check_toeplitz_certificate_manifest.py",),
        expected=("validated 95 promoted positive certificate summaries and 1 negative control",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Toeplitz/Jacobi-Trudi reindexing",
        command=("work/rh_compute/scripts/check_toeplitz_jacobi_trudi_map.py",),
        expected=("validated Toeplitz/Jacobi-Trudi reindexing: N=10, orders<=5, 124129 minors, 39094 nonzero maps, 85035 structural zeros",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Hankel sign-consistency reduction point audits",
        command=("work/rh_compute/scripts/check_hankel_sign_consistency_reduction_audit.py",),
        expected=("validated 20 reshaped Hankel sign-consistency point audits with 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Hankel sign-consistency reduction finite certificates",
        command=("work/rh_compute/scripts/check_arb_hankel_sign_consistency_reduction_manifest.py",),
        expected=("validated 689795 Arb reshaped-Hankel sign-consistency finite certificates with 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="shifted Hankel staircase finite certificates",
        command=("work/rh_compute/scripts/check_arb_shifted_hankel_staircase_manifest.py",),
        expected=("validated 3154515 shifted Arb reshaped-Hankel staircase finite certificates with 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen/Hankel bridge algebra gate",
        command=("work/rh_compute/scripts/check_jensen_hankel_bridge_algebra.py",),
        expected=("validated Jensen/Hankel bridge algebra gate: degree2 identity and degree3 finite countermodel with 0 issues",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF obligation algebra gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_obligation_algebra.py",),
        expected=("validated Jensen-window PF obligation algebra with 0 issues",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Arb Jensen-window PF obligation finite diagnostics",
        command=("work/rh_compute/scripts/check_arb_jensen_window_pf_obligation_manifest.py",),
        expected=("validated 1470 Arb Jensen-window PF obligation finite diagnostics with 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Arb Jensen-window Sturm hyperbolicity finite diagnostics",
        command=("work/rh_compute/scripts/check_arb_jensen_window_sturm_manifest.py",),
        expected=("validated 210 Arb Jensen-window Sturm hyperbolicity finite diagnostics with 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Arb Jensen-window Sturm degree-5 hyperbolicity finite diagnostics",
        command=(
            "work/rh_compute/scripts/check_arb_jensen_window_sturm_manifest.py",
            "--summary",
            "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520_summary.json",
            "--rows",
            "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520.jsonl",
            "--expected-degrees",
            "5",
        ),
        expected=("validated 105 Arb Jensen-window Sturm hyperbolicity finite diagnostics with 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Arb Jensen-window Sturm degree-6 through degree-12 hyperbolicity finite diagnostics",
        command=(
            "work/rh_compute/scripts/check_arb_jensen_window_sturm_manifest.py",
            "--summary",
            "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d6_d12_dps520_summary.json",
            "--rows",
            "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d6_d12_dps520.jsonl",
            "--expected-degrees",
            "6..12",
        ),
        expected=("validated 735 Arb Jensen-window Sturm hyperbolicity finite diagnostics with 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="finite Sturm-to-PF Jensen-window consequences",
        command=("work/rh_compute/scripts/check_jensen_window_sturm_pf_consequence.py",),
        expected=("validated 1050 finite Sturm-to-PF Jensen-window consequences with 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="signed-Hankel/Jensen bridge target specification",
        command=("work/rh_compute/scripts/check_signed_hankel_jensen_bridge_target.py",),
        expected=("validated signed-Hankel/Jensen bridge target specification with 0 issues",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF bridge target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_bridge_target.py",),
        expected=("validated Jensen-window PF bridge target with 0 issues",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF coefficient-PF equivalence gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_coefficient_pf_equivalence_gate.py",),
        expected=(
            "validated Jensen-window PF coefficient-PF equivalence gate: 13 rows, 0 issues, 3 exact coefficient identities, 4 classical/closure steps, 1 seven-way equivalence, 3 guards, 1 open structural handoff",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Edrei-Stieltjes equivalence gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_edrei_stieltjes_equivalence_gate.py",),
        expected=(
            "validated Jensen-window PF Edrei-Stieltjes equivalence gate: 14 rows, 0 issues, 12 exact indexing checks, 7 exact Hankel checks, 1 unified endpoint, 3 finite/nonpromotion guards, 1 open Xi/Phi handoff",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Edrei heat-flow boundary gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_edrei_heat_flow_boundary_gate.py",),
        expected=(
            "validated Jensen-window PF Edrei heat-flow boundary gate: 11 rows, 0 issues, 4 exact flow identities, 9 rank-one orientation checks, 2 exact heat witnesses, 1 rejected generic backward-invariance shortcut, 1 open Xi/Phi handoff",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Edrei Hankel boundary-flux gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_edrei_hankel_boundary_flux_gate.py",),
        expected=(
            "validated Edrei Hankel boundary-flux gate: 18 rows, 0 issues, 25 null-form audits, 12 determinant audits, 15 Cauchy-Binet audits, 24 orthogonal-quotient audits, 43 countermodel/collision audits, 1 open Xi/Phi determinant-growth gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Edrei raw-moment collision-resolution gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_edrei_raw_moment_collision_resolution_gate.py",),
        expected=(
            "validated Edrei raw-moment collision-resolution gate: 14 rows, 0 issues, 104 triangular transfer audits, 36 split determinant audits, 75 simple-root/flux audits, 36 repeated-flux audits, 5 conditioning audits, 24 general-r collision determinant audits, 24 general-r Cauchy-Binet audits, 24 general-r crossover-ratio audits, 4 crossover audits, 1 noncommuting-limit obstruction, 1 open Xi/Phi precision handoff",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Phi Pick-kernel target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_phi_pick_kernel_target.py",),
        expected=(
            "validated Jensen-window PF Phi Pick-kernel target: 12 rows, 0 issues, 7 exact identities, 3 independent polarization checks, 1 exact mixture guard, 2 open structural routes",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Xi Pick/Suzuki Hankel bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_xi_pick_suzuki_hankel_bridge.py",),
        expected=(
            "validated Jensen-window PF Xi Pick/Suzuki Hankel bridge: 16 rows, 0 issues, 7 exact coordinate/guard identities, 4 published operator steps, 3 exact countermodels, 1 open global gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Suzuki spectral frontier",
        command=("work/rh_compute/scripts/check_jensen_window_pf_suzuki_spectral_frontier.py",),
        expected=(
            "validated Jensen-window PF Suzuki spectral frontier: 16 rows, 0 issues, 10 exact path/reduction identities, 5 route guards, 1 open global obligation",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Suzuki determinant-only reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_suzuki_determinant_only_reduction.py",),
        expected=(
            "validated Suzuki determinant-only reduction: 16 rows, 0 issues, 12 exact bridge steps, 2 theorem/corollary candidates, 1 route guard, 1 open arithmetic gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Suzuki fixed-omega phase diagram",
        command=("work/rh_compute/scripts/check_jensen_window_pf_suzuki_fixed_omega_phase_diagram.py",),
        expected=(
            "validated Suzuki fixed-omega phase diagram: 21 rows, 0 issues, 3 exact countermodel guards, 2 published scalar targets, 1 open arithmetic gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Suzuki Jordan-totient sign scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_totient_sign_scout.py",),
        expected=(
            "validated Suzuki Jordan-totient sign scout: 6000 samples, 5 omega values, 0 negative rows, 0 issues",
        ),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF Suzuki cofinal monotonicity hierarchy",
        command=("work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.py",),
        expected=(
            "validated Suzuki cofinal monotonicity hierarchy: 21 rows, 0 issues, 3 signed/nonpromotion guards, 1 cofinal equivalence candidate, 1 open arithmetic gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Suzuki Jordan-error kernel reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_error_kernel_reduction.py",),
        expected=(
            "validated Suzuki Jordan-error kernel reduction: 16 rows, 0 issues, 8 exact identities/bounds, 2 route guards, 1 open cancellation gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Suzuki cofinal L2 hierarchy",
        command=("work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_l2_hierarchy.py",),
        expected=(
            "validated Suzuki cofinal L2 hierarchy: 22 rows, 0 issues, 3 countermodel guards, 1 literature-fit guard, 1 cofinal equivalence candidate, 1 open arithmetic energy gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Jordan-Muntz causal-energy bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_causal_energy_bridge.py",),
        expected=(
            "validated Jordan-Muntz causal-energy bridge: 32 rows, 0 issues, 5 theorem candidates, 5 proof guards, 1 open all-height mollifier gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Jordan-Muntz/Burnol Hardy intertwiner",
        command=("work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.py",),
        expected=(
            "validated Jordan-Muntz/Burnol Hardy intertwiner: 18 rows, 0 issues, 11 exact identities, 1 source-backed cofinal theorem candidate, 2 nonpromotion guards, 1 open natural-mollifier gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Burnol cell-energy/tail obstruction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_burnol_cell_energy_tail_obstruction.py",),
        expected=(
            "validated Burnol cell-energy/tail obstruction: 18 rows, 0 issues, 10 exact identities, 4 norm-reduction steps, 1 reciprocal-zeta tail obstruction, 1 open short-interval gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Burnol tail-discrepancy/dyadic reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.py",),
        expected=(
            "validated Burnol tail-discrepancy/dyadic reduction: 20 rows, 0 issues, 12 exact identities, 4 equivalence steps, 2 literature guards, 1 open dyadic square-function gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF weighted fractional autocorrelation Gram bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.py",),
        expected=(
            "validated weighted fractional autocorrelation Gram bridge: 20 rows, 0 issues, 12 exact identities, 3 proof guards, 1 open signed off-diagonal gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF weighted autocorrelation OU tail-energy reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.py",),
        expected=(
            "validated weighted-autocorrelation OU tail-energy reduction: 22 rows, 0 issues, 13 exact identities, 4 proof guards, 1 open comparison gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF OU/Mertens mean-square reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_ou_mertens_mean_square_reduction.py",),
        expected=(
            "validated OU/Mertens mean-square reduction: 25 rows, 0 issues, 16 exact reductions, 4 proof guards, 1 open anchored mean-square gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens weighted-prefix/affine-defect reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.py",),
        expected=(
            "validated Mertens weighted-prefix/affine-defect reduction: 34 rows, 0 issues, 24 exact reductions, 4 proof guards, 2 open handoff gates",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens anchor/logarithmic tail-energy reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.py",),
        expected=(
            "validated Mertens anchor/logarithmic tail-energy reduction: 26 rows, 0 issues, 17 exact reductions, 4 proof guards, 1 open logarithmic gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens dyadic cosine-mode bottleneck",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.py",),
        expected=(
            "validated Mertens dyadic cosine-mode bottleneck: 26 rows, 0 issues, 15 exact reductions, 4 proof guards, 1 open low-mode gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens local-path cosine transference",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_local_path_cosine_transference.py",),
        expected=(
            "validated Mertens local-path cosine transference: 32 rows, 0 issues, 24 exact reductions, 5 proof guards, 1 open gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens centered-bridge Vaughan handoff",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_centered_bridge_vaughan_handoff.py",),
        expected=(
            "validated Mertens centered-bridge Vaughan handoff: 40 rows, 0 issues, 34 exact reductions, 3 proof guards, 1 open gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens affine-tent/bridge handoff",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_affine_tent_bridge_handoff.py",),
        expected=(
            "validated Mertens affine-tent/bridge handoff: 53 rows, 0 issues, 44 exact reductions, 4 proof guards, 1 open gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens ordinary-vector Vaughan handoff",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.py",),
        expected=(
            "validated Mertens ordinary-vector Vaughan handoff: 34 rows, 0 issues, 28 exact reductions, 3 proof guards, 1 open gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens mixed-boundary sine/Vaughan handoff",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.py",),
        expected=(
            "validated Mertens mixed-boundary sine/Vaughan handoff: 35 rows, 0 issues, 27 exact reductions, 3 proof guards, 1 open gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens spectral anchor/high-mode reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.py",),
        expected=(
            "validated Mertens spectral anchor/high-mode reduction: 24 rows, 0 issues, 16 exact reductions, 2 proof guards, 1 open gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens truncated half-odd kernel handoff",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.py",),
        expected=(
            "validated Mertens truncated half-odd kernel handoff: 45 rows, 0 issues, 39 exact reductions, 2 proof guards, 1 open gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens shift-kernel variation handoff",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_shift_kernel_variation_handoff.py",),
        expected=(
            "validated Mertens shift-kernel variation handoff: 45 rows, 0 issues, 39 exact reductions, 2 proof guards, 1 open gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens planar Abel handoff",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_abel_handoff.py",),
        expected=(
            "validated Mertens planar Abel handoff: 36 rows, 0 issues, 28 exact reductions, 3 proof guards, 1 open gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens planar edge-Gram/Vaughan handoff",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.py",),
        expected=(
            "validated planar edge-Gram/Vaughan handoff: 38 rows, 0 issues, 948 Vaughan audits, 4720 band audits, 20 Gram audits, 87 Fourier audits, 1 open signed edge gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens planar incidence/anti-diagonal handoff",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.py",),
        expected=(
            "validated planar incidence/anti-diagonal handoff: 43 rows, 0 issues, 3561 incidence audits, 13137 cut audits, 334 Fourier audits, 528 anti-diagonal audits, 1 open projection gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens planar boundary-flux/anchor handoff",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.py",),
        expected=(
            "validated planar boundary-flux/anchor handoff: 35 rows, 0 issues, 2122 flux audits, 897 endpoint audits, 216 block audits, 6 witness audits, 70 Mobius audits, 1 open joined energy gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens planar boundary suffix-energy handoff",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.py",),
        expected=(
            "validated planar boundary suffix-energy handoff: 26 rows, 0 issues, 132 algebra audits, 492 edge audits, 906 bound audits, 42 Mobius audits, 1 open anchored energy gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens planar boundary suffix-localization gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.py",),
        expected=(
            "validated planar boundary suffix-localization gate: 26 rows, 0 issues, 41192 odd-sine audits, 10339 sliver audits, 9259 bulk audits, 306 collar audits, 49 diagnostic/countermodel audits, 1 open RH-equivalent energy gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens planar joined-energy equivalence gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.py",),
        expected=(
            "validated planar joined-energy equivalence gate: 20 rows, 0 issues, 189 arbitrary-vector audits, 63 Mobius audits, 1 open lossless energy gate",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Mertens planar curvature-energy scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_curvature_energy_scout.py",),
        expected=(
            "validated Mertens planar curvature-energy scout: 28 rows, 0 issues, K<=1024, 4 alpha values",
        ),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF fixed-shift inner/reciprocal-boundary separation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.py",),
        expected=(
            "validated fixed-shift inner/reciprocal-boundary separation gate: 18 rows, 0 issues, 15 exact reductions, 2 proof guards, 1 open target",
        ),
        category="countermodel proof safety",
    ),
    GateSpec(
        name="Jensen-window PF bridge obligation ledger",
        command=("work/rh_compute/scripts/check_jensen_window_pf_bridge_obligations.py",),
        expected=("validated Jensen-window PF bridge obligations: 16 obligations, 0 issues, 3 open obligations",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF theorem machinery fit matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py",),
        expected=("validated Jensen-window PF theorem machinery fit matrix: 11 rows, 0 issues, 0 ready-to-apply rows",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF sign-regular transfer gap matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_sign_regular_transfer_gap_matrix.py",),
        expected=(
            "validated Jensen-window PF sign-regular transfer gap matrix: 9 transfer rows, 2 countermodel gates, 3 open requirements, 3 rejected shortcuts, 0 ready-to-apply rows, 0 issues",
        ),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF factorial multiplier split audit",
        command=("work/rh_compute/scripts/check_jensen_window_pf_factorial_multiplier_split_audit.py",),
        expected=(
            "validated Jensen-window PF factorial multiplier split audit: 5 exact rows, 315 raw degree-2 anti-hyperbolic rows, 315 normalized degree-2 positive rows, 0 ready-to-apply rows, 0 issues",
        ),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF reciprocal-gamma mixture sign gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_reciprocal_gamma_mixture_sign_gate.py",),
        expected=(
            "validated Jensen-window PF reciprocal-gamma mixture sign gate: 10 rows, 0 issues, 3 exact kernel identities, 1 published all-order theorem, 1 fixed-scale theorem, 2 exact mixture countermodels, 1 tilted-variance equivalence, 1 Xi order-two composition, 1 higher-order handoff",
        ),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF reciprocal-defect compound order-three gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_reciprocal_defect_compound_order3_gate.py",),
        expected=(
            "validated Jensen-window PF reciprocal-defect compound order-three gate: 11 rows, 0 issues, 2 exact coordinate identities, 2 exact sign equivalences, 1 sufficient increment theorem, 1 exact boundary benchmark, 1 strict cone countermodel, 1 all-shift lambda=-100 entry theorem, 1 full forward propagation theorem, 1 arbitrary-column order-three theorem, 1 open order-four handoff",
        ),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda -100 compound order-three entry certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.py",),
        expected=(
            "validated Jensen-window PF negative-lambda -100 compound order-three entry certificate: 11 rows, 0 issues, 322 positive coefficients, 318 prefix compound margins, 318 prefix defect gaps, 1 exact tail theorem, 1 all-shift entry theorem, 1 open forward handoff",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-three forward-invariance certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order3_forward_invariance_certificate.py",),
        expected=(
            "validated Jensen-window PF compound order-three forward-invariance certificate: 10 rows, 0 issues, 2 exact identities, 1 cooperative flow, 1 inward boundary theorem, 1 coefficient-growth lemma, 1 infinite maximum principle, 1 full forward propagation theorem, 1 lambda=0 theorem, 1 open higher-compound handoff",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF order-three noncontiguous secant-transfer lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.py",),
        expected=(
            "validated Jensen-window PF order-three noncontiguous secant-transfer lemma: 8 rows, 0 issues, 2 exact identities, 1 secant-averaging lemma, 1 arbitrary-column order-two theorem, 1 arbitrary-column order-three theorem, 1 open order-four handoff",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four condensation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_condensation_gate.py",),
        expected=(
            "validated Jensen-window PF compound order-four condensation gate: 8 rows, 0 issues, 2 exact identities, 1 exact sign equivalence, 317 positive lambda=-100 prefix margins, 1 strict lower-order countermodel, 1 forbidden promotion, 2 open handoffs",
        ),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four first-summand curvature bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_first_summand_curvature_bridge.py",),
        expected=(
            "validated Jensen-window PF compound order-four first-summand curvature bridge: 11 rows, 0 issues, 2 exact identities, 2 exact reductions, 1 proved gap floor, 1 compact interval theorem, 1 proved full-kernel perturbation theorem, 1 open analytic ray, 8 positive finite scouts",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four localized-curvature compact certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_localized_curvature_compact_certificate.py",),
        expected=(
            "validated localized order-four compact curvature certificate: 107452 H tiles, 1073 positive blocks, 1 open analytic ray, 0 issues",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four Gaussian cumulant ray target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.py",),
        expected=(
            "validated order-four Gaussian cumulant ray target: 10 rows, 0 issues, 4 exact formal rows, 2 proved formal-corridor rows, 7 positive conditional collars, 3 open analytic rows",
        ),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four formal cumulant corridor certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.py",),
        expected=(
            "validated order-four formal cumulant corridor certificate: 1800000 blocks, 1800 chunks, 7 positive corridors, 0 issues",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four formal cumulant asymptotic ray certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.py",),
        expected=(
            "validated order-four formal cumulant asymptotic ray: 8 rows, 0 issues, 7 exact rows, 7 buffered corridors, 14 jet-remainder sign gates, 1 open exact-density remainder",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four formal cumulant next-parity certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.py",),
        expected=(
            "validated order-four formal cumulant next parity: 6 rows, 0 issues, 4 exact rows, 42 epsilon-six audits, 7 next-parity coefficients, 2 open analytic rows",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four next-parity finite certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_formal_cumulant_next_parity_finite_certificate.py",),
        expected=(
            "validated order-four next-parity finite bounds: 1800 Taylor blocks, 180 chunks, 7 signed coefficients, 0 issues",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four next-parity asymptotic ray certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate.py",),
        expected=(
            "validated order-four next-parity asymptotic ray: 8 rows, 0 issues, 7 exact rows, 7 global coefficient bounds, 14 leading buffer gates, 4 new jet-remainder sign gates, 1 open exact-density remainder",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four exact cumulant remainder budget",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.py",),
        expected=(
            "validated order-four exact cumulant remainder budget: 9 rows, 0 issues, 8 exact rows, finite epsilon-ten budget 9/1000, ray epsilon-ten budget 1/(100u), 1 open central-tail theorem",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four formal cumulant second-next parity certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate.py",),
        expected=(
            "validated order-four formal cumulant second-next parity: 5 rows, 0 issues, 3 exact rows, 56 epsilon-eight audits, 7 second-next coefficients, 2 open analytic rows",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four second-next finite certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate.py",),
        expected=(
            "validated order-four second-next finite bounds: 3600 Taylor blocks, 360 chunks, 7 signed coefficients, 0 issues",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four second-next asymptotic ray certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.py",),
        expected=(
            "validated order-four second-next asymptotic ray: 8 rows, 0 issues, 7 exact rows, 7 global coefficient bounds, 14 leading buffer gates, 4 new jet-remainder sign gates, 1 open exact-density remainder",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four exact cumulant complex-disk contract",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.py",),
        expected=(
            "validated order-four exact cumulant complex-disk contract: 10 rows, 0 issues, 9 exact rows, 70 cumulant audits, 2 sufficient partition targets, 1 open central-tail theorem",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four exact cumulant formal tails",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.py",),
        expected=(
            "validated order-four exact cumulant formal tails: 8 rows, 0 issues, 7 exact rows, 31 polynomial orders, 2 closed formal tails, 3 open exact components",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four exact cumulant exact tails",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.py",),
        expected=(
            "validated order-four exact cumulant exact tails: 9 rows, 0 issues, 8 exact rows, 2 positive-coefficient polynomials, 2 closed exact tails, 1 open central residual",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four finite partition extension",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate.py",),
        expected=(
            "validated order-four finite partition extension: 8 rows, 0 issues, 7 exact rows, 4 partition orders, 78 scalar functions, 5400 partition blocks, 5430 shifted-jet blocks, 5 new jet caps",
        ),
        category="interval certificates",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF compound order-four exact cumulant central residual",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate.py",),
        expected=(
            "validated order-four exact cumulant central residual: 11 rows, 0 issues, 11 exact rows, 222 Bell terms, 2 central regimes, 4 inherited tails, 0 open partition components",
        ),
        category="exact theorem composition",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF compound order-four exact cumulant corridor theorem",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.py",),
        expected=(
            "validated order-four exact cumulant corridors: 8 rows, 0 issues, 7 exact rows, 7 global exact corridors, 2 strict reserve regimes, 1 open localized-curvature composition",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four finite exact-corridor curvature theorem",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate.py",),
        expected=(
            "validated order-four exact-corridor finite curvature: 7 rows, 0 issues, 6 exact rows, 20700 mode blocks, 41400 t-collar gates, 20700 positive localized blocks, 1 open asymptotic ray",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four exact-corridor curvature ray theorem",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.py",),
        expected=(
            "validated order-four exact-corridor curvature ray: 8 rows, 0 issues, 8 exact rows, 7 normalized H boxes, 6 defect bounds, 5 localized gates, global corridor-to-curvature closed",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four lambda=-100 entry theorem",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_m100_entry_certificate.py",),
        expected=(
            "validated Jensen-window PF compound order-four lambda=-100 entry: 10 rows, 0 issues, 9 exact rows, 317 prefix margins, 1 analytic tail, 1 all-shift entry theorem, 1 open forward handoff",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four forward-flow reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_forward_flow_reduction.py",),
        expected=(
            "validated Jensen-window PF compound order-four forward flow: 8 rows, 0 issues, 7 exact rows, 3 exact identities, 2 cooperative flow lemmas, 1 maximum-principle reduction, 1 open spatial-tail bound",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Arb Xi lambda-zero order-four prefix certificate",
        command=("work/rh_compute/scripts/check_arb_xi_lambda0_order4_prefix_certificate.py",),
        expected=(
            "validated Arb Xi lambda-zero order-four prefix: 507 coefficients, 501 positive H4 rows, 501 positive stable margins, 0 inconclusive, 0 issues",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four lambda-zero eventual positivity",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.py",),
        expected=(
            "validated Jensen-window PF compound order-four lambda-zero eventual positivity: 8 rows, 0 issues, 7 exact rows, 7 symbolic coefficients, 501 finite prefix rows, 1 eventual theorem, 1 open effective splice",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF lambda-zero first-summand dominance transfer",
        command=("work/rh_compute/scripts/check_jensen_window_pf_lambda0_first_summand_dominance_transfer.py",),
        expected=(
            "validated lambda-zero first-summand dominance transfer: 7 rows, 0 issues, 6 exact rows, 1 lambda-zero dominance theorem, 1 order-four penalty transfer, 1 open curvature tail",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four eventual-tail forward reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.py",),
        expected=(
            "validated order-four eventual-tail forward invariance: 8 rows, 0 issues, 7 exact/input rows, 1 finite-confinement reduction, 1 conditional forward theorem, 1 open uniform tail",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four uniform-heat tail reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.py",),
        expected=(
            "validated order-four uniform-heat eventual tail: 7 rows, 0 issues, 6 exact/input rows, 7 symbolic coefficients, 1 conditional uniform transfer, 1 open heat-tilt theorem",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF uniform superpolynomial first-summand dominance",
        command=("work/rh_compute/scripts/check_jensen_window_pf_uniform_superpolynomial_first_summand_dominance.py",),
        expected=(
            "validated uniform superpolynomial first-summand dominance: 6 rows, 0 issues, 6 exact rows, 1 superpolynomial theorem, 1 local-difference theorem, 0 open rows",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF uniform first-summand heat-tilt asymptotic theorem",
        command=("work/rh_compute/scripts/check_jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.py",),
        expected=(
            "validated uniform first-summand heat tilt: 8 rows, 0 issues, 8 exact/published rows, 8 suitability coefficients, 7 Lambert derivative orders, 0 open rows",
        ),
        category="published theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-four uniform-heat forward invariance",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.py",),
        expected=(
            "validated order-four uniform-heat forward invariance: 9 rows, 0 issues, 9 ready rows, 1 uniform tail theorem, 1 forward theorem, 1 lambda-zero all-shift theorem, 0 open rows",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF order-four noncontiguous total-positivity transfer",
        command=("work/rh_compute/scripts/check_jensen_window_pf_order4_noncontiguous_total_positivity_transfer.py",),
        expected=(
            "validated order-four noncontiguous total-positivity transfer: 11 rows, 0 issues, 4 reversal orders, 240 solid-block maps, 1020 signed benchmark minors, 1 arbitrary-column order-four theorem, 1 fixed-order transfer theorem, 1 open order-five handoff",
        ),
        category="published theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-five uniform-tail and flow reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order5_uniform_tail_flow_reduction.py",),
        expected=(
            "validated order-five uniform tail and flow reduction: 12 rows, 0 issues, 12 suitability coefficients, 11 Lambert orders, 8 Newton coefficients, 120 determinant permutations, 1 uniform tail theorem, 1 cooperative flow theorem, 1 conditional forward theorem, 1 open lambda=-100 entry",
        ),
        category="exact theorem composition",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF compound order-five lambda=-100 prefix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order5_m100_prefix_certificate.py",),
        expected=(
            "validated order-five lambda=-100 prefix: 325 coefficients, 317 positive J5 rows, 317 positive relative margins, 317 positive H5 signs, 0 inconclusive, 1 open analytic tail",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-five lambda=-100 tail curvature reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order5_m100_tail_curvature_reduction.py",),
        expected=(
            "validated order-five lambda=-100 tail curvature reduction: 8 rows, 0 issues, 2 exact identities, 1 rational comparison, 1 conditional tail theorem, 1 open curvature target, cap 100/k^2 from k=321",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-five first-summand curvature bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order5_first_summand_curvature_bridge.py",),
        expected=(
            "validated order-five first-summand curvature bridge: 8 rows, 0 issues, 2 positive floors, 1 full-kernel transfer, 5 scout rows, 1 open continuous target, budgets 37+63=100",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-five nested-curvature compact certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order5_nested_curvature_compact_certificate.py",),
        expected=(
            "validated order-five nested curvature compact certificate: 5 rows, 0 issues, 36 interval blocks, 2 positive stable layers, 1 compact curvature theorem, 1 open ray, largest scaled upper 2.202668",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-five nested-curvature finite ray",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.py",),
        expected=(
            "validated order-five nested curvature finite-ray certificate: 5 rows, 0 issues, 100 extension tiles, 1850 exact-corridor blocks, 1 finite-ray theorem, 1 open asymptotic ray, largest scaled upper 11.9132",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-five nested-curvature asymptotic ray",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.py",),
        expected=(
            "validated order-five nested curvature asymptotic certificate: 6 rows, 0 issues, 1 normalized-H theorem, 1 stable-log majorant, 1 dimensionless interval theorem, 1 asymptotic curvature theorem, 0 open rows, scaled upper 9.15835",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-five lambda=-100 entry",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order5_m100_entry_certificate.py",),
        expected=(
            "validated order-five lambda=-100 entry certificate: 10 rows, 0 issues, 1 continuous curvature theorem, 1 complete scalar ceiling, 1 analytic tail theorem, 1 all-shift entry theorem, 0 open rows",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-five uniform-heat forward invariance",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.py",),
        expected=(
            "validated order-five uniform-heat forward invariance: 8 rows, 0 issues, 7 ready rows, 1 contiguous theorem, 1 arbitrary-column theorem, 1 open order-six handoff",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-six uniform-tail and flow reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order6_uniform_tail_flow_reduction.py",),
        expected=(
            "validated order-six uniform tail and flow reduction: 13 rows, 0 issues, 17 suitability coefficients, 16 Lambert orders, 10 Newton coefficients, 720 determinant permutations, 684 weighted monomials, 1 uniform signed tail theorem, 1 condensation recursion, 1 cooperative recursion, 1 conditional forward theorem, 1 open lambda=-100 entry",
        ),
        category="exact theorem composition",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF compound order-six lambda=-100 prefix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order6_m100_prefix_certificate.py",),
        expected=(
            "validated order-six lambda=-100 prefix: 327 coefficients, 317 positive relative H5 margins, 317 positive Q6 signs, 0 inconclusive, 39 repaired coefficients, 1 open analytic tail, weakest n=316 above 7/1000",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-six lambda=-100 tail curvature reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order6_m100_tail_curvature_reduction.py",),
        expected=(
            "validated order-six tail curvature reduction: 8 rows, 0 issues, 2 exact factorizations, 1 open curvature target",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF inverse-seventh-power first-summand dominance",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power7_dominance_extension.py",),
        expected=(
            "validated power-seven first-summand dominance: 6 rows, 0 issues, 11 positive gates, tail k>=316",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-six first/full curvature bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order6_first_summand_curvature_bridge.py",),
        expected=(
            "validated order-six first/full curvature bridge: 8 rows, 0 issues, 1 full transfer, 1 open continuous target",
        ),
        category="exact theorem composition",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF compound order-six high-cumulant corridor",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.py",),
        expected=(
            "validated order-six high-cumulant corridor: 2 formal orders, 109 terms, 0 issues, 2 exact corridors",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-six nested-curvature compact certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order6_nested_curvature_compact_certificate.py",),
        expected=(
            "validated order-six nested compact certificate: 38 blocks, 0 issues, 1 compact theorem, 1 open ray",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-six nested-curvature finite ray",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.py",),
        expected=(
            "validated order-six nested finite-ray certificate: 17999 ray blocks, 0 issues, 1 theorem, 1 open asymptotic ray",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-six nested-curvature asymptotic ray",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.py",),
        expected=(
            "validated order-six nested curvature asymptotic certificate: 7 rows, 0 issues, 2 normalized-H theorems, 1 dimensionless interval theorem, 1 asymptotic curvature theorem, 0 open rows, scaled upper 2.27683610696345567703247070312471784171958323101989922910682E+1",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-six lambda=-100 entry",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order6_m100_entry_certificate.py",),
        expected=(
            "validated order-six lambda=-100 entry certificate: 10 rows, 0 issues, 1 continuous curvature theorem, 1 complete scalar ceiling, 1 analytic tail theorem, 1 all-shift entry theorem, 0 open rows",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-six uniform-heat forward invariance",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate.py",),
        expected=(
            "validated order-six uniform-heat forward invariance: 8 rows, 0 issues, 1 contiguous theorem, 1 arbitrary-column theorem, 1 open order-seven handoff",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF graded-kernel all-order Vandermonde lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.py",),
        expected=(
            "validated graded-kernel all-order Vandermonde lemma: 12 order rows, 46233 permutation stress cases, 169 coefficient-valuation cells, 1 all-fixed-order eventual signed-tail theorem, 0 issues",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-seven uniform-tail and flow reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order7_uniform_tail_flow_reduction.py",),
        expected=(
            "validated order-seven uniform tail and flow reduction: 9 rows, 0 issues, 1 universal-tail specialization, 1 condensation coordinate, 1 cooperative recursion, 1 conditional forward theorem, 1 lower-cone countermodel, 1 open lambda=-100 entry",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-seven lambda=-100 prefix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order7_m100_prefix_certificate.py",),
        expected=(
            "validated order-seven lambda=-100 prefix: 327 coefficients, 317 positive Q6 values, 315 positive relative Q6 margins, 315 positive Q7 signs, 0 inconclusive, 12 repaired coefficients, 1 open analytic tail, weakest n=314 above 9/1000",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-seven lambda=-100 tail curvature reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order7_m100_tail_curvature_reduction.py",),
        expected=(
            "validated order-seven tail curvature reduction: 7 rows, 0 issues, 1 exact factorization, 1 exact curvature reduction, 1 positive comparison, 1 conditional tail theorem, 1 open curvature target",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF inverse-eighth-power rebalanced first-summand dominance",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.py",),
        expected=(
            "validated power-eight rebalanced first-summand dominance: 7 rows, 14 positive gates, tail k>=300, 0 issues",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-seven first/full curvature bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order7_first_summand_curvature_bridge.py",),
        expected=(
            "validated order-seven first/full curvature bridge: 9 rows, 0 issues, 1 fourth-gap floor theorem, 1 full transfer, 3 conditional theorems, 1 open continuous target, degree 102, 103 positive coefficients",
        ),
        category="exact theorem composition",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF compound order-seven shifted-jet compact bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.py",),
        expected=(
            "validated order-seven shifted-jet t=320..1000 certificate: 5 rows, 186 contiguous blocks, 11 shifts per anchor, 0 issues",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-seven nested-curvature compact certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order7_nested_curvature_compact_certificate.py",),
        expected=(
            "validated order-seven nested compact certificate: 82 blocks, 0 issues, 1 compact theorem, 1 open finite ray",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-seven high-cumulant coarse corridor",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order7_high_cumulant_coarse_corridor.py",),
        expected=(
            "validated order-seven high-cumulant corridor: 72 formal terms, 0 issues, 2 exact corridors",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-seven nested-curvature finite-ray certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.py",),
        expected=(
            "validated order-seven nested finite-ray certificate: 17999 ray blocks, 0 issues, 1 theorem, 1 open asymptotic ray",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-seven nested-curvature asymptotic-ray certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate.py",),
        expected=(
            "validated order-seven nested curvature asymptotic certificate: 7 rows, 0 issues, 2 normalized-H theorems, 1 dimensionless interval theorem, 1 asymptotic curvature theorem, 0 open rows",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-seven lambda=-100 entry",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order7_m100_entry_certificate.py",),
        expected=(
            "validated order-seven lambda=-100 entry certificate: 11 rows, 0 issues, 1 continuous curvature theorem, 1 complete scalar ceiling, 1 analytic tail theorem, 1 all-shift entry theorem, 0 open rows",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-seven uniform-heat forward invariance",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate.py",),
        expected=(
            "validated order-seven uniform-heat forward invariance: 8 rows, 0 issues, 1 contiguous theorem, 1 arbitrary-column theorem, 1 open all-order handoff",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eight uniform-tail and flow reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order8_uniform_tail_flow_reduction.py",),
        expected=(
            "validated order-eight uniform tail and flow reduction: 9 rows, 0 issues, 1 universal-tail specialization, 1 condensation coordinate, 1 cooperative recursion, 1 conditional forward theorem, 1 lower-cone countermodel, 1 open lambda=-100 entry",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eight lambda=-100 prefix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order8_m100_prefix_certificate.py",),
        expected=(
            "validated order-eight lambda=-100 prefix: 1257 coefficients, 1245 positive Q7 values, 1243 positive relative Q7 margins, 1243 positive Q8 signs, 0 inconclusive, 12 repaired coefficients, 1 open analytic tail, weakest n=1242 above 1/300",
        ),
        category="finite interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eight lambda=-100 tail-curvature reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order8_m100_tail_curvature_reduction.py",),
        expected=(
            "validated order-eight tail curvature reduction: 7 rows, 0 issues, 1 exact factorization, 1 exact curvature reduction, 1 positive comparison, 1 conditional tail theorem, 1 open curvature target",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF inverse-ninth-power first-summand dominance",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power9_rebalanced_dominance_extension.py",),
        expected=(
            "validated power-nine rebalanced first-summand dominance: 6 rows, 0 issues, 14 positive analytic gates, 1 dominance theorem",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eight first/full curvature bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order8_first_summand_curvature_bridge.py",),
        expected=(
            "validated order-eight first/full curvature bridge: 9 rows, 0 issues, 1 fifth-gap floor theorem, 134 positive transfer coefficients, 1 full-kernel transfer theorem, 1 open continuous target",
        ),
        category="exact theorem composition",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF compound order-eight high-cumulant coarse corridor",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order8_high_cumulant_coarse_corridor.py",),
        expected=(
            "validated order-eight high-cumulant corridor: 0 formal terms, 0 issues, 2 exact corridors",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eight shifted-jet certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate.py",),
        expected=(
            "validated order-eight shifted-jet certificate: 185 blocks, 0 issues, 1 continuous theorem, largest scaled upper 8.81335404534865833076786551392429128216127678108895244813038E+2, 1 open ray",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eight nested-curvature compact certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order8_nested_curvature_compact_certificate.py",),
        expected=(
            "validated order-eight nested compact certificate: 96 blocks, 0 issues, largest scaled upper 3.99455119863063431363867738935234812098029095541283895880362E+3, weakest U lower 1.22933158930953976688262799664500547709540427726036454298503E-4, 1 open finite ray",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eight nested-curvature finite-ray certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate.py",),
        expected=(
            "validated order-eight nested finite-ray certificate: 17999 ray blocks, 0 issues, 1 theorem, 1 open asymptotic ray",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eight nested-curvature asymptotic-ray certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate.py",),
        expected=(
            "validated order-eight nested curvature asymptotic certificate: 7 rows, 0 issues, 2 normalized-H theorems, 1 dimensionless interval theorem, 1 asymptotic curvature theorem, 0 open rows, scaled upper 1.34489839907184243202209472656240124460170043090398255726938E+2",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eight lambda=-100 entry",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order8_m100_entry_certificate.py",),
        expected=(
            "validated order-eight lambda=-100 entry certificate: 11 rows, 0 issues, 1 continuous curvature theorem, 1 analytic tail theorem, 1 all-shift entry theorem, 0 open rows",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eight uniform-heat forward invariance",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate.py",),
        expected=(
            "validated order-eight uniform-heat forward invariance: 8 rows, 0 issues, 1 contiguous theorem, 1 arbitrary-column theorem, 1 open all-order handoff",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine uniform-tail and flow reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_uniform_tail_flow_reduction.py",),
        expected=(
            "validated order-nine uniform tail and flow reduction: 9 rows, 0 issues, 1 universal-tail specialization, 1 condensation coordinate, 1 cooperative recursion, 1 conditional forward theorem, 1 lower-cone countermodel, 1 open lambda=-100 entry",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine lambda=-100 prefix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_m100_prefix_certificate.py",),
        expected=(
            "validated order-nine lambda=-100 prefix: 1257 coefficients, 1243 positive Q8 values, 1241 positive relative Q8 margins, 1241 positive Q9 signs, 0 inconclusive, 38 repaired coefficients, 1 open analytic tail, weakest n=1240 above 1/250",
        ),
        category="finite interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine lambda=-100 tail-curvature reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_m100_tail_curvature_reduction.py",),
        expected=(
            "validated order-nine tail curvature reduction: 7 rows, 0 issues, 1 exact factorization, 1 exact curvature reduction, 1 positive comparison, 1 conditional tail theorem, 1 open curvature target",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine first/full curvature bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_first_summand_curvature_bridge.py",),
        expected=(
            "validated order-nine first/full curvature bridge: 10 rows, 0 issues, 1 sixth-gap floor theorem, 169 positive transfer coefficients, 1 full-kernel transfer theorem, 1 open continuous target, 2 finite splice indices",
        ),
        category="exact theorem composition",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine high-cumulant coarse corridor",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_high_cumulant_coarse_corridor.py",),
        expected=(
            "validated order-nine high-cumulant corridor: 0 formal terms, 0 issues, 2 exact corridors",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine exact-point H0-H8 cache",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_shifted_point_h0_h8_cache.py",),
        expected=(
            "validated order-nine exact-point H0-H8 cache: 8929 rows, 0 issues",
        ),
        category="finite interval certificates",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine localized lower bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_localized_lower_bridge_certificate.py",),
        expected=(
            "validated order-nine localized lower bridge: 279 root segments, 874 accepted blocks, 0 issues, largest scaled upper 4.19424425522037111044561248832644227467685288026319528153785E+3",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine nested-curvature compact certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_nested_curvature_compact_certificate.py",),
        expected=(
            "validated order-nine nested compact certificate: 108 blocks, 0 issues, largest scaled upper 4.19918548221032914747349727361016012002522663395440335471778E+3, weakest V lower 1.43232298117533817321514966026794097632441289481887521956705E-4, 2 open handoff ranges",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine nested-curvature finite-ray certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate.py",),
        expected=(
            "validated order-nine finite-ray certificate: 17999 blocks, 0 issues, largest scaled upper 2.17882197094171936230311013389392316198614992290990750440229E+3, weakest V lower 5.39596935966714165720214930394256250853232758887469547152339E+0, 1 open asymptotic ray",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine nested-curvature asymptotic-ray certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate.py",),
        expected=(
            "validated order-nine asymptotic-ray certificate: 0 issues, scaled upper 3.24905922088046558201313018798812324136265244865916789195364E+2<500<4200, 6 positive stable coordinates, 1 global above-5700 composition",
        ),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine global first-summand curvature",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_first_summand_curvature_certificate.py",),
        expected=(
            "validated global order-nine first-summand curvature: 0 issues, largest scaled upper 4.19918548221032914747349727361016012002522663395440335471778E+3",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine finite endpoint splice",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_m100_finite_splice_certificate.py",),
        expected=(
            "validated order-nine finite splice: 1259 coefficients, 0 issues, 2 splice rows, 1243 combined positive signs",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine lambda=-100 entry",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_m100_entry_certificate.py",),
        expected=(
            "validated all-shift signed order-nine entry at lambda=-100: 0 issues",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-nine uniform-heat forward invariance",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate.py",),
        expected=(
            "validated signed contiguous and arbitrary-column order nine on -100<=lambda<=0: 0 issues",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten first/full curvature bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_first_summand_curvature_bridge.py",),
        expected=("validated order-ten first/full curvature bridge:", "scaled transfer", "0 issues"),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten high-cumulant corridor",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_high_cumulant_coarse_corridor.py",),
        expected=("validated order-ten high-cumulant corridor:", "0 issues"),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten compact H2-H24 cache",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_compact_h2_h24_unit_cache.py",),
        expected=("validated order-ten compact H2-H24 cache: 32336 contiguous unit tiles, 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten sparse exact H0-H23 cache",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache.py",),
        expected=("validated order-ten sparse exact H0-H23 cache: 4042 step-eight points, 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten propagated H0-H7 cache",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_compact_propagated_point_h0_h7_cache.py",),
        expected=("validated order-ten propagated H0-H7 cache: 32335 integer-grid points, 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten localized lower bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_localized_lower_bridge_certificate.py",),
        expected=("validated order-ten localized lower bridge: 279 segments, 9,996 contiguous blocks, 0 issues",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten compact curvature certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_nested_curvature_compact_certificate.py",),
        expected=("validated order-ten compact curvature certificate: 18310 contiguous blocks, 0 issues",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten finite-ray curvature certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_nested_curvature_finite_ray_certificate.py",),
        expected=("validated order-ten finite ray: 17999 blocks, 0 issues",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten asymptotic-ray curvature certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_nested_curvature_asymptotic_ray_certificate.py",),
        expected=("validated order-ten asymptotic ray: 0 issues",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten global first-summand curvature",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_first_summand_curvature_certificate.py",),
        expected=("validated global order-ten first-summand curvature certificate: 1 half-line first-summand theorem, 0 full-kernel claims, 0 RH claims",),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten endpoint-tail curvature reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_m100_tail_curvature_reduction.py",),
        expected=("validated order-ten endpoint-tail curvature reduction: 10 rows, 0 issues",),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten finite endpoint splice",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_m100_finite_splice_certificate.py",),
        expected=("validated order-ten finite splice:", "0 issues"),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten delayed lambda=-100 entry",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_m100_delayed_entry_certificate.py",),
        expected=("validated delayed signed order-ten entry at lambda=-100: 0 issues",),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF delayed cooperative heat-tail lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_delayed_cooperative_heat_tail_lemma.py",),
        expected=("validated delayed cooperative heat-tail lemma: 1 shifted theorem, order-ten and order-eleven specializations, 0 issues",),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten lambda-zero prefix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_lambda0_prefix_certificate.py",),
        expected=("validated order-ten lambda-zero prefix: 22 coefficients, 4 positive Q10 rows, 4 source overlaps, 0 issues",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-ten lambda-zero completion",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order10_lambda0_completion_certificate.py",),
        expected=("validated order-ten lambda-zero completion: all shifts, arbitrary columns through order ten, 0 issues",),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eleven lambda=-100 prefix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order11_m100_prefix_certificate.py",),
        expected=("validated order-eleven lambda=-100 prefix: 1263 coefficients, 1243 positive Q11 rows, 5 direct audits, 0 issues",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eleven lambda-zero prefix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order11_lambda0_prefix_certificate.py",),
        expected=("validated order-eleven lambda-zero prefix: 24 coefficients, 4 positive Q11 rows, 4 source overlaps, 0 issues",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eleven first-summand point scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order11_first_summand_point_scout.py",),
        expected=("validated order-eleven first-summand point scout: 8 points, 72 positive coordinates", "0 issues"),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eleven localized final-gap core",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order11_localized_final_gap_interval_core.py",),
        expected=("validated order-eleven localized final-gap core: H24 derivative budget, 9 positive point coordinates, 9 positive localized coordinates, 0 issues",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eleven lower sparse exact H0-H23 source",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_compound_order11_lower_sparse_point_h0_h23_cache.py",
            "--require-complete",
        ),
        expected=("validated order-eleven lower sparse exact H0-H23 cache (complete): 2233/2233 rows, 0 issues",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eleven sparse H23 propagation",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order11_sparse_h0_h23_propagation.py",),
        expected=("validated order-eleven sparse H23 propagation: 5 polynomial translations + 5 H24 remainder enclosures, 0 issues",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eleven shifted Taylor-model algebra",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order11_shifted_taylor_model_algebra.py",),
        expected=("validated order-eleven shifted Taylor-model algebra: 9 products + 81 stable-log enclosures, 0 issues",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eleven sparse H23 lower-bridge pilots",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_compound_order11_sparse_h23_lower_bridge_pilot.py",
            "work/rh_compute/results/jensen_window_pf_compound_order11_sparse_h23_lower_bridge_pilot_t1252.json",
            "work/rh_compute/results/jensen_window_pf_compound_order11_sparse_h23_lower_bridge_pilot_t1500.json",
            "work/rh_compute/results/jensen_window_pf_compound_order11_sparse_h23_lower_bridge_pilot_t2200.json",
            "work/rh_compute/results/jensen_window_pf_compound_order11_sparse_h23_lower_bridge_pilot_t3000.json",
        ),
        expected=("validated order-eleven sparse H23 lower-bridge pilots: 4 cells, 8 quarter blocks, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF compound order-eleven curvature bridge target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_compound_order11_curvature_bridge_target.py",),
        expected=("validated order-eleven curvature bridge target: 18 envelope rows, scaled transfer", "<37, 0 issues"),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF all-order endpoint-to-heat reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_all_order_endpoint_heat_reduction.py",),
        expected=(
            "validated all-order endpoint-to-heat reduction: 13 rows, 0 issues, 4 symbolic orders, 255 parity checks, 1 endpoint/interval equivalence, 1 arbitrary-column consequence, 1 rejected m>=10 endpoint hierarchy, 1 separate Jensen/PF bridge",
        ),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF endpoint deep-Schur coordinate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_endpoint_deep_schur_coordinate.py",),
        expected=(
            "validated endpoint deep-Schur coordinate: 14 rows, 0 issues, 4096 rectangle checks, 984 arbitrary-column checks, 1023 inverse checks, 1 rigorous endpoint PF_3 counterexample, 4 negative deep rectangles, 1 rejected m>=10 rectangle hierarchy, 1 separate Jensen/PF bridge",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF deep-Schur Toda/boundary gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_deep_schur_toda_boundary_gate.py",),
        expected=(
            "validated deep-Schur Toda/boundary gate: 15 rows, 0 issues, 4 symbolic Toda checks, 251 boundary checks, 138 strict-Schur checks, 1 exact full-PF Jensen counterexample, 4 negative order-ten rectangles, 1 rejected rectangle hierarchy, 1 Xi-specific bridge",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF endpoint order-ten counterexample",
        command=("work/rh_compute/scripts/check_jensen_window_pf_endpoint_order10_counterexample.py",),
        expected=(
            "validated endpoint order-ten counterexample: 8 rows, 0 issues, 5 direct checks, 4 negative deep rectangles, 1237 positive scanned rows, 0 inconclusive, 1 rejected all-order endpoint hierarchy, 1 surviving Xi/Phi bridge target",
        ),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF order-moment transport fit gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_order_moment_transport_fit_gate.py",),
        expected=(
            "validated Jensen-window PF order-moment transport fit gate: 7 rows, 0 issues, 1 exact reparametrization, 1 positive origin derivative obstruction, 1 orientation mismatch, 1 forbidden promotion, 1 open signed-transport handoff",
        ),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF structural ansatz matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py",),
        expected=("validated Jensen-window PF structural ansatz matrix: 6 ansatz rows, 0 issues, 0 ready-to-apply rows",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Schur shape contract",
        command=("work/rh_compute/scripts/check_jensen_window_pf_schur_shape_contract.py",),
        expected=("validated Jensen-window PF Schur shape contract: 4 grid rows, 0 issues, 2 frontier rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF column recurrence contract",
        command=("work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_contract.py",),
        expected=("validated Jensen-window PF column recurrence contract: 4 degree rows, 0 issues, 2 hard frontier rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF column recurrence finite coverage",
        command=("work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_finite_coverage.py",),
        expected=("validated Jensen-window PF column recurrence finite coverage: 1470 direct positive rows, 210 hard recurrence rows, 315 Sturm/PF windows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Arb Jensen-window column recurrence stress",
        command=("work/rh_compute/scripts/check_arb_jensen_window_column_recurrence_stress.py",),
        expected=("validated 12600 Arb Jensen-window column recurrence stress rows with 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF reciprocal coefficient extended stress",
        command=("work/rh_compute/scripts/check_jensen_window_pf_reciprocal_coefficient_extended_stress.py",),
        expected=("validated Jensen-window PF reciprocal coefficient extended stress: 72600 rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF reciprocal positivity route matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_reciprocal_positivity_route_matrix.py",),
        expected=("validated Jensen-window PF reciprocal positivity route matrix: 9 rows, 0 issues, 0 ready-to-apply rows",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF reciprocal fraction scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_reciprocal_fraction_scout.py",),
        expected=("validated Jensen-window PF reciprocal fraction scout: 3 symbolic rows, 735 finite rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF reciprocal signed J-fraction scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_reciprocal_signed_j_fraction_scout.py",),
        expected=("validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF reciprocal signed Jacobi beta scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_reciprocal_signed_jacobi_beta_scout.py",),
        expected=("validated Jensen-window PF reciprocal signed Jacobi beta scout: 3 symbolic rows, 3675 beta rows, 2940 positive rows, 630 negative rows, 105 terminal-zero rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF reciprocal Motzkin path obstruction scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.py",),
        expected=("validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF reciprocal Motzkin parity-lift obstruction scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.py",),
        expected=("validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF signed J-fraction theorem target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_signed_j_fraction_theorem_target.py",),
        expected=("validated Jensen-window PF signed J-fraction theorem target: 7 fit rows, 0 issues, 0 ready-to-apply rows",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF modified signed-model target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_modified_signed_model_target.py",),
        expected=("validated Jensen-window PF modified signed-model target: 9 model rows, 0 issues, 0 ready-to-apply rows, 4 live modified candidates",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF oscillatory resolvent fit matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_oscillatory_resolvent_fit_matrix.py",),
        expected=("validated Jensen-window PF oscillatory resolvent fit matrix: 8 fit rows, 0 issues, 0 ready-to-apply rows",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF positive readout theorem target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_positive_readout_theorem_target.py",),
        expected=("validated Jensen-window PF positive readout theorem target: 8 candidate rows, 0 issues, 0 ready-to-apply rows, 2 live foundational routes",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF positive spectral moment obstruction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_positive_spectral_moment_obstruction.py",),
        expected=("validated Jensen-window PF positive spectral moment obstruction: 3 symbolic rows, 735 finite Delta2 obstruction rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF nonordinary positive transform ansatz matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.py",),
        expected=("validated Jensen-window PF nonordinary positive transform ansatz matrix: 8 ansatz rows, 0 issues, 0 ready-to-apply rows, 3 live ansatz rows",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF nonpower functional low-degree scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_nonpower_functional_low_degree_scout.py",),
        expected=("validated Jensen-window PF nonpower functional low-degree scout: 7 scout rows, 0 issues, 0 ready-to-apply rows, 1 live contract rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF nonpower functional cone candidate matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_nonpower_functional_cone_candidate_matrix.py",),
        expected=("validated Jensen-window PF nonpower functional cone candidate matrix: 8 cone rows, 0 issues, 0 ready-to-apply rows, 2 live cone rows",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Cauchy-Binet cone frontier matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_cone_frontier_matrix.py",),
        expected=("validated Jensen-window PF Cauchy-Binet cone frontier matrix: 8 frontier rows, 0 issues, 0 ready-to-apply rows, 2 live frontier rows",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF monotone contraction frontier scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_frontier_scout.py",),
        expected=("validated Jensen-window PF monotone contraction frontier scout: 2 exact rows, 88 Bernstein coefficients, 210 finite zeta rows, 0 issues",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF monotone-contraction column extension scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_column_extension_scout.py",),
        expected=(
            "validated Jensen-window PF monotone-contraction column extension scout: 25 column rows, 3329 Bernstein coefficients, 3 beyond-frontier rows, 0 negative Bernstein rows, 0 issues",
        ),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF monotone-contraction sparse degree-6 scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree6_scout.py",),
        expected=(
            "validated Jensen-window PF monotone-contraction sparse degree-6 scout: 10 degree-6 rows, 63347 Bernstein coefficients, m<=10, 0 negative Bernstein rows, 0 zero Bernstein rows, 0 issues",
        ),
        category="exact theorem-search algebra",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF monotone-contraction sparse degree-7 frontier scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.py",),
        expected=(
            "validated Jensen-window PF monotone-contraction sparse degree-7 frontier scout: 9 positive rows, 1 certificate-obstruction row, 932691 Bernstein coefficients, first obstruction m=10, 126 negative Bernstein coefficients, 0 zero Bernstein coefficients, 0 issues",
        ),
        category="exact theorem-search algebra",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF monotone-contraction sparse degree-7 subdivision scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.py",),
        expected=(
            "validated Jensen-window PF monotone-contraction sparse degree-7 subdivision scout: 3 dyadic slabs, 785400 slab Bernstein coefficients, 0 negative slab coefficients, 0 zero slab coefficients, repaired m=10 obstruction, 0 issues",
        ),
        category="exact theorem-search algebra",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF monotone-contraction all-m counterexample",
        command=("work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_all_m_counterexample.py",),
        expected=(
            "validated Jensen-window PF monotone-contraction all-m counterexample: degree 7, m=11, exact full-cone witness, 6 lower walls, negative normalized value, 0 issues",
        ),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF monotone contraction theorem target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_theorem_target.py",),
        expected=("validated Jensen-window PF monotone contraction theorem target: 9 rows, 0 issues, 0 ready-to-apply rows, 2 live routes",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF heat-flow monotone closure scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_heat_flow_monotone_closure_scout.py",),
        expected=("validated Jensen-window PF heat-flow monotone closure scout: 4 exact rows, 315 threshold rows, 305 flow-bracket rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF heat-flow boundary threshold lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_heat_flow_boundary_threshold_lemma.py",),
        expected=("validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF kernel Mellin upper-wall certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_kernel_mellin_upper_wall_certificate.py",),
        expected=("validated Jensen-window PF kernel Mellin upper-wall certificate: 8 rows, 0 issues, 200 positive compact intervals, 1 positive analytic ray, 1 remaining open cone clause, 0 ready-to-apply rows",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF log-concave Mellin monotone-wall countermodel",
        command=("work/rh_compute/scripts/check_jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.py",),
        expected=("validated Jensen-window PF log-concave Mellin monotone-wall countermodel: 6 rows, 0 issues, 2 upper-wall contractions, 1 monotone-wall violation",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF T=1156 monotone-wall counterexample certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.py",),
        expected=("validated Jensen-window PF T=1156 monotone-wall counterexample certificate: 7 rows, 0 issues, 4 coefficient enclosures, 1 zeta monotone-wall violation",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda kernel summand-shift lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.py",),
        expected=("validated Jensen-window PF negative-lambda kernel summand-shift lemma: 8 rows, 0 issues, 6 exact rows, 1 compact interval row, 1 open far-tail row, 0 ready-to-apply rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda first-summand dominance certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_dominance_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda first-summand dominance certificate: 10 rows, 0 issues, 4 exact rows, 5 interval rows, 15 positive analytic gates, 1 open dominant-wall row, 0 ready-to-apply rows",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda -100 k320 collar extension certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda -100 k320 collar extension certificate: 6 rows, 0 issues, 76 positive coefficients, 74 cone rows, 73 adjacent-wall rows, 19 new extension rows, 0 ready-to-apply rows",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda -100 full cone-entry certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda -100 full cone-entry certificate: 8 rows, 0 issues, 321 positive coefficients, 319 pointwise cone rows, 318 adjacent rows, 1 analytic adjacent tail, 1 open flow handoff, 3 ready-to-apply rows",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF infinite heat-flow cone-invariance certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.py",),
        expected=("validated Jensen-window PF heat-flow infinite cone-invariance certificate: 8 rows, 0 issues, 1 infinite maximum principle, 1 full cone propagation, 1 endpoint cone theorem, 1 open all-order handoff, 3 ready-to-apply rows",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF defect complete-monotonicity scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_defect_complete_monotonicity_scout.py",),
        expected=("validated Jensen-window PF defect complete-monotonicity scout: 3284 defect positives, 3288 log positives, 838 inconclusive, both certified through order 8, 5 lambdas, 1 exact all-shape countermodel, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF multiplier complete-monotonicity frontier scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.py",),
        expected=("validated Jensen-window PF multiplier complete-monotonicity frontier scout: 7980 positive intervals, 0 inconclusive, orders 0..55, 5 lambdas, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF multiplier Hausdorff-uniqueness bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_multiplier_hausdorff_uniqueness_bridge.py",),
        expected=("validated Jensen-window PF multiplier Hausdorff-uniqueness bridge: 10 rows, 0 issues, 1 Hausdorff measure theorem, 1 unit-atomic characterization, 1 interpolation guard, 1 open recovery handoff",),
        category="exact structural lemmas",
    ),
    GateSpec(
        name="Jensen-window PF multiplier leading-atom bound certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_multiplier_leading_atom_bound_certificate.py",),
        expected=("validated Jensen-window PF multiplier leading-atom bound certificate: 8 rows, 0 issues, 56 difference orders, beta6 in (4.863538496,4.863538497), alpha_min>4.863538496, N(11/2)<=1, 1 open existence handoff",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF multiplier unit-atomic obstruction certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_multiplier_unit_atomic_obstruction_certificate.py",),
        expected=("validated Jensen-window PF multiplier unit-atomic obstruction certificate: 8 rows, 0 issues, 1 atom cutoff, 1 ratio cap, 1 unit-atomic route ruled out, 0 open requirements",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF heat-flow Jensen hierarchy lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_heat_flow_jensen_hierarchy_lemma.py",),
        expected=("validated Jensen-window PF heat-flow Jensen hierarchy lemma: 9 rows, 0 issues, 5 exact hierarchy identities, 2 cubic countermodels, 1 open higher-minor handoff, 0 ready-to-apply rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF cubic reciprocal-defect invariance lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.py",),
        expected=("validated Jensen-window PF cubic reciprocal-defect invariance lemma: 10 rows, 0 issues, 6 exact coordinate rows, 1 conditional maximum principle, 318 lambda=-100 prefix margins, 310 nonnegative-grid margins, 0 failed or inconclusive, 1 open tail handoff",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF cubic lambda=-100 tail-entry certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_cubic_m100_tail_entry_certificate.py",),
        expected=("validated Jensen-window PF cubic lambda=-100 tail-entry certificate: 10 rows, 0 issues, 4074 negative-skewness blocks, 1 analytic negative-skewness ray, 318 prefix margins, 1 all-k cubic tail, 1 full cubic entry theorem, 1 open forward-uniform tail",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF cubic forward-uniform tail certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_cubic_forward_uniform_tail_certificate.py",),
        expected=("validated Jensen-window PF cubic forward-uniform tail certificate: 10 rows, 0 issues, 3 exact flow identities, 1 weighted source cap, 1 initial weighted tail, 1 forward-uniform tail, 1 full cubic propagation theorem, 1 lambda=0 cubic theorem, 1 open higher-degree handoff",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF quartic boundary-flow obstruction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_quartic_boundary_flow_obstruction.py",),
        expected=("validated Jensen-window PF quartic boundary-flow obstruction: 10 rows, 0 issues, 4 exact quartic identities, 1 hyperbolic boundary point, 4 positive ratio margins, 3 strict cubic margins, 1 negative quartic derivative, 1 blocked promotion, 1 open coupled-invariant handoff",),
        category="countermodel guards",
    ),
    GateSpec(
        name="Jensen-window PF strong-log-concave local quartic countermodel",
        command=("work/rh_compute/scripts/check_jensen_window_pf_strong_logconcave_local_quartic_countermodel.py",),
        expected=("validated Jensen-window PF strong-log-concave local quartic countermodel: 10 rows, 0 issues, 6 positive Mellin moments, 4 local ratio-wall contractions, 3 monotone gaps, 3 strict cubic margins, 1 negative quartic frontier, 1 full-support approximation theorem, 1 Xi-specific handoff",),
        category="countermodel guards",
    ),
    GateSpec(
        name="Jensen-window PF quartic double-root threshold lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_quartic_double_root_threshold_lemma.py",),
        expected=("validated Jensen-window PF quartic double-root threshold lemma: 10 rows, 0 issues, 5 exact coordinate identities, 1 double-root splitting criterion, 1 branch-aware inward threshold, 1 triple-root equality, 1 tangent factor, 1 explained countermodel, 1 open global-invariant handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF quartic signed-Hankel branch-exclusion lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.py",),
        expected=("validated Jensen-window PF quartic signed-Hankel branch-exclusion lemma: 10 rows, 0 issues, 4 exact identities, 1 imported Xi Hankel theorem, 2 excluded boundary strata, 1 reduced outer threshold, 1 countermodel separation, 1 open outer threshold",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF quartic outer-threshold order-four nonpromotion gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.py",),
        expected=("validated quartic outer-threshold order-four nonpromotion gate: 11 rows, 0 issues, 1 exact outer quartic contact, 8 strict ratio coordinates, 7 strict scaled-defect steps, 7 reciprocal-defect increment bounds, 7 strict cubic frontiers, 120 signed order-two minors, 126 signed order-three minors, 56 signed order-four minors, 1 failed outer threshold, 1 negative quintic discriminant, 1 forbidden promotion, 1 closed downstream handoff",),
        category="countermodel guards",
    ),
    GateSpec(
        name="Jensen-window PF quartic outer-branch length-13 obstruction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_branch_length13_obstruction.py",),
        expected=("validated quartic outer-branch length-13 obstruction: 10 rows, 0 issues, 3 exact corridor parameters, 3 derivative certificates, 1516 derivative Bernstein coefficients, 8 denominator-factor Bernstein coefficients, 1 negative cube-corner maximum, 1 forbidden all-length promotion, 1 repaired handoff",),
        category="countermodel guards",
    ),
    GateSpec(
        name="Jensen-window PF quartic alternate length-13 survivor gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.py",),
        expected=("validated alternate quartic length-13 survivor gate: 14 rows, 0 issues, 8 exact corridor parameters, 12 contractions, 11 strict adjacent scalar steps, 10 positive order-three gaps, 8 positive order-four margins, 364 signed order-two minors, 715 signed order-three minors, 792 signed order-four minors, 1 positive length-13 compatibility, 1 negative fixed-tail length-14 compatibility, 0 uniform length-14 theorems, 1 repaired scope handoff",),
        category="countermodel guards",
    ),
    GateSpec(
        name="Jensen-window PF quartic one-contact uniform length-14 interval certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_branch_uniform_length14_interval_certificate.py",),
        expected=("validated one-contact uniform quartic length-14 interval certificate: 9 rows, 0 issues, 37474 certified leaves, 37473 dyadic splits, maximum depth 23, 8 tail parameters, 1 scaled wall, 1 uniform fixed-contact obstruction, 0 uniform all-contact theorems",),
        category="interval theorem certificates",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF quartic outer-contact normal-form gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_contact_normal_form_gate.py",),
        expected=("validated quartic outer-contact normal-form gate: 13 rows, 0 issues, 3 contact variables, 1 excluded outer branch, 4 normalized defects, 2 initial gap scalings, 2 first-extension regimes, 1 normalized Delta_14, 1 adjacent-quintic discriminant factorization, 1 low-q negative-discriminant collar, 0 uniform all-contact theorems",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF quartic outer-contact length-14 survivor gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_contact_length14_survivor_gate.py",),
        expected=("validated quartic outer-contact length-14 survivor gate: 15 rows, 0 issues, 3 contact variables, 9 corridor parameters, 13 contractions, 15 coefficients, 4043 positive supported order-two through order-eight minors, 1 positive Delta_14, 1 positive Delta_15, 1 nonhyperbolic adjacent quintic, 0 uniform all-contact obstructions",),
        category="countermodel guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman zero-slab degree-71 sector certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_zero_slab_degree71_sector_certificate.py",),
        expected=("validated Newman zero-slab degree-71 sector certificate: 13 rows, 0 issues, 4 directed integrals, 1 positive slab margin, 1 zero-free complex slab, 3 published inputs, all shifts through degree 71, 0 all-degree theorems",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman real-zero-band degree-361 sector certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.py",),
        expected=("validated Newman real-zero-band degree-361 sector certificate: 12 rows, 0 issues, 200 certified boundary boxes, 0 unresolved, 1 real-zero band, all shifts through degree 361, 0 all-degree theorems, weakest negative margin",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF quartic-quintic polar-contact lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_quartic_quintic_polar_contact_lemma.py",),
        expected=("validated Jensen-window PF quartic-quintic polar-contact lemma: 10 rows, 0 issues, 4 exact polar identities, 1 strict nonroot test, 1 multiplicity rule, 1 double-to-triple theorem, 1 quintic contact factorization, 1 cofactor gate, 1 open all-degree handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF cofinal-degree polar-closure lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_cofinal_degree_polar_closure_lemma.py",),
        expected=("validated Jensen-window PF cofinal-degree polar-closure lemma: 10 rows, 0 issues, 3 exact polar identities, 1 interlacing theorem, 1 multiplicity theorem, 1 finite-tower closure, 1 cofinal-degree closure, 1050 finite Sturm rows, 2875 contraction-only rows, 1 open terminal handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF cofinal scaling-limit equivalence gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_cofinal_scaling_limit_equivalence_gate.py",),
        expected=("validated Jensen-window PF cofinal scaling-limit equivalence gate: 10 rows, 0 issues, 2 exact scaling identities, 2 analytic limit steps, 1 cofinal-to-LP theorem, 1 LP-to-all-degrees theorem, 1 fixed-shift equivalence, 2 non-promotion guards, 1 open independent-route handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF polar heat-collision cascade lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_polar_heat_collision_cascade_lemma.py",),
        expected=("validated Jensen-window PF polar heat-collision cascade lemma: 10 rows, 0 issues, 3 exact local identities, 1 double-root criterion, 1 higher-multiplicity gate, 1 infinite polar cascade, 1 exponential-polynomial classification, 1 unbounded-degree escape theorem, 1 open scaled-tail handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF fixed-root shift-tangency rigidity gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_fixed_root_shift_tangency_rigidity_gate.py",),
        expected=("validated fixed-root shift-tangency rigidity gate: 32 rows, 84 independent derivatives, 21 adjacent shifts, 42 transfer jets, 33 persistent checks, 10 Hankel checks, 75 Newton identities, 1 lambda=-100 anchor, 1 heat-interval application, 0 zeta near-chain bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF first-boundary eventual-Hankel escape gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_first_boundary_eventual_hankel_escape_gate.py",),
        expected=("validated first-boundary eventual-Hankel escape gate: 15 rows, 9 sources, 11 effective-order records, 190 collision pairs, 760 propagation-index checks, 190 global rows covered, 0 global rows uncovered, 139 local all-shift rows, 6 local delayed rows, 45 local open rows, 0 uniform-order thresholds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF scaled double-zero boundary-layer lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_scaled_double_zero_boundary_layer_lemma.py",),
        expected=("validated Jensen-window PF scaled double-zero boundary-layer lemma: 10 rows, 0 issues, 3 exact scaling identities, 1 heat PDE, 1 double-zero transversality, 1 universal boundary layer, 1 D^(-3/2) gap law, 1 external-field D^(-2) collision law, 1 exact toy family, 1 threshold-exhaustion theorem, 1 open uniform handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman arbitrary-multiplicity Jensen boundary-layer gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_arbitrary_multiplicity_jensen_boundary_layer_gate.py",),
        expected=("validated arbitrary-multiplicity Jensen boundary-layer gate: 15 rows, 7 sources, 15 multiplicities, 79 universal coefficients, 15 Hermite checks, 15 imaginary-side checks, 15 orientation checks, 14 exact polynomial layers, 1 fixed shift, 0 uniform-degree bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman arbitrary-multiplicity secondary cubic layer gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_arbitrary_multiplicity_secondary_cubic_layer_gate.py",),
        expected=("validated secondary cubic Jensen layer gate: 16 rows, 4 sources, 14 multiplicities, 200 universal coefficients, 12 independent finite models, 1 cubic commutator, 1 exact m=3 threshold, 2 open handoffs, 0 Xi-specific signs, 0 degree-uniform bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman arbitrary-multiplicity tertiary universal layer gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_arbitrary_multiplicity_tertiary_universal_layer_gate.py",),
        expected=("validated tertiary universal Jensen layer gate: 15 rows, 3 sources, 14 multiplicities, 14 independent operators, 12 independent finite models, 28 Appell identities, 12 Xi-unit isolations, 2 exact m3 corrections, 2 open handoffs, 0 degree-uniform bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman arbitrary-multiplicity first Xi-jet layer gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_arbitrary_multiplicity_first_xi_jet_layer_gate.py",),
        expected=("validated first Xi-jet Jensen layer gate: 16 rows, 3 sources, 14 multiplicities, 28 independent operators, 12 independent finite models, 28 Appell identities, 12 full Xi-jet isolations, 2 exact m3 corrections, 2 sign countermodels, 2 open handoffs, 0 degree-uniform bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman regular-field first-loss nonclosure gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_regular_field_first_loss_nonclosure_gate.py",),
        expected=("validated regular-field first-loss nonclosure gate: 16 rows, 3 sources, 10 arbitrary-field models, 15 multiplicities, 10 positive-time samples, 15 forward Hermite checks, 15 backward nonreal checks, 2 open handoffs, 0 Xi field bounds, 0 degree-uniform bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Fourier-moment regular-field gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_fourier_moment_regular_field_gate.py",),
        expected=("validated Fourier-moment regular-field gate: 20 rows, 5 sources, 15 multiplicities, 15 direct moment checks, 15 score moment checks, 15 score-jet checks, 15 sinc models, 5 Gaussian models, 5 strong-curvature checks, 15 Hermite checks, 2 open handoffs, 0 Xi field bounds, 0 degree-uniform bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman modular-primitive contact nonclosure gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_modular_primitive_contact_nonclosure_gate.py",),
        expected=("validated modular-primitive contact nonclosure gate: 21 rows, 4 sources, 17 operator jets, 15 contact multiplicities, 15 sinc-Bessel models, 3 Bessel diagnostics, 8 modular odd jets, 15 Hermite checks, 2 open handoffs, 0 Xi field bounds, 0 degree-uniform bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman discrete-theta contact ladder scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_discrete_theta_contact_ladder_scout.py",),
        expected=(
            "validated discrete theta contact ladder scout: 4 rosters, 24 high-precision field witnesses, 2 unbounded-ratio guards, 2 refined target witnesses",
            "0 issues",
        ),
        category="finite theorem-search diagnostics",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman fixed-theta mass-profile dual scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_fixed_theta_mass_profile_dual_scout.py",),
        expected=(
            "validated fixed-theta mass-profile dual scout:",
            "closure TV=1.23179830387e-16",
            "0 issues",
        ),
        category="finite theorem-search diagnostics",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman fourth theta-summand tangency homotopy scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fourth_summand_tangency_homotopy_scout.py",),
        expected=(
            "validated fourth-summand tangency homotopy: 4 independent anchors",
            "0 issues",
        ),
        category="finite theorem-search diagnostics",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman fourth theta-summand tangency interval continuation certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fourth_summand_tangency_interval_continuation_certificate.py",),
        expected=(
            "checked fourth-summand interval continuation: charts=11, connectors=10, independent_anchors=3, crossing=pass, 0 issues",
        ),
        category="interval certificates",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman fourth theta-summand compact contact degree certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fourth_summand_compact_contact_degree_certificate.py",),
        expected=(
            "checked compact contact degree certificate: boundary_cells=352, winding=-1, independent_anchors=4, unique_contact=True, 0 issues",
        ),
        category="interval certificates",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman fourth theta-summand positive-time degree certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fourth_summand_positive_time_degree_certificate.py",),
        expected=(
            "checked positive-time degree certificate: boundary_cells=1248, winding=0, independent_anchors=4, no_contacts=True, 0 issues",
        ),
        category="interval certificates",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman fifth theta-summand positive-time degree certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fifth_summand_positive_time_degree_certificate.py",),
        expected=(
            "checked fifth-summand positive-time degree certificate: boundary_cells=32, winding=0, independent_anchors=4, one_slab=True, no_contacts=True, 0 issues",
        ),
        category="interval certificates",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman complete-theta local positive-time degree certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_complete_tail_positive_time_degree_certificate.py",),
        expected=(
            "checked complete-theta local degree certificate: tail_orders=2, tail_to_margin<4e-21, complete_theta=True, no_contacts=True, 0 issues",
        ),
        category="interval certificates",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman complete-theta positive-time spatial tile degree certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_complete_spatial_tile_degree_certificate.py",),
        expected=(
            "checked complete-theta spatial tile degree certificate: boundary_cells=64, tile_windings=[0, 0], independent_anchors=8, combined_band=[135,136.5], no_contacts=True, 0 issues",
        ),
        category="interval certificates",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman root external-field lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_root_external_field_lemma.py",),
        expected=("validated Jensen-window PF Newman root external-field lemma: 10 rows, 0 issues, 5 exact canonical-product identities, 1 pair-flow reduction, 1 gap-stiffness expansion, 2 sign countermodels, 1 cosine equilibrium benchmark, 1 open Xi-balance handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman classical-field balance gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_classical_field_balance_gate.py",),
        expected=("validated Jensen-window PF Newman classical-field balance gate: 10 rows, 0 issues, 3 exact field identities, 1 arithmetic equilibrium, 1 continuum -pi/8 benchmark, 1 quantile-drift match, 1 published fixed-time localization theorem, 2 exact sensitivity countermodels, 1 compactness reduction, 1 open lambda-uniform handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman local odd-count reduction lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_local_odd_count_reduction_lemma.py",),
        expected=("validated Jensen-window PF Newman local odd-count reduction lemma: 10 rows, 0 issues, 3 exact Stieltjes identities, 1 explicit outer bound, 1 log-squared localization theorem, 1 odd-count formula, 3 finite reciprocal-gap checks, 1 published uniform counting input, 1 exact classical-field birth countermodel, 1 open Xi collision-exclusion handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman boundary-energy direction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_boundary_energy_direction_gate.py",),
        expected=("validated Jensen-window PF Newman boundary-energy direction gate: 10 rows, 0 issues, 1 universal gap law, 1 exact higher-jet birth model, 1 nonintegrable collision-energy asymptotic, 1 conditional exclusion criterion, 1 published-scope audit, 1 open Xi boundary-energy handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman counterfactual birth-signature atlas",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_counterfactual_birth_signature_atlas.py",),
        expected=("validated Newman counterfactual birth-signature atlas: 12 rows, 0 issues, 1 symmetric quartet heat collision, 1 exact Jensen threshold, 1 exact Li quartet law, 1 Suzuki shift law, 1 diagonal finite-sensor escape theorem, 2 open global handoffs",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman cofinal boundary-degree transfer target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.py",),
        expected=("validated Newman cofinal boundary-degree transfer target: 12 rows, 0 issues, 2 orientation identities, 1 vector boundary-Rouche theorem, 1 sign-definite degree composition, 1 Q207 base, 2 explicit C1 proxy budgets, 1 endpoint-uniformity countermodel, 1 open cofinal Xi boundary theorem",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Q208 selected-boundary pilot",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_q208_selected_boundary_pilot.py",),
        expected=("validated Q208 selected-boundary pilot: 8 panels, 46 certified leaves, 0 unresolved, 0 global promotions",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF Newman convex phase-cell unwrapping lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_convex_phase_cell_unwrapping_lemma.py",),
        expected=("validated Newman convex phase-cell unwrapping lemma: 10 rows, 0 issues, 5 projection audits, 4 exact polygon tests, 3 rejection guards, 20 Q208 right-edge cells, 1 open bottom-edge target",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Q208 bottom phase-cell certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_q208_bottom_phase_cell_certificate.py",),
        expected=("validated Q208 bottom phase-cell certificate: 492/492 panels, 492 cells, 0 unresolved, complete=True, 0 structural issues",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman Q208 top phase-cell certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_q208_top_phase_cell_certificate.py",),
        expected=("validated Q208 top phase-cell certificate: 492/492 panels, 492 cells, 0 unresolved, complete=True, 0 structural issues",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman Q208 closed-boundary winding certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_q208_closed_boundary_winding_certificate.py",),
        expected=("validated Q208 closed-boundary winding certificate: 1005 cells, 1005 exact witnesses, winding=0, Q1..Q208 certified, 0 cofinal promotions",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF Newman cofinal phase-cell scaling diagnostics",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_cofinal_phase_cell_scaling_diagnostics.py",),
        expected=("validated cofinal phase-cell scaling diagnostics: 414 Q207 comparisons, 492+492 Q208 cells, 408/414 and 470/492 branch agreements, crossing=-19, 0 cofinal promotions",),
        category="diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF Newman adiabatic phase-cell successor lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma.py",),
        expected=("validated Newman adiabatic phase-cell successor lemma: 10 rows, 3 exact lemmas, 2 conditional compositions, 1 finite Q208 calibration, 1 open all-j Xi antecedent",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Q207-Q208 adiabatic bottom collar certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_q207_q208_adiabatic_bottom_collar_certificate.py",),
        expected=("validated Q207-Q208 adiabatic bottom collar: 380/490 panels, 379 certified, 1 failed, complete=False, 0 all-j promotions",),
        category="diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF Newman Q207-Q208 refined-tail adiabatic collar certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_q207_q208_adiabatic_bottom_collar_refined_tail_certificate.py",),
        expected=("validated hybrid Q207-Q208 adiabatic collar: 379 coarse + 222/222 refined panels, 0 failures, complete=True, 0 all-j promotions",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman Q207-Q208 forward adiabatic successor certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_q207_q208_forward_adiabatic_successor_certificate.py",),
        expected=("validated forward Q207-Q208 adiabatic successor: 414 old cells, 525/525 transport panels, max ratio<1, 1 finite successor, 0 all-j promotions",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman single-carrier adiabatic benchmark",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_single_carrier_adiabatic_benchmark.py",),
        expected=("validated Newman single-carrier adiabatic benchmark: 10 rows, 8 exact benchmark reductions, 1 four-carrier obstruction, 1 open Xi small-ball handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman corrected crossing slope-gap reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_crossing_slope_gap_reduction.py",),
        expected=("validated Newman crossing slope-gap reduction: 10 rows, 8 exact reductions, 1 exact countermodel, 4 corrected-crossing diagnostics, 1 open Xi target",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman multiplicity-compatible adaptive-jet benchmark",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_multiplicity_compatible_adaptive_jet_benchmark.py",),
        expected=("validated Newman multiplicity-compatible adaptive-jet benchmark: 10 rows, 9 exact reductions, 1 open two-regime Xi handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman time-dependent scaled-jet successor lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_time_dependent_scaled_successor_lemma.py",),
        expected=("validated Newman time-dependent scaled-jet successor lemma: 10 rows, 9 exact reductions, 1 open two-regime Xi antecedent",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman ray-aligned parabolic-frequency reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_ray_aligned_parabolic_frequency_reduction.py",),
        expected=("validated Newman ray-aligned parabolic-frequency reduction: 14 rows, 11 exact reductions, 1 asymptotic composition, 0 open new-strip antecedents, 1 open Xi descendant theorem",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Q207-Q208 parabolic-frequency relative diagnostic",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_q207_q208_parabolic_frequency_relative_diagnostic.py",),
        expected=("validated Q207-Q208 parabolic-frequency relative diagnostic: 525 panels, 4 rigorous finite scale certificates, 0 all-j promotions",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman parabolic-frequency contact-normal hierarchy gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_parabolic_frequency_contact_normal_hierarchy_gate.py",),
        expected=("validated Newman parabolic-frequency contact-normal hierarchy gate: 17 rows, 12 exact identities, 3 exact route guards, 1 literature guard, 1 open Xi target",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman parabolic-frequency energy/current gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_parabolic_frequency_energy_current_gate.py",),
        expected=("validated Newman parabolic-frequency energy/current gate: 18 rows, 8 exact identities, 3 coefficient bounds, 2 pointwise bridge rows, 3 nonpromotion guards, 1 open Xi bulk-margin target, 0 pointwise contact exclusions",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman ray-aligned energy endpoint-margin gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_ray_aligned_energy_endpoint_margin_gate.py",),
        expected=("validated ray-aligned energy endpoint-margin gate: 18 rows, 10 compact edge records, 1 left raw margin, 1 right normalized margin, 1 exact contact floor, 1 open full-collar bulk budget, 2 nonpromotion guards, 0 pointwise contact exclusions",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman ray-aligned energy cell-localization gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_ray_aligned_energy_cell_localization_gate.py",),
        expected=("validated ray-aligned energy cell-localization gate: 20 rows, 2 exact local balances, 2 exact contact floors, 2 conditional no-contact criteria, 2 outer anchors, 0 cofinal internal anchor registries, 1 finite successor instance, 2 nonpromotion guards, 0 pointwise contact exclusions",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 critical C1 endpoint peeling contract",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract.py",),
        expected=("validated Newman Polymath-15 critical C1 endpoint peeling contract: 12 rows, 5 exact coefficient identities, 2 primary-source inputs, 1 corrected-main definition, 2 conditional transfers, 1 route guard, 1 open quantitative target",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 critical Dirichlet first correction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_dirichlet_first_correction_gate.py",),
        expected=("validated Newman Polymath-15 critical Dirichlet first correction gate: 10 rows, 4 exact identities, 2 primary-source inputs, 2 analytic guards, 1 signed-main definition, 1 open uniform target",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 critical Dirichlet second-order remainder certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_dirichlet_second_order_remainder_certificate.py",),
        expected=("validated Newman critical Dirichlet second-order remainder certificate: 12 rows, C_D=400000, fixed-cell C1 constant 8000000, 1 open global splice guard",),
        category="analytic certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 critical first-order global remainder certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_first_order_global_remainder_certificate.py",),
        expected=("validated Newman critical first-order global remainder certificate: 15 rows, eta_0=100000, eta_1=200000, 1 open signed-contact target",),
        category="analytic certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 critical first-order signed-contact reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_first_order_signed_contact_reduction.py",),
        expected=("validated Newman critical first-order signed-contact reduction: 15 rows, eta_0=100000, eta_1=200000, 1 open frequency target",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order target reconciliation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_target_reconciliation_gate.py",),
        expected=("validated Newman first-order target reconciliation gate: 22 rows, 1 superseded handoff, 1 cutoff-uniform first-order remainder, 4 exact energy identities, 2 contact-box coordinates, 1 radial target, 1 box-optimal target, 1 normalized Abel target, 2 nonpromotion guards, 0 pointwise contact exclusions",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order rectangular boundary-degree reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_rectangular_boundary_degree_reduction.py",),
        expected=("validated Newman first-order rectangular boundary-degree reduction: 22 rows, 2 error-box coordinates, 2 rectangular homotopies, 1 positive-index degree transfer, 1 boundary-only Abel target, 2 route branches, 2 nonpromotion guards, 0 contact exclusions",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 centered ray-bottom logarithmic-flow reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_ray_bottom_logarithmic_flow_reduction.py",),
        expected=("validated Newman centered ray-bottom logarithmic-flow reduction: 24 rows, 1 exact cutoff partition, 8 logarithmic-flow identities, 2 phase-flux identities, 1 orientation transfer, 1 cutoff-join contract, 2 nonpromotion guards, 2 open Xi targets, 0 Abel gaps, 0 horizontal phase bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 geometric prime-power phase-monotonicity gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_geometric_prime_power_phase_monotonicity_gate.py",),
        expected=("validated geometric prime-power phase-monotonicity gate: 18 rows, 1 phase-derivative identity, 3 uniform monotone prime families, 1 complete short-chain recrossing guard, 2 nonpromotion guards, 1 open C1 perturbation target, 0 joined Abel gaps, 0 successor winding bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 prime-power logarithmic phase-flow gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_logarithmic_phase_flow_gate.py",),
        expected=("validated prime-power logarithmic phase-flow gate: 20 rows, 7 exact logarithmic-phase identities, 1 dimensionless heat profile, 1 complete dyadic counterexample, 1 rejected absolute C1 target, 1 replacement starlikeness target, 0 joined Abel gaps, 0 successor winding bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 prime-power heat-starlikeness base certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_base_certificate.py",),
        expected=("validated prime-power heat-starlikeness base certificate: 18 rows, 3 exact Bernstein base certificates, 1 analytic p>=5 two-level family, 3 actual ray-positive base families, 0 length-propagation theorems, 0 joined Abel gaps, 0 successor winding bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 prime-power heat-starlikeness length propagation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_length_propagation_gate.py",),
        expected=("validated prime-power heat-starlikeness length propagation gate: 20 rows, 1 append-one identity, 1 exact dyadic next-length certificate, 1 analytic all-length p>=5 family, 2 actual ray-positive propagated families, 2 open small-prime length families, 0 joined Abel gaps, 0 successor winding bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 ternary all-length prime-power heat-starlikeness gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_ternary_length_gate.py",),
        expected=("validated ternary all-length prime-power heat-starlikeness gate: 18 rows, 5 small-d clusters, 12 finite energies, 1 analytic terminal tail, 1 ideal all-length ternary family, 1 actual all-length ternary family, 1 open small-prime length family, 0 joined Abel gaps, 0 successor winding bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 dyadic all-length Fourier-defect localization gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_dyadic_defect_gate.py",),
        expected=("validated dyadic all-length Fourier-defect localization gate: 13 rows, 2 initial clusters, 1 shortest-length certificate, 2 low terminal anchors, 12 middle terminal clusters, 1 analytic terminal tail, 1 proved positive-offset family, 2 possible defects, 0 ideal all-length dyadic families, 0 actual all-length dyadic families, 0 joined Abel gaps, 0 successor winding bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 dyadic all-length prime-power heat-starlikeness gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_dyadic_length_gate.py",),
        expected=("validated dyadic all-length prime-power heat-starlikeness gate: 18 rows, 5 finite cores, 1 initial core, 5 finite terminal guards, 1 analytic terminal tail, 1 ideal all-length dyadic family, 1 actual all-length dyadic family, complete dyadic minimum length 8, 0 remaining defects, 0 joined Abel gaps, 0 successor winding bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 short small-prime heat-starlikeness counter-gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_starlikeness_short_family_counter_gate.py",),
        expected=("validated short small-prime heat-starlikeness counter-gate: 10 rows, 8 exact interior negative witnesses, 6 short dyadic lengths rejected, 2 short ternary lengths rejected, 0 uniform short-block positive families, 1 joined-phase handoff, 0 joined Abel gaps, 0 successor winding bounds",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 joined phase-current polarization gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_joined_phase_current_polarization_gate.py",),
        expected=("validated joined phase-current polarization gate: 10 rows, 1 exact polarization, 1 exact chain insertion, 1 cubic contact, 1 negative joined current, 3 long-family inputs, 0 joined cross-current bounds, 0 joined Abel gaps, 0 successor winding bounds",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 joined pair-kernel and p-free rejoin gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_joined_pair_kernel_pfree_rejoin_gate.py",),
        expected=("validated joined pair-kernel p-free rejoin gate: 10 rows, 1 exact pair kernel, 1 moment collapse, 1 heat telescope, 1 quadratic rejoin, 1 cross-cancellation guard, 0 cross-current signs, 0 joined Abel gaps, 0 successor winding bounds",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 corrected Mangoldt-Abel contact gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_abel_contact_gate.py",),
        expected=("validated corrected Mangoldt-Abel contact gate: 12 rows, 1 correction polynomial, 4 logarithmic moments, 3 symmetric Mangoldt moments, 8 prime-edge witnesses, 24 indefinite minors, 0 Abel gaps, 0 winding bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 contact-centered Mangoldt-Abel equivalence gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_contact_centering_gate.py",),
        expected=("validated contact-centered Mangoldt gate: 12 rows, 1 physical moment, 2 contact fibres, 1 Abel duality, 8 balanced hyperbola audits, 0 mismatches, 0 Abel gaps, 0 winding bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 centered Mangoldt adjacent-cutoff transport gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_adjacent_cutoff_transport_gate.py",),
        expected=("validated centered Mangoldt cutoff transport gate: 12 rows, 1 endpoint cancellation, 79 transitions, 71 ordinary, 8 square, 0 mismatches, 0 Abel gaps, 0 winding bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 centered Mangoldt endpoint-composed Vaughan gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_mangoldt_endpoint_composed_vaughan_gate.py",),
        expected=("validated centered Mangoldt endpoint-composed Vaughan gate: 14 rows, 1 convolution identity, 214 transitions, 5 cube transfers, 1 square/cube overlap, 2 exact countermodels, 0 signed lower bounds, 0 Abel gaps, 0 winding bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 complex endpoint source-normalization corrigendum",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_complex_endpoint_source_normalization_gate.py",),
        expected=("validated complex endpoint source-normalization gate: 16 rows, 1 nonreal source witness, 4 corrected Cartesian projections, 2 contact minors, 2 rank guards, 7 historical artifacts quarantined, 0 signed lower bounds, 0 Abel gaps, 0 winding bounds",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 endpoint-relative phase-current recurrence gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_endpoint_relative_phase_current_recurrence_gate.py",),
        expected=("validated endpoint-relative phase-current recurrence gate: 16 rows, 2 exact endpoint currents, 1 division-free relative current, 2 asymptotic sign witnesses, 1 terminal-tail recurrence, 1 contact-minor sum, 0 uniform terminal signs, 0 Abel gaps, 0 winding bounds",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 real-edge projective-current gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_real_edge_projective_current_gate.py",),
        expected=("validated real-edge projective-current gate: 14 rows, 1 division-free edge current, 4 asymptotic edge jets, 3072 Arb intervals, 2 removable charts, 1 strict curvature margin, 1 uniform q=1 leading sign, 0 finite-height edge signs, 0 Abel gaps, 0 winding bounds",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 q=1 finite-height real-edge remainder gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_q1_finite_height_real_edge_remainder_gate.py",),
        expected=("validated q=1 finite-height real-edge remainder gate: 16 rows, 4096 Arb boxes, 2 Cauchy charts, 4 finite edge jets, remainder <1/10000000, 1 retained-model finite-height sign, 0 q>1 signs, 0 Xi-level signs",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 critical-ray finite-height real-edge gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_critical_ray_finite_height_real_edge_gate.py",),
        expected=("validated critical-ray finite-height real-edge gate: 13 rows, 4096 Arb boxes, 27 normalized majorants, 4 finite edge jets, 1 full critical-ray model sign, 1 q>=1 model sign, 0 cutoff splices, 0 Xi-level signs",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 adjacent-cutoff real-edge projective splice gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_adjacent_cutoff_real_edge_projective_splice_gate.py",),
        expected=("validated adjacent-cutoff real-edge projective splice gate: 12 rows, 4 exact endpoint values, 1 exact leading wedge, 1 rational remainder budget, 1 nonvanishing affine join, 1 signed q>=1 model splice, 0 second adjacent x-jets, 0 Xi-level splices",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 residual-placement C1 cumulative handoff gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_residual_placement_c1_cumulative_handoff_gate.py",),
        expected=("validated residual-placement C1 cumulative handoff gate: 16 rows, 2 C1 residual coordinates, 1 allocation-sign guard, 0 C2 boundary requirements, 1 whole-jet homotopy, 1 cumulative polarization, 1 open cross-current/Abel target",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 edge-anchored near-terminal pair-current guard",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_edge_anchored_near_terminal_pair_current_guard.py",),
        expected=("validated edge-anchored near-terminal pair-current guard: 12 rows, 1 phase limit, 4 carrier jets, 4 edge jets, 2 Arb witnesses, 1 asymptotic sign reversal, 0 uniform termwise absorptions, 0 aggregate counterexamples",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 contiguous terminal-tail recurrence/current gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_recurrence_current_gate.py",),
        expected=("validated contiguous terminal-tail recurrence/current gate: 16 rows, 1 exact C0 tail collapse, 3 grouped jets, 1024 Arb boxes, 1 all-fixed-M clockwise theorem, 0 growing-tail estimates, 0 aggregate closures",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 contiguous terminal-tail growing-prefix finite-height gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_growing_prefix_finite_height_gate.py",),
        expected=("validated growing-prefix finite-height gate: 16 rows, 1 exact ratio, 1 differentiated ratio, 4 carrier errors, 4 aggregate errors, exponent 1/18, M=3 at L=50, current <-1/400, 0 bulk closures",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 contiguous terminal-tail anchored Abel cross-current reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_abel_cross_current_reduction.py",),
        expected=("validated terminal-tail anchored Abel cross-current reduction: 14 rows, 2 Abel identities, 2 current identities, 1 signed Phi_B target, 0 bulk closures",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 contiguous terminal-tail anchored five-current bulk closure",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_current_bulk_closure_reduction.py",),
        expected=("validated terminal-tail five-current bulk closure: 12 rows, 5 complex currents, 3 derivative identities, 1 algebraic closure, 0 signed bulk bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 contiguous terminal-tail anchored five-current Mangoldt normal form",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_current_mangoldt_normal_form_gate.py",),
        expected=("validated five-current Mangoldt normal form: 12 rows, 5 correction-free moments, 4 Mangoldt orders, 4 indefinite minor families, 0 signed Type-I/II bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 contiguous terminal-tail explicit five-moment Phi quadratic reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_explicit_phi_quadratic_reduction.py",),
        expected=("validated explicit five-moment Phi_B quadratic reduction: 13 rows, 4 real observations, rank <=4, generic inertia (2,2,6), 1 double-contact obstruction, 0 signed Phi_B bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 contiguous terminal-tail five-moment two-carrier kernel reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_two_carrier_kernel_reduction.py",),
        expected=("validated five-moment two-carrier kernel reduction: 14 rows, 1 endpoint-linear kernel, 1 Hermitian kernel, 1 transpose kernel, 1 physical chi_N bound, 0 signed Type-I/II bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 contiguous terminal-tail logarithmic-phase Turan reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_logarithmic_phase_turan_reduction.py",),
        expected=("validated logarithmic-phase Turan reduction: 15 rows, 5 Fourier jet identities, 1 Bochner kernel, 1 leading Turan identity, 1 zero-fibre sign, 1 integer-log guard, 0 signed joint bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 physical q=1 saddle-phase variance reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_q1_saddle_phase_variance_reduction.py",),
        expected=("validated physical q=1 saddle-phase variance reduction: 16 rows, 4 parameter laws, 5 saddle transfers, 1 two-carrier sign, 1 phase discriminant, 1 ordered-profile guard, 0 signed physical bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 physical q=1 saddle-phase quadratic transfer barrier",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_q1_saddle_phase_quadratic_transfer_barrier.py",),
        expected=("validated q=1 quadratic transfer barrier: 15 rows, 1 discriminant collapse, 1 signed line integral, 1 terminal amplitude lower bound, 1 absolute-mass barrier, 0 signed physical bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 endpoint-composed six-moment saddle-flow reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_saddle_flow_reduction.py",),
        expected=("validated six-moment saddle-flow reduction: 15 rows, 1 full-current flow, 2 flow symmetries, 1 leading pair identity, 6 moments, 1 new Mangoldt order, 0 signed flow bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment flow-matrix reciprocal-phase reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_flow_matrix_phase_reduction.py",),
        expected=("validated six-moment flow matrix/phase reduction: 22 rows, rank <=8, generic inertia (4,4,4), 3 reciprocal phase families, 1 Hermitian diagonal null, 0 signed flow bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Hermitian reciprocal-pairing reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_hermitian_reciprocal_pairing_reduction.py",),
        expected=("validated Hermitian reciprocal pairing reduction: 17 rows, 1 swap identity, 1 lattice-defect bound, 0 Poisson remainder bounds, 0 signed flow bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment finite-Poisson transport reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_finite_poisson_transport_reduction.py",),
        expected=("validated six-moment finite-Poisson/transport reduction: 23 rows, 6 moments, 5 transport recurrences, 2 cutoff jump laws, 0 explicit uniform remainder constants, 0 signed flow bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Morse-Fresnel endpoint reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_endpoint_reduction.py",),
        expected=("validated six-moment Morse-Fresnel endpoint reduction: 27 rows, 6 exact transition integrals, 6 transition majorants, 3 stable endpoint sums, 6 tail majorants, 0 evaluated physical constants, 0 signed flow bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Morse-Fresnel aggregate-scaling gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_aggregate_scaling_gate.py",),
        expected=("validated Morse-Fresnel aggregate-scaling gate: 24 rows, 1 transition collar, 0 enumerated modes, majorant floor 91/96800, 6 h^2 tail bounds, 0 grouped roster bounds, 0 signed flow bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Morse-Fresnel endpoint-coherence Abel gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_endpoint_coherence_abel_gate.py",),
        expected=("validated Morse-Fresnel endpoint-coherence Abel gate: 24 rows, 2 phase differences, phase sums <3sqrt(alpha), 1 h-scale terminal subblock, 2 coherent endpoint phases, 0 endpoint-composed bounds, 0 signed flow bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Morse-Fresnel endpoint-composition retention gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_endpoint_composition_retention_gate.py",),
        expected=("validated Morse-Fresnel endpoint-composition retention gate: 24 rows, 6 composed identities, 2 roster representations, 7 retained endpoint floors, 1 observation expansion, 0 oscillatory interior bounds, 0 signed flow bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Morse-Fresnel observation-image compression gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_observation_image_compression_gate.py",),
        expected=("validated Morse-Fresnel observation-image compression gate: 25 rows, 8 visible residual observations, 1 degree-five linear functional, 4 quadratic pairings, 1 tail projection template, 0 interior bounds, 0 signed flow bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Morse-Fresnel physical P_lin amplitude gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_amplitude_gate.py",),
        expected=("validated physical P_lin/Morse-amplitude gate: 24 rows, 2 correction channels, 6 centered coefficients, degree 5, 1 total-value leading factor, 2 c-prime formulas, 0 grouped interior bounds, 0 signed flow bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Morse-Fresnel physical P_lin correction-coefficient envelope gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_correction_coefficient_envelope_gate.py",),
        expected=("validated physical correction-coefficient envelope gate: 16 rows, 2 quotient bounds, 6 centered coefficients, 6 coefficient envelopes, 1 nonvanishing constant channel, 0 P_lin bounds, 0 signed flow bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Morse-Fresnel physical P_lin joined-coefficient envelope gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_joined_coefficient_envelope_gate.py",),
        expected=("validated physical P_lin joined-coefficient envelope gate: 17 rows, 5 rate bounds, 3 quadratic-rate bounds, 6 propagated p_j, 0 numerical row bounds, 0 signed flow bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Morse-Fresnel physical P_lin basis-conjugation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_basis_conjugation_gate.py",),
        expected=("validated physical P_lin basis-conjugation gate: 18 rows, 6 basis coefficients, 6 translated weights, 3 joined identities, 0 ideal-cubic bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Morse-Fresnel physical P_lin ideal-cubic nilpotent transport gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_nilpotent_transport_gate.py",),
        expected=("validated ideal-cubic nilpotent transport gate: 18 rows, nilpotence orders 3/4, 4 cubic coefficients, 5 transport identities, 2 collar bounds, 0 grouped bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Morse-Fresnel physical P_lin ideal-cubic reciprocal gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_reciprocal_gate.py",),
        expected=("validated ideal-cubic Fresnel reciprocal gate: 20 rows, 4 weighted transports, 1 Fubini reindexing, 6 reciprocal phase identities, 1 absolute-family barrier, 0 grouped bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Morse-Fresnel ideal-cubic reciprocal first-correction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_reciprocal_first_correction_gate.py",),
        expected=("validated reciprocal first-correction gate: 18 rows, 4 cell identities, 5 dual-Morse identities, 5 first-correction identities, 1 cancellation, 0 finite-cell bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment Morse-Fresnel ideal-cubic reciprocal finite-cell inversion gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_reciprocal_finite_cell_inversion_gate.py",),
        expected=("validated reciprocal finite-cell inversion gate: 23 rows, 3 cell-Poisson identities, 3 carrier inversions, 1 tie cancellation, 3 global recompositions, 0 exterior-tail bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment reciprocal terminal full-support homotopy gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_ideal_cubic_reciprocal_terminal_full_support_homotopy_gate.py",),
        expected=("validated reciprocal terminal full-support homotopy gate: 27 rows, 1 all-carrier homotopy, 3 terminal extensions, 1 terminal strip, 0 outer-complement bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment full-support outer endpoint kernel gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_outer_endpoint_kernel_gate.py",),
        expected=("validated full-support outer endpoint kernel gate: 27 rows, 1215 kernel points, 270 S1 checks, 3645 kernel bounds, 1020 removable checks, 0 complete outer bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment full-support terminal endpoint obstruction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_terminal_endpoint_obstruction_gate.py",),
        expected=("validated full-support terminal endpoint obstruction gate: 24 rows, 585 S1 checks, 1 row collapse, 1 order-one obstruction, 0 signed kernel bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment full-support terminal conditional quadrature gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_terminal_conditional_quadrature_gate.py",),
        expected=("validated full-support terminal conditional quadrature gate: 28 rows, 630 kernel checks, 2 carrier halves, 1 quadrature defect, 0 signed bulk-terminal bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment full-support terminal-bulk phase completion gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_terminal_bulk_phase_completion_gate.py",),
        expected=("validated full-support terminal-bulk phase completion gate: 28 rows, 2 phase families, 1 coefficient completion, 1 Hermitian null, 1 transpose collapse, 0 signed completed bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment full-support terminal/nonterminal relative-lift gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_terminal_nonterminal_relative_lift_gate.py",),
        expected=("validated full-support terminal/nonterminal relative-lift gate: 30 rows, 4 remainder components, 2 relative lifts, 2 degree-five functionals, 4 tie invariants, 528 independent geometry checks, 0 signed lower/interior bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment full-support joined lower/interior Abel recombination gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_joined_lower_interior_abel_recombination_gate.py",),
        expected=("validated joined lower/interior Abel recombination gate: 28 rows, 2 physical factorizations, 4 global recombinations, 4 relative joins, 256 independent physical checks, 64 independent Abel checks, 0 signed defect bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment full-support finite-band Dirichlet-cell reduction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_finite_band_dirichlet_cell_reduction_gate.py",),
        expected=("validated finite-band Dirichlet-cell reduction gate: 33 rows, 1900 independent kernel checks, 21 independent cell checks, 2 relative reductions, 8096 independent far-gap checks, 5120 independent Fourier witnesses, 0 signed physical bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment full-support ideal-cubic reciprocal collar transport gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_ideal_cubic_reciprocal_collar_transport_gate.py",),
        expected=("validated ideal-cubic reciprocal collar transport gate: 31 rows, 230 independent collar checks, 3310 independent far compressions, 2387 independent tie checks, 1519 independent nonclosure witnesses, 864 independent phase checks, 4 numerical IBP checks, 6 numerical moment checks, 0 signed collar bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment full-support ideal-cubic double-Morse rail-flux gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_ideal_cubic_double_morse_rail_flux_gate.py",),
        expected=("validated ideal-cubic double-Morse rail-flux gate: 28 rows, 90 independent endpoint checks, 90 independent interface checks, 6 independent symbolic checks, 5 tensor flux checks, 180 nonzero shear checks, 80 tie gauges, 0 signed ideal flux bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment full-support ideal-cubic double-Morse hyperbolic-transport gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_ideal_cubic_double_morse_hyperbolic_transport_gate.py",),
        expected=("validated ideal-cubic double-Morse hyperbolic-transport gate: 26 rows, 24 independent profile checks, 8 independent transport checks, 3 independent hyperbolic checks, 12 independent transverse kernels, 5 ideal checks, 10 ridge sign witnesses, 0 signed ideal transport bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment full-support ideal-cubic double-Morse ridge/carrier inversion gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_ideal_cubic_double_morse_ridge_carrier_inversion_gate.py",),
        expected=("validated double-Morse ridge/carrier inversion gate: 28 rows, 5 independent inversion checks, 3 independent ridge checks, 9 independent witness checks, 7 endpoint-weight checks, 1 nonzero exterior witness, 0 signed ridge-completed bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 six-moment full-support ideal-cubic symmetric outer-pairing gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_six_moment_full_support_ideal_cubic_symmetric_outer_pairing_gate.py",),
        expected=("validated ideal-cubic symmetric outer-pairing gate: 34 rows, 5 recombinations, 105 independent roster checks, 9 independent means, 6 independent remote pairs, 6 endpoint identities, 2 reflected-phase checks, 3 conjugation checks, 0 signed complete ideal bounds",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order cofinal boundary reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_cofinal_boundary_reduction.py",),
        expected=("validated Newman first-order cofinal boundary reduction: 13 rows, full budget 50000000000, 3 open boundary targets",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order oriented successor-winding reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_oriented_successor_winding_reduction.py",),
        expected=("validated Newman oriented successor-winding reduction: 14 rows, 1 exact chain audit, 2 exact crossing reductions, 1 one-sided integer trap, 3 open cofinal obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order Wronskian crossing reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_wronskian_crossing_reduction.py",),
        expected=("validated Newman first-order Wronskian crossing reduction: 12 rows, eta_0=100000, eta_1=200000, 2 exact crossing classes, 3 open Xi obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered real-residual reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_real_residual_reduction.py",),
        expected=("validated Newman first-order centered real-residual reduction: 15 rows, |u_a|<e^-L, |v_a|<3/x^2, core error <1e-6*e^-5L/4, 2 open Xi cell obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered complex-zero scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_complex_zero_scout.py",),
        expected=("validated Newman first-order centered complex-zero scout: N=6, |E_[1]|<1e-70, |U|>0.5, |det D_(x,t)E_[1]|>0.3, cross-precision drift <1e-50",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered adjacent-saddle recurrence",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_adjacent_saddle_recurrence.py",),
        expected=("validated Newman first-order centered adjacent-saddle recurrence: 13 rows, exact Delta A identity, vanished a^-1 coefficient, 18 selected L>=50 rows with e^-7L/4 scalar scaling, 2 open scalar obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered adjacent-chart stability certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_adjacent_chart_stability_certificate.py",),
        expected=("validated Newman first-order centered adjacent-chart stability: 14 rows, |rho-1|<11/a^2, |rho_x|<3/a^3, |Delta A|<5000*e^-7L/4, 1 open bulk Xi obligation",),
        category="asymptotic theorem-search certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered bulk pair-transfer gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_bulk_pair_transfer_gate.py",),
        expected=("validated Newman centered bulk pair-transfer gate: 14 rows, rho_n<exp(-2h_n/5), det(B_n)>=h_n^2/16, 1 exact aggregate countermodel family, 2 open Xi obligations",),
        category="asymptotic theorem-search certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered absolute-phase anchor reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_absolute_phase_anchor_reduction.py",),
        expected=("validated Newman centered absolute-phase anchor reduction: 16 rows, 1 branch-free anchor, 1 determinant collapse, 2 open arithmetic obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered direct-projection regime reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_direct_projection_regime_reduction.py",),
        expected=("validated Newman centered direct-projection regime reduction: 16 rows, 2 exact projection normal forms, 1 heat-direction countermodel, 2 open Xi regimes",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered carrier-kernel Abel-prefix reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_carrier_kernel_abel_prefix_reduction.py",),
        expected=("validated Newman centered carrier-kernel Abel reduction: 19 rows, 3 exact kernel forms, 1 all-fiber prefix scalar, 2 exact countermodels, 2 open Xi obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered normalized-prefix phase-flux reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_normalized_prefix_phase_flux_reduction.py",),
        expected=("validated Newman normalized-prefix phase-flux reduction: 15 rows, 2 certified spiral bounds, 1 exact recrossing guard, 1 dyadic completion diagnostic, 1 O(N) first-jet flux, 2 open Xi obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered prime-power heat-block composition guard",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_block_composition_guard.py",),
        expected=("validated Newman prime-power heat-block guard: 14 rows, 1 absolute-rate collar, 1 exact normalized-block zero-free/shifted-winding theorem, 1 Gaussian-mixture audit, 1 joined p-free decomposition, 1 exact two-block winding-3 countermodel, 2 open Xi obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered joined dyadic odd-prefix first-jet reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_joined_dyadic_odd_prefix_first_jet_reduction.py",),
        expected=("validated Newman joined dyadic odd-prefix first-jet reduction: 17 rows, 1 exact heat-shift factorization, 1 five-current closure, 1 phase-frozen joined zero-free theorem, 2 exact route guards, 2 open Xi obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered phase-cylinder Jacobian transport guard",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_phase_cylinder_jacobian_transport_guard.py",),
        expected=("validated Newman phase-cylinder Jacobian transport guard: 16 rows, 1 exact pair-kernel Jacobian, 1 cylinder winding-transport theorem, 3 exact route guards, 2 open Xi obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered endpoint Schur-Cohn first-jet guard",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_endpoint_schur_cohn_first_jet_guard.py",),
        expected=("validated Newman endpoint Schur-Cohn first-jet guard: 20 rows, 1 exact outside-stability theorem, 1 conditional reflection-product margin, 3 exact route guards, 2 open Xi obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered endpoint first-pivot odd small-ball guard",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_endpoint_first_pivot_odd_small_ball_guard.py",),
        expected=("validated Newman endpoint first-pivot odd small-ball guard: 18 rows, 1 exact denominator-cancelled pivot, 1 linked-logarithmic-phase countermodel, 1 small-endpoint countermodel family, 2 open Xi routes",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered endpoint odd-fibre correlation feasibility guard",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_endpoint_odd_fibre_correlation_feasibility_guard.py",),
        expected=("validated Newman endpoint odd-fibre correlation feasibility guard: 18 rows, 1 exact physical odd-fibre pivot, 1 endpoint-projection normal form, 1 five-current null-family guard, 1 linked-two-jet guard, 1 Schur route downgrade, 2 open routes",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered contact signed-transport reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_contact_signed_transport_reduction.py",),
        expected=("validated Newman contact signed-transport reduction: 20 rows, exact centered endpoint/carrier transport, slope-order Abel identity, conditional margin, Gram-rank guard, endpoint-shaped null guard, five-current first jet, and adjacent-cutoff law",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered interior projective-current gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_interior_projective_current_gate.py",),
        expected=("validated interior projective-current gate: 20 rows, exact radial-current cancellation, proved 5/64 negative-definite interior current, terminal-edge and order-swap nonpromotion guards",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered pairwise projective-alignment gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_pairwise_projective_alignment_gate.py",),
        expected=("validated pairwise projective-alignment gate: 21 rows, exact relative current and pole join, frozen-model Sturm term, moving-scale and residual-order guards, aggregate occupation identity",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered signed occupation transport reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_signed_occupation_transport_reduction.py",),
        expected=("validated signed occupation transport reduction: 25 rows, exact weak/cumulative laws, Xi source collapse, N^sigma normalization, sub-1e-7 residual budget, direct-transversality equivalence, 3 route guards",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered complete prime-power-chain occupation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_complete_prime_power_chain_occupation_gate.py",),
        expected=("validated complete prime-power-chain occupation gate: 21 rows, exact chain factorization and Euler moment, negative-zero-positive threshold guard, external-phase reversal, order-N singleton theorem, p-free telescoping, joined-boundary handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Q208-base/Q209 successor-shell coverage gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_q208_base_q209_shell_coverage_gate.py",),
        expected=("validated Q208-base/Q209 shell coverage gate: 15 rows, 5 coverage regions, 2 open shell regions, 7 oriented boundary arcs, 416 outer old half-cells, 2 new strip half-columns, Q209 not certified",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered Abel-scalar shear-flux reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_abel_scalar_shear_flux_reduction.py",),
        expected=("validated Newman Abel-scalar shear-flux reduction: 18 rows, 1 exact shear/scale homotopy, 1 reduced O(N) physical flux, 1 certified crossing-orientation handoff, 1 backward-heat many-crossing guard, 3 open obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered dominant-ray connector phase cap",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_dominant_ray_connector_phase_cap.py",),
        expected=("validated Newman dominant-ray connector phase cap: 17 rows, 1 uniform cone, 1 strict quarter-turn connector cap, 1 reduced 3*pi/2 horizontal ledger, 1 endpoint-only winding guard, 3 open obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered uniform-q>=1 degree-excision reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_uniform_q_ge_1_degree_excision_reduction.py",),
        expected=("validated Newman uniform-q>=1 degree excision: 16 rows, 1 exact interface chain, 1 conditional outer-degree excision, 1 paired many-turn cancellation guard, 1 boundary-only contact guard, 2 open Xi obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered oscillatory-spliced outer-collar reduction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_oscillatory_spliced_outer_collar_reduction.py",),
        expected=("validated Newman oscillatory-spliced outer collar: 15 rows, 1 exact current frontier, 1 epsilon-raised collar, 1 three-regime cover, 1 conditional degree excision, 2 open Xi obligations",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered reciprocal-saddle self-duality gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_saddle_self_duality_gate.py",),
        expected=("validated Newman reciprocal-saddle self-duality gate: 14 rows, 1 exact stationary map, 1 exact heat-amplitude self-duality, 1 Xi defect bound, 1 critical reciprocal-tail map, 2 nonpromotion guards, 1 open theorem",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered reciprocal normalizer-phase reinforcement guard",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_normalizer_phase_reinforcement_guard.py",),
        expected=("validated Newman reciprocal normalizer-phase reinforcement guard: 13 rows, 2 exact phase identities, 1 uniform conjugate lock, 3 nonpromotion guards, 1 rejected raw-pair route, 1 open signed-kernel target",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 first-order centered signed-handoff reciprocal-symbol guard",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_first_order_centered_signed_handoff_reciprocal_symbol_guard.py",),
        expected=("validated Newman signed-handoff reciprocal-symbol guard: 13 rows, 3 exact symbol identities, 2 first-jet identities, 2 rational exponent audits, 2 endpoint/cutoff guards, 1 retired pairwise route",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman positive-boundary attainment lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_positive_boundary_attainment_lemma.py",),
        expected=("validated Jensen-window PF Newman positive-boundary attainment lemma: 10 rows, 0 issues, 2 published compactness inputs, 1 finite-boundary attainment theorem, 1 positive-time simplicity equivalence, 1 arbitrary-multiplicity Hermite split, 9 exact Hermite checks, 1 cluster-energy blow-up, 1 open Xi endpoint handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman positive-boundary delta-localization gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_positive_boundary_delta_localization_gate.py",),
        expected=("validated Newman positive-boundary delta-localization gate: 10 rows, 0 issues, 5 exact reductions, 1 cofinal-strip criterion, 1 quadratic endpoint-uniformity countermodel, 1 delta-dependent open target",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman positive-boundary diagonal-exhaustion gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.py",),
        expected=("validated Newman positive-boundary diagonal-exhaustion gate: 10 rows, 0 issues, 5 exact exhaustion reductions, 1 compact-shell composition, 1 arbitrary-height classical-field boundary countermodel, 1 independent-rate open target",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman first-jet winding gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_first_jet_winding_gate.py",),
        expected=("validated Newman first-jet winding gate: 10 rows, 0 issues, 9 Hermite audits, 5 exact local/global index reductions, 1 signed winding theorem, 1 half-rectangle flux criterion, 1 open Xi edge target",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta first diagonal-shell interval certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_first_diagonal_shell_interval_certificate.py",),
        expected=("validated Newman theta first diagonal-shell interval certificate: 8 rows, 0 issues, 5 certified edge boxes, 0 subdivisions, 0 unresolved, minimum ratio >6/5, 1 first-stage no-contact theorem, 1 zero-winding composition, 1 open second-stage handoff",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta second diagonal-shell first-block route guard",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_second_diagonal_shell_first_block_route_guard.py",),
        expected=("validated Newman theta second diagonal-shell first-block route guard: 6 rows, 0 issues, 5 exact point audits, 4 strict two-branch failures, maximum failed ratio <24/25, 1 second-stage route guard, 1 open replacement handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta second diagonal-shell two-block interval certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_second_diagonal_shell_two_block_interval_certificate.py",),
        expected=("validated Newman theta second diagonal-shell two-block interval certificate: 10 rows, 0 issues, 25 certified slab boxes, 0 subdivisions, 0 unresolved, 25 negative-value boxes, two analytic n>=3 moment bounds, minimum value ratio >39000, 1 second-stage no-contact theorem, 1 zero-winding composition, 1 open third-stage handoff",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta third diagonal-shell two-block interval certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_third_diagonal_shell_two_block_interval_certificate.py",),
        expected=("validated Newman theta third diagonal-shell two-block interval certificate: 8 rows, 0 issues, 48 certified slab boxes, 0 subdivisions, 0 unresolved, 48 negative-value boxes, minimum value ratio >1000, 1 third-stage no-contact theorem, 1 zero-winding composition, 1 open fourth-stage handoff",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta fourth diagonal-shell two-block interval certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fourth_diagonal_shell_two_block_interval_certificate.py",),
        expected=("validated Newman theta fourth diagonal-shell two-block interval certificate: 8 rows, 0 issues, 120 certified slab boxes, 0 subdivisions, 0 unresolved, 105 negative-value boxes, 15 derivative-only boxes, minimum disjunction ratio >88000, 1 fourth-stage no-contact theorem, 1 zero-winding composition, 1 open fifth-stage handoff",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta two-block shell frontier scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_two_block_shell_frontier_scout.py",),
        expected=("validated Newman theta two-block shell frontier scout: 26 shell rows, 0 issues, j=5..30, 26 certified shells, first finite frontier=none",),
        category="bounded numerical evidence",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta fixed-block cofinal obstruction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fixed_block_cofinal_obstruction_gate.py",),
        expected=("validated Newman theta fixed-block cofinal obstruction gate: 9 rows, 0 issues, 3 exact integration-by-parts bounds, 1 cofinal method obstruction, 1 endpoint-cancellation theorem, 26 finite diagnostic shells, 1 cancellation-aware open handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta modular-blend gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_blend_gate.py",),
        expected=("validated Jensen-window PF Newman theta modular-blend gate: 12 rows, 0 issues, 1 exact positive modular partition, 1 positive-time normal series, 6 transform witnesses, 2 Jensen witnesses",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta modular-blend high-frequency scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_blend_high_frequency_scout.py",),
        expected=("validated Newman modular-blend high-frequency scout: 10 rows, 0 issues, 2 times, 5 frequencies, 3 fixed blocks, 1 cancellation non-promotion gate",),
        category="bounded numerical evidence",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta modular-blend adaptive-saddle gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.py",),
        expected=("validated Newman modular-blend adaptive-saddle gate: 8 rows, 0 issues, 2 exact saddle laws, 1 square-root transition theorem, 10 adaptive diagnostics, 10 matched monotonicity diagnostics, 1 collar guard, 1 retired monotonicity branch",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta adaptive modular C1 remainder contract",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.py",),
        expected=("validated Newman theta adaptive modular C1 remainder contract: 10 rows, 0 issues, 2 exact modular theorems, 3 exact C1 inequalities/compositions, 1 fixed-block guard, 1 adaptive saddle theorem, 10 cancellation diagnostics, 1 open quantitative target",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta modular-tail derivative envelope gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_tail_derivative_envelope_gate.py",),
        expected=("validated Newman theta modular-tail derivative envelope gate: 9 rows, 0 issues, 3 exact inequalities/identities, 3 exact tail-scale theorems/compositions, 1 fixed-collar guard, 1 open separation handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta modular-tail derivative budget scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_tail_derivative_budget_scout.py",),
        expected=("validated Newman theta modular-tail derivative-budget scout: 35 derivative rows, 64 C1 rows, 0 issues, 2 node ladders, 1 arithmetic-cap stress, 1 fixed-collar non-promotion pattern",),
        category="bounded numerical evidence",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta modular-tail Arb quadratic pilot",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_tail_arb_quadratic_pilot.py",),
        expected=("Newman theta modular-tail Arb quadratic pilot checker: ok=True rows=5 issues=0",),
        category="rigorous numerical pilot",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta forward-remainder tail gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_forward_remainder_tail_gate.py",),
        expected=("validated Newman theta forward-remainder tail gate: 6 rows, 0 issues, 10 exact derivative polynomials, 10 explicit arithmetic-tail bounds, 1 open stable-matrix handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta arbitrary-N stable-remainder gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_arbitrary_n_stable_remainder_gate.py",),
        expected=("validated Newman theta arbitrary-N stable remainder gate: 8 rows, 0 issues, 10 derivative polynomials, 100 witness bounds, 1 exact all-N split, 1 cofinal handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta switch-defect explicit-constant gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_switch_defect_explicit_constant_gate.py",),
        expected=("validated Newman theta switch-defect explicit constant gate: 9 rows, 8 all-N formula witnesses, 0 issues, 1 open cofinal retained-separation handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta adaptive gamma-scale tail gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_adaptive_gamma_scale_tail_gate.py",),
        expected=("validated Newman theta adaptive gamma-scale tail gate: 8 rows, 3 boundary witnesses, 0 issues, cofinal x>=245 omitted-tail theorem, 1 open retained-margin handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta forward six-term bridge tail gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_forward_six_term_bridge_tail_gate.py",),
        expected=("validated Newman theta forward six-term bridge tail gate: 8 rows, 3 witnesses, 0 issues, direct gamma-scale tail theorem on 38<=x<=245, 1 open retained-cover handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta forward adaptive square-root tail gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_forward_adaptive_sqrt_tail_gate.py",),
        expected=("validated Newman theta forward adaptive square-root tail gate: 8 rows, 0 issues, tunable cofinal exp(-3h)/50 value/derivative tail with saddle-scale retained count, 1 open retained-separation handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta square-root to corrected RS C1 transfer gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_forward_sqrt_to_corrected_rs_C1_transfer_gate.py",),
        expected=("validated theta square-root to corrected RS C1 transfer gate: 9 rows, 0 issues, exact dual approximation, direct ordinary threshold, and quantified certified-envelope mismatch",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman one-sided phase/moment bridge gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.py",),
        expected=("validated Newman one-sided phase/moment bridge gate: 16 rows, 0 issues, 1 zero-free complex lift, 1 score probability, 1 exact phase-contact system, 1 signed-Hankel moment bridge, 1 score/Beta Abel measure, 1 full-scale interlacing guard, 1 Abel-kernel sign guard, 1 closed score flow, 1 normalizer transfer, 1 conditioning guard, 1 shape countermodel, 3 precision-stable diagnostics, 1 open Xi joint-avoidance target",),
        category="exact theorem-search algebra",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman score-Abel radial dimension-lift gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_score_abel_radial_dimension_lift_gate.py",),
        expected=("validated score-Abel radial dimension-lift gate: 15 rows, 0 issues, 1 imported Abel coordinate, 1 radial probability ladder, 1 planar marginal, 1 Bessel derivative identity, 1 own-dimension positive-definiteness theorem, 1 dimension walk, 1 Jensen identification, 1 Schoenberg scope guard, 3 exact Gaussian countermodel rows, 2 imported Xi low-degree closures, 1 imported degree-361 closure, 1 nonpromotion gate",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta stable-remainder outer-tail gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate.py",),
        expected=("validated Newman theta stable-remainder outer-tail gate: 7 rows, 203 directed entries, 7 retained counts, 0 issues",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta stable-remainder Arb quadratic matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_stable_remainder_arb_quadratic_matrix.py",),
        expected=("validated Newman theta stable-remainder Arb matrix: 223/223 cache rows, 7 complete retained counts, 0 issues",),
        category="rigorous numerical certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta full derivative-budget certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_full_derivative_budget_certificate.py",),
        expected=("validated Newman theta full derivative-budget certificate: 7 rows, 7 full d0/d1 budgets, 0 issues, 1 open cofinal first-jet handoff",),
        category="rigorous numerical certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta modular-retained Q31 interval certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_modular_retained_q31_interval_certificate.py",),
        expected=("validated Newman theta modular-retained Q31 interval certificate: 1240/1240 initial boxes, 1892 Taylor boxes, 0 unresolved, 0 issues",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta forward six-term finite-bridge interval certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_forward_six_term_finite_bridge_interval_certificate.py",),
        expected=("validated Newman theta forward six-term finite-bridge interval certificate: 7 rows, 414/414 panels, 7102 certified leaves, 0 unresolved panels, 0 issues",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta Q31 margin-geometry audit",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_q31_margin_geometry_audit.py",),
        expected=("validated Q31 margin geometry audit: 1240 cache records, 1566 leaves, 7 rows, 0 issues, 1 saddle-normalized cofinal handoff",),
        category="bounded numerical evidence",
    ),
    GateSpec(
        name="Jensen-window PF Newman strict-Laguerre correlation target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_strict_laguerre_correlation_target.py",),
        expected=("validated Jensen-window PF Newman strict-Laguerre correlation target: 10 rows, 0 issues, 1 strict-Laguerre equivalence, 1 exact correlation identity, 1 Wiener-density equivalence, 1 RH-equivalent density target, 1 exact strict-log-concavity/positive-definiteness countermodel, 2 non-promotion gates, 1 open Xi handoff",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta-curvature probability/operator gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_curvature_probability_operator_gate.py",),
        expected=("validated Newman theta-curvature probability/operator gate: 20 rows, 0 issues, 3 theta-primitive identities, 2 monotone-convex inequalities, 1 probability law, 1 fixed-weight theta mixture, 1 uniform C2 dominant-block budget, 1 dominant-mass endpoint-jet guard, 2 transform identities, 1 oscillator factorization, 1 Doob diffusion identity, 1 Sturm nodal-loss guard, 1 asymptotic drift-curvature guard, 1 characteristic contact reduction, 1 componentwise contact decomposition, 1 explicit C1 tail disjunction, 1 endpoint Laplace comparison, 1 generic double-contact guard, 1 open Xi transversality handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta compact-transversality scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_compact_transversality_scout.py",),
        expected=("validated Newman theta compact-transversality scout: 12 rows, 0 issues, 3 exact moment/kernel inequalities, 1 uniform near-origin no-contact theorem, 61951 compact diagnostic points, 0 compact grid failures, 1 independent quadrature crosscheck, 1 outer nonpromotion guard, 2 high-frequency composition theorems, 2 open certification targets",),
        category="bounded numerical evidence",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta compact-transversality interval certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_compact_transversality_interval_certificate.py",),
        expected=("validated Newman theta compact-transversality interval certificate: 8 rows, 0 issues, 4 exact identities/inequalities, 1900 certified boxes, 10 adaptive subdivisions, 0 unresolved boxes, 1 compact no-contact theorem, 1 origin composition, 1 open high-frequency handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman backward-Pick collision-bridge audit",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_backward_pick_collision_bridge_audit.py",),
        expected=("validated Newman backward-Pick collision-bridge audit: 12 rows, 0 issues, 3 exact flow identities, 2 Pick-sign regions, 1 square-root speed blowup, 1 cutoff-hiding theorem, 1 noncommuting-limit obstruction, 1 preprint gap, 1 open Xi repair",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 oscillatory zeta handoff theorem",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.py",),
        expected=("validated Newman Polymath-15 oscillatory zeta handoff theorem: 12 rows, 0 issues, 11 exponent pairs, 10 exact transitions, threshold 4911678521/1933561194, 1 zeta-jet handoff, 1 exact-H asymptotic theorem",),
        category="asymptotic theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 cancellation/zero-free wall gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.py",),
        expected=("validated Newman Polymath-15 cancellation/zero-free wall gate: 10 rows, 0 issues, 11 frontier points, 1 exact c_* maximum, 1 conditional c=2 handoff, 1 inner-wall nonpromotion gate, 1 open Wronskian handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 Gaussian/Legendre duality gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate.py",),
        expected=("validated Newman Polymath-15 Gaussian/Legendre duality gate: 10 rows, 0 issues, 1 exact Gaussian identity, 1 exact Legendre equivalence, 1 c_* equality point, 1 c=2 deficit, 1 nonpromotion gate",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 ANTEDB beta-frontier audit",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit.py",),
        expected=("validated Newman Polymath-15 ANTEDB beta-frontier audit: 10 rows, 0 issues, 6 post-2023 pairs, 1 exact current-hull maximum, 1 direct-beta contact, 1 unchanged c=2 deficit, 1 nonpromotion gate",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 Lambda 0.1965 provenance audit",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_lambda01965_provenance_audit.py",),
        expected=("validated Newman Polymath-15 Lambda 0.1965 provenance audit: 10 rows, 0 issues, 1 published interval, 22 portable passes, 3 compiled-package boundaries, 1 candidate nonpromotion gate",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 critical scaled coercivity target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target.py",),
        expected=("validated Newman Polymath-15 critical scaled coercivity target: 10 rows, 0 issues, 2 exact curvature identities, 1 published endpoint correction, 1 refined remainder, 1 open coercivity target",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman Polymath-15 critical component Wronskian gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_component_wronskian_gate.py",),
        expected=("validated Newman Polymath-15 critical component Wronskian gate: 10 rows, 0 issues, 5 exact identities/reductions, 1 exact ordered-speed countermodel, 4 corrected diagnostics, 1 open arithmetic small-ball target",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman correlation hierarchy Gaussian-mixture gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.py",),
        expected=("validated Jensen-window PF Newman correlation hierarchy Gaussian-mixture gate: 11 rows, 0 issues, 3 exact hierarchy identities, 1 universal boundary-contact signature, 1 Gaussian-mixture sufficient theorem, 2 numerical diagnostics, 1 exact super-Gaussian tail theorem, 2 non-promotion gates, 1 tail-compatible handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman positive-time strong-log-concavity gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_positive_time_strong_logconcavity_gate.py",),
        expected=("validated Jensen-window PF Newman positive-time strong-log-concavity gate: 9 rows, 0 issues, 1 published input, 2 exact curvature/admissibility theorems, 1 Arb threshold certificate, 1 Prekopa correlation theorem, 1 square-transform identity, 1 Xi-specific nonpromotion gate, 1 target-window margin, 1 weighted-hierarchy handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman weighted strong-log-concavity countermodel gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.py",),
        expected=("validated Jensen-window PF Newman weighted strong-log-concavity countermodel gate: 9 rows, 0 issues, 1 exact strong-curvature bound, 1 exact root-variable concavity theorem, 1 explicit theta-tail admissible kernel, 1 Gaussian endpoint identity, 1 exact endpoint witness, 1 Arb theta-tail witness, 1 weighted-correlation countermodel, 1 Xi-specific handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman strict-Laguerre monotonicity scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_strict_laguerre_monotonicity_scout.py",),
        expected=("validated Jensen-window PF Newman strict-Laguerre monotonicity scout: 9 rows, 0 issues, 5 exact target identities, 2 exact classical-route obstructions, 6 dense time rows, 20 high-frequency rows, 1 theta-tail nonpromotion guard, 1 Arb Xi monotonicity rejection",),
        category="non-promotion guards",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman theta-summand spectral-square gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_summand_spectral_square_gate.py",),
        expected=("validated Jensen-window PF Newman theta-summand spectral-square gate: 12 rows, 0 issues, 7 exact transform identities, 1 xi-reconstruction non-promotion gate, 2 exact finite-truncation theorems, 1 numerical sign diagnostic, 1 infinite-cancellation handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta/Bessel higher-shift regularization gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.py",),
        expected=("validated Jensen-window PF Newman theta/Bessel higher-shift regularization gate: 9 rows, 0 issues, 3 exact expansion identities, 1 coefficient sign theorem, 1 fixed-block Bessel theorem, 3 spectral non-promotion gates, 1 coupled modular handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman theta cell-renormalization gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_theta_cell_renormalization_gate.py",),
        expected=("validated Jensen-window PF Newman theta cell-renormalization gate: 10 rows, 0 issues, 4 exact kernel/transform identities, 3 exact convergence/sign theorems, 1 coupled Laguerre identity, 1 positive-time obstruction, 1 modular handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman Gasper fake-Xi remainder gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_gasper_fake_xi_remainder_gate.py",),
        expected=("validated Jensen-window PF Newman Gasper fake-Xi remainder gate: 8 rows, 0 issues, 2 exact transform identities, 1 established real-zero benchmark, 1 scalar algebra theorem, 2 interval scalar witnesses, 2 high-precision cross-checks, 1 exact positive-convolution obstruction, 1 sign-aware handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman Gasper residual two-block gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_gasper_residual_two_block_gate.py",),
        expected=("validated Jensen-window PF Newman Gasper residual two-block gate: 8 rows, 0 issues, 2 exact kernel theorems, 1 exact Laguerre quadratic, 2 Acb derivative certificates, 2 beta intervals covered, 1 exhaustive positive-residual obstruction, 1 multiplier guard, 1 signed handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman classical three-block residual gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_classical_three_block_residual_gate.py",),
        expected=("validated Jensen-window PF Newman classical three-block residual gate: 10 rows, 0 issues, 2 established classical real-zero benchmarks, 2 positive-kernel residual theorems, 1 compact parameter reduction, 1 exact bivariate Laguerre identity, 3 Acb spectral certificates, 64908 parameter boxes covered, 2 classical residual obstructions, 1 Gasper square-scope guard, 1 coupled handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman signed universal-factor residual gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_signed_universal_factor_residual_gate.py",),
        expected=("validated Jensen-window PF Newman signed universal-factor residual gate: 8 rows, 0 issues, 1 exact multiplier reduction, 1 rational parameter rectangle, 2 discriminant exclusions, 2 Acb spectral certificates, 4094 adaptive leaves, 3416 base boxes, maximum depth 6, 1 exhaustive signed universal-factor obstruction, 1 coupled handoff",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Laguerre scale-mixture gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_laguerre_scale_mixture_gate.py",),
        expected=("validated Jensen-window PF Laguerre scale-mixture gate: 10 rows, 0 issues, 3 exact kernel identities, 1 individual-kernel hyperbolicity theorem, 1 positive-mixture countermodel, 1 Gamma reduction, 1 half-integer all-degree theorem, 1 log-concavity countermodel, 1 open Xi-specific handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF rank-two boundary-family lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_rank_two_boundary_family_lemma.py",),
        expected=("validated Jensen-window PF rank-two boundary-family lemma: 11 rows, 0 issues, 4 exact identities, 1 all-degree factorization, 1 integer-product closure, 2 exact countermodels, 1 open structural handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF multiplier counting-measure target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_multiplier_counting_measure_target.py",),
        expected=("validated Jensen-window PF multiplier counting-measure target: 10 rows, 0 issues, 4 exact rows, 1 finite evidence row, 3 countermodel rows, 1 live route, 0 ready-to-apply rows",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Mellin multiplier power-sum obstruction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_mellin_multiplier_power_sum_obstruction.py",),
        expected=("validated Jensen-window PF Mellin multiplier power-sum obstruction: 11 log moments, 9 power sums, 6 Hankel determinants, 3 negative, 0 inconclusive, 1 continuous route ruled out, 0 discrete routes ruled out, 0 issues",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda first-summand saddle-wall closure",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.py",),
        expected=("validated Jensen-window PF negative-lambda first-summand saddle-wall closure: 9 rows, 0 issues, 9 positive samples, 9 quarter-k2 samples, 9 bracketed saddles, 0 open requirements, 2 ready-to-apply rows",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda first-summand cumulant bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.py",),
        expected=("validated Jensen-window PF negative-lambda first-summand cumulant bridge: 8 rows, 0 issues, 4 exact identities, 1 conditional bridge, 9 positive samples, 0 open requirements, 3 ready-to-apply rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda first-summand leading-saddle certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda first-summand leading-saddle certificate: 8 rows, 0 issues, 40740 positive leading intervals, 40740 positive cubic-correction intervals, 40740 positive fifth-correction intervals, 3 positive analytic ray gates, 9 positive seventh-remainder samples, 1 open remainder, 0 ready-to-apply rows",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda first-summand paired-remainder certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda first-summand paired-remainder certificate: 8 rows, 0 issues, 40736 eighth-envelope intervals, 4074 compact remainder blocks, 1 open asymptotic ray, 0 ready-to-apply rows",),
        category="interval theorem certificates",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda first-summand paired-remainder ray certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda first-summand paired-remainder ray certificate: 9 rows, 0 issues, 1 analytic ray theorem, 1 global remainder closure, 0 open rays, 2 ready-to-apply rows",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF heat-flow ratio cone invariance lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.py",),
        expected=("validated Jensen-window PF heat-flow ratio cone invariance lemma: 6 exact rows, 315 lower rows, 315 upper rows, 310 monotone rows, 0 issues",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF heat-flow cone-entry asymptotic target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_heat_flow_cone_entry_asymptotic_target.py",),
        expected=("validated Jensen-window PF heat-flow cone-entry asymptotic target: 8 rows, 0 issues, 1 ready-to-apply rows, 0 live routes",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Phi Taylor cone-entry sign scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_phi_taylor_cone_entry_sign_scout.py",),
        expected=("validated Jensen-window PF Phi Taylor cone-entry sign scout: 4 coefficient balls, 2 certified signs, 0 ready-to-apply rows, 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda cone-entry prefix scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py",),
        expected=("validated Jensen-window PF negative-lambda cone-entry prefix scout: 69 coefficient rows, 63 lower-wall rows, 63 upper-wall rows, 60 monotone-gap rows, 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda cone-entry prefix k30 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_cone_entry_prefix_k30_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k30_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda cone-entry prefix scout: 93 coefficient rows, 87 lower-wall rows, 87 upper-wall rows, 84 monotone-gap rows, 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda cone-entry prefix k50 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_cone_entry_prefix_k50_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k50_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda cone-entry prefix scout: 153 coefficient rows, 147 lower-wall rows, 147 upper-wall rows, 144 monotone-gap rows, 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda cone-entry prefix k60 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_cone_entry_prefix_k60_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k60_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda cone-entry prefix scout: 183 coefficient rows, 177 lower-wall rows, 177 upper-wall rows, 174 monotone-gap rows, 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda cone-entry prefix k80 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_cone_entry_prefix_k80_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k80_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda cone-entry prefix scout: 243 coefficient rows, 237 lower-wall rows, 237 upper-wall rows, 234 monotone-gap rows, 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda cone-entry prefix k100 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_cone_entry_prefix_k100_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k100_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda cone-entry prefix scout: 303 coefficient rows, 297 lower-wall rows, 297 upper-wall rows, 294 monotone-gap rows, 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda cone-entry prefix k150 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_cone_entry_prefix_k150_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k150_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda cone-entry prefix scout: 453 coefficient rows, 447 lower-wall rows, 447 upper-wall rows, 444 monotone-gap rows, 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda cone-entry prefix k200 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_cone_entry_prefix_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_cone_entry_prefix_k200_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k200_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda cone-entry prefix scout: 603 coefficient rows, 597 lower-wall rows, 597 upper-wall rows, 594 monotone-gap rows, 0 issues",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda finite-collar contract",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py",),
        expected=("validated Jensen-window PF negative-lambda finite-collar contract: active depth K=19, 57 active lower rows, 57 active upper rows, 57 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda finite-collar k30 contract",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_k30_contract.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_finite_collar_k30_contract.md",
        ),
        expected=("validated Jensen-window PF negative-lambda finite-collar contract: active depth K=27, 81 active lower rows, 81 active upper rows, 81 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda finite-collar k50 contract",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_k50_contract.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_finite_collar_k50_contract.md",
        ),
        expected=("validated Jensen-window PF negative-lambda finite-collar contract: active depth K=47, 141 active lower rows, 141 active upper rows, 141 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda finite-collar k60 contract",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_k60_contract.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_finite_collar_k60_contract.md",
        ),
        expected=("validated Jensen-window PF negative-lambda finite-collar contract: active depth K=57, 171 active lower rows, 171 active upper rows, 171 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda finite-collar k80 contract",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_k80_contract.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_finite_collar_k80_contract.md",
        ),
        expected=("validated Jensen-window PF negative-lambda finite-collar contract: active depth K=77, 231 active lower rows, 231 active upper rows, 231 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda finite-collar k100 contract",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_k100_contract.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_finite_collar_k100_contract.md",
        ),
        expected=("validated Jensen-window PF negative-lambda finite-collar contract: active depth K=97, 291 active lower rows, 291 active upper rows, 291 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda finite-collar k150 contract",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_k150_contract.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_finite_collar_k150_contract.md",
        ),
        expected=("validated Jensen-window PF negative-lambda finite-collar contract: active depth K=147, 441 active lower rows, 441 active upper rows, 441 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda finite-collar k200 contract",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_finite_collar_contract.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_finite_collar_k200_contract.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_finite_collar_k200_contract.md",
        ),
        expected=("validated Jensen-window PF negative-lambda finite-collar contract: active depth K=197, 591 active lower rows, 591 active upper rows, 591 active monotone rows, 6 collar lower/upper rows, 3 collar monotone rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda tail-barrier scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py",),
        expected=("validated Jensen-window PF negative-lambda tail-barrier scout: 63 cone-buffer rows, 60 defect-monotone rows, 63 one-third-buffer rows, 60 scaled-defect increase rows, 1 rejected candidate, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda tail-barrier k30 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_k30_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_tail_barrier_k30_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda tail-barrier scout: 87 cone-buffer rows, 84 defect-monotone rows, 87 one-third-buffer rows, 84 scaled-defect increase rows, 1 rejected candidate, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda tail-barrier k50 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_k50_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_tail_barrier_k50_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda tail-barrier scout: 147 cone-buffer rows, 144 defect-monotone rows, 128 one-third-buffer rows, 144 scaled-defect increase rows, 1 rejected candidate, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda tail-barrier k60 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_k60_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_tail_barrier_k60_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda tail-barrier scout: 177 cone-buffer rows, 174 defect-monotone rows, 139 one-third-buffer rows, 174 scaled-defect increase rows, 1 rejected candidate, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda tail-barrier k80 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_k80_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_tail_barrier_k80_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda tail-barrier scout: 237 cone-buffer rows, 234 defect-monotone rows, 159 one-third-buffer rows, 234 scaled-defect increase rows, 1 rejected candidate, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda tail-barrier k100 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_k100_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_tail_barrier_k100_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda tail-barrier scout: 297 cone-buffer rows, 294 defect-monotone rows, 179 one-third-buffer rows, 294 scaled-defect increase rows, 1 rejected candidate, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda tail-barrier k150 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_k150_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_tail_barrier_k150_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda tail-barrier scout: 447 cone-buffer rows, 444 defect-monotone rows, 179 one-third-buffer rows, 444 scaled-defect increase rows, 1 rejected candidate, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda tail-barrier k200 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_k200_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_tail_barrier_k200_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda tail-barrier scout: 597 cone-buffer rows, 594 defect-monotone rows, 179 one-third-buffer rows, 594 scaled-defect increase rows, 1 rejected candidate, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda scaled-defect frontier k50 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k50_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k50_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda scaled-defect frontier scout: 147 scaled rows, 147 cone rows, 147 half-width rows, 128 one-third rows, 19 one-third failures, 144 scaled-increase rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda scaled-defect frontier k60 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k60_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k60_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda scaled-defect frontier scout: 177 scaled rows, 177 cone rows, 177 half-width rows, 139 one-third rows, 38 one-third failures, 174 scaled-increase rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda scaled-defect frontier k80 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k80_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k80_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda scaled-defect frontier scout: 237 scaled rows, 237 cone rows, 237 half-width rows, 159 one-third rows, 78 one-third failures, 234 scaled-increase rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda scaled-defect frontier k100 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k100_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k100_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda scaled-defect frontier scout: 297 scaled rows, 297 cone rows, 297 half-width rows, 179 one-third rows, 118 one-third failures, 294 scaled-increase rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda scaled-defect frontier k150 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda scaled-defect frontier scout: 447 scaled rows, 447 cone rows, 430 half-width rows, 179 one-third rows, 268 one-third failures, 444 scaled-increase rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda scaled-defect frontier k200 scout",
        command=(
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py",
            "--target",
            "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.json",
            "--note",
            "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.md",
        ),
        expected=("validated Jensen-window PF negative-lambda scaled-defect frontier scout: 597 scaled rows, 597 cone rows, 521 half-width rows, 179 one-third rows, 418 one-third failures, 594 scaled-increase rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda defect-recurrence scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_recurrence_scout.py",),
        expected=("validated Jensen-window PF negative-lambda defect-recurrence scout: 63 buffered rows, 60 defect-monotone rows, 60 width-recurrence rejections, 1 live sufficient routes, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda log-curvature bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_log_curvature_bridge.py",),
        expected=("validated Jensen-window PF negative-lambda log-curvature bridge: 63 simple log-buffer rows, 63 exact defect-buffer rows, 60 curvature-monotone rows, 5 bridge rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda bounded log-curvature target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_bounded_log_curvature_target.py",),
        expected=("validated Jensen-window PF negative-lambda bounded log-curvature target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 live routes, 63 raw-threshold rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda bounded log-curvature k300 obstruction",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.py",),
        expected=("validated Jensen-window PF negative-lambda bounded log-curvature k300 obstruction: 7 rows, 0 issues, 718 two-thirds failures, 894 scaled-curvature increase rows, 0 ready-to-apply rows",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda Gaussian curvature matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_gaussian_curvature_matrix.py",),
        expected=("validated Jensen-window PF negative-lambda Gaussian curvature matrix: 7 matrix rows, 63 positive-deficit rows, 63 bounded-deficit rows, 63 raw-threshold rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda signed Gaussian perturbation matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.py",),
        expected=("validated Jensen-window PF negative-lambda signed Gaussian perturbation matrix: 8 matrix rows, 2 certified Taylor signs, 1 fixed-k activation estimates, 0 ready-to-apply rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda uniform remainder target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_uniform_remainder_target.py",),
        expected=("validated Jensen-window PF negative-lambda uniform remainder target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 open requirements, 3 leading-scale rows",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda Taylor moment budget",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_taylor_moment_budget.py",),
        expected=("validated Jensen-window PF negative-lambda Taylor moment budget: 9 budget rows, 7 tail-start samples, 4 invalid truncation rows, 2 bounded truncation rows, 0 ready-to-apply rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda high-order Taylor scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_high_order_taylor_scout.py",),
        expected=("validated Jensen-window PF negative-lambda high-order Taylor scout: 8 coefficient rows, 35 truncation rows, 9 invalid normalizers, 2 upper-wall violations, 3 overbound rows, 0 ready-to-apply rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda defect-tail theorem target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_tail_theorem_target.py",),
        expected=("validated Jensen-window PF negative-lambda defect-tail theorem target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 live routes",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda half-width tail target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_half_width_tail_target.py",),
        expected=("validated Jensen-window PF negative-lambda half-width tail target: 9 rows, 0 issues, 0 ready-to-apply rows, 0 live routes, 430 half-width rows, 17 half-width failures",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda adaptive scaled-defect target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.py",),
        expected=("validated Jensen-window PF negative-lambda adaptive scaled-defect target: 8 rows, 0 issues, 2 live routes, 597 exact-cone rows, 76 half-width failures",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda adaptive envelope matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_matrix.py",),
        expected=("validated Jensen-window PF negative-lambda adaptive envelope matrix: 7 matrix rows, 0 issues, 594 k-increase rows, 398 lambda-order rows, 76 half-width failures",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda adaptive envelope obligations",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_obligations.py",),
        expected=("validated Jensen-window PF negative-lambda adaptive envelope obligations: 9 obligation rows, 0 issues, 3 exact rows, 3 open requirements, 1 rejected routes",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda raw-moment bridge matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.py",),
        expected=("validated Jensen-window PF negative-lambda raw-moment bridge matrix: 8 matrix rows, 0 issues, 597 raw-cone rows, 594 corridor rows, 76 half-width failures",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda raw-ratio decrement-corridor scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_ratio_decrement_corridor_scout.py",),
        expected=("validated Jensen-window PF negative-lambda raw-ratio decrement-corridor scout: 9 rows, 0 issues, 594 decrement-corridor rows, 591 theta-k-monotone rows, 2 exact counterexamples, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda k300 precision-repair audit",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_k300_precision_repair_audit.py",),
        expected=("validated Jensen-window PF negative-lambda k300 precision-repair audit: 7 rows, 0 issues, 894 repaired decrement-corridor rows, 891 repaired theta-k-monotone rows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda raw-log decrement bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_log_decrement_bridge.py",),
        expected=("validated Jensen-window PF negative-lambda raw-log decrement bridge: 8 rows, 0 issues, 894 log-corridor rows, 894 log-decrease rows, 2 exact counterexamples, 0 ready-to-apply rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda coefficient-curvature corridor bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.py",),
        expected=("validated Jensen-window PF negative-lambda coefficient-curvature corridor bridge: 9 rows, 0 issues, 894 curvature-corridor rows, 894 monotone-curvature rows, 2 exact counterexamples, 0 ready-to-apply rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda linear curvature-barrier scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.py",),
        expected=("validated Jensen-window PF negative-lambda linear curvature-barrier scout: 8 rows, 0 issues, 894 linear-barrier rows, 894 monotone-curvature rows, 2 exact counterexamples, 0 ready-to-apply rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda scaled-curvature monotonicity target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.py",),
        expected=("validated Jensen-window PF negative-lambda scaled-curvature monotonicity target: 10 rows, 0 issues, 2 live routes, 894 scaled-curvature increase rows, 0 ready-to-apply rows",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda scaled-curvature continuous bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.py",),
        expected=("validated Jensen-window PF negative-lambda scaled-curvature continuous bridge: 10 rows, 0 issues, 16074 compact blocks, 318 prefix gaps, 1 analytic ray, 1 all-k scaled-curvature theorem, 0 open requirements",),
        category="interval theorem certificates",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda scaled-curvature log-ceiling bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.py",),
        expected=("validated Jensen-window PF negative-lambda scaled-curvature log-ceiling bridge: 8 rows, 0 issues, 894 scaled-ceiling rows, 894 scaled-log-corridor rows, 894 ceiling-dominance rows, 0 ready-to-apply rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian curvature bridge",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian curvature bridge: 8 rows, 0 issues, 897 B-positive rows, 894 B-decrease rows, 894 C-increase rows, 598 C-lambda-order rows, 0 ready-to-apply rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian Taylor stencil scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian Taylor stencil scout: 8 rows, 0 issues, 3 certified leading-sign rows, 35 truncation rows, 4 all-positive stencil rows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian stencil remainder obligations",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian stencil remainder obligations: 9 rows, 0 issues, 4 positive baseline rows, 31 blocked baseline rows, 4 exact stencil rows, 0 ready-to-apply rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian pointwise tail budget",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian pointwise tail budget: 9 rows, 0 issues, 4 positive baseline rows, 4 budget rows, 0 ready-to-apply rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian next-increment stencil stress",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian next-increment stencil stress: 8 rows, 0 issues, 2 tested next-increment rows, 2 pointwise budget failures, 2 stencil-sign-preserving rows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian degree-16 stencil continuation",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian degree-16 stencil continuation: 7 rows, 0 issues, 4 tested continuation rows, 3 stencil-sign-preserving rows, 1 stencil-sign-failure rows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian degree-16 collar scan",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian degree-16 collar scan: 7 rows, 0 issues, 1301 scan rows, 1045 continuation-positive rows, 718 half-safety rows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian degree-16 real-T collar scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian degree-16 real-T collar scout: 8 rows, 0 issues, 4 positive normalizer rows, 3 certified surrogate stencil rows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian degree-16 Arb real-T collar certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian degree-16 Arb real-T collar certificate: 8 rows, 0 issues, 4 positive normalizer rows, 3 certified stencil rows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian degree-40 Arb collar ladder stress",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian degree-40 Arb collar ladder stress: 8 rows, 0 issues, 13 degree levels, max degree 40, 0 failed Bernstein rows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian degree-40 residual tail budget",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian degree-40 residual tail budget: 8 rows, 0 issues, 5 budget inequalities, 4 finite tail profile rows, 0 ready-to-apply rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian formal-tail obstruction scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian formal-tail obstruction scout: 8 rows, 0 issues, 4 profile rows, 4 formal-tail turnaround rows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian asymptotic remainder target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian asymptotic remainder target: 6 rows, 0 issues, 4 first-omitted rows, 4 optimized-window rows, 0 ready-to-apply rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian actual endpoint remainder scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian actual endpoint remainder scout: 6 rows, 0 issues, 4 endpoint rows, 5 quadrature orders, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian cancellation-reduced remainder grid scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian cancellation-reduced remainder grid scout: 6 rows, 0 issues, 20 grid rows, 5 T values, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian intervalization target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian intervalization target: 6 rows, 0 issues, 8 obligations, 5 open requirements, 0 ready-to-apply rows",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian Phi tail bound scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian Phi tail bound scout: 6 rows, 0 issues, 3 tail bounds below 1e-1000, 2 conditional requirements, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian node-c0 range certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian node-c0 range certificate: 5 rows, 0 issues, 16 Laguerre bound rows, 2 certified side conditions, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian Phi-tail grid certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian Phi-tail grid certificate: 6 rows, 0 issues, 3 certified tail sources, 2 certified side conditions, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian quadrature ladder scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian quadrature ladder scout: 5 rows, 0 issues, 7 ladder rows, 320 reference order, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian quadrature-remainder route matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian quadrature-remainder route matrix: 7 rows, 0 issues, derivative order 640, 2 derivative-sup caps, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian worst-row far-tail split certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian worst-row far-tail split certificate: 7 rows, 0 issues, split y=200, 2 tail ratios below quadrature cap, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian worst-row compact-interval integration scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian worst-row compact-interval integration scout: 7 rows, 0 issues, 6 panels, plain interval Riemann rejected, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian worst-row Chebyshev panel-moment scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian worst-row Chebyshev panel-moment scout: 7 rows, 0 issues, 5 degrees, 4 Cauchy pairs, 3 cap-safe pairs, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian worst-row Arb Chebyshev interpolant-moment scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian worst-row Arb Chebyshev interpolant-moment scout: 7 rows, 0 issues, 4 degrees, 3 Cauchy pairs, 3 cap-safe pairs, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian worst-row interpolation-remainder route matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian worst-row interpolation-remainder route matrix: 8 rows, 0 issues, 6 panel masses, 20 Bernstein budgets, 16 minimal-degree rows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian endpoint parity-repair matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian endpoint parity-repair matrix: 7 rows, 0 issues, 8 odd Taylor rows, order 15, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian endpoint x-panel route matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian endpoint x-panel route matrix: 7 rows, 0 issues, x interval 0<=x<=0.01, 18 Bernstein budgets, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian endpoint x-moment Taylor certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian endpoint x-moment Taylor certificate: 7 rows, 0 issues, 65 exact-moment rows, 1 certified first panel, 5 open compact panels, 0 ready-to-apply rows",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian worst-row compact x-moment Taylor certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian worst-row compact x-moment Taylor certificate: 7 rows, 0 issues, 129 exact-moment rows, 1 certified compact interval, 0 open compact panels, 0 ready-to-apply rows",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian worst-row full-expectation certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian worst-row full-expectation certificate: 8 rows, 0 issues, 3 composed sources, 2 below-one full ratios, 0 open worst-row integral sources, 0 ready-to-apply rows",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian all-row direct expectation certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian all-row direct expectation certificate: 8 rows, 0 issues, 20 grid rows, 20 negative values, 20 negative derivatives, 0 open recorded-grid integral sources, 0 ready-to-apply rows",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian recorded-grid stencil composition certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian recorded-grid stencil composition certificate: 8 rows, 0 issues, 20 residual rows, 4 positive perturbation margins, 5 certified T systems, 0 open recorded-grid stencil sources, 0 ready-to-apply rows",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian finite-collar-segment stencil certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian finite-collar-segment stencil certificate: 10 rows, 0 issues, 4 uniform residual rows, 2 T regimes, 4 positive perturbation margins, 0 open finite-segment stencil sources, 0 ready-to-apply rows",),
        category="promoted interval evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian full-kernel evenness/Cauchy lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian full-kernel evenness/Cauchy lemma: 9 rows, 0 issues, 3 symbolic identities, order>=42 residual zero, 0 ready-to-apply rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian full-real-T fixed-k stencil certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian full-real-T fixed-k stencil certificate: 10 rows, 0 issues, 4 ray residual rows, 4 full-kernel n-tail channels, 4 positive perturbation margins, 0 open full-T fixed-k stencil sources, 0 ready-to-apply rows",),
        category="promoted interval evidence",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian first-omitted denominator certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian first-omitted denominator certificate: 6 rows, 0 issues, 20 denominator rows, 2 ratio-cap rows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian coefficient-core propagation certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian coefficient-core propagation certificate: 6 rows, 0 issues, 22 coefficient rows, 20 propagation rows, 2 intervalization rows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian worst-row Laguerre root-bracket certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian worst-row Laguerre root-bracket certificate: 6 rows, 0 issues, 320 root brackets, 30 zero floating weights, 2 intervalization rows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight midpoint scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight midpoint scout: 6 rows, 0 issues, 320 midpoint weights, 30 repaired floating underflows, 320 direct interval obstructions, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight interval certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight interval certificate: 6 rows, 0 issues, 320 interval weights, 0 Taylor denominator obstructions, 30 repaired floating underflows, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian worst-row finite-part weighted-sum interval certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian worst-row finite-part weighted-sum interval certificate: 6 rows, 0 issues, 320 refined nodes, 320 interval weights, 2 below-one ratios, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda relative-Gaussian worst-row finite-plus-tail budget certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda relative-Gaussian worst-row finite-plus-tail budget certificate: 6 rows, 0 issues, 2 composed ratios, 3 tail sources, 0 ready-to-apply rows",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda raw-moment obstruction matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.py",),
        expected=("validated Jensen-window PF negative-lambda raw-moment obstruction matrix: 7 matrix rows, 0 issues, 3 exact counterexamples, 0 ready-to-apply rows",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda zeta-specific raw-corridor target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.py",),
        expected=("validated Jensen-window PF negative-lambda zeta-specific raw-corridor target: 9 rows, 0 issues, 2 live routes, 2 rejected shortcuts, 0 ready-to-apply rows",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda -100 raw-corridor certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_raw_corridor_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda -100 raw-corridor certificate: 6 rows, 0 issues, 2 theorem inputs, 1 raw-corridor theorem, 0 open requirements",),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF negative-lambda -100 adaptive-defect certificate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.py",),
        expected=("validated Jensen-window PF negative-lambda -100 adaptive-defect certificate: 8 rows, 0 issues, 2 theorem inputs, 4 defect conclusions, 0 open requirements",),
        category="exact theorem composition",
    ),
    GateSpec(
        name="Jensen-window PF monotone contraction stress",
        command=("work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_stress.py",),
        expected=("validated Jensen-window PF monotone contraction stress: 2875 rows, 2875 positive rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF state-space sign-lift obstruction scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_state_space_sign_lift_obstruction_scout.py",),
        expected=("validated Jensen-window PF state-space sign-lift obstruction scout: 3 symbolic rows, 735 mu2 sign-lift obstruction rows, 0 issues",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF Cauchy-Binet low-degree scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py",),
        expected=("validated Jensen-window PF Cauchy-Binet low-degree scout: 15 formula rows, 0 issues, 0 kernel identities found",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF log-concavity frontier scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_log_concavity_frontier_scout.py",),
        expected=("validated Jensen-window PF log-concavity frontier scout: 14 contiguous rows, 0 issues",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF ratio-condition scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_ratio_condition_scout.py",),
        expected=("validated Jensen-window PF ratio-condition scout: 7 candidate rows, 0 issues, 4 rejected by countermodel, 1 rejected by construction",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF contraction-log-concavity scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_contraction_log_concavity_scout.py",),
        expected=("validated Jensen-window PF contraction-log-concavity scout: 1 rejected by construction, 0 issues, 2 negative frontier rows",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="sign-regularity theorem fit matrix",
        command=("work/rh_compute/scripts/check_sign_regularity_theorem_fit_matrix.py",),
        expected=("validated sign-regularity theorem fit matrix with 0 issues",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="positive Schur-specialization target note",
        command=("work/rh_compute/scripts/check_positive_schur_specialization_target.py",),
        expected=("validated positive Schur-specialization target note with 0 issues",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Edrei-log sign diagnostics",
        command=("work/rh_compute/scripts/check_edrei_log_sign_manifest.py",),
        expected=("validated 320 finite Edrei-log sign diagnostics",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Edrei power-Hankel diagnostics",
        command=("work/rh_compute/scripts/check_edrei_power_hankel_manifest.py",),
        expected=("validated 4205 finite Edrei power-Hankel diagnostics",),
        category="promoted finite evidence",
        slow=True,
    ),
    GateSpec(
        name="Edrei midpoint frontier non-promotion guard",
        command=("work/rh_compute/scripts/check_edrei_power_hankel_frontier_manifest.py",),
        expected=("validated 5 non-rigorous Edrei midpoint frontier scouts",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Edrei power-Hankel boundary repair manifest",
        command=("work/rh_compute/scripts/check_edrei_power_hankel_boundary_manifest.py",),
        expected=("validated 2 retired inconclusive blocker rows and 3 repaired positive boundary rows",),
        category="promoted finite evidence",
    ),
    GateSpec(
        name="Edrei moment-recurrence scout manifest",
        command=("work/rh_compute/scripts/check_edrei_quadrature_scout_manifest.py",),
        expected=("validated 1 positive Arb recurrence scout and 1 inconclusive frontier scout",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF Newman outer-frequency frontier consolidation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_outer_frequency_frontier_consolidation_gate.py",),
        expected=("validated Newman outer-frequency frontier consolidation gate: 8 rows, 0 issues, 1900 compact boxes, Q208 base, Q209 first unresolved, 2 open shell regions, no compact rescout",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 remainder supersession frontier gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_remainder_supersession_frontier_gate.py",),
        expected=("validated Newman C1 remainder supersession frontier gate: 10 rows, 0 issues, eta_0=100000, eta_1=200000, 0 C2 boundary requirements, 1 live carrier-plus-near target",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 carrier-near odd-harmonic anchor gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_carrier_near_odd_harmonic_anchor_gate.py",),
        expected=("validated Newman C1 carrier-near odd-harmonic anchor gate: 18 rows, 0 issues, 96 half-cell identities, q=1 D_N>L/4, 2 explicit raw anchors, 1 phase-rotation guard, 1 live stationary residual",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 reciprocal-stationary disk geometry gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_reciprocal_stationary_disk_geometry_gate.py",),
        expected=("validated Newman C1 reciprocal-stationary disk geometry gate: 18 rows, 0 issues, 96 disk samples, 5849 block counts, 2 combined residuals, 2 absolute-budget obstructions, 1 live signed stationary target",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 terminal-phase winding and covariant two-channel gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_terminal_phase_winding_covariant_two_channel_gate.py",),
        expected=("validated Newman C1 terminal-phase winding/covariant two-channel gate: 16 rows, 0 issues, 63 circle-coverage samples, 288 matrix samples, 164 criterion samples, 1 full-cell winding theorem, 1 phase-uniform criterion, 1 live covariant current target",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source saddle-phase lock and centre-symmetry gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_saddle_phase_lock_center_symmetry_gate.py",),
        expected=("validated Newman C1 actual-source saddle-phase lock/centre-symmetry gate: 19 rows, 0 issues, 9 decrement samples, 4 rational audits, 32 centre samples, 5 mismatch samples, 1 reverse full-turn theorem, 1 live joint-phase target",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source raw-centre sign-obstruction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_raw_center_sign_obstruction_gate.py",),
        expected=("validated Newman C1 actual-source raw-centre sign-obstruction gate: 13 rows, 0 issues, 9 scale samples, 6 sign-corner audits, 2 strict sign witnesses, 1 actual-source centre obstruction, 1 macroscopic compensation condition, 1 live complete-current target",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source complete-compensation ledger gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_complete_compensation_ledger_gate.py",),
        expected=("validated Newman C1 actual-source complete-compensation ledger gate: 17 rows, 0 issues, 8 disjoint scalar packages, 147 tone audits, 7 budget audits, 12 nonpromotion audits, 1 sharp 98/100 adverse threshold, 0 common-unit secondary bounds, 1 live signed finite-package target",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source terminal-secondary budget gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_terminal_secondary_budget_gate.py",),
        expected=("validated Newman C1 actual-source terminal-secondary budget gate: 12 rows, 0 issues, 6 exponent audits, 9 rational audits, rho<(L+9)/6, 2 terminal packages closed, joint ratio 1/17160, 5 unbounded secondary packages, 1 live common-unit secondary target",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source full-phase remote budget gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_full_phase_remote_budget_gate.py",),
        expected=("validated Newman C1 actual-source full-phase remote budget gate: 15 rows, 0 issues, 2 symbolic audits, 6 rational audits, 4 constant audits, 1 paired denominator cancellation, 2 ideal remote channels closed, |E_R|<rho*A_T/6, 4 unbounded secondary packages, 1 live common-unit four-package target",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source physical-correction budget gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_physical_correction_budget_gate.py",),
        expected=("validated Newman C1 actual-source physical-correction budget gate: 16 rows, 0 issues, 4 symbolic audits, 9 rational audits, 6 kernel audits, 2 terminal coefficient channels closed, 1 contiguous-kernel L1 theorem, |E_Delta|<rho*A_T/2000, 3 unbounded secondary packages, 1 live common-unit three-package target",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source edge-affine ownership gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_edge_affine_ownership_gate.py",),
        expected=("validated Newman C1 actual-source edge-affine ownership gate: 15 rows, 0 issues, 4 matrix audits, 3 polynomial audits, 2 source audits, 6 rational audits, 1 edge-only ownership correction, 0 absolute affine budgets, 1 signed-main reclassification, 2 remaining absolute secondary packages, 1 live finite-affine target",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source edge-affine remote budget gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_edge_affine_remote_budget_gate.py",),
        expected=("validated Newman C1 actual-source edge-affine remote budget gate: 16 rows, 0 issues, 4 symbolic audits, 16 rational audits, 7 source audits, B_0T_0/rho<61/100, |E_aff^C|<4rho*A_T/25, 2 remaining absolute secondary packages, 1 live signed near-affine target",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source moving-tail budget gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_moving_tail_budget_gate.py",),
        expected=("validated Newman C1 actual-source moving-tail budget gate: 16 rows, 0 issues, 5 symbolic audits, 12 rational audits, 8 source audits, |E_move|<rho*A_T/1000000, 1 remaining absolute secondary package, 1 live signed near-affine target",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source current-transfer quadratic-primitive gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_current_transfer_quadratic_primitive_gate.py",),
        expected=("validated Newman C1 current/transfer quadratic-primitive gate: 17 rows, 0 issues, 6 symbolic audits, 6 rational audits, 6 source audits, pointwise |E_move^cur|<rho*A_T/100, 1 determinant primitive, 1 signed-main reclassification, 1 open theorem row",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source quadratic-residual two-carrier kernel gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_quadratic_residual_two_carrier_kernel_gate.py",),
        expected=("validated Newman C1 quadratic residual two-carrier kernel gate: 16 rows, 0 issues, 8 symbolic audits, 4 finite-functional audits, 4 source audits, 2 exact kernel families, 1 open signed theorem row",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source nonterminal determinant-completion gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_nonterminal_determinant_completion_gate.py",),
        expected=("validated Newman C1 nonterminal determinant-completion gate: 17 rows, 0 issues, 8 symbolic audits, 5 exact-rational audits, 7 source audits, 1 missing linear package identified, 1 corrected determinant completion, 1 open arithmetic row",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source complete outer-determinant recombination gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_complete_outer_determinant_recombination_gate.py",),
        expected=("validated Newman C1 complete outer-determinant recombination gate: 18 rows, 0 issues, 9 symbolic audits, 6 exact-rational audits, 11 source audits, 8 packages collapsed, 2 pointwise packages remain, 1 finite-cell determinant form, 1 open arithmetic row",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 actual-source finite determinant reciprocal-block symmetry gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_actual_source_finite_determinant_reciprocal_block_symmetry_gate.py",),
        expected=("validated Newman C1 finite determinant reciprocal-block symmetry gate: 20 rows, 0 issues, 10 symbolic audits, 6 exact-rational audits, 3 source audits, 3 reciprocal block classes, 1 determinant symmetry, 1 singular-base guard, 1 open arithmetic row",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 calibrated roster-interior phase-root scout",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_calibrated_roster_interior_phase_root_scout.py",),
        expected=("validated Newman C1 calibrated roster-interior phase-root scout: 12 rows, 0 issues, 6 source audits, 3 precision ladders, 16 independent numerical audits, 3 physical q=1 cells, 3 calibrated B=-49/100 roots, 3 fixed-roster margins, 0 determinant rows, 1 open evaluator row",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 fixed-cell low-N determinant evaluator gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_fixed_cell_low_n_determinant_evaluator_gate.py",),
        expected=("validated fixed-cell low-N determinant evaluator gate: 20 rows, 3 fixtures, 184 finite modes, maximum alternate-basis error 4.903e-12, maximum alternate-current error 1.734e-13, 0 physical signed bounds",),
        category="open theorem target hygiene",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 physical-observation polynomial compiler gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_physical_observation_polynomial_compiler_gate.py",),
        expected=(
            "independent physical-observation compiler check passed: 16 exact symbolic rows, degree-5 joined fixture, 3 low-N fixtures",
            "max lift identity error",
        ),
        category="exact theorem-search algebra",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 physical scalar and genuine-edge adapter gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_physical_scalar_edge_adapter_gate.py",),
        expected=(
            "validated Newman C1 physical scalar and genuine-edge adapter gate: 15 rows, 6 source audits, 3 physical roots, 10 scalar inputs per root, 2 independent chi routes, 4 direct edge-derivative checks per root, 12 base and 12 tangent rows, 3 degree-five edge-affine polynomials",
            "0 retained-observation rows and 0 signed determinant rows",
        ),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 shifted-Hardy diagnostic bridge gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_shifted_hardy_diagnostic_bridge_gate.py",),
        expected=("validated Newman C1 shifted-Hardy diagnostic bridge gate: 15 rows, 0 issues, 16 source/hash audits, 9 independent fixture audits, 6 kernel audits, 24 physical kernel fits, shift_step 0.01 retained, 0.02 and 0.04 rejected, 0 physical values and 1 open resumable-evaluator row",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy resumable-evaluator equivalence gate",
        command=("work/rh_compute/scripts/check_hardy_resumable_equivalence_fixture.py",),
        expected=("validated resumable Hardy equivalence fixtures: 2 heights, 192 journal records, exact stop/resume state, 6 parked invocations, 2 fail-closed resumes",),
        category="finite computational reproducibility",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy joint correlated-error ledger gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_joint_correlated_error_ledger_gate.py",),
        expected=("validated Newman C1 Hardy joint correlated-error ledger gate: 2 residual vectors, 15 shifts, 24 physical rows, 21 constant-annihilating rows, 21/24 projected errors increase, 6 cutoff replays, scalar extrapolation retired, 0 physical values and 1 open modewise-error theorem",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy modewise derivative handoff gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_modewise_derivative_handoff_gate.py",),
        expected=("validated Newman C1 Hardy modewise derivative handoff gate: 12 rows, 6 exact low modes, 7 derivative orders, 1 exact rough-tail constant, 4 primary-paper locations, 7 source markers, 2 conditional routes and 1 open external error theorem",),
        category="exact conditional reduction",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy exact low-mode projected-kernel gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_exact_low_mode_projected_kernel_gate.py",),
        expected=("validated Newman C1 Hardy exact low-mode projected-kernel gate: 24 rows, 84 exact annihilations, 6 low modes, max held-out error 9.32193472449e-11, 21/24 projected errors still increase, 0 physical values and 1 open derivative theorem",),
        category="finite diagnostic nonpromotion",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy local-defect obligation ledger gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_local_defect_obligation_ledger_gate.py",),
        expected=("validated Newman C1 Hardy local-defect obligation ledger gate: 16 rows, 6 paper-to-source mappings, 13 error components, 2 exact defect identities, 1 exact-shell bypass, 2 conditional derivative routes and 0 external error bounds",),
        category="exact conditional reduction",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy chain-telemetry source contract",
        command=("work/rh_compute/scripts/check_hardy_chain_telemetry_source_contract.py",),
        expected=("validated Hardy telemetry source boundary: accepted hash pinned, 16 observer blocks exact, direct parent/child defects instrumented",),
        category="finite computational reproducibility",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy chain-telemetry scout",
        command=("work/rh_compute/scripts/check_hardy_chain_telemetry_scout.py",),
        expected=("validated Hardy chain telemetry scout: 64 chains, 64 direct defects, 2 branches, exact evaluator state, 3 diagnostic routes and 0 uniform bounds",),
        category="finite diagnostic nonpromotion",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy exact-binary128 point-ball defect gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_point_ball_defect_gate.py",),
        expected=("validated Hardy point-ball defect gate: 64 exact binary128 chains, 64 zero exclusions, 2 Arb precisions, 1 pinned runtime and 0 uniform bounds",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy initial-level parity adapter gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_initial_level_adapter_gate.py",),
        expected=("validated Hardy initial-level adapter gate: 64 exact source maps, 6720 saved parity checks, 36 synthetic sign/parity fixtures and 1 open binary phase budget",),
        category="exact conditional reduction",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy cubic transformed-child Legendre adapter gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_cubic_child_legendre_adapter_gate.py",),
        expected=("validated Hardy cubic child Legendre adapter gate: 64 chains, 131 parity checks, 131 denominator checks, 20 formal series cases",),
        category="exact conditional reduction",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy local q selector-cell gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_q_selector_cell_gate.py",),
        expected=("validated Hardy q selector-cell gate: 64 cells, maximum 11 halvings, 6 PSI and 2 ERF paths per cell, 0 physical coverage claims",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy coefficient-to-physical-RAE roster gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_rae_coefficient_roster_gate.py",),
        expected=("validated Hardy coefficient-to-RAE roster gate: 32 pivots, 64 branches, spacing 420, continuous coverage open",),
        category="exact conditional reduction",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy q-cell physical-spacing obstruction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_q_cell_physical_spacing_gate.py",),
        expected=("validated Hardy q-cell physical-spacing gate: 62/62 phase projections disjoint, discrete-index route selected",),
        category="countermodel nonpromotion",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy full block-20 telemetry fixture",
        command=("work/rh_compute/scripts/check_hardy_block20_full_telemetry_fixture.py",),
        expected=("validated Hardy block-20 full telemetry fixture: 424 chains, 374 recurrences, exact state/displayed values",),
        category="finite computational reproducibility",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 discrete selector atlas gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.py",),
        expected=("validated Hardy block-20 discrete selector atlas: 424 classified calls, 374 q cells, 50 direct kernels, 0 unclassified",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 source special-function probe",
        command=("work/rh_compute/scripts/check_hardy_block20_special_function_probe.py",),
        expected=("validated Hardy block-20 special-function probe: 374 calls, 2992 source values, 1122 exact t1/t2/t4 bit replays",),
        category="finite computational reproducibility",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 rigorous special-function residual gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_special_function_residual_gate.py",),
        expected=("validated Hardy block-20 special-function residual gate: 2244 PSI, 748 ERF, max correlated q shift",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 rigorous direct-kernel gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_direct_kernel_gate.py",),
        expected=("validated Hardy block-20 direct-kernel gate: 50 kernels, max true/logged gap",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 corrected recurrence-defect gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_corrected_recurrence_defect_gate.py",),
        expected=("validated Hardy block-20 corrected recurrence-defect gate: 374 defects, 374 corrected nonzero, 232 improved/142 worsened",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 source t5 component probe",
        command=("work/rh_compute/scripts/check_hardy_block20_t5_component_probe.py",),
        expected=("validated Hardy block-20 t5 component probe: 374 calls, 913 intrinsic erfc values, 374 exact t5 bit replays",),
        category="finite computational reproducibility",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 rigorous t5 real-erfc residual gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_t5_real_erfc_residual_gate.py",),
        expected=("validated Hardy block-20 t5 real-erfc residual gate: 913 calls, max residual",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 t5-erfc-corrected recurrence-defect gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_t5_erfc_corrected_recurrence_defect_gate.py",),
        expected=("validated Hardy block-20 t5-erfc-corrected recurrence-defect gate: 374 nonzero defects, scale separation",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 native-q required-compensation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate.py",),
        expected=("validated Hardy block-20 native-q compensation gate: 374 nonzero targets, range",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 ecor redundancy gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_ecor_redundancy_gate.py",),
        expected=("validated Hardy block-20 ecor redundancy gate: exact quartic match, no executable phase use",),
        category="exact conditional reduction",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 displayed-W1 finite-ip completion gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_w1_finite_ip_completion_gate.py",),
        expected=("validated Hardy block-20 displayed-W1 finite-ip completion gate: 204 fully covered, 175 omitted terms, 374 residuals exclude zero",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 t5 sqrt2 literal-kind gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_t5_sqrt2_literal_kind_gate.py",),
        expected=("validated Hardy block-20 t5 sqrt2 literal-kind gate: 913 calls, 374 stacked residuals exclude zero",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 exact equation-(69) saddle-family gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_exact_eq69_saddle_family_gate.py",),
        expected=("validated Hardy block-20 exact equation-(69) saddle-family gate: 544 integrals, 229 improved, 374 residuals exclude zero",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 t2 L+1 source/paper mismatch gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_t2_lplus1_mismatch_gate.py",),
        expected=("validated Hardy block-20 t2 L+1 mismatch gate: 374 improved, max residual ratio below 0.02, 374 residuals exclude zero",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 W2--W5 source/paper correspondence gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_w2_w5_correspondence_gate.py",),
        expected=("validated Hardy block-20 W2--W5 correspondence gate: 165 W3, 209 W4, 374/374 residuals exclude zero",),
        category="finite interval certificate",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 cubic vertical-contour admissibility gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_cubic_vertical_contour_admissibility_gate.py",),
        expected=("validated cubic contour admissibility gate: 187 failing 67b, 187 failing 67c",),
        category="exact conditional reduction",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 cubic sector-contour repair gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_cubic_sector_contour_repair_gate.py",),
        expected=("validated cubic sector-contour repair gate: 374 strict-convexity calls, 187+187 sign-aware slanted rays",),
        category="exact conditional reduction",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 logarithmic endpoint-ray reduction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_reduction_gate.py",),
        expected=("validated logarithmic endpoint-ray reduction gate: 374 calls, 4 diagnostic branches",),
        category="exact conditional reduction",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 logarithmic endpoint-ray Arb pilot gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate.py",),
        expected=("validated logarithmic endpoint-ray Arb pilot: 4/4 sign branches, both exact identities enclose zero",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 logarithmic endpoint-ray full Arb gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate.py",),
        expected=("validated logarithmic endpoint-ray full Arb gate: 374/374 calls, both exact identities enclose zero",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 endpoint-linearization decomposition pilot gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_pilot_gate.py",),
        expected=("validated endpoint-linearization decomposition pilot: 4/4 branches, both exact identities enclose zero",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 endpoint-linearization decomposition full gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_endpoint_linearization_decomposition_full_gate.py",),
        expected=("validated endpoint-linearization decomposition full gate: 374/374 calls, both exact identities enclose zero",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 paired exceptional-cancellation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_paired_exceptional_cancellation_gate.py",),
        expected=("validated paired exceptional-cancellation gate: 374/374 identities, worst paired ratio below 211",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 exceptional-mode recombination pilot gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_pilot_gate.py",),
        expected=("validated exceptional-mode recombination pilot: 4/4 branches, W2--W4 identities and direct modes certified",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 exceptional-mode recombination full gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_full_gate.py",),
        expected=("validated exceptional-mode recombination full gate: 374/374 model and recombination identities",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 shared-contour exceptional majorant gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_shared_contour_exceptional_majorant_gate.py",),
        expected=("validated shared-contour exceptional majorants: 374/374 W2 and 165/165 W3 bounds",),
        category="exact conditional reduction",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 W4 contour endpoint Hurwitz majorant gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_w4_contour_endpoint_hurwitz_majorant_gate.py",),
        expected=("validated W4/Hurwitz endpoint majorants: 209/209 W4 contours and 374/374 complete bounds",),
        category="exact conditional reduction",
    ),
    GateSpec(
        name="Hardy block-20 exact accumulation-weight fixture gate",
        command=("work/rh_compute/scripts/check_hardy_block20_accumulation_weight_fixture.py",),
        expected=("validated block-20 accumulation weights: 3180 rows, exact evaluator equivalence",),
        category="source equivalence",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 transported endpoint budget gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate.py",),
        expected=("validated transported endpoint budget: 424 calls, 15 outputs, absolute majorant exceeds 0.005 on every output",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 retained-decay endpoint majorant gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_retained_decay_endpoint_majorant_gate.py",),
        expected=("validated retained-decay endpoint majorants: 374 calls, 15/15 transported outputs below 0.005",),
        category="exact conditional reduction",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 coefficient-neighborhood transport gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate.py",),
        expected=("validated coefficient-neighborhood transport: 374 cells, 15/15 endpoint outputs below 0.005, other columns open",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 complete cubic Legendre-tail budget gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_complete_cubic_legendre_tail_budget_gate.py",),
        expected=("validated complete cubic Legendre-tail budget: 374 calls, 753 child indices, 15/15 combined outputs below 0.005",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 coefficient-cell complete cubic Legendre-tail budget gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_complete_cubic_legendre_tail_budget_gate.py",),
        expected=("validated coefficient-cell complete cubic Legendre-tail budget: 374 cells, 753 child indices, 15/15 combined outputs below 0.005",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 coefficient-cell joint recurrence identity gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_joint_recurrence_identity_gate.py",),
        expected=("validated coefficient-cell joint recurrence identity: 374 cells, 15/15 outputs below 0.005",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy block-20 source-to-corrected-model bridge gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_source_to_corrected_model_bridge_gate.py",),
        expected=("validated source-to-corrected-model bridge: 374 calls, 15/15 outputs below 0.005",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 cross-block structural atlas gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate.py",),
        expected=("validated t=1e10 cross-block structural atlas: 6784 calls, 1414 recursive, blocks 29-35 all direct",),
        category="finite theorem-search diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 cross-block retained endpoint gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_crossblock_retained_endpoint_gate.py",),
        expected=("validated t=1e10 cross-block retained endpoint gate: 1414 calls, blocks 20-28, max 4.87965760106136006892983729083960062751062901236036E-3",),
        category="interval certificates",
    ),
    GateSpec(
        name="Hardy t=1e10 cross-block accumulation-weight fixture",
        command=("work/rh_compute/scripts/check_hardy_crossblock_accumulation_weight_fixture.py",),
        expected=("validated cross-block accumulation weights: 50880 rows, exact evaluator equivalence",),
        category="source-equivalent fixture",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 cross-block recursive output budget gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate.py",),
        expected=("validated t=1e10 cross-block recursive output budget: 1414 calls, max 5.80580088833578459106118035839597594671279405744097E-3, 0/15 below 0.005",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 cross-block direct-kernel output gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_crossblock_direct_kernel_output_gate.py",),
        expected=("validated t=1e10 cross-block direct-kernel output gate: 5370 direct, combined max 5.80580088833578459106118039690461177653815320721018E-3, 0/15 below 0.005",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 cross-block hybrid exact output gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_crossblock_hybrid_exact_output_gate.py",),
        expected=("validated t=1e10 cross-block hybrid exact output gate: 6784 calls, max 3.80729575575206586213444515889870616774671274440825E-3, 15/15 below 0.005",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 block-21 signed bridge pilot gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_pilot_gate.py",),
        expected=("validated block-21 signed bridge pilot: witnesses=4, max_formula_gap=3.57018059931571603293765093667389010079205036163330E-14",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 complete block-21 signed bridge gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_full_gate.py",),
        expected=("validated complete block-21 signed bridge: calls=307, promoted_max=2.47142888476766613852340788184259027151181401602216E-3",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 blocks-22--28 geometry bridge pilot gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_pilot_gate.py",),
        expected=("validated blocks22--28 bridge geometry pilot: witnesses=18, max_gap=2.33303854962671323482192864418038880103267729282379E-14",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 all-recursive exact-point bridge gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_full_gate.py",),
        expected=("validated all-recursive exact-point bridge: recursive=1414, max=4.43988874036384893269534544300356252002080618765265E-4",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 later coefficient-cell signed bridge pilot gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_pilot_gate.py",),
        expected=("validated later coefficient-cell signed bridge pilot: 2/2 nonzero cells, 9 algebra checks, finite margin preserved",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 complete later coefficient-cell signed bridge gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_full_gate.py",),
        expected=("validated complete later coefficient-cell signed bridge: 1040/1040 cells, 15/15 outputs, max=4.48798733005244497141947844440854195181195374505537E-4",),
        category="interval certificates",
    ),
    GateSpec(
        name="Hardy pinned binary128 rounding-mode probe",
        command=("work/rh_compute/scripts/check_hardy_binary128_rounding_mode_probe.py",),
        expected=("validated binary128 rounding-mode probe: radix=2, digits=113, mode=nearest",),
        category="source reproduction",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 later complete-recurrence stress-cell pilot gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_stress_cell_pilot_gate.py",),
        expected=("validated later complete-recurrence stress-cell pilot: 4/4 cells, 15/15 outputs, max=4.56639665774875056167646386825368796568692156096651E-4",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 complete later-cell recurrence-tail gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_cell_full_gate.py",),
        expected=("validated complete later recurrence-tail cell transport: 1040/1040 cells, 15/15 outputs, max=1.18534368330532488984495230980059347380722157622465E-3",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 rigorous Arb source-q dependency-cut gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_rigorous_arb_source_q_dependency_cut_gate.py",),
        expected=("validated rigorous Arb source-q dependency cut: 1040 cells, 4 stripped replays, 15 outputs, max=1.18534368330532563659717863467183140923914676684614E-3",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 outer-Hardy representation dependency-ledger gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_outer_hardy_representation_dependency_ledger_gate.py",),
        expected=("validated outer-Hardy dependency audit: 15 outputs, exact=4, enclosed=1, open=13, Hardy-certified=0",),
        category="diagnostics",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 exact Hardy Arb calibration gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_exact_hardy_arb_calibration_gate.py",),
        expected=("validated exact Hardy Arb calibration: 15/15 reject 0.005",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 lower-RS versus hybrid-upper component-split gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_lower_rs_hybrid_upper_component_split_gate.py",),
        expected=("validated lower-RS/hybrid-upper split: 15 identities, lower-dominant=0, upper-dominant=15, upper-rejects-0.005=15",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 classical upper-main remainder-split gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_classical_upper_main_remainder_split_gate.py",),
        expected=("validated classical upper-main split: 15 identities, remainder-subordinate=15, hybrid-main-rejects-0.005=15",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 equation-124 classical block-partition gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate.py",),
        expected=("validated equation-124 block partition: 540 rows, block0-rejects-0.005=15, later-net-cancels=15",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 block-zero source-formula identity gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_block_zero_source_formula_identity_gate.py",),
        expected=("validated block-zero source identity: 1515 terms, transition-zero=15, source-model<1e-20=15, eq62-rejects-0.005=15",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 equation-62 cell-partition gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_block_zero_equation62_cell_partition_gate.py",),
        expected=("validated equation-(62) cell partition independently: roots=1530, prefixes=1515, final-midpoint<0.005=15/15",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 pole-free radius five-saddle collar gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation62_pole_free_radius_saddle_collar_gate.py",),
        expected=("validated pole-free radius saddle collar independently: outputs=15, crossed=[37941, 37942, 37943, 37944, 37945], pole-roster-unchanged=15",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 five-saddle D-transport gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation62_five_saddle_D_transport_gate.py",),
        expected=("validated five-saddle D transport independently: integrals=75, gap>0.005=15/15",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 signed midpoint-remainder decomposition gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation62_signed_midpoint_remainder_decomposition_gate.py",),
        expected=("validated signed midpoint remainder decomposition independently: outputs=15, opposite-sign=15, signed<0.005=15, triangle>0.005=15",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 source-aligned hybrid-error ledger gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_source_aligned_hybrid_error_ledger_gate.py",),
        expected=("validated source-aligned hybrid error ledger independently: outputs=15, blocks=540, open-first=source_aligned_upper_block_sum",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 equation-10 global theta-current reduction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation10_global_theta_current_reduction_gate.py",),
        expected=("validated equation-(10) global theta-current reduction independently: blocks=36, alpha=2481423, pairs=7034",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 finite-Poisson portcullis saddle reduction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.py",),
        expected=("validated finite-Poisson portcullis saddle reduction independently: interior=39231, turning=84, endpoint-classical=42",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 symmetric finite-Poisson interchange gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.py",),
        expected=("validated symmetric finite-Poisson interchange independently: alpha=2481423, uniform_tail=True, quantitative=False",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 Fresnel endpoint normal-form gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_fresnel_endpoint_normal_form_gate.py",),
        expected=("validated Fresnel endpoint normal form independently: modes=84, chart<0.044, ranges=6",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 universal logistic Morse characteristic-fold reduction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.py",),
        expected=("validated universal logistic Morse fold reduction independently: modes=84, split=42+42, closest<6.27e-6",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 lower-endpoint logit Airy-core reduction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_endpoint_logit_airy_core_reduction_gate.py",),
        expected=("validated lower-endpoint logit Airy core independently: lambda=5.1099, roots=2, phase<0.001",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 symmetric endpoint-tail reassembly gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_symmetric_endpoint_tail_reassembly_gate.py",),
        expected=("checked exact symmetric endpoint-tail residual and 39273-mode partition",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 half-Kummer reflection branch-reduction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_reflection_branch_reduction_gate.py",),
        expected=("checked odd-alpha reflection and one-branch half-domain saddle ledger",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 grouped tangent-B profile barrier gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_kummer_B_crossing_tangent_profile_barrier_gate.py",),
        expected=("checked B tangent profile in reverse order at 140 digits using direct erfc tails",),
        category="countermodel gates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 midpoint interpolation guard",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_boundary_linear_roster_pair_anchor_gate.py",),
        expected=("checked midpoint interpolation guard at 130 digits using direct erfc Fresnel values",),
        category="countermodel gates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 paired half-domain target-residual gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_paired_target_residual_normal_form_gate.py",),
        expected=("checked paired half-domain target residual by independent real-coordinate replay",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 signed-pair half-domain Morse-transform gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_signed_pair_logistic_morse_transform_gate.py",),
        expected=("checked signed-pair half-domain Morse transform by signed-variable replay",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 endpoint-driven pair-transport gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_endpoint_transport_volterra_gate.py",),
        expected=("checked endpoint-driven pair transport by independent Fourier/PDE replay",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 incomplete-Gamma outer-kernel gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_volterra_incomplete_gamma_outer_kernel_gate.py",),
        expected=("checked outer Volterra kernel by alternate Jacobian and 100-digit tail replay",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 separable pair-triangle geometry gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate.py",),
        expected=("checked pair-triangle geometry by independent gap and face-sign replay",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 exact bi-Morse face-fold gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_exact_bimorse_face_fold_gate.py",),
        expected=("checked exact bi-Morse face-fold chart by hyperbolic replay",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 face-trace lattice deghosting gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_face_trace_lattice_deghosting_gate.py",),
        expected=("checked face-trace lattice geometry and mode-257 deghosting at 100 digits",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 exact B step-tail normal-form gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_positive_B_crossing_step_tail_normal_form_gate.py",),
        expected=("checked exact B step-tail normal form by truth-table and partial-fraction replay",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 Abel-theta modular dual-roster gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_modular_dual_roster_gate.py",),
        expected=("checked Abel-theta modular roster in independent face-distance coordinates",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 B face-window Fresnel-remainder gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate.py",),
        expected=("checked B face window and Fresnel remainder independently",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 pole-subtracted B cotangent-current gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_pole_subtracted_cotangent_current_gate.py",),
        expected=("checked pole-subtracted B cotangent current by direct-series replay",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 discrete B-face joint-gradient gap gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_discrete_joint_gradient_gap_gate.py",),
        expected=("checked discrete B-face joint-gradient gap by minimax replay",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 exact local B quadratic-domain gap gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_exact_quadratic_domain_gap_gate.py",),
        expected=("checked exact local B quadratic-domain gaps by Lambert-free replay",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 dual-continuation stationary-cancellation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_continuation_stationary_cancellation_gate.py",),
        expected=("checked dual-continuation cancellation in face-distance coordinates with rational pi bounds",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 fixed-label Weber-completion gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_weber_completion_cancellation_gate.py",),
        expected=("checked Weber completion in ratio coordinates and the convergent half-plane",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 dual-cutoff shift-defect gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_cutoff_shift_defect_vanishing_gate.py",),
        expected=("checked dual-cutoff shift-defect decay in the complementary s coordinate",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 Weber source-roster tail-reassembly gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_weber_source_roster_tail_reassembly_gate.py",),
        expected=("checked source normalization and theta-tail reassembly in reciprocal-tail coordinates",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 double-Weber source-reconstruction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_double_weber_exact_source_reconstruction_gate.py",),
        expected=("checked double-Weber cancellation and endpoint bookkeeping independently",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 local B normal-safe Fresnel-remainder gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_normal_safe_fresnel_remainder_gate.py",),
        expected=("independently checked local B normal-safe remainders with 3/2 distance panels",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 signed pole-subtracted B cotangent-window gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_pole_subtracted_cotangent_signed_window_gate.py",),
        expected=("independently validated signed B cotangent window projection",),
        category="interval certificates",
        slow=True,
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 higher B boundary-current window gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_higher_boundary_current_absolute_window_gate.py",),
        expected=("independently checked higher B face currents",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 B-face Fresnel/boundary truncation-dictionary gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_fresnel_boundary_truncation_dictionary_gate.py",),
        expected=("independently checked Fresnel/boundary truncation dictionary correction",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 B-face window partial-budget gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_window_partial_budget_gate.py",),
        expected=("independently checked partial B-window budget",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 exact local B-window current gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_outer_safe_exact_current_gate.py",),
        expected=("independent local B current replay",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 complete signed B trace-package window gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_complete_gaussian_window_signed_gate.py",),
        expected=("independently checked complete B trace-package window",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 global residual B-window extraction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_B_window_extraction_target_gate.py",),
        expected=("independently checked B trace-package allocation, window extraction, and grouped target arithmetic",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 B trace-package bulk-step scope guard",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_trace_package_bulk_step_scope_guard_gate.py",),
        expected=("independently checked B trace-package scope and nonlocal bulk-step ownership",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 endpoint-safe exterior B trace tangential-cutoff gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_tangential_exterior_cutoff_gate.py",),
        expected=("independently checked C2 corner cutoff, exterior tangential gaps, and global allocation",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 endpoint/source-hybrid nonidentification guard",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_endpoint_residual_source_hybrid_nonidentification_guard_gate.py",),
        expected=("checked endpoint/source nonidentification independently",),
        category="non-promotion guards",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 beta^-4 canonical Airy-derivative reduction gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_airy_derivative_reduction_gate.py",),
        expected=("independently checked beta^-4 Airy reduction",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 exact positive-X Hankel branch-factorization gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_canonical_exact_hankel_branch_factorization_gate.py",),
        expected=("independently checked exact Airy/Hankel split",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 all-event Hankel carrier-selection gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_hankel_carrier_event_selection_gate.py",),
        expected=("independently checked Hankel event selection",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 scaled Airy-branch amplitude-envelope gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_scaled_airy_branch_amplitude_envelope_gate.py",),
        expected=("independently checked scaled Airy amplitude envelopes",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 three-roster branch-reflection pairing gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_opposite_branch_roster_reflection_pairing_gate.py",),
        expected=("independently checked 399-mode branch pairing: 3 rosters",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 completed beta^-4 selector-projection gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.py",),
        expected=("independently validated beta^-4 completed selector branch projection",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 selected/logistic stationary-action bridge gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_logistic_stationary_action_bridge_gate.py",),
        expected=("independently validated stationary-action bridge on 399 modes",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 event-zero top-corridor transport-extension gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_top_corridor_transport_extension_gate.py",),
        expected=("independently validated event-0 theta<=2 transport extension",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 event-ordered 398-term pairing multiplier gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_event_ordered_398_pairing_multiplier_gate.py",),
        expected=("independently checked event-ordered 398-term pairing and leading multipliers",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 selector-centred event-zero transport gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_selector_center_recentered_transport_gate.py",),
        expected=("independently validated selector-centered event-zero transport",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 selected beta^-4 endpoint-coherent detuning-ODE gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_endpoint_coherent_detuning_ode_gate.py",),
        expected=("independently checked endpoint-coherent beta^-4 detuning ODE",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 selected beta^-4 detuning Airy Green-kernel gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_detuning_airy_green_kernel_gate.py",),
        expected=("independently checked selected detuning Airy Green kernel",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 detuning Airy-lattice cubic-resonance gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_airy_lattice_cubic_resonance_gate.py",),
        expected=("independently checked Airy-lattice cubic resonance normal form",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 detuning Airy-lattice grouped Green-packet gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_airy_lattice_grouped_green_packet_gate.py",),
        expected=("independently checked grouped Airy phase and Green packets",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 leading-defect weighted Airy Green-packet gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_leading_defect_weighted_airy_green_packet_gate.py",),
        expected=("independently checked leading-defect weighted Airy Green packets",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 weighted variation-of-constants roster-reassembly gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_variation_of_constants_roster_reassembly_gate.py",),
        expected=("independently checked weighted variation-of-constants roster reassembly",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 weighted-completion partial-sum obligation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_completion_partial_sum_obligation_gate.py",),
        expected=("independently checked weighted-completion partial-sum obligation",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 completed cutoff-family uniform Dirichlet gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_completed_cutoff_family_uniform_dirichlet_gate.py",),
        expected=("independently checked uniform completed cutoff-family Dirichlet bound",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 weighted opposite-branch Dirichlet gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_opposite_branch_dirichlet_gate.py",),
        expected=("independently checked weighted opposite-branch Dirichlet bound",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 weighted selected-branch grouped-integral gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_selected_branch_grouped_integral_gate.py",),
        expected=("independently checked grouped selected-branch leading-defect integral",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 exact completed-point weighted-profile obligation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_completed_point_weighted_profile_obligation_gate.py",),
        expected=("independently checked exact completed-point versus weighted-profile obligation",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 exact finite-t common-profile operator gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.py",),
        expected=("independently checked exact finite-t common-profile Fourier operator",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 beta^-6 common-profile leading-coefficient gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta6_common_profile_leading_coefficient_gate.py",),
        expected=("independently checked beta^-6 common-profile leading coefficient",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 exact common-profile initial-slope gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_profile_initial_slope_gate.py",),
        expected=("independently generated degree-17/radius-11 initial-slope recomputation",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 exact Kummer-ODE whole-profile closure gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_kummer_ode_profile_closure_gate.py",),
        expected=("independent exact algebra and 16384-panel Kummer-ODE profile closure check",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 exact corrected selected-defect splice gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_corrected_selected_defect_splice_gate.py",),
        expected=("independent 23-mode algebra and corrected selected-defect sign arithmetic",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 corrected selected source-projection orientation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_corrected_selected_source_projection_orientation_gate.py",),
        expected=("independent parity, normalization, and source-projected corrected sign check",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 fold-owned target-allocation embedding-obligation gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_residual_target_allocation_embedding_obligation_gate.py",),
        expected=("independent inverse-map and four-term fold-allocation check; embedding remains open",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 selector-increment paired-sector embedding gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_increment_paired_sector_embedding_gate.py",),
        expected=("independent selector-cocycle and paired-sector increment check; remainder remains open",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 weighted-defect raw-strip companion decomposition gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_defect_raw_strip_companion_decomposition_gate.py",),
        expected=("independent weighted-defect/raw-strip companion and selector-increment check; joint remainder remains open",),
        category="exact theorem-search analysis",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 inner/outer corrected-defect projection gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_inner_outer_corrected_defect_projection_gate.py",),
        expected=("alternate-grid inner/outer corrected-defect projection replay; paired companion remains open",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 paired-increment two-packet positive-barrier gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_increment_two_packet_positive_barrier_gate.py",),
        expected=("alternate-grid two-packet positive selector-increment barrier; fixed-state residual remains open",),
        category="interval certificates",
    ),
    GateSpec(
        name="Jensen-window PF Newman C1 Hardy t=1e10 selector completion-sector nonlocality gate",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_selector_completion_sector_nonlocality_gate.py",),
        expected=("independent selector completion-sector barrier; global fixed-state endpoint residual remains open",),
        category="interval certificates",
    ),
)


def command_for(spec: GateSpec) -> list[str]:
    return [sys.executable, *spec.command]


def tail(text: str, lines: int = 12) -> str:
    parts = text.splitlines()
    return "\n".join(parts[-lines:])


def run_gate(spec: GateSpec, timeout: int) -> dict:
    start = perf_counter()
    try:
        completed = subprocess.run(
            command_for(spec),
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        elapsed = perf_counter() - start
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode(errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors="replace")
        return {
            "name": spec.name,
            "category": spec.category,
            "command": " ".join(command_for(spec)),
            "returncode": None,
            "elapsed_seconds": round(elapsed, 3),
            "ok": False,
            "timed_out": True,
            "missing_expected": list(spec.expected),
            "stdout_tail": tail(stdout),
            "stderr_tail": tail(stderr),
        }
    elapsed = perf_counter() - start
    combined = completed.stdout + "\n" + completed.stderr
    missing = [needle for needle in spec.expected if needle not in combined]
    ok = completed.returncode == 0 and not missing
    return {
        "name": spec.name,
        "category": spec.category,
        "command": " ".join(command_for(spec)),
        "returncode": completed.returncode,
        "elapsed_seconds": round(elapsed, 3),
        "ok": ok,
        "timed_out": False,
        "missing_expected": missing,
        "stdout_tail": tail(completed.stdout),
        "stderr_tail": tail(completed.stderr),
    }


def resolved_control_path(path: Path | None) -> Path | None:
    if path is None or path.is_absolute():
        return path
    return REPO_ROOT / path


def registry_signature(specs: list[GateSpec]) -> str:
    rows = []
    for spec in specs:
        script = REPO_ROOT / spec.command[0]
        script_sha256 = None
        if script.is_file():
            script_sha256 = hashlib.sha256(script.read_bytes()).hexdigest()
        rows.append(
            {
                "name": spec.name,
                "command": list(spec.command),
                "expected": list(spec.expected),
                "category": spec.category,
                "slow": spec.slow,
                "script_sha256": script_sha256,
            }
        )
    encoded = json.dumps(
        rows, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def checkpoint_payload(
    *,
    specs: list[GateSpec],
    results: list[dict],
    signature: str,
    timeout: int,
    skip_slow: bool,
    parked: bool,
) -> dict:
    return {
        "schema_version": 1,
        "kind": "core_proof_programme_replay_checkpoint",
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "registry_signature": signature,
        "skip_slow": skip_slow,
        "per_gate_timeout_seconds": timeout,
        "total_gates": len(specs),
        "completed_gates": len(results),
        "next_gate": (
            specs[len(results)].name
            if len(results) < len(specs)
            else None
        ),
        "parked": parked,
        "complete": len(results) == len(specs) and not parked,
        "all_completed_ok": all(result["ok"] for result in results),
        "results": results,
    }


def write_checkpoint(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def load_checkpoint(
    path: Path,
    *,
    specs: list[GateSpec],
    signature: str,
    skip_slow: bool,
) -> list[dict]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("kind") != "core_proof_programme_replay_checkpoint":
        raise ValueError("checkpoint kind does not match this runner")
    if payload.get("registry_signature") != signature:
        raise ValueError(
            "checkpoint registry signature does not match current gates"
        )
    if payload.get("skip_slow") is not skip_slow:
        raise ValueError("checkpoint --skip-slow mode does not match")
    results = payload.get("results")
    if not isinstance(results, list):
        raise ValueError("checkpoint results are not a list")
    if len(results) > len(specs):
        raise ValueError("checkpoint contains too many gate results")
    expected_names = [spec.name for spec in specs[: len(results)]]
    actual_names = [result.get("name") for result in results]
    if actual_names != expected_names:
        raise ValueError("checkpoint results are not a registry prefix")
    if not all(result.get("ok") is True for result in results):
        raise ValueError(
            "checkpoint contains a failed gate; start a fresh replay "
            "after fixing it"
        )
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=int, default=600, help="Per-gate timeout in seconds.")
    parser.add_argument("--skip-slow", action="store_true", help="Skip gates marked slow.")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON summary.")
    parser.add_argument(
        "--checkpoint",
        type=Path,
        help=(
            "Atomically save each completed gate so a parked serial replay "
            "can resume."
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume the validated prefix in --checkpoint.",
    )
    parser.add_argument(
        "--stop-file",
        type=Path,
        help=(
            "Park after the current gate when this file exists. "
            "A parked run exits with status 75."
        ),
    )
    return parser


def main() -> int:
    set_below_normal_priority()
    args = build_parser().parse_args()
    if args.resume and args.checkpoint is None:
        raise SystemExit("--resume requires --checkpoint")
    specs = [spec for spec in GATES if not (args.skip_slow and spec.slow)]
    signature = registry_signature(specs)
    checkpoint = resolved_control_path(args.checkpoint)
    stop_file = resolved_control_path(args.stop_file)
    results: list[dict] = []
    if args.resume:
        if checkpoint is None or not checkpoint.is_file():
            raise SystemExit("resume checkpoint does not exist")
        try:
            results = load_checkpoint(
                checkpoint,
                specs=specs,
                signature=signature,
                skip_slow=args.skip_slow,
            )
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            raise SystemExit(f"cannot resume checkpoint: {exc}") from exc

    parked = False
    if checkpoint is not None:
        write_checkpoint(
            checkpoint,
            checkpoint_payload(
                specs=specs,
                results=results,
                signature=signature,
                timeout=args.timeout,
                skip_slow=args.skip_slow,
                parked=False,
            ),
        )

    for spec in specs[len(results) :]:
        if stop_file is not None and stop_file.exists():
            parked = True
            break
        results.append(run_gate(spec, args.timeout))
        if checkpoint is not None:
            write_checkpoint(
                checkpoint,
                checkpoint_payload(
                    specs=specs,
                    results=results,
                    signature=signature,
                    timeout=args.timeout,
                    skip_slow=args.skip_slow,
                    parked=False,
                ),
            )
        if stop_file is not None and stop_file.exists():
            parked = True
            break

    ok = (
        not parked
        and len(results) == len(specs)
        and all(result["ok"] for result in results)
    )
    if checkpoint is not None:
        write_checkpoint(
            checkpoint,
            checkpoint_payload(
                specs=specs,
                results=results,
                signature=signature,
                timeout=args.timeout,
                skip_slow=args.skip_slow,
                parked=parked,
            ),
        )

    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "parked": parked,
                    "completed_gates": len(results),
                    "total_gates": len(specs),
                    "gates": results,
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for result in results:
            status = "OK" if result["ok"] else "FAIL"
            print(
                f"{status} core gate: {result['name']} "
                f"({result['category']}, {result['elapsed_seconds']}s)"
            )
            if not result["ok"]:
                if result["timed_out"]:
                    print(f"  timed out after {args.timeout}s")
                if result["missing_expected"]:
                    print(f"  missing expected: {result['missing_expected']}")
                if result["stdout_tail"]:
                    print("  stdout tail:")
                    print(result["stdout_tail"])
                if result["stderr_tail"]:
                    print("  stderr tail:")
                    print(result["stderr_tail"])
        if parked:
            print(
                "parked core proof-programme replay after "
                f"{len(results)}/{len(specs)} gates"
            )
        else:
            print(
                "validated "
                f"{sum(1 for result in results if result['ok'])}/"
                f"{len(specs)} core proof-programme gates"
            )

    if parked:
        return 75
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
