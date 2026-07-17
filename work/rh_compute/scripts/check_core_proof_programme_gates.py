#!/usr/bin/env python3
"""Run the core proof-programme reproducibility gates.

This is a compact smoke/ledger runner for the RH dynamical-rigidity corpus.
It does not prove PF-infinity, Laguerre-Polya membership, RH, or Lambda <= 0.
It checks that the promoted finite-evidence manifests and executable
countermodel guards still match the advertised proof-programme status.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import subprocess
import sys
from time import perf_counter


REPO_ROOT = Path(__file__).resolve().parents[3]


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
        expected=("validated proof-claim ledger: 324 claims, 0 issues, 9 open theorem targets",),
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
        name="Jensen-window PF bridge obligation ledger",
        command=("work/rh_compute/scripts/check_jensen_window_pf_bridge_obligations.py",),
        expected=("validated Jensen-window PF bridge obligations: 11 obligations, 0 issues, 3 open obligations",),
        category="open theorem target hygiene",
    ),
    GateSpec(
        name="Jensen-window PF theorem machinery fit matrix",
        command=("work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py",),
        expected=("validated Jensen-window PF theorem machinery fit matrix: 7 rows, 0 issues, 0 ready-to-apply rows",),
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
        name="Jensen-window PF quartic double-root threshold lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_quartic_double_root_threshold_lemma.py",),
        expected=("validated Jensen-window PF quartic double-root threshold lemma: 10 rows, 0 issues, 5 exact coordinate identities, 1 double-root splitting criterion, 1 branch-aware inward threshold, 1 triple-root equality, 1 tangent factor, 1 explained countermodel, 1 open global-invariant handoff",),
        category="exact theorem-search algebra",
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
        name="Jensen-window PF scaled double-zero boundary-layer lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_scaled_double_zero_boundary_layer_lemma.py",),
        expected=("validated Jensen-window PF scaled double-zero boundary-layer lemma: 10 rows, 0 issues, 3 exact scaling identities, 1 heat PDE, 1 double-zero transversality, 1 universal boundary layer, 1 D^(-3/2) gap law, 1 external-field D^(-2) collision law, 1 exact toy family, 1 threshold-exhaustion theorem, 1 open uniform handoff",),
        category="exact theorem-search algebra",
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
        name="Jensen-window PF Newman positive-boundary attainment lemma",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_positive_boundary_attainment_lemma.py",),
        expected=("validated Jensen-window PF Newman positive-boundary attainment lemma: 10 rows, 0 issues, 2 published compactness inputs, 1 finite-boundary attainment theorem, 1 positive-time simplicity equivalence, 1 arbitrary-multiplicity Hermite split, 9 exact Hermite checks, 1 cluster-energy blow-up, 1 open Xi endpoint handoff",),
        category="exact theorem-search algebra",
    ),
    GateSpec(
        name="Jensen-window PF Newman strict-Laguerre correlation target",
        command=("work/rh_compute/scripts/check_jensen_window_pf_newman_strict_laguerre_correlation_target.py",),
        expected=("validated Jensen-window PF Newman strict-Laguerre correlation target: 10 rows, 0 issues, 1 strict-Laguerre equivalence, 1 exact correlation identity, 1 Wiener-density equivalence, 1 RH-equivalent density target, 1 exact strict-log-concavity/positive-definiteness countermodel, 2 non-promotion gates, 1 open Xi handoff",),
        category="open theorem target hygiene",
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
)


def command_for(spec: GateSpec) -> list[str]:
    return [sys.executable, *spec.command]


def tail(text: str, lines: int = 12) -> str:
    parts = text.splitlines()
    return "\n".join(parts[-lines:])


def run_gate(spec: GateSpec, timeout: int) -> dict:
    start = perf_counter()
    completed = subprocess.run(
        command_for(spec),
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
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
        "missing_expected": missing,
        "stdout_tail": tail(completed.stdout),
        "stderr_tail": tail(completed.stderr),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=int, default=600, help="Per-gate timeout in seconds.")
    parser.add_argument("--skip-slow", action="store_true", help="Skip gates marked slow.")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON summary.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    specs = [spec for spec in GATES if not (args.skip_slow and spec.slow)]
    results = [run_gate(spec, args.timeout) for spec in specs]
    ok = all(result["ok"] for result in results)

    if args.json:
        print(json.dumps({"ok": ok, "gates": results}, indent=2, sort_keys=True))
    else:
        for result in results:
            status = "OK" if result["ok"] else "FAIL"
            print(
                f"{status} core gate: {result['name']} "
                f"({result['category']}, {result['elapsed_seconds']}s)"
            )
            if not result["ok"]:
                if result["missing_expected"]:
                    print(f"  missing expected: {result['missing_expected']}")
                if result["stdout_tail"]:
                    print("  stdout tail:")
                    print(result["stdout_tail"])
                if result["stderr_tail"]:
                    print("  stderr tail:")
                    print(result["stderr_tail"])
        print(f"validated {sum(1 for result in results if result['ok'])}/{len(results)} core proof-programme gates")

    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
