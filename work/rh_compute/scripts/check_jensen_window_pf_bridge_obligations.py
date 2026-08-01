#!/usr/bin/env python3
"""Validate the Jensen-window PF bridge obligation ledger."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_LEDGER = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_bridge_obligations.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_bridge_obligations.md"


ALLOWED_STATUSES = {
    "available_exact",
    "finite_evidence",
    "guard_validated",
    "open_obligation",
    "conditional_consequence",
    "rejected_by_countermodel",
    "route_mismatch",
}

ALLOWED_ROLES = {
    "exact_equivalence",
    "exact_contact",
    "algebraic_obstruction",
    "finite_evidence",
    "open_antecedent",
    "open_bridge_theorem",
    "open_uniformity_obligation",
    "conditional_consequence",
    "rejection_test",
    "route_separation",
}

ROLE_STATUSES = {
    "exact_equivalence": {"available_exact"},
    "exact_contact": {"guard_validated"},
    "algebraic_obstruction": {"guard_validated"},
    "finite_evidence": {"finite_evidence"},
    "open_antecedent": {"open_obligation"},
    "open_bridge_theorem": {"open_obligation"},
    "open_uniformity_obligation": {"open_obligation"},
    "conditional_consequence": {"conditional_consequence"},
    "rejection_test": {"rejected_by_countermodel"},
    "route_separation": {"route_mismatch"},
}

NON_CLOSING_ROLES = {
    "exact_equivalence",
    "exact_contact",
    "algebraic_obstruction",
    "finite_evidence",
    "conditional_consequence",
    "rejection_test",
    "route_separation",
}

REQUIRED_IDS = {
    "jwpf_01_window_pf_jensen_equivalence",
    "jwpf_02_degree2_signed_hankel_contact",
    "jwpf_03_low_degree_extra_toeplitz_obligations",
    "jwpf_04_current_finite_pf_sturm_evidence",
    "jwpf_05_all_order_shifted_sign_consistency",
    "jwpf_05b_weaker_xi_specific_antecedent",
    "jwpf_06_sign_regular_to_jensen_pf_conversion",
    "jwpf_07_binomial_weight_and_shift_uniformity",
    "jwpf_08_jensen_to_laguerre_polya_limit",
    "jwpf_09_finite_rectangle_promotion_rejected",
    "jwpf_10_ordinary_coefficient_pf_route_separated",
    "jwpf_11_edrei_log_stieltjes_endpoint_equivalence",
    "jwpf_12_phi_pick_kernel_endpoint_equivalence",
    "jwpf_13_generic_backward_stieltjes_invariance_rejected",
    "jwpf_14_suzuki_arithmetic_hankel_equivalence",
    "jwpf_15_suzuki_determinant_only_reduction",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Bridge Obligations",
    "Status: theorem-obligation ledger",
    "This is not a proof of PF-infinity",
    "work/rh_compute/results/jensen_window_pf_bridge_obligations.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_bridge_obligations.py",
    "validated Jensen-window PF bridge obligations: 16 obligations, 0 issues, 3 open obligations",
    "jwpf_06_sign_regular_to_jensen_pf_conversion",
    "jwpf_05b_weaker_xi_specific_antecedent",
    "outputs/jensen_window_pf_quartic_outer_contact_normal_form_gate.md",
    "outputs/jensen_window_pf_quartic_outer_contact_length14_survivor_gate.md",
    "outputs/jensen_window_pf_quartic_quintic_polar_contact_lemma.md",
    "outputs/jensen_window_pf_newman_zero_slab_degree71_sector_certificate.md",
    "work/rh_compute/results/jensen_window_pf_newman_zero_slab_degree71_sector_certificate.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_newman_zero_slab_degree71_sector_certificate.py",
    "outputs/jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.md",
    "work/rh_compute/results/jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.py",
    "through degree 361 uniformly on 0<=t<=1/5",
    "does not prove degree 362",
    "does not supply the all-degree antecedent",
    "negative-discriminant outward collar",
    "force the contact",
    "Q_(10,n)(-100)<0",
    "outputs/jensen_window_pf_endpoint_order10_counterexample.md",
    "would_close_target=true",
    "outputs/jensen_window_pf_bridge_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
    "outputs/jensen_window_pf_theorem_machinery_fit_matrix.md",
    "work/rh_compute/results/proof_claim_ledger.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py",
    "validated Jensen-window PF theorem machinery fit matrix: 11 rows, 0 issues, 0 ready-to-apply rows",
    "jwpf_11_edrei_log_stieltjes_endpoint_equivalence",
    "outputs/jensen_window_pf_edrei_stieltjes_equivalence_gate.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_edrei_stieltjes_equivalence_gate.py",
    "jwpf_12_phi_pick_kernel_endpoint_equivalence",
    "outputs/jensen_window_pf_phi_pick_kernel_target.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_phi_pick_kernel_target.py",
    "jwpf_13_generic_backward_stieltjes_invariance_rejected",
    "outputs/jensen_window_pf_edrei_heat_flow_boundary_gate.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_edrei_heat_flow_boundary_gate.py",
    "outputs/jensen_window_pf_edrei_hankel_boundary_flux_gate.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_edrei_hankel_boundary_flux_gate.py",
    "orthogonal-polynomial Schur quotient",
    "shifted determinant-flux criterion is exactly equivalent",
    "outputs/jensen_window_pf_edrei_raw_moment_collision_resolution_gate.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_edrei_raw_moment_collision_resolution_gate.py",
    "condition number of order",
    "gap-adapted exponentially",
    "jwpf_14_suzuki_arithmetic_hankel_equivalence",
    "jwpf_15_suzuki_determinant_only_reduction",
    "outputs/jensen_window_pf_xi_pick_suzuki_hankel_bridge.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_xi_pick_suzuki_hankel_bridge.py",
    "outputs/jensen_window_pf_suzuki_determinant_only_reduction.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_determinant_only_reduction.py",
    "outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_fixed_omega_phase_diagram.py",
    "outputs/jensen_window_pf_suzuki_jordan_totient_sign_scout.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_totient_sign_scout.py",
    "outputs/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.py",
    "outputs/jensen_window_pf_suzuki_jordan_error_kernel_reduction.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_error_kernel_reduction.py",
    "outputs/jensen_window_pf_suzuki_cofinal_l2_hierarchy.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_l2_hierarchy.py",
    "H_(omega,k)(exp(t))-P_(omega,k)(t) in L2(0,infinity)",
    "outputs/jensen_window_pf_jordan_muntz_causal_energy_bridge.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_causal_energy_bridge.py",
    "sup_N ||E_(omega,N)||_H",
    "outputs/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.py",
    "(1-2omega)||t^(-omega)f_(2omega,N)||_2",
    "outputs/jensen_window_pf_burnol_cell_energy_tail_obstruction.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_burnol_cell_energy_tail_obstruction.py",
    "Q_(omega,N)",
    "N^(-1/2-omega)",
    "short-multiplicative-interval",
    "outputs/jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.py",
    "Q_(omega,infinity)<infinity",
    "V_(omega,N)(2^j*N)",
    "(k/(q+1),k/q]",
    "outputs/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.md",
    "work/rh_compute/results/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.py",
    "stationary Gram",
    "signed Mobius off-diagonal",
    "kernel positivity",
    "outputs/jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.md",
    "work/rh_compute/results/jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.py",
    "Ornstein-Uhlenbeck",
    "reciprocal-Mobius tail energy",
    "sign-indefinite remainder",
    "outputs/jensen_window_pf_ou_mertens_mean_square_reduction.md",
    "work/rh_compute/results/jensen_window_pf_ou_mertens_mean_square_reduction.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_ou_mertens_mean_square_reduction.py",
    "dyadic Mertens mean square",
    "origin-anchored",
    "averaged-Chowla",
    "outputs/jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.md",
    "work/rh_compute/results/jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.py",
    "constant anchor mode",
    "B_1(X)-TI_X+TII_X",
    "-TI_X+TII_X",
    "endogenous",
    "outputs/jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.md",
    "work/rh_compute/results/jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.py",
    "sum_(N>=1)P_(alpha/2,N)/N<infinity",
    "one logarithm short",
    "outputs/jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.md",
    "work/rh_compute/results/jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.py",
    "r>=ceil(sqrt(K))",
    "K^(1-alpha)",
    "outputs/jensen_window_pf_mertens_local_path_cosine_transference.md",
    "work/rh_compute/results/jensen_window_pf_mertens_local_path_cosine_transference.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_local_path_cosine_transference.py",
    "sqrt(K)r_K-tau_(K,0)",
    "zero dyadic endpoint tails",
    "outputs/jensen_window_pf_mertens_centered_bridge_vaughan_handoff.md",
    "work/rh_compute/results/jensen_window_pf_mertens_centered_bridge_vaughan_handoff.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_centered_bridge_vaughan_handoff.py",
    "min(i,j)-ij/K",
    "C_K=-TI_K+TII_K",
    "outputs/jensen_window_pf_mertens_affine_tent_bridge_handoff.md",
    "work/rh_compute/results/jensen_window_pf_mertens_affine_tent_bridge_handoff.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_affine_tent_bridge_handoff.py",
    "u_K-(1/2)u_(2K)",
    "O_(loc,K)=-TI_(loc,K)+TII_(loc,K)",
    "beta_K",
    "H_K=O(K^(2+epsilon))",
    "outputs/jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.md",
    "work/rh_compute/results/jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.py",
    "B_K=-V_(I,K)+V_(II,K)",
    "outputs/jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.md",
    "work/rh_compute/results/jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.py",
    "X_(K,r)=-T_(I,K,r)+T_(II,K,r)",
    "outputs/jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.md",
    "work/rh_compute/results/jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.py",
    "R_(alpha,K)=ceil(K^(1-alpha/2))",
    "2K^(-1)+8K^(-alpha/2)",
    "outputs/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.md",
    "work/rh_compute/results/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.py",
    "7pi^2/24",
    "-2<u_I,G u_II>",
    "outputs/jensen_window_pf_mertens_shift_kernel_variation_handoff.md",
    "work/rh_compute/results/jensen_window_pf_mertens_shift_kernel_variation_handoff.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_shift_kernel_variation_handoff.py",
    "Var_(n in I_(K,h))W_(K,h)(n)<6rho_K",
    "joint cancellation across shifts",
    "outputs/jensen_window_pf_mertens_planar_abel_handoff.md",
    "work/rh_compute/results/jensen_window_pf_mertens_planar_abel_handoff.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_abel_handoff.py",
    "3*pi^2*(1+log(2R))/K",
    "E_(alpha,K)=O_epsilon(K^(1+epsilon))",
    "outputs/jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.md",
    "work/rh_compute/results/jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.py",
    "(27/2)*pi^2*K*(1+log(2R))",
    "|C_(alpha,K)|=O_epsilon(K^(1+epsilon))",
    "one-vertex collisions",
    "outputs/jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.md",
    "work/rh_compute/results/jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.py",
    "C_3=210",
    "O_(CC,int,K)",
    "sum_(r<=R)|Y_(K,r)|^2",
    "arXiv:2607.15574",
    "outputs/jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.md",
    "work/rh_compute/results/jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.py",
    "B_CC-A_CF",
    "sum_r Z_(K,r)/q_r^2",
    "E_perp+c(beta+gamma/c)^2",
    "outputs/jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.md",
    "work/rh_compute/results/jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.py",
    "R_B=(E_B-D_B)/2+C_delta",
    "D_B<5pi^2/12",
    "anchored `E_B`",
    "E_B<=pi(1/2+2pi)K^(-2)H_K",
    "outputs/jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.md",
    "work/rh_compute/results/jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.py",
    "O(K^(alpha/2))",
    "K^2 E_B/H_K->0",
    "RH-equivalent",
    "O_(CC,int)+E_B/2",
    "outputs/jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.md",
    "work/rh_compute/results/jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.py",
    "J_B=L/2+(D_B-D_L)/2-C_delta",
    "returns the old lossless gate",
    "outputs/jensen_window_pf_mertens_planar_curvature_energy_scout.md",
    "work/rh_compute/results/jensen_window_pf_mertens_planar_curvature_energy_scout.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_curvature_energy_scout.py",
    "max `K*V_K<6.571`",
    "do not prove a uniform bound",
    "outputs/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.md",
    "work/rh_compute/results/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.py",
    "inner fixed-shift quotient",
    "universal same-shift",
    "outputs/jensen_window_pf_structural_ansatz_matrix.md",
    "work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py",
    "validated Jensen-window PF structural ansatz matrix: 6 ansatz rows, 0 issues, 0 ready-to-apply rows",
)


@dataclass(frozen=True)
class ObligationIssue:
    obligation_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(obligation_id: str, ref: str) -> list[ObligationIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [ObligationIssue(obligation_id, "missing-ref", ref)]


def has_boundary(text: str) -> bool:
    lowered = text.lower()
    return any(
        marker in lowered
        for marker in (
            "not ",
            "not a ",
            "only",
            "open",
            "finite",
            "conditional",
            "rejected",
            "separate",
        )
    )


def validate_obligation(row: dict) -> list[ObligationIssue]:
    issues: list[ObligationIssue] = []
    obligation_id = str(row.get("id", "<missing-id>"))
    for key in ("id", "role", "status", "refs", "claim", "needed_upgrade", "would_close_target", "proof_boundary"):
        if key not in row:
            issues.append(ObligationIssue(obligation_id, "missing-field", key))

    role = row.get("role")
    status = row.get("status")
    if role not in ALLOWED_ROLES:
        issues.append(ObligationIssue(obligation_id, "bad-role", repr(role)))
    if status not in ALLOWED_STATUSES:
        issues.append(ObligationIssue(obligation_id, "bad-status", repr(status)))
    if role in ROLE_STATUSES and status not in ROLE_STATUSES[role]:
        issues.append(ObligationIssue(obligation_id, "role-status-mismatch", f"{role!r} -> {status!r}"))

    would_close = row.get("would_close_target")
    if not isinstance(would_close, bool):
        issues.append(ObligationIssue(obligation_id, "bad-would-close-target", repr(would_close)))
    if role in NON_CLOSING_ROLES and would_close is True:
        issues.append(ObligationIssue(obligation_id, "nonclosing-row-closes-target", role))
    if would_close is True and status != "open_obligation":
        issues.append(ObligationIssue(obligation_id, "closing-row-not-open", repr(status)))

    refs = row.get("refs")
    if not isinstance(refs, list) or not refs:
        issues.append(ObligationIssue(obligation_id, "missing-refs", "refs must be a nonempty list"))
    else:
        for ref in refs:
            if not isinstance(ref, str):
                issues.append(ObligationIssue(obligation_id, "bad-ref", repr(ref)))
            else:
                issues.extend(validate_ref(obligation_id, ref))

    boundary = str(row.get("proof_boundary", ""))
    if not has_boundary(boundary):
        issues.append(ObligationIssue(obligation_id, "weak-proof-boundary", boundary))

    text = f"{row.get('claim', '')} {row.get('needed_upgrade', '')} {boundary}".lower()
    if row.get("status") in {"finite_evidence", "guard_validated", "rejected_by_countermodel", "route_mismatch"}:
        for forbidden in ("therefore rh", "therefore lambda <= 0", "proves lambda <= 0"):
            if forbidden in text:
                issues.append(ObligationIssue(obligation_id, "forbidden-overclaim", forbidden))
    return issues


def validate_note(path: Path) -> list[ObligationIssue]:
    if not path.exists():
        return [ObligationIssue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ObligationIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(ObligationIssue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "the bridge is proved"):
        if forbidden in lowered:
            issues.append(ObligationIssue("note", "forbidden-text", forbidden))
    return issues


def validate(ledger_path: Path, note_path: Path) -> tuple[list[ObligationIssue], int, int]:
    ledger = load_json(ledger_path)
    issues: list[ObligationIssue] = []
    if ledger.get("kind") != "jensen_window_pf_bridge_obligation_ledger":
        issues.append(ObligationIssue("<ledger>", "bad-kind", repr(ledger.get("kind"))))
    if ledger.get("target") != "target_jensen_window_pf_bridge":
        issues.append(ObligationIssue("<ledger>", "bad-target", repr(ledger.get("target"))))
    boundary = str(ledger.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(ObligationIssue("<ledger>", "weak-proof-boundary", ledger.get("proof_boundary", "")))

    obligations = ledger.get("obligations", [])
    if not isinstance(obligations, list) or not obligations:
        issues.append(ObligationIssue("<ledger>", "missing-obligations", "obligations must be a nonempty list"))
        obligations = []

    seen: set[str] = set()
    open_count = 0
    closing_count = 0
    for row in obligations:
        if not isinstance(row, dict):
            issues.append(ObligationIssue("<ledger>", "bad-obligation", repr(row)))
            continue
        obligation_id = str(row.get("id", "<missing-id>"))
        if obligation_id in seen:
            issues.append(ObligationIssue(obligation_id, "duplicate-id", obligation_id))
        seen.add(obligation_id)
        if row.get("status") == "open_obligation":
            open_count += 1
        if row.get("would_close_target") is True:
            closing_count += 1
        issues.extend(validate_obligation(row))

    for missing in sorted(REQUIRED_IDS - seen):
        issues.append(ObligationIssue(missing, "missing-required-obligation", missing))
    if open_count < 3:
        issues.append(ObligationIssue("<ledger>", "too-few-open-obligations", str(open_count)))
    if closing_count != 1:
        issues.append(ObligationIssue("<ledger>", "bad-closing-obligation-count", str(closing_count)))
    if "jwpf_04_current_finite_pf_sturm_evidence" in seen:
        finite_row = next(row for row in obligations if isinstance(row, dict) and row.get("id") == "jwpf_04_current_finite_pf_sturm_evidence")
        boundary_text = f"{finite_row.get('claim', '')} {finite_row.get('proof_boundary', '')}".lower()
        if "not all-minor jensen-window pf-infinity" not in boundary_text:
            issues.append(
                ObligationIssue(
                    "jwpf_04_current_finite_pf_sturm_evidence",
                    "missing-finite-boundary",
                    "must say not all-minor Jensen-window PF-infinity",
                )
            )

    issues.extend(validate_note(note_path))
    return issues, len(obligations), open_count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, obligation_count, open_count = validate(args.ledger, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "obligations": obligation_count,
                    "open_obligations": open_count,
                    "issues": [asdict(issue) for issue in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JENSEN-WINDOW-PF-OBLIGATION {item.obligation_id} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF bridge obligations: "
            f"{obligation_count} obligations, {len(issues)} issues, {open_count} open obligations"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
