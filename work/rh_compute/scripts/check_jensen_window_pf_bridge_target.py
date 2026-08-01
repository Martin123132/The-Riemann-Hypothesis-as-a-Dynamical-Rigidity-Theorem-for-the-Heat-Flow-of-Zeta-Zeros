#!/usr/bin/env python3
"""Validate the Jensen-window PF bridge target specification."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_bridge_target.md"


@dataclass(frozen=True)
class TargetIssue:
    section: str
    issue: str
    detail: str


REQUIRED_STRINGS = (
    "Status: theorem target artifact",
    "This is not a proof of PF-infinity",
    "B^{d,n,lambda}_j = binom(d,j) A_{n+j}(lambda)",
    "for every d >= 1 and n >= 0",
    "the finite sequence B^{d,n,0}_0,...,B^{d,n,0}_d is PF-infinity",
    "Toeplitz total-positivity language",
    "infinite banded Toeplitz matrix",
    "T_{d,n}(r,c) = B^{d,n,0}_{c-r}",
    "Exact Relation To Jensen",
    "Contact With Signed Hankel",
    "d = 2",
    "A_{n+1}^2 - A_n A_{n+2} >= 0",
    "m = 1",
    "For `d >= 3`",
    "outputs/jensen_window_pf_obligation_algebra.md",
    "work/rh_compute/results/jensen_window_pf_obligation_algebra.json",
    "python work/rh_compute/scripts/jensen_window_pf_obligation_algebra.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_obligation_algebra.py",
    "selected degree-3 and degree-4",
    "outputs/arb_jensen_window_pf_obligation_diagnostic.md",
    "work/rh_compute/results/arb_jensen_window_pf_obligations_lamgrid_n0_n20_d3_m1_m8_d4_m1_m6_dps520_summary.json",
    "work/rh_compute/results/arb_jensen_window_pf_obligations_lamgrid_n0_n20_d3_m1_m8_d4_m1_m6_dps520.jsonl",
    "python work/rh_compute/scripts/arb_jensen_window_pf_obligation_probe.py",
    "python work/rh_compute/scripts/check_arb_jensen_window_pf_obligation_manifest.py",
    "validates `1470/1470` selected Arb interval determinants",
    "not all-minor Jensen-window PF-infinity",
    "outputs/arb_jensen_window_sturm_hyperbolicity_diagnostic.md",
    "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520_summary.json",
    "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520.jsonl",
    "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520_summary.json",
    "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520.jsonl",
    "python work/rh_compute/scripts/arb_jensen_window_sturm_hyperbolicity_probe.py",
    "python work/rh_compute/scripts/check_arb_jensen_window_sturm_manifest.py",
    "validates `210/210` degree-3/4 and `105/105` degree-5 Jensen-window",
    "not all-degree or all-shift Jensen hyperbolicity",
    "outputs/jensen_window_sturm_pf_consequence.md",
    "python work/rh_compute/scripts/check_jensen_window_sturm_pf_consequence.py",
    "finite Polya-frequency consequence",
    "`315/315`",
    "not all-minor Jensen-window PF-infinity as an infinite theorem",
    "outputs/jensen_hankel_bridge_algebra.md",
    "work/rh_compute/results/jensen_hankel_bridge_algebra.json",
    "python work/rh_compute/scripts/check_jensen_hankel_bridge_algebra.py",
    "outputs/jensen_window_pf_deep_schur_toda_boundary_gate.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_deep_schur_toda_boundary_gate.py",
    "H(z)=exp(z/100)/((1-z)(1-2z))",
    "coefficient specialization is strictly Schur-positive",
    "-222484532394597/2000000000000<0",
    "additional Xi/Phi-specific analytic structure",
    "outputs/jensen_window_pf_endpoint_order10_counterexample.md",
    "Q_(10,n)(-100)<0",
    "s_((N^10))(h)<0",
    "proposed all-shift signed-Hankel/deep-Schur cone",
    "outputs/jensen_window_pf_newman_zero_slab_degree71_sector_certificate.md",
    "work/rh_compute/results/jensen_window_pf_newman_zero_slab_degree71_sector_certificate.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_newman_zero_slab_degree71_sector_certificate.py",
    "Re H_t(x+iy)>0",
    "|x|<=84/5",
    "sin(delta)=840/7081",
    "through degree `71` hyperbolic",
    "degree `72`",
    "d<=9.36*10^20",
    "neither bounded result supplies the required all-degree",
    "Required Bridge Theorem",
    "an Xi/Phi-specific kernel, determinant, or variation-diminishing condition",
    "every binomially weighted Jensen window B^{d,n,0} is PF-infinity",
    "Obligation Ledger",
    "outputs/jensen_window_pf_bridge_obligations.md",
    "work/rh_compute/results/jensen_window_pf_bridge_obligations.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_bridge_obligations.py",
    "16 obligations",
    "3 open obligations",
    "jwpf_06_sign_regular_to_jensen_pf_conversion is the central open bridge theorem",
    "finite evidence rows have would_close_target=false",
    "jwpf_05_all_order_shifted_sign_consistency is rejected by counterexample",
    "jwpf_05b_weaker_xi_specific_antecedent is the open replacement antecedent",
    "outputs/jensen_window_pf_theorem_machinery_fit_matrix.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py",
    "11` source-anchored theorem-family rows",
    "ready-to-apply rows",
    "outputs/jensen_window_pf_structural_ansatz_matrix.md",
    "work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py",
    "`6` candidate, blocked, and rejected ansatz rows",
    "finite countermodel kill gate",
    "python work/rh_compute/scripts/check_arb_shifted_hankel_staircase_manifest.py",
    "3,154,515/3,154,515 finite shifted reshaped-Hankel minors",
    "outputs/countermodel_library.md",
    "python work/rh_compute/scripts/countermodel_gate_examples.py",
    "A_0..A_25",
    "positive `A_26`",
    "breaks the next degree-2 Jensen discriminant at shift `24`",
    "c_k(lambda) = A_k(lambda)/k!",
    "c is PF-infinity <=> every B^{d,n}_j is finite PF-infinity for all d,n",
    "outputs/jensen_window_pf_coefficient_pf_equivalence_gate.md",
    "work/rh_compute/results/jensen_window_pf_coefficient_pf_equivalence_gate.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_coefficient_pf_equivalence_gate.py",
    "a_r=p_(r+1)=(-1)^r[z^r]H'(z)/H(z)",
    "Sokal's logarithmic-derivative criterion",
    "a is a Stieltjes moment sequence",
    "outputs/jensen_window_pf_edrei_stieltjes_equivalence_gate.md",
    "work/rh_compute/results/jensen_window_pf_edrei_stieltjes_equivalence_gate.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_edrei_stieltjes_equivalence_gate.py",
    "Finite rows do not enter the equivalence",
    "P_Phi(z)=-Im(F_0'(z)*conj(F_0(z)))",
    "F_0 in LP+ <=> P_Phi(z)>0 for every Im(z)>0",
    "outputs/jensen_window_pf_phi_pick_kernel_target.md",
    "work/rh_compute/results/jensen_window_pf_phi_pick_kernel_target.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_phi_pick_kernel_target.py",
    "M(w)=xi((1+w)/2)/4",
    "hyperbolic directional derivative",
    "outputs/jensen_window_pf_xi_pick_suzuki_hankel_bridge.md",
    "work/rh_compute/results/jensen_window_pf_xi_pick_suzuki_hankel_bridge.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_xi_pick_suzuki_hankel_bridge.py",
    "outputs/jensen_window_pf_suzuki_spectral_frontier.md",
    "work/rh_compute/results/jensen_window_pf_suzuki_spectral_frontier.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_spectral_frontier.py",
    "outputs/jensen_window_pf_suzuki_determinant_only_reduction.md",
    "work/rh_compute/results/jensen_window_pf_suzuki_determinant_only_reduction.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_determinant_only_reduction.py",
    "outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md",
    "work/rh_compute/results/jensen_window_pf_suzuki_fixed_omega_phase_diagram.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_fixed_omega_phase_diagram.py",
    "outputs/jensen_window_pf_suzuki_jordan_totient_sign_scout.md",
    "work/rh_compute/results/jensen_window_pf_suzuki_jordan_totient_sign_scout.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_totient_sign_scout.py",
    "outputs/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.md",
    "work/rh_compute/results/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.py",
    "outputs/jensen_window_pf_suzuki_jordan_error_kernel_reduction.md",
    "work/rh_compute/results/jensen_window_pf_suzuki_jordan_error_kernel_reduction.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_error_kernel_reduction.py",
    "outputs/jensen_window_pf_suzuki_cofinal_l2_hierarchy.md",
    "work/rh_compute/results/jensen_window_pf_suzuki_cofinal_l2_hierarchy.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_cofinal_l2_hierarchy.py",
    "outputs/jensen_window_pf_jordan_muntz_causal_energy_bridge.md",
    "work/rh_compute/results/jensen_window_pf_jordan_muntz_causal_energy_bridge.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_causal_energy_bridge.py",
    "sup_N ||E_(omega_j,N)||_H < infinity",
    "outputs/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.md",
    "work/rh_compute/results/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.py",
    "(1-2omega)||t^(-omega)f_(2omega,N)||_2",
    "outputs/jensen_window_pf_burnol_cell_energy_tail_obstruction.md",
    "work/rh_compute/results/jensen_window_pf_burnol_cell_energy_tail_obstruction.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_burnol_cell_energy_tail_obstruction.py",
    "sup_N Q_(omega,N)<infinity",
    "N^(-1/2-omega)",
    "short-multiplicative-interval",
    "outputs/jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.md",
    "work/rh_compute/results/jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.py",
    "Q_(omega,infinity)<infinity",
    "V_(omega,N)(2^j*N)",
    "(k/(q+1),k/q]",
    "outputs/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.md",
    "work/rh_compute/results/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.py",
    "stationary Gram form",
    "signed Mobius off-diagonal",
    "Positive definiteness",
    "outputs/jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.md",
    "work/rh_compute/results/jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.py",
    "Ornstein-Uhlenbeck",
    "reciprocal-Mobius tail energy",
    "sign-indefinite",
    "cofinal sequence",
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
    "P_(K,infinity)(i,j)",
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
    "outputs/jensen_window_pf_mertens_planar_curvature_energy_scout.md",
    "work/rh_compute/results/jensen_window_pf_mertens_planar_curvature_energy_scout.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_curvature_energy_scout.py",
    "max `K*V_K` is",
    "non-falsification only",
    "outputs/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.md",
    "work/rh_compute/results/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.py",
    "rational inner fixed-shift quotient",
    "universal same-shift",
    "all-height natural-mollifier",
    "Brownian max-kernel total nonnegativity",
    "dilation orbit Bessel",
    "det(I+/-K_(omega_n,nu_n)[t])!=0 for every finite t",
    "distinct shifted zeros accumulating at itself",
    "RH iff D(omega_n) holds for some sequence omega_n->0",
    "x^(-1/2)*1_(1,infinity)(x)-h_omega^<1>(x) belongs to L2(1,infinity)",
    "all 6000 sampled",
    "arbitrary selected integers k_j>=1",
    "k=2",
    "positive Jordan-totient main term vanishes identically",
    "signed Mobius-error cancellation theorem",
    "H_(omega,k)(exp(t))-P_(omega,k)(t) in L2(0,infinity)",
    "finite Suzuki residual-energy cutoff",
    "generic backward Stieltjes-cone invariance is false",
    "outputs/jensen_window_pf_edrei_heat_flow_boundary_gate.md",
    "work/rh_compute/results/jensen_window_pf_edrei_heat_flow_boundary_gate.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_edrei_heat_flow_boundary_gate.py",
    "uses positive mixing to assert the polarized Phi Pick kernel is nonnegative",
    "promotes horizontal xi-modulus monotonicity to the hyperbolic Pick direction",
    "promotes total nonnegativity of a Hankel kernel without ruling out a finite +/-1 crossing",
    "uses |Theta|=1 on the real boundary to move a high contour across possible poles",
    "uses Hilbert-Schmidt domination as an all-t contraction certificate",
    "repeats Suzuki's M_v^2 exp(4vt)<1 estimate beyond its finite time ceiling",
    "demands one positive contraction gap uniform in t",
    "promotes Suzuki's local interval or a finite t grid to global Fredholm nonvanishing",
    "claims one fixed determinant family directly proves its terminal limit",
    "promotes finite positivity of Suzuki's summatory kernel to eventual sign or an L2 tail",
    "promotes finite positivity of a logarithmic smoothing to eventual sign",
    "uses the positive Jordan-totient residue main term as a dominant contribution after the Suzuki kernel annihilates it",
    "uses generic backward invariance of the Stieltjes moment cone",
    "Kill Gates",
    "The Jensen-window PF bridge is open",
    "former all-order signed-Hankel",
)


FORBIDDEN_STRINGS = (
    "therefore RH",
    "therefore `Lambda <= 0`",
    "we have proved Lambda <= 0",
    "the bridge is proved",
    "This proves the Jensen-window PF bridge",
)


def validate_note(path: Path) -> list[TargetIssue]:
    issues: list[TargetIssue] = []
    if not path.exists():
        return [TargetIssue("<file>", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    for required in REQUIRED_STRINGS:
        if required not in text:
            issues.append(TargetIssue("<required>", "missing-text", required))
    lowered = text.lower()
    for forbidden in FORBIDDEN_STRINGS:
        if forbidden.lower() in lowered:
            issues.append(TargetIssue("<forbidden>", "forbidden-text", forbidden))

    kill_gate_count = 0
    in_kill_gates = False
    for line in text.splitlines():
        if line.startswith("## Kill Gates"):
            in_kill_gates = True
            continue
        if in_kill_gates and line.startswith("## "):
            break
        if in_kill_gates and line.strip().endswith(";"):
            kill_gate_count += 1
    if kill_gate_count < 6:
        issues.append(TargetIssue("Kill Gates", "too-few-kill-gates", str(kill_gate_count)))

    required_refs = (
        "outputs/signed_hankel_jensen_bridge_target.md",
        "outputs/sign_regularity_theorem_fit_matrix.md",
        "work/rh_compute/scripts/check_signed_hankel_jensen_bridge_target.py",
        "work/rh_compute/scripts/check_sign_regularity_theorem_fit_matrix.py",
    )
    for ref in required_refs:
        if ref not in text:
            issues.append(TargetIssue("Integration Points", "missing-ref", ref))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues = validate_note(args.note)
    if args.json:
        print(json.dumps({"ok": not issues, "issues": [asdict(issue) for issue in issues]}, indent=2, sort_keys=True))
    else:
        for issue in issues:
            print(f"JENSEN-WINDOW-PF-TARGET {issue.section} [{issue.issue}] {issue.detail}")
        print(f"validated Jensen-window PF bridge target with {len(issues)} issues")
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
