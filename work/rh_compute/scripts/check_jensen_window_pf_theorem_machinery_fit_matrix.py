#!/usr/bin/env python3
"""Validate the Jensen-window PF theorem-machinery fit matrix."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MATRIX = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_theorem_machinery_fit_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_theorem_machinery_fit_matrix.md"


ALLOWED_VERDICTS = {
    "endpoint_equivalence_only",
    "possible_if_new_kernel_representation",
    "possible_if_preserver_hypotheses_are_proved",
    "route_mismatch",
    "conditional_downstream_only",
    "rejected_circular",
}

REQUIRED_IDS = {
    "tm_01_asw_edrei_pf_sequence_characterization",
    "tm_08_sokal_log_derivative_stieltjes_criterion",
    "tm_09_krein_stieltjes_pick_phi_kernel",
    "tm_10_edrei_radial_heat_backward_boundary",
    "tm_11_suzuki_arithmetic_hankel_canonical_system",
    "tm_02_schoenberg_pf_functions_variation_diminishing",
    "tm_03_karlin_basic_composition_cauchy_binet",
    "tm_04_polya_schur_multiplier_preservers",
    "tm_05_gantmacher_krein_sign_regular_matrices",
    "tm_06_laguerre_polya_jensen_limit",
    "tm_07_finite_grid_or_rh_assuming_shortcuts",
}

NON_CLOSING_VERDICTS = {
    "endpoint_equivalence_only",
    "route_mismatch",
    "conditional_downstream_only",
    "rejected_circular",
}

REQUIRED_FEATURES = {
    "outputs every binomially weighted Jensen-window Toeplitz matrix for all d,n",
    "handles binomial weights binom(d,j)",
    "handles all shifts n uniformly",
    "starts from proved noncircular hypotheses about the actual A_k(0)",
    "does not assume Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Theorem Machinery Fit Matrix",
    "Status: theorem-search fit matrix",
    "This is not a proof of Jensen-window",
    "jwpf_06_sign_regular_to_jensen_pf_conversion",
    "work/rh_compute/results/jensen_window_pf_theorem_machinery_fit_matrix.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_theorem_machinery_fit_matrix.py",
    "validated Jensen-window PF theorem machinery fit matrix: 11 rows, 0 issues, 0 ready-to-apply rows",
    "tm_01_asw_edrei_pf_sequence_characterization",
    "tm_08_sokal_log_derivative_stieltjes_criterion",
    "tm_09_krein_stieltjes_pick_phi_kernel",
    "tm_10_edrei_radial_heat_backward_boundary",
    "outputs/jensen_window_pf_edrei_hankel_boundary_flux_gate.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_edrei_hankel_boundary_flux_gate.py",
    "outputs/jensen_window_pf_edrei_raw_moment_collision_resolution_gate.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_edrei_raw_moment_collision_resolution_gate.py",
    "partial_lambda log D_(r,s)",
    "subexponential",
    "gap-adapted",
    "tm_11_suzuki_arithmetic_hankel_canonical_system",
    "outputs/jensen_window_pf_suzuki_determinant_only_reduction.md",
    "outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md",
    "outputs/jensen_window_pf_suzuki_jordan_totient_sign_scout.md",
    "outputs/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.md",
    "outputs/jensen_window_pf_suzuki_jordan_error_kernel_reduction.md",
    "outputs/jensen_window_pf_suzuki_cofinal_l2_hierarchy.md",
    "outputs/jensen_window_pf_jordan_muntz_causal_energy_bridge.md",
    "outputs/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.md",
    "outputs/jensen_window_pf_burnol_cell_energy_tail_obstruction.md",
    "outputs/jensen_window_pf_burnol_tail_discrepancy_dyadic_reduction.md",
    "outputs/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.md",
    "H_(omega,k)(exp(t))-P_(omega,k)(t)",
    "sup_N ||E_(omega,N)||_H<infinity",
    "(1-2omega)||t^(-omega)f_(2omega,N)||_2",
    "Q_(omega,N)",
    "N^(-1/2-omega)",
    "short-multiplicative",
    "Q_(omega,infinity)<infinity",
    "V_(omega,N)(2^j*N)",
    "(k/(q+1),k/q]",
    "stationary Gram",
    "signed Mobius off-diagonal",
    "bounded Gram diagonal",
    "outputs/jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.md",
    "Ornstein-Uhlenbeck",
    "reciprocal-Mobius tail energy",
    "sign-indefinite remainder",
    "N^(1-alpha)",
    "outputs/jensen_window_pf_ou_mertens_mean_square_reduction.md",
    "dyadic Mertens mean square",
    "origin-anchored",
    "averaged-Chowla",
    "outputs/jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.md",
    "constant anchor mode",
    "B_1(X)-TI_X+TII_X",
    "-TI_X+TII_X",
    "endogenous",
    "outputs/jensen_window_pf_mertens_anchor_logarithmic_tail_energy_reduction.md",
    "sum_N P_(alpha/2,N)/N",
    "one logarithm short",
    "outputs/jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.md",
    "ceil(sqrt(K))",
    "Davenport's arbitrary logarithmic",
    "outputs/jensen_window_pf_mertens_local_path_cosine_transference.md",
    "sqrt(K)r_K-tau_(K,0)",
    "zero dyadic endpoint tails",
    "outputs/jensen_window_pf_mertens_centered_bridge_vaughan_handoff.md",
    "min(i,j)-ij/K",
    "C_K=-TI_K+TII_K",
    "outputs/jensen_window_pf_mertens_affine_tent_bridge_handoff.md",
    "u_K-(1/2)u_(2K)",
    "O_(loc,K)=-TI_(loc,K)+TII_(loc,K)",
    "H_K=O(K^(2+epsilon))",
    "outputs/jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.md",
    "B_K=-V_(I,K)+V_(II,K)",
    "outputs/jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.md",
    "X_(K,r)=-T_(I,K,r)+T_(II,K,r)",
    "outputs/jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.md",
    "R_(alpha,K)=ceil(K^(1-alpha/2))",
    "2K^(-1)+8K^(-alpha/2)",
    "outputs/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.md",
    "P_(K,infinity)(i,j)",
    "7pi^2/24",
    "-2<u_I,G u_II>",
    "outputs/jensen_window_pf_mertens_shift_kernel_variation_handoff.md",
    "Var_(n in I_(K,h))W_(K,h)(n)<6rho_K",
    "joint cancellation across shifts",
    "outputs/jensen_window_pf_mertens_planar_abel_handoff.md",
    "3*pi^2*(1+log(2R))/K",
    "E_(alpha,K)=O_epsilon(K^(1+epsilon))",
    "outputs/jensen_window_pf_mertens_planar_edge_gram_vaughan_handoff.md",
    "(27/2)*pi^2*K*(1+log(2R))",
    "signed edge-pair",
    "outputs/jensen_window_pf_mertens_planar_incidence_antidiagonal_handoff.md",
    "C_3=210, C_4=-210, C=0",
    "O_(CC,int,K)",
    "Y_(K,r)",
    "outputs/jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.md",
    "O=O_(CC,int)+sum_(r<=R)Z_(K,r)/q_r^2",
    "phi_r(K-1)U_r-e_r V_r",
    "E_perp+c(beta+gamma/c)^2",
    "joined `Z_(K,r)`",
    "outputs/jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.md",
    "R_B=(E_B-D_B)/2+C_delta",
    "D_B<5pi^2/12",
    "|R_B-E_B/2|=O(1+log R)",
    "E_B<=pi(1/2+2pi)K^(-2)H_K",
    "outputs/jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.md",
    "O(K^(alpha/2))",
    "K^2 E_B/H_K->0",
    "RH-equivalent",
    "O_(CC,int)+E_B/2",
    "outputs/jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.md",
    "J_B=O_(CC,int)+E_B/2=L/2+(D_B-D_L)/2-C_delta",
    "join returns",
    "arXiv:2607.15574",
    "Lewko-Lewko",
    "iterated-log",
    "outputs/jensen_window_pf_mertens_planar_curvature_energy_scout.md",
    "K*V_K<6.571",
    "non-falsification only",
    "outputs/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.md",
    "fixed-shift rational innerness",
    "reciprocal square",
    "Burnol's weighted natural",
    "weighted all-height",
    "Brownian max kernel",
    "orbit Bessel",
    "shifted-zero accumulation",
    "tm_02_schoenberg_pf_functions_variation_diminishing",
    "tm_03_karlin_basic_composition_cauchy_binet",
    "tm_04_polya_schur_multiplier_preservers",
    "tm_05_gantmacher_krein_sign_regular_matrices",
    "tm_06_laguerre_polya_jensen_limit",
    "tm_07_finite_grid_or_rh_assuming_shortcuts",
    "no `ready_to_apply` row",
    "degree-3 countermodel",
    "Structural Ansatz Workbench",
    "outputs/jensen_window_pf_structural_ansatz_matrix.md",
    "work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py",
    "validated Jensen-window PF structural ansatz matrix: 6 ansatz rows, 0 issues, 0 ready-to-apply rows",
    "outputs/jensen_window_pf_schur_shape_contract.md",
    "work/rh_compute/results/jensen_window_pf_schur_shape_contract.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_schur_shape_contract.py",
    "finite-band shape obligations",
    "outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md",
    "work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py",
    "`15` formula rows with nonnegative Bernstein coefficients",
    "`0` kernel identities found",
    "outputs/jensen_window_pf_edrei_stieltjes_equivalence_gate.md",
    "H_0'/H_0",
    "outputs/jensen_window_pf_phi_pick_kernel_target.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_phi_pick_kernel_target.py",
    "outputs/jensen_window_pf_edrei_heat_flow_boundary_gate.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_edrei_heat_flow_boundary_gate.py",
    "outputs/jensen_window_pf_xi_pick_suzuki_hankel_bridge.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_xi_pick_suzuki_hankel_bridge.py",
    "det(I+/-K",
    "terminal premise is redundant",
    "L2 Jordan-totient summatory",
    "6000-point positive",
    "logarithmic smoothing order",
    "k=2",
    "positive residue main term cancels",
    "signed Mobius-error convolution",
    "P_Phi(z)",
)


@dataclass(frozen=True)
class MatrixIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_source(row_id: str, source: str) -> list[MatrixIssue]:
    if source.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / source).exists():
        return []
    return [MatrixIssue(row_id, "missing-source", source)]


def validate_row(row: dict) -> list[MatrixIssue]:
    issues: list[MatrixIssue] = []
    row_id = str(row.get("id", "<missing-id>"))
    required_fields = (
        "id",
        "theorem_family",
        "primary_sources",
        "fit_to_jwpf06",
        "required_hypotheses",
        "current_evidence",
        "fatal_gap",
        "next_action",
        "verdict",
        "hypotheses_verified",
        "would_close_jwpf06_if_verified",
        "proof_boundary",
    )
    for key in required_fields:
        if key not in row:
            issues.append(MatrixIssue(row_id, "missing-field", key))

    verdict = row.get("verdict")
    if verdict not in ALLOWED_VERDICTS:
        issues.append(MatrixIssue(row_id, "bad-verdict", repr(verdict)))
    if row.get("hypotheses_verified") is not False:
        issues.append(MatrixIssue(row_id, "hypotheses-should-not-be-verified", repr(row.get("hypotheses_verified"))))
    if verdict in NON_CLOSING_VERDICTS and row.get("would_close_jwpf06_if_verified") is True:
        issues.append(MatrixIssue(row_id, "nonclosing-verdict-closes", str(verdict)))
    if row.get("would_close_jwpf06_if_verified") is True and verdict not in {
        "possible_if_new_kernel_representation",
        "possible_if_preserver_hypotheses_are_proved",
    }:
        issues.append(MatrixIssue(row_id, "bad-closing-verdict", str(verdict)))

    for key in ("primary_sources", "required_hypotheses", "current_evidence"):
        value = row.get(key)
        if not isinstance(value, list) or not value:
            issues.append(MatrixIssue(row_id, f"bad-{key}", repr(value)))
    for source in row.get("primary_sources", []):
        if isinstance(source, str):
            issues.extend(validate_source(row_id, source))
        else:
            issues.append(MatrixIssue(row_id, "bad-source", repr(source)))

    if not str(row.get("fatal_gap", "")).strip():
        issues.append(MatrixIssue(row_id, "missing-fatal-gap", row_id))
    if not str(row.get("next_action", "")).strip():
        issues.append(MatrixIssue(row_id, "missing-next-action", row_id))
    boundary = str(row.get("proof_boundary", "")).lower()
    if not any(marker in boundary for marker in ("not ", "only", "missing", "rejected", "conditional", "candidate")):
        issues.append(MatrixIssue(row_id, "weak-proof-boundary", row.get("proof_boundary", "")))

    combined = " ".join(str(row.get(key, "")) for key in ("fit_to_jwpf06", "fatal_gap", "proof_boundary")).lower()
    for forbidden in ("therefore rh", "proves lambda <= 0", "bridge is proved", "ready_to_apply"):
        if forbidden in combined:
            issues.append(MatrixIssue(row_id, "forbidden-overclaim", forbidden))
    return issues


def validate_note(path: Path) -> list[MatrixIssue]:
    if not path.exists():
        return [MatrixIssue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[MatrixIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(MatrixIssue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "the bridge is proved"):
        if forbidden in lowered:
            issues.append(MatrixIssue("note", "forbidden-text", forbidden))
    return issues


def validate(matrix_path: Path, note_path: Path) -> tuple[list[MatrixIssue], int, int]:
    matrix = load_json(matrix_path)
    issues: list[MatrixIssue] = []
    if matrix.get("kind") != "jensen_window_pf_theorem_machinery_fit_matrix":
        issues.append(MatrixIssue("<matrix>", "bad-kind", repr(matrix.get("kind"))))
    if matrix.get("target_obligation") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(MatrixIssue("<matrix>", "bad-target-obligation", repr(matrix.get("target_obligation"))))

    boundary = str(matrix.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(MatrixIssue("<matrix>", "weak-proof-boundary", matrix.get("proof_boundary", "")))

    features = set(matrix.get("required_bridge_features", []))
    missing_features = REQUIRED_FEATURES - features
    if missing_features:
        issues.append(MatrixIssue("<matrix>", "missing-required-features", repr(sorted(missing_features))))

    allowed = set(matrix.get("allowed_verdicts", []))
    if not ALLOWED_VERDICTS.issubset(allowed):
        issues.append(MatrixIssue("<matrix>", "missing-allowed-verdicts", repr(sorted(ALLOWED_VERDICTS - allowed))))

    rows = matrix.get("rows", [])
    if not isinstance(rows, list) or not rows:
        issues.append(MatrixIssue("<matrix>", "missing-rows", "rows must be a nonempty list"))
        rows = []

    seen: set[str] = set()
    ready_count = 0
    closing_candidate_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(MatrixIssue("<matrix>", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        if row_id in seen:
            issues.append(MatrixIssue(row_id, "duplicate-id", row_id))
        seen.add(row_id)
        if row.get("verdict") == "ready_to_apply":
            ready_count += 1
        if row.get("would_close_jwpf06_if_verified") is True:
            closing_candidate_count += 1
        issues.extend(validate_row(row))

    for missing in sorted(REQUIRED_IDS - seen):
        issues.append(MatrixIssue(missing, "missing-required-row", missing))
    if ready_count != 0:
        issues.append(MatrixIssue("<matrix>", "ready-to-apply-row-present", str(ready_count)))
    if closing_candidate_count < 2:
        issues.append(MatrixIssue("<matrix>", "too-few-structural-candidates", str(closing_candidate_count)))
    if "tm_07_finite_grid_or_rh_assuming_shortcuts" in seen:
        shortcut = next(row for row in rows if isinstance(row, dict) and row.get("id") == "tm_07_finite_grid_or_rh_assuming_shortcuts")
        if shortcut.get("verdict") != "rejected_circular":
            issues.append(MatrixIssue("tm_07_finite_grid_or_rh_assuming_shortcuts", "shortcut-not-rejected", str(shortcut.get("verdict"))))

    issues.extend(validate_note(note_path))
    return issues, len(rows), ready_count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, row_count, ready_count = validate(args.matrix, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "rows": row_count,
                    "ready_to_apply_rows": ready_count,
                    "issues": [asdict(issue) for issue in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for issue in issues:
            print(f"JWPF-THEOREM-MACHINERY {issue.row_id} [{issue.issue}] {issue.detail}")
        print(
            "validated Jensen-window PF theorem machinery fit matrix: "
            f"{row_count} rows, {len(issues)} issues, {ready_count} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
