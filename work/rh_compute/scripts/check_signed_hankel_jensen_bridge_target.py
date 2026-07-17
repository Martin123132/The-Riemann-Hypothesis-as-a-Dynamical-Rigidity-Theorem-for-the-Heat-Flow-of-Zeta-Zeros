#!/usr/bin/env python3
"""Validate the signed-Hankel/Jensen bridge target specification."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_NOTE = REPO_ROOT / "outputs/signed_hankel_jensen_bridge_target.md"


@dataclass(frozen=True)
class TargetIssue:
    section: str
    issue: str
    detail: str


REQUIRED_STRINGS = (
    "Status: theorem target",
    "This is not a proof of RH or `Lambda <= 0`",
    "for every d >= 1 and n >= 0",
    "outputs/jensen_window_pf_bridge_target.md",
    "outputs/jensen_window_pf_obligation_algebra.md",
    "outputs/arb_jensen_window_pf_obligation_diagnostic.md",
    "outputs/arb_jensen_window_sturm_hyperbolicity_diagnostic.md",
    "python work/rh_compute/scripts/check_jensen_window_pf_bridge_target.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_obligation_algebra.py",
    "python work/rh_compute/scripts/check_arb_jensen_window_pf_obligation_manifest.py",
    "python work/rh_compute/scripts/check_arb_jensen_window_sturm_manifest.py",
    "B^{d,n,0}_j = binom(d,j) A_{n+j}(0)",
    "finite PF-infinity sequence",
    "degree 3 and degree 4 introduce additional",
    "1470/1470",
    "210/210",
    "105/105",
    "positive-root counts for `Q_{d,n,lambda}(y)=P_{d,n,lambda}(-y)`",
    "R_{k,n}(j_1,...,j_k)",
    "(-1)^(k(k-1)/2)",
    "n = 0",
    "k = 2..7",
    "689,795/689,795 finite minors positive",
    "n = 0..20",
    "k = 2..5",
    "1,322,685/1,322,685 finite minors positive",
    "k = 6",
    "840,840/840,840 finite minors positive",
    "k = 7",
    "675,675/675,675 finite minors positive",
    "k = 8",
    "315,315/315,315 finite minors positive",
    "3,154,515/3,154,515 finite shifted minors positive",
    "All certificates remain finite",
    "Candidate Theorem B-Star",
    "Degree 2 is exact",
    "outputs/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.md",
    "Karlin's 1964 reciprocal-gamma theorem",
    "sign-changing integrand",
    "CV_n^2<=2/(2n+1)",
    "higher compound concentration",
    "outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md",
    "C_n=q_(n+1)*q_(n+3)-q_(n+2)^2+1>0",
    "unit-buffered log-concavity",
    "wrong positive",
    "outputs/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.md",
    "s_319>251/500",
    "C_n(-100)>57613471/66107054971>0, n>=318",
    "outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md",
    "C_n'/r_(n+2)=alpha_n*C_(n+1)+beta_n*C_n, alpha_n>0",
    "C_n(lambda)>0 for every n>=0 and finite lambda>=-100",
    "D_(3,n)(0)<0 for every n>=0",
    "outputs/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.md",
    "R_(3,n)(j_1,j_2,j_3)<0",
    "reshaped-Hankel orders two and three are complete",
    "Compound order four is",
    "outputs/jensen_window_pf_compound_order4_condensation_gate.md",
    "G_(n+1)^2>x_(n+3)^3*G_n*G_(n+2)",
    "all 317 available margins",
    "P_n<=4/(n+3)^2, n>=317",
    "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
    "K_1(t)=(log g_1(t))''<=7/(2t^2), t>=319",
    "J_1(t)>=1/(7t)",
    "|P_n-P_n^(1)|<=2/(5k^2)",
    "s=(3/10,2/5,41/100,49/100,11/20)",
    "H_(4,0)<0",
    "order four is a genuinely new compound inequality",
    "outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md",
    "H_(4,n)(lambda)>0 for every integer n>=0 and every lambda in [-100,0]",
    "outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md",
    "R_(4,n)(j_1,j_2,j_3,j_4)>0",
    "outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md",
    "294912*G_2^10",
    "a_n=(4*n+34)*H_(4,n)/H_(4,n+1)>0",
    "target_compound_order5_m100_entry",
    "H_(5,n)(-100)>0 for every n>=0",
    "outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md",
    "H_(5,n)(-100)=W_n*J_n(-100)>0 for every 0<=n<=316",
    "relative_316>0.006269",
    "analytic tail `n>=317` remained",
    "outputs/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.md",
    "C_n=Delta^2 log(F_n)-Delta^2 log(d_(n+3))",
    "J_n>0 iff C_n<-4log(x_k)",
    "C_n<=100/(n+4)^2",
    "outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md",
    "|C_n-C_n^(1)|<=37/k^2",
    "q_1''(t)<=60/t^2 for every real t>=320",
    "outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md",
    "outputs/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.md",
    "outputs/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.md",
    "t^2*q_1''(t)<10<60 for every mode u>=20",
    "outputs/jensen_window_pf_compound_order5_m100_entry_certificate.md",
    "discharges `target_compound_order5_m100_entry`",
    "outputs/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.md",
    "H_(5,n)(lambda)>0 for every n>=0 and -100<=lambda<=0",
    "R_(5,n)(j_1,j_2,j_3,j_4,j_5;lambda)>0",
    "outputs/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.md",
    "outputs/jensen_window_pf_compound_order6_m100_prefix_certificate.md",
    "H_(5,n)=A_(n+4)^5*exp(p(n+4))",
    "P_k<=320/k^2",
    "p_1''(t)<=200/t^2 for every real t>=321",
    "outputs/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.md",
    "outputs/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.md",
    "outputs/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.md",
    "t^2*p_1''(t)<22.769<200 for every mode u>=20",
    "outputs/jensen_window_pf_compound_order6_m100_entry_certificate.md",
    "Q_(6,n)(-100)=-H_(6,n)(-100)>0 for every n>=0",
    "outputs/jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate.md",
    "epsilon_6*R_(6,n)(j_1,...,j_6;lambda)>0",
    "outputs/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.md",
    "for every fixed m there exists N_m such that Q_(m,n)(lambda)>0 for n>=N_m",
    "outputs/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.md",
    "52183852646400*G_2^21*h^21",
    "Q_(7,n)*H_(5,n+2)=Q_(6,n+1)^2-Q_(6,n)*Q_(6,n+2)",
    "outputs/jensen_window_pf_compound_order7_m100_prefix_certificate.md",
    "Q_(7,n)(-100)>0 for every 0<=n<=314",
    "L_314>9/1000",
    "outputs/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.md",
    "R_k<=900/k^2 for every k>=321",
    "outputs/jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension.md",
    "0<=delta_k<2/k^8 for every k>=300",
    "outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md",
    "min(T_j,T_j^(1))>=3/(2*j), j>=320",
    "|R_k-R_k^(1)|<262/k^2, k>=321",
    "r_1''(t)<=600/t^2 => R_k<863/k^2<900/k^2",
    "continuous fourth-nested curvature theorem",
    "All-shift order seven remains open",
    "The historical order-seven frontier above is now superseded",
    "outputs/jensen_window_pf_all_order_endpoint_heat_reduction.md",
    "outputs/jensen_window_pf_endpoint_deep_schur_coordinate.md",
    "Q_(m,n)(-100)=A_0(-100)^m s_((n+m-1)^m)(h)",
    "s_((N^m))(h)>0 for every m>=10 and N>=m-1",
    "outputs/jensen_window_pf_endpoint_order10_counterexample.md",
    "s_((N^10))(h)<0 for N=9,10,11,12",
    "The four negative shapes lie inside the required deep cone",
    "all-order endpoint antecedent cannot be supplied",
    "s_(1,1,1)(h)<0",
    "outputs/jensen_window_pf_deep_schur_toda_boundary_gate.md",
    "tau_(m+1,N) tau_(m-1,N)",
    "h_0 h_2/h_1^2>0",
    "H(z)=exp(z/100)/((1-z)(1-2z))",
    "-222484532394597/2000000000000<0",
    "strict full unweighted",
    "weaker properties specific to the actual",
    "Assumption 2 is false",
    "weaker than all-shift signed-Hankel positivity",
    "former all-order signed-Hankel route to",
    "closed by counterexample",
    "Degree 3 is the first nontrivial bridge obstruction",
    "low-order finite sign checks cannot be promoted into Jensen hyperbolicity",
    "Proof Obligations",
    "Kill Gates",
    "This theorem target is open",
)


FORBIDDEN_STRINGS = (
    "therefore RH",
    "therefore `Lambda <= 0`",
    "we have proved Lambda <= 0",
    "the bridge is proved",
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

    proof_obligation_count = sum(
        1 for line in text.splitlines() if line.startswith(("1. ", "2. ", "3. ", "4. ", "5. ", "6. "))
    )
    if proof_obligation_count < 6:
        issues.append(TargetIssue("Proof Obligations", "too-few-obligations", str(proof_obligation_count)))

    if "outputs/jensen_hankel_bridge_algebra.md" not in text:
        issues.append(TargetIssue("Degree 3", "missing-algebra-note-ref", "outputs/jensen_hankel_bridge_algebra.md"))
    if "work/rh_compute/results/jensen_hankel_bridge_algebra.json" not in text:
        issues.append(TargetIssue("Degree 3", "missing-algebra-json-ref", "work/rh_compute/results/jensen_hankel_bridge_algebra.json"))
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
            print(f"SIGNED-HANKEL-JENSEN-TARGET {issue.section} [{issue.issue}] {issue.detail}")
        print(f"validated signed-Hankel/Jensen bridge target specification with {len(issues)} issues")
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
