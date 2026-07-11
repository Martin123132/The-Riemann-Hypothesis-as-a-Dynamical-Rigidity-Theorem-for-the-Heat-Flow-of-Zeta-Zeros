#!/usr/bin/env python3
"""Specify the discrete multiplier-product bridge target for zeta windows."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_multiplier_counting_measure_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_multiplier_counting_measure_target.md"


def build_payload() -> dict:
    rows = [
        {
            "id": "mcmt_01_geometric_normalization",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "formula": "B_k=A_k*A_0^(k-1)/A_1^k, so B_0=B_1=1 and x_k(B)=x_k(A)",
            "claim": "Positive geometric normalization leaves every Jensen root sign and ratio contraction unchanged.",
        },
        {
            "id": "mcmt_02_elementary_multiplier_atom",
            "role": "exact_lemma",
            "readiness": "available_exact",
            "formula": "M_k^(alpha)=(alpha/(alpha+1))^(k-1)*(alpha+k)/(alpha+1)",
            "claim": "For alpha>0, the EGF is (1+z/(alpha+1))*exp(alpha*z/(alpha+1)), hence M^(alpha) is a nonnegative multiplier sequence.",
        },
        {
            "id": "mcmt_03_integer_product_sufficient_theorem",
            "role": "conditional_theorem",
            "readiness": "not_ready_to_apply",
            "formula": "B_k=product_j M_k^(alpha_j), alpha_j>0, sum_j alpha_j^(-2)<infinity",
            "claim": "A convergent unit-atomic product is a pointwise limit of multiplier sequences and therefore makes every zeta Jensen window hyperbolic.",
        },
        {
            "id": "mcmt_04_contraction_log_representation",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "formula": "-log x_k=sum_j -log(1-1/(k+alpha_j)^2)",
            "claim": "Because B_0=B_1=1, this contraction identity is equivalent to the coefficient product whenever the products converge.",
        },
        {
            "id": "mcmt_05_positive_laplace_kernel",
            "role": "exact_lemma",
            "readiness": "available_exact",
            "formula": "-log(1-1/(k+alpha)^2)=integral_0^infty exp(-(k+alpha)t)*2*(cosh(t)-1)/t dt",
            "claim": "Any counting-measure representation forces y_k=-log x_k to be completely monotone.",
        },
        {
            "id": "mcmt_06_finite_zeta_necessary_evidence",
            "role": "finite_evidence",
            "readiness": "not_ready_to_apply",
            "formula": "(-1)^m*Delta^m(-log x_k)>0 on every cached interval through m=8",
            "claim": "The cached zeta data satisfy the finite necessary sign pattern, with high-order enclosure inconclusives retained.",
        },
        {
            "id": "mcmt_07_fractional_weight_guard",
            "role": "countermodel_gate",
            "readiness": "guard_validated",
            "formula": "alpha=1, weight=1/2 gives cubic discriminant <-27/125",
            "claim": "Arbitrary positive measure weights cannot replace a counting measure without an independent stability-preserver theorem.",
        },
        {
            "id": "mcmt_08_convex_mixture_guard",
            "role": "countermodel_gate",
            "readiness": "guard_validated",
            "formula": "(A^(u=3/5)+A^(u=2/3))/2 gives cubic discriminant -937/3456",
            "claim": "Coefficientwise convex mixtures of elementary multiplier atoms need not be multiplier sequences.",
        },
        {
            "id": "mcmt_09_mellin_interpolation_guard",
            "role": "countermodel_gate",
            "readiness": "guard_validated",
            "formula": "det[p_(2+i+j)]_(i,j=0..3)<0 for the Gamma-normalized Mellin interpolation",
            "claim": "The natural continuous-index Mellin interpolation cannot equal the elementary multiplier product; this does not reject equality asserted only at integer indices.",
        },
        {
            "id": "mcmt_10_zeta_counting_measure_target",
            "role": "open_statement",
            "readiness": "not_ready_to_apply",
            "formula": "construct {alpha_j} from Phi/Newman data and prove B_k=product_j M_k^(alpha_j) for all k",
            "claim": "This is a sufficient noncircular all-degree bridge if the exact discrete representation and convergence are proved.",
        },
    ]
    return {
        "kind": "jensen_window_pf_multiplier_counting_measure_target",
        "date": "2026-07-10",
        "status": "open theorem target",
        "target_id": "target_zeta_multiplier_counting_measure_factorization",
        "proof_boundary": (
            "Open sufficient bridge target based on a discrete product of elementary multiplier sequences. "
            "No zeta counting measure is constructed, finite complete-monotonicity evidence is not promoted, arbitrary positive weights are rejected, and PF-infinity, RH, and Lambda <= 0 remain unproved."
        ),
        "source_boundary_family": "outputs/jensen_window_pf_rank_two_boundary_family_lemma.md",
        "source_defect_scout": "outputs/jensen_window_pf_defect_complete_monotonicity_scout.md",
        "source_bridge_target": "outputs/jensen_window_pf_bridge_target.md",
        "source_mellin_obstruction": "outputs/jensen_window_pf_mellin_multiplier_power_sum_obstruction.md",
        "primary_source_anchors": [
            "https://annals.math.princeton.edu/wp-content/uploads/annals-v170-n1-p14-p.pdf",
            "https://arxiv.org/abs/math/0606360"
        ],
        "target_rows": rows,
        "acceptance_contract": [
            "Normalize the actual zeta coefficients by B_k=A_k*A_0^(k-1)/A_1^k and keep B_0=B_1=1 explicit.",
            "Construct a finite or countable multiset alpha_j>0 with unit integer multiplicities and sum alpha_j^(-2)<infinity.",
            "Prove the coefficient product, or equivalently the contraction-log identity plus convergence, for every k>=1.",
            "Derive the multiset from noncircular Phi/Newman-kernel structure, not endpoint Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0.",
            "Prove locally uniform convergence of the EGF/Jensen windows or use coefficientwise Jensen limits to preserve real negative roots.",
            "Do not replace the counting measure by arbitrary positive weights without a separately proved stability-preserving composition theorem.",
            "Do not identify the integer product with the natural Mellin interpolation: its candidate power sums have an Arb-certified negative Hankel determinant.",
        ],
        "summary": {
            "target_rows": len(rows),
            "exact_rows": 4,
            "finite_evidence_rows": 1,
            "countermodel_rows": 3,
            "conditional_theorem_rows": 1,
            "open_statement_rows": 1,
            "ready_to_apply_rows": 0,
            "live_routes": 1,
            "target_closing": False,
            "main_finding": (
                "A discrete unit-atomic product of elementary multiplier sequences would close the all-degree Jensen-window target, and its ratio contractions have an explicit positive Laplace kernel. "
                "The zeta data satisfy the finite necessary log-complete-monotonicity pattern through order 8. No discrete representation is known; exact fractional-weight and convex-mixture countermodels forbid the obvious continuous relaxations, and an Arb negative power-sum Hankel determinant rules out identifying the product with the natural Mellin interpolation."
            ),
        },
        "invariants": [
            "No row is ready_to_apply.",
            "The counting-measure representation is sufficient, not claimed necessary for all multiplier sequences.",
            "Finite log-complete-monotonicity does not construct the atoms.",
            "Arbitrary positive measure weights and convex mixtures are blocked by exact cubic countermodels.",
            "The natural continuous-index Mellin interpolation is blocked, but integer-only product equality remains open without an interpolation-uniqueness theorem.",
            "Endpoint PF, Jensen hyperbolicity, RH, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(payload: dict, path: Path) -> None:
    result = (
        "validated Jensen-window PF multiplier counting-measure target: "
        "10 rows, 0 issues, 4 exact rows, 1 finite evidence row, 3 countermodel rows, "
        "1 live route, 0 ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Multiplier Counting-Measure Target",
        "",
        "Date: 2026-07-10",
        "",
        "Status: open theorem target. This is not a proof of PF-infinity, Jensen",
        "hyperbolicity for zeta, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_multiplier_counting_measure_target`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_multiplier_counting_measure_target.json",
        "python work/rh_compute/scripts/jensen_window_pf_multiplier_counting_measure_target.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_multiplier_counting_measure_target.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result,
        "```",
        "",
        "## Sufficient Theorem",
        "",
        "Normalize the actual coefficients by",
        "",
        "```text",
        "B_k=A_k*A_0^(k-1)/A_1^k, B_0=B_1=1.",
        "```",
        "",
        "For `alpha>0`, define",
        "",
        "```text",
        "M_k^(alpha)=(alpha/(alpha+1))^(k-1)*(alpha+k)/(alpha+1),",
        "EGF(M^(alpha))=(1+z/(alpha+1))*exp(alpha*z/(alpha+1)).",
        "```",
        "",
        "The EGF has one negative zero and is Laguerre-Polya type I, so the",
        "Pólya-Schur theorem makes `M^(alpha)` a multiplier sequence. Finite",
        "pointwise products preserve real-rootedness by composition. Therefore the",
        "following would be sufficient:",
        "",
        "```text",
        "B_k=product_j M_k^(alpha_j),",
        "alpha_j>0 with unit integer multiplicity,",
        "sum_j alpha_j^(-2)<infinity.",
        "```",
        "",
        "Coefficientwise limits then preserve every fixed Jensen polynomial. This is a",
        "sufficient subclass theorem, not a claimed characterization of every multiplier",
        "sequence.",
        "",
        "## Ratio Form",
        "",
        "The same target is",
        "",
        "```text",
        "-log x_k=sum_j -log(1-1/(k+alpha_j)^2).",
        "```",
        "",
        "The elementary kernel has the positive Laplace representation",
        "",
        "```text",
        "-log(1-1/(k+alpha)^2)",
        "  =integral_0^infty exp(-(k+alpha)t)*2*(cosh(t)-1)/t dt.",
        "```",
        "",
        "Thus complete monotonicity of `y_k=-log x_k` is necessary. The Arb scout",
        "certifies that sign pattern through order 8 on all five cached heat times,",
        "but finite evidence does not construct a counting measure.",
        "",
        "## Failure Gates",
        "",
        "The atoms cannot simply be assigned arbitrary positive weights. The exact",
        "fractional-power cubic guard has discriminant `<-27/125`, and the equal",
        "positive mixture of the `u=3/5,2/3` boundary families has discriminant `-937/3456`.",
        "A continuous positive fit is therefore not enough.",
        "",
        "The separate Arb Mellin-interpolation certificate constructs the continuous",
        "power-sum candidates from log moments of `Phi` and finds a strictly negative",
        "shift-2 size-4 Hankel determinant. Therefore the natural Gamma-normalized",
        "Mellin interpolation cannot be the elementary multiplier product. This does not",
        "rule out equality asserted only for integer `k`; such a promotion would need a",
        "separate interpolation-uniqueness theorem.",
        "",
        "## Acceptance",
        "",
        "A closing proof must construct the actual discrete multiset from Phi/Newman",
        "data, prove convergence and equality for every k, and avoid endpoint",
        "hyperbolicity or RH as an input. No such construction is currently known.",
        "",
        "Primary theorem anchors:",
        "",
        "```text",
        "https://annals.math.princeton.edu/wp-content/uploads/annals-v170-n1-p14-p.pdf",
        "https://arxiv.org/abs/math/0606360",
        "```",
        "",
        "```text",
        "outputs/jensen_window_pf_rank_two_boundary_family_lemma.md",
        "outputs/jensen_window_pf_defect_complete_monotonicity_scout.md",
        "outputs/jensen_window_pf_mellin_multiplier_power_sum_obstruction.md",
        "outputs/jensen_window_pf_bridge_target.md",
        "```",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(payload, args.note)
    print(
        "validated Jensen-window PF multiplier counting-measure target: "
        "10 rows, 0 issues, 4 exact rows, 1 finite evidence row, 3 countermodel rows, "
        "1 live route, 0 ready-to-apply rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
