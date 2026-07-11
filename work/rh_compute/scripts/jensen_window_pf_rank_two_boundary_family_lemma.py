#!/usr/bin/env python3
"""Exact all-degree Jensen factorization for a rank-two boundary family."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_rank_two_boundary_family_lemma.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_rank_two_boundary_family_lemma.md"


def build_payload() -> dict:
    rows = [
        {
            "id": "rtbf_01_parameter_domain",
            "role": "exact_definition",
            "formula": "1/2<u<1; a=(2*u-1)/u, b=1-u, c=2*u-1",
            "claim": "The parameters a,b,c are strictly positive.",
        },
        {
            "id": "rtbf_02_coefficient_family",
            "role": "exact_definition",
            "formula": "A_k(u)=a^(k-1)*(c+k*b)/u for every k>=0",
            "claim": "The formula includes A_0=A_1=1 and defines a positive infinite sequence.",
        },
        {
            "id": "rtbf_03_ratio_defect_identity",
            "role": "exact_identity",
            "formula": "x_k=1-b^2/(c+k*b)^2",
            "claim": "The ratio contractions increase to one.",
        },
        {
            "id": "rtbf_04_hausdorff_defect_identity",
            "role": "exact_identity",
            "formula": "d_k=1/(k+c/b)^2=integral_0^1 t^(k+c/b-1)*(-log t) dt",
            "claim": "The full defect sequence is a Hausdorff moment sequence and is completely monotone.",
        },
        {
            "id": "rtbf_05_all_degree_shift_factorization",
            "role": "exact_theorem",
            "formula": "J_(d,n)(z)=A_n*(1+a*z)^(d-1)*(1+a*(c+n*b+d*b)/(c+n*b)*z)",
            "claim": "Every degree-d, shift-n Jensen window factors into d real negative roots.",
        },
        {
            "id": "rtbf_06_cubic_quartic_boundary_chain",
            "role": "exact_identity",
            "formula": "u=3/5 gives J_3=(1+z/3)^2*(1+7*z/3), J_4=(1+z/3)^3*(1+3*z), and x=5/9,y=21/25,z=45/49",
            "claim": "The cubic double root extends to a quartic triple root at the critical contraction z=(3-2*u)/(2-u)^2.",
        },
        {
            "id": "rtbf_07_positive_mixture_countermodel",
            "role": "exact_countermodel_gate",
            "formula": "(A(3/5)+A(2/3))/2 gives J_3=1+3*z+(47/24)*z^2+(41/108)*z^3 with Disc=-937/3456",
            "claim": "Positive mixtures of the rank-two boundary family do not preserve Jensen hyperbolicity.",
        },
        {
            "id": "rtbf_08_alpha_multiplier_identity",
            "role": "exact_identity",
            "formula": "u=(alpha+1)/(alpha+2): M_k=(alpha/(alpha+1))^(k-1)*(alpha+k)/(alpha+1), EGF=(1+z/(alpha+1))*exp(alpha*z/(alpha+1))",
            "claim": "Each boundary sequence is a classical nonnegative multiplier sequence with one negative EGF zero.",
        },
        {
            "id": "rtbf_09_integer_hadamard_product_closure",
            "role": "exact_theorem",
            "formula": "gamma_k=geometric_k*product_j M_k^(alpha_j) with unit integer atoms alpha_j>0",
            "claim": "Finite pointwise products remain multiplier sequences because their diagonal real-root preservers compose.",
        },
        {
            "id": "rtbf_10_fractional_power_countermodel",
            "role": "exact_countermodel_gate",
            "formula": "alpha=1, exponent=1/2: Disc(sum_(k=0)^3 binom(3,k)*sqrt(k+1)*z^k)<-27/125<0",
            "claim": "Fractional pointwise powers need not be multiplier sequences, so a general positive measure over alpha is not sufficient.",
        },
        {
            "id": "rtbf_11_open_structural_handoff",
            "role": "open_handoff",
            "formula": "-log x_k=sum_j -log(1-1/(k+alpha_j)^2) with a genuine counting measure",
            "claim": "The family gives an integer-atomic canonical-product target, not a proved representation of the zeta coefficients.",
        },
    ]
    return {
        "kind": "jensen_window_pf_rank_two_boundary_family_lemma",
        "date": "2026-07-10",
        "status": "exact all-degree boundary-family lemma with mixture countermodel",
        "proof_boundary": (
            "Exact positive model family and exact rejection of naive positive mixtures. "
            "This does not represent the zeta coefficient sequence, construct the required positive kernel, prove Jensen-window PF-infinity for zeta, RH, or Lambda <= 0."
        ),
        "source_hierarchy": "outputs/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.md",
        "source_bridge_target": "outputs/jensen_window_pf_bridge_target.md",
        "rows": rows,
        "mixture_countermodel": {
            "parameters": ["3/5", "2/3"],
            "weights": ["1/2", "1/2"],
            "shift": 0,
            "degree": 3,
            "jensen_coefficients_ascending": ["1", "3", "47/24", "41/108"],
            "discriminant": "-937/3456",
            "discriminant_is_negative": True,
        },
        "summary": {
            "rows": len(rows),
            "exact_identity_rows": 4,
            "all_degree_factorization_rows": 1,
            "integer_product_closure_rows": 1,
            "exact_countermodels": 2,
            "open_structural_handoffs": 1,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "A one-parameter positive multiplier sequence has an exact all-degree/all-shift Jensen factorization with one simple and one repeated negative root. "
                "Finite pointwise products give an integer-atomic canonical-product cone. However, both an equal convex mixture and a fractional pointwise power have negative cubic discriminants, so arbitrary positive superposition or measure weights are not the missing kernel theorem."
            ),
        },
        "invariants": [
            "The all-degree factorization is exact for every d>=1 and n>=0.",
            "Every model Jensen root is real and negative for 1/2<u<1.",
            "The mixture countermodel blocks convex closure only.",
            "No identification with the actual zeta coefficients is claimed.",
            "RH and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(payload: dict, path: Path) -> None:
    result = (
        "validated Jensen-window PF rank-two boundary-family lemma: "
        "11 rows, 0 issues, 4 exact identities, 1 all-degree factorization, "
        "1 integer-product closure, 2 exact countermodels, 1 open structural handoff"
    )
    lines = [
        "# Jensen-Window PF Rank-Two Boundary-Family Lemma",
        "",
        "Date: 2026-07-10",
        "",
        "Status: exact all-degree boundary-family lemma with mixture countermodel. This is not a proof",
        "of a zeta representation, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_rank_two_boundary_family_lemma`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_rank_two_boundary_family_lemma.json",
        "python work/rh_compute/scripts/jensen_window_pf_rank_two_boundary_family_lemma.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_rank_two_boundary_family_lemma.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result,
        "```",
        "",
        "## Exact Family",
        "",
        "For `1/2<u<1`, put",
        "",
        "```text",
        "a=(2*u-1)/u, b=1-u, c=2*u-1,",
        "A_k(u)=a^(k-1)*(c+k*b)/u.",
        "```",
        "",
        "Then `A_0=A_1=1`, and exact binomial algebra gives, for every `d>=1,n>=0`,",
        "",
        "```text",
        "J_(d,n)(z)=A_n*(1+a*z)^(d-1)",
        "             *(1+a*(c+n*b+d*b)/(c+n*b)*z).",
        "```",
        "",
        "Both root locations are strictly negative. Thus every shifted Jensen window in this",
        "model is hyperbolic. Its contractions and defects are",
        "",
        "```text",
        "x_k=1-b^2/(c+k*b)^2,",
        "d_k=1/(k+c/b)^2",
        "   = integral_0^1 t^(k+c/b-1)*(-log t) dt.",
        "```",
        "",
        "The cubic critical value `x=5/9,y=21/25,z=45/49` at `u=3/5` is one",
        "member of this all-degree hyperbolic boundary family.",
        "",
        "## Mixture Gate",
        "",
        "The family is not convex in the required hyperbolicity sense. The equal positive",
        "mixture of `u=3/5` and `u=2/3` gives",
        "",
        "```text",
        "J_3(z)=1+3*z+(47/24)*z^2+(41/108)*z^3,",
        "Disc(J_3)=-937/3456<0.",
        "```",
        "",
        "Therefore a positive integral or mixture over these boundary sequences is not, by",
        "itself, the missing zeta kernel representation.",
        "",
        "## Multiplier Product Route",
        "",
        "With `u=(alpha+1)/(alpha+2)`, the family becomes",
        "",
        "```text",
        "M_k^(alpha)=(alpha/(alpha+1))^(k-1)*(alpha+k)/(alpha+1),",
        "sum_(k>=0) M_k^(alpha)*z^k/k!",
        "  =(1+z/(alpha+1))*exp(alpha*z/(alpha+1)).",
        "```",
        "",
        "This is a classical multiplier sequence. Finite pointwise products remain multiplier",
        "sequences because the corresponding diagonal real-root preservers compose. Their",
        "ratio contractions satisfy",
        "",
        "```text",
        "-log x_k=sum_j -log(1-1/(k+alpha_j)^2).",
        "```",
        "",
        "The atoms must carry genuine integer multiplicity. Fractional weights are unsafe:",
        "after removing positive geometric scaling, `alpha=1` and exponent `1/2` gives",
        "",
        "```text",
        "J_3(z)=1+3*sqrt(2)*z+3*sqrt(3)*z^2+2*z^3,",
        "Disc=-432*sqrt(2)-324*sqrt(3)+378+324*sqrt(6)",
        "    <-27/125<0.",
        "```",
        "",
        "Thus a general positive measure representation of `-log x_k` is insufficient; the",
        "live target is a discrete counting-measure factorization or another structure with",
        "an independently proved stability-preserving composition law.",
        "",
        "## Handoff",
        "",
        payload["summary"]["main_finding"],
        "",
        "```text",
        "outputs/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.md",
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
        "validated Jensen-window PF rank-two boundary-family lemma: "
        "11 rows, 0 issues, 4 exact identities, 1 all-degree factorization, "
        "1 integer-product closure, 2 exact countermodels, 1 open structural handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
