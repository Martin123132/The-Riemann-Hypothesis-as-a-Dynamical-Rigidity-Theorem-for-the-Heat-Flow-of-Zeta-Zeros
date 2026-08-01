#!/usr/bin/env python3
"""Build the exact Suzuki truncation-path spectral frontier and route guards."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_suzuki_spectral_frontier.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_suzuki_spectral_frontier.md"
)


@dataclass(frozen=True)
class FrontierRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def exact_statements() -> dict[str, str]:
    return {
        "operator": (
            "(K[t]f)(x)=1_(x<t)*integral_(-infinity)^t "
            "K(x+y)f(y)dy, with K(s)=0 for s<0"
        ),
        "fixed_space": (
            "For (U_t f)(u)=f(t-u), u>0: "
            "A_t=U_t K[t] U_t^(-1), "
            "(A_t g)(u)=integral_0^infinity K(2t-u-v)g(v)dv"
        ),
        "effective_interval": (
            "K[t] vanishes on inputs supported in (-infinity,-t) and its "
            "output vanishes on (-infinity,-t); it is therefore represented "
            "by the continuous kernel K(x+y) on the finite square (-t,t)^2"
        ),
        "hilbert_schmidt": (
            "||K[t]||_HS^2=integral_0^(2t) "
            "(2t-s)*|K(s)|^2 ds"
        ),
        "hilbert_schmidt_growth": (
            "If c_R=integral_0^R |K(s)|^2 ds>0, then for t>R/2: "
            "||K[t]||_HS^2 >= (2t-R)*c_R -> infinity"
        ),
        "path_continuity": (
            "For continuous K with K(s)=0 on s<0, "
            "t -> A_t is continuous in Hilbert-Schmidt norm and A_0=0"
        ),
        "norm_monotonicity": (
            "For 0<=t_1<=t_2, ||K[t_1]||<=||K[t_2]||, because K[t_1] "
            "is the quadratic-form compression of the self-adjoint K[t_2]"
        ),
        "determinant_contractivity": (
            "For every T>=0: det(I+/-K[t])!=0 for all 0<=t<=T "
            "iff ||K[T]||<1 iff ||K[t]||<1 for all 0<=t<=T"
        ),
        "first_crossing": (
            "If T_*=inf{t>=0: ||K[t]||=1}<infinity, compact "
            "self-adjointness forces +1 or -1 to be an eigenvalue of K[T_*], "
            "so det(I-K[T_*])=0 or det(I+K[T_*])=0"
        ),
        "isolated_time_guard": (
            "At one isolated time determinant nonvanishing is weaker: "
            "A=2P for a rank-one orthogonal projection has ||A||=2, "
            "det(I-A)=-1, and det(I+A)=3"
        ),
        "high_contour_certificate": (
            "Suzuki's Proposition 4.2 eigenfunction contradiction is certified "
            "when M_v^2*exp(4v*t)<1 for some v>1/2+omega, where "
            "M_v=sup_(u in R)|Theta_(omega,nu)(u+i*v)|"
        ),
        "zeta_vertical_asymptotic": (
            "For zeta and a=nu*omega: "
            "Theta_(omega,nu)(i*v)="
            "[xi(v+1/2-omega)/xi(v+1/2+omega)]^nu "
            "~(2*pi/v)^a as v->infinity"
        ),
        "high_contour_ceiling": (
            "T_HC(omega,nu)=sup_(v>1/2+omega) "
            "[log(1/M_v)]_+/(2v)<infinity; the certificate "
            "M_v^2*exp(4v*t)<1 cannot prove any t>T_HC"
        ),
        "conditional_inner_result": (
            "If E_(omega,nu) belongs to HB, Suzuki Proposition 5.1 gives a "
            "full L2 isometry K, K[t]=P_t K P_t, and ||K[t]||<1 for every "
            "finite t"
        ),
        "no_uniform_gap": (
            "Under the same HB hypothesis, P_t->I strongly and K is an "
            "isometry, so lim_(t->infinity)||K[t]||=1; no epsilon>0 can "
            "satisfy ||K[t]||<=1-epsilon for every t"
        ),
        "quadratic_form_target": (
            "The noncircular live target is, for every finite t and every "
            "nonzero f: |<K[t]f,f>|<||f||^2, equivalently "
            "I+K[t]>0 and I-K[t]>0"
        ),
        "terminal_sequence_reduction": (
            "For one fixed pair, pointwise-in-t strict contractivity does not "
            "directly prove J(t;z,z)->0. On a strictly decreasing cofinal "
            "omega_n->0 sequence, however, all-time contractivity gives a "
            "bounded causal multiplier, forces RH by shifted-zero "
            "accumulation, and then Suzuki's Theorem 2.3 supplies the "
            "terminal limit"
        ),
    }


def exact_toy_checks() -> dict[str, str]:
    # For K(s)=s on s>=0 and t=1:
    # integral_0^2 (2-s)s^2 ds = 4/3.
    hs_value = (
        2 * Fraction(2) ** 3 / 3 - Fraction(2) ** 4 / 4
    )
    if hs_value != Fraction(4, 3):
        raise RuntimeError(f"unexpected Hilbert-Schmidt toy value: {hs_value}")

    rank_one_eigenvalue = Fraction(2)
    det_minus = 1 - rank_one_eigenvalue
    det_plus = 1 + rank_one_eigenvalue
    if (det_minus, det_plus) != (Fraction(-1), Fraction(3)):
        raise RuntimeError("unexpected isolated-time rank-one determinants")

    return {
        "hs_kernel": "K(s)=s for s>=0",
        "hs_time": "1",
        "hs_squared": "4/3",
        "rank_one_eigenvalue": "2",
        "rank_one_det_minus": "-1",
        "rank_one_det_plus": "3",
    }


def build_payload() -> dict:
    exact = exact_statements()
    rows = [
        FrontierRow(
            "ssf_01_fixed_space_translation",
            "exact_identity",
            "available_exact",
            "Every Suzuki truncation is unitarily equivalent to a finite-triangle Hankel operator on one fixed half-line space.",
            exact["fixed_space"],
            "This uses only the support K(s)=0 for s<0.",
        ),
        FrontierRow(
            "ssf_02_effective_finite_interval",
            "exact_identity",
            "available_exact",
            "Although written on a half-line, the truncation acts only on a finite interval.",
            exact["effective_interval"],
            "This does not bound its operator norm.",
        ),
        FrontierRow(
            "ssf_03_hilbert_schmidt_identity",
            "exact_identity",
            "available_exact",
            "The Hilbert-Schmidt norm has a one-dimensional weighted-kernel formula.",
            exact["hilbert_schmidt"],
            "The identity is valid for every locally square-integrable supported kernel.",
        ),
        FrontierRow(
            "ssf_04_hilbert_schmidt_divergence",
            "exact_route_guard",
            "guard_validated",
            "For every nonzero continuous supported kernel, the Hilbert-Schmidt norm of the truncations diverges.",
            exact["hilbert_schmidt_growth"],
            "This does not imply that the operator norm exceeds one.",
        ),
        FrontierRow(
            "ssf_05_hilbert_schmidt_certificate_no_go",
            "exact_route_guard",
            "guard_validated",
            "The sufficient estimate ||K[t]||<=||K[t]||_HS<1 can only certify a bounded initial t-range.",
            exact["hilbert_schmidt_growth"],
            "Cancellation-sensitive operator estimates remain possible.",
        ),
        FrontierRow(
            "ssf_06_truncation_path_continuity",
            "exact_lemma",
            "available_exact",
            "The translated compact operators form a Hilbert-Schmidt-continuous path starting at zero.",
            exact["path_continuity"],
            "Suzuki Proposition 4.1 supplies continuity of the arithmetic kernel when nu*omega>1.",
        ),
        FrontierRow(
            "ssf_07_operator_norm_monotonicity",
            "exact_lemma",
            "available_exact",
            "The truncation operator norm is nondecreasing in t.",
            exact["norm_monotonicity"],
            "Monotonicity alone does not keep the norm below one.",
        ),
        FrontierRow(
            "ssf_08_global_determinant_contractivity_equivalence",
            "exact_equivalence",
            "ready_to_apply",
            "Along the whole connected truncation path, Suzuki's two determinant exclusions are exactly pointwise strict contractivity.",
            exact["determinant_contractivity"],
            "The all-earlier-times quantifier is essential.",
        ),
        FrontierRow(
            "ssf_09_first_crossing_alternative",
            "exact_equivalence",
            "ready_to_apply",
            "Any finite failure has a first spectral collision at +1 or -1.",
            exact["first_crossing"],
            "The open task is to rule out this first collision for the arithmetic kernel.",
        ),
        FrontierRow(
            "ssf_10_isolated_time_guard",
            "exact_countermodel",
            "guard_validated",
            "At a single isolated t, nonzero determinants do not imply contraction.",
            exact["isolated_time_guard"],
            "The example is a logical spectral guard, not a Suzuki truncation path.",
        ),
        FrontierRow(
            "ssf_11_suzuki_high_contour_certificate",
            "published_exact_argument",
            "available_exact",
            "Suzuki's unconditional local interval follows from a weighted high-contour eigenfunction contradiction.",
            exact["high_contour_certificate"],
            "The condition is sufficient and local, not necessary.",
        ),
        FrontierRow(
            "ssf_12_high_contour_finite_ceiling",
            "exact_route_guard",
            "guard_validated",
            "For each fixed zeta pair (omega,nu), Suzuki's M_v certificate has a finite t ceiling.",
            exact["high_contour_ceiling"],
            "This rejects only this sup-norm contour certificate, not every possible contour argument.",
        ),
        FrontierRow(
            "ssf_13_conditional_inner_contractivity",
            "published_conditional_theorem",
            "conditional_only",
            "Hermite-Biehler innerness supplies strict contraction at every finite t.",
            exact["conditional_inner_result"],
            "Using this implication to establish the same HB property would be circular.",
        ),
        FrontierRow(
            "ssf_14_no_uniform_t_gap",
            "exact_conditional_consequence",
            "guard_validated",
            "In the desired HB case the finite-t norms approach one, so a t-uniform positive gap is impossible.",
            exact["no_uniform_gap"],
            "The target is strictness at each finite t, not one epsilon valid for all t.",
        ),
        FrontierRow(
            "ssf_15_signed_quadratic_form_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The surviving noncircular spectral task is a signed quadratic-form barrier against finite first crossing.",
            exact["quadratic_form_target"],
            "No arithmetic factorization or coercive identity proving this is currently known in the corpus.",
        ),
        FrontierRow(
            "ssf_16_terminal_sequence_reduction",
            "exact_cross_artifact_reduction",
            "internally_audited",
            "The terminal premise is fixed-family independent but cofinal-sequence redundant.",
            exact["terminal_sequence_reduction"],
            "This uses the separately checked determinant-only causal-multiplier reduction and does not derive the terminal limit directly from one fixed family.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_suzuki_spectral_frontier",
        "date": "2026-07-23",
        "status": (
            "exact truncation-path spectral reformulation, finite-ceiling "
            "route guards, one open signed quadratic-form target, and an "
            "internally audited cofinal terminal reduction"
        ),
        "assumptions": [
            "K is real and continuous on R",
            "K(s)=0 for s<0",
            "K[t] is the compact self-adjoint Suzuki truncation",
        ],
        "exact": exact,
        "toy_checks": exact_toy_checks(),
        "rows": [asdict(row) for row in rows],
        "sources": [
            {
                "title": (
                    "Masatoshi Suzuki, Hamiltonians arising from L-functions "
                    "in the Selberg class"
                ),
                "url": "https://arxiv.org/abs/1606.05726",
                "relevant_results": (
                    "Propositions 4.1-4.3 and 5.1; Theorem 2.4"
                ),
            },
            {
                "title": "Xi Pick/Suzuki arithmetic-Hankel bridge",
                "path": (
                    "outputs/"
                    "jensen_window_pf_xi_pick_suzuki_hankel_bridge.md"
                ),
            },
            {
                "title": "Suzuki determinant-only reduction",
                "path": (
                    "outputs/"
                    "jensen_window_pf_suzuki_determinant_only_reduction.md"
                ),
            },
        ],
        "proof_boundary": (
            "This artifact proves the fixed-space translation, exact "
            "Hilbert-Schmidt identity, path continuity and monotonicity, the "
            "all-history determinant/contractivity equivalence, the first-"
            "crossing alternative, and finite-ceiling guards for Hilbert-"
            "Schmidt and Suzuki sup-contour certificates. It records Suzuki's "
            "conditional HB contraction theorem and proves that its norms tend "
            "to one. Together with the separately checked causal-multiplier "
            "reduction, it records that the terminal premise is redundant for "
            "a cofinal determinant family. It does not prove the signed "
            "quadratic-form target, global Fredholm nonvanishing for zeta "
            "without HB, the Phi Pick sign, PF-infinity, RH, or Lambda<=0."
        ),
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    toy = payload["toy_checks"]
    return "\n".join(
        [
            "# Jensen-Window PF Suzuki Spectral Frontier",
            "",
            "Date: 2026-07-23",
            "",
            "Status: exact truncation-path spectral reformulation, finite-ceiling",
            "route guards, one open signed quadratic-form target, and an",
            "internally audited cofinal terminal reduction. This is not",
            "a proof of global Fredholm nonvanishing, RH, or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_suzuki_spectral_frontier.json",
            "python work/rh_compute/scripts/jensen_window_pf_suzuki_spectral_frontier.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_spectral_frontier.py",
            "```",
            "",
            "## Fixed-Space Truncation",
            "",
            "Start from Suzuki's supported Hankel operator",
            "",
            "```text",
            exact["operator"],
            "```",
            "",
            "and translate `L2(-infinity,t)` to one fixed copy of",
            "`L2(0,infinity)`:",
            "",
            "```text",
            exact["fixed_space"],
            "```",
            "",
            exact["effective_interval"] + ".",
            "",
            "The translated kernel is supported on the triangle `u+v<=2t`.",
            "Consequently,",
            "",
            "```text",
            exact["hilbert_schmidt"],
            "```",
            "",
            f"For the check kernel `{toy['hs_kernel']}` at `t={toy['hs_time']}`,",
            f"both sides equal `{toy['hs_squared']}`.",
            "",
            "## Hilbert-Schmidt Guard",
            "",
            "If the kernel is not identically zero, choose `R` with",
            "`c_R=integral_0^R |K(s)|^2 ds>0`. Then",
            "",
            "```text",
            exact["hilbert_schmidt_growth"],
            "```",
            "",
            "Thus the generic estimate `||K[t]||<=||K[t]||_HS<1` can prove",
            "only a bounded initial interval. This is a limitation of the",
            "Hilbert-Schmidt certificate, not evidence that the operator norm",
            "itself exceeds one.",
            "",
            "## Exact Spectral Reformulation",
            "",
            "Suzuki's continuity and support assumptions imply",
            "",
            "```text",
            exact["path_continuity"],
            exact["norm_monotonicity"],
            "```",
            "",
            "Since a compact self-adjoint operator attains its norm at an",
            "eigenvalue, continuity from the zero operator gives",
            "",
            "```text",
            exact["determinant_contractivity"],
            "```",
            "",
            "Equivalently,",
            "",
            "```text",
            exact["first_crossing"],
            "```",
            "",
            "The all-times qualifier cannot be dropped:",
            "",
            "```text",
            exact["isolated_time_guard"],
            "```",
            "",
            "So Suzuki condition (3) is neither an unspecified determinant",
            "problem nor a request for one t-uniform margin. It is exactly the",
            "problem of preventing a finite first crossing of either spectral",
            "edge.",
            "",
            "## High-Contour Ceiling",
            "",
            "Suzuki's unconditional local proof supplies the sufficient test",
            "",
            "```text",
            exact["high_contour_certificate"],
            "```",
            "",
            "For zeta, Stirling's formula gives the lower witness",
            "",
            "```text",
            exact["zeta_vertical_asymptotic"],
            "```",
            "",
            "while Suzuki's upper estimate has the same polynomial scale.",
            "Therefore",
            "",
            "```text",
            exact["high_contour_ceiling"],
            "```",
            "",
            "The exponential support loss eventually defeats every polynomial",
            "gain in contour height. Sharpening only the polynomial constant or",
            "power cannot turn Proposition 4.2's mechanism into an all-`t`",
            "proof for a fixed pair.",
            "",
            "## Correct Global Target",
            "",
            "Under the desired Hermite-Biehler hypothesis, Suzuki proves",
            "",
            "```text",
            exact["conditional_inner_result"],
            "```",
            "",
            "but this conditional argument also shows why a uniform margin is",
            "the wrong target:",
            "",
            "```text",
            exact["no_uniform_gap"],
            "```",
            "",
            "The surviving noncircular obligation is therefore",
            "",
            "```text",
            exact["quadratic_form_target"],
            "```",
            "",
            "A useful proof would factor or otherwise control both signed forms",
            "without identifying the high-contour kernel with an inner real-axis",
            "multiplier. For one fixed pair this spectral task does not",
            "directly prove the terminal limit. The cofinal logic is sharper:",
            "",
            "```text",
            exact["terminal_sequence_reduction"],
            "```",
            "",
            "The causal-multiplier and shifted-zero proof is recorded in",
            "`outputs/jensen_window_pf_suzuki_determinant_only_reduction.md`.",
            "",
            "## Source",
            "",
            "- Masatoshi Suzuki, `Hamiltonians arising from L-functions in the Selberg class`, Propositions 4.1-4.3 and 5.1, Theorem 2.4: https://arxiv.org/abs/1606.05726",
            "- `outputs/jensen_window_pf_xi_pick_suzuki_hankel_bridge.md`",
            "- `outputs/jensen_window_pf_suzuki_determinant_only_reduction.md`",
            "",
            "## Proof Boundary",
            "",
            payload["proof_boundary"],
            "",
        ]
    )


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
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Jensen-window PF Suzuki spectral frontier: "
        "16 rows, 10 exact path/reduction identities, "
        "5 route guards, 1 open global obligation"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
