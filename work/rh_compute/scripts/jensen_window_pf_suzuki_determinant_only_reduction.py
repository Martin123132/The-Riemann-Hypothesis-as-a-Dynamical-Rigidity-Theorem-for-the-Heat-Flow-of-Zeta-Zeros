#!/usr/bin/env python3
"""Build the audited determinant-only Suzuki reduction for RH."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_determinant_only_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_suzuki_determinant_only_reduction.md"
)


@dataclass(frozen=True)
class ReductionRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def exact_statements() -> dict[str, str]:
    return {
        "determinant_to_contraction": (
            "det(I+/-K[t])!=0 for every t>=0 "
            "implies ||K[t]||<1 for every finite t"
        ),
        "dense_form": (
            "D=L2_c(R)={compactly supported L2 functions}; "
            "B(f,g)=integral_R integral_R "
            "K(x+y)f(y)conj(g(x)) dy dx"
        ),
        "uniform_form_bound": (
            "For f,g in D and t above both supports: "
            "B(f,g)=<K[t]f,g>, hence "
            "|B(f,g)|<=||f||_2||g||_2"
        ),
        "global_hankel_extension": (
            "B extends uniquely to a bounded Hermitian form on L2(R), "
            "represented by a self-adjoint H with ||H||<=1"
        ),
        "reflection_convolution": (
            "Jf(x)=f(-x), G=H*J; for f in D, "
            "(Gf)(x)=integral_R K(x-y)f(y)dy=(K*f)(x)"
        ),
        "translation_causality": (
            "G commutes with every translation and "
            "P_a G=P_a G P_a for P_a=1_(-infinity,a); "
            "thus G is bounded, translation-invariant, and causal"
        ),
        "causal_multiplier": (
            "There is M in H-infinity({Re(s)>0}) with "
            "||M||_infinity=||G||<=1 and "
            "L(Gf)(s)=M(s)L(f)(s) for causal f"
        ),
        "high_strip_identification": (
            "For Re(s)>c: M(s)=integral_0^infinity "
            "K(x)e^(-s*x)dx=Theta_(omega,nu)(i*s)"
        ),
        "probe_nonvanishing": (
            "For q=1_[0,1], L(q)(s)=(1-e^(-s))/s!=0 "
            "when Re(s)>0, so the high-strip multiplier identity "
            "can be divided by L(q)(s)"
        ),
        "theta_extension": (
            "Theta_tilde(z)=M(-i*z) belongs to H-infinity(C+) with "
            "||Theta_tilde||_infinity<=1 and agrees with "
            "[xi(1/2-omega-i*z)/xi(1/2+omega-i*z)]^nu "
            "on Im(z)>c"
        ),
        "pole_removal": (
            "The meromorphic identity theorem makes every apparent pole "
            "of the xi quotient in C+ removable"
        ),
        "shifted_zero": (
            "If xi(rho)=0 and Re(rho)>1/2+omega, pole removal at "
            "z=-Im(rho)+i*(Re(rho)-1/2-omega) forces "
            "xi(rho-2*omega)=0"
        ),
        "cofinal_accumulation": (
            "If omega_n decreases strictly to 0, an off-line zero rho "
            "would force distinct zeros rho-2*omega_n -> rho for all "
            "large n, contradicting isolation of zeros of nonzero entire xi"
        ),
        "determinant_only_equivalence": (
            "RH iff there exist strictly decreasing omega_n->0 and "
            "positive integers nu_n with nu_n*omega_n>1 such that "
            "det(I+/-K_(omega_n,nu_n)[t])!=0 for every n and every t>=0"
        ),
        "terminal_redundancy": (
            "Within Suzuki's cofinal criterion, the terminal condition "
            "J_(omega_n,nu_n)(t;z,z)->0 follows after the determinant "
            "conditions imply RH; it is not an additional logical premise"
        ),
        "single_pair_guard": (
            "For one fixed omega, holomorphic continuation permits a "
            "denominator zero to be canceled by a numerator zero at "
            "rho-2*omega; cofinal distinct shifts are essential"
        ),
        "open_arithmetic_gate": (
            "The reduction does not establish the antecedent: one still "
            "must prove ||K_(omega_n,nu_n)[t]||<1 for every finite t "
            "on a cofinal omega_n sequence without assuming RH or HB"
        ),
    }


def exact_coordinate_check() -> dict[str, str]:
    beta = Fraction(3, 4)
    gamma = Fraction(7)
    omega = Fraction(1, 8)
    imaginary_z = beta - Fraction(1, 2) - omega
    denominator_real = Fraction(1, 2) + omega + imaginary_z
    numerator_real = Fraction(1, 2) - omega + imaginary_z
    if denominator_real != beta:
        raise RuntimeError("denominator zero coordinate did not recover beta")
    if numerator_real != beta - 2 * omega:
        raise RuntimeError("numerator zero coordinate did not shift by 2 omega")
    if imaginary_z <= 0:
        raise RuntimeError("test zero did not map into the upper half-plane")
    return {
        "rho": "3/4+7i",
        "omega": "1/8",
        "z": "-7+i/8",
        "denominator_argument": "3/4+7i",
        "numerator_argument": "1/2+7i",
        "shift": "1/4",
    }


def build_payload() -> dict:
    exact = exact_statements()
    rows = [
        ReductionRow(
            "sdr_01_all_history_contraction",
            "exact_input_reformulation",
            "available_exact",
            "Suzuki's all-time determinant condition supplies a uniform norm-one bound on every compactly supported test pair.",
            exact["determinant_to_contraction"],
            "This uses the continuous compact self-adjoint truncation path from K[0]=0.",
        ),
        ReductionRow(
            "sdr_02_dense_hankel_form",
            "exact_lemma",
            "available_exact",
            "The local truncations define one compatible Hankel form on compactly supported functions.",
            exact["dense_form"],
            "For each test pair, any t above both supports gives the same value.",
        ),
        ReductionRow(
            "sdr_03_uniform_form_bound",
            "exact_lemma",
            "available_exact",
            "All-time strict contractivity bounds the dense Hankel form by one independently of support size.",
            exact["uniform_form_bound"],
            "Strict constants may approach one; only the uniform non-strict bound is used.",
        ),
        ReductionRow(
            "sdr_04_global_hankel_extension",
            "exact_lemma",
            "available_exact",
            "The dense form has a unique bounded global Hankel-operator extension.",
            exact["global_hankel_extension"],
            "This is a Riesz-representation argument and does not assume HB or RH.",
        ),
        ReductionRow(
            "sdr_05_reflection_to_convolution",
            "exact_identity",
            "available_exact",
            "Reflection turns the global Hankel extension into convolution by Suzuki's supported kernel.",
            exact["reflection_convolution"],
            "Equality is first distributional on compact tests and then holds in L2.",
        ),
        ReductionRow(
            "sdr_06_causal_translation_invariance",
            "exact_lemma",
            "available_exact",
            "The reflected operator is a bounded causal translation-invariant L2 operator.",
            exact["translation_causality"],
            "Causality uses only K(x)=0 for x<0.",
        ),
        ReductionRow(
            "sdr_07_hardy_multiplier",
            "published_functional_analysis",
            "ready_to_apply",
            "The causal L2 multiplier theorem gives a bounded holomorphic transfer function on the right half-plane.",
            exact["causal_multiplier"],
            "This is the standard Paley-Wiener multiplier theorem, not an arithmetic assumption.",
        ),
        ReductionRow(
            "sdr_08_high_strip_identification",
            "exact_identity",
            "available_exact",
            "The transfer function is Suzuki's absolutely convergent high-strip Laplace transform.",
            exact["high_strip_identification"],
            exact["probe_nonvanishing"],
        ),
        ReductionRow(
            "sdr_09_theta_hardy_extension",
            "exact_consequence",
            "available_exact",
            "All-time determinant nonvanishing forces the meromorphic xi quotient to have a contractive H-infinity extension to C+.",
            exact["theta_extension"],
            "The equality starts on Im(z)>c and extends meromorphically by identity.",
        ),
        ReductionRow(
            "sdr_10_pole_removal",
            "exact_consequence",
            "available_exact",
            "Every upper-half-plane denominator zero is canceled by a numerator zero of at least the same multiplicity.",
            exact["pole_removal"],
            "Cancellation need not make the underlying E function Hermite-Biehler for one fixed omega.",
        ),
        ReductionRow(
            "sdr_11_shifted_zero_rule",
            "exact_arithmetic_consequence",
            "available_exact",
            "A zero to the right of the shifted critical line propagates left by exactly 2 omega.",
            exact["shifted_zero"],
            "Only existence, not multiplicity, is needed in the cofinal argument.",
        ),
        ReductionRow(
            "sdr_12_cofinal_accumulation",
            "exact_arithmetic_consequence",
            "available_exact",
            "A cofinal sequence of distinct shifts rules out every zero to the right of the critical line.",
            exact["cofinal_accumulation"],
            "The functional equation then rules out zeros to the left.",
        ),
        ReductionRow(
            "sdr_13_determinant_only_equivalence",
            "theorem_candidate",
            "internally_audited",
            "Suzuki's cofinal RH criterion can be sharpened by deleting its terminal-kernel premise.",
            exact["determinant_only_equivalence"],
            "The implication from RH is Suzuki's conditional all-time contraction theorem; the reverse implication is rows sdr_01 through sdr_12.",
        ),
        ReductionRow(
            "sdr_14_terminal_redundancy",
            "exact_corollary_candidate",
            "internally_audited",
            "For a cofinal family satisfying the determinant premise, Suzuki's terminal limit follows after RH and is logically redundant.",
            exact["terminal_redundancy"],
            "This does not derive the terminal limit directly from one fixed determinant family.",
        ),
        ReductionRow(
            "sdr_15_single_pair_cancellation_guard",
            "exact_route_guard",
            "guard_validated",
            "A single omega is insufficient because common-zero cancellation can survive.",
            exact["single_pair_guard"],
            "Do not promote one fixed contractive family to a zero-free half-plane without an extra no-common-zero theorem.",
        ),
        ReductionRow(
            "sdr_16_open_all_time_gate",
            "open_theorem_target",
            "not_ready_to_apply",
            "The determinant-only equivalence isolates one remaining Suzuki-side arithmetic gate.",
            exact["open_arithmetic_gate"],
            "No proof of this all-time signed contraction gate is presently in the corpus.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_suzuki_determinant_only_reduction",
        "date": "2026-07-23",
        "status": (
            "internally audited exact reduction and theorem candidate; "
            "external expert review still required before publication"
        ),
        "assumptions": [
            "L is the Riemann zeta function and xi is its nonzero entire completion",
            "omega>0, nu is a positive integer, and nu*omega>1",
            "Suzuki's K_(omega,nu) is real, continuous, supported in [0,infinity), and has an absolutely convergent high-strip transform",
            "both Fredholm determinants are nonzero for every finite truncation time",
        ],
        "exact": exact,
        "coordinate_check": exact_coordinate_check(),
        "rows": [asdict(row) for row in rows],
        "sources": [
            {
                "title": (
                    "Masatoshi Suzuki, Hamiltonians arising from L-functions "
                    "in the Selberg class"
                ),
                "url": "https://arxiv.org/abs/1606.05726",
                "relevant_results": (
                    "Theorems 2.3-2.4; Propositions 4.1-4.4 in arXiv v3"
                ),
            },
            {
                "title": (
                    "Guiver, Logemann, and Opmeer, Operator-valued "
                    "multiplier theorems for causal translation-invariant "
                    "operators"
                ),
                "url": (
                    "https://link.springer.com/article/"
                    "10.1007/s00498-024-00387-4"
                ),
                "relevant_results": (
                    "classical L2 causal multiplier equivalence, equations "
                    "(1.1) and (1.3)"
                ),
            },
            {
                "title": "Suzuki spectral frontier",
                "path": (
                    "outputs/"
                    "jensen_window_pf_suzuki_spectral_frontier.md"
                ),
            },
        ],
        "proof_boundary": (
            "This artifact proves, modulo the cited standard causal-multiplier "
            "theorem, that a cofinal sequence of Suzuki all-time determinant "
            "conditions implies RH and hence makes Suzuki's terminal J-limit "
            "premise redundant in that cofinal criterion. It does not prove "
            "that the determinant conditions hold for zeta, does not rule out "
            "common-zero cancellation for one fixed omega, and does not prove "
            "RH or Lambda<=0. The determinant-only equivalence is a new "
            "corpus theorem candidate and should receive independent expert "
            "review before being represented as a published result."
        ),
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    check = payload["coordinate_check"]
    return "\n".join(
        [
            "# Jensen-Window PF Suzuki Determinant-Only Reduction",
            "",
            "Date: 2026-07-23",
            "",
            "Status: internally audited exact reduction and theorem candidate.",
            "It is not a proof of RH or `Lambda <= 0`.",
            "It sharpens the logical form of Suzuki's criterion but does not establish the",
            "still-open all-time determinant premise. Independent expert review remains",
            "required before publication.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_suzuki_determinant_only_reduction.json",
            "python work/rh_compute/scripts/jensen_window_pf_suzuki_determinant_only_reduction.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_suzuki_determinant_only_reduction.py",
            "```",
            "",
            "## Statement",
            "",
            "For zeta, Suzuki's Theorem 2.4 remains equivalent to RH after",
            "condition (4), the terminal canonical-kernel limit, is deleted:",
            "",
            "```text",
            exact["determinant_only_equivalence"],
            "```",
            "",
            "The reverse implication is not a proof that these determinants",
            "are nonzero. It is a reduction showing that this one all-time",
            "operator condition would already be enough.",
            "",
            "## 1. From Every Truncation To One Global Form",
            "",
            "The spectral-frontier lemma gives",
            "",
            "```text",
            exact["determinant_to_contraction"],
            "```",
            "",
            "On the dense subspace of compactly supported L2 functions define",
            "",
            "```text",
            exact["dense_form"],
            "```",
            "",
            "For any fixed `f,g`, choose `t` above both supports. Then",
            "",
            "```text",
            exact["uniform_form_bound"],
            "```",
            "",
            "The bound is independent of the chosen support cutoff. Therefore",
            "",
            "```text",
            exact["global_hankel_extension"],
            "```",
            "",
            "No limit of the pointwise strict constants is claimed; their",
            "common non-strict upper bound one is exactly what the extension",
            "requires.",
            "",
            "## 2. Reflection Produces A Causal Multiplier",
            "",
            "Let `J` be reflection. A change of variables gives",
            "",
            "```text",
            exact["reflection_convolution"],
            "```",
            "",
            "The equality holds first against compact tests. Since both sides",
            "represent the same distribution and `G` is bounded on L2, it",
            "identifies Suzuki's locally defined kernel with the global",
            "convolution distribution. Its support now supplies",
            "",
            "```text",
            exact["translation_causality"],
            "```",
            "",
            "The standard L2 Paley-Wiener multiplier theorem therefore gives",
            "",
            "```text",
            exact["causal_multiplier"],
            "```",
            "",
            "This is the step that must not be replaced by an unjustified",
            "real-axis contour shift: boundedness is obtained from all finite",
            "forms before the multiplier theorem is invoked.",
            "",
            "## 3. Identify The Arithmetic Transfer Function",
            "",
            "Suzuki proves absolute convergence of the kernel transform in a",
            "high half-plane. There, ordinary Fubini calculation gives",
            "",
            "```text",
            exact["high_strip_identification"],
            "```",
            "",
            "A concrete division guard is",
            "",
            "```text",
            exact["probe_nonvanishing"],
            "```",
            "",
            "Rotating the transfer function to the upper half-plane yields",
            "",
            "```text",
            exact["theta_extension"],
            "```",
            "",
            "The two meromorphic functions agree on a nonempty open strip, so",
            "",
            "```text",
            exact["pole_removal"],
            "```",
            "",
            "This conclusion derives inner-type analyticity from the all-time",
            "finite-section bound; it does not assume HB, RH, or a contour",
            "shift across unknown poles.",
            "",
            "## 4. Pole Cancellation Becomes A Zero Shift",
            "",
            "Let `rho=beta+i*gamma` be a zero of xi with",
            "`beta>1/2+omega`. The corresponding denominator zero occurs at",
            "",
            "```text",
            "z=-gamma+i*(beta-1/2-omega) in C+.",
            "```",
            "",
            "Pole removal then gives",
            "",
            "```text",
            exact["shifted_zero"],
            "```",
            "",
            f"The exact coordinate check uses `rho={check['rho']}`,",
            f"`omega={check['omega']}`, and `z={check['z']}`. The denominator",
            f"argument is `{check['denominator_argument']}` and the numerator",
            f"argument is `{check['numerator_argument']}`, a real shift of",
            f"`{check['shift']}=2*omega`.",
            "",
            "For one omega this is only cancellation, not a contradiction.",
            "For a strictly decreasing cofinal sequence, however,",
            "",
            "```text",
            exact["cofinal_accumulation"],
            "```",
            "",
            "Thus xi has no zeros with real part greater than one half. Its",
            "functional equation excludes the reflected left half, proving",
            "the reverse implication in the displayed equivalence. Conversely,",
            "RH gives Suzuki's HB hypothesis and hence all-time strict",
            "contractivity and determinant nonvanishing.",
            "",
            "## Terminal Condition",
            "",
            "The resulting logical consequence is",
            "",
            "```text",
            exact["terminal_redundancy"],
            "```",
            "",
            "This is sequence-level redundancy. It is not the stronger and",
            "generally unjustified claim that one fixed all-time determinant",
            "family directly forces its terminal limit. The guard is",
            "",
            "```text",
            exact["single_pair_guard"],
            "```",
            "",
            "## Surviving Proof Obligation",
            "",
            "The sharpened Suzuki route now has one arithmetic gate:",
            "",
            "```text",
            exact["open_arithmetic_gate"],
            "```",
            "",
            "Finite grids, the unconditional local interval, Hilbert-Schmidt",
            "domination, and an RH/HB-derived isometry do not close this gate.",
            "",
            "## Sources",
            "",
            "- Masatoshi Suzuki, `Hamiltonians arising from L-functions in the Selberg class`, Theorems 2.3-2.4 and Propositions 4.1-4.4 in arXiv v3: https://arxiv.org/abs/1606.05726",
            "- Chris Guiver, Hartmut Logemann, and Mark R. Opmeer, `Operator-valued multiplier theorems for causal translation-invariant operators with applications to control theoretic input-output stability`, classical L2 equivalence in equations (1.1) and (1.3): https://link.springer.com/article/10.1007/s00498-024-00387-4",
            "- `outputs/jensen_window_pf_suzuki_spectral_frontier.md`",
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
        "wrote Suzuki determinant-only reduction: "
        "16 rows, 12 exact bridge steps, 2 theorem/corollary candidates, "
        "1 route guard, 1 open arithmetic gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
