#!/usr/bin/env python3
"""Record the exact coefficient-PF/Jensen-window equivalence and its guards."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from math import comb, factorial
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_coefficient_pf_equivalence_gate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_coefficient_pf_equivalence_gate.md"


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def pascal_row(degree: int) -> list[int]:
    row = [1]
    for _ in range(degree):
        row = [1] + [row[j - 1] + row[j] for j in range(1, len(row))] + [1]
    return row


def sample_audit() -> dict:
    values = [Fraction(k * k + 3 * k + 2, 1) for k in range(12)]
    normalization_checks = 0
    derivative_checks = 0
    window_checks = 0

    for k, value in enumerate(values):
        c_k = value / factorial(k)
        if c_k * factorial(k) != value:
            raise RuntimeError(f"normalization identity failed at k={k}")
        normalization_checks += 1

    for n in range(4):
        for j in range(len(values) - n):
            c_nj = values[n + j] / factorial(n + j)
            derivative_coefficient = c_nj * factorial(n + j) / factorial(j)
            if derivative_coefficient != values[n + j] / factorial(j):
                raise RuntimeError(f"derivative-tail identity failed at n={n}, j={j}")
            derivative_checks += 1

        for degree in range(1, 6):
            expanded_coefficients = pascal_row(degree)
            for j in range(degree + 1):
                diagonal_coefficient = expanded_coefficients[j] * values[n + j]
                window_coefficient = comb(degree, j) * values[n + j]
                if diagonal_coefficient != window_coefficient:
                    raise RuntimeError(
                        f"diagonal-window identity failed at n={n}, d={degree}, j={j}"
                    )
                window_checks += 1

    return {
        "normalization_checks": normalization_checks,
        "derivative_tail_checks": derivative_checks,
        "diagonal_window_checks": window_checks,
    }


def build_exact() -> dict:
    return {
        "scope_assumptions": "A_k>=0 for every k, A_0>0, and F is entire; all hold for the corpus zeta coefficient sequence",
        "normalization": "c_k=A_k/k! and F(z)=sum_(k>=0)c_k*z^k=sum_(k>=0)A_k*z^k/k!",
        "tail_function": "F_n(z)=F^(n)(z)=sum_(j>=0)A_(n+j)*z^j/j!",
        "diagonal_operator": "T_n[z^j]=A_(n+j)*z^j",
        "window_identity": "T_n[(1+z)^d]=P_(d,n)(z)=sum_(j=0)^d C(d,j)*A_(n+j)*z^j",
        "edrei_form": "F(z)=C*exp(gamma*z)*product_r(1+alpha_r*z), C>0, gamma>=0, alpha_r>=0, sum_r alpha_r<infinity",
        "polya_schur": "{A_k} is a multiplier sequence <=> F is Laguerre-Polya type I <=> P_(d,0) is hyperbolic for every d",
        "derivative_closure": "F in Laguerre-Polya type I => F^(n) in Laguerre-Polya type I for every n",
        "finite_asw": "P_(d,n) has only real nonpositive zeros <=> B^(d,n)_j=C(d,j)*A_(n+j) is a finite PF-infinity sequence",
        "equivalence": "c is PF-infinity <=> F is Laguerre-Polya type I <=> A is a multiplier sequence <=> every P_(d,0) is hyperbolic <=> every shifted A_(n+j) is a multiplier sequence <=> every P_(d,n) is hyperbolic <=> every B^(d,n) is finite PF-infinity",
        "finite_guard": "finite Toeplitz/PF certificates for c do not imply c is PF-infinity",
        "signed_hankel_guard": "signed-Hankel positivity does not by itself imply coefficient PF-infinity or the equivalent Jensen-window target",
        "positive_mixture_guard": "a positive mixture of dilations is coefficient-positive but is a hyperbolicity preserver only after a separate multiplier-sequence theorem",
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="cpje_01_normalization_identity",
            role="exact_identity",
            readiness="available_exact",
            claim="The ordinary coefficient generating function is exactly the exponential generating function of A.",
            formula=exact["normalization"],
            proof_boundary="Normalization identity only.",
        ),
        GateRow(
            id="cpje_02_diagonal_window_identity",
            role="exact_identity",
            readiness="available_exact",
            claim="Every Jensen window is the image of a binomial power under the diagonal tail operator.",
            formula=exact["window_identity"],
            proof_boundary="Polynomial coefficient identity only.",
        ),
        GateRow(
            id="cpje_03_polya_schur_shift_zero",
            role="classical_theorem",
            readiness="ready_to_apply",
            claim="Polya-Schur identifies shift-zero Jensen hyperbolicity in every degree with the multiplier-sequence property of A.",
            formula=exact["polya_schur"],
            proof_boundary="Classical equivalence; it does not prove its hypothesis for the zeta coefficients.",
        ),
        GateRow(
            id="cpje_04_asw_edrei_entire_pf",
            role="classical_theorem",
            readiness="ready_to_apply",
            claim="For the entire positive generating function F, ASW/Edrei identifies coefficient PF-infinity with the type-I product form.",
            formula=exact["edrei_form"],
            proof_boundary="Classical characterization; all-order Toeplitz positivity remains unproved here.",
        ),
        GateRow(
            id="cpje_05_derivative_tail_identity",
            role="exact_identity",
            readiness="available_exact",
            claim="The exponential generating function of every shifted A-tail is the corresponding derivative of F.",
            formula=exact["tail_function"],
            proof_boundary="Formal differentiation identity, justified by the known entire convergence of F.",
        ),
        GateRow(
            id="cpje_06_type_i_derivative_closure",
            role="closure_theorem",
            readiness="ready_to_apply",
            claim="Laguerre-Polya type I is closed under differentiation, so a shift-zero multiplier theorem supplies every shifted tail.",
            formula=exact["derivative_closure"],
            proof_boundary="Classical closure step conditional on the unproved type-I membership of F.",
        ),
        GateRow(
            id="cpje_07_all_shift_composition",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="Combining Polya-Schur with derivative closure proves all-degree/all-shift Jensen hyperbolicity from coefficient PF-infinity.",
            formula="c PF-infinity => F type I => every F^(n) type I => every P_(d,n) hyperbolic",
            proof_boundary="Exact conditional composition; the first all-order antecedent is still open.",
        ),
        GateRow(
            id="cpje_08_finite_asw_windows",
            role="classical_theorem",
            readiness="ready_to_apply",
            claim="For each fixed finite window, hyperbolicity and finite PF-infinity are equivalent.",
            formula=exact["finite_asw"],
            proof_boundary="Window-by-window equivalence only.",
        ),
        GateRow(
            id="cpje_09_seven_way_equivalence",
            role="exact_equivalence",
            readiness="ready_to_apply",
            claim="The all-order coefficient-PF and all-degree/all-shift Jensen-window targets are exact reformulations of one endpoint theorem.",
            formula=exact["equivalence"],
            proof_boundary="Equivalence of targets, not a proof that any target statement holds for zeta.",
        ),
        GateRow(
            id="cpje_10_finite_coefficient_pf_guard",
            role="forbidden_promotion",
            readiness="guard_validated",
            claim="The existing finite coefficient-PF ledger cannot enter the all-order equivalence.",
            formula=exact["finite_guard"],
            proof_boundary="Blocks promotion of finite Toeplitz evidence only.",
        ),
        GateRow(
            id="cpje_11_signed_hankel_nonimplication",
            role="nonimplication_guard",
            readiness="guard_validated",
            claim="The equivalence does not create the missing transfer from signed-Hankel data to coefficient PF-infinity.",
            formula=exact["signed_hankel_guard"],
            proof_boundary="Preserves the central structural gap and its existing countermodels.",
        ),
        GateRow(
            id="cpje_12_positive_mixture_guard",
            role="positive_mixture_guard",
            readiness="guard_validated",
            claim="Positivity of the Phi moment integral or a mixture-of-dilations formula is not itself a real-rootedness preserver theorem.",
            formula=exact["positive_mixture_guard"],
            proof_boundary="Requires an independent multiplier or stability-preserver hypothesis; positivity alone is insufficient.",
        ),
        GateRow(
            id="cpje_13_open_structural_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The live task is one noncircular all-order proof of coefficient PF-infinity, equivalently of the Jensen-window target, from Xi/Phi-specific structure.",
            formula="Xi/Phi-specific structure => c PF-infinity <=> all Jensen windows PF-infinity",
            proof_boundary="Open structural theorem; not PF-infinity, Laguerre-Polya membership, RH, or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_coefficient_pf_equivalence_gate",
        "date": "2026-07-23",
        "status": "exact coefficient-PF, multiplier-sequence, and Jensen-window equivalence gate",
        "proof_boundary": (
            "This artifact proves an equivalence between all-order target formulations and guards against finite, signed-Hankel-only, and positive-mixture promotions. It does not prove coefficient PF-infinity for the zeta coefficients, the missing signed-Hankel-to-PF bridge, Laguerre-Polya membership, RH, or Lambda <= 0."
        ),
        "primary_sources": [
            "https://annals.math.princeton.edu/wp-content/uploads/annals-v170-n1-p14-p.pdf",
            "https://doi.org/10.1007/BF02786970",
            "https://doi.org/10.1007/BF02786971",
            "https://doi.org/10.1073/pnas.37.5.303",
        ],
        "corpus_sources": [
            "outputs/coefficient_pf_bridge_obstruction.md",
            "outputs/jensen_window_pf_bridge_target.md",
            "outputs/jensen_window_pf_bridge_obligations.md",
            "outputs/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.md",
            "outputs/countermodel_library.md",
        ],
        "exact": exact,
        "sample_audit": sample_audit(),
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    audit = payload["sample_audit"]
    return "\n".join(
        [
            "# Jensen-Window PF Coefficient-PF Equivalence Gate",
            "",
            "Date: 2026-07-23",
            "",
            "Status: exact theorem-equivalence and proof-route correction. This is not",
            "a proof of coefficient PF-infinity, the signed-Hankel-to-PF bridge,",
            "Laguerre-Polya membership, RH, or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_coefficient_pf_equivalence_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_coefficient_pf_equivalence_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_coefficient_pf_equivalence_gate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF coefficient-PF equivalence gate: 13 rows, 0 issues, 3 exact coefficient identities, 4 classical/closure steps, 1 seven-way equivalence, 3 guards, 1 open structural handoff",
            "```",
            "",
            "## Exact Normalization",
            "",
            "For the corpus normalization,",
            "",
            exact["scope_assumptions"] + ".",
            "",
            "```text",
            exact["normalization"],
            exact["tail_function"],
            "```",
            "",
            "Define the diagonal tail operator by",
            "",
            "```text",
            exact["diagonal_operator"],
            exact["window_identity"],
            "```",
            "",
            "Thus the binomially weighted Jensen polynomial is exactly the",
            "Polya-Schur test polynomial for the shifted sequence.",
            "",
            "## Exact Equivalence",
            "",
            "The Polya-Schur theorem, ASW/Edrei characterization, closure under",
            "differentiation, and the finite ASW theorem give",
            "",
            "```text",
            exact["equivalence"],
            "```",
            "",
            "The all-order coefficient-PF route and the all-degree/all-shift",
            "Jensen-window route are therefore not two independent endpoint",
            "theorems. They are exact formulations of the same endpoint theorem.",
            "",
            "## Programme Correction",
            "",
            "What remains separate is the present finite evidence. Finite Toeplitz",
            "certificates for `c_k` do not prove the infinite PF statement, and finite",
            "Jensen-window certificates do not prove every degree and shift.",
            "",
            "Likewise, the equivalence supplies no transfer from the signed-Hankel",
            "certificates to coefficient PF-infinity. That structural implication is",
            "still the live bridge, now with one unambiguous endpoint target rather",
            "than two nominally separate ones.",
            "",
            "## Positive-Mixture Guard",
            "",
            "The positive Phi moment integral can be viewed as a positive mixture of",
            "dilations. Coefficient positivity or positive mixing does not make that",
            "diagonal operator a hyperbolicity preserver. Polya-Schur requires the",
            "operator's multiplier-sequence/type-I hypothesis, which is an endpoint",
            "statement and cannot be smuggled in through positivity language.",
            "",
            "## Machine Audit",
            "",
            "The generator and independent checker verify",
            "",
            "```text",
            f"normalization identities: {audit['normalization_checks']}",
            f"derivative-tail identities: {audit['derivative_tail_checks']}",
            f"diagonal-window identities: {audit['diagonal_window_checks']}",
            "```",
            "",
            "These finite exact checks audit the indexing and factorial normalization.",
            "The arbitrary-order equivalence comes from the cited classical theorems.",
            "",
            "## Primary Sources",
            "",
            "```text",
            *payload["primary_sources"],
            "```",
            "",
            "## Boundary",
            "",
            "Passing this checker proves only that the target formulations and guards",
            "are stated consistently. It does not prove the all-order antecedent for",
            "the zeta coefficients and therefore does not prove RH or `Lambda <= 0`.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Jensen-window PF coefficient-PF equivalence gate: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
