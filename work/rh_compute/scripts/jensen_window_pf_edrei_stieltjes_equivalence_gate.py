#!/usr/bin/env python3
"""Record the exact Edrei-log/Stieltjes endpoint equivalence and its guards."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from math import factorial
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_edrei_stieltjes_equivalence_gate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_edrei_stieltjes_equivalence_gate.md"


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def multiply(left: list[Fraction], right: list[Fraction], degree: int) -> list[Fraction]:
    out = [Fraction(0) for _ in range(degree + 1)]
    for i, a_i in enumerate(left):
        for j, b_j in enumerate(right):
            if i + j <= degree:
                out[i + j] += a_i * b_j
    return out


def reciprocal(series: list[Fraction], degree: int) -> list[Fraction]:
    if series[0] != 1:
        raise ValueError("reciprocal audit expects unit constant term")
    out = [Fraction(0) for _ in range(degree + 1)]
    out[0] = Fraction(1)
    for n in range(1, degree + 1):
        out[n] = -sum(series[k] * out[n - k] for k in range(1, n + 1))
    return out


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    work = [row[:] for row in matrix]
    sign = 1
    det = Fraction(1)
    for col in range(len(work)):
        pivot = next((row for row in range(col, len(work)) if work[row][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            sign *= -1
        pivot_value = work[col][col]
        det *= pivot_value
        for row in range(col + 1, len(work)):
            factor = work[row][col] / pivot_value
            for j in range(col + 1, len(work)):
                work[row][j] -= factor * work[col][j]
    return det * sign


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def sample_audit() -> dict:
    degree = 12
    gamma = Fraction(2, 5)
    atoms = [(Fraction(1, 7), 2), (Fraction(2, 9), 1)]

    h = [gamma**n / factorial(n) for n in range(degree + 1)]
    for beta, multiplicity in atoms:
        for _ in range(multiplicity):
            h = multiply(h, [Fraction(1), beta], degree)

    derivative = [Fraction(n + 1) * h[n + 1] for n in range(degree)]
    inverse = reciprocal(h, degree - 1)
    log_derivative = multiply(derivative, inverse, degree - 1)

    moments = []
    for n in range(degree):
        value = (gamma if n == 0 else Fraction(0)) + sum(
            multiplicity * beta ** (n + 1) for beta, multiplicity in atoms
        )
        moments.append(value)
        if log_derivative[n] != (-1) ** n * value:
            raise RuntimeError(f"log-derivative/moment identity failed at n={n}")

    h0 = []
    for size in range(1, 5):
        matrix = [[moments[i + j] for j in range(size)] for i in range(size)]
        h0.append(determinant(matrix))
    h1 = []
    for size in range(1, 4):
        matrix = [[moments[i + j + 1] for j in range(size)] for i in range(size)]
        h1.append(determinant(matrix))

    if not all(value > 0 for value in h0[:3]) or h0[3] != 0:
        raise RuntimeError("unexpected unshifted finite-support Hankel pattern")
    if not all(value > 0 for value in h1[:2]) or h1[2] != 0:
        raise RuntimeError("unexpected shifted finite-support Hankel pattern")

    return {
        "degree": degree,
        "gamma": fraction_text(gamma),
        "atoms": [
            {"beta": fraction_text(beta), "multiplicity": multiplicity}
            for beta, multiplicity in atoms
        ],
        "h_coefficients": [fraction_text(value) for value in h],
        "log_derivative_coefficients": [
            fraction_text(value) for value in log_derivative
        ],
        "shifted_log_moments": [fraction_text(value) for value in moments],
        "hankel_s1_determinants_sizes_1_to_4": [
            fraction_text(value) for value in h0
        ],
        "hankel_s2_determinants_sizes_1_to_3": [
            fraction_text(value) for value in h1
        ],
        "indexing_checks": degree,
        "hankel_checks": len(h0) + len(h1),
    }


def build_exact() -> dict:
    return {
        "corpus_normalization": "H_lambda(z)=F_lambda(z)/F_lambda(0)=sum_(k>=0)d_k(lambda)z^k, H_lambda(0)=1, and H_lambda is entire",
        "log_coordinates": "log H(z)=sum_(n>=1)ell_n*z^n, q_n=n*ell_n, p_n=(-1)^(n-1)*q_n",
        "indexing_identity": "a_r:=(-1)^r*[z^r](H'(z)/H(z))=p_(r+1)",
        "sokal_criterion": "For entire H with H(0)!=0: H in LP+ <=> (a_r)_(r>=0) is a Stieltjes moment sequence",
        "unique_measure": "a_r=integral_[0,infinity)x^r*dnu(x), nu=gamma*delta_0+sum_j m_j*beta_j*delta_(beta_j), gamma>=0, beta_j>0, m_j in Z_(>=1), sum_j m_j*beta_j<infinity",
        "stieltjes_transform": "H'(z)/H(z)=integral_[0,infinity)dnu(x)/(1+x*z)=gamma+sum_j m_j*beta_j/(1+beta_j*z)",
        "hankel_psd": "a is Stieltjes <=> (a_(i+j))_(i,j>=0) and (a_(i+j+1))_(i,j>=0) are positive semidefinite on every finite block",
        "strict_two_column_target": "D_(m,1)=det(p_(i+j+1))_(i,j=0)^m>0 and D_(m,2)=det(p_(i+j+2))_(i,j=0)^m>0 for every m>=0 => a is Stieltjes",
        "s_fraction": "a is Stieltjes <=> sum_(r>=0)a_r*t^r has a formal Stieltjes continued fraction with nonnegative coefficients, including terminating cases",
        "endpoint_equivalence": "c is PF-infinity <=> H is LP+ <=> a_r=p_(r+1) is Stieltjes <=> every shifted Jensen window is finite PF-infinity",
        "phi_integral": "F_lambda(z)=integral_R exp(lambda*u^2)*Phi(u)*cosh(u*sqrt(z))*du",
        "phi_log_derivative": "H_lambda'(z)/H_lambda(z)=[integral_R exp(lambda*u^2)*Phi(u)*u*sinh(u*sqrt(z))/(2*sqrt(z))*du]/[integral_R exp(lambda*u^2)*Phi(u)*cosh(u*sqrt(z))*du]",
        "finite_guard": "finitely many positive D_(m,s), log signs, or recurrence coefficients do not imply the all-order Stieltjes moment property",
        "object_separation": "Hankel matrices of a_r=p_(r+1) are nonlinear Edrei-log objects and are not the original signed-Hankel matrices of A_k",
        "mixture_guard": "positivity of the Phi integral or fixed-scale LP+ kernels does not imply that their positive mixture has Stieltjes logarithmic derivative",
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            "esee_01_corpus_normalization",
            "exact_identity",
            "available_exact",
            "The normalized coefficient generating function is entire and has unit constant term.",
            exact["corpus_normalization"],
            "Known analytic normalization only; no zero-location conclusion.",
        ),
        GateRow(
            "esee_02_log_derivative_indexing",
            "exact_identity",
            "available_exact",
            "The shifted signed Edrei-log powers are exactly the Taylor coefficients of the signed logarithmic derivative.",
            exact["indexing_identity"],
            "Formal identity near zero; audited independently at exact rational coefficients.",
        ),
        GateRow(
            "esee_03_sokal_log_derivative_criterion",
            "classical_theorem",
            "ready_to_apply",
            "Sokal's logarithmic-derivative criterion makes the all-order shifted Stieltjes property equivalent to LP+ for an entire function.",
            exact["sokal_criterion"],
            "The criterion identifies an endpoint-equivalent target; it does not prove the Stieltjes antecedent for zeta.",
        ),
        GateRow(
            "esee_04_unique_edrei_measure",
            "classical_theorem",
            "ready_to_apply",
            "Entireness forces any Stieltjes representing measure for the logarithmic derivative into the discrete Edrei zero-measure form with integer residues.",
            exact["unique_measure"],
            "Conditional on the all-order Stieltjes property.",
        ),
        GateRow(
            "esee_05_stieltjes_transform_identity",
            "exact_equivalence",
            "ready_to_apply",
            "The representing measure is equivalently a Stieltjes-transform representation of the logarithmic derivative.",
            exact["stieltjes_transform"],
            "An exact target representation, not an established representation for H_0.",
        ),
        GateRow(
            "esee_06_hankel_psd_criterion",
            "classical_theorem",
            "ready_to_apply",
            "The Stieltjes moment property is equivalent to positivity of both unshifted and once-shifted Hankel quadratic forms at every finite size.",
            exact["hankel_psd"],
            "All sizes are required; finite blocks remain diagnostics.",
        ),
        GateRow(
            "esee_07_strict_two_column_target",
            "sufficient_theorem_target",
            "not_ready_to_apply",
            "Strict positivity of the two leading power-Hankel columns s=1 and s=2 at every order is a clean sufficient endpoint theorem.",
            exact["strict_two_column_target"],
            "Open for all m; currently checked only on a finite staircase.",
        ),
        GateRow(
            "esee_08_stieltjes_s_fraction",
            "classical_theorem",
            "ready_to_apply",
            "A nonnegative Stieltjes continued fraction is an exact alternative all-order certificate for the same shifted moment sequence.",
            exact["s_fraction"],
            "Requires every continued-fraction coefficient or a valid terminating certificate.",
        ),
        GateRow(
            "esee_09_unified_endpoint",
            "exact_equivalence",
            "ready_to_apply",
            "The Edrei-log Stieltjes target joins the already unified coefficient-PF and Jensen-window endpoint.",
            exact["endpoint_equivalence"],
            "Equivalence of targets only; no endpoint statement is proved for zeta.",
        ),
        GateRow(
            "esee_10_finite_evidence_alignment",
            "finite_evidence_map",
            "diagnostic_validated",
            "The existing log-sign, power-Hankel, and recurrence artifacts test finite pieces of this exact Stieltjes target.",
            "320 log-sign rows; 4205 power-Hankel rows; 55 positive recurrence rows through order 12",
            "Finite evidence only; higher recurrence orders are inconclusive and all-order positivity is open.",
        ),
        GateRow(
            "esee_11_original_hankel_separation",
            "object_separation_guard",
            "guard_validated",
            "The viable Edrei-log Hankel target is not invalidated by the order-ten failure of the original A_k signed-Hankel antecedent.",
            exact["object_separation"],
            "Prevents conflation of two different Hankel constructions.",
        ),
        GateRow(
            "esee_12_finite_prefix_guard",
            "forbidden_promotion",
            "guard_validated",
            "No finite determinant or recurrence prefix enters the all-order equivalence.",
            exact["finite_guard"],
            "Blocks finite-to-infinite extrapolation.",
        ),
        GateRow(
            "esee_13_phi_ratio_target",
            "xi_phi_exact_handoff",
            "not_ready_to_apply",
            "The Phi moment integral gives an exact ratio whose proof as a Stieltjes function would establish the endpoint noncircularly.",
            exact["phi_log_derivative"],
            "The ratio identity is exact; its Stieltjes property is open and does not follow from positivity of the integrands.",
        ),
        GateRow(
            "esee_14_positive_mixture_guard",
            "positive_mixture_guard",
            "guard_validated",
            "Positive mixing of fixed-scale LP+ kernels is not a preservation theorem for the logarithmic-derivative Stieltjes class.",
            exact["mixture_guard"],
            "The existing exact scale-mixture countermodel remains active.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_edrei_stieltjes_equivalence_gate",
        "date": "2026-07-23",
        "status": "exact Edrei-log/Stieltjes/LP+/coefficient-PF/Jensen endpoint equivalence gate",
        "proof_boundary": (
            "This artifact proves an exact equivalence of all-order target formulations, identifies a strict two-column Hankel sufficient target, and records an exact Xi/Phi logarithmic-derivative handoff. It does not prove the all-order Stieltjes property, coefficient PF-infinity, Jensen hyperbolicity for zeta, RH, or Lambda <= 0."
        ),
        "primary_sources": [
            "https://doi.org/10.1016/j.jmaa.2022.126432",
            "https://discovery.ucl.ac.uk/10150527/1/Sokal_1-s2.0-S0022247X22004462-main.pdf",
            "https://doi.org/10.1007/BF02786970",
            "https://doi.org/10.1007/BF02786971",
        ],
        "corpus_sources": [
            "outputs/jensen_window_pf_coefficient_pf_equivalence_gate.md",
            "outputs/edrei_log_sign_diagnostic.md",
            "outputs/edrei_moment_quadrature_scout.md",
            "outputs/jensen_window_pf_laguerre_scale_mixture_gate.md",
            "outputs/missing_coefficient_pf_theorem.md",
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
            "# Jensen-Window PF Edrei-Stieltjes Equivalence Gate",
            "",
            "Date: 2026-07-23",
            "",
            "Status: exact endpoint-equivalence and theorem-route correction. This is not",
            "a proof of the all-order Stieltjes property, coefficient PF-infinity, Jensen",
            "hyperbolicity for zeta, RH, or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_edrei_stieltjes_equivalence_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_edrei_stieltjes_equivalence_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_edrei_stieltjes_equivalence_gate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF Edrei-Stieltjes equivalence gate: 14 rows, 0 issues, 12 exact indexing checks, 7 exact Hankel checks, 1 unified endpoint, 3 finite/nonpromotion guards, 1 open Xi/Phi handoff",
            "```",
            "",
            "## Exact Indexing",
            "",
            "```text",
            exact["corpus_normalization"],
            exact["log_coordinates"],
            exact["indexing_identity"],
            "```",
            "",
            "The last line is the important off-by-one check: the sequence tested by the",
            "existing moment-recurrence scout, `a_r=p_(r+1)`, is exactly Sokal's signed",
            "Taylor-coefficient sequence for `H'/H`.",
            "",
            "## Classical Equivalence",
            "",
            "Sokal's Proposition 6 applies because `H` is entire and `H(0)=1`:",
            "",
            "```text",
            exact["sokal_criterion"],
            exact["stieltjes_transform"],
            "```",
            "",
            "Entireness is decisive. It forces the representing measure to be the unique",
            "discrete zero measure",
            "",
            "```text",
            exact["unique_measure"],
            "```",
            "",
            "so no extra integrality or product-reconstruction conjecture is needed after",
            "the all-order Stieltjes property has actually been proved.",
            "",
            "Combining this with the coefficient-PF/Jensen gate gives",
            "",
            "```text",
            exact["endpoint_equivalence"],
            "```",
            "",
            "This is an equivalence of endpoint targets, not evidence that the endpoint",
            "holds.",
            "",
            "## Exact Hankel Target",
            "",
            "The Stieltjes moment criterion is",
            "",
            "```text",
            exact["hankel_psd"],
            "```",
            "",
            "A clean strict sufficient target in the corpus notation is therefore",
            "",
            "```text",
            exact["strict_two_column_target"],
            "```",
            "",
            "Only the two columns `s=1,2` are structurally required. The wider finite",
            "staircase remains useful stress evidence, but it is not an all-order proof.",
            "The exact continued-fraction alternative is",
            "",
            "```text",
            exact["s_fraction"],
            "```",
            "",
            "## Existing Finite Evidence",
            "",
            "```text",
            "320 finite Edrei-log sign rows pass",
            "4205 finite power-Hankel rows pass",
            "55 finite recurrence rows through order 12 are positive",
            "orders 13..20 are interval-inconclusive, not negative",
            "```",
            "",
            "These are finite pieces of an endpoint-equivalent Stieltjes target. They remain",
            "diagnostics and cannot be promoted by extrapolation.",
            "",
            "## Xi/Phi Handoff",
            "",
            "The moment integral gives",
            "",
            "```text",
            exact["phi_integral"],
            exact["phi_log_derivative"],
            "```",
            "",
            "A noncircular proof that this ratio is a Stieltjes function at `lambda=0`",
            "would prove the common endpoint. Positivity of `Phi`, or the fact that every",
            "fixed-scale `cosh(u*sqrt(z))` lies in `LP+`, is insufficient: positive",
            "mixtures need not preserve the required zero or logarithmic-derivative class.",
            "",
            "## Object-Separation Guard",
            "",
            "```text",
            exact["object_separation"],
            "```",
            "",
            "Thus the known order-ten failure in the original signed-Hankel hierarchy does",
            "not refute this nonlinear Edrei-log Hankel target. It only forbids conflating",
            "the two objects.",
            "",
            "## Machine Audit",
            "",
            "The generator and independent checker use a rational finite Edrei product to",
            "verify the indexing and both Hankel columns:",
            "",
            "```text",
            f"indexing checks: {audit['indexing_checks']}",
            f"Hankel determinant checks: {audit['hankel_checks']}",
            f"s=1 signs by sizes 1..4: {audit['hankel_s1_determinants_sizes_1_to_4']}",
            f"s=2 signs by sizes 1..3: {audit['hankel_s2_determinants_sizes_1_to_3']}",
            "```",
            "",
            "The terminal zeros are the expected finite-support rank termination, not",
            "failures of the Stieltjes property.",
            "",
            "## Primary Source",
            "",
            "Sokal, Proposition 6 and Proposition 7:",
            "",
            "```text",
            "https://doi.org/10.1016/j.jmaa.2022.126432",
            "https://discovery.ucl.ac.uk/10150527/1/Sokal_1-s2.0-S0022247X22004462-main.pdf",
            "```",
            "",
            "## Boundary",
            "",
            payload["proof_boundary"],
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
        "wrote Jensen-window PF Edrei-Stieltjes equivalence gate: "
        f"{args.out.relative_to(REPO_ROOT)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
