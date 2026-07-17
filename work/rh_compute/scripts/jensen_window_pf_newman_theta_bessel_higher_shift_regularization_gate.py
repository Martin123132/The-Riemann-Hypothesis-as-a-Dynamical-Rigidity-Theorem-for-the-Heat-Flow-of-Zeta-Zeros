#!/usr/bin/env python3
"""Expose the exact theta/Bessel higher-shift expansion and its summation gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def coefficient_sign_rows() -> list[dict]:
    rows: list[dict] = []
    for n_value in (1, 2, 3):
        threshold = 2 * sp.pi**2 * n_value**4 / 3
        last_positive = int(sp.floor(sp.N(threshold, 80)))
        rows.append(
            {
                "n": n_value,
                "threshold": str(threshold),
                "last_positive_m": last_positive,
                "first_negative_m": last_positive + 1,
                "sample_signs": [
                    {
                        "m": m_value,
                        "sign": (
                            "positive"
                            if sp.N(2 * sp.pi**2 * n_value**4 - 3 * m_value, 80) > 0
                            else "negative"
                        ),
                    }
                    for m_value in sorted({0, last_positive, last_positive + 1})
                ],
            }
        )
    return rows


def build_exact() -> dict:
    n, m = sp.symbols("n m", integer=True, positive=True)
    p = sp.pi * n**2
    numerator = 2 * sp.pi**2 * n**4 * p**m - 3 * sp.pi * n**2 * m * p ** (m - 1)
    expected = p**m * (2 * sp.pi**2 * n**4 - 3 * m)
    if sp.simplify(numerator - expected) != 0:
        raise RuntimeError("higher-shift coefficient identity failed")

    f1, f2, f1p, f2p, f1pp, f2pp, c1, c2 = sp.symbols(
        "f1 f2 f1p f2p f1pp f2pp c1 c2", real=True
    )
    direct = (c1 * f1p + c2 * f2p) ** 2 - (c1 * f1 + c2 * f2) * (
        c1 * f1pp + c2 * f2pp
    )
    matrix = (
        c1**2 * (f1p**2 - f1 * f1pp)
        + 2
        * c1
        * c2
        * (f1p * f2p - (f1 * f2pp + f2 * f1pp) / 2)
        + c2**2 * (f2p**2 - f2 * f2pp)
    )
    if sp.expand(direct - matrix) != 0:
        raise RuntimeError("mixed Laguerre matrix identity failed")

    return {
        "theta_summand": {
            "definition": (
                "phi_n(u)=(2*pi^2*n^4*exp(9u)-3*pi*n^2*exp(5u))*"
                "exp(-pi*n^2*exp(4u))"
            ),
            "shift": "phi_n(u)=n^(-1/2)*phi_1(u+(log n)/2)",
            "full_evenness": "Phi(u)=sum_(n>=1)phi_n(u)=Phi(-u)",
            "symmetrized_block": "g_n(u)=(phi_n(u)+phi_n(-u))/2",
            "pointwise_sum": "Phi(u)=sum_(n>=1)g_n(u)",
        },
        "higher_shift_expansion": {
            "coefficient": (
                "c_(n,m)=pi^m*n^(2m)*(2*pi^2*n^4-3*m)/m!"
            ),
            "identity": (
                "g_n(u)=exp(-2*pi*n^2*cosh(4u))*sum_(m>=0) "
                "c_(n,m)*cosh((9-4*m)*u)"
            ),
            "derivation": (
                "Symmetrize phi_n first, factor exp(-pi*n^2*(exp(4u)+exp(-4u))), "
                "expand the remaining exponentials, and pair the exp(9u) term at "
                "index m with the exp(5u) term at index m-1."
            ),
            "convergence": (
                "For each fixed n the m-series and all u-derivatives converge "
                "locally uniformly; the absolute m-majorant is integrable on [0,infinity)."
            ),
            "absolute_majorant": (
                "(1/2)*(2*pi^2*n^4)*(exp(9u-pi*n^2*exp(4u))+"
                "exp(-9u-pi*n^2*exp(-4u))) + "
                "(1/2)*(3*pi*n^2)*(exp(5u-pi*n^2*exp(4u))+"
                "exp(-5u-pi*n^2*exp(-4u)))"
            ),
        },
        "coefficient_sign": {
            "law": "sign(c_(n,m))=sign(2*pi^2*n^4-3*m)",
            "threshold": "c_(n,m)>0 iff m<2*pi^2*n^4/3",
            "no_zero": "c_(n,m) is never zero because pi^2 is irrational",
            "first_block": (
                "For n=1, c_(1,m)>0 for 0<=m<=6 and c_(1,m)<0 for m>=7; "
                "the bounds 3<pi<22/7 place 2*pi^2/3 strictly between 6 and 7."
            ),
            "rows": coefficient_sign_rows(),
        },
        "fixed_block_transform": {
            "definition": (
                "P_(n,m)(x)=1/8*(K_(9/4-m+i*x/4)(2*pi*n^2)+"
                "K_(9/4-m-i*x/4)(2*pi*n^2))"
            ),
            "bessel_integral": (
                "integral_0^infinity exp(-2*pi*n^2*cosh(4u))*"
                "cosh((9-4m)u)*cos(xu)du=P_(n,m)(x)"
            ),
            "termwise_transform": (
                "I_n(x):=integral_0^infinity g_n(u)cos(xu)du="
                "sum_(m>=0)c_(n,m)*P_(n,m)(x)"
            ),
            "scope": (
                "The identity is an ordinary absolutely justified transform for each fixed n."
            ),
        },
        "spectral_summation_obstruction": {
            "profile_transform": (
                "Qhat(x)=-(1+x^2)/32*pi^(-(1+i*x)/4)*Gamma((1+i*x)/4)"
            ),
            "shift_transform": (
                "Fourier[phi_n](x)=n^(-1/2-i*x/2)*Qhat(x)"
            ),
            "symmetrized_cosine_transform": (
                "I_n(x)=1/4*(n^(-1/2-i*x/2)*Qhat(x)+"
                "n^(-1/2+i*x/2)*Qhat(-x))"
            ),
            "zero_frequency": (
                "I_n(0)=Qhat(0)/(2*sqrt(n)), "
                "Qhat(0)=-Gamma(1/4)/(32*pi^(1/4))<0"
            ),
            "divergence": "sum_(n>=1)I_n(0)=-infinity",
            "finite_target": (
                "H_0(0)=integral_0^infinity Phi(u)du is finite, so "
                "H_0 is not the ordinary termwise sum of the Bessel-block transforms."
            ),
            "fubini_boundary": (
                "The n-sum is pointwise and locally uniform in the kernel variable, "
                "but it is not interchangeable with the half-line Fourier integral."
            ),
        },
        "analytic_continuation_guard": {
            "formal_assembly": (
                "Qhat(x)*sum_(n>=1)n^(-1/2-i*x/2)"
            ),
            "continued_identity": "Qhat(x)*zeta((1+i*x)/2)=xi((1+i*x)/2)/4",
            "decision": (
                "The Dirichlet series is outside its convergence half-plane on the "
                "real spectral axis. Analytic continuation reconstructs xi itself; "
                "it is not a positive or zero-free block summation theorem."
            ),
        },
        "coupled_matrix_handoff": {
            "matrix_entry": (
                "M_(i,j)=F_i'*F_j'-(F_i*F_j''+F_j*F_i'')/2"
            ),
            "quadratic_identity": (
                "L[sum_j a_j*F_j]=sum_(i,j)a_i*a_j*M_(i,j)"
            ),
            "required_upgrade": (
                "Group the theta summands through the modular identity before the "
                "spectral transform, then prove positivity of the resulting coupled "
                "mixed-term quadratic form. Alternatively supply a noncircular "
                "renormalized summation theorem that preserves the needed sign."
            ),
            "subsequent_endpoint_result": (
                "The theta cell-renormalization gate supplies an ordinary normally "
                "convergent coupled matrix at t=0, but proves that every individual "
                "cell block retains an exp(-5u) tail and cannot be deformed to t>0."
            ),
        },
    }


def build_artifact() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="ntbhsr_01_theta_symmetrization",
            role="exact_identity",
            readiness="available_exact",
            claim="The even Xi kernel is the pointwise sum of symmetrized arithmetic theta blocks.",
            formula=exact["theta_summand"]["pointwise_sum"],
            proof_boundary="Exact kernel identity; no spectral interchange is asserted.",
            diagnostics=exact["theta_summand"],
        ),
        GateRow(
            id="ntbhsr_02_higher_shift_expansion",
            role="exact_identity",
            readiness="available_exact",
            claim="Each symmetrized theta block has an exact infinite expansion in higher Bessel-shift kernels.",
            formula=exact["higher_shift_expansion"]["identity"],
            proof_boundary="Exact for each fixed n with locally uniform kernel convergence.",
            diagnostics=exact["higher_shift_expansion"],
        ),
        GateRow(
            id="ntbhsr_03_coefficient_formula",
            role="exact_identity",
            readiness="available_exact",
            claim="The two neighboring exponential series combine into one closed coefficient formula.",
            formula=exact["higher_shift_expansion"]["coefficient"],
            proof_boundary="Exact coefficient algebra only.",
        ),
        GateRow(
            id="ntbhsr_04_coefficient_sign",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Every arithmetic block has a sharp positive-to-negative coefficient threshold.",
            formula=exact["coefficient_sign"]["threshold"],
            proof_boundary="Signs the expansion coefficients; it does not sign mixed Laguerre terms.",
            diagnostics=exact["coefficient_sign"],
        ),
        GateRow(
            id="ntbhsr_05_fixed_block_bessel_transform",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="For each fixed theta index, the higher-shift series may be Fourier transformed term by term.",
            formula=exact["fixed_block_transform"]["termwise_transform"],
            proof_boundary="Fixed n only; summation over n is deliberately excluded.",
            diagnostics=exact["fixed_block_transform"],
        ),
        GateRow(
            id="ntbhsr_06_zero_frequency_divergence",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The ordinary sum of the fixed-block transforms diverges already at zero frequency.",
            formula=exact["spectral_summation_obstruction"]["divergence"],
            proof_boundary="Rejects Fubini promotion of this block expansion, not the exact kernel identity.",
            diagnostics=exact["spectral_summation_obstruction"],
        ),
        GateRow(
            id="ntbhsr_07_fubini_boundary",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Pointwise theta convergence does not supply an ordinary Bessel-transform decomposition of Xi.",
            formula=exact["spectral_summation_obstruction"]["fubini_boundary"],
            proof_boundary="A modularly grouped or renormalized transform may still exist.",
        ),
        GateRow(
            id="ntbhsr_08_continuation_guard",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Analytic continuation of the divergent arithmetic transform reconstructs xi rather than a new positive factor.",
            formula=exact["analytic_continuation_guard"]["continued_identity"],
            proof_boundary="No zero information or strict Laguerre positivity is obtained.",
            diagnostics=exact["analytic_continuation_guard"],
        ),
        GateRow(
            id="ntbhsr_09_coupled_matrix_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The surviving route is a modularly grouped mixed-term matrix identity or a sign-preserving renormalized summation theorem.",
            formula=exact["coupled_matrix_handoff"]["quadratic_identity"],
            proof_boundary="Open; this is not Wiener density, RH, or Lambda<=0.",
            diagnostics=exact["coupled_matrix_handoff"],
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate",
        "date": "2026-07-11",
        "status": "exact higher-shift expansion with spectral summation obstruction",
        "proof_boundary": (
            "This artifact proves an exact higher-shift expansion for every "
            "symmetrized theta summand, its coefficient sign law, and a fixed-index "
            "Bessel transform theorem. It also proves that summing those transforms "
            "over the arithmetic index is invalid in the ordinary sense, with an "
            "explicit divergence at zero frequency. It does not reject modularly "
            "grouped or renormalized coupled identities and does not prove strict "
            "Laguerre positivity, Wiener density, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_theta_summand_spectral_square_gate.md",
            "outputs/jensen_window_pf_newman_signed_universal_factor_residual_gate.md",
            "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md",
            "outputs/formal_core.md",
            "https://dlmf.nist.gov/10.32",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Theta/Bessel Higher-Shift Regularization Gate",
            "",
            "Date: 2026-07-11",
            "",
            "Status: exact higher-shift expansion with spectral summation obstruction.",
            "This is not a proof or disproof of RH or `Lambda <= 0`.",
            "",
            "Artifact kind: `jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_bessel_higher_shift_regularization_gate.py",
            "```",
            "",
            "Builds on `outputs/jensen_window_pf_newman_theta_summand_spectral_square_gate.md`.",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF Newman theta/Bessel higher-shift regularization gate: 9 rows, 0 issues, 3 exact expansion identities, 1 coefficient sign theorem, 1 fixed-block Bessel theorem, 3 spectral non-promotion gates, 1 coupled modular handoff",
            "```",
            "",
            "## Exact Higher Shifts",
            "",
            "Set",
            "",
            "```text",
            exact["theta_summand"]["definition"],
            exact["theta_summand"]["symmetrized_block"],
            exact["theta_summand"]["pointwise_sum"],
            "```",
            "",
            "Symmetrizing before expanding gives, for every fixed `n>=1`,",
            "",
            "```text",
            exact["higher_shift_expansion"]["coefficient"],
            exact["higher_shift_expansion"]["identity"],
            "```",
            "",
            "The series and all kernel derivatives converge locally uniformly. Its",
            "absolute `m`-majorant is integrable on the positive half-line for each",
            "fixed `n`, so the individual block transform is ordinary, not formal.",
            "",
            "## Sign Law",
            "",
            "```text",
            exact["coefficient_sign"]["law"],
            exact["coefficient_sign"]["threshold"],
            exact["coefficient_sign"]["first_block"],
            "```",
            "",
            "Thus higher shifts really do enter, but with infinitely many signed",
            "coefficients. The sign change is intrinsic and cannot be removed by",
            "treating every shifted block as an independent positive residual.",
            "",
            "## Fixed-Block Transform",
            "",
            "The standard Bessel integral gives",
            "",
            "```text",
            exact["fixed_block_transform"]["definition"],
            exact["fixed_block_transform"]["bessel_integral"],
            exact["fixed_block_transform"]["termwise_transform"],
            "```",
            "",
            "This theorem is deliberately indexed by one fixed `n`. The arithmetic",
            "sum cannot then be moved through the Fourier integral.",
            "",
            "## Spectral Obstruction",
            "",
            "The translate identity gives exactly",
            "",
            "```text",
            exact["spectral_summation_obstruction"]["shift_transform"],
            exact["spectral_summation_obstruction"]["zero_frequency"],
            exact["spectral_summation_obstruction"]["divergence"],
            "```",
            "",
            "But `H_0(0)=integral_0^infinity Phi(u)du` is finite. Therefore the",
            "exact pointwise kernel expansion is not an ordinary termwise Bessel",
            "decomposition of the Xi transform. This is a concrete failure of",
            "Fubini/Tonelli, not a failure of the theta identity.",
            "",
            "Analytic continuation does not repair the proof route:",
            "",
            "```text",
            exact["analytic_continuation_guard"]["formal_assembly"],
            exact["analytic_continuation_guard"]["continued_identity"],
            "```",
            "",
            "It reconstructs Xi itself and supplies no positive block summation.",
            "",
            "## Live Handoff",
            "",
            "```text",
            exact["coupled_matrix_handoff"]["matrix_entry"],
            exact["coupled_matrix_handoff"]["quadratic_identity"],
            "```",
            "",
            exact["coupled_matrix_handoff"]["required_upgrade"],
            "",
            exact["coupled_matrix_handoff"]["subsequent_endpoint_result"],
            "See `outputs/jensen_window_pf_newman_theta_cell_renormalization_gate.md`.",
            "",
            "Naive ordinary higher-shift summation is now closed. A surviving proof",
            "must retain the modular cancellation inside the coupled mixed terms",
            "before taking the spectral transform.",
            "",
            "Reference: https://dlmf.nist.gov/10.32.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(f"wrote Newman theta/Bessel higher-shift regularization gate: {len(artifact['rows'])} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
