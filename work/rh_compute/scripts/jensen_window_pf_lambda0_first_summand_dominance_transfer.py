#!/usr/bin/env python3
"""Transfer first-summand dominance from lambda=-100 through lambda=0."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_lambda0_first_summand_dominance_transfer.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_lambda0_first_summand_dominance_transfer.md"
)
SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_dominance_certificate.json"
)

TAIL_START_K = 300
LAMBDA0_SPLICE_K = 504


@dataclass(frozen=True)
class TransferRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_source() -> dict:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    if source.get("kind") != (
        "jensen_window_pf_negative_lambda_first_summand_dominance_certificate"
    ):
        raise RuntimeError("first-summand dominance source has the wrong kind")
    if source.get("summary", {}).get("tail_start_k") != TAIL_START_K:
        raise RuntimeError("first-summand dominance source has the wrong tail start")
    matches = [
        row for row in source.get("rows", []) if row.get("id") == "fsdc_08_full_moment_dominance"
    ]
    if len(matches) != 1 or matches[0].get("readiness") != "interval_validated":
        raise RuntimeError("full-moment dominance source row is not validated")
    if matches[0].get("formula") != (
        "M_k=M_k^(1)*(1+delta_k), 0<=delta_k<=2/k^6"
    ):
        raise RuntimeError("full-moment dominance source formula changed")
    return source


def exact_arithmetic() -> dict:
    n = sp.symbols("n", nonnegative=True, integer=True)
    k = n + LAMBDA0_SPLICE_K
    polynomial = sp.Poly(
        sp.expand((k - 3) ** 6 - 25280 * k**2 * (k + 1) ** 2), n
    )
    coefficients = [int(value) for value in polynomial.all_coeffs()]
    if any(value <= 0 for value in coefficients):
        raise RuntimeError("lambda-zero perturbation polynomial is not coefficient-positive")

    first_index_error = 10112 * Fraction(
        (LAMBDA0_SPLICE_K + 1) ** 2,
        (LAMBDA0_SPLICE_K - 3) ** 6,
    )
    first_index_cap = Fraction(2, 5 * LAMBDA0_SPLICE_K**2)
    if not first_index_error < first_index_cap:
        raise RuntimeError("lambda-zero perturbation endpoint comparison failed")

    return {
        "heat_parameter": "T=-lambda, 0<=T<=100",
        "first_summand_measure": (
            "dmu_(k,T)(u)=Z_(k,T)^(-1)*u^(2k)*exp(-T*u^2)*Phi_1(u)du"
        ),
        "tail_ratio": "epsilon(u)=sum_(m>=2) Phi_m(u)/Phi_1(u)",
        "relative_tail": "delta_k(T)=E_(mu_(k,T))[epsilon(U)]",
        "covariance_derivative": "delta_k'(T)=-Cov_(mu_(k,T))(epsilon(U),U^2)",
        "two_copy_identity": (
            "2*Cov(f(U),g(U))=E[(f(U)-f(V))*(g(U)-g(V))]"
        ),
        "opposite_monotonicity": (
            "epsilon is decreasing and u^2 is increasing, hence Cov(epsilon(U),U^2)<=0"
        ),
        "heat_monotonicity": "delta_k'(T)>=0 on 0<=T<=100",
        "uniform_dominance": (
            "0<=delta_k(T)<=delta_k(100)<=2/k^6, k>=300, 0<=T<=100"
        ),
        "lambda_zero_dominance": "0<=delta_k(0)<=2/k^6, every integer k>=300",
        "moment_error": (
            "a_j=2*((j-1)^(-6)+2*j^(-6)+(j+1)^(-6)); "
            "|B_j-B_j^(1)|<=a_j"
        ),
        "penalty_error": (
            "|P_n-P_n^(1)|<=10112*(k+1)^2/(k-3)^6<=2/(5*k^2), "
            "k=n+3>=504"
        ),
        "splice": "lambda=0 prefix ends at n=500; tail starts at n=501, k=n+3=504",
        "perturbation_polynomial": str(polynomial.as_expr()),
        "perturbation_coefficients_descending": coefficients,
        "endpoint_error": str(first_index_error),
        "endpoint_cap": str(first_index_cap),
        "endpoint_margin": str(first_index_cap - first_index_error),
    }


def build_artifact() -> dict:
    source = load_source()
    exact = exact_arithmetic()
    rows = [
        TransferRow(
            id="l0fsdt_01_tail_expectation",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The complete-to-first-summand relative moment tail is an expectation of the pointwise theta-tail ratio.",
            formula=exact["relative_tail"],
            proof_boundary="Exact first-summand probability-coordinate identity.",
        ),
        TransferRow(
            id="l0fsdt_02_covariance_derivative",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Differentiation in the nonnegative heat-penalty parameter is a covariance identity.",
            formula=exact["covariance_derivative"],
            proof_boundary="Exact differentiation identity for finite first-summand moments.",
        ),
        TransferRow(
            id="l0fsdt_03_opposite_monotonicity",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="The decreasing theta-tail ratio and increasing square coordinate force the covariance to be nonpositive.",
            formula=exact["two_copy_identity"] + "; " + exact["opposite_monotonicity"],
            proof_boundary="Uses the already-proved pointwise monotonicity of every higher theta-summand ratio.",
        ),
        TransferRow(
            id="l0fsdt_04_uniform_heat_transfer",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="The lambda=-100 dominance theorem transfers uniformly through the complete heat interval to lambda zero.",
            formula=exact["uniform_dominance"],
            proof_boundary="Relative moment dominance only; no Hankel sign or curvature conclusion.",
            diagnostics={
                "source_status": source.get("status"),
                "source_row": "fsdc_08_full_moment_dominance",
            },
        ),
        TransferRow(
            id="l0fsdt_05_lambda_zero_dominance",
            role="exact_interval_composition",
            readiness="ready_to_apply",
            claim="At lambda zero, every complete moment from index 300 onward differs from its first theta summand by at most two inverse sixth powers relatively.",
            formula=exact["lambda_zero_dominance"],
            proof_boundary="All-index lambda-zero coefficient dominance; not a determinant-sign theorem.",
        ),
        TransferRow(
            id="l0fsdt_06_order4_penalty_transfer",
            role="exact_perturbation_theorem",
            readiness="ready_to_apply",
            claim="At the first unproved lambda-zero order-four index and beyond, the complete-kernel logarithmic penalty differs from the first-summand penalty by at most two-fifths inverse square.",
            formula=exact["penalty_error"],
            proof_boundary="Exact Lipschitz transfer only; the first-summand curvature ceiling remains open.",
            diagnostics={
                "splice": exact["splice"],
                "shifted_coefficients_descending": exact[
                    "perturbation_coefficients_descending"
                ],
                "endpoint_margin": exact["endpoint_margin"],
            },
        ),
        TransferRow(
            id="l0fsdt_07_curvature_handoff",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove a uniform lambda-zero first-summand stable-gap curvature ceiling from the real tent boundary t=503 onward.",
            formula="K_1(t)<=(7/2)/t^2 for every real t>=503",
            proof_boundary="Open first-summand analytic theorem; no lambda-zero all-shift order-four conclusion is asserted.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_lambda0_first_summand_dominance_transfer",
        "date": "2026-07-13",
        "status": "exact heat-parameter dominance transfer with open lambda-zero curvature tail",
        "proof_boundary": (
            "This artifact proves first-summand dominance uniformly for -100<=lambda<=0 "
            "and the resulting lambda-zero order-four perturbation budget. It does not "
            "prove the first-summand curvature ceiling, all-shift order-four positivity, "
            "forward order-four invariance, PF-infinity, RH, or Lambda<=0."
        ),
        "source": str(SOURCE.relative_to(REPO_ROOT)).replace("\\", "/"),
        "source_sha256": sha256(SOURCE),
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_lambda0_first_summand_dominance_transfer.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_lambda0_first_summand_dominance_transfer.py"
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": len(rows) - 1,
            "ready_to_apply_rows": len(rows) - 1,
            "open_analytic_rows": 1,
            "uniform_heat_intervals": 1,
            "lambda_zero_dominance_theorems": 1,
            "order4_penalty_transfers": 1,
            "tail_start_k": TAIL_START_K,
            "lambda0_splice_k": LAMBDA0_SPLICE_K,
        },
        "remaining_target": "K_1(t)<=(7/2)/t^2 for every real t>=503 at lambda=0",
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Lambda-Zero First-Summand Dominance Transfer",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact heat-parameter dominance transfer with an open lambda-zero",
        "first-summand curvature tail. This is not a proof of all-shift order-four",
        "positivity, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_lambda0_first_summand_dominance_transfer.json",
        "python work/rh_compute/scripts/jensen_window_pf_lambda0_first_summand_dominance_transfer.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_lambda0_first_summand_dominance_transfer.py",
        "```",
        "",
        "## Monotone Heat Transfer",
        "",
        "Write `T=-lambda`, so `0<=T<=100`, and normalize the first theta",
        "summand to the probability measure",
        "",
        "```text",
        exact["first_summand_measure"],
        "```",
        "",
        "For the decreasing pointwise higher-summand ratio `epsilon`,",
        "",
        "```text",
        exact["relative_tail"],
        exact["covariance_derivative"],
        exact["two_copy_identity"],
        "```",
        "",
        "Since `epsilon(u)` decreases while `u^2` increases, the two-copy",
        "integrand is nonpositive. Therefore",
        "",
        "```text",
        exact["heat_monotonicity"],
        exact["uniform_dominance"],
        "```",
        "",
        "The last inequality uses the certified all-index theorem at `T=100`.",
        "It follows in particular that",
        "",
        "```text",
        exact["lambda_zero_dominance"],
        "```",
        "",
        "## Order-Four Perturbation",
        "",
        "The direct Arb prefix ends at `n=500`; the first unproved tail index is",
        "`n=501`, hence `k=n+3=504`. The same exact Lipschitz calculation now gives",
        "",
        "```text",
        exact["moment_error"],
        exact["penalty_error"],
        "```",
        "",
        "After writing `k=504+n`, the final rational comparison has the",
        "coefficient-positive numerator",
        "",
        "```text",
        exact["perturbation_polynomial"] + ".",
        "```",
        "",
        "The remaining analytic target is now isolated as",
        "",
        "```text",
        artifact["remaining_target"] + ".",
        "```",
        "",
        "This transfer does not supply that curvature theorem or promote the",
        "lambda-zero finite prefix to an all-shift order-four theorem.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
        "outputs/arb_xi_lambda0_order4_prefix_certificate.md",
        "outputs/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "derived lambda-zero first-summand dominance transfer: "
        f"{summary['rows']} rows, {summary['exact_rows']} exact rows, "
        f"{summary['lambda_zero_dominance_theorems']} lambda-zero dominance theorem, "
        f"{summary['order4_penalty_transfers']} order-four penalty transfer, "
        f"{summary['open_analytic_rows']} open curvature tail"
    )


if __name__ == "__main__":
    main()
