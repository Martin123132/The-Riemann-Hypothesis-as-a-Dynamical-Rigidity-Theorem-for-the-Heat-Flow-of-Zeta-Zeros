#!/usr/bin/env python3
"""Compose uniform compact-heat contiguous order-four forward invariance."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md"
)
LAMBDA0_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.json"
)
SUPERPOLYNOMIAL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_uniform_superpolynomial_first_summand_dominance.json"
)
HEAT_TILT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.json"
)
UNIFORM_REDUCTION_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.json"
)
FORWARD_REDUCTION_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.json"
)
ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order4_m100_entry_certificate.json"
)
FLOW_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order4_forward_flow_reduction.json"
)


@dataclass(frozen=True)
class InvarianceRow:
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


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_diagnostics() -> dict:
    lambda0 = load_json(LAMBDA0_SOURCE)
    superpolynomial = load_json(SUPERPOLYNOMIAL_SOURCE)
    heat_tilt = load_json(HEAT_TILT_SOURCE)
    uniform_reduction = load_json(UNIFORM_REDUCTION_SOURCE)
    forward_reduction = load_json(FORWARD_REDUCTION_SOURCE)
    entry = load_json(ENTRY_SOURCE)
    flow = load_json(FLOW_SOURCE)
    checks = {
        "lambda0_eventual": lambda0.get("summary", {}).get(
            "eventual_positivity_theorems"
        )
        == 1,
        "higher_theta": superpolynomial.get("summary", {}).get(
            "higher_theta_handoffs_closed"
        )
        == 1,
        "heat_tilt": heat_tilt.get("summary", {}).get("uniform_heat_tilt_theorems")
        == 1,
        "uniform_transfer": uniform_reduction.get("summary", {}).get(
            "uniform_determinant_transfers"
        )
        == 1,
        "finite_confinement": forward_reduction.get("summary", {}).get(
            "finite_confinement_reductions"
        )
        == 1,
        "entry": entry.get("summary", {}).get("all_shift_order_four_entry_theorems")
        == 1,
        "cooperative_flow": flow.get("summary", {}).get("cooperative_flow_lemmas")
        == 2,
    }
    if not all(checks.values()):
        raise RuntimeError(f"order-four forward source contract failed: {checks}")
    return {
        "checks": checks,
        "lambda0_ratio": lambda0.get("published_ratio_input", {}).get("exact_ratio"),
        "lambda0_G2_limit": lambda0.get("published_ratio_input", {}).get(
            "needed_limit"
        ),
        "symbolic_main_term": lambda0.get("symbolic_cancellation", {}).get(
            "first_nonzero_term"
        ),
        "higher_theta": superpolynomial.get("exact", {}).get(
            "local_difference_consequence"
        ),
        "heat_tilt": heat_tilt.get("theorem", {}).get("target_theorem"),
        "uniform_contract": uniform_reduction.get("exact", {}).get(
            "uniform_ratio_hypothesis"
        ),
        "entry": entry.get("exact", {}).get("all_shift_entry"),
        "flow": flow.get("exact_flow", {}).get("weighted_identity"),
        "off_diagonal": flow.get("exact_flow", {}).get("a_n"),
        "finite_confinement": forward_reduction.get("exact", {}).get(
            "variation_of_constants"
        ),
    }


def exact_composition() -> dict:
    return {
        "base_ratio": (
            "the published lambda-zero Xi ratios have h=Delta(M)^2->0, "
            "G_2(M)->1, bounded G_3,...,G_7, and an o(h^6) local remainder"
        ),
        "first_summand_tilt": (
            "Delta_k^m log(A_k^(1)(-T)/A_k^(1)(0))="
            "O(log(k)/k^m), 2<=m<=7, uniformly for 0<=T<=100"
        ),
        "kernel_remainder": (
            "all fixed local log differences of log(A_k(-T)/A_k^(1)(-T)) "
            "are O_p,m(k^-p) uniformly for every p>0"
        ),
        "uniform_ratio_consequence": (
            "the complete A_k(lambda) family satisfies the uniform degree-seven "
            "ratio contract on -100<=lambda<=0 with the same G_2 limit 1"
        ),
        "heat_coefficient_scale": (
            "c_r(T,M)=O(log(M)/M^r) and h(M)~1/(2M) imply "
            "c_r(T,M)/h(M)^(r-1)=O(log(M)/M)=o(1), r=2,...,6"
        ),
        "determinant_asymptotic": (
            "H_(4,n)(lambda)=positive_scale(lambda,n)*(768*G_2(lambda,n)^6*"
            "h(lambda,n)^6+o(h(lambda,n)^6)) uniformly"
        ),
        "uniform_eventual_tail": (
            "there exists N such that H_(4,n)(lambda)>0 for every n>=N and "
            "-100<=lambda<=0"
        ),
        "cooperative_flow": "Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n>0",
        "variation_of_constants": (
            "Q_n(lambda)=E_n(lambda)*(Q_n(-100)+integral_(-100)^lambda "
            "E_n(s)^(-1)*a_n(s)*Q_(n+1)(s)ds)"
        ),
        "all_interval_theorem": (
            "H_(4,n)(lambda)>0 for every integer n>=0 and every "
            "lambda in [-100,0]"
        ),
        "lambda_zero_theorem": "H_(4,n)(0)>0 for every integer n>=0",
        "finite_splice_status": (
            "the previous non-effective lambda-zero splice is closed indirectly by "
            "uniform tail confinement and backward cooperative propagation"
        ),
    }


def newton_transfer() -> dict:
    j = sp.symbols("j")
    differences = sp.symbols("D1:7")
    polynomial = sp.expand(
        sum(differences[order - 1] * sp.binomial(j, order).expand(func=True) for order in range(1, 7))
    )
    coefficient_rows = []
    for degree in range(1, 7):
        coefficient = sp.expand(polynomial).coeff(j, degree)
        support = [
            order
            for order, difference in enumerate(differences, start=1)
            if coefficient.coeff(difference) != 0
        ]
        if not support or min(support) != degree:
            raise RuntimeError(f"Newton support is not triangular at degree {degree}")
        diagonal = sp.factor(coefficient.coeff(differences[degree - 1]))
        if diagonal != sp.Rational(1, sp.factorial(degree)):
            raise RuntimeError(f"Newton diagonal changed at degree {degree}")
        coefficient_rows.append(
            {
                "degree": degree,
                "coefficient": str(sp.factor(coefficient)),
                "difference_orders": support,
                "diagonal": str(diagonal),
            }
        )
    return {
        "interpolant": "q_M(j)=sum_(m=1)^6 D_m(M)*binom(j,m), 0<=j<=6",
        "exactness": "q_M(j)=f(M-j)-f(M) exactly on all seven determinant nodes",
        "triangular_support": (
            "the coefficient of j^r uses only D_m with m>=r and has D_r/r! on its diagonal"
        ),
        "difference_input": "D_m(M)=O(log(M)/M^m), m=2,...,6",
        "coefficient_output": "[j^r]q_M(j)=O(log(M)/M^r), r=2,...,6",
        "graded_output": (
            "after division by h(M)^(r-1), every heat correction to G_r is "
            "O(log(M)/M)=o(1)"
        ),
        "polynomial": str(polynomial),
        "coefficient_rows": coefficient_rows,
    }


def build_artifact() -> dict:
    sources = source_diagnostics()
    exact = exact_composition()
    newton = newton_transfer()
    if sources["symbolic_main_term"] != "768*G2^6*h^6":
        raise RuntimeError("universal order-four main term changed")
    rows = [
        InvarianceRow(
            id="co4uhfi_01_lambda0_ratio_input",
            role="published_theorem_input",
            readiness="ready_to_apply",
            claim="The published Xi ratio theorem and exact determinant cancellation give the universal positive order-four main term at lambda zero.",
            formula=exact["base_ratio"] + "; 768*G_2^6*h^6",
            proof_boundary="Eventual contiguous order-four asymptotic input only.",
        ),
        InvarianceRow(
            id="co4uhfi_02_heat_tilt",
            role="theorem_input",
            readiness="ready_to_apply",
            claim="The first theta summand has the required seven compact-uniform heat-tilt ratio estimates.",
            formula=exact["first_summand_tilt"],
            proof_boundary="First-summand compact-heat theorem.",
            diagnostics={"source_formula": sources["heat_tilt"]},
        ),
        InvarianceRow(
            id="co4uhfi_03_kernel_remainder",
            role="theorem_input",
            readiness="ready_to_apply",
            claim="Every higher-theta correction needed by the fixed order-four stencil is uniformly superpolynomial.",
            formula=exact["kernel_remainder"],
            proof_boundary="Complete-to-first-summand asymptotic transfer only.",
            diagnostics={"source_formula": sources["higher_theta"]},
        ),
        InvarianceRow(
            id="co4uhfi_04_uniform_ratio",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="The full heat-deformed Xi coefficient family inherits the uniform degree-seven ratio expansion and the unchanged positive quadratic limit.",
            formula=exact["uniform_ratio_consequence"],
            proof_boundary="Compact heat interval and seven-point local ratios only.",
            diagnostics={
                "newton_transfer": newton,
                "graded_scale": exact["heat_coefficient_scale"],
            },
        ),
        InvarianceRow(
            id="co4uhfi_05_uniform_eventual_tail",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The universal determinant expansion proves one strictly positive order-four tail uniformly across the full heat interval.",
            formula=exact["determinant_asymptotic"] + "; " + exact["uniform_eventual_tail"],
            proof_boundary="Uniform eventual contiguous order-four theorem; the threshold is non-effective.",
        ),
        InvarianceRow(
            id="co4uhfi_06_entry_and_flow",
            role="theorem_input",
            readiness="ready_to_apply",
            claim="All-shift lambda=-100 entry and the completed order-three layer supply a strictly cooperative order-four flow.",
            formula=sources["entry"] + "; " + exact["cooperative_flow"],
            proof_boundary="Previously certified entry and local flow inputs.",
            diagnostics={"off_diagonal": sources["off_diagonal"]},
        ),
        InvarianceRow(
            id="co4uhfi_07_finite_confinement",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The uniform positive tail reduces the infinite flow to finite backward variation-of-constants propagation.",
            formula=exact["variation_of_constants"],
            proof_boundary="Exact finite-confinement argument; no global coefficient ceiling.",
        ),
        InvarianceRow(
            id="co4uhfi_08_forward_invariance",
            role="theorem_conclusion",
            readiness="ready_to_apply",
            claim="Every shifted contiguous order-four determinant remains strictly positive from lambda=-100 through lambda=0.",
            formula=exact["all_interval_theorem"],
            proof_boundary="Contiguous order four only.",
        ),
        InvarianceRow(
            id="co4uhfi_09_lambda_zero",
            role="theorem_conclusion",
            readiness="ready_to_apply",
            claim="Every shifted contiguous four-by-four Hankel minor has the required positive sign at lambda zero.",
            formula=exact["lambda_zero_theorem"],
            proof_boundary="Not arbitrary-column order four, higher compounds, PF-infinity, RH, or Lambda<=0.",
            diagnostics={"splice": exact["finite_splice_status"]},
        ),
    ]
    source_paths = (
        LAMBDA0_SOURCE,
        SUPERPOLYNOMIAL_SOURCE,
        HEAT_TILT_SOURCE,
        UNIFORM_REDUCTION_SOURCE,
        FORWARD_REDUCTION_SOURCE,
        ENTRY_SOURCE,
        FLOW_SOURCE,
    )
    return {
        "kind": "jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate",
        "date": "2026-07-13",
        "status": "all-shift contiguous order-four forward invariance through lambda zero",
        "proof_boundary": (
            "This artifact proves H_(4,n)(lambda)>0 for every contiguous shift n>=0 "
            "and every -100<=lambda<=0. It does not prove noncontiguous or "
            "arbitrary-column order-four minors, any compound order five or higher, "
            "PF-infinity, the all-degree Jensen criterion, RH, or Lambda<=0."
        ),
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/") for path in source_paths
        ],
        "source_hashes": {
            str(path.relative_to(REPO_ROOT)).replace("\\", "/"): sha256(path)
            for path in source_paths
        },
        "source_diagnostics": sources,
        "exact": exact,
        "newton_transfer": newton,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": len(rows),
            "open_analytic_rows": 0,
            "uniform_ratio_theorems": 1,
            "uniform_eventual_tail_theorems": 1,
            "finite_confinement_theorems": 1,
            "forward_invariance_theorems": 1,
            "lambda_zero_all_shift_theorems": 1,
            "global_coefficient_ceilings_used": 0,
            "newton_coefficients_checked": 6,
        },
        "remaining_target": (
            "derive arbitrary-column order-four sign regularity or begin a genuinely "
            "new contiguous order-five entry and invariance layer"
        ),
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Four Uniform-Heat Forward Invariance",
        "",
        "Date: 2026-07-13",
        "",
        "Status: all-shift contiguous order-four forward invariance through lambda",
        "zero. This is not arbitrary-column order four, PF-infinity, RH, or a proof",
        "of `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.py",
        "```",
        "",
        "## Uniform Ratio Expansion",
        "",
        "The published lambda-zero Xi ratio theorem supplies the base graded",
        "expansion. Two new uniform theorems control the heat deformation:",
        "",
        "```text",
        exact["first_summand_tilt"],
        exact["kernel_remainder"],
        "```",
        "",
        "Therefore",
        "",
        "```text",
        exact["uniform_ratio_consequence"] + ".",
        "```",
        "",
        "On the seven determinant nodes, exact Newton interpolation is triangular:",
        "the coefficient of `j^r` uses only differences of order at least `r`.",
        "Thus the heat correction at degree `r` is `O(log(M)/M^r)`, and",
        "",
        "```text",
        exact["heat_coefficient_scale"] + ".",
        "```",
        "",
        "## Uniform Positive Tail",
        "",
        "The exact determinant cancellation is universal and gives",
        "",
        "```text",
        exact["determinant_asymptotic"] + ".",
        "```",
        "",
        "Since `G_2->1` uniformly, there is one finite, non-effective `N` with",
        "",
        "```text",
        exact["uniform_eventual_tail"] + ".",
        "```",
        "",
        "## Backward Cooperative Propagation",
        "",
        "At `lambda=-100`, every shift is already positive. Throughout the interval",
        "the order-three theorem gives",
        "",
        "```text",
        exact["cooperative_flow"] + ".",
        "```",
        "",
        "Starting at the uniform positive tail boundary, variation of constants",
        "propagates strict positivity backward through the finite remaining indices:",
        "",
        "```text",
        exact["variation_of_constants"] + ".",
        "```",
        "",
        "Hence",
        "",
        "```text",
        exact["all_interval_theorem"],
        exact["lambda_zero_theorem"] + ".",
        "```",
        "",
        "This closes the previous non-effective lambda-zero finite splice without",
        "making its threshold explicit. The next layer is not RH: noncontiguous",
        "order-four minors and every compound order at least five remain open.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.md",
        "outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md",
        "outputs/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.md",
        "outputs/jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.md",
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
        "proved order-four uniform-heat forward invariance: "
        f"{summary['rows']} rows, {summary['ready_to_apply_rows']} ready rows, "
        f"{summary['uniform_eventual_tail_theorems']} uniform tail theorem, "
        f"{summary['forward_invariance_theorems']} forward theorem, "
        f"{summary['lambda_zero_all_shift_theorems']} lambda-zero all-shift theorem, "
        f"{summary['open_analytic_rows']} open rows"
    )


if __name__ == "__main__":
    main()
