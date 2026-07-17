#!/usr/bin/env python3
"""Reduce the compact-heat order-four tail to a uniform ratio expansion."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate import (  # noqa: E402
    determinant_cancellation,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.md"
)
LAMBDA0_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.json"
)
FORWARD_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.json"
)
DOMINANCE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_lambda0_first_summand_dominance_transfer.json"
)
SUPERPOLYNOMIAL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_uniform_superpolynomial_first_summand_dominance.json"
)


@dataclass(frozen=True)
class UniformRow:
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


def source_diagnostics() -> dict:
    lambda0 = json.loads(LAMBDA0_SOURCE.read_text(encoding="utf-8"))
    forward = json.loads(FORWARD_SOURCE.read_text(encoding="utf-8"))
    dominance = json.loads(DOMINANCE_SOURCE.read_text(encoding="utf-8"))
    superpolynomial = json.loads(
        SUPERPOLYNOMIAL_SOURCE.read_text(encoding="utf-8")
    )
    if lambda0.get("summary", {}).get("eventual_positivity_theorems") != 1:
        raise RuntimeError("lambda-zero eventual theorem source is not closed")
    if forward.get("summary", {}).get("finite_confinement_reductions") != 1:
        raise RuntimeError("finite-confinement source is not closed")
    if dominance.get("summary", {}).get("lambda_zero_dominance_theorems") != 1:
        raise RuntimeError("first-summand dominance transfer source is not closed")
    if superpolynomial.get("summary", {}).get("higher_theta_handoffs_closed") != 1:
        raise RuntimeError("superpolynomial higher-theta source is not closed")
    return {
        "lambda0_status": lambda0.get("status"),
        "lambda0_ratio_source": lambda0.get("published_ratio_input", {}).get("source"),
        "lambda0_eventual": lambda0.get("eventual_theorem", {}).get("strict_sign"),
        "forward_status": forward.get("status"),
        "forward_remaining_target": forward.get("remaining_target"),
        "dominance_status": dominance.get("status"),
        "dominance_uniform_interval": dominance.get("exact", {}).get(
            "uniform_dominance"
        ),
        "superpolynomial_status": superpolynomial.get("status"),
        "superpolynomial_local_differences": superpolynomial.get("exact", {}).get(
            "local_difference_consequence"
        ),
    }


def exact_reduction(cancellation: dict) -> dict:
    if cancellation.get("first_nonzero_term") != "768*G2^6*h^6":
        raise RuntimeError("universal determinant cancellation changed")
    return {
        "parameter_interval": "I=[-100,0]",
        "indexing": "M=n+6; backward shifts j=6-i-j_column in {0,...,6}",
        "uniform_ratio_hypothesis": (
            "log(A_(M-j)(lambda)/A_M(lambda))="
            "-G_1(lambda,M)*j-sum_(m=2)^7 G_m(lambda,M)*h(lambda,M)^(m-1)*j^m"
            "+r_(lambda,M,j)"
        ),
        "uniform_small_parameter": (
            "sup_(lambda in I) h(lambda,M)->0 and "
            "sup_(lambda in I)|G_2(lambda,M)-1|->0"
        ),
        "uniform_coefficient_control": (
            "sup_(lambda in I)|G_m(lambda,M)|=O(1), m=3,...,7"
        ),
        "uniform_remainder": (
            "max_(lambda in I,0<=j<=6)|r_(lambda,M,j)|=o(h(lambda,M)^6) uniformly"
        ),
        "normalized_matrix": (
            "K_(i,j)=exp(-sum_(m=2)^7 G_m*h^(m-1)*(6-i-j)^m+r_(6-i-j))"
        ),
        "determinant_expansion": (
            "det K=768*G_2(lambda,M)^6*h(lambda,M)^6+o(h(lambda,M)^6) uniformly"
        ),
        "uniform_tail_consequence": (
            "there exists N such that H_(4,n)(lambda)>0 for every n>=N and "
            "lambda in [-100,0]"
        ),
        "forward_consequence": (
            "the uniform tail plus lambda=-100 entry and the cooperative flow imply "
            "H_(4,n)(lambda)>0 for all n>=0 and lambda in [-100,0]"
        ),
        "heat_tilt_factor": (
            "R_T^(1)(k)=A_k^(1)(-T)/A_k^(1)(0), 0<=T<=100"
        ),
        "sufficient_heat_tilt_estimate": (
            "Delta_k^m log R_T^(1)(k)=O(log(k)/k^m) uniformly for 0<=T<=100, "
            "m=2,...,7"
        ),
        "why_sufficient": (
            "O(log(k)/k^m)=o(k^(-(m-1))) for m>=2, so the bounded heat tilt "
            "does not change the universal G_2 limit or the graded ratio scales"
        ),
        "higher_theta_handoff": (
            "closed: every fixed local log difference of the complete-to-first-"
            "summand correction is uniformly O_p,m(k^-p) for every p>0"
        ),
        "symbolic_coefficients": cancellation["coefficients_h0_through_h6"],
        "permutations_checked": cancellation["permutations_checked"],
    }


def build_artifact() -> dict:
    sources = source_diagnostics()
    cancellation = determinant_cancellation()
    exact = exact_reduction(cancellation)
    rows = [
        UniformRow(
            id="co4uhet_01_uniform_ratio_contract",
            role="exact_hypothesis_contract",
            readiness="ready_to_apply",
            claim="A compact-parameter degree-seven backward-ratio expansion isolates all analytic input needed by the order-four determinant.",
            formula=exact["uniform_ratio_hypothesis"],
            proof_boundary="Sufficient hypothesis contract only; it is not yet proved for the heat-deformed Xi coefficients.",
        ),
        UniformRow(
            id="co4uhet_02_symbolic_cancellation",
            role="exact_symbolic_lemma",
            readiness="ready_to_apply",
            claim="Uniformity in the heat parameter does not alter the universal cancellation through order h^5 or its positive first coefficient.",
            formula="[h^0,...,h^6] det K=[0,0,0,0,0,0,768*G_2^6]",
            proof_boundary="Exact 24-permutation formal determinant computation.",
            diagnostics={
                "coefficients": exact["symbolic_coefficients"],
                "permutations_checked": exact["permutations_checked"],
            },
        ),
        UniformRow(
            id="co4uhet_03_uniform_determinant_transfer",
            role="exact_conditional_theorem",
            readiness="conditional",
            claim="The uniform ratio contract forces a strictly positive contiguous order-four tail simultaneously across the complete compact heat interval.",
            formula=exact["determinant_expansion"] + "; " + exact["uniform_tail_consequence"],
            proof_boundary="Conditional on co4uhet_01 for the actual heat-deformed Xi coefficients.",
        ),
        UniformRow(
            id="co4uhet_04_finite_confinement_composition",
            role="exact_conditional_theorem",
            readiness="conditional",
            claim="The uniform determinant tail composes with entry and cooperativity to close every finite order-four coordinate.",
            formula=exact["forward_consequence"],
            proof_boundary="Conditional contiguous order-four theorem only.",
            diagnostics={"source": sources["forward_status"]},
        ),
        UniformRow(
            id="co4uhet_05_heat_tilt_reduction",
            role="exact_asymptotic_reduction",
            readiness="ready_to_apply",
            claim="It is sufficient to control seven local finite differences of the bounded heat-tilt multiplier relative to the published lambda-zero sequence.",
            formula=exact["sufficient_heat_tilt_estimate"],
            proof_boundary="Sufficient asymptotic estimate; no estimate is asserted here.",
            diagnostics={
                "tilt": exact["heat_tilt_factor"],
                "scale_comparison": exact["why_sufficient"],
            },
        ),
        UniformRow(
            id="co4uhet_06_first_summand_input",
            role="theorem_input",
            readiness="ready_to_apply",
            claim="The complete-to-first-summand correction and every fixed local logarithmic difference are uniformly superpolynomial throughout the heat interval.",
            formula=sources["superpolynomial_local_differences"],
            proof_boundary="Closed higher-theta handoff; only the first-summand heat tilt remains.",
            diagnostics={
                "baseline_power_six": sources["dominance_uniform_interval"],
                "status": sources["superpolynomial_status"],
            },
        ),
        UniformRow(
            id="co4uhet_07_uniform_heat_tilt_target",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove the uniform first-summand heat-tilt finite-difference estimate.",
            formula=exact["sufficient_heat_tilt_estimate"],
            proof_boundary="Sole Xi-specific asymptotic target in this reduction; no forward invariance or RH conclusion is asserted.",
            diagnostics={"higher_theta": exact["higher_theta_handoff"]},
        ),
    ]
    source_paths = (
        LAMBDA0_SOURCE,
        FORWARD_SOURCE,
        DOMINANCE_SOURCE,
        SUPERPOLYNOMIAL_SOURCE,
    )
    return {
        "kind": "jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction",
        "date": "2026-07-13",
        "status": "exact uniform-asymptotic order-four reduction with one open heat-tilt theorem",
        "proof_boundary": (
            "This artifact proves the universal determinant and finite-confinement "
            "implications of a uniform compact-heat ratio expansion. It does not prove "
            "that Xi-specific expansion, a uniform eventual tail, unconditional "
            "order-four invariance, arbitrary-column order four, PF-infinity, RH, or "
            "Lambda<=0."
        ),
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/") for path in source_paths
        ],
        "source_hashes": {
            str(path.relative_to(REPO_ROOT)).replace("\\", "/"): sha256(path)
            for path in source_paths
        },
        "source_diagnostics": sources,
        "symbolic_cancellation": cancellation,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_or_input_rows": len(rows) - 1,
            "ready_to_apply_rows": 4,
            "conditional_rows": 2,
            "open_analytic_rows": 1,
            "symbolic_coefficients_checked": 7,
            "uniform_determinant_transfers": 1,
            "finite_confinement_compositions": 1,
            "open_heat_tilt_theorems": 1,
            "higher_theta_handoffs_closed": 1,
        },
        "remaining_target": exact["sufficient_heat_tilt_estimate"],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Four Uniform-Heat Eventual-Tail Reduction",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact uniform-asymptotic order-four reduction with one open",
        "heat-tilt theorem. This is not a proof of unconditional order-four",
        "invariance, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.py",
        "```",
        "",
        "## Uniform Ratio Contract",
        "",
        "For `I=[-100,0]`, `M=n+6`, and backward shifts `0<=j<=6`, it is",
        "sufficient to prove uniformly in `lambda` that",
        "",
        "```text",
        exact["uniform_ratio_hypothesis"],
        exact["uniform_small_parameter"],
        exact["uniform_coefficient_control"],
        exact["uniform_remainder"],
        "```",
        "",
        "## Universal Determinant",
        "",
        "The exact 24-permutation calculation is parameter-blind:",
        "",
        "```text",
        "[h^0,...,h^6] det K=[0,0,0,0,0,0,768*G_2^6].",
        exact["determinant_expansion"] + ".",
        "```",
        "",
        "Since `G_2->1` uniformly, this would prove",
        "",
        "```text",
        exact["uniform_tail_consequence"] + ".",
        "```",
        "",
        "The finite-confinement theorem would then give",
        "",
        "```text",
        exact["forward_consequence"] + ".",
        "```",
        "",
        "## Xi-Specific Target",
        "",
        "Relative to the published lambda-zero coefficient sequence, set",
        "",
        "```text",
        exact["heat_tilt_factor"] + ".",
        "```",
        "",
        "A sufficient compact-heat estimate is",
        "",
        "```text",
        exact["sufficient_heat_tilt_estimate"] + ".",
        "```",
        "",
        exact["why_sufficient"] + ".",
        "The superpolynomial dominance theorem now closes the complete-kernel",
        "correction and all seven required local log differences. The only remaining",
        "asymptotic input is the displayed first-summand heat-tilt estimate.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.md",
        "outputs/jensen_window_pf_lambda0_first_summand_dominance_transfer.md",
        "outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md",
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
        "derived order-four uniform-heat eventual-tail reduction: "
        f"{summary['rows']} rows, {summary['exact_or_input_rows']} exact/input rows, "
        f"{summary['symbolic_coefficients_checked']} symbolic coefficients, "
        f"{summary['uniform_determinant_transfers']} conditional uniform transfer, "
        f"{summary['open_heat_tilt_theorems']} open heat-tilt theorem"
    )


if __name__ == "__main__":
    main()
