#!/usr/bin/env python3
"""Reduce order-four forward invariance to a uniform eventual positive tail."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.md"
)
ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order4_m100_entry_certificate.json"
)
FLOW_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order4_forward_flow_reduction.json"
)
ORDER3_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_compound_order3_forward_invariance_certificate.json"
)
LAMBDA0_EVENTUAL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.json"
)


@dataclass(frozen=True)
class ReductionRow:
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
    entry = load_json(ENTRY_SOURCE)
    flow = load_json(FLOW_SOURCE)
    order3 = load_json(ORDER3_SOURCE)
    eventual = load_json(LAMBDA0_EVENTUAL_SOURCE)
    if entry.get("summary", {}).get("all_shift_order_four_entry_theorems") != 1:
        raise RuntimeError("lambda=-100 order-four entry source is not closed")
    if flow.get("summary", {}).get("cooperative_flow_lemmas") != 2:
        raise RuntimeError("order-four cooperative flow source is not closed")
    if order3.get("summary", {}).get("full_forward_propagation_rows") != 1:
        raise RuntimeError("order-three forward source is not closed")
    if eventual.get("summary", {}).get("eventual_positivity_theorems") != 1:
        raise RuntimeError("lambda-zero eventual source is not closed")
    return {
        "entry_status": entry.get("status"),
        "entry_theorem": "H_(4,n)(-100)>0 for every n>=0",
        "flow_status": flow.get("status"),
        "flow_identity": flow.get("exact_flow", {}).get("weighted_identity"),
        "off_diagonal": flow.get("exact_flow", {}).get("a_n"),
        "order3_status": order3.get("status"),
        "order3_theorem": "H_(3,n)(lambda)<0 on -100<=lambda<=0",
        "lambda0_eventual_status": eventual.get("status"),
        "lambda0_eventual_theorem": "exists N_H4: H_(4,n)(0)>0 for n>=N_H4",
    }


def exact_reduction() -> dict:
    return {
        "cooperative_flow": "Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n>0",
        "uniform_tail_hypothesis": (
            "exists N such that Q_n(lambda)>0 for every n>=N and "
            "-100<=lambda<=0"
        ),
        "integrating_factor": (
            "E_n(lambda)=exp(integral_(-100)^lambda b_n(s)ds)>0"
        ),
        "variation_of_constants": (
            "Q_n(lambda)=E_n(lambda)*(Q_n(-100)+"
            "integral_(-100)^lambda E_n(s)^(-1)*a_n(s)*Q_(n+1)(s)ds)"
        ),
        "backward_induction": (
            "Q_N>0 on the interval and Q_n(-100)>0 imply successively "
            "Q_(N-1)>0,...,Q_0>0"
        ),
        "finite_confinement": (
            "uniform eventual positivity confines every possible first crossing "
            "to the finite set 0<=n<N"
        ),
        "conditional_forward_theorem": (
            "uniform eventual order-four positivity on [-100,0] implies "
            "H_(4,n)(lambda)>0 for every n>=0 and -100<=lambda<=0"
        ),
        "lambda_zero_consequence": (
            "uniform eventual order-four positivity on [-100,0] implies "
            "H_(4,n)(0)>0 for every n>=0"
        ),
        "uniform_asymptotic_target": (
            "prove one N works simultaneously for every lambda in [-100,0]"
        ),
        "endpoint_warning": (
            "eventual positivity at lambda=0 alone does not provide the uniform "
            "interior tail hypothesis"
        ),
    }


def build_artifact() -> dict:
    sources = source_diagnostics()
    exact = exact_reduction()
    rows = [
        ReductionRow(
            id="co4etfi_01_entry",
            role="theorem_input",
            readiness="ready_to_apply",
            claim="Every contiguous order-four determinant is strictly positive at the initial heat parameter.",
            formula=sources["entry_theorem"],
            proof_boundary="Previously certified all-shift entry at lambda=-100.",
        ),
        ReductionRow(
            id="co4etfi_02_cooperative_flow",
            role="exact_flow_input",
            readiness="ready_to_apply",
            claim="The completed order-three layer makes the order-four flow one-sided cooperative throughout the target heat interval.",
            formula=exact["cooperative_flow"],
            proof_boundary="Contiguous order-four flow only; the spatial tail is separate.",
            diagnostics={
                "off_diagonal": sources["off_diagonal"],
                "order3_theorem": sources["order3_theorem"],
            },
        ),
        ReductionRow(
            id="co4etfi_03_uniform_tail_target",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove that one eventual order-four tail is positive uniformly across the compact heat interval.",
            formula=exact["uniform_tail_hypothesis"],
            proof_boundary="Open uniform asymptotic theorem; endpoint eventual positivity is insufficient.",
            diagnostics={
                "available_endpoint": sources["lambda0_eventual_theorem"],
                "warning": exact["endpoint_warning"],
            },
        ),
        ReductionRow(
            id="co4etfi_04_finite_confinement",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="A uniform positive tail confines any possible loss of positivity to finitely many order-four coordinates.",
            formula=exact["finite_confinement"],
            proof_boundary="Exact consequence of the open uniform-tail hypothesis.",
        ),
        ReductionRow(
            id="co4etfi_05_variation_of_constants",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Each finite coordinate has an exact positive integrating-factor representation driven only by its next neighbor.",
            formula=exact["variation_of_constants"],
            proof_boundary="Finite-index linear ODE identity.",
            diagnostics={"integrating_factor": exact["integrating_factor"]},
        ),
        ReductionRow(
            id="co4etfi_06_backward_induction",
            role="exact_conditional_theorem",
            readiness="conditional",
            claim="Starting at the uniform positive tail boundary, variation of constants propagates strict positivity backward through every remaining finite index.",
            formula=exact["backward_induction"],
            proof_boundary="Conditional only on co4etfi_03; no coefficient ceiling is required.",
        ),
        ReductionRow(
            id="co4etfi_07_forward_invariance",
            role="exact_conditional_theorem",
            readiness="conditional",
            claim="Uniform eventual positivity closes contiguous order-four forward invariance from lambda=-100 through lambda=0.",
            formula=exact["conditional_forward_theorem"],
            proof_boundary="Conditional on the single uniform-tail theorem; contiguous order four only.",
        ),
        ReductionRow(
            id="co4etfi_08_lambda_zero_handoff",
            role="exact_conditional_consequence",
            readiness="conditional",
            claim="The same hypothesis promotes the existing finite lambda-zero prefix to an all-shift contiguous order-four theorem.",
            formula=exact["lambda_zero_consequence"],
            proof_boundary="Conditional; not arbitrary-column order four, PF-infinity, RH, or Lambda<=0.",
        ),
    ]
    source_paths = (ENTRY_SOURCE, FLOW_SOURCE, ORDER3_SOURCE, LAMBDA0_EVENTUAL_SOURCE)
    return {
        "kind": "jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction",
        "date": "2026-07-13",
        "status": "exact finite-confinement forward-invariance reduction with one open uniform tail",
        "proof_boundary": (
            "This artifact proves that a uniform eventual positive order-four tail on "
            "[-100,0] is sufficient for full contiguous order-four forward invariance, "
            "without a global coefficient ceiling. It does not prove that uniform tail, "
            "unconditional order-four invariance, arbitrary-column order four, "
            "PF-infinity, RH, or Lambda<=0."
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
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_or_input_rows": len(rows) - 1,
            "ready_to_apply_rows": 4,
            "conditional_rows": 3,
            "open_analytic_rows": 1,
            "finite_confinement_reductions": 1,
            "variation_of_constants_identities": 1,
            "conditional_forward_theorems": 1,
            "global_coefficient_ceilings_required": 0,
        },
        "remaining_target": exact["uniform_tail_hypothesis"],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    sources = artifact["source_diagnostics"]
    lines = [
        "# Jensen-Window PF Compound Order-Four Eventual-Tail Forward-Invariance Reduction",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact finite-confinement forward-invariance reduction with one",
        "open uniform tail. This is not a proof of unconditional order-four",
        "invariance, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.py",
        "```",
        "",
        "## Available Inputs",
        "",
        "The completed entry and order-three propagation theorems give",
        "",
        "```text",
        sources["entry_theorem"],
        sources["order3_theorem"],
        exact["cooperative_flow"],
        "```",
        "",
        "## Alternative Tail Hypothesis",
        "",
        "Instead of bounding every effective diagonal coefficient, assume only",
        "",
        "```text",
        exact["uniform_tail_hypothesis"] + ".",
        "```",
        "",
        "The endpoint theorem currently proves only",
        "",
        "```text",
        sources["lambda0_eventual_theorem"] + ".",
        "```",
        "",
        "That endpoint statement does not yet control the interior heat interval.",
        "",
        "## Finite Confinement",
        "",
        "For each fixed `n`, put",
        "",
        "```text",
        exact["integrating_factor"],
        "```",
        "",
        "Variation of constants gives",
        "",
        "```text",
        exact["variation_of_constants"] + ".",
        "```",
        "",
        "Under the uniform tail hypothesis, `Q_N` is positive throughout the",
        "interval. Since `Q_n(-100)>0`, `E_n>0`, and `a_n>0`, backward induction",
        "from `N-1` to `0` proves",
        "",
        "```text",
        exact["conditional_forward_theorem"] + ".",
        "```",
        "",
        "This route needs no supremum bound on `beta_n+alpha_n*(n+2)/(n+1)`.",
        "Its sole open input is a compact-heat, uniform version of the available",
        "lambda-zero eventual asymptotic theorem.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_m100_entry_certificate.md",
        "outputs/jensen_window_pf_compound_order4_forward_flow_reduction.md",
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
        "derived order-four eventual-tail forward-invariance reduction: "
        f"{summary['rows']} rows, {summary['exact_or_input_rows']} exact/input rows, "
        f"{summary['finite_confinement_reductions']} finite-confinement reduction, "
        f"{summary['conditional_forward_theorems']} conditional forward theorem, "
        f"{summary['open_analytic_rows']} open uniform tail"
    )


if __name__ == "__main__":
    main()
