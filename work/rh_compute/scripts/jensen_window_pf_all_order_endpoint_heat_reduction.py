#!/usr/bin/env python3
"""Prove the all-order endpoint-to-heat reduction for signed Hankel layers."""

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
    / "work/rh_compute/results/jensen_window_pf_all_order_endpoint_heat_reduction.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_all_order_endpoint_heat_reduction.md"

TAIL_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.json"
)
FLOW_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_uniform_tail_flow_reduction.json"
)
ORDER9_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate.json"
)
TRANSFER_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_order4_noncontiguous_total_positivity_transfer.json"
)
ORDER10_COUNTEREXAMPLE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_endpoint_order10_counterexample.json"
)
SOURCE_PATHS = (
    TAIL_SOURCE,
    FLOW_SOURCE,
    ORDER9_SOURCE,
    TRANSFER_SOURCE,
    ORDER10_COUNTEREXAMPLE_SOURCE,
)

SYMBOLIC_ORDERS = tuple(range(2, 6))
PARITY_MAX_ORDER = 256


@dataclass(frozen=True)
class ReductionRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_record(path: Path, payload: dict) -> dict:
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "sha256": sha256(path),
        "kind": payload.get("kind"),
        "status": payload.get("status"),
    }


def validate_sources() -> tuple[list[dict], dict]:
    tail = load_json(TAIL_SOURCE)
    flow = load_json(FLOW_SOURCE)
    order9 = load_json(ORDER9_SOURCE)
    transfer = load_json(TRANSFER_SOURCE)
    order10 = load_json(ORDER10_COUNTEREXAMPLE_SOURCE)

    if tail.get("summary", {}).get("all_fixed_order_eventual_tail_theorems") != 1:
        raise RuntimeError("all-fixed-order eventual-tail source is not closed")
    if flow.get("summary", {}).get("signed_condensation_recursions") != 1:
        raise RuntimeError("general signed-condensation source is not closed")
    if flow.get("summary", {}).get("cooperative_flow_recursions") != 1:
        raise RuntimeError("general cooperative-flow source is not closed")
    if order9.get("summary", {}).get("contiguous_order_nine_theorems") != 1:
        raise RuntimeError("contiguous order-nine source is not closed")
    if order9.get("summary", {}).get("arbitrary_column_order_nine_theorems") != 1:
        raise RuntimeError("arbitrary-column order-nine source is not closed")
    if transfer.get("summary", {}).get("fixed_order_transfer_theorems") != 1:
        raise RuntimeError("fixed-order initial-minor transfer source is not closed")
    if order10.get("summary", {}).get("negative_deep_rectangles") != 4:
        raise RuntimeError("order-ten endpoint counterexample source changed")
    if order10.get("summary", {}).get("rejected_all_order_endpoint_hierarchies") != 1:
        raise RuntimeError("all-order endpoint hierarchy is not marked rejected")

    flow_exact = flow.get("exact_flow", {})
    required_flow = {
        "signed_condensation_recursion": (
            "Q_(m,n)*Q_(m-2,n+2)=Q_(m-1,n+1)^2-"
            "Q_(m-1,n)*Q_(m-1,n+2)"
        ),
        "general_affine_derivative": (
            "Q_(m,n)'=c_(m,n)*delta(Q_(m,n)), c_(m,n)=4*n+8*m-6"
        ),
        "general_plucker": (
            "Q_(m-1,n+1)*delta(Q_(m,n))=Q_(m-1,n)*Q_(m,n+1)+"
            "delta(Q_(m-1,n+1))*Q_(m,n)"
        ),
    }
    for key, expected in required_flow.items():
        if flow_exact.get(key) != expected:
            raise RuntimeError(f"general flow source changed at {key}")

    payloads = (tail, flow, order9, transfer, order10)
    records = [
        source_record(path, payload) for path, payload in zip(SOURCE_PATHS, payloads)
    ]
    contract = {
        "fixed_order_tail": (
            "for every fixed m there exists N_m such that Q_(m,n)(lambda)>0 "
            "for n>=N_m and -100<=lambda<=0"
        ),
        "general_condensation": required_flow["signed_condensation_recursion"],
        "general_affine_derivative": required_flow["general_affine_derivative"],
        "general_flag_plucker": required_flow["general_plucker"],
        "completed_base": (
            "Q_(m,n)(lambda)>0 for 1<=m<=9, n>=0, and -100<=lambda<=0"
        ),
        "fixed_order_column_transfer": (
            "completed signed contiguous layers through each fixed m imply all "
            "consecutive-row arbitrary-column signed minors through m"
        ),
        "order10_counterexample": (
            "Q_(10,n)(-100)<0 for n=0,1,2,3"
        ),
    }
    return records, contract


def hankel(values: tuple[sp.Symbol, ...], order: int, shift: int) -> sp.Expr:
    if order == 0:
        return sp.Integer(1)
    return sp.Matrix(
        [[values[shift + i + j] for j in range(order)] for i in range(order)]
    ).det(method="domain-ge")


def delta_derivative(expression: sp.Expr, values: tuple[sp.Symbol, ...]) -> sp.Expr:
    return sp.expand(
        sum(
            sp.diff(expression, values[index]) * values[index + 1]
            for index in range(len(values) - 1)
        )
    )


def heat_derivative(
    expression: sp.Expr, values: tuple[sp.Symbol, ...], n: sp.Symbol
) -> sp.Expr:
    return sp.expand(
        sum(
            sp.diff(expression, values[index])
            * (4 * (n + index) + 2)
            * values[index + 1]
            for index in range(len(values) - 1)
        )
    )


def symbolic_identity_audit() -> dict:
    n = sp.symbols("n", integer=True, nonnegative=True)
    rows = []
    for order in SYMBOLIC_ORDERS:
        values = sp.symbols(f"a0:{2 * order + 2}")
        hm = hankel(values, order, 0)
        hm_next = hankel(values, order, 1)
        lower = hankel(values, order - 1, 0)
        lower_mid = hankel(values, order - 1, 1)
        lower_next = hankel(values, order - 1, 2)
        lower_two = hankel(values, order - 2, 2)

        condensation = sp.factor(
            hm * lower_two - (lower * lower_next - lower_mid**2)
        )
        plucker = sp.factor(
            lower_mid * delta_derivative(hm, values)
            - lower * hm_next
            - delta_derivative(lower_mid, values) * hm
        )
        affine = sp.factor(
            heat_derivative(hm, values, n)
            - (4 * n + 8 * order - 6) * delta_derivative(hm, values)
        )
        residuals = {
            "condensation": str(condensation),
            "flag_plucker": str(plucker),
            "affine_heat": str(affine),
        }
        if any(value != "0" for value in residuals.values()):
            raise RuntimeError(f"symbolic identity audit failed at order {order}")
        rows.append({"order": order, "residuals": residuals})
    return {
        "orders": list(SYMBOLIC_ORDERS),
        "rows": rows,
        "all_residuals_zero": True,
        "scope": (
            "Independent finite symbolic specializations; arbitrary-order validity "
            "comes from the determinant cancellation, Desnanot-Jacobi, and flag "
            "Plucker derivations recorded in the theorem rows."
        ),
    }


def parity_and_coefficient_audit() -> dict:
    parity_rows = []
    for order in range(2, PARITY_MAX_ORDER + 1):
        exponent_residual = (
            order * (order - 1) // 2
            + (order - 2) * (order - 3) // 2
            - 2 * ((order - 1) * (order - 2) // 2)
        )
        if exponent_residual != 1:
            raise RuntimeError(f"orientation parity failed at order {order}")
        parity_rows.append({"order": order, "exponent_residual": exponent_residual})

    m, n = sp.symbols("m n", integer=True, positive=True)
    current = 4 * n + 8 * m - 6
    lower_shifted = 4 * (n + 1) + 8 * (m - 1) - 6
    coefficient_residual = sp.expand(lower_shifted - (current - 4))
    if coefficient_residual != 0:
        raise RuntimeError("cooperative coefficient shift failed")
    return {
        "parity_formula": (
            "binom(m,2)+binom(m-2,2)-2*binom(m-1,2)=1"
        ),
        "parity_checks": len(parity_rows),
        "parity_min_order": 2,
        "parity_max_order": PARITY_MAX_ORDER,
        "parity_residual_values": sorted(
            {row["exponent_residual"] for row in parity_rows}
        ),
        "coefficient_formula": "c_(m-1,n+1)=c_(m,n)-4",
        "coefficient_residual": str(coefficient_residual),
    }


def exact_statements() -> dict:
    return {
        "definitions": (
            "H_(m,n)=det[A_(n+i+j)]_(0<=i,j<m), H_(0,n)=1, "
            "epsilon_m=(-1)^binom(m,2), Q_(m,n)=epsilon_m*H_(m,n)"
        ),
        "coefficient_heat_flow": "A_n'=(4*n+2)*A_(n+1)",
        "delta_border": (
            "delta(H_(m,n))=det[C_0,...,C_(m-2),C_m], "
            "C_j=(A_(n+i+j))_(0<=i<m)"
        ),
        "affine_heat": (
            "Q_(m,n)'=c_(m,n)*delta(Q_(m,n)), c_(m,n)=4*n+8*m-6"
        ),
        "affine_cancellation": (
            "column part=[4*(n+m-1)+2]*delta(H); "
            "row part=4*(m-1)*delta(H)"
        ),
        "signed_condensation": (
            "Q_(m,n)*Q_(m-2,n+2)=Q_(m-1,n+1)^2-"
            "Q_(m-1,n)*Q_(m-1,n+2)"
        ),
        "endpoint_coordinate": (
            "Q_(m,n)>0 iff Q_(m-1,n+1)^2>"
            "Q_(m-1,n)*Q_(m-1,n+2), provided Q_(m-2,n+2)>0"
        ),
        "flag_plucker": (
            "Q_(m-1,n+1)*delta(Q_(m,n))=Q_(m-1,n)*Q_(m,n+1)+"
            "delta(Q_(m-1,n+1))*Q_(m,n)"
        ),
        "cooperative_flow": (
            "Q_(m,n)'=a_(m,n)*Q_(m,n+1)+b_(m,n)*Q_(m,n), "
            "a_(m,n)=c_(m,n)*Q_(m-1,n)/Q_(m-1,n+1)>0, "
            "b_(m,n)=c_(m,n)/(c_(m,n)-4)*(log Q_(m-1,n+1))'"
        ),
        "fixed_order_tail": (
            "for every fixed m exists N_m such that for every n>=N_m and "
            "-100<=lambda<=0, Q_(m,n)(lambda)>0"
        ),
        "variation_of_constants": (
            "Q_(m,n)(lambda)=E_(m,n)(lambda)*(Q_(m,n)(-100)+"
            "integral_(-100)^lambda E_(m,n)(s)^(-1)*a_(m,n)(s)*"
            "Q_(m,n+1)(s)ds), E_(m,n)>0"
        ),
        "single_layer_implication": (
            "[Q_(k,n)>0 on the heat interval for every k<m,n and "
            "Q_(m,n)(-100)>0 for every n] => "
            "[Q_(m,n)(lambda)>0 for every n and -100<=lambda<=0]"
        ),
        "known_base": (
            "Q_(m,n)(lambda)>0 for 1<=m<=9, every n>=0, and "
            "-100<=lambda<=0"
        ),
        "candidate_endpoint": (
            "Q_(m,n)(-100)>0 for every integer m>=10 and n>=0"
        ),
        "order10_counterexample": "Q_(10,n)(-100)<0 for n=0,1,2,3",
        "all_order_interval": (
            "Q_(m,n)(lambda)>0 for every integer m>=1,n>=0 and "
            "-100<=lambda<=0"
        ),
        "endpoint_interval_equivalence": (
            "[Q_(m,n)(-100)>0 for every m>=10,n>=0] iff "
            "[Q_(m,n)(lambda)>0 for every m>=10,n>=0,-100<=lambda<=0]"
        ),
        "arbitrary_column_consequence": (
            "all-order contiguous interval positivity => all-order "
            "consecutive-row arbitrary-column signed-Hankel positivity"
        ),
        "bridge_boundary": (
            "all-order shifted signed-Hankel positivity does not by itself prove "
            "PF-infinity of every binomially weighted Jensen window"
        ),
    }


def build_artifact() -> dict:
    sources, source_contract = validate_sources()
    symbolic = symbolic_identity_audit()
    arithmetic = parity_and_coefficient_audit()
    exact = exact_statements()

    rows = [
        ReductionRow(
            "jwpfaoehr_01_coordinate",
            "exact_definition",
            "ready_to_apply",
            "Use one signed coordinate at every compound order, with the empty minor normalized to one.",
            exact["definitions"],
            "Definition only; no sign is asserted.",
        ),
        ReductionRow(
            "jwpfaoehr_02_fixed_order_tail",
            "theorem_input",
            "ready_to_apply",
            "Every fixed order has a heat-uniform eventual positive signed tail, although its threshold may depend on the order.",
            exact["fixed_order_tail"],
            "The quantifier is for every m there exists N_m, not one N uniform in m.",
            {"source": sources[0]},
        ),
        ReductionRow(
            "jwpfaoehr_03_signed_condensation",
            "exact_identity",
            "ready_to_apply",
            "Desnanot-Jacobi converts each new endpoint layer into strict log-concavity of the preceding positive layer.",
            exact["signed_condensation"] + "; " + exact["endpoint_coordinate"],
            "The sign follows from one exact orientation-parity flip.",
            {
                "parity_formula": arithmetic["parity_formula"],
                "parity_checks": arithmetic["parity_checks"],
            },
        ),
        ReductionRow(
            "jwpfaoehr_04_affine_heat",
            "exact_all_order_identity",
            "ready_to_apply",
            "Multilinearity and duplicate-column/row cancellation give the affine Hankel heat derivative at arbitrary order.",
            exact["coefficient_heat_flow"] + "; " + exact["affine_heat"],
            "Only the final shifted column and final shifted row survive.",
            {"cancellation": exact["affine_cancellation"]},
        ),
        ReductionRow(
            "jwpfaoehr_05_flag_plucker",
            "exact_all_order_identity",
            "ready_to_apply",
            "The adjacent flag-Plucker relation converts the affine derivative into a nearest-neighbor system.",
            exact["flag_plucker"],
            "Apply the three-term flag relation to C_0,...,C_m and their first m-1 rows.",
            symbolic,
        ),
        ReductionRow(
            "jwpfaoehr_06_cooperative_flow",
            "exact_all_order_lemma",
            "ready_to_apply",
            "Whenever the preceding signed layer is positive, the order-m flow is cooperative in the forward shift direction.",
            exact["cooperative_flow"],
            "For m>=2; c_(m,n)>0 and the lower-layer quotient is positive.",
            {
                "coefficient_formula": arithmetic["coefficient_formula"],
                "coefficient_residual": arithmetic["coefficient_residual"],
            },
        ),
        ReductionRow(
            "jwpfaoehr_07_single_layer_propagation",
            "exact_quantifier_lemma",
            "ready_to_apply",
            "At one fixed order, endpoint entry and the order-dependent eventual tail propagate positivity to the whole heat interval.",
            exact["single_layer_implication"] + "; " + exact["variation_of_constants"],
            "Fix m first, choose N_m, and descend only through the finite shifts N_m-1,...,0.",
        ),
        ReductionRow(
            "jwpfaoehr_08_completed_base",
            "theorem_input",
            "ready_to_apply",
            "The rigorous fixed-order programme supplies the induction base through order nine.",
            exact["known_base"],
            "No order at least ten is imported.",
            {"source": sources[2]},
        ),
        ReductionRow(
            "jwpfaoehr_09_endpoint_equivalence",
            "exact_all_order_reduction",
            "ready_to_apply",
            "Ordinary induction over m makes the remaining all-order heat theorem equivalent to one static endpoint hierarchy.",
            exact["endpoint_interval_equivalence"],
            "The reverse implication is restriction to lambda=-100; the forward implication iterates the single-layer lemma.",
        ),
        ReductionRow(
            "jwpfaoehr_10_arbitrary_columns",
            "exact_conditional_consequence",
            "ready_to_apply",
            "After contiguous completion at every finite order, the initial-minor theorem transfers the result to arbitrary columns pointwise in lambda.",
            exact["arbitrary_column_consequence"],
            "For each finite m separately; no infinite determinant is invoked.",
            {"source": sources[3]},
        ),
        ReductionRow(
            "jwpfaoehr_11_static_hierarchy",
            "exact_target_coordinate",
            "ready_to_apply",
            "The conditional all-order theorem is equivalently an infinite strict log-concavity hierarchy in the positive signed coordinates.",
            exact["endpoint_coordinate"],
            "This is an exact reformulation of the candidate antecedent; the next row records that the antecedent is false.",
        ),
        ReductionRow(
            "jwpfaoehr_12_rejected_endpoint",
            "countermodel_gate",
            "rejected_by_counterexample",
            "The candidate all-order signed endpoint hierarchy fails at its first uncompleted order.",
            exact["candidate_endpoint"] + "; but " + exact["order10_counterexample"],
            "Rejects the all-order signed-Hankel antecedent for the actual endpoint sequence, not the exact conditional heat algebra or RH.",
            {"source": sources[4]},
        ),
        ReductionRow(
            "jwpfaoehr_13_bridge_boundary",
            "route_separation_guard",
            "separate_open_obligation",
            "The Jensen-window PF problem remains separate and now requires an antecedent weaker than all-shift all-order signed-Hankel positivity.",
            exact["bridge_boundary"],
            "Blocks promotion of signed-Hankel sign regularity directly to RH.",
        ),
    ]

    return {
        "kind": "jensen_window_pf_all_order_endpoint_heat_reduction",
        "date": "2026-07-16",
        "status": (
            "exact all-order endpoint-to-heat equivalence with a completed base "
            "through order nine and the candidate endpoint hierarchy rejected "
            "by the order-ten counterexample"
        ),
        "proof_boundary": (
            "This artifact proves the arbitrary-order determinant identities, "
            "single-layer tail propagation, induction over compound order, and "
            "the endpoint/heat-interval equivalence. The all-order positivity "
            "antecedent is false at order ten. This does not prove the "
            "Jensen-window PF bridge, PF-infinity, Jensen hyperbolicity, RH, or "
            "Lambda<=0."
        ),
        "sources": sources,
        "source_contract": source_contract,
        "exact": exact,
        "symbolic_identity_audit": symbolic,
        "arithmetic_audit": arithmetic,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": sum(
                row.readiness == "ready_to_apply" for row in rows
            ),
            "open_endpoint_rows": sum(
                row.readiness == "not_ready_to_apply" for row in rows
            ),
            "rejected_endpoint_rows": sum(
                row.readiness == "rejected_by_counterexample" for row in rows
            ),
            "separate_bridge_rows": sum(
                row.readiness == "separate_open_obligation" for row in rows
            ),
            "symbolic_specialization_orders": len(SYMBOLIC_ORDERS),
            "orientation_parity_checks": arithmetic["parity_checks"],
            "all_fixed_order_tail_theorems": 1,
            "all_order_affine_identities": 1,
            "all_order_flag_plucker_identities": 1,
            "all_order_cooperative_recursions": 1,
            "single_layer_propagation_lemmas": 1,
            "endpoint_interval_equivalences": 1,
            "arbitrary_column_consequences": 1,
            "completed_base_order": 9,
            "first_failed_order": 10,
            "negative_order10_endpoint_rows": 4,
            "open_endpoint_hierarchies": 0,
            "rejected_endpoint_hierarchies": 1,
            "separate_jensen_pf_bridges": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_all_order_endpoint_heat_reduction.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_all_order_endpoint_heat_reduction.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF All-Order Endpoint-To-Heat Reduction",
        "",
        "Date: 2026-07-16",
        "",
        "Status: exact all-order endpoint-to-heat equivalence, with the signed",
        "base completed through order nine and the candidate continuation",
        "rejected by a rigorous order-ten endpoint counterexample. This is not",
        "a proof of the separate Jensen-window PF bridge, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_all_order_endpoint_heat_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_all_order_endpoint_heat_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_all_order_endpoint_heat_reduction.py",
        "```",
        "",
        "## Coordinates And Eventual Tail",
        "",
        "Put",
        "",
        "```text",
        exact["definitions"],
        "```",
        "",
        "The graded-kernel Vandermonde theorem has the correctly ordered",
        "quantifiers",
        "",
        "```text",
        exact["fixed_order_tail"],
        "```",
        "",
        "The threshold may depend on `m`. That is sufficient below because the",
        "argument fixes one finite `m` before choosing and using `N_m`.",
        "",
        "## Arbitrary-Order Heat Algebra",
        "",
        "Let `C_j=(A_(n+i+j))_(0<=i<m)`. Under the sequence-shift derivation,",
        "multilinearity gives",
        "",
        "```text",
        exact["delta_border"],
        "```",
        "",
        "Every earlier shifted column duplicates its right neighbor. Splitting",
        "`4(n+i+j)+2` into a column part and a row part leaves only the final",
        "shifted column and final shifted row:",
        "",
        "```text",
        exact["affine_cancellation"],
        exact["affine_heat"],
        "```",
        "",
        "For the `m` by `(m+1)` matrix with columns `C_0,...,C_m`, apply the",
        "three-term flag-Plucker relation to its maximal minors and the maximal",
        "minors of its first `m-1` rows. With the Hankel index map this is",
        "",
        "```text",
        exact["flag_plucker"],
        "```",
        "",
        "Since `c_(m-1,n+1)=c_(m,n)-4`, division by the positive preceding",
        "layer yields",
        "",
        "```text",
        exact["cooperative_flow"],
        "```",
        "",
        "Thus the order-`m` system is one-sided cooperative whenever order",
        "`m-1` has already been completed.",
        "",
        "## Tail-To-Endpoint Propagation",
        "",
        "Fix `m`. Choose its own `N_m` from the eventual-tail theorem. For",
        "`n=N_m-1,...,0`, variation of constants reads",
        "",
        "```text",
        exact["variation_of_constants"],
        "```",
        "",
        "The endpoint term is strictly positive. Backward induction in `n`",
        "makes the integrand nonnegative because `a_(m,n)>0` and the next shift",
        "is positive throughout the interval. Hence",
        "",
        "```text",
        exact["single_layer_implication"],
        "```",
        "",
        "There is no exchange of `forall m` and `exists N_m` in this proof.",
        "Ordinary induction completes one order before moving to the next.",
        "",
        "## Exact All-Order Reduction",
        "",
        "The fixed-order programme now supplies",
        "",
        "```text",
        exact["known_base"],
        "```",
        "",
        "Iterating the single-layer lemma from `m=10` proves the exact",
        "equivalence",
        "",
        "```text",
        exact["endpoint_interval_equivalence"],
        "```",
        "",
        "The reverse implication is simply restriction to `lambda=-100`. Once",
        "the contiguous hierarchy is complete, the fixed-order initial-minor",
        "criterion applies separately at each finite `m` and each `lambda`,",
        "giving",
        "",
        "```text",
        exact["arbitrary_column_consequence"],
        "```",
        "",
        "## Rejected Endpoint Antecedent",
        "",
        "Desnanot-Jacobi and the orientation identity",
        "`binom(m,2)+binom(m-2,2)-2binom(m-1,2)=1` give",
        "",
        "```text",
        exact["signed_condensation"],
        exact["endpoint_coordinate"],
        "```",
        "",
        "The conditional all-order route would require the static statement",
        "",
        "```text",
        exact["candidate_endpoint"],
        "```",
        "",
        "equivalently the continuation from order nine of the strict",
        "log-concavity hierarchy in the positive signed coordinates. The",
        "order-ten counterexample now proves that this antecedent is false:",
        "",
        "```text",
        exact["order10_counterexample"],
        "```",
        "",
        "Thus the endpoint/interval equivalence remains exact, but it cannot",
        "serve as a proof route for the actual sequence.",
        "",
        "## Separate Bridge",
        "",
        "Abstract all-order shifted signed-Hankel positivity would not",
        "automatically establish Jensen hyperbolicity:",
        "",
        "```text",
        exact["bridge_boundary"],
        "```",
        "",
        "The binomially weighted Jensen-window problem remains open, but its",
        "replacement Xi/Phi-specific theorem cannot assume the rejected",
        "all-shift signed-Hankel hierarchy.",
        "",
        "## Machine Audit",
        "",
        f"The artifact records `{summary['orientation_parity_checks']}` exact orientation-parity checks",
        f"and independent symbolic determinant specializations at `{summary['symbolic_specialization_orders']}`",
        "orders. Those finite checks audit the implementation; arbitrary-order",
        "validity comes from the determinant cancellation and the two general",
        "determinant identities proved above.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "wrote all-order endpoint-to-heat reduction: "
        f"{artifact['summary']['rows']} rows, "
        f"base through order {artifact['summary']['completed_base_order']}, "
        f"first failed order {artifact['summary']['first_failed_order']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
