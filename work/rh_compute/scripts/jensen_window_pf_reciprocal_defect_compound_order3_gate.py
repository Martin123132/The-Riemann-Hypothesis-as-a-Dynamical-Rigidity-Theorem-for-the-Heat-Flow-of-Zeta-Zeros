#!/usr/bin/env python3
"""Build the reciprocal-defect coordinate for contiguous order-three minors."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_defect_compound_order3_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_reciprocal_defect_compound_order3_gate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def build_countermodel() -> dict:
    reciprocal_defects = (sp.Integer(10), sp.Rational(109, 10), sp.Rational(58, 5))
    defects = tuple(sp.factor(1 / q**2) for q in reciprocal_defects)
    contractions = tuple(sp.factor(1 - defect) for defect in defects)
    increments = tuple(
        sp.factor(reciprocal_defects[index + 1] - reciprocal_defects[index])
        for index in range(2)
    )
    scaled_defects = tuple(
        sp.factor(sp.Rational(2 * index + 1, 2) * defects[index - 1])
        for index in range(1, 4)
    )

    ratio_rows = [sp.Integer(1)]
    for contraction in contractions:
        ratio_rows.append(sp.factor(ratio_rows[-1] * contraction))
    coefficients = [sp.Integer(1)]
    for ratio in ratio_rows:
        coefficients.append(sp.factor(coefficients[-1] * ratio))
    coefficients = coefficients[:5]

    determinant = sp.factor(
        sp.Matrix(
            [[coefficients[i + j] for j in range(3)] for i in range(3)]
        ).det()
    )
    compound_margin = sp.factor(
        reciprocal_defects[0] * reciprocal_defects[2]
        - reciprocal_defects[1] ** 2
        + 1
    )
    cubic_values = tuple(
        sp.factor(
            contractions[index] ** 2 * contractions[index + 1] ** 2
            - 6 * contractions[index] * contractions[index + 1]
            + 4 * contractions[index]
            + 4 * contractions[index + 1]
            - 3
        )
        for index in range(2)
    )

    if not (
        defects[0] > defects[1] > defects[2] > 0
        and scaled_defects[0] < scaled_defects[1] < scaled_defects[2] < 1
        and all(0 < increment < 1 for increment in increments)
        and all(value < 0 for value in cubic_values)
        and determinant > 0
        and compound_margin < 0
    ):
        raise RuntimeError("strict order-three countermodel failed")
    lower_walls = (sp.Rational(1, 3), sp.Rational(3, 5), sp.Rational(5, 7))
    if not all(x > wall for x, wall in zip(contractions, lower_walls)):
        raise RuntimeError("countermodel lower-wall check failed")
    if not contractions[0] < contractions[1] < contractions[2] < 1:
        raise RuntimeError("countermodel contraction monotonicity failed")

    return {
        "reciprocal_defects": [str(value) for value in reciprocal_defects],
        "defects": [str(value) for value in defects],
        "contractions": [str(value) for value in contractions],
        "increments": [str(value) for value in increments],
        "scaled_defects": [str(value) for value in scaled_defects],
        "coefficients": [str(value) for value in coefficients],
        "cubic_frontier_values": [str(value) for value in cubic_values],
        "compound_margin": str(compound_margin),
        "hankel_determinant": str(determinant),
        "verified_properties": [
            "strict full ratio cone at indices 1,2,3",
            "strictly decreasing defects",
            "strictly increasing scaled defects below one",
            "both reciprocal-defect increments strictly between zero and one",
            "both neighboring cubic Jensen frontier polynomials strictly negative",
            "contiguous order-three signed-Hankel determinant has the wrong positive sign",
        ],
    }


def build_exact() -> dict:
    a0, ratio, x1, x2, x3 = sp.symbols("a0 ratio x1 x2 x3", positive=True)
    ratios = (ratio, ratio * x1, ratio * x1 * x2, ratio * x1 * x2 * x3)
    coefficients = [a0]
    for item in ratios:
        coefficients.append(sp.expand(coefficients[-1] * item))
    determinant = sp.factor(
        sp.Matrix(
            [[coefficients[i + j] for j in range(3)] for i in range(3)]
        ).det()
    )
    frontier = x1 * x2**2 * x3 - x1 * x2**2 - x2**2 * x3 + 2 * x2 - 1
    expected = a0**3 * ratio**6 * x1**3 * frontier
    if sp.factor(determinant - expected) != 0:
        raise RuntimeError("order-three contraction coordinate failed")

    d1, d2, d3 = sp.symbols("d1 d2 d3", positive=True)
    defect_frontier = sp.factor(
        frontier.subs({x1: 1 - d1, x2: 1 - d2, x3: 1 - d3})
    )
    expected_defect = d1 * d3 * (1 - d2) ** 2 - d2**2
    if sp.expand(defect_frontier - expected_defect) != 0:
        raise RuntimeError("order-three defect coordinate failed")

    q1, q2, q3 = sp.symbols("q1 q2 q3", positive=True)
    compound_margin = q1 * q3 - q2**2 + 1
    compound_gap = sp.factor(
        (d2**2 - (1 - d2) ** 2 * d1 * d3).subs(
            {d1: q1**-2, d2: q2**-2, d3: q3**-2}
        )
    )
    positive_factorization = (
        compound_margin
        * (q1 * q3 + q2**2 - 1)
        / (q1**2 * q2**4 * q3**2)
    )
    if sp.factor(compound_gap - positive_factorization) != 0:
        raise RuntimeError("reciprocal-defect sign factorization failed")

    left_increment, right_increment = sp.symbols(
        "left_increment right_increment", nonnegative=True
    )
    increment_form = sp.expand(
        compound_margin.subs(
            {
                q1: q2 - left_increment,
                q3: q2 + right_increment,
            },
            simultaneous=True,
        )
    )
    expected_increment = (
        1
        - left_increment * right_increment
        + q2 * (right_increment - left_increment)
    )
    if sp.expand(increment_form - expected_increment) != 0:
        raise RuntimeError("order-three increment budget failed")

    countermodel = build_countermodel()
    return {
        "determinant_coordinate": {
            "definitions": (
                "rho_k=A_(k+1)/A_k, x_k=rho_k/rho_(k-1), "
                "d_k=1-x_k, q_k=d_k^(-1/2)"
            ),
            "minor": (
                "D_(3,n)=det[A_(n+i+j)]_(i,j=0..2)"
            ),
            "factorization": (
                "D_(3,n)=A_n^3*rho_n^6*x_(n+1)^3*F_n"
            ),
            "frontier": (
                "F_n=x_(n+1)*x_(n+2)^2*x_(n+3)-"
                "x_(n+1)*x_(n+2)^2-x_(n+2)^2*x_(n+3)+"
                "2*x_(n+2)-1"
            ),
            "expected_sign": "The order-three signed-Hankel condition is D_(3,n)<0.",
        },
        "defect_coordinate": {
            "identity": (
                "F_n=d_(n+1)*d_(n+3)*x_(n+2)^2-d_(n+2)^2"
            ),
            "condition": (
                "D_(3,n)<0 iff d_(n+2)^2>"
                "x_(n+2)^2*d_(n+1)*d_(n+3)"
            ),
            "interpretation": (
                "The defect sequence may be log-convex, but its centered "
                "log-curvature must stay below the explicit buffer -2*log x_(n+2)."
            ),
        },
        "reciprocal_defect_coordinate": {
            "condition": (
                "D_(3,n)<0 iff C_n:=q_(n+1)*q_(n+3)-q_(n+2)^2+1>0"
            ),
            "positive_factorization": (
                "d_(n+2)^2-x_(n+2)^2*d_(n+1)*d_(n+3)="
                "C_n*(q_(n+1)*q_(n+3)+q_(n+2)^2-1)/"
                "(q_(n+1)^2*q_(n+2)^4*q_(n+3)^2)"
            ),
            "interpretation": (
                "This is unit-buffered log-concavity of the reciprocal defects."
            ),
        },
        "increment_budget": {
            "definitions": (
                "a_n=q_(n+2)-q_(n+1), b_n=q_(n+3)-q_(n+2)"
            ),
            "identity": (
                "C_n=1-a_n*b_n+q_(n+2)*(b_n-a_n)"
            ),
            "exact_target": (
                "If b_n<a_n, prove q_(n+2)*(a_n-b_n)+a_n*b_n<1."
            ),
        },
        "sufficient_increment_theorem": {
            "hypotheses": "0<=a_n<=b_n<=1",
            "conclusion": "C_n>=1-a_n*b_n>=0",
            "strictness": (
                "C_n>0 unless a_n=b_n=1; the strict cubic cone excludes that equality."
            ),
            "boundary_family": (
                "For q_k=alpha+beta*k with 0<=beta<1, C_n=1-beta^2>0."
            ),
            "scope": (
                "This is sufficient only. Actual reciprocal-defect increments may decrease."
            ),
        },
        "countermodel": countermodel,
        "m100_entry_theorem": {
            "source": (
                "outputs/jensen_window_pf_negative_lambda_m100_"
                "compound_order3_entry_certificate.md"
            ),
            "statement": (
                "A repaired Arb prefix through n=317 and an exact scaled-defect "
                "tail prove C_n(-100)>0 for every n>=0."
            ),
            "uniform_tail_lower": "57613471/66107054971",
            "boundary": (
                "All shifts at lambda=-100, but only the contiguous order-three family."
            ),
        },
        "forward_invariance_theorem": {
            "source": (
                "outputs/jensen_window_pf_compound_order3_"
                "forward_invariance_certificate.md"
            ),
            "statement": (
                "The exact cooperative system C_n'/r_(n+2)=alpha_n*C_(n+1)"
                "+beta_n*C_n, the cubic forward-uniform tail, and a weighted "
                "maximum principle propagate C_n>0 through lambda=0."
            ),
            "conclusion": "C_n(lambda)>0 for every n>=0 and finite lambda>=-100",
            "boundary": (
                "Complete shifted contiguous order three only; noncontiguous "
                "order-three minors are not included."
            ),
        },
        "noncontiguous_transfer_theorem": {
            "source": (
                "outputs/jensen_window_pf_order3_noncontiguous_"
                "secant_transfer_lemma.md"
            ),
            "statement": (
                "Positive column scaling embeds the three-row Hankel columns as "
                "planar points. Contiguous negativity makes their edge slopes "
                "strictly decrease, and secant averaging transfers the sign to "
                "every strictly increasing column triple."
            ),
            "conclusion": (
                "R_(3,n)(j_1,j_2,j_3)<0 for every n>=0 and "
                "0<=j_1<j_2<j_3 at lambda=0"
            ),
            "boundary": "Complete reshaped-Hankel order three, not order four.",
        },
        "live_handoff": {
            "target": (
                "Find coordinates, entry theorems, and forward-invariance laws for "
                "contiguous order four; then transfer to arbitrary columns and "
                "continue to higher compound orders."
            ),
            "forbidden_promotion": (
                "The full ratio cone, decreasing defect, increasing scaled defect, "
                "and strict cubic Jensen increment bounds do not imply C_n>0 "
                "without the lambda=-100 scaled-defect anchor used by the entry theorem."
            ),
        },
    }


def build_artifact() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="rdco3_01_determinant_coordinate",
            role="exact_identity",
            readiness="available_exact",
            claim="The contiguous order-three signed-Hankel minor has an exact contraction-coordinate factorization.",
            formula=exact["determinant_coordinate"]["factorization"],
            proof_boundary="Exact contiguous 3x3 identity only.",
            diagnostics=exact["determinant_coordinate"],
        ),
        GateRow(
            id="rdco3_02_defect_coordinate",
            role="exact_equivalence",
            readiness="ready_to_apply",
            claim="The required minor sign is exactly a buffered defect log-curvature inequality.",
            formula=exact["defect_coordinate"]["condition"],
            proof_boundary="Equivalent to contiguous order three, not arbitrary order-three column sets.",
            diagnostics=exact["defect_coordinate"],
        ),
        GateRow(
            id="rdco3_03_reciprocal_coordinate",
            role="exact_equivalence",
            readiness="ready_to_apply",
            claim="The same sign is unit-buffered log-concavity of reciprocal defects.",
            formula=exact["reciprocal_defect_coordinate"]["condition"],
            proof_boundary="Exact scalar coordinate for one contiguous minor family.",
            diagnostics=exact["reciprocal_defect_coordinate"],
        ),
        GateRow(
            id="rdco3_04_increment_budget",
            role="exact_identity",
            readiness="available_exact",
            claim="The compound margin is an exact first-increment product plus second-increment curvature budget.",
            formula=exact["increment_budget"]["identity"],
            proof_boundary="Exact coordinate; no sign is assumed.",
            diagnostics=exact["increment_budget"],
        ),
        GateRow(
            id="rdco3_05_sufficient_increment",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Nondecreasing reciprocal-defect increments inside the cubic unit cone are sufficient for the order-three compound sign.",
            formula=exact["sufficient_increment_theorem"]["conclusion"],
            proof_boundary="Sufficient theorem only; the monotonicity hypothesis is not proved for Xi.",
            diagnostics=exact["sufficient_increment_theorem"],
        ),
        GateRow(
            id="rdco3_06_boundary_family",
            role="exact_benchmark",
            readiness="available_exact",
            claim="The reciprocal-square rank-two boundary family lies strictly inside the compound cone when its increment is below one.",
            formula=exact["sufficient_increment_theorem"]["boundary_family"],
            proof_boundary="Exact benchmark family, not the Xi trajectory.",
        ),
        GateRow(
            id="rdco3_07_strict_countermodel",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="A positive rational sequence satisfies the proved ratio, adaptive-defect, and strict cubic Jensen cones while reversing the contiguous order-three Hankel sign.",
            formula="q_1=10, q_2=109/10, q_3=58/5, C_0=-181/100",
            proof_boundary="Abstract local countermodel, not the Xi coefficient sequence.",
            diagnostics=exact["countermodel"],
        ),
        GateRow(
            id="rdco3_08_m100_entry",
            role="interval_analytic_theorem",
            readiness="ready_to_apply",
            claim="The actual Xi coefficient sequence has the required contiguous order-three sign at every shift at lambda=-100.",
            formula="C_n(-100)>0 and D_(3,n)(-100)<0 for every n>=0",
            proof_boundary="All-shift entry at one heat parameter; no forward propagation or noncontiguous minors.",
            diagnostics=exact["m100_entry_theorem"],
        ),
        GateRow(
            id="rdco3_09_forward_invariance",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="The contiguous order-three compound margin propagates strictly through lambda=0 at every shift.",
            formula="C_n(lambda)>0 for every n>=0 and finite lambda>=-100",
            proof_boundary="Complete shifted contiguous order three only.",
            diagnostics=exact["forward_invariance_theorem"],
        ),
        GateRow(
            id="rdco3_10_noncontiguous_transfer",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="Planar secant geometry transfers contiguous order-three signs to every increasing column triple.",
            formula=exact["noncontiguous_transfer_theorem"]["conclusion"],
            proof_boundary="Complete arbitrary-column reshaped-Hankel order three only.",
            diagnostics=exact["noncontiguous_transfer_theorem"],
        ),
        GateRow(
            id="rdco3_11_live_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The next compound target is contiguous order four, followed by its arbitrary-column transfer and higher orders.",
            formula=exact["live_handoff"]["target"],
            proof_boundary="Open; not all order-three column sets, the all-order antecedent, RH, or Lambda<=0.",
            diagnostics=exact["live_handoff"],
        ),
    ]
    return {
        "kind": "jensen_window_pf_reciprocal_defect_compound_order3_gate",
        "date": "2026-07-12",
        "status": (
            "complete all-column order-three theorem with strict abstract "
            "countermodel gate"
        ),
        "proof_boundary": (
            "This artifact factors the contiguous 3x3 signed-Hankel minor into "
            "contraction, defect, reciprocal-defect, and increment coordinates; "
            "proves one sufficient increment theorem; and supplies a strict rational "
            "countermodel to promotion from the already proved ratio, adaptive-defect, "
            "and cubic Jensen cones; and composes the all-shift entry, forward-"
            "invariance, and noncontiguous secant-transfer theorems. It does not "
            "prove compound order four or higher, the all-order signed-Hankel "
            "antecedent, the Jensen bridge, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.md",
            "outputs/signed_hankel_jensen_bridge_target.md",
            "outputs/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.md",
            "outputs/jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate.md",
            "outputs/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.md",
            "outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md",
            "outputs/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.md",
            "outputs/formal_core.md",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    countermodel = exact["countermodel"]
    return "\n".join(
        [
            "# Jensen-Window PF Reciprocal-Defect Compound Order-Three Gate",
            "",
            "Date: 2026-07-12",
            "",
            "Status: complete all-column order-three theorem with strict abstract",
            "countermodel gate. This is not a proof of RH or `Lambda <= 0`.",
            "",
            "Artifact kind: `jensen_window_pf_reciprocal_defect_compound_order3_gate`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_reciprocal_defect_compound_order3_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_reciprocal_defect_compound_order3_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_defect_compound_order3_gate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF reciprocal-defect compound order-three gate: 11 rows, 0 issues, 2 exact coordinate identities, 2 exact sign equivalences, 1 sufficient increment theorem, 1 exact boundary benchmark, 1 strict cone countermodel, 1 all-shift lambda=-100 entry theorem, 1 full forward propagation theorem, 1 arbitrary-column order-three theorem, 1 open order-four handoff",
            "```",
            "",
            "## Exact Coordinate",
            "",
            "For the contiguous minor",
            "",
            "```text",
            exact["determinant_coordinate"]["minor"],
            exact["determinant_coordinate"]["factorization"],
            exact["determinant_coordinate"]["frontier"],
            "```",
            "",
            "all factors outside `F_n` are positive. In defect variables this",
            "collapses to",
            "",
            "```text",
            exact["defect_coordinate"]["identity"],
            exact["defect_coordinate"]["condition"],
            "```",
            "",
            "and in reciprocal defects it becomes",
            "",
            "```text",
            exact["reciprocal_defect_coordinate"]["condition"],
            exact["reciprocal_defect_coordinate"]["positive_factorization"],
            "```",
            "",
            "This is the first explicit higher compound concentration coordinate.",
            "",
            "## Increment Budget",
            "",
            "Set consecutive reciprocal-defect increments `a_n,b_n`. Then",
            "",
            "```text",
            exact["increment_budget"]["definitions"],
            exact["increment_budget"]["identity"],
            exact["increment_budget"]["exact_target"],
            "```",
            "",
            "One sufficient theorem is immediate:",
            "",
            "```text",
            exact["sufficient_increment_theorem"]["hypotheses"],
            exact["sufficient_increment_theorem"]["conclusion"],
            exact["sufficient_increment_theorem"]["boundary_family"],
            "```",
            "",
            "Actual Xi increments need not be nondecreasing, so the exact curvature",
            "budget remains the live statement.",
            "",
            "## Strict Countermodel",
            "",
            "Take",
            "",
            "```text",
            f"q_1,q_2,q_3={countermodel['reciprocal_defects']}",
            f"d_1,d_2,d_3={countermodel['defects']}",
            f"x_1,x_2,x_3={countermodel['contractions']}",
            f"s_1,s_2,s_3={countermodel['scaled_defects']}",
            f"q increments={countermodel['increments']}",
            "```",
            "",
            "The resulting positive coefficient prefix is",
            "",
            "```text",
            f"A_0,...,A_4={countermodel['coefficients']}",
            f"two cubic frontier values={countermodel['cubic_frontier_values']}",
            f"C_0={countermodel['compound_margin']}",
            f"det[A_(i+j)]_(i,j=0..2)={countermodel['hankel_determinant']}>0",
            "```",
            "",
            "It lies strictly inside the full ratio cone, has decreasing defects,",
            "increasing scaled defects, reciprocal-defect increments below one, and",
            "strictly hyperbolic neighboring cubic Jensen windows. Nevertheless the",
            "order-three signed-Hankel minor has the wrong sign. This blocks promotion",
            "from every currently proved scalar/cubic cone to the new compound margin.",
            "",
            "## Lambda=-100 Entry",
            "",
            exact["m100_entry_theorem"]["statement"],
            "",
            "```text",
            "C_n(-100)>0 and D_(3,n)(-100)<0 for every n>=0,",
            "uniform analytic tail lower bound=57613471/66107054971.",
            "```",
            "",
            "The finite prefix is spliced to an exact scaled-defect tail theorem;",
            "this is no longer a finite-only observation.",
            "",
            "## Forward Invariance",
            "",
            exact["forward_invariance_theorem"]["statement"],
            "",
            "```text",
            exact["forward_invariance_theorem"]["conclusion"],
            "D_(3,n)(0)<0 for every n>=0.",
            "```",
            "",
            "## Arbitrary Columns",
            "",
            exact["noncontiguous_transfer_theorem"]["statement"],
            "",
            "```text",
            exact["noncontiguous_transfer_theorem"]["conclusion"],
            "```",
            "",
            "## Live Handoff",
            "",
            exact["live_handoff"]["target"],
            "",
            "All-column reshaped-Hankel order three is closed. Compound order four",
            "and every higher order remain open.",
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
    print(f"wrote reciprocal-defect compound order-three gate: {len(artifact['rows'])} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
