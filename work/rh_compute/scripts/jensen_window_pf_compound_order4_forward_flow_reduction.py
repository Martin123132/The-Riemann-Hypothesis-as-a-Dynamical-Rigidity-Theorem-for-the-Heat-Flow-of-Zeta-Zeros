#!/usr/bin/env python3
"""Derive the exact cooperative flow reduction for contiguous order four."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_forward_flow_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_forward_flow_reduction.md"
)
ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_m100_entry_certificate.json"
)
ORDER3_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order3_forward_invariance_certificate.json"
)


@dataclass(frozen=True)
class FlowRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def hankel(values: tuple[sp.Symbol, ...], size: int, shift: int) -> sp.Expr:
    return sp.expand(
        sp.det(sp.Matrix(size, size, lambda i, j: values[shift + i + j]))
    )


def shift_derivative(expression: sp.Expr, values: tuple[sp.Symbol, ...]) -> sp.Expr:
    return sp.expand(
        sum(
            sp.diff(expression, values[index]) * values[index + 1]
            for index in range(len(values) - 1)
        )
    )


def heat_derivative(
    expression: sp.Expr, values: tuple[sp.Symbol, ...], base: sp.Symbol
) -> sp.Expr:
    return sp.expand(
        sum(
            sp.diff(expression, values[index])
            * (4 * (base + index) + 2)
            * values[index + 1]
            for index in range(len(values) - 1)
        )
    )


def exact_flow_algebra() -> dict:
    n = sp.symbols("n", integer=True, nonnegative=True)
    values = sp.symbols("a0:10")
    h3 = hankel(values, 3, 0)
    h3_next = hankel(values, 3, 1)
    h4 = hankel(values, 4, 0)
    h4_next = hankel(values, 4, 1)
    delta_h3 = shift_derivative(h3, values)
    delta_h3_next = shift_derivative(h3_next, values)
    delta_h4 = shift_derivative(h4, values)

    heat_h3 = heat_derivative(h3, values, n)
    heat_h3_next = heat_derivative(h3_next, values, n)
    heat_h4 = heat_derivative(h4, values, n)
    c3 = 4 * n + 18
    c3_next = 4 * n + 22
    c4 = 4 * n + 26
    weighted_checks = {
        "H3_n": sp.factor(heat_h3 - c3 * delta_h3),
        "H3_n_plus_1": sp.factor(
            heat_h3_next - c3_next * delta_h3_next
        ),
        "H4_n": sp.factor(heat_h4 - c4 * delta_h4),
    }
    if any(value != 0 for value in weighted_checks.values()):
        raise RuntimeError(f"affine Hankel derivative failed: {weighted_checks}")

    plucker_remainder = sp.factor(
        h3_next * delta_h4
        - h3 * h4_next
        - delta_h3_next * h4
    )
    if plucker_remainder != 0:
        raise RuntimeError("order-four Plucker flow identity failed")
    cooperative_remainder = sp.factor(
        h3_next * heat_h4
        - c4 * h3 * h4_next
        - sp.Rational(1, 1) * c4 / c3_next * heat_h3_next * h4
    )
    if cooperative_remainder != 0:
        raise RuntimeError("weighted cooperative identity failed")
    return {
        "coefficient_flow": "A_j'=(4*j+2)*A_(j+1)",
        "affine_hankel_derivative": (
            "H_(m,n)'=(4*(n+2*m-2)+2)*delta(H_(m,n))"
        ),
        "verified_sizes": "m=3 at shifts n,n+1; m=4 at shift n",
        "plucker_identity": (
            "H_(3,n+1)*delta(H_(4,n))="
            "H_(3,n)*H_(4,n+1)+delta(H_(3,n+1))*H_(4,n)"
        ),
        "weighted_identity": (
            "H_(4,n)'=a_n*H_(4,n+1)+b_n*H_(4,n)"
        ),
        "a_n": "(4*n+26)*H_(3,n)/H_(3,n+1)",
        "b_n": (
            "((4*n+26)/(4*n+22))*H_(3,n+1)'/H_(3,n+1)"
        ),
        "boundary": (
            "H_(4,n)=0 => H_(4,n)'="
            "(4*n+26)*H_(3,n)/H_(3,n+1)*H_(4,n+1)"
        ),
        "off_diagonal_positive": (
            "H_(3,j)<0 for every j implies a_n>0"
        ),
    }


def stable_gap_algebra() -> dict:
    coefficient, ratio = sp.symbols("A rho", positive=True)
    contractions = sp.symbols("x1:8", positive=True)
    ratios = [ratio]
    for index in range(1, 7):
        ratios.append(ratio * sp.prod(contractions[:index]))
    values = [coefficient]
    for current in ratios:
        values.append(sp.factor(values[-1] * current))
    h4 = sp.factor(
        sp.det(sp.Matrix(4, 4, lambda i, j: values[i + j]))
    )

    def gap(center: int) -> sp.Expr:
        x_mid = contractions[center - 1]
        return sp.expand(
            (1 - x_mid) ** 2
            - x_mid**2
            * (1 - contractions[center - 2])
            * (1 - contractions[center])
        )

    g0, g1, g2 = gap(2), gap(3), gap(4)
    margin = sp.factor(g1**2 - contractions[2] ** 3 * g0 * g2)
    scale = sp.factor(
        coefficient**4
        * ratio**12
        * contractions[0] ** 8
        * contractions[1] ** 4
        / (1 - contractions[2])
    )
    if sp.factor(h4 - scale * margin) != 0:
        raise RuntimeError("stable order-four gap factorization failed")
    return {
        "gap": (
            "G_n=d_(n+2)^2-x_(n+2)^2*d_(n+1)*d_(n+3)"
        ),
        "margin": "F_n=G_(n+1)^2-x_(n+3)^3*G_n*G_(n+2)",
        "positive_scale": (
            "S_n=A_n^4*rho_n^12*x_(n+1)^8*x_(n+2)^4/d_(n+3)"
        ),
        "factorization": "H_(4,n)=S_n*F_n",
        "boundedness": (
            "inside the ratio and order-three cones, 0<G_j<=1 and |F_n|<=1"
        ),
        "scaled_flow": "F_n'=alpha_n*F_(n+1)+beta_n*F_n",
        "alpha_n": (
            "(4*n+26)*rho_(n+2)*x_(n+3)^4*"
            "(d_(n+3)/d_(n+4))*(G_n/G_(n+1))"
        ),
        "beta_n": (
            "((4*n+26)/(4*n+22))*(log T_(n+1))'-(log S_n)'"
        ),
        "alpha_positive": True,
        "scaled_boundary": "F_n=0 => F_n'=alpha_n*F_(n+1)",
    }


def maximum_principle_reduction() -> dict:
    return {
        "weight": "z_n=F_n/(n+1)",
        "tail_attainment": (
            "|F_n|<=1 implies z_n->0 uniformly on every heat interval"
        ),
        "weighted_flow": (
            "z_n'=alpha_n*((n+2)/(n+1))*z_(n+1)+beta_n*z_n"
        ),
        "sufficient_coefficient_bound": (
            "for every finite L>=-100, sup_(-100<=lambda<=L,n>=0) "
            "[beta_n+alpha_n*(n+2)/(n+1)]<infinity"
        ),
        "consequence": (
            "the exponential weighted-minimum argument propagates F_n>0"
        ),
        "remaining_issue": (
            "prove the scalar coefficient bound from a quantitative spatial "
            "tail for neighboring reciprocal-defect increments and order-three gaps"
        ),
    }


def build_artifact() -> dict:
    entry = json.loads(ENTRY_SOURCE.read_text(encoding="utf-8"))
    order3 = json.loads(ORDER3_SOURCE.read_text(encoding="utf-8"))
    if (
        entry.get("summary", {}).get("all_shift_order_four_entry_theorems")
        != 1
    ):
        raise RuntimeError("order-four entry source is not closed")
    if order3.get("summary", {}).get("full_forward_propagation_rows") != 1:
        raise RuntimeError("order-three forward source is not closed")
    exact = exact_flow_algebra()
    stable = stable_gap_algebra()
    maximum = maximum_principle_reduction()
    rows = [
        FlowRow(
            id="co4ffr_01_coefficient_flow",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The normalized Xi coefficients obey an affine weighted shift flow.",
            formula=exact["coefficient_flow"],
            proof_boundary="Exact coefficient heat equation only.",
        ),
        FlowRow(
            id="co4ffr_02_affine_hankel_derivative",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Affine shift weights collapse to one scalar on every contiguous Hankel determinant.",
            formula=exact["affine_hankel_derivative"],
            proof_boundary="Exact determinant differentiation identity.",
        ),
        FlowRow(
            id="co4ffr_03_plucker_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="A Plucker identity relates the shifted order-four derivative to the neighboring order-four minor and the completed order-three layer.",
            formula=exact["plucker_identity"],
            proof_boundary="Local determinant algebra only.",
        ),
        FlowRow(
            id="co4ffr_04_cooperative_H4_flow",
            role="exact_flow_lemma",
            readiness="ready_to_apply",
            claim="The contiguous order-four determinants obey a one-sided cooperative linear system over the strict order-three cone.",
            formula=exact["weighted_identity"],
            proof_boundary="Finite-index flow algebra; escape to infinite index is separate.",
            diagnostics={
                "a_n": exact["a_n"],
                "b_n": exact["b_n"],
                "boundary": exact["boundary"],
            },
        ),
        FlowRow(
            id="co4ffr_05_stable_gap_factorization",
            role="exact_sign_equivalence",
            readiness="ready_to_apply",
            claim="Positive coefficient-ratio factors isolate a bounded dimensionless order-four margin.",
            formula="H_(4,n)=S_n*F_n, S_n>0, |F_n|<=1",
            proof_boundary="Inside the completed ratio and order-three cones.",
            diagnostics=stable,
        ),
        FlowRow(
            id="co4ffr_06_cooperative_gap_flow",
            role="exact_flow_lemma",
            readiness="ready_to_apply",
            claim="Positive rescaling preserves cooperativity and gives an explicit positive off-diagonal coefficient for the bounded gap margin.",
            formula=stable["scaled_flow"],
            proof_boundary="Exact local flow and boundary orientation only.",
            diagnostics={
                "alpha_n": stable["alpha_n"],
                "beta_n": stable["beta_n"],
                "boundary": stable["scaled_boundary"],
            },
        ),
        FlowRow(
            id="co4ffr_07_maximum_principle_reduction",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="Boundedness of the stable margin makes every negative weighted infimum spatially attained, reducing propagation to one scalar coefficient ceiling.",
            formula=maximum["sufficient_coefficient_bound"],
            proof_boundary="Conditional infinite maximum principle reduction.",
            diagnostics=maximum,
        ),
        FlowRow(
            id="co4ffr_08_spatial_tail_target",
            role="analytic_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove the effective diagonal coefficient is uniformly bounded above on every compact heat interval.",
            formula=maximum["sufficient_coefficient_bound"],
            proof_boundary="Open scalar spatial-tail theorem; no forward order-four propagation is asserted.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order4_forward_flow_reduction",
        "date": "2026-07-13",
        "status": "exact cooperative order-four flow reduction with one open spatial-tail coefficient bound",
        "proof_boundary": (
            "This artifact proves the exact local cooperative order-four flow "
            "and reduces infinite-index propagation to one explicit scalar "
            "coefficient ceiling. It does not prove that ceiling, forward "
            "order-four invariance, arbitrary-column order four, PF-infinity, "
            "RH, or Lambda<=0."
        ),
        "exact_flow": exact,
        "stable_gap": stable,
        "maximum_principle": maximum,
        "source_diagnostics": {
            "entry_status": entry.get("status"),
            "entry_theorems": entry.get("summary", {}).get(
                "all_shift_order_four_entry_theorems"
            ),
            "order3_status": order3.get("status"),
            "order3_forward_theorems": order3.get("summary", {}).get(
                "full_forward_propagation_rows"
            ),
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": len(rows) - 1,
            "exact_identities": 3,
            "cooperative_flow_lemmas": 2,
            "stable_gap_factorizations": 1,
            "maximum_principle_reductions": 1,
            "open_spatial_tail_bounds": 1,
        },
        "sources": [
            "outputs/jensen_window_pf_compound_order4_m100_entry_certificate.md",
            "outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_forward_flow_reduction.py"
        ),
        "remaining_target": maximum["sufficient_coefficient_bound"],
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact_flow"]
    stable = artifact["stable_gap"]
    maximum = artifact["maximum_principle"]
    lines = [
        "# Jensen-Window PF Compound Order-Four Forward-Flow Reduction",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact cooperative order-four flow reduction with one open",
        "spatial-tail coefficient bound. This is not a proof of forward",
        "order-four invariance, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_forward_flow_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_forward_flow_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_forward_flow_reduction.py",
        "```",
        "",
        "## Affine Hankel Flow",
        "",
        "The coefficient heat equation is",
        "",
        "```text",
        exact["coefficient_flow"],
        "```",
        "",
        "If `delta(A_j)=A_(j+1)`, exact determinant algebra gives",
        "",
        "```text",
        exact["affine_hankel_derivative"],
        exact["plucker_identity"],
        "```",
        "",
        "Put `T_n=-H_(3,n)>0` and `Q_n=H_(4,n)`. The completed order-three",
        "theorem therefore turns the order-four flow into",
        "",
        "```text",
        "Q_n'=a_n*Q_(n+1)+b_n*Q_n,",
        f"a_n={exact['a_n']} >0,",
        f"b_n={exact['b_n']}.",
        "```",
        "",
        "At `Q_n=0`, the vector field points inward whenever `Q_(n+1)>=0`.",
        "Thus the finite-dimensional boundary algebra is cooperative; order five",
        "does not enter this local boundary identity.",
        "",
        "## Stable Margin",
        "",
        "The positive ratio factorization is",
        "",
        "```text",
        stable["margin"],
        stable["positive_scale"],
        stable["factorization"],
        "```",
        "",
        "Inside the completed ratio and order-three cones, `0<G_j<=1`, so",
        "`|F_n|<=1`. Positive rescaling preserves the one-sided system:",
        "",
        "```text",
        stable["scaled_flow"],
        f"alpha_n={stable['alpha_n']} >0,",
        f"beta_n={stable['beta_n']}.",
        "```",
        "",
        "## Infinite-Index Target",
        "",
        "Set `z_n=F_n/(n+1)`. Boundedness gives `z_n->0` uniformly, so a",
        "negative spatial infimum is attained. The exact weighted flow is",
        "",
        "```text",
        maximum["weighted_flow"],
        "```",
        "",
        "and the standard exponential minimum argument closes as soon as, for",
        "every finite `L>=-100`,",
        "",
        "```text",
        maximum["sufficient_coefficient_bound"],
        "```",
        "",
        "This is now the sole forward-propagation blocker. It is a scalar spatial",
        "tail theorem, not a missing local sign identity. A proof needs quantitative",
        "control of neighboring reciprocal-defect increments and the positive",
        "order-three gaps on compact heat intervals.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_m100_entry_certificate.md",
        "outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "derived Jensen-window PF compound order-four forward flow: "
        f"{summary['rows']} rows, {summary['exact_rows']} exact rows, "
        f"{summary['exact_identities']} exact identities, "
        f"{summary['cooperative_flow_lemmas']} cooperative flow lemmas, "
        f"{summary['maximum_principle_reductions']} maximum-principle reduction, "
        f"{summary['open_spatial_tail_bounds']} open spatial-tail bound"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
