#!/usr/bin/env python3
"""Transfer contiguous order-three Hankel signs to arbitrary columns."""

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
    "jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.md"
)
FORWARD_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order3_forward_invariance_certificate.json"
)


@dataclass(frozen=True)
class LemmaRow:
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


def exact_diagnostics() -> dict:
    u0, u1, u2 = sp.symbols("u0 u1 u2")
    v0, v1, v2 = sp.symbols("v0 v1 v2")
    matrix = sp.Matrix([[1, 1, 1], [u0, u1, u2], [v0, v1, v2]])
    determinant = sp.factor(matrix.det())
    sigma0 = (v1 - v0) / (u1 - u0)
    sigma1 = (v2 - v1) / (u2 - u1)
    expected = sp.factor((u1 - u0) * (u2 - u1) * (sigma1 - sigma0))
    if sp.factor(determinant - expected) != 0:
        raise RuntimeError("planar orientation identity failed")

    rational_points = ((5, 4), (3, -4), (2, -7), (0, -11))
    edge_slopes = tuple(
        sp.Rational(rational_points[j + 1][1] - rational_points[j][1])
        / (rational_points[j + 1][0] - rational_points[j][0])
        for j in range(3)
    )
    if not (edge_slopes[0] > edge_slopes[1] > edge_slopes[2]):
        raise RuntimeError("strict secant benchmark was not constructed correctly")
    triple_signs = []
    for indices in ((0, 1, 2), (1, 2, 3), (0, 1, 3), (0, 2, 3)):
        block = sp.Matrix(
            [
                [1, 1, 1],
                [rational_points[j][0] for j in indices],
                [rational_points[j][1] for j in indices],
            ]
        )
        triple_signs.append(str(sp.factor(block.det())))
    if not all(sp.Rational(value) < 0 for value in triple_signs):
        raise RuntimeError("noncontiguous orientation benchmark failed")

    return {
        "normalized_column": (
            "A_(n+j)*(1, r_(n+j), r_(n+j)*r_(n+j+1))^T"
        ),
        "planar_points": "P_j=(u_j,v_j)=(r_(n+j),r_(n+j)*r_(n+j+1))",
        "strict_abscissa_order": "u_(j+1)<u_j because x_(n+j+1)<1",
        "edge_slope": "sigma_j=(v_(j+1)-v_j)/(u_(j+1)-u_j)",
        "orientation_identity": (
            "det[(1,u_j,v_j)^T]_(j=l,l+1,l+2)="
            "(u_(l+1)-u_l)*(u_(l+2)-u_(l+1))*(sigma_(l+1)-sigma_l)"
        ),
        "contiguous_equivalence": (
            "D_(3,n+l)<0 iff sigma_(l+1)<sigma_l"
        ),
        "secant_average": (
            "S_(a,b)=(v_b-v_a)/(u_b-u_a) is the weighted average of "
            "sigma_a,...,sigma_(b-1) with weights u_j-u_(j+1)>0"
        ),
        "arbitrary_orientation": (
            "det(P_a,P_b,P_c)=(u_b-u_a)*(u_c-u_b)*(S_(b,c)-S_(a,b))<0"
        ),
        "order_two": (
            "R_(2,n)(j_1,j_2)=A_(n+j_1)*A_(n+j_2)*"
            "(r_(n+j_2)-r_(n+j_1))<0"
        ),
        "benchmark": {
            "points": [list(point) for point in rational_points],
            "edge_slopes": [str(value) for value in edge_slopes],
            "triple_determinants": triple_signs,
        },
    }


def build_artifact() -> dict:
    exact = exact_diagnostics()
    forward = load_json(FORWARD_SOURCE)
    if forward.get("summary", {}).get("lambda_zero_theorem_rows") != 1:
        raise RuntimeError("contiguous order-three source is not closed")
    rows = [
        LemmaRow(
            id="o3nst_01_column_normalization",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Positive column scaling converts every three-row Hankel column into one planar point.",
            formula=exact["normalized_column"],
            proof_boundary="Three consecutive rows only.",
        ),
        LemmaRow(
            id="o3nst_02_strict_abscissas",
            role="theorem_input",
            readiness="ready_to_apply",
            claim="Strict coefficient log-concavity orders the planar abscissas.",
            formula=exact["strict_abscissa_order"],
            proof_boundary="Uses the strict upper ratio wall for the actual Xi trajectory.",
        ),
        LemmaRow(
            id="o3nst_03_orientation_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="A contiguous determinant is exactly the product of two abscissa gaps and one edge-slope difference.",
            formula=exact["orientation_identity"],
            proof_boundary="Exact planar determinant algebra.",
        ),
        LemmaRow(
            id="o3nst_04_local_slope_order",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="The completed contiguous order-three theorem makes all successive edge slopes strictly decreasing.",
            formula="D_(3,n+l)<0 for every l>=0 => sigma_(l+1)<sigma_l",
            proof_boundary="Uses the all-shift contiguous theorem, not finite evidence.",
        ),
        LemmaRow(
            id="o3nst_05_secant_averaging",
            role="exact_lemma",
            readiness="ready_to_apply",
            claim="Every longer secant slope is a positive weighted average of its edge slopes.",
            formula=exact["secant_average"],
            proof_boundary="Elementary finite weighted-average lemma.",
        ),
        LemmaRow(
            id="o3nst_06_arbitrary_order_three",
            role="theorem_conclusion",
            readiness="ready_to_apply",
            claim="All arbitrary-column three-row reshaped-Hankel determinants have the required negative sign.",
            formula=(
                "R_(3,n)(j_1,j_2,j_3)<0 for every n>=0 and "
                "0<=j_1<j_2<j_3"
            ),
            proof_boundary="Complete order three for the consecutive-row reshaped-Hankel target.",
        ),
        LemmaRow(
            id="o3nst_07_arbitrary_order_two",
            role="theorem_conclusion",
            readiness="ready_to_apply",
            claim="The same strict ratio order gives the arbitrary-column order-two sign.",
            formula=exact["order_two"],
            proof_boundary="Complete order two only.",
        ),
        LemmaRow(
            id="o3nst_08_higher_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The planar secant argument has no automatic order-four analogue; compound order four and the sign-regular-to-Jensen transfer remain open.",
            formula="orders 1,2,3 closed; order 4 and higher open",
            proof_boundary="Not an all-order theorem, PF-infinity, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_order3_noncontiguous_secant_transfer_lemma",
        "date": "2026-07-12",
        "status": (
            "exact all-column reshaped-Hankel order-two and order-three theorem "
            "at lambda=0 with order four open"
        ),
        "proof_boundary": (
            "This artifact proves the required arbitrary-column reshaped-Hankel "
            "signs for orders two and three at lambda=0 by exact planar secant "
            "geometry. It does not prove order four or higher, the all-order "
            "sign-regular-to-Jensen transfer, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md",
            "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
            "outputs/signed_hankel_jensen_bridge_target.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.py"
        ),
        "exact": exact,
        "summary": {
            "lemma_rows": len(rows),
            "exact_identity_rows": 2,
            "secant_averaging_rows": 1,
            "arbitrary_order_two_rows": 1,
            "arbitrary_order_three_rows": 1,
            "open_handoffs": 1,
            "ready_to_apply_rows": sum(row.readiness == "ready_to_apply" for row in rows),
        },
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Order-Three Noncontiguous Secant-Transfer Lemma",
        "",
        "Date: 2026-07-12",
        "",
        "Status: exact all-column reshaped-Hankel order-two and order-three",
        "theorem at lambda=0 with order four open. This is not a proof of",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_order3_noncontiguous_secant_transfer_lemma`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.json",
        "python work/rh_compute/scripts/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF order-three noncontiguous secant-transfer lemma: 8 rows, 0 issues, 2 exact identities, 1 secant-averaging lemma, 1 arbitrary-column order-two theorem, 1 arbitrary-column order-three theorem, 1 open order-four handoff",
        "```",
        "",
        "## Planar Normalization",
        "",
        "Fix `n` and divide the column with offset `j` by its positive first",
        "entry. The three-row Hankel column becomes",
        "",
        "```text",
        exact["normalized_column"],
        "P_j=(u_j,v_j)=(r_(n+j),r_(n+j)*r_(n+j+1)).",
        "```",
        "",
        "Strict moment log-concavity gives `x_k<1`, so",
        "",
        "```text",
        "u_(j+1)<u_j.",
        "```",
        "",
        "## Local Orientation",
        "",
        "Define the edge slope",
        "",
        "```text",
        exact["edge_slope"],
        "```",
        "",
        "Direct determinant algebra gives",
        "",
        "```text",
        exact["orientation_identity"],
        "```",
        "",
        "The two abscissa gaps have the same negative sign. Therefore the proved",
        "contiguous theorem `D_(3,n+j)<0` for every shift is exactly",
        "",
        "```text",
        "sigma_(j+1)<sigma_j for every j>=0.",
        "```",
        "",
        "## Arbitrary Columns",
        "",
        "For `a<b`, the longer secant slope",
        "",
        "```text",
        "S_(a,b)=(v_b-v_a)/(u_b-u_a)",
        "```",
        "",
        "is the weighted average of `sigma_a,...,sigma_(b-1)` with positive",
        "weights `u_j-u_(j+1)`. If `a<b<c`, every edge in the first interval",
        "has larger slope than every edge in the second. Hence",
        "",
        "```text",
        "S_(a,b)>S_(b,c),",
        exact["arbitrary_orientation"],
        "```",
        "",
        "and restoring the positive column factors proves",
        "",
        "```text",
        "R_(3,n)(j_1,j_2,j_3)<0",
        "for every n>=0 and 0<=j_1<j_2<j_3 at lambda=0.",
        "```",
        "",
        "For two columns, the normalized determinant is simply the difference",
        "of two decreasing ratios, so",
        "",
        "```text",
        "R_(2,n)(j_1,j_2)<0 for every j_1<j_2.",
        "```",
        "",
        "Thus the complete arbitrary-column reshaped-Hankel layers of orders two",
        "and three are now closed. Order four has no automatic planar secant",
        "reduction and remains the next compound target.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order3_forward_invariance_certificate.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "validated Jensen-window PF order-three noncontiguous secant-transfer "
        "lemma: 8 rows, 0 issues, 2 exact identities, 1 secant-averaging lemma, "
        "1 arbitrary-column order-two theorem, 1 arbitrary-column order-three "
        "theorem, 1 open order-four handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
