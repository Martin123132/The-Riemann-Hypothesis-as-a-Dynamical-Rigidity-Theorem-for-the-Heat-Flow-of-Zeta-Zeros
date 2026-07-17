#!/usr/bin/env python3
"""Transfer contiguous signed-Hankel signs to arbitrary columns through total positivity."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from itertools import combinations
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_order4_noncontiguous_total_positivity_transfer.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md"
)
ORDER3_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.json"
)
ORDER4_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.json"
)


@dataclass(frozen=True)
class TransferRow:
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


def determinant_reversal_audits() -> list[dict]:
    audits = []
    for order in range(1, 5):
        symbols = sp.symbols(f"z0:{2 * order - 1}")
        hankel = sp.Matrix(
            [[symbols[i + j] for j in range(order)] for i in range(order)]
        )
        reversed_columns = hankel[:, ::-1]
        orientation = (-1) ** (order * (order - 1) // 2)
        residual = sp.expand(reversed_columns.det() - orientation * hankel.det())
        if residual != 0:
            raise RuntimeError(f"column-reversal determinant audit failed at order {order}")
        audits.append(
            {
                "order": order,
                "orientation": orientation,
                "identity": (
                    f"det(reverse_columns(H_{order}))="
                    f"(-1)^binom({order},2)*det(H_{order})"
                ),
                "residual": str(residual),
            }
        )
    return audits


def index_mapping_audits() -> dict:
    audits = 0
    N = 8
    for n in range(3):
        for order in range(1, 5):
            for row_start in range(5 - order):
                for column_start in range(N - order + 2):
                    shift = n + row_start + N - column_start - order + 1
                    if shift < 0:
                        raise RuntimeError("negative mapped Hankel shift")
                    for i in range(order):
                        for j in range(order):
                            block_index = n + row_start + i + N - (column_start + j)
                            reversed_hankel_index = shift + i + (order - 1 - j)
                            if block_index != reversed_hankel_index:
                                raise RuntimeError("solid-minor index map failed")
                    audits += 1
    return {
        "finite_column_bound": N,
        "solid_blocks_checked": audits,
        "formula": (
            "det B[r:r+k,c:c+k]=(-1)^binom(k,2)*H_(k,n+r+N-c-k+1)"
        ),
    }


def reciprocal_factorial_benchmark() -> dict:
    signed_minor_rows = 0
    initial_minor_rows = 0
    minimum_signed_minor: sp.Rational | None = None
    N = 8

    for n in range(4):
        for order in range(1, 5):
            orientation = (-1) ** (order * (order - 1) // 2)
            for columns in combinations(range(N + 1), order):
                determinant = sp.det(
                    sp.Matrix(
                        [
                            [
                                sp.Rational(1, sp.factorial(n + i + j))
                                for j in columns
                            ]
                            for i in range(order)
                        ]
                    )
                )
                signed = orientation * determinant
                if signed <= 0:
                    raise RuntimeError("reciprocal-factorial signed-minor benchmark failed")
                minimum_signed_minor = (
                    signed
                    if minimum_signed_minor is None
                    else min(minimum_signed_minor, signed)
                )
                signed_minor_rows += 1

    n = 0
    for order in range(1, 5):
        for row_start in range(5 - order):
            for column_start in range(N - order + 2):
                if row_start != 0 and column_start != 0:
                    continue
                block = sp.Matrix(
                    [
                        [
                            sp.Rational(
                                1,
                                sp.factorial(
                                    n + row_start + i + N - column_start - j
                                ),
                            )
                            for j in range(order)
                        ]
                        for i in range(order)
                    ]
                )
                if block.det() <= 0:
                    raise RuntimeError("reciprocal-factorial initial minor was not positive")
                initial_minor_rows += 1

    if initial_minor_rows != 4 * (N + 1):
        raise RuntimeError("initial-minor count did not equal rows times columns")

    return {
        "sequence": "a_k=1/k!",
        "shifts": "n=0,1,2,3",
        "column_offsets": f"0<=j<={N}",
        "strict_signed_minors": signed_minor_rows,
        "positive_initial_minors_for_4_by_9_block": initial_minor_rows,
        "minimum_signed_minor": str(minimum_signed_minor),
    }


def nonpromotion_countermodel() -> dict:
    values = tuple(map(sp.Integer, (10, 9, 29, 18, 21, 25, 3, 16)))

    def determinant(n: int, columns: tuple[int, ...]) -> sp.Integer:
        order = len(columns)
        return sp.det(
            sp.Matrix([[values[n + i + j] for j in columns] for i in range(order)])
        )

    contiguous_0 = determinant(0, (0, 1, 2, 3))
    contiguous_1 = determinant(1, (0, 1, 2, 3))
    noncontiguous = determinant(0, (0, 1, 3, 4))
    lower_order = determinant(0, (0, 1))
    if not (contiguous_0 > 0 and contiguous_1 > 0 and noncontiguous < 0):
        raise RuntimeError("order-four-only nonpromotion countermodel failed")
    if lower_order <= 0:
        raise RuntimeError("countermodel did not visibly violate the lower signed layer")
    return {
        "sequence": [int(value) for value in values],
        "H4_n0": str(contiguous_0),
        "H4_n1": str(contiguous_1),
        "R4_columns_0_1_3_4": str(noncontiguous),
        "H2_n0_wrong_sign": str(lower_order),
        "conclusion": (
            "contiguous order four alone does not imply arbitrary-column order four; "
            "the lower initial-minor signs are essential"
        ),
    }


def exact_diagnostics() -> dict:
    return {
        "orientation": "epsilon_k=(-1)^binom(k,2)",
        "reversed_block": "B_(i,q)^(n,N)=A_(n+i+N-q)(0)",
        "solid_minor_identity": (
            "det B[r:r+k,c:c+k]=epsilon_k*H_(k,n+r+N-c-k+1)(0)"
        ),
        "arbitrary_column_identity": (
            "det B[0:k|N-j_k,...,N-j_1]=epsilon_k*"
            "R_(k,n)(j_1,...,j_k)"
        ),
        "determinant_reversal_audits": determinant_reversal_audits(),
        "index_mapping": index_mapping_audits(),
        "reciprocal_factorial_benchmark": reciprocal_factorial_benchmark(),
        "nonpromotion_countermodel": nonpromotion_countermodel(),
    }


def validate_sources() -> dict:
    order3 = load_json(ORDER3_SOURCE)
    order4 = load_json(ORDER4_SOURCE)
    order3_summary = order3.get("summary", {})
    order4_summary = order4.get("summary", {})
    if order3_summary.get("arbitrary_order_two_rows") != 1:
        raise RuntimeError("all-shift order-two source is not closed")
    if order3_summary.get("arbitrary_order_three_rows") != 1:
        raise RuntimeError("all-shift order-three source is not closed")
    if order4_summary.get("lambda_zero_all_shift_theorems") != 1:
        raise RuntimeError("all-shift contiguous order-four source is not closed")
    return {
        "positive_coefficients": "A_s(0)>0 for every s>=0",
        "signed_contiguous_orders": (
            "epsilon_k*H_(k,s)(0)>0 for every s>=0 and k=1,2,3,4"
        ),
        "order3_source_status": order3.get("status"),
        "order4_source_status": order4.get("status"),
    }


def build_artifact() -> dict:
    exact = exact_diagnostics()
    sources = validate_sources()
    rows = [
        TransferRow(
            id="o4ntp_01_signed_contiguous_inputs",
            role="theorem_input",
            readiness="ready_to_apply",
            claim="The actual lambda-zero Xi sequence has the strict signed contiguous signs through order four at every shift.",
            formula=sources["signed_contiguous_orders"],
            proof_boundary="Uses completed all-shift theorems, not finite evidence.",
        ),
        TransferRow(
            id="o4ntp_02_reversed_finite_block",
            role="exact_definition",
            readiness="ready_to_apply",
            claim="Reverse a finite block of Hankel columns to convert the signed orientation into ordinary total positivity.",
            formula=exact["reversed_block"],
            proof_boundary="Finite 4 by (N+1) block for arbitrary N.",
        ),
        TransferRow(
            id="o4ntp_03_solid_minor_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Every solid minor of the reversed block is a signed contiguous Hankel minor.",
            formula=exact["solid_minor_identity"],
            proof_boundary="Exact index map and column-reversal sign.",
            diagnostics=exact["index_mapping"],
        ),
        TransferRow(
            id="o4ntp_04_positive_initial_minors",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="All initial minors of every reversed finite block are strictly positive.",
            formula="all solid minors are positive, hence all initial minors are positive",
            proof_boundary="Requires every lower signed contiguous layer through the block row count.",
        ),
        TransferRow(
            id="o4ntp_05_gasca_pena_criterion",
            role="published_theorem_input",
            readiness="ready_to_apply",
            claim="The rectangular Gasca-Pena criterion promotes positive initial minors to strict total positivity.",
            formula="all initial minors of an m by p real matrix are positive iff the matrix is strictly totally positive",
            proof_boundary="Published finite-dimensional total-positivity criterion.",
        ),
        TransferRow(
            id="o4ntp_06_finite_block_total_positivity",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="Every reversed 4 by (N+1) Xi Hankel block is strictly totally positive.",
            formula="B^(n,N) is strictly totally positive for every n,N>=0",
            proof_boundary="Only four consecutive Hankel rows are used.",
        ),
        TransferRow(
            id="o4ntp_07_arbitrary_column_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="An arbitrary increasing Hankel column set is the reverse of an increasing column set in the finite block.",
            formula=exact["arbitrary_column_identity"],
            proof_boundary="Choose N at least the largest requested column offset.",
        ),
        TransferRow(
            id="o4ntp_08_arbitrary_order_four",
            role="theorem_conclusion",
            readiness="ready_to_apply",
            claim="Every arbitrary-column four-row reshaped-Hankel determinant has the required positive sign at lambda zero.",
            formula=(
                "R_(4,n)(j_1,j_2,j_3,j_4)>0 for every n>=0 and "
                "0<=j_1<j_2<j_3<j_4"
            ),
            proof_boundary="Complete arbitrary-column order four for consecutive rows.",
        ),
        TransferRow(
            id="o4ntp_09_fixed_order_transfer",
            role="exact_structural_theorem",
            readiness="ready_to_apply",
            claim="At every fixed order m, strict signed contiguous signs through m imply all arbitrary-column signed signs through m.",
            formula=(
                "[epsilon_k H_(k,s)>0 for 1<=k<=m, all s] => "
                "[epsilon_k R_(k,n)(j_1,...,j_k)>0 for 1<=k<=m]"
            ),
            proof_boundary="A finite-dimensional theorem for each fixed m; it supplies no new contiguous signs.",
        ),
        TransferRow(
            id="o4ntp_10_lower_layer_kill_gate",
            role="non_promotion_guard",
            readiness="ready_to_apply",
            claim="Contiguous order four alone cannot replace the complete family of lower initial-minor signs.",
            formula="H4_0=288076>0, H4_1=264875>0, but R4_(0,1,3,4)=-231169<0",
            proof_boundary="Exact rational countermodel; it does not concern the actual Xi sequence.",
            diagnostics=exact["nonpromotion_countermodel"],
        ),
        TransferRow(
            id="o4ntp_11_order_five_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The first new signed-Hankel layer is now contiguous order five; once it is proved, arbitrary columns follow automatically from the fixed-order transfer.",
            formula="next target: epsilon_5*H_(5,n)(0)=H_(5,n)(0)>0 for every n>=0",
            proof_boundary="Not an order-five theorem, PF-infinity, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_order4_noncontiguous_total_positivity_transfer",
        "date": "2026-07-13",
        "status": (
            "exact arbitrary-column reshaped-Hankel order-four theorem at lambda zero"
        ),
        "proof_boundary": (
            "This artifact proves arbitrary-column signed-Hankel signs through "
            "order four from the completed contiguous signs and a published "
            "strict-total-positivity criterion. It does not prove contiguous "
            "order five, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.md",
            "outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md",
            "https://doi.org/10.1016/0024-3795(92)90226-Z",
            "https://arxiv.org/abs/1207.3613",
            "outputs/signed_hankel_jensen_bridge_target.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_order4_noncontiguous_total_positivity_transfer.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_order4_noncontiguous_total_positivity_transfer.py"
        ),
        "source_contract": sources,
        "exact": exact,
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": sum(
                row.readiness == "ready_to_apply" for row in rows
            ),
            "determinant_reversal_orders": len(
                exact["determinant_reversal_audits"]
            ),
            "index_mapping_blocks": exact["index_mapping"]["solid_blocks_checked"],
            "reciprocal_factorial_signed_minors": exact[
                "reciprocal_factorial_benchmark"
            ]["strict_signed_minors"],
            "published_theorem_rows": 1,
            "arbitrary_order_four_theorems": 1,
            "fixed_order_transfer_theorems": 1,
            "nonpromotion_guards": 1,
            "open_handoffs": 1,
        },
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    benchmark = exact["reciprocal_factorial_benchmark"]
    countermodel = exact["nonpromotion_countermodel"]
    lines = [
        "# Jensen-Window PF Order-Four Noncontiguous Total-Positivity Transfer",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact arbitrary-column reshaped-Hankel order-four theorem at",
        "lambda zero. This is not a proof of contiguous order five, PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_order4_noncontiguous_total_positivity_transfer`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.json",
        "python work/rh_compute/scripts/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_order4_noncontiguous_total_positivity_transfer.py",
        "```",
        "",
        "## Published Criterion",
        "",
        "Gasca and Pena, `Total positivity and Neville elimination`, Linear",
        "Algebra and its Applications 165 (1992), 25-44, prove the rectangular",
        "initial-minor criterion for strict total positivity:",
        "",
        "```text",
        "an m by p real matrix is strictly totally positive if and only if",
        "all its initial minors are strictly positive.",
        "```",
        "",
        "An initial minor uses consecutive rows and consecutive columns and",
        "touches the first row or first column. A primary-source restatement is",
        "given in Launois and Lenagan, arXiv:1207.3613, lines 15-17 of the",
        "introduction.",
        "",
        "## Reversed Finite Block",
        "",
        "Fix `n>=0` and `N>=0`, and define the finite reversed Hankel block",
        "",
        "```text",
        exact["reversed_block"],
        "0<=i<=3, 0<=q<=N.",
        "```",
        "",
        "For a solid minor of order `1<=k<=4`, direct index substitution and",
        "column reversal give",
        "",
        "```text",
        exact["solid_minor_identity"],
        "```",
        "",
        "where the mapped shift is nonnegative. The completed lambda-zero",
        "theorems give",
        "",
        "```text",
        "A_s(0)>0, -H_(2,s)(0)>0, -H_(3,s)(0)>0, H_(4,s)(0)>0",
        "for every integer s>=0.",
        "```",
        "",
        "Thus every solid minor, hence every initial minor, of `B^(n,N)` is",
        "positive. The published criterion makes `B^(n,N)` strictly totally",
        "positive.",
        "",
        "## Arbitrary Columns",
        "",
        "Given `0<=j_1<...<j_k<=N`, use the increasing columns",
        "`N-j_k<...<N-j_1` of `B`. Exact reversal gives",
        "",
        "```text",
        exact["arbitrary_column_identity"],
        "```",
        "",
        "and strict total positivity makes the left side positive. Since `N`",
        "can be the largest requested offset,",
        "",
        "```text",
        "R_(4,n)(j_1,j_2,j_3,j_4)>0",
        "for every n>=0 and 0<=j_1<j_2<j_3<j_4 at lambda=0.",
        "```",
        "",
        "More generally, for every fixed `m`,",
        "",
        "```text",
        "epsilon_k=(-1)^binom(k,2)",
        "[epsilon_k H_(k,s)>0 for 1<=k<=m and every s]",
        "  => [epsilon_k R_(k,n)(j_1,...,j_k)>0 for 1<=k<=m].",
        "```",
        "",
        "Therefore arbitrary columns require no new analytic theorem once the",
        "contiguous layers are complete. The first new layer is contiguous",
        "order five.",
        "",
        "## Exact Audits",
        "",
        f"The checker verifies the reversal identity through order 4 and {exact['index_mapping']['solid_blocks_checked']} solid-block index maps.",
        f"The rational benchmark `{benchmark['sequence']}` has {benchmark['strict_signed_minors']} strictly signed arbitrary minors across four shifts and orders one through four.",
        "",
        "The lower layers are essential. The exact positive sequence",
        "",
        "```text",
        str(countermodel["sequence"]),
        "```",
        "",
        "has",
        "",
        "```text",
        f"H_(4,0)={countermodel['H4_n0']}>0, H_(4,1)={countermodel['H4_n1']}>0,",
        f"R_(4,0)(0,1,3,4)={countermodel['R4_columns_0_1_3_4']}<0,",
        f"H_(2,0)={countermodel['H2_n0_wrong_sign']}>0 (wrong signed layer).",
        "```",
        "",
        "So contiguous order four alone cannot be promoted; the proof genuinely",
        "uses every initial-minor sign through order four.",
        "",
        "```text",
        "outputs/jensen_window_pf_order3_noncontiguous_secant_transfer_lemma.md",
        "outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
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
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_note(args.note, artifact)
    print(
        "wrote order-four noncontiguous total-positivity transfer: "
        f"{artifact['summary']['rows']} rows, "
        f"{artifact['summary']['reciprocal_factorial_signed_minors']} benchmark minors"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
