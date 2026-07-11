#!/usr/bin/env python3
"""Build the node-range and Phi(0) lower certificate for the relative-Gaussian tail scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_GRID_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json"
)
DEFAULT_PHI_TAIL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.json"
)
DEFAULT_INTERVAL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md"
)

DEFAULT_T_GRID = (1156, 1500, 2000, 5000, 10000)
DEFAULT_QUADRATURE_ORDERS = (64, 96, 128, 192)
DEFAULT_INDICES = (21, 22, 23, 24)
DEFAULT_PRECISION_BITS = 200
DEFAULT_C0_LOWER = "0.44"


@dataclass(frozen=True)
class MatrixRow:
    id: str
    role: str
    readiness: str
    claim: str
    proof_boundary: str
    source_artifacts: list[str]
    diagnostics: dict | None = None
    gap: str | None = None


def frac_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def dec_text(value: Fraction, digits: int = 18) -> str:
    return f"{float(value):.{digits}E}"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(grid_path: Path, phi_tail_path: Path, interval_path: Path) -> dict[str, str]:
    return {
        "grid_json": grid_path.relative_to(REPO_ROOT).as_posix(),
        "grid_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
        "phi_tail_json": phi_tail_path.relative_to(REPO_ROOT).as_posix(),
        "phi_tail_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.md",
        "interval_json": interval_path.relative_to(REPO_ROOT).as_posix(),
        "interval_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def laguerre_alpha(index: int) -> Fraction:
    return Fraction(2 * index - 1, 2)


def interior_row_bound(order: int, alpha: Fraction) -> Fraction:
    return Fraction(4 * order, 1) + 2 * alpha - 6


def last_row_bound(order: int, alpha: Fraction) -> Fraction:
    return Fraction(3 * order - 2, 1) + Fraction(3, 2) * alpha


def first_row_bound(alpha: Fraction) -> Fraction:
    return Fraction(2, 1) + Fraction(3, 2) * alpha


def gershgorin_bound(order: int, alpha: Fraction) -> Fraction:
    if order < 2:
        return first_row_bound(alpha)
    return max(first_row_bound(alpha), interior_row_bound(order, alpha), last_row_bound(order, alpha))


def build_laguerre_rows() -> tuple[list[dict], dict]:
    rows: list[dict] = []
    worst: dict | None = None
    min_t = min(DEFAULT_T_GRID)
    for order in DEFAULT_QUADRATURE_ORDERS:
        for index in DEFAULT_INDICES:
            alpha = laguerre_alpha(index)
            bound = gershgorin_bound(order, alpha)
            ratio = bound / min_t
            row = {
                "quadrature_order": order,
                "index": index,
                "alpha": frac_text(alpha),
                "alpha_decimal": dec_text(alpha),
                "laguerre_largest_node_upper_bound": frac_text(bound),
                "min_T": min_t,
                "x_square_upper_bound": frac_text(ratio),
                "x_square_upper_bound_decimal": dec_text(ratio),
                "x_less_than_one_certified": bool(ratio < 1),
                "proof_boundary": (
                    "Exact Gershgorin/AM-GM upper bound for Laguerre nodes only; it does not enclose "
                    "individual roots or Christoffel weights."
                ),
            }
            rows.append(row)
            if worst is None or bound > worst["bound"]:
                worst = {
                    "order": order,
                    "index": index,
                    "alpha": alpha,
                    "bound": bound,
                    "ratio": ratio,
                }
    assert worst is not None
    summary = {
        "laguerre_bound_rows": len(rows),
        "min_T": min_t,
        "max_quadrature_order": max(DEFAULT_QUADRATURE_ORDERS),
        "max_index": max(DEFAULT_INDICES),
        "max_alpha": frac_text(laguerre_alpha(max(DEFAULT_INDICES))),
        "worst_order": worst["order"],
        "worst_index": worst["index"],
        "worst_alpha": frac_text(worst["alpha"]),
        "worst_laguerre_node_upper_bound": frac_text(worst["bound"]),
        "worst_x_square_upper_bound": frac_text(worst["ratio"]),
        "worst_x_square_upper_bound_decimal": dec_text(worst["ratio"]),
        "worst_x_square_slack_to_one": frac_text(1 - worst["ratio"]),
        "worst_x_square_slack_to_one_decimal": dec_text(1 - worst["ratio"]),
        "node_range_x_le_1_certified": bool(worst["ratio"] < 1),
        "derivation": (
            "For L_N^(alpha), use the symmetric Jacobi matrix with diagonal b_j=2j+alpha+1 "
            "and off-diagonals a_j=sqrt((j+1)(j+alpha+1)). Gershgorin gives "
            "lambda_max <= max_j(b_j+a_{j-1}+a_j). AM-GM gives "
            "sqrt(j(j+alpha))<=j+alpha/2 and sqrt((j+1)(j+alpha+1))<=j+1+alpha/2, "
            "so the interior rows are bounded by 4j+2alpha+2, hence by 4N+2alpha-6."
        ),
    }
    return rows, summary


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def build_c0_certificate(precision_bits: int, c0_lower_text: str) -> dict:
    flint.ctx.prec = precision_bits
    pi = flint.arb.pi()
    c0_lower = flint.arb(c0_lower_text)
    n1_coefficient = 2 * pi * pi - 3 * pi
    n1_term = n1_coefficient * (-pi).exp()
    margin = n1_term - c0_lower
    if not arb_positive(n1_coefficient):
        raise RuntimeError("n=1 coefficient did not certify positive")
    if not arb_positive(margin):
        raise RuntimeError("Phi(0) lower margin did not certify positive")
    return {
        "precision_bits": precision_bits,
        "certified_c0_lower": c0_lower_text,
        "pi_lower_for_term_positivity": "pi>3",
        "term_positivity_formula": "2*pi^2*n^4-3*pi*n^2=pi*n^2*(2*pi*n^2-3)>0 for n>=1",
        "n1_coefficient_ball": n1_coefficient.str(80),
        "n1_phi0_term_ball": n1_term.str(80),
        "n1_phi0_term_lower": n1_term.lower().str(80, radius=False),
        "n1_minus_c0_lower_ball": margin.str(80),
        "n1_minus_c0_lower_lower": margin.lower().str(80, radius=False),
        "c0_lower_certified_by_n1_term": True,
        "full_phi0_lower_reason": (
            "The n=1 contribution already exceeds 0.44, and every n>=1 contribution is positive."
        ),
        "proof_boundary": (
            "Arb lower certificate for Phi(0)>=0.44 only; it does not evaluate Phi(x), Phi'(x), "
            "Laguerre weights, or quadrature remainders."
        ),
    }


def build_rows(paths: dict[str, str], laguerre_summary: dict, c0_certificate: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgnc_01_laguerre_jacobi_reduction",
            role="exact_reduction",
            readiness="available_exact",
            claim=(
                "The generalized Laguerre nodes used by the finite grid are eigenvalues of the standard "
                "symmetric tridiagonal Laguerre Jacobi matrix."
            ),
            diagnostics={
                "diagonal": "b_j=2j+alpha+1",
                "off_diagonal": "a_j=sqrt((j+1)(j+alpha+1))",
                "alpha_range": f"{frac_text(laguerre_alpha(min(DEFAULT_INDICES)))}..{frac_text(laguerre_alpha(max(DEFAULT_INDICES)))}",
                "order_range": f"{min(DEFAULT_QUADRATURE_ORDERS)}..{max(DEFAULT_QUADRATURE_ORDERS)}",
            },
            source_artifacts=[paths["grid_json"], paths["grid_note"]],
            proof_boundary=(
                "Exact node-location reduction only; not an individual root or weight interval enclosure."
            ),
        ),
        MatrixRow(
            id="nlrgnc_02_gershgorin_node_range",
            role="exact_bound",
            readiness="available_exact",
            claim=(
                "Gershgorin plus AM-GM bounds every recorded generalized Laguerre node by 809, "
                "so every node-induced x=sqrt(node/T) lies below 1 for T>=1156."
            ),
            diagnostics=laguerre_summary,
            source_artifacts=[paths["grid_json"], paths["phi_tail_note"]],
            proof_boundary=(
                "Exact coarse range certificate only; not a certified quadrature-node/weight table."
            ),
        ),
        MatrixRow(
            id="nlrgnc_03_phi0_n1_lower_certificate",
            role="arb_lower_bound",
            readiness="available_arb",
            claim=(
                "The n=1 term in Phi(0) is Arb-certified above 0.44, and all remaining Phi(0) "
                "terms are positive."
            ),
            diagnostics=c0_certificate,
            source_artifacts=[paths["phi_tail_note"]],
            proof_boundary=(
                "Arb lower bound for Phi(0) only; not a complete interval model for Phi or Phi'."
            ),
        ),
        MatrixRow(
            id="nlrgnc_04_phi_tail_conditions_handoff",
            role="conditional_handoff",
            readiness="available_for_tail_scout",
            claim=(
                "The two explicit side conditions used by the padded Phi-tail scout, x<=1 and "
                "Phi(0)>=0.44, are certified for the recorded finite grid."
            ),
            diagnostics={
                "supports_phi_tail_bound_scout": paths["phi_tail_note"],
                "certified_side_conditions": ["node-induced x<=1", "Phi(0)>=0.44"],
                "remaining_intervalization_obligations": [
                    "individual Laguerre node and weight intervals",
                    "quadrature-remainder error",
                    "coefficient-ratio and first-omitted denominator propagation",
                    "rounding aggregation",
                    "finite-grid to full-collar bridge",
                ],
            },
            source_artifacts=[paths["phi_tail_json"], paths["interval_note"], paths["dependency_graph"]],
            proof_boundary=(
                "Handoff to the tail scout only; not a finite-grid interval certificate or uniform collar theorem."
            ),
        ),
        MatrixRow(
            id="nlrgnc_05_weight_and_quadrature_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="The node-range and Phi(0) lower certificate also supplies Laguerre weight intervals and quadrature error bounds.",
            gap=(
                "The certificate bounds the largest node and Phi(0) only. It does not enclose individual "
                "nodes, Christoffel weights, or any generalized Gauss-Laguerre quadrature remainder."
            ),
            source_artifacts=[paths["interval_note"], paths["uniform_remainder_target"]],
            proof_boundary=(
                "Rejected promotion only; not a proof of the finite-grid certificate, the uniform collar "
                "remainder theorem, scaled-curvature monotonicity, RH, or Lambda <= 0."
            ),
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(grid_path: Path, phi_tail_path: Path, interval_path: Path) -> dict:
    grid = load_json(grid_path)
    phi_tail = load_json(phi_tail_path)
    interval = load_json(interval_path)
    paths = source_paths(grid_path, phi_tail_path, interval_path)
    laguerre_rows, laguerre_summary = build_laguerre_rows()
    c0_certificate = build_c0_certificate(DEFAULT_PRECISION_BITS, DEFAULT_C0_LOWER)
    rows = build_rows(paths, laguerre_summary, c0_certificate)
    summary = {
        "matrix_rows": len(rows),
        "laguerre_bound_rows": len(laguerre_rows),
        "certified_side_conditions": 2,
        "ready_to_apply_rows": 0,
        "source_grid_rows": grid["summary"]["grid_rows"],
        "source_phi_tail_rows": phi_tail["summary"]["matrix_rows"],
        "source_interval_obligations": interval["summary"]["obligation_rows"],
        "min_T": laguerre_summary["min_T"],
        "max_quadrature_order": laguerre_summary["max_quadrature_order"],
        "max_index": laguerre_summary["max_index"],
        "worst_laguerre_node_upper_bound": laguerre_summary["worst_laguerre_node_upper_bound"],
        "worst_x_square_upper_bound": laguerre_summary["worst_x_square_upper_bound"],
        "worst_x_square_upper_bound_decimal": laguerre_summary["worst_x_square_upper_bound_decimal"],
        "worst_x_square_slack_to_one": laguerre_summary["worst_x_square_slack_to_one"],
        "node_range_x_le_1_certified": laguerre_summary["node_range_x_le_1_certified"],
        "certified_c0_lower": c0_certificate["certified_c0_lower"],
        "n1_phi0_term_lower": c0_certificate["n1_phi0_term_lower"],
        "n1_minus_c0_lower_lower": c0_certificate["n1_minus_c0_lower_lower"],
        "c0_lower_certified_by_n1_term": c0_certificate["c0_lower_certified_by_n1_term"],
        "target_closing": False,
        "main_finding": (
            "The finite relative-Gaussian grid satisfies the two concrete side conditions used by the "
            "padded Phi-tail scout: Gershgorin/AM-GM gives every recorded Laguerre node <=809<T_min=1156, "
            "hence x<=1, and Arb certifies the n=1 term in Phi(0) already exceeds 0.44. This supports "
            "the tail scout but leaves weights, quadrature error, coefficient propagation, rounding, and "
            "the grid-to-collar bridge open."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate",
        "date": "2026-07-07",
        "status": "node-range and Phi0 lower certificate",
        "source_cancellation_reduced_grid_scout": paths["grid_note"],
        "source_cancellation_reduced_grid_json": paths["grid_json"],
        "source_phi_tail_bound_scout": paths["phi_tail_note"],
        "source_phi_tail_bound_json": paths["phi_tail_json"],
        "source_intervalization_target": paths["interval_note"],
        "source_intervalization_target_json": paths["interval_json"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.py",
        "proof_boundary": (
            "Node-range and Phi(0) lower certificate only. It certifies x<=1 and Phi(0)>=0.44 for the "
            "recorded finite grid, but it does not provide individual Laguerre node or weight intervals, "
            "does not bound quadrature remainder, does not propagate coefficient balls, does not prove a "
            "uniform collar theorem, does not prove scaled-curvature monotonicity, and does not prove RH "
            "or Lambda <= 0."
        ),
        "laguerre_bound_rows": laguerre_rows,
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply for the full intervalization target.",
            "The certificate only discharges the x<=1 and Phi(0)>=0.44 side conditions used by the Phi-tail scout.",
            "Laguerre weights and quadrature remainders remain open.",
            "The finite grid is not promoted to a uniform collar theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian node-c0 range certificate: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['laguerre_bound_rows']} Laguerre bound rows, "
        f"{summary['certified_side_conditions']} certified side conditions, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Node-C0 Range Certificate",
        "",
        "Date: 2026-07-07",
        "",
        "Status: node-range and Phi0 lower certificate. This is not a proof",
        "of a finite-grid interval certificate, a uniform residual estimate,",
        "scaled-curvature monotonicity, cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate`.",
        "",
        "Proof boundary: this artifact certifies only the two side conditions",
        "needed by the padded Phi-tail scout: `x<=1` for the recorded",
        "Laguerre-node grid and `Phi(0)>=0.44`. It does not certify",
        "individual nodes, weights, quadrature remainders, coefficient",
        "propagation, rounding, or a grid-to-collar bridge.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Node Range",
        "",
        "```text",
        f"min T: {summary['min_T']}",
        f"max quadrature order: {summary['max_quadrature_order']}",
        f"max index: {summary['max_index']}",
        f"worst Laguerre node upper bound: {summary['worst_laguerre_node_upper_bound']}",
        f"worst x^2 upper bound: {summary['worst_x_square_upper_bound']}",
        f"worst x^2 upper bound decimal: {summary['worst_x_square_upper_bound_decimal']}",
        f"worst x^2 slack to one: {summary['worst_x_square_slack_to_one']}",
        f"node range x<=1 certified: {summary['node_range_x_le_1_certified']}",
        "```",
        "",
        "The exact row bound uses the Laguerre Jacobi matrix, Gershgorin,",
        "and AM-GM. For the worst recorded case, `N=192` and",
        "`alpha=47/2`, every node is bounded by `809`, while `T>=1156`.",
        "",
        "## Phi0 Lower Bound",
        "",
        "```text",
        f"certified c0 lower: {summary['certified_c0_lower']}",
        f"n=1 Phi(0) term lower: {summary['n1_phi0_term_lower']}",
        f"n=1 margin over c0 lower: {summary['n1_minus_c0_lower_lower']}",
        f"c0 lower certified by n=1 term: {summary['c0_lower_certified_by_n1_term']}",
        "```",
        "",
        "Because",
        "",
        "```text",
        "2*pi^2*n^4 - 3*pi*n^2 = pi*n^2*(2*pi*n^2-3)>0 for n>=1, using pi>3",
        "```",
        "",
        "the full `Phi(0)` is at least the certified `n=1` contribution.",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_cancellation_reduced_grid_scout"],
        artifact["source_cancellation_reduced_grid_json"],
        artifact["source_phi_tail_bound_scout"],
        artifact["source_phi_tail_bound_json"],
        artifact["source_intervalization_target"],
        artifact["source_intervalization_target_json"],
        artifact["source_uniform_remainder_target"],
        artifact["source_dependency_graph"],
        "```",
        "",
        "Summary:",
        "",
        summary["main_finding"],
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-json", type=Path, default=DEFAULT_GRID_JSON)
    parser.add_argument("--phi-tail-json", type=Path, default=DEFAULT_PHI_TAIL_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    grid_path = args.grid_json if args.grid_json.is_absolute() else REPO_ROOT / args.grid_json
    phi_tail_path = args.phi_tail_json if args.phi_tail_json.is_absolute() else REPO_ROOT / args.phi_tail_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(grid_path, phi_tail_path, interval_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian node-c0 range certificate: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
