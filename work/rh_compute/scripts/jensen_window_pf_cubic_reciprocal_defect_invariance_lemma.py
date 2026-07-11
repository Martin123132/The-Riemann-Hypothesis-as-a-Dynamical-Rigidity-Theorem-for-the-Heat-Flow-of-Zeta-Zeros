#!/usr/bin/env python3
"""Build the reciprocal-defect form of cubic Jensen heat invariance."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))
if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

import flint  # noqa: E402

import jensen_window_pf_monotone_contraction_stress as stress  # noqa: E402
import jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate as m100  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.md"
DEFAULT_PRECISION_BITS = 1024


@dataclass(frozen=True)
class LemmaRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def lower_text(value: flint.arb, digits: int = 50) -> str:
    return value.lower().str(digits, radius=False)


def reciprocal_defect(x: flint.arb) -> flint.arb:
    defect = 1 - x
    if not arb_positive(defect):
        raise RuntimeError(f"contraction defect is not strictly positive: {defect}")
    return flint.arb(1) / defect.sqrt()


def finite_m100_diagnostics() -> dict:
    values, sources = m100.merged_coefficients()
    contractions = {
        k: values[k + 1] * values[k - 1] / values[k] ** 2 for k in range(1, 320)
    }
    q = {k: reciprocal_defect(contractions[k]) for k in contractions}
    margins = {k: 1 - (q[k + 1] - q[k]) for k in range(1, 319)}
    failed = [k for k, value in margins.items() if not arb_positive(value)]
    if failed:
        raise RuntimeError(f"lambda=-100 cubic reciprocal-defect margins failed: {failed[:10]}")
    minimum_k = min(margins, key=lambda k: margins[k].lower())
    selected_k = (1, 2, 100, 219, 220, 244, 245, 299, 300, 318)
    return {
        "lambda": "-100",
        "coefficient_range": [0, 320],
        "increment_k_range": [1, 318],
        "certified_positive_margins": len(margins),
        "failed_or_inconclusive_margins": len(failed),
        "minimum_margin_at_k": minimum_k,
        "minimum_margin_ball": str(margins[minimum_k]),
        "minimum_margin_lower": lower_text(margins[minimum_k]),
        "selected_margins": [
            {"k": k, "margin_ball": str(margins[k]), "margin_lower": lower_text(margins[k])}
            for k in selected_k
        ],
        "merged_sources": sources,
    }


def finite_nonnegative_diagnostics() -> dict:
    balls, samples, labels = stress.load_enclosures(list(stress.DEFAULT_ENCLOSURE_JSONL))
    rows = []
    global_minimum: tuple[flint.arb, str, int] | None = None
    total = 0
    for lam in sorted(labels):
        contractions = {}
        for k in range(1, 64):
            a_prev = balls[(lam, k - 1)]
            a_cur = balls[(lam, k)]
            a_next = balls[(lam, k + 1)]
            contractions[k] = a_next * a_prev / a_cur**2
        q = {k: reciprocal_defect(contractions[k]) for k in contractions}
        margins = {k: 1 - (q[k + 1] - q[k]) for k in range(1, 63)}
        failed = [k for k, value in margins.items() if not arb_positive(value)]
        if failed:
            raise RuntimeError(f"lambda={labels[lam]} reciprocal-defect margins failed: {failed[:10]}")
        minimum_k = min(margins, key=lambda k: margins[k].lower())
        candidate = (margins[minimum_k], labels[lam], minimum_k)
        if global_minimum is None or candidate[0].lower() < global_minimum[0].lower():
            global_minimum = candidate
        total += len(margins)
        rows.append(
            {
                "lambda": labels[lam],
                "increment_rows": len(margins),
                "certified_positive_margins": len(margins),
                "failed_or_inconclusive_margins": 0,
                "minimum_margin_at_k": minimum_k,
                "minimum_margin_ball": str(margins[minimum_k]),
                "minimum_margin_lower": lower_text(margins[minimum_k]),
            }
        )
    assert global_minimum is not None
    return {
        "lambda_count": len(rows),
        "coefficient_max_index": 64,
        "increment_rows": total,
        "certified_positive_margins": total,
        "failed_or_inconclusive_margins": 0,
        "global_minimum": {
            "lambda": global_minimum[1],
            "k": global_minimum[2],
            "margin_ball": str(global_minimum[0]),
            "margin_lower": lower_text(global_minimum[0]),
        },
        "lambda_rows": rows,
        "source_enclosures": [path.relative_to(REPO_ROOT).as_posix() for path in stress.DEFAULT_ENCLOSURE_JSONL],
    }


def build_payload() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    m100_diag = finite_m100_diagnostics()
    nonnegative_diag = finite_nonnegative_diagnostics()
    rows = [
        LemmaRow(
            id="crdi_01_reciprocal_defect_coordinate",
            role="exact_definition",
            readiness="available_exact",
            claim="Strict upper contraction walls permit reciprocal-square-root defect coordinates.",
            formula="d_k=1-x_k, q_k=d_k^(-1/2)",
            proof_boundary="Coordinate definition only; x_k<1 must be supplied by the coefficient theorem.",
        ),
        LemmaRow(
            id="crdi_02_cubic_discriminant_factorization",
            role="exact_identity",
            readiness="available_exact",
            claim="The normalized cubic frontier factors into four linear terms in consecutive reciprocal defects.",
            formula="q_k^4*q_(k+1)^4*F=(q_k-q_(k+1)-1)*(q_k-q_(k+1)+1)*(q_k+q_(k+1)-1)*(q_k+q_(k+1)+1)",
            proof_boundary="Exact cubic discriminant algebra only.",
        ),
        LemmaRow(
            id="crdi_03_cubic_hyperbolicity_equivalence",
            role="exact_equivalence",
            readiness="available_exact",
            claim="Inside the increasing-contraction cone, cubic Jensen hyperbolicity is exactly the unit-Lipschitz upper bound on q.",
            formula="Disc(J_(3,n))>=0 iff F<=0 iff 0<=q_(n+2)-q_(n+1)<=1",
            proof_boundary="Degree-3 equivalence only; not an all-degree criterion.",
        ),
        LemmaRow(
            id="crdi_04_frontier_parameterization",
            role="exact_identity",
            readiness="available_exact",
            claim="The nontrivial cubic frontier is the rank-two reciprocal-square defect family.",
            formula="t=sqrt(1-x): x=1-t^2, y=(1+2*t)/(1+t)^2, d_j=t^2/(1+(j-1)*t)^2",
            proof_boundary="Exact frontier parameterization for 0<t<1 only.",
        ),
        LemmaRow(
            id="crdi_05_heat_frontier_factorization",
            role="exact_identity",
            readiness="available_exact",
            claim="Modulo the cubic frontier, the heat derivative has a shift-independent threshold on the next contraction.",
            formula="partial_lambda F/r_(n+1)=C_n(t)*(z-Z(t)), C_n(t)>0, Z(t)=(1+t)*(1+3*t)/(1+2*t)^2",
            proof_boundary="Exact boundary vector-field factorization only.",
        ),
        LemmaRow(
            id="crdi_06_neighbor_inward_condition",
            role="exact_boundary_condition",
            readiness="available_exact",
            claim="At a saturated cubic increment, the vector field points inward exactly when the next cubic increment is not saturated outward.",
            formula="q_(k+1)-q_k=1 => d/dlambda(q_(k+1)-q_k)<=0 iff q_(k+2)-q_(k+1)<=1",
            proof_boundary="One-boundary inward condition; it assumes the next coordinate exists.",
        ),
        LemmaRow(
            id="crdi_07_conditional_infinite_invariance",
            role="conditional_exact_lemma",
            readiness="ready_to_apply",
            claim="A first-crossing maximum principle propagates all cubic inequalities once the full ratio cone, a cubic entry time, and uniform reciprocal-defect increment tails are supplied.",
            formula="sup_k(q_(k+1)-q_k-1)<=0 is forward invariant if q_(k+1)-q_k->0 uniformly on compact lambda intervals",
            proof_boundary="Conditional infinite maximum principle only; the reciprocal-defect tail and all-k entry theorem remain separate obligations.",
        ),
        LemmaRow(
            id="crdi_08_lambda_m100_prefix_certificate",
            role="interval_certificate",
            readiness="finite_validated",
            claim="The repaired lambda=-100 coefficients satisfy every reciprocal-defect cubic margin through k=318.",
            formula="1-(q_(k+1)-q_k)>0, 1<=k<=318",
            proof_boundary="Finite prefix with A_0..A_320 only; no analytic k>=319 tail.",
            diagnostics=m100_diag,
        ),
        LemmaRow(
            id="crdi_09_nonnegative_grid_certificate",
            role="interval_certificate",
            readiness="finite_validated",
            claim="All five cached nonnegative heat times satisfy every available reciprocal-defect cubic margin.",
            formula="1-(q_(k+1)-q_k)>0, 1<=k<=62",
            proof_boundary="Finite five-lambda, k<=62 interval evidence only.",
            diagnostics=nonnegative_diag,
        ),
        LemmaRow(
            id="crdi_10_open_tail_and_higher_degree_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The next analytic task is an all-k reciprocal-defect increment tail at lambda=-100 and on compact forward heat intervals; higher degrees still require stronger minor coordinates.",
            formula="prove q_(k+1)-q_k->0 uniformly and 1-(q_(k+1)-q_k)>0 for k>=319 at lambda=-100",
            proof_boundary="Open cubic-tail theorem and all-degree handoff; not PF-infinity, RH, or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_cubic_reciprocal_defect_invariance_lemma",
        "date": "2026-07-10",
        "status": "exact cubic reciprocal-defect invariance reduction with finite entry certificates",
        "proof_boundary": (
            "This artifact gives an exact degree-3 coordinate equivalence, boundary factorization, and conditional infinite maximum principle, plus finite Arb entry evidence. "
            "It does not prove the required all-k reciprocal-defect tail, all-degree Jensen hyperbolicity, PF-infinity, RH, or Lambda <= 0."
        ),
        "source_hierarchy": "outputs/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.md",
        "source_flow_cone": "outputs/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.md",
        "source_m100_entry": "outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md",
        "source_rank_two_family": "outputs/jensen_window_pf_rank_two_boundary_family_lemma.md",
        "diagnostics": {
            "precision_bits": DEFAULT_PRECISION_BITS,
            "lambda_minus_100": m100_diag,
            "nonnegative_grid": nonnegative_diag,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_coordinate_rows": 6,
            "conditional_invariance_rows": 1,
            "finite_certificate_rows": 2,
            "open_handoff_rows": 1,
            "lambda_minus_100_positive_margins": m100_diag["certified_positive_margins"],
            "nonnegative_grid_positive_margins": nonnegative_diag["certified_positive_margins"],
            "failed_or_inconclusive_margins": (
                m100_diag["failed_or_inconclusive_margins"]
                + nonnegative_diag["failed_or_inconclusive_margins"]
            ),
            "ready_to_apply_rows": 1,
            "target_closing": False,
            "main_finding": (
                "Cubic Jensen hyperbolicity becomes the unit-increment cone q_(k+1)-q_k<=1 for reciprocal square-root defects. "
                "At a saturated boundary the exact heat flow points inward precisely when the next shift obeys the same inequality, enabling a conditional infinite maximum principle. "
                "Arb certifies 318 lambda=-100 prefix margins and 310 nonnegative-grid margins, while the all-k reciprocal-defect tail remains open."
            ),
        },
        "invariants": [
            "The reciprocal-defect coordinate is used only where 0<x_k<1.",
            "The cubic equivalence assumes increasing contractions, supplied separately by the full ratio cone.",
            "Finite prefix evidence is not promoted to an all-k entry theorem.",
            "The maximum principle is conditional on a uniform spatial tail for reciprocal-defect increments.",
            "Degree-3 invariance is not an all-degree Jensen theorem.",
        ],
    }


def result_line(payload: dict) -> str:
    summary = payload["summary"]
    return (
        "validated Jensen-window PF cubic reciprocal-defect invariance lemma: "
        f"{summary['rows']} rows, 0 issues, {summary['exact_coordinate_rows']} exact coordinate rows, "
        f"{summary['conditional_invariance_rows']} conditional maximum principle, "
        f"{summary['lambda_minus_100_positive_margins']} lambda=-100 prefix margins, "
        f"{summary['nonnegative_grid_positive_margins']} nonnegative-grid margins, "
        f"{summary['failed_or_inconclusive_margins']} failed or inconclusive, "
        f"{summary['open_handoff_rows']} open tail handoff"
    )


def write_note(payload: dict, path: Path) -> None:
    m100_diag = payload["diagnostics"]["lambda_minus_100"]
    grid_diag = payload["diagnostics"]["nonnegative_grid"]
    lines = [
        "# Jensen-Window PF Cubic Reciprocal-Defect Invariance Lemma",
        "",
        "Date: 2026-07-10",
        "",
        "Status: exact cubic reciprocal-defect invariance reduction with finite entry",
        "certificates. This is not a proof of all-degree Jensen hyperbolicity, PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.json",
        "python work/rh_compute/scripts/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(payload),
        "```",
        "",
        "## Exact Cubic Coordinate",
        "",
        "Set `d_k=1-x_k` and `q_k=d_k^(-1/2)`. For the normalized cubic",
        "frontier `F=x^2*y^2-6*x*y+4*x+4*y-3`, exact factorization gives",
        "",
        "```text",
        "q_k^4*q_(k+1)^4*F",
        " =(q_k-q_(k+1)-1)*(q_k-q_(k+1)+1)",
        "  *(q_k+q_(k+1)-1)*(q_k+q_(k+1)+1).",
        "```",
        "",
        "Inside the full ratio cone, `q_(k+1)>=q_k>=1`. Since the cubic",
        "discriminant is `-27*x^2*F`, degree-3 hyperbolicity is exactly",
        "",
        "```text",
        "0 <= q_(k+1)-q_k <= 1.",
        "```",
        "",
        "## Boundary Flow",
        "",
        "At the nontrivial boundary, write `t=sqrt(1-x)`. Then",
        "",
        "```text",
        "x=1-t^2,",
        "y=(1+2*t)/(1+t)^2,",
        "d_j=t^2/(1+(j-1)*t)^2.",
        "```",
        "",
        "This is exactly the reciprocal-square defect law from the rank-two boundary",
        "family. Substitution into the heat ratio ODE gives",
        "",
        "```text",
        "partial_lambda F/r_(n+1)=C_n(t)*(z-Z(t)),",
        "C_n(t)=8*t^3*(2*n+7)*(1-t)*(1+2*t)^2/(1+t)^3 > 0,",
        "Z(t)=(1+t)*(1+3*t)/(1+2*t)^2,",
        "1-Z(t)=t^2/(1+2*t)^2.",
        "```",
        "",
        "Therefore, when `q_(k+1)-q_k=1`, the cubic frontier points inward",
        "exactly when `q_(k+2)-q_(k+1)<=1`. The threshold is independent of the",
        "shift. If the increments tend uniformly to zero in the spatial tail, a",
        "finite-active-set first-crossing argument propagates the whole cubic cone.",
        "",
        "## Finite Arb Entry",
        "",
        f"At `lambda=-100`, all `{m100_diag['certified_positive_margins']}` repaired-prefix margins",
        "through `k=318` are strictly positive. The minimum is",
        "",
        "```text",
        f"k={m100_diag['minimum_margin_at_k']}: {m100_diag['minimum_margin_ball']}",
        "```",
        "",
        f"Across the five cached nonnegative heat times, all `{grid_diag['certified_positive_margins']}`",
        "available margins are strictly positive. The global finite minimum is",
        "",
        "```text",
        f"lambda={grid_diag['global_minimum']['lambda']}, k={grid_diag['global_minimum']['k']}:",
        grid_diag["global_minimum"]["margin_ball"],
        "```",
        "",
        "## Remaining Theorem",
        "",
        "The cubic propagation route now needs one specific analytic input:",
        "",
        "```text",
        "q_(k+1)-q_k -> 0 uniformly on compact lambda intervals,",
        "and q_(k+1)-q_k<1 for every k>=319 at lambda=-100.",
        "```",
        "",
        "That would promote the existing finite entry and exact boundary algebra to all",
        "degree-3 shifts. Higher degrees still need additional minor coordinates, so this",
        "is a controlled cubic advance rather than the missing all-order bridge.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(payload, args.note)
    print(result_line(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
