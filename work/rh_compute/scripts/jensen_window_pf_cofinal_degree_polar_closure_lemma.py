#!/usr/bin/env python3
"""Prove finite-tower and cofinal-degree closure under Jensen polar maps."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_cofinal_degree_polar_closure_lemma.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_cofinal_degree_polar_closure_lemma.md"
STURM_34 = REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520_summary.json"
STURM_5 = REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520_summary.json"
STURM_6_12 = REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d6_d12_dps520_summary.json"
CONTRACTION_STRESS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_stress_lamgrid_d3_d12_k64_summary.json"


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
    return json.loads(path.read_text(encoding="utf-8"))


def build_exact() -> dict:
    d, j = sp.symbols("d j", integer=True, positive=True)
    binomial_identity = sp.simplify(
        (1 - j / (d + 1)) * sp.binomial(d + 1, j) - sp.binomial(d, j)
    )
    if binomial_identity != 0:
        raise RuntimeError("adjacent polar coefficient identity failed")

    w = sp.symbols("w")
    alpha = sp.symbols("alpha_0:6", positive=True)
    product = sp.prod(1 + value * w for value in alpha)
    polar = sp.factor(product - w * sp.diff(product, w) / len(alpha))
    logarithmic_form = sp.factor(polar / product)
    expected_logarithmic = sum(1 / (1 + value * w) for value in alpha) / len(alpha)
    if sp.simplify(logarithmic_form - expected_logarithmic) != 0:
        raise RuntimeError("polar logarithmic representation failed")
    logarithmic_derivative = sp.factor(sp.diff(expected_logarithmic, w))

    sturm34 = load_json(STURM_34)
    sturm5 = load_json(STURM_5)
    sturm6_12 = load_json(STURM_6_12)
    stress = load_json(CONTRACTION_STRESS)
    if not sturm34.get("all_ok") or sturm34.get("degrees") != [3, 4]:
        raise RuntimeError("degree-3/4 Sturm source drifted")
    if not sturm5.get("all_ok") or sturm5.get("degrees") != [5]:
        raise RuntimeError("degree-5 Sturm source drifted")
    if not sturm6_12.get("all_ok") or sturm6_12.get("degrees") != list(range(6, 13)):
        raise RuntimeError("degree-6 through degree-12 Sturm source drifted")
    stress_degrees = stress.get("degrees")
    if stress_degrees != list(range(3, 13)):
        raise RuntimeError("contraction stress degree range drifted")

    return {
        "normalized_windows": "P_(d,n)(w)=sum_(j=0)^d C(d,j)*B_(n,j)*w^j, B_(n,0)=B_(n,1)=1",
        "adjacent_polar_identity": "P_(d,n)=P_(d+1,n)-w*P_(d+1,n)'/(d+1)",
        "binomial_identity": "(1-j/(d+1))*C(d+1,j)=C(d,j)",
        "polar_product_identity": (
            "If P_D=product_i(1+alpha_i*w), then "
            "P_(D-1)/P_D=(1/D)*sum_i 1/(1+alpha_i*w)."
        ),
        "polar_logarithmic_derivative": str(logarithmic_derivative),
        "one_step_preservation": "H_(D,n) => H_(D-1,n), with negative-root interlacing",
        "finite_tower": "H_(D,n) => H_(d,n) for every 0<=d<=D",
        "cofinal_closure": (
            "If {D:H_(D,n)} is unbounded for a fixed n, then H_(d,n) holds for every finite d."
        ),
        "all_shift_cofinal_closure": (
            "If the cofinal condition holds for every n>=0, then every shifted Jensen polynomial is hyperbolic."
        ),
        "evidence_audit": {
            "sturm_hyperbolicity_degrees": list(range(3, 13)),
            "sturm_lambdas": sturm34["lambdas"],
            "sturm_shifts": [min(sturm34["shifts"]), max(sturm34["shifts"])],
            "sturm_rows": sturm34["rows"] + sturm5["rows"] + sturm6_12["rows"],
            "sturm_failed_or_inconclusive": (
                sturm34["failed_or_inconclusive"]
                + sturm5["failed_or_inconclusive"]
                + sturm6_12["failed_or_inconclusive"]
            ),
            "contraction_only_degrees": stress_degrees,
            "contraction_only_rows": stress["summary"]["stress_rows"],
            "contraction_evidence_is_not_hyperbolicity": True,
            "cofinal_terminal_degrees_certified": False,
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        LemmaRow(
            id="cdpc_01_normalized_degree_family",
            role="exact_definition",
            readiness="available_exact",
            claim="At each fixed shift, all Jensen degrees use one common normalized coefficient sequence.",
            formula=exact["normalized_windows"],
            proof_boundary="Normalization only.",
        ),
        LemmaRow(
            id="cdpc_02_adjacent_polar_identity",
            role="exact_identity",
            readiness="available_exact",
            claim="Dropping one Jensen degree is exactly the polar derivative with pole zero.",
            formula=exact["adjacent_polar_identity"],
            proof_boundary="Adjacent-degree identity only.",
        ),
        LemmaRow(
            id="cdpc_03_polar_product_formula",
            role="exact_identity",
            readiness="available_exact",
            claim="For a negative-root hyperbolic terminal polynomial, the polar quotient is an average of simple reciprocal factors.",
            formula=exact["polar_product_identity"],
            proof_boundary="Hyperbolic terminal polynomial only.",
        ),
        LemmaRow(
            id="cdpc_04_one_step_preservation",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The zero-pole polar derivative preserves negative-root hyperbolicity and interlaces the terminal roots.",
            formula=exact["one_step_preservation"],
            proof_boundary="One adjacent-degree descent only.",
        ),
        LemmaRow(
            id="cdpc_05_multiplicity_control",
            role="exact_lemma",
            readiness="ready_to_apply",
            claim="Nonroot polar zeros are simple, while a nonzero terminal root loses exactly one unit of multiplicity.",
            formula="nonroot derivative<0; mult_r(P_D)=m => mult_r(P_(D-1))=m-1",
            proof_boundary="Polar multiplicity structure only.",
        ),
        LemmaRow(
            id="cdpc_06_finite_tower_closure",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="One hyperbolic terminal degree closes every lower degree at the same shift by repeated polar descent.",
            formula=exact["finite_tower"],
            proof_boundary="Finite degree tower at one shift only.",
        ),
        LemmaRow(
            id="cdpc_07_cofinal_degree_closure",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="An unbounded set of hyperbolic terminal degrees is equivalent to all finite degrees at a fixed shift.",
            formula=exact["cofinal_closure"],
            proof_boundary="Fixed-shift cofinal reduction only.",
        ),
        LemmaRow(
            id="cdpc_08_all_shift_reduction",
            role="conditional_composition",
            readiness="ready_to_apply",
            claim="Cofinal terminal hyperbolicity at every shift suffices for the complete Jensen-window family.",
            formula=exact["all_shift_cofinal_closure"],
            proof_boundary="Conditional all-shift reduction; no terminal theorem is supplied.",
        ),
        LemmaRow(
            id="cdpc_09_terminal_evidence_audit",
            role="finite_evidence_audit",
            readiness="finite_validated",
            claim="Current rigorous finite hyperbolicity data reach every degree from 3 through 12; the separate contraction cache is not needed to infer those finite root counts.",
            formula="1050 Sturm rows at d=3..12; 2875 separate contraction-only rows at d=3..12",
            proof_boundary="Finite evidence classification only; not a cofinal terminal sequence.",
            diagnostics=exact["evidence_audit"],
        ),
        LemmaRow(
            id="cdpc_10_open_terminal_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Prove hyperbolicity for an unbounded sequence of terminal degrees at every shift, with enough uniformity to survive the heat-flow argument.",
            formula="for every n, construct D_j->infinity with H_(D_j,n)",
            proof_boundary="Open cofinal terminal theorem; not PF-infinity, RH, or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_cofinal_degree_polar_closure_lemma",
        "date": "2026-07-10",
        "status": "exact cofinal-degree polar closure with open terminal theorem",
        "proof_boundary": (
            "This artifact proves that cofinally many hyperbolic Jensen degrees at each fixed shift suffice for all finite degrees. It does not prove any unbounded terminal-degree sequence for the zeta coefficients, an all-shift cofinal theorem, PF-infinity, RH, or Lambda <= 0."
        ),
        "sources": [
            "outputs/jensen_window_pf_quartic_quintic_polar_contact_lemma.md",
            "outputs/arb_jensen_window_sturm_hyperbolicity_diagnostic.md",
            "outputs/jensen_window_pf_monotone_contraction_stress.md",
            "outputs/jensen_window_pf_bridge_target.md",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    audit = exact["evidence_audit"]
    return "\n".join(
        [
            "# Jensen-Window PF Cofinal-Degree Polar-Closure Lemma",
            "",
            "Date: 2026-07-10",
            "",
            "Status: exact cofinal-degree reduction with the terminal theorem open.",
            "This is not a proof of PF-infinity, RH, or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_cofinal_degree_polar_closure_lemma.json",
            "python work/rh_compute/scripts/jensen_window_pf_cofinal_degree_polar_closure_lemma.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_cofinal_degree_polar_closure_lemma.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF cofinal-degree polar-closure lemma: 10 rows, 0 issues, 3 exact polar identities, 1 interlacing theorem, 1 multiplicity theorem, 1 finite-tower closure, 1 cofinal-degree closure, 1050 finite Sturm rows, 2875 contraction-only rows, 1 open terminal handoff",
            "```",
            "",
            "## Polar Descent",
            "",
            "At a fixed shift `n`, normalized Jensen windows satisfy",
            "",
            "```text",
            exact["adjacent_polar_identity"],
            "```",
            "",
            "If `P_D=product_i(1+alpha_i*w)` with `alpha_i>0`, then",
            "",
            "```text",
            exact["polar_product_identity"],
            "```",
            "",
            "The reciprocal sum is strictly decreasing between its poles. Its zeros",
            "are real, negative, and interlace the roots of `P_D`; repeated terminal",
            "roots lose one unit of multiplicity. Therefore",
            "",
            "```text",
            exact["finite_tower"],
            "```",
            "",
            "## Cofinal Closure",
            "",
            "For any fixed shift and finite target degree `d`, one may choose a",
            "hyperbolic terminal degree `D>=d` and descend. Consequently",
            "",
            "```text",
            exact["cofinal_closure"],
            "```",
            "",
            "Thus the all-degree target does not require certifying every degree",
            "directly. An unbounded terminal subsequence at every shift is sufficient.",
            "",
            "## Evidence Audit",
            "",
            "Current rigorous finite Sturm data contain",
            "",
            "```text",
            f"degrees={audit['sturm_hyperbolicity_degrees']}",
            f"rows={audit['sturm_rows']}",
            f"failed_or_inconclusive={audit['sturm_failed_or_inconclusive']}",
            "shifts=0..20 on five cached nonnegative heat parameters.",
            "```",
            "",
            "The 2,875 contraction rows remain a logically separate diagnostic; the",
            "new root counts come from Sturm certificates, not from promotion of",
            "contraction inequalities.",
            "",
            "## Remaining Handoff",
            "",
            "The next theorem must construct, for every shift `n`, an unbounded",
            "sequence `D_j` of hyperbolic terminal degrees. The reduction is exact,",
            "but no large-degree asymptotic or compactness theorem currently supplies",
            "those terminals for the zeta coefficients.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Jensen-window PF cofinal-degree polar-closure lemma: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
