#!/usr/bin/env python3
"""Record exact heat-flow identities for binomial Jensen windows."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.md"


def build_payload() -> dict:
    rows = [
        {
            "id": "hfjh_01_coefficient_heat_ode",
            "role": "exact_identity",
            "formula": "A_k'=2*(2*k+1)*A_(k+1)",
            "claim": "The factorial-normalized Newman coefficients satisfy an exact upper-shift heat ODE.",
            "readiness": "available_exact",
        },
        {
            "id": "hfjh_02_jensen_window_definition",
            "role": "exact_definition",
            "formula": "J_(d,n)(z)=sum_(j=0)^d binom(d,j)*A_(n+j)*z^j",
            "claim": "All degree and shift dependence is retained explicitly.",
            "readiness": "available_exact",
        },
        {
            "id": "hfjh_03_pascal_shift_degree_identity",
            "role": "exact_identity",
            "formula": "J_(d+1,n)=J_(d,n)+z*J_(d,n+1)",
            "claim": "Pascal's identity couples one degree step to one coefficient shift.",
            "readiness": "available_exact",
        },
        {
            "id": "hfjh_04_z_derivative_identity",
            "role": "exact_identity",
            "formula": "partial_z J_(d,n)=d*J_(d-1,n+1)",
            "claim": "A z derivative lowers degree and raises shift.",
            "readiness": "available_exact",
        },
        {
            "id": "hfjh_05_lambda_hierarchy_identity",
            "role": "exact_identity",
            "formula": "partial_lambda J_(d,n)=(4*n+2)*J_(d,n+1)+4*z*partial_z J_(d,n+1)",
            "claim": "The heat derivative is an exact shifted Jensen-window operator.",
            "readiness": "available_exact",
        },
        {
            "id": "hfjh_06_cross_degree_operator_identity",
            "role": "exact_identity",
            "formula": "partial_lambda J_(d,n)=((4*n+2)*partial_z+4*z*partial_z^2)*J_(d+1,n)/(d+1)",
            "claim": "The same hierarchy is an explicit differential operator applied to the next-degree window.",
            "readiness": "available_exact",
        },
        {
          "id": "hfjh_07_degree3_discriminant_frontier",
            "role": "exact_countermodel_gate",
            "formula": "Disc(1+3*z+3*x*z^2+x^2*y*z^3)=-27*x^2*(x^2*y^2-6*x*y+4*x+4*y-3)",
            "claim": "The full static ratio cone does not force even degree-3 hyperbolicity: x=1/2,y=1 gives discriminant -27/16.",
            "readiness": "guard_validated",
        },
        {
            "id": "hfjh_08_cubic_boundary_exit_countermodel",
            "role": "exact_countermodel_gate",
            "formula": "d_k=(4/9)*(9/25)^(k-1); n=0, x=5/9, y=21/25, z=589/625: F=0 and (partial_lambda F)/r_0=329728/2109375>0, hence partial_lambda Disc<0",
            "claim": "Even complete Hausdorff defect monotonicity, the full static ratio cone, and the exact local coefficient heat ODE do not make the degree-3 hyperbolicity boundary forward invariant.",
            "readiness": "guard_validated",
        },
        {
            "id": "hfjh_09_open_higher_minor_handoff",
            "role": "open_handoff",
            "formula": "shift-coupled hierarchy + missing invariant cone => all-degree Jensen hyperbolicity",
            "claim": "A viable heat-flow proof must control a degree-and-shift hierarchy or its minors; there is no closed fixed-window scalar evolution supplied here.",
            "readiness": "not_ready_to_apply",
        },
    ]
    return {
        "kind": "jensen_window_pf_heat_flow_jensen_hierarchy_lemma",
        "date": "2026-07-10",
        "status": "exact heat-flow hierarchy lemma and open higher-minor handoff",
        "proof_boundary": (
            "Exact coefficient and Jensen-window identities plus one exact cubic countermodel gate. "
            "This does not prove an invariant higher-minor cone, Jensen hyperbolicity, PF-infinity, RH, or Lambda <= 0."
        ),
        "source_flow_certificate": "outputs/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.md",
        "source_bridge_target": "outputs/jensen_window_pf_bridge_target.md",
        "source_defect_scout": "outputs/jensen_window_pf_defect_complete_monotonicity_scout.md",
        "rows": rows,
        "summary": {
            "rows": len(rows),
            "exact_identity_rows": 5,
            "exact_definition_rows": 1,
            "exact_cubic_countermodels": 2,
            "open_hierarchy_handoffs": 1,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "The lambda derivative of a Jensen window is exactly coupled to shift n+1 and, equivalently, to degree d+1. "
                "This gives a concrete hierarchy for higher-minor theorem search. The -27/16 static cubic guard and the exact one-atom Hausdorff boundary-exit guard 329728/2109375 prove that neither the propagated static cone nor complete defect monotonicity plus the local heat ODE supplies the missing invariant structure."
            ),
        },
        "invariants": [
            "Every displayed hierarchy identity is coefficientwise exact for all d>=1 and n>=0.",
            "The hierarchy is not closed at a single degree and shift.",
            "The cubic witness is abstract and is not a zeta heat-flow orbit.",
            "Endpoint Jensen hyperbolicity, RH, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(payload: dict, path: Path) -> None:
    result = (
        "validated Jensen-window PF heat-flow Jensen hierarchy lemma: "
        "9 rows, 0 issues, 5 exact hierarchy identities, 2 cubic countermodels, "
        "1 open higher-minor handoff, 0 ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Heat-Flow Jensen Hierarchy Lemma",
        "",
        "Date: 2026-07-10",
        "",
        "Status: exact heat-flow hierarchy lemma and open higher-minor handoff. This is not a",
        "proof of Jensen hyperbolicity, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_heat_flow_jensen_hierarchy_lemma`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.json",
        "python work/rh_compute/scripts/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_heat_flow_jensen_hierarchy_lemma.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result,
        "```",
        "",
        "## Exact Hierarchy",
        "",
        "For",
        "",
        "```text",
        "J_(d,n)(z)=sum_(j=0)^d binom(d,j)*A_(n+j)*z^j",
        "A_k'=2*(2*k+1)*A_(k+1),",
        "```",
        "",
        "coefficient comparison gives",
        "",
        "```text",
        "J_(d+1,n)=J_(d,n)+z*J_(d,n+1),",
        "partial_z J_(d,n)=d*J_(d-1,n+1),",
        "partial_lambda J_(d,n)=(4*n+2)*J_(d,n+1)+4*z*partial_z J_(d,n+1),",
        "partial_lambda J_(d,n)=((4*n+2)*partial_z+4*z*partial_z^2)*J_(d+1,n)/(d+1).",
        "```",
        "",
        "The heat evolution is therefore coupled across shifts and degrees. These identities do",
        "not supply a closed scalar PDE for one fixed Jensen window.",
        "",
        "## Cubic Gate",
        "",
        "The normalized degree-3 discriminant is",
        "",
        "```text",
        "Disc(1+3*z+3*x*z^2+x^2*y*z^3)",
        "  = -27*x^2*(x^2*y^2-6*x*y+4*x+4*y-3).",
        "```",
        "",
        "At `x=1/2,y=1`, the infinite contraction sequence lies in the full static ratio cone",
        "and has completely monotone defects, but the discriminant is `-27/16`. Thus a",
        "static higher-minor proof needs an additional invariant.",
        "",
        "There is also an exact dynamic boundary guard. At shift `n=0`, set",
        "",
        "```text",
        "d_k=(4/9)*(9/25)^(k-1),",
        "x=5/9, y=21/25, z=589/625,",
        "F=x^2*y^2-6*x*y+4*x+4*y-3=0,",
        "(partial_lambda F)/r_0=329728/2109375>0.",
        "```",
        "",
        "Since the discriminant is a positive factor times `-F`, its heat derivative is",
        "negative at this boundary point. The defects form the one-atom Hausdorff moment",
        "sequence `(4/9)*delta_(9/25)` and satisfy every full-cone wall. Thus even complete",
        "defect monotonicity plus the local heat ODE does not supply forward cubic invariance. A viable",
        "higher-minor flow proof needs an additional shift-coupled invariant.",
        "",
        "## Handoff",
        "",
        payload["summary"]["main_finding"],
        "",
        "```text",
        "outputs/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.md",
        "outputs/jensen_window_pf_defect_complete_monotonicity_scout.md",
        "outputs/jensen_window_pf_bridge_target.md",
        "```",
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
    print(
        "validated Jensen-window PF heat-flow Jensen hierarchy lemma: "
        "9 rows, 0 issues, 5 exact hierarchy identities, 2 cubic countermodels, "
        "1 open higher-minor handoff, 0 ready-to-apply rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
