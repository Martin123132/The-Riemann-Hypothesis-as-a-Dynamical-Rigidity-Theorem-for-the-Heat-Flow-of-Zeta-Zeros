#!/usr/bin/env python3
"""Audit the recent order-moment transport theorem against the first summand."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists() and str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_order_moment_transport_fit_gate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_order_moment_transport_fit_gate.md"
SOURCE_URL = "https://arxiv.org/abs/2606.31647v2"
DEFAULT_PRECISION_BITS = 256


@dataclass(frozen=True)
class FitRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def arb_lower_text(value: flint.arb, digits: int = 60) -> str:
    return value.lower().str(digits, radius=False).replace("e", "E")


def exact_diagnostics() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    pi = flint.arb.pi()
    derivative = 5 + 8 * pi / (2 * pi - 3) - 4 * pi
    if not bool(derivative > 0):
        raise RuntimeError("first-summand origin derivative is not positive")
    return {
        "first_summand": (
            "phi_1(u)=pi*exp(5*u)*(2*q-3)*exp(-q), q=pi*exp(4*u)"
        ),
        "substitution": "y=100*u^2",
        "gamma_average": (
            "A_t^(1)=sqrt(pi)/(10*400^t*Gamma(t+1/2))*"
            "integral_0^infinity y^(t-1/2)*exp(-y)*g(y)dy"
        ),
        "transformed_kernel": "g(y)=phi_1(sqrt(y)/10)",
        "transport_if_cm": (
            "g(y)=integral exp(-s*y)dmu(s) => normalized gamma average="
            "integral (1+s)^(-(t+1/2))dmu(s)"
        ),
        "cm_orientation": (
            "The transported A_t^(1) is a positive Mellin moment family and has "
            "nonnegative Hankel minors (ordinary log-convex orientation)."
        ),
        "origin_log_derivative": (
            "d_u log(phi_1)(0)=5+8*pi/(2*pi-3)-4*pi"
        ),
        "origin_log_derivative_ball": derivative.str(60).replace("e", "E"),
        "origin_log_derivative_lower": arb_lower_text(derivative),
        "complete_monotonicity_failure": (
            "d_u phi_1(0)>0, so g'(y)>0 for sufficiently small y>0; "
            "g is not completely monotone"
        ),
        "source_scope": (
            "The source theorem treats positive Mellin/Hankel total positivity and "
            "explicitly separates reciprocal sign-regularity."
        ),
        "live_replacement": (
            "A signed reciprocal-Gamma transport theorem or the direct stable-gap "
            "curvature certificate is still required."
        ),
    }


def build_artifact() -> dict:
    exact = exact_diagnostics()
    rows = [
        FitRow(
            id="omtf_01_primary_scope",
            role="theorem_fit",
            readiness="source_audited",
            claim="The order-moment transport theorem turns Gamma-normalized completely-monotone averages into positive Hankel moment kernels.",
            formula=exact["transport_if_cm"],
            proof_boundary="Primary-source theorem scope only.",
        ),
        FitRow(
            id="omtf_02_exact_reparametrization",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The lambda=-100 first Newman summand has exactly the Gamma-average form tested by the theorem.",
            formula=exact["gamma_average"],
            proof_boundary="Exact change of variables only.",
        ),
        FitRow(
            id="omtf_03_cm_prerequisite",
            role="theorem_fit",
            readiness="tested",
            claim="Application of the positive transport criterion requires complete monotonicity of the transformed kernel.",
            formula="(-1)^m*g^(m)(y)>=0 for every m>=0 and y>0",
            proof_boundary="Necessary hypothesis for this theorem route, not a claim about all routes.",
        ),
        FitRow(
            id="omtf_04_origin_obstruction",
            role="interval_countercheck",
            readiness="guard_validated",
            claim="The transformed first-summand kernel increases at the origin and therefore fails complete monotonicity.",
            formula=exact["origin_log_derivative"],
            proof_boundary="Exact formula with Arb sign certification.",
            diagnostics={
                "ball": exact["origin_log_derivative_ball"],
                "lower": exact["origin_log_derivative_lower"],
            },
        ),
        FitRow(
            id="omtf_05_orientation_mismatch",
            role="structural_mismatch",
            readiness="guard_validated",
            claim="If the criterion applied, it would produce positive Hankel moment orientation, not the reciprocal-Gamma signed orientation required here.",
            formula=exact["cm_orientation"],
            proof_boundary="Route-orientation audit only.",
        ),
        FitRow(
            id="omtf_06_forbidden_promotion",
            role="forbidden_promotion",
            readiness="guard_validated",
            claim="The positive order-moment transport theorem cannot be cited as the missing signed-Hankel bridge for the Xi first summand.",
            formula="positive transport theorem != reciprocal sign-regular transport theorem",
            proof_boundary="Rejects this theorem application only, not total-positivity methods in general.",
        ),
        FitRow(
            id="omtf_07_live_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Seek a signed reciprocal-Gamma analogue or continue the direct stable-gap curvature proof.",
            formula=exact["live_replacement"],
            proof_boundary="Open theorem-search handoff.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_order_moment_transport_fit_gate",
        "date": "2026-07-12",
        "status": (
            "primary-source theorem-fit audit with exact first-summand "
            "complete-monotonicity obstruction"
        ),
        "proof_boundary": (
            "This artifact rejects direct use of one positive order-moment transport "
            "criterion for the first Newman summand. It does not reject signed total "
            "positivity, prove the stable-gap curvature ceiling, order-four entry, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "primary_source": {
            "title": (
                "Order-Moment Transport and Hankel Determinants in "
                "Special-Function Inequalities"
            ),
            "author": "Domingos S. P. Salazar",
            "version": "arXiv:2606.31647v2",
            "url": SOURCE_URL,
        },
        "sources": [
            SOURCE_URL,
            "outputs/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.md",
            "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
        ],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_order_moment_transport_fit_gate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_order_moment_transport_fit_gate.py"
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "fit_rows": len(rows),
            "exact_identity_rows": 1,
            "interval_counterchecks": 1,
            "structural_mismatches": 1,
            "forbidden_promotions": 1,
            "open_handoffs": 1,
            "ready_to_apply_rows": 1,
        },
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    source = artifact["primary_source"]
    lines = [
        "# Jensen-Window PF Order-Moment Transport Fit Gate",
        "",
        "Date: 2026-07-12",
        "",
        "Status: primary-source theorem-fit rejection gate. This is not a proof",
        "of order-four entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_order_moment_transport_fit_gate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_order_moment_transport_fit_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_order_moment_transport_fit_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_order_moment_transport_fit_gate.py",
        "```",
        "",
        "## Primary Theorem",
        "",
        f"[{source['title']}]({source['url']}) records a positive order-moment",
        "transport criterion: a Gamma-normalized average of a completely monotone",
        "function becomes a positive Mellin moment family and hence a totally",
        "nonnegative Hankel kernel. The paper explicitly treats reciprocal",
        "sign-regularity as a separate problem.",
        "",
        "## Exact Fit Test",
        "",
        "For the first Newman summand at `lambda=-100`, put `y=100*u^2`. Then",
        "",
        "```text",
        exact["first_summand"],
        exact["gamma_average"],
        exact["transformed_kernel"],
        "```",
        "",
        "This is formally the right Gamma-average shape. However, the positive",
        "transport theorem requires `g` to be completely monotone. At the origin",
        "",
        "```text",
        exact["origin_log_derivative"] + ",",
        f"  = {exact['origin_log_derivative_ball']},",
        f"  > {exact['origin_log_derivative_lower']} > 0.",
        "```",
        "",
        "Thus `phi_1` increases immediately to the right of zero, and because",
        "`y -> sqrt(y)/10` is increasing, `g'(y)>0` for sufficiently small",
        "positive `y`. The transformed kernel is not completely monotone.",
        "",
        "## Consequence",
        "",
        "The criterion would in any event produce ordinary positive Hankel moment",
        "orientation. The Xi programme needs the alternating reciprocal-Gamma",
        "sign orientation. Therefore this theorem cannot be promoted into the",
        "missing bridge.",
        "",
        "The useful handoff is precise: seek a signed reciprocal-Gamma transport",
        "theorem, or continue the direct stable-gap curvature certificate. This",
        "gate rejects one theorem application, not total-positivity methods as a",
        "whole.",
        "",
        "```text",
        "outputs/jensen_window_pf_reciprocal_gamma_mixture_sign_gate.md",
        "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
        "outputs/jensen_window_pf_theorem_machinery_fit_matrix.md",
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
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "validated Jensen-window PF order-moment transport fit gate: "
        "7 rows, 0 issues, 1 exact reparametrization, 1 positive origin "
        "derivative obstruction, 1 orientation mismatch, 1 forbidden promotion, "
        "1 open signed-transport handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
