#!/usr/bin/env python3
"""Build the cancellation/zero-free wall for the Newman scaled layer."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem import (
    build_envelope,
    phase_exponent,
    required_scaled_time,
    transition_radius,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.md"
)

POLYMATH_SOURCE_URL = "https://arxiv.org/abs/1904.12438"
EXPONENT_PAIR_SOURCE_URL = "https://arxiv.org/abs/2306.05599"
BOURGAIN_SOURCE_URL = "https://arxiv.org/abs/1408.5794"
ZERO_FREE_SOURCE_URL = "https://arxiv.org/abs/2212.06867"
ONE_LINE_SOURCE_URL = "https://arxiv.org/abs/2312.09412"
RECIPROCAL_SOURCE_URL = "https://arxiv.org/abs/2405.04869"

C_STAR = Fraction(4_911_678_521, 1_933_561_194)
C_ZETA_WALL = Fraction(2, 1)
R_STAR = Fraction(125_662, 155_153)
KNOWN_BOURGAIN_THETA = Fraction(13, 84)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def fraction_record(value: Fraction) -> dict:
    return {
        "exact": f"{value.numerator}/{value.denominator}",
        "decimal": f"{float(value):.15f}",
    }


def target_phase(radius: Fraction, scaled_time: Fraction) -> Fraction:
    """Largest phase-sum exponent compatible with decay at this radius."""
    return radius / 2 + scaled_time * radius * (2 - radius) / 8


def build_frontier() -> dict:
    envelope = build_envelope()
    critical_time = Fraction(envelope["critical_scaled_time"]["exact"])
    critical_radius = Fraction(envelope["critical_radius"]["exact"])
    if critical_time != C_STAR or critical_radius != R_STAR:
        raise RuntimeError("known cancellation threshold drifted")

    transition_rows: list[dict] = []
    for upper_index in range(10, 0, -1):
        lower_index = upper_index - 1
        radius = transition_radius(upper_index, lower_index)
        known = phase_exponent(upper_index, radius)
        if known != phase_exponent(lower_index, radius):
            raise RuntimeError("transition phase identity failed")
        required = required_scaled_time(upper_index, radius)
        c2_target = target_phase(radius, C_ZETA_WALL)
        transition_rows.append(
            {
                "upper_index": upper_index,
                "lower_index": lower_index,
                "radius": fraction_record(radius),
                "known_phase_exponent": fraction_record(known),
                "required_scaled_time": fraction_record(required),
                "c2_target_phase_exponent": fraction_record(c2_target),
                "c2_phase_excess": fraction_record(known - c2_target),
                "is_current_maximum": required == C_STAR,
            }
        )

    endpoint_radius = Fraction(1, 1)
    endpoint_known = phase_exponent(0, endpoint_radius)
    endpoint_required = required_scaled_time(0, endpoint_radius)
    transition_rows.append(
        {
            "upper_index": 0,
            "lower_index": None,
            "radius": fraction_record(endpoint_radius),
            "known_phase_exponent": fraction_record(endpoint_known),
            "required_scaled_time": fraction_record(endpoint_required),
            "c2_target_phase_exponent": fraction_record(
                target_phase(endpoint_radius, C_ZETA_WALL)
            ),
            "c2_phase_excess": fraction_record(
                endpoint_known - target_phase(endpoint_radius, C_ZETA_WALL)
            ),
            "is_current_maximum": endpoint_required == C_STAR,
        }
    )

    required_values = [
        Fraction(row["required_scaled_time"]["exact"])
        for row in transition_rows
    ]
    if max(required_values) != C_STAR:
        raise RuntimeError("frontier no longer reproduces c_star")
    if sum(row["is_current_maximum"] for row in transition_rows) != 1:
        raise RuntimeError("frontier maximum is not unique")
    if any(Fraction(row["c2_phase_excess"]["exact"]) <= 0 for row in transition_rows):
        raise RuntimeError("known frontier unexpectedly reaches c=2")

    critical_known = phase_exponent(2, R_STAR)
    critical_c2_target = target_phase(R_STAR, C_ZETA_WALL)
    critical_deficit = critical_known - critical_c2_target
    symmetric_theta_cap = R_STAR * (2 - R_STAR) / 8
    if critical_known != Fraction(220_633, 310_306):
        raise RuntimeError("critical phase exponent drifted")

    return {
        "transition_rows": transition_rows,
        "critical_radius": fraction_record(R_STAR),
        "critical_known_phase_exponent": fraction_record(critical_known),
        "critical_c2_target_phase_exponent": fraction_record(critical_c2_target),
        "critical_c2_phase_deficit": fraction_record(critical_deficit),
        "critical_symmetric_theta_cap": fraction_record(symmetric_theta_cap),
        "known_bourgain_theta": fraction_record(KNOWN_BOURGAIN_THETA),
        "known_scaled_threshold": fraction_record(C_STAR),
        "zeta_line_scaled_wall": fraction_record(C_ZETA_WALL),
    }


def build_exact() -> dict:
    frontier = build_frontier()
    return {
        "scaled_geometry": (
            "L=log(x/(4*pi)), c=t*L, N=exp(L/2+o(1)), "
            "sigma=Re(s_*)=1/2+c/4+o(1), tau=|Im(s_*)|=2*pi*N^2*(1+o(1))"
        ),
        "weighted_block": (
            "If M=N^r and sup_I|sum_(n in I)n^(-i*tau)|<<N^(Phi(r)+eta), "
            "then the weighted dyadic block is "
            "<<N^(Phi(r)-r/2-c*r*(2-r)/8+eta)"
        ),
        "cancellation_frontier": (
            "Decay at scaled time c requires "
            "Phi(r)<r/2+c*r*(2-r)/8 on every untreated radius r"
        ),
        "known_threshold": (
            "The current published exponent-pair hull gives "
            f"c_*={C_STAR.numerator}/{C_STAR.denominator}="
            f"{float(C_STAR):.15f}..."
        ),
        "outer_strip": (
            "For every fixed c>2, sigma>1 and the Euler product supplies a "
            "zeta floor; therefore 2<c<=c_* is a cancellation gap, not a "
            "zero-free gap"
        ),
        "one_line": (
            "At c=2, sigma=1+o(1). Published zero-free estimates reach "
            "sigma>=1-1/(13*log(tau)); in that strip "
            "|1/zeta|<=4904*log(tau)^(11/12), while published 1-line "
            "logarithmic-derivative bounds are O(log(tau)/loglog(tau))"
        ),
        "conditional_c2": (
            "A cancellation envelope satisfying "
            "Phi(r)<=r-r^2/4-delta*r for some delta>0 transfers "
            "D_0,D_1,D_2 to zeta at c=2; the resulting M^(-delta) "
            "dyadic majorant handles growing blocks, while each fixed block "
            "vanishes because t=2/L; the "
            "published zero-free/1-line bounds then leave phase speed "
            "-L/4+o(L) and a polylogarithmic nonzero amplitude floor"
        ),
        "inner_wall": (
            "For fixed c<2, sigma=1/2+c/4<1. Existing zero-free regions "
            "approach 1 and do not supply a uniform zeta floor on that line; "
            "improved exponential-sum cancellation alone does not close the "
            "phase-critical-value target"
        ),
        "critical_deficit": (
            "At r_*=125662/155153, Phi_known=220633/310306 but the c=2 "
            "frontier is 15549101725/24072453409; the exact phase-exponent "
            "excess is 3133668399/48144906818"
        ),
        "symmetric_benchmark": (
            "At r_*, a symmetric pair (theta,1/2+theta) can meet the local "
            "c=2 frontier only if theta<=2900341791/24072453409="
            "0.120483846898448..., versus the published theta=13/84"
        ),
        "wronskian_handoff": (
            "The outer improvement would shrink the live thin-collar "
            "Wronskian problem from 0<c<=c_*+o(1) to fixed c<2 (plus a "
            "vanishing transition strip); it would not prove the inner phase theorem"
        ),
        "frontier": frontier,
    }


def build_artifact() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="np15czfw_01_scaled_geometry",
            role="exact_asymptotic_geometry",
            readiness="ready_to_apply",
            claim="The scaled Newman saddle maps c=2 exactly to the zeta 1-line at leading order.",
            formula=exact["scaled_geometry"],
            proof_boundary="Exact leading asymptotic coordinate map only.",
        ),
        GateRow(
            id="np15czfw_02_weighted_frontier",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="Every prospective cancellation estimate has an exact weighted decay frontier.",
            formula=exact["weighted_block"] + "; " + exact["cancellation_frontier"],
            proof_boundary="Abel-summation exponent bookkeeping; it does not supply a new exponential-sum estimate.",
        ),
        GateRow(
            id="np15czfw_03_known_threshold",
            role="validated_published_composition",
            readiness="ready_to_apply",
            claim="The known exponent-pair hull reproduces the exact current threshold c_*.",
            formula=exact["known_threshold"],
            proof_boundary="Uses the previously validated eleven-pair envelope.",
            diagnostics=exact["frontier"]["transition_rows"],
        ),
        GateRow(
            id="np15czfw_04_outer_strip",
            role="theorem_search_split",
            readiness="ready_to_apply",
            claim="The still-open fixed strip above c=2 is purely a cancellation deficit within this route.",
            formula=exact["outer_strip"],
            proof_boundary="Classifies the missing input; it does not close 2<c<=c_*.",
        ),
        GateRow(
            id="np15czfw_05_one_line_input",
            role="published_zero_free_input",
            readiness="available_published",
            claim="The c=2 line retains a quantitative nonzero zeta floor and subleading logarithmic phase jets.",
            formula=exact["one_line"],
            proof_boundary="These inputs become useful only after a stronger weighted cancellation theorem.",
        ),
        GateRow(
            id="np15czfw_06_conditional_c2",
            role="conditional_handoff",
            readiness="not_ready_to_apply",
            claim="An all-radius c=2 cancellation frontier with radius-proportional slack would close the zeta phase-amplitude handoff on the 1-line.",
            formula=exact["conditional_c2"],
            proof_boundary="Conditional on an exponential-sum estimate stronger than the current published hull.",
        ),
        GateRow(
            id="np15czfw_07_inner_wall",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Cancellation alone cannot be promoted through a fixed line Re(s_*)<1 using the present zeta-floor argument.",
            formula=exact["inner_wall"],
            proof_boundary="Route boundary only; it does not prove that another inner-strip method is impossible.",
        ),
        GateRow(
            id="np15czfw_08_critical_deficit",
            role="exact_frontier_target",
            readiness="ready_to_apply",
            claim="The worst known transition carries an explicit local phase-exponent deficit relative to c=2.",
            formula=exact["critical_deficit"],
            proof_boundary="One necessary local benchmark, not a sufficient global exponent-pair construction.",
        ),
        GateRow(
            id="np15czfw_09_symmetric_benchmark",
            role="theorem_search_benchmark",
            readiness="ready_to_apply",
            claim="A central symmetric exponent pair would need a substantial improvement even at the current worst radius.",
            formula=exact["symmetric_benchmark"],
            proof_boundary="Local necessary benchmark only; all radii still require coverage.",
        ),
        GateRow(
            id="np15czfw_10_wronskian_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The cancellation programme can at best hand the proof to an inner phase-critical Wronskian theorem.",
            formula=exact["wronskian_handoff"],
            proof_boundary="The inner theorem remains RH-level and open; this is not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate",
        "date": "2026-07-17",
        "status": (
            "exact cancellation/zero-free wall decomposition for the Newman "
            "scaled layer; conditional c=2 handoff and open inner phase target, "
            "not Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "This artifact proves exponent bookkeeping and classifies the known "
            "published inputs. It does not improve an exponential-sum estimate, "
            "close 2<c<=c_*, prove the conditional c=2 handoff, supply a zeta "
            "floor for fixed c<2, prove thin-collar Wronskian separation, prove "
            "Lambda<=0, or prove RH."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "sources": [
            POLYMATH_SOURCE_URL,
            EXPONENT_PAIR_SOURCE_URL,
            BOURGAIN_SOURCE_URL,
            ZERO_FREE_SOURCE_URL,
            ONE_LINE_SOURCE_URL,
            RECIPROCAL_SOURCE_URL,
            "outputs/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    frontier = exact["frontier"]
    lines = [
        "# Jensen-Window PF Newman Polymath-15 Cancellation/Zero-Free Wall Gate",
        "",
        "Date: 2026-07-17",
        "",
        "Status: exact theorem-search gate separating the cancellation and",
        "zero-free walls. This is not a proof of `Lambda <= 0` or RH.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.py",
        "```",
        "",
        "## Exact Frontier",
        "",
        "```text",
        exact["scaled_geometry"],
        exact["weighted_block"],
        exact["cancellation_frontier"],
        "```",
        "",
        "The eleven-pair published envelope gives",
        "",
        "```text",
        exact["known_threshold"],
        "```",
        "",
        "Its maximum occurs at the pair-2/pair-1 transition:",
        "",
        "```text",
        exact["critical_deficit"],
        exact["symmetric_benchmark"],
        "```",
        "",
        "The symmetric-pair number is only a necessary local benchmark. A clean",
        "sufficient condition for `c=2` is the complete phase curve",
        "`Phi(r)<=r-r^2/4-delta*r` for some `delta>0`. The radius-proportional",
        "margin becomes an `M^(-delta)` dyadic gain even as `r` tends to zero.",
        "",
        "## Frontier Table",
        "",
        "| active transition | r | required c | excess above c=2 curve |",
        "|---:|---:|---:|---:|",
    ]
    for row in frontier["transition_rows"]:
        if row["lower_index"] is None:
            label = "pair 0 endpoint"
        else:
            label = f"{row['upper_index']} -> {row['lower_index']}"
        lines.append(
            "| {label} | {r} | {c} | {gap} |".format(
                label=label,
                r=row["radius"]["decimal"],
                c=row["required_scaled_time"]["decimal"],
                gap=row["c2_phase_excess"]["decimal"],
            )
        )
    lines.extend(
        [
            "",
            "## Two Different Walls",
            "",
            "```text",
            exact["outer_strip"],
            "```",
            "",
            "Thus `2<c<=c_*` is a concrete analytic-number-theory target: improve",
            "the weighted log-phase cancellation while retaining the existing Euler",
            "product floor. It is not yet closed.",
            "",
            "At the endpoint `c=2`, the leading saddle lies on the zeta 1-line.",
            "The relevant published results are:",
            "",
            f"- [Mossinghoff-Trudgian-Yang zero-free region]({ZERO_FREE_SOURCE_URL})",
            f"- [Cully-Hugill-Leong 1-line estimates]({ONE_LINE_SOURCE_URL})",
            f"- [Leong reciprocal estimates]({RECIPROCAL_SOURCE_URL})",
            "",
            "```text",
            exact["one_line"],
            exact["conditional_c2"],
            "```",
            "",
            "The `o(L)` logarithmic phase correction cannot cancel the normalizer's",
            "`-L/4` phase speed, and the polylogarithmic zeta floor is still vastly",
            "larger than the exponentially small Polymath-15 remainder. This closes",
            "the phase-amplitude step only conditionally on the missing strict",
            "`c=2` cancellation frontier.",
            "",
            "For every fixed `c<2`, however, the saddle enters a fixed vertical line",
            "inside the critical strip:",
            "",
            "```text",
            exact["inner_wall"],
            "```",
            "",
            "This is a route boundary, not an impossibility theorem. A direct",
            "Wronskian, zero-dynamical, or other Xi-specific argument could still",
            "cross it; ordinary exponent-pair improvement by itself cannot provide",
            "the missing nonzero phase amplitude.",
            "",
            "## Proof Handoff",
            "",
            "```text",
            exact["wronskian_handoff"],
            "```",
            "",
            "So the next analytic-number-theory milestone is unambiguous: lower the",
            "known cancellation frontier from `c_*` toward `2`. The genuinely inner",
            "problem remains quantitative phase-critical-value avoidance for the",
            "corrected Riemann-Siegel main.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 cancellation/zero-free wall gate: "
        f"{len(artifact['rows'])} rows, "
        f"{len(artifact['exact']['frontier']['transition_rows'])} frontier points, "
        f"c_*={C_STAR.numerator}/{C_STAR.denominator}, 1 conditional c=2 handoff, "
        "1 open inner phase target"
    )


if __name__ == "__main__":
    main()
