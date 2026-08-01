#!/usr/bin/env python3
"""Validate the first-order oriented successor-winding reduction."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_oriented_successor_winding_reduction as builder


HISTORICAL_RAY_SHA256 = (
    "2f2c7e69957446d461354e673ff27af66c346863e9264cf2d573be0b470c10e5"
)


EXPECTED_IDS = [
    "nfosw_01_orientation",
    "nfosw_02_scaled_pruefer_form",
    "nfosw_03_signed_ray_count",
    "nfosw_04_horizontal_crossings",
    "nfosw_05_vertical_crossings",
    "nfosw_06_successor_chain",
    "nfosw_07_base_degree",
    "nfosw_08_outer_strip_degree",
    "nfosw_09_successor_recurrence",
    "nfosw_10_one_sided_integer_trap",
    "nfosw_11_proxy_transfer",
    "nfosw_12_frequency_margin",
    "nfosw_13_parabolic_and_shoulder_margin",
    "nfosw_14_one_sided_phase_target",
]


def independent_phase_audit(issues: list[str]) -> None:
    value, value_x, value_xx = sp.symbols("J J_x J_xx", real=True)
    value_t, value_xt = sp.symbols("J_t J_xt", real=True)
    scale = sp.symbols("ell", positive=True)
    scale_x, scale_t = sp.symbols("ell_x ell_t", real=True)
    norm_sq = value**2 + (value_x / scale) ** 2
    phase_x = sp.factor(
        (
            value
            * (value_xx / scale - value_x * scale_x / scale**2)
            - value_x**2 / scale
        )
        / norm_sq
    )
    phase_t = sp.factor(
        (
            value
            * (value_xt / scale - value_x * scale_t / scale**2)
            - value_x * value_t / scale
        )
        / norm_sq
    )
    expected_x = (
        scale * value * value_xx
        - value * value_x * scale_x
        - scale * value_x**2
    ) / (scale**2 * norm_sq)
    expected_t = (
        scale * value * value_xt
        - value * value_x * scale_t
        - scale * value_x * value_t
    ) / (scale**2 * norm_sq)
    if sp.simplify(phase_x - expected_x) != 0:
        issues.append("independent horizontal phase identity failed")
    if sp.simplify(phase_t - expected_t) != 0:
        issues.append("independent vertical phase identity failed")
    if sp.simplify(phase_x.subs(value, 0) + scale) != 0:
        issues.append("independent simple-zero phase sign failed")


def independent_chain_audit(issues: list[str]) -> None:
    p_j = Counter(
        bottom_old=1,
        right_upper=1,
        top_old=1,
        axis_upper=1,
    )
    descendant = Counter(
        bottom_new_inner=1,
        right_lower=1,
        bottom_old=-1,
        axis_lower=1,
    )
    outer = Counter(
        bottom_new_outer=1,
        right_next=1,
        top_outer=1,
        right_lower=-1,
        right_upper=-1,
    )
    expected = Counter(
        bottom_new_inner=1,
        bottom_new_outer=1,
        right_next=1,
        top_outer=1,
        top_old=1,
        axis_upper=1,
        axis_lower=1,
    )
    if p_j + descendant + outer != expected:
        issues.append("independent successor chain cancellation failed")


def independent_countermodel_audit(issues: list[str]) -> None:
    t, x = sp.symbols("t x", real=True)
    flow = x**2 - 2 * t
    if sp.diff(flow, t) + sp.diff(flow, x, 2) != 0:
        issues.append("quadratic countermodel heat equation failed")
    jacobian = sp.Matrix(
        [
            [sp.diff(flow, x), sp.diff(flow, t)],
            [sp.diff(flow, x, 2), sp.diff(flow, x, t)],
        ]
    ).det().subs({t: 0, x: 0})
    if jacobian != 4:
        issues.append("quadratic countermodel standard degree sign failed")
    # F+iF_x can vanish only at x=0 and t=0, which is interior.
    if sp.solve(
        [flow, sp.diff(flow, x)],
        [t, x],
        dict=True,
    ) != [{t: 0, x: 0}]:
        issues.append("quadratic countermodel common-zero audit failed")


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != builder.STEM:
        issues.append("artifact kind mismatch")
    stored_sources = artifact.get("source_sha256", {})
    current_sources = builder.source_hashes()
    normalized_sources = dict(current_sources)
    stored_ray = stored_sources.get("ray_schedule")
    if stored_ray not in {
        current_sources.get("ray_schedule"),
        HISTORICAL_RAY_SHA256,
    }:
        issues.append("unrecognized ray source snapshot")
    normalized_sources["ray_schedule"] = stored_ray
    if stored_sources != normalized_sources:
        issues.append("source hash chain drifted")
    if artifact.get("source_audit") != builder.source_audit():
        issues.append("source audit drifted")
    if artifact.get("exact") != builder.build_exact():
        issues.append("exact payload drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    open_rows = [
        row for row in rows if row.get("role") == "open_theorem_target"
    ]
    if len(open_rows) != 3:
        issues.append(f"open target count drifted: {len(open_rows)}")
    if any(row.get("readiness") != "not_ready_to_apply" for row in open_rows):
        issues.append("an open target was promoted")

    independent_phase_audit(issues)
    independent_chain_audit(issues)
    independent_countermodel_audit(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "ind_(x,t)(H,H_x/ell)=+floor(m/2)>0",
        "partial_x arg(W_J)",
        "sum_(u=0,v>0) sign(-d_s u)",
        "partial P_(j+1)=partial P_j+partial D_j+partial S_j",
        "wind(V_H(partial S_j),0)=0",
        "kappa_j=sum_(p in D_j) floor(m_p/2)",
        "wind(proxy_j)<1",
        "13/500<1",
        "q=2t_jL^2>=1",
        "q<1",
        "Boundary nonvanishing alone does not imply",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    proof_boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove",
        "bottom margin",
        "finite shoulders",
        "one-sided phase bound",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in proof_boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Orientation",
        "Scaled Pruefer Form",
        "Signed Crossings",
        "Successor Chain",
        "One-Sided Integer Trap",
        "First-Order Transfer",
        "Remaining Xi Input",
        "Scope Guard",
        "Live Handoff",
        "not by itself make its open-path phase contribution zero",
        "Computing exact zero winding is unnecessary",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=builder.DEFAULT_OUT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Newman oriented successor-winding reduction: "
        "14 rows, 1 exact chain audit, 2 exact crossing reductions, "
        "1 one-sided integer trap, 3 open cofinal obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
