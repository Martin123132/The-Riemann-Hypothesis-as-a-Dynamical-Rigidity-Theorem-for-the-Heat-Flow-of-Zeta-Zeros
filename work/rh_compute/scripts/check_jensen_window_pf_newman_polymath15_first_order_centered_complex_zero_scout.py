#!/usr/bin/env python3
"""Validate the selected corrected-complex-main zero scout."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

import jensen_window_pf_newman_polymath15_first_order_centered_complex_zero_scout as builder


EXPECTED_IDS = [
    "nfocczs_01_selected_lift",
    "nfocczs_02_complex_zero",
    "nfocczs_03_transverse_complex_zero",
    "nfocczs_04_simple_real_crossing",
    "nfocczs_05_wronskian_degeneracy",
    "nfocczs_06_centered_core",
    "nfocczs_07_precision",
    "nfocczs_08_interval_handoff",
]


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != builder.STEM:
        issues.append("artifact kind mismatch")
    if artifact.get("source_sha256") != builder.source_hashes():
        issues.append("source hash chain drifted")
    if artifact.get("configuration") != {
        "N_fixed": 6,
        "low_dps": 70,
        "high_dps": 95,
        "initial_pair_1": ["519.38", "0.403"],
        "initial_pair_2": ["519.40", "0.405"],
    }:
        issues.append("configuration drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if rows[-1].get("readiness") != "not_ready_to_apply":
        issues.append("interval handoff was promoted")

    selected = artifact.get("selected_zero", {})
    mp.mp.dps = 100
    try:
        if int(selected.get("N", -1)) != 6:
            issues.append("selected cutoff drifted")
        if not mp.mpf("0") < mp.mpf(selected["t"]) < mp.mpf("0.5"):
            issues.append("selected time left 0<t<1/2")
        if not mp.mpf(selected["L"]) < 50:
            issues.append("moderate-height proof boundary drifted")
        if not mp.mpf(selected["q"]) > 1:
            issues.append("selected point left q>1")
        if not mp.mpf(selected["E_abs"]) < mp.mpf("1e-70"):
            issues.append("complex-main residual exceeded 1e-70")
        if not abs(mp.mpf(selected["E_x_real_U"])) > mp.mpf("0.5"):
            issues.append("real crossing slope lost its margin")
        if not abs(mp.mpf(selected["complex_zero_jacobian_xt"])) > mp.mpf("0.3"):
            issues.append("two-variable Jacobian lost its margin")
        if not abs(mp.mpf(selected["wronskian"])) < mp.mpf("1e-70"):
            issues.append("Wronskian residual exceeded 1e-70")
        if not abs(
            mp.mpf(selected["centered_scalar_S"])
            - mp.mpf(selected["E_x_real_U"])
        ) < mp.mpf("1e-70"):
            issues.append("centered scalar no longer matches U at E=0")
        if not abs(
            mp.mpf(selected["E_x_real_U"])
            - mp.mpf(selected["core_scalar_A"])
            - mp.mpf(selected["real_D1x"])
        ) < mp.mpf("1e-65"):
            issues.append("core plus D1x reconstruction failed")
    except (KeyError, ValueError) as exc:
        issues.append(f"selected numeric payload malformed: {exc}")

    replay = builder.solve_selected_zero(105)
    for key in ("x", "t", "E_x_real_U", "complex_zero_jacobian_xt"):
        if abs(mp.mpf(replay[key]) - mp.mpf(selected[key])) >= mp.mpf("1e-50"):
            issues.append(f"independent high-precision replay drifted: {key}")

    proof_boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove",
        "interval arithmetic",
        "L>=50",
        "q>=1",
        "contact exclusion",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in proof_boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Selected Point",
        "Precision Replay",
        "Route Consequence",
        "not an interval certificate",
        "crossing slope is about `0.556`",
        "optional next",
        "not a proof of RH",
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
        "validated Newman first-order centered complex-zero scout: "
        "N=6, |E_[1]|<1e-70, |U|>0.5, "
        "|det D_(x,t)E_[1]|>0.3, cross-precision drift <1e-50"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
