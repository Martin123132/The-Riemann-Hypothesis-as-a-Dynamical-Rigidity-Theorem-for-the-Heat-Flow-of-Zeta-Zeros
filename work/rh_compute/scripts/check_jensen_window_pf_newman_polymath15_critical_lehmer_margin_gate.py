#!/usr/bin/env python3
"""Validate the critical Lehmer-margin stress gate independently."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import mpmath as mp

import jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.md"
GENERATOR = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.py"
)
EXPECTED_IDS = {
    "np15clmg_01_exact_pair_identity",
    "np15clmg_02_heat_transversality",
    "np15clmg_03_quadratic_model",
    "np15clmg_04_arb_exact_point",
    "np15clmg_05_fixed_floor_rejected",
    "np15clmg_06_corrected_scout",
    "np15clmg_07_target_refinement",
    "np15clmg_08_nonpromotion",
    "np15clmg_09_monotonicity_counterexample",
}


def load_artifact(rebuild: bool) -> dict:
    if rebuild:
        subprocess.run([sys.executable, str(GENERATOR)], cwd=REPO_ROOT, check=True)
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def validate(artifact: dict) -> list[str]:
    issues: list[str] = []
    if artifact.get("kind") != "jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate":
        issues.append("kind drifted")
    rows = artifact.get("rows", [])
    if {row.get("id") for row in rows} != EXPECTED_IDS:
        issues.append("gate row ids drifted")
    if len(rows) != 9:
        issues.append("expected nine gate rows")

    exact = gate.build_exact()
    if artifact.get("exact") != exact:
        issues.append("exact identity payload drifted")
    recomputed_arb = gate.arb_exact_point()
    if artifact.get("arb") != recomputed_arb:
        issues.append("Arb exact-point payload drifted")
    if recomputed_arb.get("certified_bounds") != "0<C[H_0/A_0]/L^2<1/1000<3/40":
        issues.append("small-margin certificate drifted")
    if recomputed_arb.get("certified_monotonicity_bounds") != (
        "M_0(x)=H_0(x)H_0'''(x)-H_0'(x)H_0''(x)=-L_0'(x)<0"
    ):
        issues.append("strict-monotonicity counter-certificate drifted")

    scout = artifact.get("scout", {})
    point_rows = scout.get("rows", [])
    if [row.get("c") for row in point_rows] != list(gate.C_VALUES):
        issues.append("scaled-time grid drifted")
    if len(point_rows) != 4:
        issues.append("expected four Lehmer diagnostic rows")
    for row in point_rows:
        try:
            curvature = mp.mpf(row["corrected_curvature"])
            scaled = mp.mpf(row["corrected_curvature_over_L2"])
            transversality = mp.mpf(row["transversality_norm"])
        except (KeyError, ValueError):
            issues.append("malformed Lehmer diagnostic row")
            continue
        if curvature <= 0 or scaled <= 0 or transversality <= 0:
            issues.append(f"nonpositive Lehmer diagnostic at c={row.get('c')}")
    try:
        precision_delta = mp.mpf(scout["max_relative_corrected_curvature_delta"])
        if precision_delta >= mp.mpf("1e-25"):
            issues.append("coarse/fine precision threshold failed")
    except (KeyError, ValueError):
        issues.append("missing precision diagnostic")

    note = NOTE.read_text(encoding="utf-8")
    required = (
        "0<C[H_0/A_0]/L^2<1/1000<3/40",
        "first-jet norm",
        "Strict Monotonicity Counter-Gate",
        "M_0(x)=H_0(x)H_0'''(x)-H_0'(x)H_0''(x)=-L_0'(x)<0",
        "dominated convergence",
        "not a proof of `Lambda <= 0` or RH",
    )
    for marker in required:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    if "not a proof" not in artifact.get("status", ""):
        issues.append("status nonpromotion language missing")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()
    artifact = load_artifact(args.rebuild)
    issues = validate(artifact)
    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        raise SystemExit(1)
    print(
        "validated Newman Polymath-15 critical Lehmer margin gate: "
        "9 rows, 0 issues, 1 Arb small-margin point, 3 exact local identities, "
        "4 corrected diagnostics, 1 fixed-floor rejection, "
        "1 Arb strict-monotonicity rejection"
    )


if __name__ == "__main__":
    main()
