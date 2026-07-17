#!/usr/bin/env python3
"""Validate the corrected critical transversality target."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import jensen_window_pf_newman_polymath15_critical_transversality_target as target


REPO_ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_transversality_target.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_polymath15_critical_transversality_target.md"
GENERATOR = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_transversality_target.py"
)
EXPECTED_IDS = {
    "np15ctt_01_normalized_double_zero",
    "np15ctt_02_corrected_split",
    "np15ctt_03_first_jet_identity",
    "np15ctt_04_sum_of_squares",
    "np15ctt_05_c1_collision_budget",
    "np15ctt_06_derivative_saving",
    "np15ctt_07_live_target",
    "np15ctt_08_global_partition",
    "np15ctt_09_conditional_endgame",
    "np15ctt_10_nonpromotion",
}


def load_artifact(rebuild: bool) -> dict:
    if rebuild:
        subprocess.run([sys.executable, str(GENERATOR)], cwd=REPO_ROOT, check=True)
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def validate(artifact: dict) -> list[str]:
    issues: list[str] = []
    if artifact.get("kind") != "jensen_window_pf_newman_polymath15_critical_transversality_target":
        issues.append("kind drifted")
    rows = artifact.get("rows", [])
    if len(rows) != 10:
        issues.append("expected ten rows")
    if {row.get("id") for row in rows} != EXPECTED_IDS:
        issues.append("row ids drifted")
    exact = target.build_exact()
    if artifact.get("exact") != exact:
        issues.append("exact payload drifted")
    if sum(row.get("role") == "open_theorem_target" for row in rows) != 1:
        issues.append("expected one open theorem target")
    if sum(row.get("role") == "exact_exclusion_lemma" for row in rows) != 1:
        issues.append("expected one exact exclusion lemma")
    live = next(
        (row for row in rows if row.get("id") == "np15ctt_07_live_target"),
        {},
    )
    if live.get("readiness") != "not_ready_to_apply":
        issues.append("live target was promoted")
    if "epsilon_0^2+epsilon_1^2" not in exact.get("strict_target", ""):
        issues.append("strict C1 budget drifted")
    if "no second-derivative" not in exact.get("cauchy_gain", ""):
        issues.append("derivative-saving statement drifted")
    if "Lambda<=0" not in exact.get("conditional_endgame", ""):
        issues.append("conditional endgame drifted")
    if "c_*=4911678521/1933561194" not in exact.get(
        "global_composition", ""
    ):
        issues.append("current asymptotic partition drifted")
    if "bounded-L remainder" not in exact.get("global_composition", ""):
        issues.append("existential-threshold compact handoff missing")

    note = NOTE.read_text(encoding="utf-8")
    required = (
        "H_t=H_t'=0 if and only if Z_t=Z_t'=0",
        "T_L[J]",
        "C1 Exclusion Lemma",
        "c_*=4911678521/1933561194",
        "remains open",
        "not a proof of `Lambda <= 0` or RH",
    )
    for marker in required:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    status = artifact.get("status", "")
    if "not a proof" not in status or "RH" not in status:
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
        "validated Newman Polymath-15 critical transversality target: "
        "10 rows, 0 issues, 3 exact identities, 1 C1 exclusion lemma, "
        "1 derivative saving, 1 open first-jet target"
    )


if __name__ == "__main__":
    main()
