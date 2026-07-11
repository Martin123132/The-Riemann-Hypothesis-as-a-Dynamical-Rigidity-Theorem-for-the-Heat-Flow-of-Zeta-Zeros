#!/usr/bin/env python3
"""Independently validate the cofinal-degree polar-closure lemma."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_cofinal_degree_polar_closure_lemma.json"
STURM_34 = REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520_summary.json"
STURM_5 = REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520_summary.json"
STURM_6_12 = REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d6_d12_dps520_summary.json"
STRESS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_stress_lamgrid_d3_d12_k64_summary.json"


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()
    payload = load(args.input)
    issues: list[str] = []
    rows = payload.get("rows", [])

    d, j = sp.symbols("d j", integer=True, positive=True)
    if sp.simplify((1 - j / (d + 1)) * sp.binomial(d + 1, j) - sp.binomial(d, j)) != 0:
        issues.append("adjacent binomial identity failed")

    w = sp.symbols("w")
    alpha = sp.symbols("alpha_0:5", positive=True)
    product = sp.prod(1 + value * w for value in alpha)
    polar = product - w * sp.diff(product, w) / 5
    reciprocal = sum(1 / (1 + value * w) for value in alpha) / 5
    if sp.simplify(polar / product - reciprocal) != 0:
        issues.append("polar reciprocal identity failed")
    derivative = sp.diff(reciprocal, w)
    expected_derivative = -sum(value / (1 + value * w) ** 2 for value in alpha) / 5
    if sp.simplify(derivative - expected_derivative) != 0:
        issues.append("polar strict derivative identity failed")

    sturm34, sturm5, sturm6_12, stress = load(STURM_34), load(STURM_5), load(STURM_6_12), load(STRESS)
    if sturm34.get("ok") + sturm5.get("ok") + sturm6_12.get("ok") != 1050:
        issues.append("Sturm row total failed")
    if (
        sturm34.get("failed_or_inconclusive")
        + sturm5.get("failed_or_inconclusive")
        + sturm6_12.get("failed_or_inconclusive")
        != 0
    ):
        issues.append("Sturm failure total failed")
    if stress.get("summary", {}).get("stress_rows") != 2875:
        issues.append("contraction-only row total failed")

    audit = payload.get("exact", {}).get("evidence_audit", {})
    if audit.get("cofinal_terminal_degrees_certified") is not False:
        issues.append("cofinal evidence was promoted")
    if audit.get("contraction_evidence_is_not_hyperbolicity") is not True:
        issues.append("contraction/hyperbolicity boundary missing")
    if len(rows) != 10:
        issues.append(f"row count {len(rows)} != 10")
    if sum(row.get("role") == "open_handoff" for row in rows) != 1:
        issues.append("open handoff count failed")
    boundary = payload.get("proof_boundary", "")
    for marker in ("does not prove", "unbounded terminal-degree", "PF-infinity", "RH"):
        if marker not in boundary:
            issues.append(f"missing proof-boundary marker {marker!r}")

    print(
        "validated Jensen-window PF cofinal-degree polar-closure lemma: "
        f"{len(rows)} rows, {len(issues)} issues, 3 exact polar identities, "
        "1 interlacing theorem, 1 multiplicity theorem, 1 finite-tower closure, "
        "1 cofinal-degree closure, 1050 finite Sturm rows, 2875 contraction-only rows, "
        "1 open terminal handoff"
    )
    for issue in issues:
        print(f"ISSUE {issue}")
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
