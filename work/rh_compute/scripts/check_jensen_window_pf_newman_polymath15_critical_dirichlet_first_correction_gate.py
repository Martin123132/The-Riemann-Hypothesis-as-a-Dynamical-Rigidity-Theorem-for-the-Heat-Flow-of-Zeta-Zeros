#!/usr/bin/env python3
"""Validate the critical Dirichlet first-correction gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_critical_dirichlet_first_correction_gate import (
    DEFAULT_OUT,
    build_exact,
)


EXPECTED_IDS = [
    f"np15cdfc_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "heat_shift"),
        (2, "scaled_gamma"),
        (3, "exact_relative_ratio"),
        (4, "gaussian_moment"),
        (5, "first_correction"),
        (6, "remainder_decomposition"),
        (7, "gamma_guard"),
        (8, "central_tail_guard"),
        (9, "signed_main_update"),
        (10, "uniform_target"),
    )
]

EXPECTED_SUMMARY = {
    "rows": 10,
    "exact_identities": 4,
    "primary_source_inputs": 2,
    "analytic_guards": 2,
    "signed_main_definitions": 1,
    "open_uniform_targets": 1,
}


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_polymath15_critical_"
        "dirichlet_first_correction_gate"
    ):
        issues.append("artifact kind mismatch")

    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")
    if artifact.get("summary") != EXPECTED_SUMMARY:
        issues.append(f"summary mismatch: {artifact.get('summary')}")

    exact = build_exact()
    if artifact.get("exact") != exact:
        issues.append("saved exact dictionary differs from regenerated identities")
    for marker in (
        "t/2+t^2*alpha_n^2/4",
        "1/(6s)+alpha'(s)",
        "C_D/T^2",
        "exp(-7L/4)",
        "r_[0]=DeltaJ+r_[1]",
        "exp(-5L/4)",
    ):
        if not any(marker in str(value) for value in exact.values()):
            issues.append(f"exact marker missing: {marker}")

    by_id = {row.get("id"): row for row in rows}
    for row_id in (
        "np15cdfc_07_gamma_guard",
        "np15cdfc_08_central_tail_guard",
    ):
        if by_id.get(row_id, {}).get("readiness") != "conditional_ready":
            issues.append(f"analytic guard was promoted: {row_id}")
    if (
        by_id.get("np15cdfc_10_uniform_target", {}).get("readiness")
        != "not_ready_to_apply"
    ):
        issues.append("uniform target was promoted")

    sources = artifact.get("sources", [])
    for marker in (
        "arxiv.org/abs/1904.12438",
        "dlmf.nist.gov/5.11",
        "critical_RS_C1_endpoint_peeling_contract",
        "contact_normal_hierarchy_gate",
    ):
        if not any(marker in source for source in sources):
            issues.append(f"source marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "uniform C_D/T^2 remainder",
        "exp(-5L/4) C^1 constants",
        "cutoff collar",
        "strict signed contact inequality",
        "contact exclusion",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman Polymath-15 critical Dirichlet first correction gate: "
        f"{len(artifact['rows'])} rows, 4 exact identities, "
        "2 primary-source inputs, 2 analytic guards, "
        "1 signed-main definition, 1 open uniform target"
    )


if __name__ == "__main__":
    main()
