#!/usr/bin/env python3
"""Validate the critical C1 endpoint peeling contract."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract import (
    DEFAULT_OUT,
    build_exact,
)


EXPECTED_IDS = [
    f"np15c1ep_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "published_expansion"),
        (2, "flint_coefficient_assembly"),
        (3, "first_recurrence"),
        (4, "first_coefficient"),
        (5, "critical_specialization"),
        (6, "centered_heat_cancellation"),
        (7, "corrected_endpoint_definition"),
        (8, "second_order_handoff"),
        (9, "signed_contact_handoff"),
        (10, "critical_cutoff_parity"),
        (11, "cutoff_guard"),
        (12, "uniform_quantitative_target"),
    )
]

EXPECTED_SUMMARY = {
    "rows": 12,
    "exact_coefficient_identities": 5,
    "primary_source_inputs": 2,
    "corrected_main_definitions": 1,
    "conditional_transfers": 2,
    "route_guards": 1,
    "open_quantitative_targets": 1,
}


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_polymath15_critical_RS_C1_"
        "endpoint_peeling_contract"
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
        "d_0^(1)=1/12",
        "F'''(p)/(12*pi^2)",
        "integral C_1",
        "r_[0]=-DeltaQ+r_[1]",
        "a^(-2)=2*pi/T_0",
        "C_1(-p,1/2)=-C_1(p,1/2)",
    ):
        if not any(marker in str(value) for value in exact.values()):
            issues.append(f"exact marker missing: {marker}")

    by_id = {row.get("id"): row for row in rows}
    if (
        by_id.get("np15c1ep_08_second_order_handoff", {}).get("readiness")
        != "conditional_ready"
    ):
        issues.append("second-order handoff was promoted")
    if (
        by_id.get("np15c1ep_11_cutoff_guard", {}).get("readiness")
        != "guard_validated"
    ):
        issues.append("cutoff guard missing")
    if (
        by_id.get("np15c1ep_12_uniform_quantitative_target", {}).get(
            "readiness"
        )
        != "not_ready_to_apply"
    ):
        issues.append("uniform quantitative target was promoted")

    sources = artifact.get("sources", [])
    for marker in (
        "arxiv.org/abs/1904.12438",
        "flintlib.org/doc/acb_dirichlet.html",
        "zeta_rs_d_coeffs.c",
        "zeta_rs_f_coeffs.c",
        "zeta_rs_r.c",
    ):
        if not any(marker in source for source in sources):
            issues.append(f"primary-source marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "uniform second-order value/derivative bounds",
        "adjacent-cutoff splice",
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
        "validated Newman Polymath-15 critical C1 endpoint peeling contract: "
        f"{len(artifact['rows'])} rows, 5 exact coefficient identities, "
        "2 primary-source inputs, 1 corrected-main definition, "
        "2 conditional transfers, 1 route guard, "
        "1 open quantitative target"
    )


if __name__ == "__main__":
    main()
