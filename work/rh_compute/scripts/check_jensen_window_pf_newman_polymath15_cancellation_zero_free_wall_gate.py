#!/usr/bin/env python3
"""Validate the Newman cancellation/zero-free wall gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate import (
    C_STAR,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    R_STAR,
    build_exact,
)


EXPECTED_IDS = [
    f"np15czfw_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "scaled_geometry"),
        (2, "weighted_frontier"),
        (3, "known_threshold"),
        (4, "outer_strip"),
        (5, "one_line_input"),
        (6, "conditional_c2"),
        (7, "inner_wall"),
        (8, "critical_deficit"),
        (9, "symmetric_benchmark"),
        (10, "wronskian_handoff"),
    )
]


def validate(path: Path, note_path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate"
    ):
        issues.append("artifact kind mismatch")
    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")

    exact = build_exact()
    frontier = exact["frontier"]
    rows = frontier["transition_rows"]
    if len(rows) != 11:
        issues.append("frontier point count drifted")
    required = [Fraction(row["required_scaled_time"]["exact"]) for row in rows]
    if max(required) != C_STAR:
        issues.append("known threshold is not the frontier maximum")
    maxima = [row for row in rows if row["is_current_maximum"]]
    if len(maxima) != 1:
        issues.append("frontier maximum is not unique")
    elif Fraction(maxima[0]["radius"]["exact"]) != R_STAR:
        issues.append("frontier maximum radius drifted")
    if any(Fraction(row["c2_phase_excess"]["exact"]) <= 0 for row in rows):
        issues.append("known frontier was promoted to c=2")
    if frontier["critical_c2_phase_deficit"]["exact"] != (
        "3133668399/48144906818"
    ):
        issues.append("critical c=2 deficit drifted")
    if frontier["critical_symmetric_theta_cap"]["exact"] != (
        "2900341791/24072453409"
    ):
        issues.append("symmetric theta benchmark drifted")

    by_id = {row.get("id"): row for row in artifact.get("rows", [])}
    if by_id.get("np15czfw_06_conditional_c2", {}).get("readiness") != (
        "not_ready_to_apply"
    ):
        issues.append("conditional c=2 handoff was promoted")
    if by_id.get("np15czfw_07_inner_wall", {}).get("readiness") != (
        "guard_validated"
    ):
        issues.append("inner-wall nonpromotion gate missing")
    if by_id.get("np15czfw_10_wronskian_handoff", {}).get("readiness") != (
        "not_ready_to_apply"
    ):
        issues.append("inner Wronskian target was promoted")

    sources = set(artifact.get("sources", []))
    for source in (
        "https://arxiv.org/abs/1904.12438",
        "https://arxiv.org/abs/2306.05599",
        "https://arxiv.org/abs/2212.06867",
        "https://arxiv.org/abs/2312.09412",
        "https://arxiv.org/abs/2405.04869",
    ):
        if source not in sources:
            issues.append(f"source missing: {source}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not improve an exponential-sum estimate",
        "fixed c<2",
        "Wronskian separation",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = note_path.read_text(encoding="utf-8")
    for marker in (
        "Cancellation/Zero-Free Wall Gate",
        "c_*=4911678521/1933561194",
        "Phi(r)<=r-r^2/4-delta*r",
        "2<c<=c_*",
        "fixed `c<2`",
        "not an impossibility theorem",
        "This is not a proof of `Lambda <= 0` or RH",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    issues = validate(args.artifact, args.note)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman Polymath-15 cancellation/zero-free wall gate: "
        f"{len(artifact['rows'])} rows, 0 issues, 11 frontier points, "
        "1 exact c_* maximum, 1 conditional c=2 handoff, "
        "1 inner-wall nonpromotion gate, 1 open Wronskian handoff"
    )


if __name__ == "__main__":
    main()
