#!/usr/bin/env python3
"""Validate the Newman Gaussian/Legendre duality gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate import (
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_exact,
)


EXPECTED_IDS = [
    f"np15gld_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "gaussian_identity"),
        (2, "moment_identity"),
        (3, "scaled_cost"),
        (4, "partial_sum_profile"),
        (5, "saddle_identity"),
        (6, "dual_equivalence"),
        (7, "cstar_contact"),
        (8, "c2_deficit"),
        (9, "nonpromotion"),
        (10, "handoff"),
    )
]


def validate(path: Path, note_path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_polymath15_gaussian_legendre_duality_gate"
    ):
        issues.append("artifact kind mismatch")
    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")
    exact = build_exact()
    if artifact.get("exact") != exact:
        issues.append("exact payload drifted")

    star = exact["critical_star"]
    if star["line"]["exact"] != "4800718975/7734244776":
        issues.append("critical dual line drifted")
    if star["gaussian_shift"]["exact"] != "31657/61548":
        issues.append("critical Gaussian shift drifted")
    if star["partial_sum_exponent"]["exact"] != (
        "1989040967/9549356844"
    ):
        issues.append("critical partial-sum exponent drifted")
    if star["partial_sum_exponent"] != star["gaussian_cost"]:
        issues.append("c_star equality contact failed")
    if Fraction(star["net_exponent"]["exact"]) != 0:
        issues.append("c_star net exponent is not zero")

    c_two = exact["critical_c2"]
    if c_two["line"]["exact"] != "92322/155153":
        issues.append("c=2 dual line drifted")
    if c_two["net_exponent"]["exact"] != (
        "3133668399/48144906818"
    ):
        issues.append("c=2 exact deficit drifted")

    by_id = {row.get("id"): row for row in artifact.get("rows", [])}
    if by_id.get("np15gld_09_nonpromotion", {}).get("readiness") != (
        "guard_validated"
    ):
        issues.append("Gaussian nonpromotion gate missing")
    if by_id.get("np15gld_10_handoff", {}).get("readiness") != (
        "not_ready_to_apply"
    ):
        issues.append("analytic handoff was promoted")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not improve a partial-sum estimate",
        "strict decay at c_*",
        "fixed c<2",
        "Wronskian separation",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = note_path.read_text(encoding="utf-8")
    for marker in (
        "Gaussian/Legendre Duality Gate",
        "q_*=4800718975/7734244776",
        "E(q_*)=1989040967/9549356844",
        "deficit=3133668399/48144906818",
        "not an artifact",
        "does not move the",
        "This is not a proof of `Lambda <= 0` or RH",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
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
        "validated Newman Polymath-15 Gaussian/Legendre duality gate: "
        f"{len(artifact['rows'])} rows, 0 issues, 1 exact Gaussian identity, "
        "1 exact Legendre equivalence, 1 c_* equality point, "
        "1 c=2 deficit, 1 nonpromotion gate"
    )


if __name__ == "__main__":
    main()
