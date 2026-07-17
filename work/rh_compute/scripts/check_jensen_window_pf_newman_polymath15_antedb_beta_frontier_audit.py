#!/usr/bin/env python3
"""Validate the current-ANTEDB beta-frontier audit."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit import (
    ALPHA_STAR,
    ANTEDB_COMMIT,
    BETA_STAR,
    C_STAR,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    R_STAR,
    beta_required_time,
    build_exact,
)


EXPECTED_IDS = [
    f"np15abf_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "coordinate_map"),
        (2, "source_snapshot"),
        (3, "post2023_pairs"),
        (4, "pair_hull"),
        (5, "beta_pipeline"),
        (6, "direct_contact"),
        (7, "c2_deficit"),
        (8, "stale_table_gate"),
        (9, "recursive_boundary"),
        (10, "handoff"),
    )
]


def validate(path: Path, note_path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_polymath15_antedb_beta_frontier_audit"
    ):
        issues.append("artifact kind mismatch")
    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")

    exact = build_exact()
    if artifact.get("exact") != exact:
        issues.append("exact payload drifted")
    source = exact["source_snapshot"]
    if source.get("antedb_commit") != ANTEDB_COMMIT:
        issues.append("ANTEDB source commit drifted")

    contact = exact["critical_contact"]
    if Fraction(contact["alpha_star"]["exact"]) != ALPHA_STAR:
        issues.append("alpha_star drifted")
    if Fraction(contact["radius_star"]["exact"]) != R_STAR:
        issues.append("radius_star drifted")
    if Fraction(contact["beta_star"]["exact"]) != BETA_STAR:
        issues.append("beta_star drifted")
    if beta_required_time(ALPHA_STAR, BETA_STAR) != C_STAR:
        issues.append("direct-beta c_star identity failed")

    hull = exact["pair_hull"]
    maximum = hull["maximum"]
    if maximum["required_scaled_time"]["exact"] != (
        "4911678521/1933561194"
    ):
        issues.append("current pair-hull maximum drifted")
    if maximum["active"] != ["TY1", "TY2"]:
        issues.append("current pair-hull active contact drifted")
    if Fraction(hull["tail_trivial_time_ceiling"]["exact"]) >= C_STAR:
        issues.append("Heath-Brown tail guard failed")

    post_rows = exact["post_2023_pairs"]
    if len(post_rows) != 6:
        issues.append("post-2023 pair count drifted")
    if any(
        Fraction(row["excess_over_beta_star"]["exact"]) < 0
        for row in post_rows
    ):
        issues.append("post-2023 pair lowered the recorded contact")

    pipeline = exact["antedb_pipeline"]
    if pipeline.get("final_piece_count_in_two_seeded_runs") != 68:
        issues.append("seeded final piece count drifted")
    if "not promoted" not in pipeline.get("reproducibility_boundary", ""):
        issues.append("optimizer tie-partition boundary missing")
    iteration = exact["beta_only_iteration_audit"]
    if iteration.get("passes_completed") != 12:
        issues.append("beta-only pass count drifted")
    if len(iteration.get("piece_counts", [])) != 13:
        issues.append("beta-only piece-count trace drifted")
    if iteration.get("contact_beta_on_every_pass", {}).get("exact") != (
        "220633/620612"
    ):
        issues.append("beta-only contact drifted")
    if iteration.get("global_threshold_on_every_pass", {}).get("exact") != (
        "4911678521/1933561194"
    ):
        issues.append("beta-only global threshold drifted")
    if "infinitely many" not in iteration.get("proof_boundary", ""):
        issues.append("finite-iteration proof boundary missing")

    by_id = {row.get("id"): row for row in artifact.get("rows", [])}
    if by_id.get("np15abf_08_stale_table_gate", {}).get("readiness") != (
        "guard_validated"
    ):
        issues.append("stale-table nonpromotion gate missing")
    if by_id.get("np15abf_10_handoff", {}).get("readiness") != (
        "not_ready_to_apply"
    ):
        issues.append("beta handoff was promoted")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove the external database complete",
        "recursive closure",
        "improve beta at the critical contact",
        "close c=2",
        "fixed c<2",
        "Wronskian separation",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = note_path.read_text(encoding="utf-8")
    for marker in (
        "ANTEDB Beta-Frontier Audit",
        ANTEDB_COMMIT,
        "alpha_*=62831/155153",
        "beta_*=220633/620612",
        "c_req(alpha_*)=4911678521/1933561194",
        "PYTHONHASHSEED=0",
        "twelve times",
        "contact beta on every pass=220633/620612",
        "global threshold on every pass=4911678521/1933561194",
        "broader exploratory beta/exponent-pair fixed-point run was terminated",
        "not a proof of cancellation",
        "This is not a proof of",
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
        "validated Newman Polymath-15 ANTEDB beta-frontier audit: "
        f"{len(artifact['rows'])} rows, 0 issues, 6 post-2023 pairs, "
        "1 exact current-hull maximum, 1 direct-beta contact, "
        "1 unchanged c=2 deficit, 1 nonpromotion gate"
    )


if __name__ == "__main__":
    main()
