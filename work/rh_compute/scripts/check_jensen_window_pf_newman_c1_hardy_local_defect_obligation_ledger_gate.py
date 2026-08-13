#!/usr/bin/env python3
"""Validate the Hardy MGS local-defect obligation ledger and bypass handoff."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from jensen_window_pf_newman_c1_hardy_local_defect_obligation_ledger_gate import (
    ARXIV_ID,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    KIND,
    SOURCE,
    SOURCE_MARKERS,
    build_artifact,
    error_components,
    exact_accumulation_fixture,
    file_hash,
    source_marker_locations,
)


EXPECTED_IDS = [
    f"held_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "source_pin"),
        (2, "exact_until_69"),
        (3, "w1_anatomy"),
        (4, "fortran_mapping"),
        (5, "cutoff_conflict"),
        (6, "special_functions"),
        (7, "scheme_kernel"),
        (8, "exact_defect"),
        (9, "accumulation"),
        (10, "exact_shell"),
        (11, "termwise_route"),
        (12, "real_C6_route"),
        (13, "disk_route"),
        (14, "representation_split"),
        (15, "next_experiment"),
        (16, "physical_block"),
    )
]


def rational_complex(record: dict) -> tuple[Fraction, Fraction]:
    return Fraction(record["real"]), Fraction(record["imag"])


def validate(path: Path, note_path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != KIND:
        issues.append("artifact kind mismatch")
    if artifact.get("status") != "exact_defect_reduction_with_componentwise_open_bounds":
        issues.append("artifact status mismatch")
    if artifact != build_artifact():
        issues.append("artifact differs from a deterministic rebuild")

    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")

    paper = artifact.get("paper_audit", {})
    if paper.get("arxiv_id") != ARXIV_ID:
        issues.append("paper version drifted")
    locations = paper.get("audited_locations", [])
    if len(locations) != 6:
        issues.append("paper location count drifted")
    equations = {
        equation
        for location in locations
        for equation in location.get("equations", [])
    }
    for equation in ("69", "81", "96", "120", "121", "122", "123"):
        if equation not in equations:
            issues.append(f"paper equation missing: {equation}")

    source = artifact.get("source_audit", {})
    if source.get("sha256") != file_hash(SOURCE):
        issues.append("source hash drifted")
    rebuilt_locations = source_marker_locations()
    if source.get("marker_lines") != rebuilt_locations:
        issues.append("source marker locations drifted")
    if set(source.get("marker_lines", {})) != set(SOURCE_MARKERS):
        issues.append("source marker roster drifted")
    mappings = source.get("paper_to_fortran_map", [])
    expected_pairs = [
        ("W1", "t5"),
        ("W2", "t2"),
        ("W3_or_W4", "t4"),
        ("W5", "t1"),
        ("endpoint_half_sum", "t3"),
        ("CW_assembly_and_recurrence", "qq_and_csum"),
    ]
    actual_pairs = [
        (row.get("paper_component"), row.get("fortran_component"))
        for row in mappings
    ]
    if actual_pairs != expected_pairs:
        issues.append(f"paper/source mapping drifted: {actual_pairs}")
    if any(not row.get("source_lines") for row in mappings):
        issues.append("paper/source mapping lacks source lines")

    components = artifact.get("error_components", [])
    if components != error_components():
        issues.append("error-component ledger drifted")
    if len(components) != 13 or len({row.get("id") for row in components}) != 13:
        issues.append("error-component count or uniqueness drifted")
    required_scopes = {
        "local_MGS",
        "local_MGS_numeric",
        "local_MGS_model",
        "whole_Hardy_evaluator",
        "whole_shifted_batch",
    }
    if not required_scopes.issubset({row.get("scope") for row in components}):
        issues.append("error-component scopes are incomplete")
    for component in components:
        for key in ("origin", "remainder", "readiness", "certificate"):
            if not component.get(key):
                issues.append(f"component {component.get('id')} lacks {key}")

    routes = artifact.get("exact_defect_routes", {})
    fixture = routes.get("accumulation_fixture", {})
    if fixture != exact_accumulation_fixture():
        issues.append("exact accumulation fixture drifted")
    if rational_complex(fixture.get("recursive_root_error", {})) != rational_complex(
        fixture.get("expanded_root_error", {})
    ):
        issues.append("recursive and expanded exact root errors differ")
    if fixture.get("direct_shell_level") != 2:
        issues.append("direct shell level drifted")
    for marker in (
        "finite level sums",
        "skipped local defects vanish",
        "No complexity claim",
        "denominator-weighted kernel",
    ):
        if not any(marker in str(value) for value in routes.values()):
            issues.append(f"exact-defect route marker missing: {marker}")

    by_id = {row.get("id"): row for row in artifact.get("rows", [])}
    if by_id.get("held_08_exact_defect", {}).get("readiness") != "available_exact":
        issues.append("exact local-defect identity was weakened")
    if by_id.get("held_10_exact_shell", {}).get("readiness") != "available_conditional":
        issues.append("exact-shell implication readiness drifted")
    if by_id.get("held_16_physical_block", {}).get("readiness") != "not_ready_to_apply":
        issues.append("physical application was promoted")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove any W1 tail",
        "transformed-level adapter",
        "evaluator C6 or disk theorem",
        "outer Hardy representation bound",
        "Lambda<=0",
        "RH",
        "prize-level conclusion",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = note_path.read_text(encoding="utf-8")
    for marker in (
        "Hardy Local-Defect Obligation Ledger Gate",
        ARXIV_ID,
        "exact through equation (69)",
        "`P=3`",
        "`ip=3`",
        "`ip~20`",
        "epsilon_l = S_l - [a_l S_(l+1) + q_l]",
        "directly summing the first parent above the kernel",
        "termwise `W1` bounds",
        "0 external error bounds",
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
    print(
        "validated Newman C1 Hardy local-defect obligation ledger gate: "
        "16 rows, 6 paper-to-source mappings, 13 error components, 2 exact "
        "defect identities, 1 exact-shell bypass, 2 conditional derivative "
        "routes and 0 external error bounds"
    )


if __name__ == "__main__":
    main()
