#!/usr/bin/env python3
"""Validate the shifted-Hardy derivative-to-mode handoff and source audit."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from jensen_window_pf_newman_c1_hardy_modewise_derivative_handoff_gate import (
    ARXIV_ID,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    KIND,
    RESUMABLE_SOURCE,
    SOURCE_MARKERS,
    UPSTREAM_COMMIT,
    build_artifact,
    build_exact_mode_handoff,
    exact_orthogonal_modes,
    file_hash,
)


EXPECTED_IDS = [
    f"hmdh_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "primary_source"),
        (2, "paper_tolerance"),
        (3, "mgs_reduction"),
        (4, "missing_iteration_constant"),
        (5, "hybrid_error"),
        (6, "sample_boundary"),
        (7, "multi_shift_structure"),
        (8, "exact_grid_modes"),
        (9, "derivative_handoff"),
        (10, "analytic_disk"),
        (11, "transition_guard"),
        (12, "open_target"),
    )
]


def validate(path: Path, note_path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != KIND:
        issues.append("artifact kind mismatch")
    if artifact.get("status") != (
        "exact_conditional_handoff_with_open_external_error_theorem"
    ):
        issues.append("artifact status mismatch")
    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")

    rebuilt = build_artifact()
    if artifact != rebuilt:
        issues.append("artifact differs from an independent deterministic rebuild")

    handoff = artifact.get("exact_mode_handoff", {})
    if handoff != build_exact_mode_handoff():
        issues.append("exact mode handoff drifted")
    if handoff.get("shift_grid") != [
        f"{Fraction(index, 100).numerator}/{Fraction(index, 100).denominator}"
        for index in range(-7, 8)
    ]:
        issues.append("shift grid drifted")

    modes = exact_orthogonal_modes()
    artifact_modes = handoff.get("orthogonal_modes", [])
    if len(artifact_modes) != 6:
        issues.append("orthogonal mode count drifted")
    for degree, mode in enumerate(modes):
        if degree >= len(artifact_modes):
            break
        row = artifact_modes[degree]
        if row.get("degree") != degree or row.get("integer_vector") != mode:
            issues.append(f"mode {degree} vector drifted")
            continue
        norm_squared = sum(value * value for value in mode)
        if row.get("norm_squared") != str(norm_squared):
            issues.append(f"mode {degree} norm drifted")
        if math.gcd(*[abs(value) for value in mode]) != 1:
            issues.append(f"mode {degree} is not primitive")
        weights = row.get("derivative_weights", {})
        for lower_degree in range(degree):
            moment = sum(
                Fraction(value) * Fraction(index, 100) ** lower_degree
                for index, value in zip(range(-7, 8), mode, strict=True)
            )
            if moment != 0:
                issues.append(f"mode {degree} fails degree {lower_degree} orthogonality")
        expected_orders = []
        for derivative in range(degree, 6):
            moment = sum(
                Fraction(value) * Fraction(index, 100) ** derivative
                for index, value in zip(range(-7, 8), mode, strict=True)
            )
            if moment == 0:
                continue
            expected_orders.append(str(derivative))
            record = weights.get(str(derivative), {})
            expected_numerator = abs(moment) / math.factorial(derivative)
            if Fraction(record.get("numerator", "0/1")) != expected_numerator:
                issues.append(f"mode {degree} derivative {derivative} weight drifted")
            if record.get("norm_squared") != str(norm_squared):
                issues.append(f"mode {degree} derivative {derivative} norm drifted")
        if list(weights) != expected_orders:
            issues.append(f"mode {degree} active derivative orders drifted")
        remainder_numerator = sum(
            Fraction(abs(value)) * abs(Fraction(index, 100)) ** 6
            for index, value in zip(range(-7, 8), mode, strict=True)
        ) / math.factorial(6)
        remainder = row.get("sixth_derivative_remainder_weight", {})
        if Fraction(remainder.get("numerator", "0/1")) != remainder_numerator:
            issues.append(f"mode {degree} sixth-order remainder weight drifted")
        if remainder.get("norm_squared") != str(norm_squared):
            issues.append(f"mode {degree} sixth-order remainder norm drifted")
    tail_squared = Fraction(
        handoff.get("rough_tail_weight", {}).get("exact_squared", "0/1")
    )
    expected_tail_squared = sum(
        Fraction(index, 100) ** 12 for index in range(-7, 8)
    ) / 720**2
    if tail_squared != expected_tail_squared:
        issues.append("rough-tail exact constant drifted")

    source = artifact.get("source_audit", {})
    paper = source.get("paper", {})
    if paper.get("arxiv_id") != ARXIV_ID:
        issues.append("paper version drifted")
    if len(paper.get("audited_locations", [])) != 4:
        issues.append("primary-paper audit location count drifted")
    for marker in (
        "neither an explicit higher-order local-iteration constant",
        "uniform six-derivative",
        "fifteen-value shifted evaluator",
    ):
        if marker not in paper.get("audit_conclusion", ""):
            issues.append(f"paper audit conclusion marker missing: {marker}")

    code = source.get("code", {})
    if code.get("upstream_commit") != UPSTREAM_COMMIT:
        issues.append("upstream commit drifted")
    if code.get("audited_derivative_sha256") != file_hash(RESUMABLE_SOURCE):
        issues.append("audited derivative source hash drifted")
    source_text = RESUMABLE_SOURCE.read_text(encoding="utf-8")
    for marker in SOURCE_MARKERS:
        if marker not in source_text:
            issues.append(f"source marker missing: {marker}")

    by_id = {row.get("id"): row for row in artifact.get("rows", [])}
    if by_id.get("hmdh_04_missing_iteration_constant", {}).get("readiness") != "open":
        issues.append("missing local-iteration constant was promoted")
    if by_id.get("hmdh_12_open_target", {}).get("readiness") != "not_ready_to_apply":
        issues.append("physical-height handoff was promoted")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove any evaluator derivative bound",
        "complex-disk bound",
        "transition exclusion",
        "physical carrier value",
        "Lambda<=0",
        "RH",
        "prize-level conclusion",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = note_path.read_text(encoding="utf-8")
    for marker in (
        "Hardy Modewise Derivative Handoff Gate",
        ARXIV_ID,
        "equations (120)--(123)",
        "does not formulate one",
        "no hard-and-fast rule",
        "sup_(|u-T|<=7/100)",
        "Cauchy's estimate",
        "transition-aware piecewise substitute",
        "scalar two-height route remains retired",
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
        "validated Newman C1 Hardy modewise derivative handoff gate: "
        "12 rows, 6 exact low modes, 7 derivative orders, 1 exact rough-tail "
        "constant, 4 primary-paper locations, 7 source markers, 2 conditional "
        "routes and 1 open external error theorem"
    )


if __name__ == "__main__":
    main()
