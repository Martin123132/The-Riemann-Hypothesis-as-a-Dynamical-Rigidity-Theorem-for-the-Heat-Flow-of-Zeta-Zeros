#!/usr/bin/env python3
"""Validate the stored Q208 bottom phase-cell certificate and cache."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import json
from pathlib import Path

import jensen_window_pf_newman_q208_bottom_phase_cell_certificate as cert


def independent_open_crossing(
    vertices: list[tuple[Fraction, Fraction]],
) -> int:
    total = 0
    for left, right in zip(vertices, vertices[1:]):
        determinant = left[0] * right[1] - left[1] * right[0]
        dot = left[0] * right[0] + left[1] * right[1]
        if (
            left == (0, 0)
            or right == (0, 0)
            or (determinant == 0 and dot <= 0)
        ):
            raise ValueError("witness edge meets the origin")
        upward = left[1] <= 0 and right[1] > 0
        downward = right[1] <= 0 and left[1] > 0
        if upward and determinant > 0:
            total += 1
        if downward and determinant < 0:
            total -= 1
    return total


def validate_leaf(
    leaf: dict,
    issues: list[str],
    label: str,
    e0,
    e1,
) -> None:
    arb = cert.compact.arb
    try:
        x_low = Fraction(leaf["x_low"])
        x_high = Fraction(leaf["x_high"])
    except (KeyError, ValueError, ZeroDivisionError) as exc:
        issues.append(f"{label} has invalid rational bounds: {exc}")
        return
    if x_low >= x_high:
        issues.append(f"{label} has nonpositive width")
        return

    try:
        h = arb(leaf["h_retained"])
        h_prime = arb(leaf["h_prime_retained"])
        stored_f_retained = arb(leaf["f_retained"])
        stored_f_prime_retained = arb(leaf["f_prime_retained"])
        stored_tail_f = arb(leaf["tail_f_upper"])
        stored_tail_f_prime = arb(leaf["tail_f_prime_upper"])
        full_f = arb(leaf["full_f"])
        full_f_prime = arb(leaf["full_f_prime"])
    except (KeyError, ValueError) as exc:
        issues.append(f"{label} has an invalid Arb enclosure: {exc}")
        return

    x_box = cert.bridge.interval_ball(x_low, x_high)
    expected_f_retained = 16 * (1 + x_box**4) * h
    expected_f_prime_retained = (
        64 * x_box**3 * h
        + 16 * (1 + x_box**4) * h_prime
    )
    try:
        stored_f_retained.intersection(expected_f_retained)
    except ValueError:
        issues.append(f"{label} retained F formula does not overlap")
    try:
        stored_f_prime_retained.intersection(
            expected_f_prime_retained
        )
    except ValueError:
        issues.append(f"{label} retained F' formula does not overlap")

    x_ceiling = arb(cert.compact.fraction_decimal(x_high))
    expected_tail_f = (
        16 * (1 + x_ceiling**4) * e0
    ).upper()
    expected_tail_f_prime = (
        64 * x_ceiling**3 * e0
        + 16 * (1 + x_ceiling**4) * e1
    ).upper()
    if not stored_tail_f.contains(expected_tail_f):
        issues.append(f"{label} F tail is too small")
    if not stored_tail_f_prime.contains(expected_tail_f_prime):
        issues.append(f"{label} F' tail is too small")

    if not full_f.contains(stored_f_retained):
        issues.append(f"{label} full F does not contain retained F")
    if not full_f_prime.contains(stored_f_prime_retained):
        issues.append(f"{label} full F' does not contain retained F'")

    value_separated = full_f.abs_lower() > 0
    derivative_separated = full_f_prime.abs_lower() > 0
    if not value_separated and not derivative_separated:
        issues.append(f"{label} cell contains the origin")
    if leaf.get("certified") is not True:
        issues.append(f"{label} is not marked certified")
    expected_branches = {
        "value" if value_separated else None,
        "derivative" if derivative_separated else None,
    }
    if leaf.get("branch") not in expected_branches:
        issues.append(f"{label} branch is inconsistent")


def validate_complete_chain(
    artifact: dict,
    records: list[dict],
    issues: list[str],
) -> None:
    leaves = [
        leaf
        for record in records
        for leaf in record["result"]["certified_leaves"]
    ]
    leaves.sort(key=lambda row: Fraction(row["x_low"]))
    if not leaves:
        issues.append("complete artifact has no phase cells")
        return
    if Fraction(leaves[0]["x_low"]) != cert.X_LOWER:
        issues.append("complete chain does not start at x=0")
    if Fraction(leaves[-1]["x_high"]) != cert.X_UPPER:
        issues.append("complete chain does not end at x=246")
    for index, (left, right) in enumerate(zip(leaves, leaves[1:])):
        if Fraction(left["x_high"]) != Fraction(right["x_low"]):
            issues.append(f"cell cover gap at join {index + 1}")

    chain = artifact.get("summary", {}).get("open_phase_chain")
    if not isinstance(chain, dict):
        issues.append("complete artifact has no open phase chain")
        return
    raw_witnesses = chain.get("witnesses", [])
    try:
        witnesses = [
            (Fraction(point[0]), Fraction(point[1]))
            for point in raw_witnesses
        ]
    except (IndexError, ValueError, ZeroDivisionError) as exc:
        issues.append(f"phase witnesses are not exact rationals: {exc}")
        return
    if len(witnesses) != len(leaves) + 1:
        issues.append("phase witness count mismatch")
        return

    cells = [
        (
            cert.compact.arb(leaf["full_f"]),
            cert.compact.arb(leaf["full_f_prime"]),
        )
        for leaf in leaves
    ]
    for index, witness in enumerate(witnesses):
        if index > 0 and not cert.point_in_cell(
            witness,
            cells[index - 1],
        ):
            issues.append(f"witness {index} misses its previous cell")
        if index < len(cells) and not cert.point_in_cell(
            witness,
            cells[index],
        ):
            issues.append(f"witness {index} misses its next cell")
    if not (
        witnesses[0][0] > 0
        and witnesses[0][1] == 0
        and chain.get("axis_witness_positive_real") is True
    ):
        issues.append("axis-corner witness is not positive real")
    try:
        crossing = independent_open_crossing(witnesses)
    except ValueError as exc:
        issues.append(str(exc))
    else:
        if crossing != chain.get("open_positive_ray_crossing_count"):
            issues.append("open phase-chain crossing count drifted")


def validate() -> list[str]:
    issues: list[str] = []
    cert.compact.flint.ctx.prec = cert.bridge.PRECISION_BITS
    if not cert.DEFAULT_OUT.exists():
        return ["stored result is missing"]
    if not cert.DEFAULT_NOTE.exists():
        issues.append("rendered note is missing")
    artifact = json.loads(cert.DEFAULT_OUT.read_text(encoding="utf-8"))
    if artifact.get("kind") != cert.STEM:
        issues.append("artifact kind drifted")
    if artifact.get("contract") != cert.contract_payload():
        issues.append("contract payload drifted")
    if artifact.get("contract_sha256") != cert.contract_hash():
        issues.append("contract hash drifted")
    if artifact.get("builder_sha256") != cert.file_hash(
        Path(cert.__file__).resolve()
    ):
        issues.append("builder hash drifted")
    if artifact.get("contract", {}).get("source_sha256") != (
        cert.source_hashes()
    ):
        issues.append("source hashes drifted")

    tasks = cert.panel_tasks()
    try:
        cache_records = cert.load_cache(cert.DEFAULT_CACHE, tasks)
    except (RuntimeError, ValueError, json.JSONDecodeError) as exc:
        issues.append(f"cache validation failed: {exc}")
        cache_records = []
    stored_records = artifact.get("records", [])
    if stored_records != cache_records:
        issues.append("artifact records differ from hash-chained cache")
    records = cache_records
    expected_cache_hash = (
        cert.file_hash(cert.DEFAULT_CACHE)
        if cert.DEFAULT_CACHE.exists()
        else None
    )
    if artifact.get("cache_sha256") != expected_cache_hash:
        issues.append("cache file hash drifted")
    if len(records) > len(tasks):
        issues.append("too many deterministic panel records")

    leaf_count = 0
    unresolved = 0
    branches: Counter = Counter()
    tail_source = json.loads(
        cert.bridge.TAIL_SOURCE.read_text(encoding="utf-8")
    )
    tail_witness = next(
        row
        for row in tail_source["witnesses"]
        if int(row["first_omitted"]) == cert.bridge.FIRST_OMITTED
    )
    e0 = cert.compact.arb(tail_witness["E0"]["enclosure"])
    e1 = cert.compact.arb(tail_witness["E1"]["enclosure"])
    for index, record in enumerate(records):
        result = record.get("result", {})
        leaves = result.get("certified_leaves", [])
        leaf_count += len(leaves)
        unresolved += int(result.get("unresolved_boxes", 0))
        for leaf_index, leaf in enumerate(leaves):
            validate_leaf(
                leaf,
                issues,
                f"panel {index + 1} leaf {leaf_index + 1}",
                e0,
                e1,
            )
            branches[leaf.get("branch")] += 1
        if len(leaves) != result.get("certified_leaf_boxes"):
            issues.append(f"panel {index + 1} leaf count mismatch")
        if result.get("status") == "certified" and result.get(
            "unresolved_boxes"
        ) != 0:
            issues.append(f"panel {index + 1} status is inconsistent")

    summary = artifact.get("summary", {})
    resource = summary.get("resource", {})
    if resource.get("worker_count") != 1:
        issues.append("resource record does not use one worker")
    if resource.get("mode") != "day":
        issues.append("resource record is not daytime mode")
    if len(resource.get("baseline_samples", [])) != cert.BASELINE_SECONDS:
        issues.append("baseline CPU sample count drifted")
    if summary.get("tasks_total") != len(tasks):
        issues.append("summary task total drifted")
    if summary.get("tasks_completed") != len(records):
        issues.append("summary completed count drifted")
    if summary.get("certified_leaf_cells") != leaf_count:
        issues.append("summary leaf count drifted")
    if summary.get("unresolved_boxes") != unresolved:
        issues.append("summary unresolved count drifted")
    if summary.get("branch_counts") != dict(branches):
        issues.append("summary branch counts drifted")
    expected_complete = len(records) == len(tasks) and unresolved == 0
    if summary.get("complete_bottom_edge_certificate") != expected_complete:
        issues.append("completion flag drifted")
    expected_status = (
        "rigorous complete Q208 bottom-edge phase-cell certificate"
        if expected_complete
        else "resumable partial Q208 bottom-edge phase-cell computation"
    )
    if artifact.get("status") != expected_status:
        issues.append("artifact status drifted")
    if expected_complete:
        validate_complete_chain(artifact, records, issues)
    elif summary.get("open_phase_chain") is not None:
        issues.append("partial artifact promoted an open phase chain")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "certifies only the fixed-time bottom path",
        "top edge",
        "cyclic winding",
        "Q208",
        "Lambda<=0",
        "RH",
        "Clay prize",
        "remain open",
    ):
        if marker not in boundary:
            issues.append(f"proof boundary missing: {marker}")
    if cert.DEFAULT_NOTE.exists():
        note = cert.DEFAULT_NOTE.read_text(encoding="utf-8")
        for marker in (
            "Contract",
            "Progress",
            "not a proof of Q208",
            "top phase chain",
            "closed polygon winding remain open",
        ):
            if marker not in note:
                issues.append(f"note marker missing: {marker}")
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    artifact = json.loads(cert.DEFAULT_OUT.read_text(encoding="utf-8"))
    summary = artifact["summary"]
    print(
        "validated Q208 bottom phase-cell certificate: "
        f"{summary['tasks_completed']}/{summary['tasks_total']} panels, "
        f"{summary['certified_leaf_cells']} cells, "
        f"{summary['unresolved_boxes']} unresolved, "
        f"complete={summary['complete_bottom_edge_certificate']}, "
        "0 structural issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
