#!/usr/bin/env python3
"""Validate the stored Q208 top phase-cell certificate and cache."""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

import check_jensen_window_pf_newman_q208_bottom_phase_cell_certificate as bottom_check
import jensen_window_pf_newman_q208_top_phase_cell_certificate as cert


def validate() -> list[str]:
    issues: list[str] = []
    core = cert.core
    core.compact.flint.ctx.prec = core.bridge.PRECISION_BITS
    if not cert.DEFAULT_OUT.exists():
        return ["stored top result is missing"]
    if not cert.DEFAULT_NOTE.exists():
        issues.append("rendered top note is missing")
    artifact = json.loads(cert.DEFAULT_OUT.read_text(encoding="utf-8"))
    if artifact.get("kind") != cert.STEM:
        issues.append("top artifact kind drifted")
    if artifact.get("contract") != cert.contract_payload():
        issues.append("top contract payload drifted")
    if artifact.get("contract_sha256") != cert.contract_hash():
        issues.append("top contract hash drifted")
    if artifact.get("builder_sha256") != core.file_hash(
        Path(cert.__file__).resolve()
    ):
        issues.append("top builder hash drifted")
    if artifact.get("contract", {}).get("source_sha256") != (
        cert.source_hashes()
    ):
        issues.append("top source hashes drifted")

    tasks = core.panel_tasks()
    try:
        cache_records = cert.load_cache(cert.DEFAULT_CACHE, tasks)
    except (RuntimeError, ValueError, json.JSONDecodeError) as exc:
        issues.append(f"top cache validation failed: {exc}")
        cache_records = []
    if artifact.get("records", []) != cache_records:
        issues.append("top artifact records differ from cache")
    expected_cache_hash = (
        core.file_hash(cert.DEFAULT_CACHE)
        if cert.DEFAULT_CACHE.exists()
        else None
    )
    if artifact.get("cache_sha256") != expected_cache_hash:
        issues.append("top cache file hash drifted")

    tail_source = json.loads(
        core.bridge.TAIL_SOURCE.read_text(encoding="utf-8")
    )
    tail_witness = next(
        row
        for row in tail_source["witnesses"]
        if int(row["first_omitted"]) == core.bridge.FIRST_OMITTED
    )
    e0 = core.compact.arb(tail_witness["E0"]["enclosure"])
    e1 = core.compact.arb(tail_witness["E1"]["enclosure"])

    records = cache_records
    leaf_count = 0
    unresolved = 0
    branches: Counter = Counter()
    for index, record in enumerate(records):
        result = record.get("result", {})
        leaves = result.get("certified_leaves", [])
        leaf_count += len(leaves)
        unresolved += int(result.get("unresolved_boxes", 0))
        for leaf_index, leaf in enumerate(leaves):
            bottom_check.validate_leaf(
                leaf,
                issues,
                f"top panel {index + 1} leaf {leaf_index + 1}",
                e0,
                e1,
            )
            branches[leaf.get("branch")] += 1
        if len(leaves) != result.get("certified_leaf_boxes"):
            issues.append(f"top panel {index + 1} leaf count mismatch")
        if result.get("status") == "certified" and result.get(
            "unresolved_boxes"
        ) != 0:
            issues.append(f"top panel {index + 1} status is inconsistent")

    summary = artifact.get("summary", {})
    resource = summary.get("resource", {})
    if resource.get("worker_count") != 1:
        issues.append("top resource record does not use one worker")
    if resource.get("mode") != "day":
        issues.append("top resource record is not daytime mode")
    if len(resource.get("baseline_samples", [])) != core.BASELINE_SECONDS:
        issues.append("top baseline CPU sample count drifted")
    if summary.get("tasks_total") != len(tasks):
        issues.append("top task total drifted")
    if summary.get("tasks_completed") != len(records):
        issues.append("top completed count drifted")
    if summary.get("certified_leaf_cells") != leaf_count:
        issues.append("top leaf count drifted")
    if summary.get("unresolved_boxes") != unresolved:
        issues.append("top unresolved count drifted")
    if summary.get("branch_counts") != dict(branches):
        issues.append("top branch counts drifted")

    complete = len(records) == len(tasks) and unresolved == 0
    if summary.get("complete_top_edge_certificate") != complete:
        issues.append("top completion flag drifted")
    expected_status = (
        "rigorous complete Q208 top-edge phase-cell certificate"
        if complete
        else "resumable partial Q208 top-edge phase-cell computation"
    )
    if artifact.get("status") != expected_status:
        issues.append("top status drifted")
    if complete:
        bottom_check.validate_complete_chain(
            artifact,
            records,
            issues,
        )
    elif summary.get("open_phase_chain") is not None:
        issues.append("partial top artifact promoted a phase chain")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "certifies only the fixed-time top path",
        "bottom chain",
        "cyclic winding",
        "Q208",
        "Lambda<=0",
        "RH",
        "Clay prize",
        "remain open",
    ):
        if marker not in boundary:
            issues.append(f"top proof boundary missing: {marker}")
    if cert.DEFAULT_NOTE.exists():
        note = cert.DEFAULT_NOTE.read_text(encoding="utf-8")
        for marker in (
            "Contract",
            "Progress",
            "not a proof of Q208",
            "bottom chain",
            "cyclic exact",
        ):
            if marker not in note:
                issues.append(f"top note marker missing: {marker}")
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
        "validated Q208 top phase-cell certificate: "
        f"{summary['tasks_completed']}/{summary['tasks_total']} panels, "
        f"{summary['certified_leaf_cells']} cells, "
        f"{summary['unresolved_boxes']} unresolved, "
        f"complete={summary['complete_top_edge_certificate']}, "
        "0 structural issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
