#!/usr/bin/env python3
"""Validate the rigorous endpoint order-ten counterexample artifact."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_endpoint_order10_counterexample import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    SOURCE_PATHS,
    build_artifact,
    sha256,
)


@dataclass(frozen=True)
class CounterexampleIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> CounterexampleIssue:
    return CounterexampleIssue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[CounterexampleIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[CounterexampleIssue] = []
    if artifact.get("kind") != "jensen_window_pf_endpoint_order10_counterexample":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    status = str(artifact.get("status", ""))
    for marker in ("first-open-order", "counterexample", "rejecting"):
        if marker not in status:
            issues.append(issue("artifact", "bad-status", marker))
    boundary = str(artifact.get("proof_boundary", ""))
    for marker in (
        "four negative",
        "does not disprove Jensen hyperbolicity",
        "RH",
        "Lambda<=0",
        "replacement Xi/Phi bridge",
    ):
        if marker not in boundary:
            issues.append(issue("artifact", "bad-proof-boundary", marker))

    expected_summary = {
        "rows": 8,
        "direct_determinant_checks": 5,
        "negative_deep_rectangles": 4,
        "first_positive_shift": 4,
        "scanned_order10_rows": 1241,
        "positive_scanned_rows": 1237,
        "negative_scanned_rows": 4,
        "inconclusive_scanned_rows": 0,
        "rejected_all_order_endpoint_hierarchies": 1,
        "surviving_xi_specific_bridge_targets": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    expected_ids = {f"eoc10_{index:02d}_{suffix}" for index, suffix in (
        (1, "orientation"),
        (2, "interval_input"),
        (3, "toda_coordinate"),
        (4, "direct_counterexample"),
        (5, "schur_counterexample"),
        (6, "finite_sign_map"),
        (7, "rejected_hierarchy"),
        (8, "rerouted_bridge"),
    )}
    rows = artifact.get("rows", [])
    row_map = {row.get("id"): row for row in rows}
    if set(row_map) != expected_ids or len(rows) != len(expected_ids):
        issues.append(issue("rows", "bad-id-set", sorted(str(key) for key in row_map)))
    if row_map.get("eoc10_07_rejected_hierarchy", {}).get("role") != "forbidden_promotion":
        issues.append(issue("rows", "hierarchy-not-rejected", row_map.get("eoc10_07_rejected_hierarchy")))
    if row_map.get("eoc10_08_rerouted_bridge", {}).get("readiness") != "not_ready_to_apply":
        issues.append(issue("rows", "bridge-overpromoted", row_map.get("eoc10_08_rerouted_bridge")))

    exact = artifact.get("exact", {})
    expected_exact = {
        "orientation": "epsilon_10=(-1)^45=-1, Q_(10,n)=-H_(10,n)",
        "counterexample": (
            "Q_(10,n)(-100)<0 and s_(((n+9)^10))(h)<0 for n=0,1,2,3"
        ),
        "finite_positive_range": (
            "Q_(10,n)(-100)>0 for every 4<=n<=1240"
        ),
        "rejected_hierarchy": (
            "s_((N^m))(h)>0 for every m>=10,N>=m-1 is false"
        ),
    }
    for key, expected in expected_exact.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))

    direct = artifact.get("direct_audit", [])
    if [row.get("n") for row in direct] != [0, 1, 2, 3, 4]:
        issues.append(issue("direct", "bad-index-set", [row.get("n") for row in direct]))
    if [row.get("expected_sign") for row in direct] != [
        "negative", "negative", "negative", "negative", "positive"
    ]:
        issues.append(issue("direct", "bad-sign-pattern", [row.get("expected_sign") for row in direct]))
    for row in direct:
        checks = row.get("checks", {})
        if not checks or not all(value is True for value in checks.values()):
            issues.append(issue("direct", "failed-check", {"n": row.get("n"), "checks": checks}))

    finite = artifact.get("finite_scan", {})
    if finite.get("negative_indices") != [0, 1, 2, 3]:
        issues.append(issue("finite", "bad-negative-indices", finite.get("negative_indices")))
    if finite.get("positive_range") != [4, 1240]:
        issues.append(issue("finite", "bad-positive-range", finite.get("positive_range")))
    if finite.get("inconclusive_rows") != 0:
        issues.append(issue("finite", "inconclusive-rows", finite.get("inconclusive_rows")))

    source_rows = artifact.get("sources", [])
    source_map = {row.get("path"): row for row in source_rows}
    if len(source_rows) != len(SOURCE_PATHS):
        issues.append(issue("sources", "bad-count", len(source_rows)))
    for source in SOURCE_PATHS:
        relative = source.relative_to(source.parents[3]).as_posix()
        row = source_map.get(relative)
        if row is None:
            issues.append(issue("sources", "missing-source", relative))
        elif row.get("sha256") != sha256(source):
            issues.append(issue("sources", "hash-mismatch", relative))

    try:
        canonical = build_artifact()
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if canonical != artifact:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[CounterexampleIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous first-open-order counterexample",
        "counterexample to RH",
        "epsilon_10=(-1)^45=-1",
        "Q_(10,n)*Q_(8,n+2)",
        "(9^10),(10^10),(11^10),(12^10)",
        "positive: 4<=n<=1240",
        "s_((N^m))(h)>0 for every m>=10,N>=m-1 is false",
        "all-order positivity antecedent is",
        "weaker Xi/Phi-specific route",
        "route correction",
    )
    issues: list[CounterexampleIssue] = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh is false",
        "rh is disproved",
        "jensen hyperbolicity is false",
        "the replacement bridge is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ENDPOINT-ORDER10 {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"endpoint order-ten counterexample: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated endpoint order-ten counterexample: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['direct_determinant_checks']} direct checks, "
        f"{summary['negative_deep_rectangles']} negative deep rectangles, "
        f"{summary['positive_scanned_rows']} positive scanned rows, "
        f"{summary['inconclusive_scanned_rows']} inconclusive, "
        "1 rejected all-order endpoint hierarchy, "
        "1 surviving Xi/Phi bridge target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
