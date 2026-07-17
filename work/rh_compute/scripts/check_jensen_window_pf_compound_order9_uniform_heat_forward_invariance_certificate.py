#!/usr/bin/env python3
"""Validate signed order-nine forward invariance through lambda zero."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate as forward9  # noqa: E402


@dataclass(frozen=True)
class ForwardIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> ForwardIssue:
    return ForwardIssue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[ForwardIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[ForwardIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    expected_status = (
        "signed contiguous and arbitrary-column order-nine theorem on "
        "-100<=lambda<=0"
    )
    if artifact.get("status") != expected_status:
        issues.append(issue("artifact", "bad-status", artifact.get("status")))
    expected_summary = {
        "rows": 8,
        "ready_rows": 7,
        "open_rows": 1,
        "uniform_tail_theorems": 1,
        "cooperative_flow_theorems": 1,
        "contiguous_order_nine_theorems": 1,
        "arbitrary_column_order_nine_theorems": 1,
        "open_all_order_handoffs": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    sources = artifact.get("source_contract", {})
    expected_sources = {
        "endpoint_entry": forward9.ENDPOINT_ENTRY,
        "lower_layers": (
            "epsilon_k*H_(k,n)(lambda)>0 for k=1,2,3,4,5,6,7,8, all n "
            "and -100<=lambda<=0"
        ),
    }
    for key, expected in expected_sources.items():
        if sources.get(key) != expected:
            issues.append(issue("sources", f"bad-{key}", sources.get(key)))
    if "a_n=(4*n+66)" not in str(sources.get("cooperative_flow", "")):
        issues.append(issue("sources", "bad-cooperative-flow", sources))
    if "there exists N_9" not in str(sources.get("uniform_tail", "")):
        issues.append(issue("sources", "bad-uniform-tail", sources))
    bound_sources = sources.get("sources", [])
    if len(bound_sources) != 4:
        issues.append(issue("sources", "bad-count", len(bound_sources)))
    for source in bound_sources:
        try:
            source_path = forward9.REPO_ROOT / source["path"]
        except Exception as exc:
            issues.append(issue("sources", "bad-path", exc))
            continue
        if not source_path.exists():
            issues.append(issue("sources", "missing-source", source_path))
        elif forward9.sha256(source_path) != source.get("sha256"):
            issues.append(issue("sources", "hash-mismatch", source_path))

    theorem = artifact.get("conclusions", {})
    expected_theorems = forward9.conclusions()
    if theorem != expected_theorems:
        issues.append(issue("theorem", "conclusion-drift", theorem))
    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    if len(rows) != 8 or len(ids) != len(set(ids)):
        issues.append(issue("rows", "count-or-uniqueness", ids))
    if sum(row.get("readiness") == "ready_to_apply" for row in rows) != 7:
        issues.append(issue("rows", "bad-ready-count", rows))
    if sum(row.get("readiness") == "not_ready_to_apply" for row in rows) != 1:
        issues.append(issue("rows", "bad-open-count", rows))
    try:
        canonical = forward9.build_artifact()
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[ForwardIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    theorem = forward9.conclusions()
    required = (
        "Status: signed contiguous and arbitrary-column order-nine theorem",
        "This is not a proof",
        forward9.ENDPOINT_ENTRY,
        "there exists N_9",
        "a_n=(4*n+66)",
        theorem["contiguous_order_nine"],
        theorem["arbitrary_order_nine"],
        "not PF-infinity",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "pf-infinity follows",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=forward9.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=forward9.DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER9-FORWARD {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-nine heat-forward invariance: {len(issues)} issues")
        return 1
    print(
        "validated signed contiguous and arbitrary-column order nine on "
        "-100<=lambda<=0: 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
