#!/usr/bin/env python3
"""Validate signed order-seven forward invariance through lambda zero."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    ENTRY_SOURCE,
    FLOW_SOURCE,
    LOWER_SOURCE,
    TRANSFER_SOURCE,
    build_artifact,
)


@dataclass(frozen=True)
class ForwardIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> ForwardIssue:
    return ForwardIssue(section, code, str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_artifact(path: Path) -> tuple[dict, list[ForwardIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = load_json(path)
    issues: list[ForwardIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != (
        "signed contiguous and arbitrary-column order-seven theorem on "
        "-100<=lambda<=0"
    ):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

    expected_summary = {
        "rows": 8,
        "ready_rows": 7,
        "open_rows": 1,
        "uniform_tail_theorems": 1,
        "cooperative_flow_theorems": 1,
        "contiguous_order_seven_theorems": 1,
        "arbitrary_column_order_seven_theorems": 1,
        "open_all_order_handoffs": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    sources = artifact.get("source_contract", {})
    expected_sources = {
        "endpoint_entry": (
            "Q_(7,n)(-100)=-H_(7,n)(-100)>0 for every integer n>=0"
        ),
        "lower_layers": (
            "epsilon_k*H_(k,n)(lambda)>0 for k=1,2,3,4,5,6, all n and "
            "-100<=lambda<=0"
        ),
    }
    for key, expected in expected_sources.items():
        if sources.get(key) != expected:
            issues.append(issue("sources", f"bad-{key}", sources.get(key)))
    if "a_n=(4*n+50)" not in str(sources.get("cooperative_flow", "")):
        issues.append(issue("sources", "bad-cooperative-flow", sources))
    if "there exists N_7" not in str(sources.get("uniform_tail", "")):
        issues.append(issue("sources", "bad-uniform-tail", sources))

    theorem = artifact.get("conclusions", {})
    expected_theorems = {
        "contiguous_order_seven": (
            "Q_(7,n)(lambda)=-H_(7,n)(lambda)>0 for every n>=0 and "
            "-100<=lambda<=0"
        ),
        "signed_layers_through_seven": (
            "epsilon_k*H_(k,n)(lambda)>0 for 1<=k<=7, every n>=0 and "
            "-100<=lambda<=0"
        ),
    }
    for key, expected in expected_theorems.items():
        if theorem.get(key) != expected:
            issues.append(issue("theorem", f"bad-{key}", theorem.get(key)))
    if "epsilon_7*R_(7,n)" not in str(theorem.get("arbitrary_order_seven", "")):
        issues.append(issue("theorem", "bad-arbitrary-order-seven", theorem))

    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    if len(rows) != 8 or len(ids) != len(set(ids)):
        issues.append(issue("rows", "count-or-uniqueness", ids))
    if sum(row.get("readiness") == "ready_to_apply" for row in rows) != 7:
        issues.append(issue("rows", "bad-ready-count", rows))
    if sum(row.get("readiness") == "not_ready_to_apply" for row in rows) != 1:
        issues.append(issue("rows", "bad-open-count", rows))
    for source in (ENTRY_SOURCE, FLOW_SOURCE, LOWER_SOURCE, TRANSFER_SOURCE):
        if not source.exists():
            issues.append(issue("sources", "missing-source", source))

    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[ForwardIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: signed contiguous and arbitrary-column order-seven theorem",
        "not PF-infinity",
        "Q_(7,n)(-100)=-H_(7,n)(-100)>0 for every integer n>=0",
        "there exists N_7",
        "Q_n'=a_n*Q_(n+1)+b_n*Q_n",
        "Q_(7,n)(lambda)=-H_(7,n)(lambda)>0 for every n>=0",
        "epsilon_7*R_(7,n)",
        "contiguous and arbitrary-column order seven are both closed",
        "uniform-in-order mechanism",
        "outputs/jensen_window_pf_compound_order7_m100_entry_certificate.md",
        "outputs/formal_core.md",
    )
    issues: list[ForwardIssue] = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "pf-infinity is proved",
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
                f"ORDER7-UNIFORM-FORWARD {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-seven uniform-heat forward certificate: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-seven uniform-heat forward invariance: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['contiguous_order_seven_theorems']} contiguous theorem, "
        f"{summary['arbitrary_column_order_seven_theorems']} arbitrary-column theorem, "
        f"{summary['open_all_order_handoffs']} open all-order handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
