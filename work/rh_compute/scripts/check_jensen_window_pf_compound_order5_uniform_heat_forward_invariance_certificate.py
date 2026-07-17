#!/usr/bin/env python3
"""Validate order-five forward invariance and arbitrary-column transfer."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate import (  # noqa: E402
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
    if artifact.get("kind") != "jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "contiguous and arbitrary-column order-five theorem on -100<=lambda<=0":
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

    expected_summary = {
        "rows": 8,
        "ready_rows": 7,
        "open_rows": 1,
        "uniform_tail_theorems": 1,
        "cooperative_flow_theorems": 1,
        "contiguous_order_five_theorems": 1,
        "arbitrary_column_order_five_theorems": 1,
        "open_order_six_handoffs": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    sources = artifact.get("source_contract", {})
    if sources.get("endpoint_entry") != "H_(5,n)(-100)>0 for every integer n>=0":
        issues.append(issue("sources", "bad-endpoint", sources))
    if "a_n=(4*n+34)*H_(4,n)/H_(4,n+1)>0" not in sources.get("cooperative_flow", ""):
        issues.append(issue("sources", "bad-flow", sources))
    if "E_n(lambda)" not in sources.get("variation_of_constants", ""):
        issues.append(issue("sources", "bad-variation", sources))
    if "epsilon_k R_(k,n)" not in sources.get("fixed_order_transfer", ""):
        issues.append(issue("sources", "bad-transfer", sources))

    conclusions = artifact.get("conclusions", {})
    expected = {
        "contiguous_order_five": "H_(5,n)(lambda)>0 for every n>=0 and -100<=lambda<=0",
        "signed_layers_through_five": "epsilon_k*H_(k,n)(lambda)>0 for 1<=k<=5, every n>=0 and -100<=lambda<=0",
        "arbitrary_order_five": "R_(5,n)(j_1,j_2,j_3,j_4,j_5;lambda)>0 for every n>=0, 0<=j_1<j_2<j_3<j_4<j_5, and -100<=lambda<=0",
    }
    for key, value in expected.items():
        if conclusions.get(key) != value:
            issues.append(issue("conclusions", f"bad-{key}", conclusions.get(key)))

    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    if len(rows) != 8 or len(ids) != len(set(ids)):
        issues.append(issue("rows", "count-or-uniqueness", ids))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or "order-six" not in open_rows[0].get("claim", ""):
        issues.append(issue("rows", "bad-open-handoff", open_rows))
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
        "Status: contiguous and arbitrary-column order-five theorem",
        "not PF-infinity",
        "H_(5,n)(-100)>0 for every integer n>=0",
        "there exists N_5",
        "Q_n'=a_n*Q_(n+1)+b_n*Q_n",
        "Variation of constants",
        "H_(5,n)(lambda)>0 for every n>=0 and -100<=lambda<=0",
        "R_(5,n)(j_1,j_2,j_3,j_4,j_5;lambda)>0",
        "order-six",
        "outputs/jensen_window_pf_compound_order5_m100_entry_certificate.md",
        "outputs/formal_core.md",
    )
    issues: list[ForwardIssue] = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "pf-infinity is proved"):
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
            print(f"ORDER5-FORWARD {finding.section} [{finding.code}] {finding.detail}")
        print(f"order-five uniform-heat forward invariance: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print("validated order-five uniform-heat forward invariance: "
          f"{summary['rows']} rows, 0 issues, {summary['ready_rows']} ready rows, "
          f"{summary['contiguous_order_five_theorems']} contiguous theorem, "
          f"{summary['arbitrary_column_order_five_theorems']} arbitrary-column theorem, "
          f"{summary['open_order_six_handoffs']} open order-six handoff")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
