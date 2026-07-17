#!/usr/bin/env python3
"""Validate the all-shift signed order-seven entry theorem at lambda=-100."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order7_m100_entry_certificate import (  # noqa: E402
    ASYMPTOTIC_SOURCE,
    BRIDGE_SOURCE,
    COMPACT_SOURCE,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    FINITE_SOURCE,
    PREFIX_SOURCE,
    SHIFTED_SOURCE,
    TAIL_SOURCE,
    build_artifact,
)


@dataclass(frozen=True)
class EntryIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> EntryIssue:
    return EntryIssue(section, code, str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_artifact(path: Path) -> tuple[dict, list[EntryIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = load_json(path)
    issues: list[EntryIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order7_m100_entry_certificate":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != (
        "all-shift signed contiguous order-seven entry theorem at lambda=-100"
    ):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

    expected_summary = {
        "rows": 11,
        "ready_rows": 11,
        "continuous_curvature_theorems": 1,
        "complete_scalar_ceiling_theorems": 1,
        "analytic_tail_theorems": 1,
        "all_shift_m100_entry_theorems": 1,
        "open_rows": 0,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    sources = artifact.get("source_contract", {})
    expected_sources = {
        "prefix": "Q_(7,n)(-100)>0 for every 0<=n<=314",
        "shifted": "r_1''(t)<=600/t^2 for 320<=t<=1000",
        "compact": "r_1''(t)<=600/t^2 for every real 1000<=t<=V'(2)",
        "finite_ray": "r_1''(t)<=600/t^2 for every saddle mode 2<=u<=20",
        "asymptotic_ray": "t^2*r_1''(t)<100<600 for u>=20",
        "complete_ceiling": "R_k<=900/k^2 for every real/integer k>=321",
    }
    for key, expected in expected_sources.items():
        if sources.get(key) != expected:
            issues.append(issue("sources", f"bad-{key}", sources.get(key)))
    try:
        asymptotic_upper = Decimal(sources.get("asymptotic_scaled_upper", "nan"))
    except Exception as exc:
        issues.append(issue("sources", "bad-asymptotic-upper", exc))
    else:
        if not asymptotic_upper.is_finite() or asymptotic_upper >= 100:
            issues.append(issue("sources", "asymptotic-upper-failed", asymptotic_upper))

    exact = artifact.get("exact", {})
    expected_exact = {
        "global_first_curvature": "r_1''(t)<=600/t^2 for every real t>=320",
        "tent_transfer": "R_k^(1)<600*[-log(1-1/k^2)]<601/k^2, k>=321",
        "full_ceiling": (
            "R_k<R_k^(1)+|R_k-R_k^(1)|<601/k^2+262/k^2="
            "863/k^2<900/k^2, k>=321"
        ),
        "tail_sign": "R_k<900/k^2 => Q_(7,n)(-100)>0, k=n+6, n>=315",
        "all_shift_entry": (
            "Q_(7,n)(-100)=-H_(7,n)(-100)>0 for every integer n>=0"
        ),
    }
    for key, expected in expected_exact.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))
    if "covers every real t>=320" not in str(exact.get("coverage", "")):
        issues.append(issue("exact", "bad-coverage", exact.get("coverage")))

    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    if len(rows) != 11 or len(ids) != len(set(ids)):
        issues.append(issue("rows", "count-or-uniqueness", ids))
    if any(row.get("readiness") != "ready_to_apply" for row in rows):
        issues.append(issue("rows", "non-ready-row", rows))
    for source in (
        PREFIX_SOURCE,
        TAIL_SOURCE,
        BRIDGE_SOURCE,
        SHIFTED_SOURCE,
        COMPACT_SOURCE,
        FINITE_SOURCE,
        ASYMPTOTIC_SOURCE,
    ):
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


def validate_note(path: Path) -> list[EntryIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: all-shift signed contiguous order-seven entry theorem",
        "This certificate is not a proof",
        "r_1''(t)<=600/t^2 for 320<=t<=1000",
        "r_1''(t)<=600/t^2 for every real 1000<=t<=V'(2)",
        "r_1''(t)<=600/t^2 for every saddle mode 2<=u<=20",
        "t^2*r_1''(t)<100<600 for u>=20",
        "r_1''(t)<=600/t^2 for every real t>=320",
        "<601/k^2",
        "<262/k^2",
        "863/k^2<900/k^2",
        "Q_(7,n)(-100)>0 for every 0<=n<=314",
        "Q_(7,n)(-100)=-H_(7,n)(-100)>0 for every integer n>=0",
        "outputs/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.md",
        "outputs/formal_core.md",
    )
    issues: list[EntryIssue] = []
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
                f"ORDER7-M100-ENTRY {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-seven lambda=-100 entry certificate: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-seven lambda=-100 entry certificate: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['continuous_curvature_theorems']} continuous curvature theorem, "
        f"{summary['complete_scalar_ceiling_theorems']} complete scalar ceiling, "
        f"{summary['analytic_tail_theorems']} analytic tail theorem, "
        f"{summary['all_shift_m100_entry_theorems']} all-shift entry theorem, "
        f"{summary['open_rows']} open rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
