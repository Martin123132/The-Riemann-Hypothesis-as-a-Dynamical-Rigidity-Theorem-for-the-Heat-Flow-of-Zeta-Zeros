#!/usr/bin/env python3
"""Validate the all-shift order-eight lambda=-100 entry certificate."""

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

from jensen_window_pf_compound_order8_m100_entry_certificate import (  # noqa: E402
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
class Issue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> Issue:
    return Issue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[Issue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[Issue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order8_m100_entry_certificate":
        issues.append(issue("artifact", "kind", artifact.get("kind")))
    if artifact.get("status") != (
        "all-shift signed contiguous order-eight entry theorem at lambda=-100"
    ):
        issues.append(issue("artifact", "status", artifact.get("status")))

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
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(issue("summary", key, summary.get(key)))

    sources = artifact.get("source_contract", {})
    expected_sources = {
        "prefix": "Q_(8,n)(-100)>0 for every 0<=n<=1242",
        "shifted": "s_1''(t)<=2000/t^2 for every 699<=t<=999",
        "compact": "s_1''(t)<=4000/t^2 for every real 999<=t<=V'(2)",
        "finite_ray": "s_1''(t)<=4000/t^2 for every saddle mode 2<=u<=20",
        "asymptotic_ray": "t^2*s_1''(t)<200<4000 for u>=20",
        "coverage": "[999<=t<=V'(2)] union [2<=u<=20] union [u>=20]",
        "tent_transfer": (
            "s_1''(t)<=4000/t^2 => W_k^(1)<=4000*[-log(1-1/k^2)]"
            "<4001/k^2, k>=1250"
        ),
        "full_transfer": (
            "|W_k-W_k^(1)|<=C_(k-1)+2*C_k+C_(k+1)<190/k^2, k>=1250"
        ),
        "complete_ceiling": "W_k<=4300/k^2 for every real/integer k>=1250",
    }
    for key, value in expected_sources.items():
        if sources.get(key) != value:
            issues.append(issue("sources", key, sources.get(key)))
    try:
        asymptotic_upper = Decimal(sources.get("asymptotic_scaled_upper", "nan"))
    except Exception as exc:
        issues.append(issue("sources", "asymptotic-upper", exc))
    else:
        if not asymptotic_upper.is_finite() or asymptotic_upper >= 200:
            issues.append(issue("sources", "asymptotic-upper", asymptotic_upper))

    exact = artifact.get("exact", {})
    expected_exact = {
        "global_first_curvature": "s_1''(t)<=4000/t^2 for every real t>=699",
        "bridge_curvature": "s_1''(t)<=4000/t^2 for every real t>=999",
        "full_ceiling": (
            "W_k<W_k^(1)+|W_k-W_k^(1)|<4001/k^2+190/k^2="
            "4191/k^2<4300/k^2, k>=1250"
        ),
        "tail_sign": (
            "W_k<4300/k^2 => Q_(8,n)(-100)>0, k=n+7, n>=1243"
        ),
        "all_shift_entry": (
            "Q_(8,n)(-100)=H_(8,n)(-100)>0 for every integer n>=0"
        ),
    }
    for key, value in expected_exact.items():
        if exact.get(key) != value:
            issues.append(issue("exact", key, exact.get(key)))
    if not 4001 + 190 == 4191 < 4300:
        issues.append(issue("exact", "scalar-arithmetic", "failed"))
    if set(range(0, 1243)).intersection(range(1243, 1250)):
        issues.append(issue("exact", "prefix-tail", "unexpected overlap"))
    if max(range(0, 1243)) + 1 != 1243:
        issues.append(issue("exact", "prefix-tail", "index gap"))

    rows = artifact.get("rows", [])
    ids = [row.get("id") for row in rows]
    if len(rows) != 11 or len(ids) != len(set(ids)):
        issues.append(issue("rows", "count-or-unique", ids))
    if any(row.get("readiness") != "ready_to_apply" for row in rows):
        issues.append(issue("rows", "not-ready", rows))
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
            issues.append(issue("sources", "missing", source))
    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover
        issues.append(issue("rebuild", "exception", exc))
    else:
        if canonical != artifact:
            issues.append(issue("rebuild", "drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[Issue]:
    if not path.exists():
        return [issue("note", "missing", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: all-shift signed contiguous order-eight entry theorem at",
        "This certificate is not a proof",
        "s_1''(t)<=4000/t^2 for every real t>=699",
        "17999 rational mode blocks",
        "4001/k^2+190/k^2=4191/k^2<4300/k^2",
        "n>=1243",
        "0<=n<=1242",
        "Q_(8,n)(-100)=H_(8,n)(-100)>0 for every integer n>=0",
        "outputs/jensen_window_pf_compound_order8_uniform_tail_flow_reduction.md",
        "outputs/formal_core.md",
    )
    return [
        issue("note", "marker", marker)
        for marker in required
        if marker not in text
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    for item in issues:
        print(f"ORDER8-M100-ENTRY {item.section} [{item.code}] {item.detail}")
    if issues:
        print(f"order-eight lambda=-100 entry certificate: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-eight lambda=-100 entry certificate: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['continuous_curvature_theorems']} continuous curvature theorem, "
        f"{summary['analytic_tail_theorems']} analytic tail theorem, "
        f"{summary['all_shift_m100_entry_theorems']} all-shift entry theorem, "
        f"{summary['open_rows']} open rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
