#!/usr/bin/env python3
"""Validate the rigorous lambda=-100 signed order-nine prefix."""

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

from jensen_window_pf_compound_order9_m100_prefix_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    MAX_COEFFICIENT_INDEX,
    PREFIX_LAST_N,
    build_artifact,
)


@dataclass(frozen=True)
class PrefixIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> PrefixIssue:
    return PrefixIssue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[PrefixIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[PrefixIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order9_m100_prefix_certificate":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "prefix through n=1240" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

    expected_summary = {
        "rows": 8,
        "ready_to_apply_rows": 7,
        "coefficients": 1257,
        "positive_Q8_rows": 1243,
        "positive_relative_Q8_rows": 1241,
        "positive_Q9_rows": 1241,
        "inconclusive_rows": 0,
        "precision_repair_rows": 38,
        "open_analytic_tails": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    exact = artifact.get("exact", {})
    expected_exact = {
        "signed_condensation": (
            "Q_(9,n)*Q_(7,n+2)=Q_(8,n+1)^2-Q_(8,n)*Q_(8,n+2)"
        ),
        "stable_coordinate": (
            "M_n=Q_(8,n+1)^2/(Q_(8,n)*Q_(8,n+2))-1"
        ),
        "prefix": "Q_(9,n)(-100)>0 for every 0<=n<=1240",
        "remaining_tail": "Q_(9,n)(-100)>0 for every n>=1241",
    }
    for key, expected in expected_exact.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))

    finite = artifact.get("finite", {})
    if finite.get("n_range") != [0, PREFIX_LAST_N]:
        issues.append(issue("finite", "bad-n-range", finite.get("n_range")))
    if finite.get("coefficient_range") != [0, MAX_COEFFICIENT_INDEX]:
        issues.append(issue("finite", "bad-coefficient-range", finite.get("coefficient_range")))
    rows = finite.get("rows", [])
    if len(rows) != 1241 or [row.get("n") for row in rows] != list(range(1241)):
        issues.append(issue("finite", "bad-row-cover", len(rows)))
    if any(row.get("Q9_sign") != "positive_by_signed_condensation" for row in rows):
        issues.append(issue("finite", "bad-q9-sign", rows))
    if finite.get("minimum_relative_n") != 1240:
        issues.append(issue("finite", "bad-minimum-index", finite.get("minimum_relative_n")))
    try:
        minimum_lower = Decimal(finite.get("minimum_relative_lower", "nan"))
    except Exception as exc:
        issues.append(issue("finite", "bad-minimum-lower", exc))
    else:
        if not minimum_lower.is_finite() or minimum_lower <= Decimal(1) / Decimal(250):
            issues.append(issue("finite", "minimum-not-above-one-over-250", minimum_lower))

    artifact_rows = artifact.get("rows", [])
    if len(artifact_rows) != 8 or len({row.get("id") for row in artifact_rows}) != 8:
        issues.append(issue("rows", "bad-rows", artifact_rows))
    try:
        canonical = build_artifact()
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[PrefixIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous signed order-nine endpoint prefix through `n=1240`",
        "This is not a proof",
        "Q_(9,n)*Q_(7,n+2)=Q_(8,n+1)^2-Q_(8,n)*Q_(8,n+2)",
        "M_n=Q_(8,n+1)^2/(Q_(8,n)*Q_(8,n+2))-1",
        "A_153,...,A_178",
        "All 38 repaired rows",
        "M_n(-100)>0 for every 0<=n<=1240",
        "Q_(9,n)(-100)>0 for every 0<=n<=1240",
        "n=1240",
        ">1/250",
        "Q_(9,n)(-100)>0 for every n>=1241",
        "outputs/jensen_window_pf_compound_order9_uniform_tail_flow_reduction.md",
        "outputs/formal_core.md",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "all-shift order nine is proved",
        "pf-infinity follows",
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
                f"ORDER9-M100-PREFIX {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-nine lambda=-100 prefix certificate: {len(issues)} issues")
        return 1
    finite = artifact["finite"]
    summary = artifact["summary"]
    print(
        "validated order-nine lambda=-100 prefix: "
        f"{summary['coefficients']} coefficients, "
        f"{summary['positive_Q8_rows']} positive Q8 values, "
        f"{summary['positive_relative_Q8_rows']} positive relative Q8 margins, "
        f"{summary['positive_Q9_rows']} positive Q9 signs, "
        f"{summary['inconclusive_rows']} inconclusive, "
        f"{summary['precision_repair_rows']} repaired coefficients, "
        f"{summary['open_analytic_tails']} open analytic tail, "
        f"weakest n={finite['minimum_relative_n']} above 1/250"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
