#!/usr/bin/env python3
"""Validate the lambda=-100 order-ten finite splice certificate."""

from __future__ import annotations

import argparse
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from jensen_window_pf_compound_order10_m100_finite_splice_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_artifact,
)


def add(issues: list[str], condition: bool, message: str) -> None:
    if not condition:
        issues.append(message)


def validate_artifact(path: Path) -> tuple[dict, list[str]]:
    if not path.exists():
        return {}, [f"missing artifact: {path}"]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    add(
        issues,
        artifact.get("kind")
        == "jensen_window_pf_compound_order10_m100_finite_splice_certificate",
        f"bad kind: {artifact.get('kind')!r}",
    )
    add(
        issues,
        "two-index order-ten finite splice" in str(artifact.get("status", "")),
        f"bad status: {artifact.get('status')!r}",
    )
    expected_summary = {
        "rows": 5,
        "ready_rows": 5,
        "open_rows": 0,
        "coefficients": 1261,
        "new_coefficient_rows": 2,
        "direct_Q10_checks": 2,
        "direct_schur_checks": 2,
        "toda_consistency_checks": 2,
        "finite_splice_rows": 2,
        "combined_positive_Q10_rows": 1239,
        "preserved_negative_endpoint_shifts": 4,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        add(
            issues,
            summary.get(key) == expected,
            f"bad summary {key}: {summary.get(key)!r}",
        )
    exact = artifact.get("exact", {})
    required_exact = {
        "splice": "Q_(10,n)(-100)>0 for n=1241,1242",
        "combined_positive_block": (
            "Q_(10,n)(-100)>0 for every 4<=n<=1242"
        ),
        "preserved_negative_prefix": (
            "Q_(10,n)(-100)<0 for n=0,1,2,3"
        ),
        "next_analytic_index": "n>=1243, equivalently k=n+9>=1252",
    }
    for key, expected in required_exact.items():
        add(issues, exact.get(key) == expected, f"bad exact {key}: {exact.get(key)!r}")
    finite = artifact.get("finite", {})
    add(
        issues,
        finite.get("coefficient_range") == [0, 1260],
        f"bad coefficient range: {finite.get('coefficient_range')!r}",
    )
    add(
        issues,
        finite.get("combined_positive_range") == [4, 1242],
        f"bad positive range: {finite.get('combined_positive_range')!r}",
    )
    add(
        issues,
        finite.get("preserved_negative_indices") == [0, 1, 2, 3],
        f"bad negative indices: {finite.get('preserved_negative_indices')!r}",
    )
    splice_rows = finite.get("splice_rows", [])
    add(
        issues,
        [row.get("n") for row in splice_rows] == [1241, 1242],
        f"bad splice indices: {splice_rows!r}",
    )
    for row in splice_rows:
        checks = row.get("checks", {})
        add(
            issues,
            len(checks) == 11 and all(checks.values()),
            f"failed direct checks at n={row.get('n')}: {checks!r}",
        )
        try:
            lower = Decimal(row.get("relative_Q9_margin_lower", "nan"))
        except Exception as exc:
            issues.append(f"bad relative lower at n={row.get('n')}: {exc}")
        else:
            add(
                issues,
                Decimal("0.0045") < lower < Decimal("0.0046"),
                f"unexpected relative lower at n={row.get('n')}: {lower}",
            )
    rows = artifact.get("rows", [])
    add(issues, len(rows) == 5, f"bad row count: {len(rows)}")
    add(
        issues,
        len({row.get("id") for row in rows}) == 5,
        "duplicate or missing row ids",
    )
    add(
        issues,
        all(row.get("readiness") == "ready_to_apply" for row in rows),
        "non-ready certificate row",
    )
    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover
        issues.append(f"canonical rebuild failed: {exc}")
    else:
        add(issues, artifact == canonical, f"artifact drift: {path}")
    return artifact, issues


def validate_note(path: Path) -> list[str]:
    if not path.exists():
        return [f"missing note: {path}"]
    text = path.read_text(encoding="utf-8")
    issues = []
    for marker in (
        "Status: rigorous two-index endpoint splice",
        "This is not a proof of",
        "A_1259(-100)",
        "A_1260(-100)",
        "Q_(10,n)(-100)>0 for n=1241,1242",
        "Q_(10,n)(-100)>0 for every 4<=n<=1242",
        "Q_(10,n)(-100)<0 for n=0,1,2,3",
        "n>=1243, equivalently k=n+9>=1252",
    ):
        add(issues, marker in text, f"missing note marker: {marker}")
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "order-ten entry is proved",
        "pf-infinity follows",
    ):
        add(issues, forbidden not in lowered, f"forbidden overclaim: {forbidden}")
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
            print(f"ORDER10-FINITE-SPLICE {finding}")
        print(f"order-ten finite splice: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-ten finite splice: "
        f"{summary['coefficients']} coefficients, 0 issues, "
        f"{summary['finite_splice_rows']} splice rows, "
        f"{summary['combined_positive_Q10_rows']} combined positive signs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
