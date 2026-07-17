#!/usr/bin/env python3
"""Validate the order-six lambda=-100 tail curvature reduction."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order6_m100_tail_curvature_reduction import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_artifact,
)


@dataclass(frozen=True)
class Issue:
    section: str
    code: str
    detail: str


def validate(artifact_path: Path, note_path: Path) -> tuple[dict, list[Issue]]:
    issues: list[Issue] = []
    if not artifact_path.exists():
        return {}, [Issue("artifact", "missing", str(artifact_path))]
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    if artifact.get("kind") != "jensen_window_pf_compound_order6_m100_tail_curvature_reduction":
        issues.append(Issue("artifact", "kind", str(artifact.get("kind"))))
    expected = {
        "rows": 8,
        "ready_to_apply_rows": 7,
        "exact_factorizations": 2,
        "exact_curvature_reductions": 1,
        "coefficient_positive_comparisons": 1,
        "conditional_tail_theorems": 1,
        "open_curvature_targets": 1,
    }
    for key, value in expected.items():
        if artifact.get("summary", {}).get(key) != value:
            issues.append(Issue("summary", key, str(artifact.get("summary", {}).get(key))))
    exact = artifact.get("exact", {})
    if exact.get("shifted_coefficients") != ["251", "129644", "15704684"]:
        issues.append(Issue("exact", "shifted-coefficients", str(exact.get("shifted_coefficients"))))
    if exact.get("prefactor_residual") != "0" or exact.get("stable_factor_residual") != "0":
        issues.append(Issue("exact", "residual", str(exact)))
    rows = artifact.get("rows", [])
    if sum(row.get("readiness") == "not_ready_to_apply" for row in rows) != 1:
        issues.append(Issue("rows", "open-count", str(rows)))
    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover
        issues.append(Issue("rebuild", "exception", repr(exc)))
    else:
        if artifact != canonical:
            issues.append(Issue("rebuild", "drift", str(artifact_path)))
    if not note_path.exists():
        issues.append(Issue("note", "missing", str(note_path)))
    else:
        text = note_path.read_text(encoding="utf-8")
        for marker in (
            "Status: exact endpoint-tail reduction",
            "H_(5,n)=A_(n+4)^5*exp(p(n+4))",
            "P_k<=320/k^2 for every k>=322",
            "251*m**2 + 129644*m + 15704684>0",
            "This is not a proof",
            "outputs/formal_core.md",
        ):
            if marker not in text:
                issues.append(Issue("note", "marker", marker))
    return artifact, issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate(args.artifact, args.note)
    for item in issues:
        print(f"ORDER6-TAIL {item.section} [{item.code}] {item.detail}")
    if issues:
        print(f"order-six tail curvature reduction: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-six tail curvature reduction: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['exact_factorizations']} exact factorizations, "
        f"{summary['open_curvature_targets']} open curvature target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
