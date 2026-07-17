#!/usr/bin/env python3
"""Validate the rigorous lambda=-100 signed order-seven prefix."""

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

from jensen_window_pf_compound_order7_m100_prefix_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    MAX_COEFFICIENT_INDEX,
    PREFIX_LAST_N,
    PRECISION_BITS,
    build_artifact,
)


@dataclass(frozen=True)
class PrefixIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> PrefixIssue:
    return PrefixIssue(section, code, str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_artifact(path: Path) -> tuple[dict, list[PrefixIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = load_json(path)
    issues: list[PrefixIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order7_m100_prefix_certificate":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    status = str(artifact.get("status", ""))
    if "one open analytic tail" not in status:
        issues.append(issue("artifact", "bad-status", status))

    expected_summary = {
        "rows": 8,
        "ready_to_apply_rows": 7,
        "coefficients": MAX_COEFFICIENT_INDEX + 1,
        "positive_Q6_rows": PREFIX_LAST_N + 3,
        "positive_relative_Q6_rows": PREFIX_LAST_N + 1,
        "positive_Q7_rows": PREFIX_LAST_N + 1,
        "inconclusive_rows": 0,
        "precision_repair_rows": 12,
        "open_analytic_tails": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    exact = artifact.get("exact", {})
    required_exact = {
        "signed_condensation": (
            "Q_(7,n)*H_(5,n+2)=Q_(6,n+1)^2-Q_(6,n)*Q_(6,n+2)"
        ),
        "stable_coordinate": (
            "L_n=Q_(6,n+1)^2/(Q_(6,n)*Q_(6,n+2))-1"
        ),
        "positive_scale": (
            "Q_(7,n)=Q_(6,n)*Q_(6,n+2)*L_n/H_(5,n+2)"
        ),
        "prefix": "Q_(7,n)(-100)>0 for every 0<=n<=314",
        "remaining_tail": "Q_(7,n)(-100)>0 for every n>=315",
    }
    for key, expected in required_exact.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))

    factorizations = artifact.get("factorizations", {})
    if factorizations.get("H4", {}).get("symbolic_residual") != "0":
        issues.append(issue("factorization", "bad-H4-residual", factorizations))
    if factorizations.get("H5", {}).get("symbolic_residual") != "0":
        issues.append(issue("factorization", "bad-H5-residual", factorizations))

    finite = artifact.get("finite", {})
    if finite.get("n_range") != [0, PREFIX_LAST_N]:
        issues.append(issue("finite", "bad-n-range", finite.get("n_range")))
    if finite.get("coefficient_range") != [0, MAX_COEFFICIENT_INDEX]:
        issues.append(issue("finite", "bad-coefficient-range", finite.get("coefficient_range")))
    if finite.get("precision_bits") != PRECISION_BITS:
        issues.append(issue("finite", "bad-precision", finite.get("precision_bits")))
    rows = finite.get("rows", [])
    if len(rows) != PREFIX_LAST_N + 1:
        issues.append(issue("finite", "bad-row-count", len(rows)))
    if [row.get("n") for row in rows] != list(range(PREFIX_LAST_N + 1)):
        issues.append(issue("finite", "bad-row-indices", rows))
    if any(row.get("Q7_sign") != "positive_by_signed_condensation" for row in rows):
        issues.append(issue("finite", "bad-q7-sign", rows))
    for row in rows:
        try:
            lower = Decimal(str(row.get("relative_Q6_margin_lower")))
        except Exception:
            issues.append(issue("finite", "bad-lower", row))
            continue
        if lower <= 0:
            issues.append(issue("finite", "nonpositive-lower", row))
    if finite.get("minimum_relative_n") != PREFIX_LAST_N:
        issues.append(issue("finite", "bad-minimum-index", finite.get("minimum_relative_n")))
    try:
        minimum_lower = Decimal(str(finite.get("minimum_relative_lower")))
    except Exception:
        minimum_lower = Decimal(0)
        issues.append(issue("finite", "bad-minimum-lower", finite.get("minimum_relative_lower")))
    if minimum_lower <= Decimal("0.009"):
        issues.append(issue("finite", "weak-minimum-lower", minimum_lower))

    artifact_rows = artifact.get("rows", [])
    if len(artifact_rows) != 8 or len({row.get("id") for row in artifact_rows}) != 8:
        issues.append(issue("rows", "bad-rows", artifact_rows))
    if sum(row.get("readiness") == "ready_to_apply" for row in artifact_rows) != 7:
        issues.append(issue("rows", "bad-ready-count", artifact_rows))
    if sum(row.get("readiness") == "not_ready_to_apply" for row in artifact_rows) != 1:
        issues.append(issue("rows", "bad-open-count", artifact_rows))

    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover - diagnostic path
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
        "Status: rigorous signed order-seven endpoint prefix",
        "This is not a proof of all-shift order",
        "Q_(7,n)*H_(5,n+2)=Q_(6,n+1)^2-Q_(6,n)*Q_(6,n+2)",
        "L_n=Q_(6,n+1)^2/(Q_(6,n)*Q_(6,n+2))-1",
        "never subtracts a raw seven-by-seven determinant",
        "A_179,...,A_190",
        "no midpoint sign is promoted",
        "L_n(-100)>0 for every 0<=n<=314",
        "Q_(7,n)(-100)>0 for every 0<=n<=314",
        "L_314 lower=",
        ">9/1000",
        "Q_(7,n)(-100)>0 for every n>=315",
        "outputs/signed_hankel_jensen_bridge_target.md",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "all-shift order seven is proved",
        "pf-infinity follows",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER7-PREFIX {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-seven lambda=-100 prefix: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    finite = artifact["finite"]
    print(
        "validated order-seven lambda=-100 prefix: "
        f"{summary['coefficients']} coefficients, "
        f"{summary['positive_Q6_rows']} positive Q6 values, "
        f"{summary['positive_relative_Q6_rows']} positive relative Q6 margins, "
        f"{summary['positive_Q7_rows']} positive Q7 signs, "
        f"{summary['inconclusive_rows']} inconclusive, "
        f"{summary['precision_repair_rows']} repaired coefficients, "
        f"{summary['open_analytic_tails']} open analytic tail, "
        f"weakest n={finite['minimum_relative_n']} above 9/1000"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
