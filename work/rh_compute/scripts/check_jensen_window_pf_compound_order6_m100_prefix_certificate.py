#!/usr/bin/env python3
"""Validate the lambda=-100 signed order-six prefix certificate."""

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

from jensen_window_pf_compound_order6_m100_prefix_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    MAX_COEFFICIENT_INDEX,
    PRECISION_BITS,
    PRECISION_REPAIR,
    PRECISION_REPAIR_SUMMARY,
    PREFIX_EXTENSION,
    PREFIX_EXTENSION_SUMMARY,
    PREFIX_LAST_N,
    SOURCE_PATHS,
    build_artifact,
    sha256,
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
    if artifact.get("kind") != "jensen_window_pf_compound_order6_m100_prefix_certificate":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "through n=316" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

    expected_summary = {
        "rows": 7,
        "ready_to_apply_rows": 6,
        "coefficients": 327,
        "positive_H5_rows": 319,
        "positive_relative_H5_rows": 317,
        "positive_Q6_rows": 317,
        "inconclusive_rows": 0,
        "precision_repair_rows": 39,
        "prefix_extension_rows": 2,
        "open_analytic_tails": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    rows = artifact.get("rows", [])
    if len(rows) != 7:
        issues.append(issue("rows", "bad-count", len(rows)))
    if sum(row.get("readiness") == "ready_to_apply" for row in rows) != 6:
        issues.append(issue("rows", "bad-ready-count", rows))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or "n=317" not in open_rows[0].get("claim", ""):
        issues.append(issue("rows", "bad-open-tail", open_rows))

    finite = artifact.get("finite", {})
    if finite.get("n_range") != [0, PREFIX_LAST_N]:
        issues.append(issue("finite", "bad-n-range", finite.get("n_range")))
    if finite.get("coefficient_range") != [0, MAX_COEFFICIENT_INDEX]:
        issues.append(issue("finite", "bad-coefficient-range", finite.get("coefficient_range")))
    if finite.get("precision_bits") != PRECISION_BITS:
        issues.append(issue("finite", "bad-precision", finite.get("precision_bits")))
    finite_rows = finite.get("rows", [])
    if len(finite_rows) != PREFIX_LAST_N + 1:
        issues.append(issue("finite", "bad-row-count", len(finite_rows)))
    if [row.get("n") for row in finite_rows] != list(range(PREFIX_LAST_N + 1)):
        issues.append(issue("finite", "bad-row-indices", finite_rows[:3]))
    if not finite.get("all_H5_positive"):
        issues.append(issue("finite", "H5-not-positive", finite))
    if not finite.get("all_relative_H5_margins_positive"):
        issues.append(issue("finite", "relative-not-positive", finite))
    if not finite.get("all_Q6_positive"):
        issues.append(issue("finite", "Q6-not-positive", finite))
    if finite.get("minimum_relative_n") != 316:
        issues.append(issue("finite", "bad-minimum-index", finite.get("minimum_relative_n")))
    if Decimal(str(finite.get("minimum_relative_lower", "0"))) <= Decimal(7) / Decimal(1000):
        issues.append(issue("finite", "weak-minimum-lower", finite.get("minimum_relative_lower")))
    for row in finite_rows:
        if Decimal(str(row.get("relative_H5_margin_lower", "0"))) <= 0:
            issues.append(issue("finite", "nonpositive-relative-lower", row))
            break
        if row.get("Q6_sign") != "positive_by_signed_condensation":
            issues.append(issue("finite", "bad-Q6-sign", row))
            break

    diagnostics = artifact.get("source_diagnostics", [])
    if [row.get("source") for row in diagnostics] != [
        source.relative_to(SOURCE_PATHS[0].parents[3]).as_posix()
        for source in SOURCE_PATHS
    ]:
        issues.append(issue("sources", "bad-source-order", diagnostics))
    for source, row in zip(SOURCE_PATHS, diagnostics):
        if not source.exists():
            issues.append(issue("sources", "missing-source", source))
        elif row.get("sha256") != sha256(source):
            issues.append(issue("sources", "hash-mismatch", source))
    repair = load_json(PRECISION_REPAIR_SUMMARY)
    if repair.get("rows") != 39 or repair.get("k_min") != 191 or repair.get("k_max") != 229:
        issues.append(issue("repair", "bad-summary", repair))
    if not PRECISION_REPAIR.exists():
        issues.append(issue("repair", "missing-jsonl", PRECISION_REPAIR))
    extension = load_json(PREFIX_EXTENSION_SUMMARY)
    if (
        extension.get("rows") != 2
        or extension.get("k_min") != 325
        or extension.get("k_max") != 326
    ):
        issues.append(issue("extension", "bad-summary", extension))
    if not PREFIX_EXTENSION.exists():
        issues.append(issue("extension", "missing-jsonl", PREFIX_EXTENSION))

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
        "Status: rigorous signed order-six endpoint prefix",
        "This is not a proof",
        "K_n=H_(5,n+1)^2/(H_(5,n)*H_(5,n+2))-1",
        "A_k(-100)>0 for every 0<=k<=326",
        "A_191,...,A_229",
        "24 interval-indeterminate rows",
        "Q_(6,n)(-100)>0 for every 0<=n<=316",
        "K_316 lower=",
        ">7/1000",
        "Q_(6,n)(-100)>0 for every n>=317",
        "open analytic tail",
        "outputs/formal_core.md",
    )
    return [
        issue("note", "missing-marker", marker)
        for marker in required
        if marker not in text
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    for finding in issues:
        print(f"ORDER6-M100-PREFIX {finding.section} [{finding.code}] {finding.detail}")
    if issues:
        print(f"order-six lambda=-100 prefix: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    finite = artifact["finite"]
    print(
        "validated order-six lambda=-100 prefix: "
        f"{summary['coefficients']} coefficients, "
        f"{summary['positive_relative_H5_rows']} positive relative H5 margins, "
        f"{summary['positive_Q6_rows']} positive Q6 signs, "
        f"{summary['inconclusive_rows']} inconclusive, "
        f"{summary['precision_repair_rows']} repaired coefficients, "
        f"{summary['open_analytic_tails']} open analytic tail, "
        f"weakest n={finite['minimum_relative_n']} above 7/1000"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
