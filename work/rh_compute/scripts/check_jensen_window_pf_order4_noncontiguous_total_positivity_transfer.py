#!/usr/bin/env python3
"""Validate the arbitrary-column order-four total-positivity transfer."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_order4_noncontiguous_total_positivity_transfer import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    ORDER3_SOURCE,
    ORDER4_SOURCE,
    build_artifact,
)


@dataclass(frozen=True)
class TransferIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> TransferIssue:
    return TransferIssue(section, code, str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_exact_independently(artifact: dict) -> list[TransferIssue]:
    issues: list[TransferIssue] = []
    exact = artifact.get("exact", {})

    audits = exact.get("determinant_reversal_audits", [])
    if [row.get("order") for row in audits] != [1, 2, 3, 4]:
        issues.append(issue("exact", "bad-reversal-orders", audits))
    for order in range(1, 5):
        entries = sp.symbols(f"x0:{order * order}")
        matrix = sp.Matrix(order, order, entries)
        orientation = (-1) ** (order * (order - 1) // 2)
        residual = sp.expand(matrix[:, ::-1].det() - orientation * matrix.det())
        if residual != 0:
            issues.append(issue("exact", "reversal-identity", order))

    index_mapping = exact.get("index_mapping", {})
    if index_mapping.get("solid_blocks_checked") != 240:
        issues.append(
            issue(
                "exact",
                "bad-index-map-count",
                index_mapping.get("solid_blocks_checked"),
            )
        )
    for n in range(3):
        N = 8
        for order in range(1, 5):
            for row_start in range(5 - order):
                for column_start in range(N - order + 2):
                    shift = n + row_start + N - column_start - order + 1
                    for i in range(order):
                        for j in range(order):
                            left = n + row_start + i + N - column_start - j
                            right = shift + i + order - 1 - j
                            if left != right:
                                issues.append(
                                    issue(
                                        "exact",
                                        "index-map",
                                        (n, order, row_start, column_start, i, j),
                                    )
                                )

    benchmark = exact.get("reciprocal_factorial_benchmark", {})
    if benchmark.get("strict_signed_minors") != 1020:
        issues.append(
            issue(
                "benchmark",
                "bad-signed-minor-count",
                benchmark.get("strict_signed_minors"),
            )
        )
    if benchmark.get("positive_initial_minors_for_4_by_9_block") != 36:
        issues.append(
            issue(
                "benchmark",
                "bad-initial-minor-count",
                benchmark.get("positive_initial_minors_for_4_by_9_block"),
            )
        )
    minimum = sp.Rational(str(benchmark.get("minimum_signed_minor", "0")))
    if minimum <= 0:
        issues.append(issue("benchmark", "nonpositive-minimum", minimum))

    countermodel = exact.get("nonpromotion_countermodel", {})
    expected_countermodel = {
        "H4_n0": sp.Integer(288076),
        "H4_n1": sp.Integer(264875),
        "R4_columns_0_1_3_4": sp.Integer(-231169),
        "H2_n0_wrong_sign": sp.Integer(209),
    }
    for key, expected in expected_countermodel.items():
        value = sp.Integer(str(countermodel.get(key, "0")))
        if value != expected:
            issues.append(issue("countermodel", f"bad-{key}", value))

    return issues


def validate_artifact(path: Path) -> tuple[dict, list[TransferIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = load_json(path)
    issues: list[TransferIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_order4_noncontiguous_total_positivity_transfer"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != (
        "exact arbitrary-column reshaped-Hankel order-four theorem at lambda zero"
    ):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))

    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 11,
        "ready_to_apply_rows": 10,
        "determinant_reversal_orders": 4,
        "index_mapping_blocks": 240,
        "reciprocal_factorial_signed_minors": 1020,
        "published_theorem_rows": 1,
        "arbitrary_order_four_theorems": 1,
        "fixed_order_transfer_theorems": 1,
        "nonpromotion_guards": 1,
        "open_handoffs": 1,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    rows = artifact.get("rows", [])
    if len(rows) != 11:
        issues.append(issue("rows", "bad-count", len(rows)))
    ids = [row.get("id") for row in rows]
    if len(ids) != len(set(ids)):
        issues.append(issue("rows", "duplicate-id", ids))
    if sum(row.get("readiness") == "ready_to_apply" for row in rows) != 10:
        issues.append(issue("rows", "bad-ready-count", rows))
    conclusions = [row for row in rows if row.get("role") == "theorem_conclusion"]
    if len(conclusions) != 1 or "R_(4,n)" not in conclusions[0].get("formula", ""):
        issues.append(issue("rows", "missing-order-four-conclusion", conclusions))
    structural = [
        row for row in rows if row.get("role") == "exact_structural_theorem"
    ]
    if len(structural) != 1 or "1<=k<=m" not in structural[0].get("formula", ""):
        issues.append(issue("rows", "missing-fixed-order-transfer", structural))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or "order five" not in open_rows[0].get("claim", ""):
        issues.append(issue("rows", "bad-open-handoff", open_rows))

    for source in (ORDER3_SOURCE, ORDER4_SOURCE):
        if not source.exists():
            issues.append(issue("sources", "missing-source", source))

    issues.extend(validate_exact_independently(artifact))

    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover - diagnostic path
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[TransferIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: exact arbitrary-column reshaped-Hankel order-four theorem",
        "This is not a proof of contiguous order five",
        "Gasca and Pena",
        "all its initial minors are strictly positive",
        "B_(i,q)^(n,N)=A_(n+i+N-q)(0)",
        "det B[r:r+k,c:c+k]=epsilon_k*H_(k,n+r+N-c-k+1)(0)",
        "R_(4,n)(j_1,j_2,j_3,j_4)>0",
        "for every n>=0 and 0<=j_1<j_2<j_3<j_4 at lambda=0",
        "The first new layer is contiguous",
        "order five",
        "H_(4,0)=288076>0",
        "R_(4,0)(0,1,3,4)=-231169<0",
        "contiguous order four alone cannot be promoted",
        "outputs/jensen_window_pf_compound_order4_uniform_heat_forward_invariance_certificate.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
    )
    issues: list[TransferIssue] = []
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
                f"ORDER4-NONCONTIGUOUS-TP {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(
            "order-four noncontiguous total-positivity transfer: "
            f"{len(issues)} issues"
        )
        return 1

    summary = artifact["summary"]
    print(
        "validated order-four noncontiguous total-positivity transfer: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['determinant_reversal_orders']} reversal orders, "
        f"{summary['index_mapping_blocks']} solid-block maps, "
        f"{summary['reciprocal_factorial_signed_minors']} signed benchmark minors, "
        f"{summary['arbitrary_order_four_theorems']} arbitrary-column order-four theorem, "
        f"{summary['fixed_order_transfer_theorems']} fixed-order transfer theorem, "
        f"{summary['open_handoffs']} open order-five handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
