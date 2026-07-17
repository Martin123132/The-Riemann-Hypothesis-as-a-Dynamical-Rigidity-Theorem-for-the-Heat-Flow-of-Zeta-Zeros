#!/usr/bin/env python3
"""Validate the normalized deep-Schur endpoint coordinate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import re
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_endpoint_deep_schur_coordinate import (  # noqa: E402
    ARBITRARY_COLUMN_BOUND,
    ARBITRARY_MAX_ORDER,
    ARBITRARY_SHIFTS,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    ENDPOINT_BALL_SOURCE,
    INVERSE_MAX_ORDER,
    INVERSE_MAX_PART,
    RECTANGLE_MAX_ORDER,
    RECTANGLE_MAX_SHIFT,
    REPO_ROOT,
    SOURCE_PATHS,
    SYMBOLIC_ORDERS,
    build_artifact,
    sha256,
)


@dataclass(frozen=True)
class CoordinateIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> CoordinateIssue:
    return CoordinateIssue(section, code, str(detail))


def direct_partition(n: int, columns: tuple[int, ...]) -> tuple[int, ...]:
    order = len(columns)
    return tuple(n + columns[order - 1 - q] + q for q in range(order))


def direct_inverse(partition: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    order = len(partition)
    n = partition[-1] - order + 1
    columns = tuple(
        partition[order - 1 - ell] - partition[-1] + ell
        for ell in range(order)
    )
    return n, columns


def direct_coordinate_audit() -> list[CoordinateIssue]:
    issues: list[CoordinateIssue] = []
    for order in range(1, RECTANGLE_MAX_ORDER + 1):
        for n in range(RECTANGLE_MAX_SHIFT + 1):
            partition = direct_partition(n, tuple(range(order)))
            if partition != (n + order - 1,) * order:
                issues.append(issue("direct-map", "rectangle-failure", (order, n)))
                return issues

    for order in range(1, ARBITRARY_MAX_ORDER + 1):
        for n in ARBITRARY_SHIFTS:
            for columns in combinations(range(ARBITRARY_COLUMN_BOUND), order):
                partition = direct_partition(n, columns)
                if any(
                    partition[index] < partition[index + 1]
                    for index in range(order - 1)
                ):
                    issues.append(
                        issue("direct-map", "not-a-partition", (n, columns, partition))
                    )
                    return issues
                if partition[-1] < order - 1:
                    issues.append(
                        issue("direct-map", "not-deep", (n, columns, partition))
                    )
                    return issues
                for row in range(order):
                    for column in range(order):
                        left = n + column + columns[order - 1 - row]
                        right = partition[row] - row + column
                        if left != right:
                            issues.append(
                                issue(
                                    "direct-map",
                                    "entry-mismatch",
                                    (n, columns, row, column, left, right),
                                )
                            )
                            return issues
                n_back, columns_back = direct_inverse(partition)
                expected_columns = tuple(value - columns[0] for value in columns)
                if n_back != n + columns[0] or columns_back != expected_columns:
                    issues.append(
                        issue(
                            "direct-map",
                            "canonical-inverse-failure",
                            (n, columns, n_back, columns_back),
                        )
                    )
                    return issues
    return issues


def direct_endpoint_pf_audit() -> list[CoordinateIssue]:
    pattern = re.compile(r"^\[([^ ]+) \+/- ([^\]]+)\]$")
    bounds: dict[int, tuple[Fraction, Fraction]] = {}
    with ENDPOINT_BALL_SOURCE.open("r", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if row.get("lam") != "-100.0" or row.get("k") not in range(4):
                continue
            match = pattern.fullmatch(str(row.get("A_ball", "")))
            if match is None:
                return [issue("direct-pf", "bad-ball", row.get("A_ball"))]
            center = Fraction(Decimal(match.group(1)))
            radius = max(
                Fraction(Decimal(match.group(2))),
                Fraction(Decimal(str(row["A_rad"]))),
            )
            bounds[int(row["k"])] = (center - radius, center + radius)
    if set(bounds) != set(range(4)):
        return [issue("direct-pf", "missing-coefficients", sorted(bounds))]
    if any(lower <= 0 for lower, _ in bounds.values()):
        return [issue("direct-pf", "coefficient-not-positive", bounds)]

    l0, u0 = bounds[0]
    l1, u1 = bounds[1]
    l2, u2 = bounds[2]
    l3, u3 = bounds[3]
    determinant_upper = u1**3 - 2 * l0 * l1 * l2 + u0**2 * u3
    determinant_lower = l1**3 - 2 * u0 * u1 * u2 + l0**2 * l3
    issues: list[CoordinateIssue] = []
    if determinant_lower > determinant_upper:
        issues.append(
            issue("direct-pf", "reversed-determinant-bounds", determinant_lower)
        )
    if determinant_upper >= 0:
        issues.append(
            issue("direct-pf", "upper-bound-not-negative", determinant_upper)
        )
    return issues


def validate_artifact(path: Path) -> tuple[dict, list[CoordinateIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[CoordinateIssue] = []

    if artifact.get("kind") != "jensen_window_pf_endpoint_deep_schur_coordinate":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    status = str(artifact.get("status", ""))
    for marker in ("deep-Schur", "endpoint PF", "rejected"):
        if marker not in status:
            issues.append(issue("artifact", "bad-status", marker))
    boundary = str(artifact.get("proof_boundary", ""))
    for marker in (
        "rejects the deep endpoint hierarchy",
        "does not settle",
        "Jensen-window PF bridge",
        "RH",
    ):
        if marker not in boundary:
            issues.append(issue("artifact", "bad-proof-boundary", marker))

    expected_summary = {
        "rows": 14,
        "ready_to_apply_rows": 12,
        "open_endpoint_rows": 0,
        "rejected_endpoint_rows": 1,
        "separate_bridge_rows": 1,
        "symbolic_orders": 5,
        "rectangle_index_checks": 4096,
        "arbitrary_column_checks": 984,
        "inverse_partition_checks": 1023,
        "deep_cone_bijections": 1,
        "endpoint_deep_equivalences": 1,
        "rigorous_endpoint_pf_counterexamples": 1,
        "direct_literature_closing_theorems": 0,
        "completed_base_order": 9,
        "first_failed_order": 10,
        "negative_deep_rectangles": 4,
        "open_rectangle_hierarchies": 0,
        "rejected_rectangle_hierarchies": 1,
        "separate_jensen_pf_bridges": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    expected_ids = {
        f"jwpfedsc_{index:02d}_{suffix}"
        for index, suffix in (
            (1, "normalization"),
            (2, "column_reversal"),
            (3, "rectangles"),
            (4, "arbitrary_columns"),
            (5, "inverse_map"),
            (6, "deep_support_boundary"),
            (7, "endpoint_equivalence"),
            (8, "minimal_rectangle_target"),
            (9, "endpoint_pf_counterexample"),
            (10, "strictly_weaker_than_pf"),
            (11, "literature_fit"),
            (12, "dynamic_composition"),
            (13, "rejected_deep_target"),
            (14, "jensen_bridge_boundary"),
        )
    }
    rows = artifact.get("rows", [])
    row_map = {row.get("id"): row for row in rows}
    if set(row_map) != expected_ids or len(rows) != len(expected_ids):
        issues.append(issue("rows", "bad-id-set", sorted(str(key) for key in row_map)))
    if row_map.get("jwpfedsc_13_rejected_deep_target", {}).get("readiness") != "rejected_by_counterexample":
        issues.append(issue("rows", "deep-target-not-rejected", row_map.get("jwpfedsc_13_rejected_deep_target")))
    if row_map.get("jwpfedsc_14_jensen_bridge_boundary", {}).get("readiness") != "separate_open_obligation":
        issues.append(issue("rows", "bridge-not-separated", row_map.get("jwpfedsc_14_jensen_bridge_boundary")))

    expected_exact = {
        "normalization": (
            "h_k=A_k(-100)/A_0(-100) for k>=0, h_k=0 for k<0, h_0=1"
        ),
        "rectangle_identity": (
            "Q_(m,n)(-100)=A_0(-100)^m*s_((n+m-1)^m)(h)"
        ),
        "partition_map": (
            "lambda_(q+1)=n+j_(m-1-q)+q for 0<=q<m"
        ),
        "inverse_map": (
            "n=lambda_m-(m-1), j_l=lambda_(m-l)-lambda_m+l for 0<=l<m"
        ),
        "endpoint_deep_equivalence": (
            "[Q_(m,n)(-100)>0 for every m>=10,n>=0] iff "
            "[s_lambda(h)>0 for every m>=10,lambda in D_m]"
        ),
        "pf_failure": "s_(1,1,1)(h)<0 at lambda=-100",
        "deep_failure": "s_((N^10))(h)<0 for N=9,10,11,12",
        "rejected_rectangle_target": (
            "s_((N^m))(h)>0 for every m>=10 and N>=m-1 is false"
        ),
    }
    exact = artifact.get("exact", {})
    for key, expected in expected_exact.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))

    source_rows = artifact.get("sources", [])
    source_map = {row.get("path"): row for row in source_rows}
    if len(source_rows) != len(SOURCE_PATHS):
        issues.append(issue("sources", "bad-count", len(source_rows)))
    for source_path in SOURCE_PATHS:
        relative = source_path.relative_to(REPO_ROOT).as_posix()
        row = source_map.get(relative)
        if row is None:
            issues.append(issue("sources", "missing-source", relative))
        elif row.get("sha256") != sha256(source_path):
            issues.append(issue("sources", "hash-mismatch", relative))

    symbolic = artifact.get("symbolic_rectangle_audit", {})
    if symbolic.get("orders") != list(SYMBOLIC_ORDERS):
        issues.append(issue("symbolic", "bad-orders", symbolic.get("orders")))
    if symbolic.get("all_residuals_zero") is not True:
        issues.append(issue("symbolic", "stored-failure", symbolic))
    for row in symbolic.get("rows", []):
        if row.get("determinant_residual") != "0":
            issues.append(issue("symbolic", "nonzero-residual", row))

    counterexample = artifact.get("endpoint_pf_counterexample", {})
    if counterexample.get("strictly_negative") is not True:
        issues.append(issue("counterexample", "not-negative", counterexample))
    if counterexample.get("excluded_from_deep_order_three") is not True:
        issues.append(issue("counterexample", "not-excluded", counterexample))
    if counterexample.get("partition") != [1, 1, 1]:
        issues.append(issue("counterexample", "wrong-partition", counterexample))
    try:
        stored_upper = Decimal(counterexample["normalized_interval"]["upper"])
    except Exception as exc:
        issues.append(issue("counterexample", "bad-stored-upper", exc))
    else:
        if stored_upper >= 0:
            issues.append(issue("counterexample", "stored-upper-not-negative", stored_upper))

    literature = artifact.get("primary_literature_audit", [])
    expected_literature = {
        "gasca_pena_initial_minors": "applicable_transfer_only",
        "edrei_pf_classification": "strictly_too_strong",
        "pena_hankel_toeplitz_orientation": "structural_orientation_support",
        "kushel_eventual_total_positivity": "coordinate_mismatch",
    }
    literature_map = {row.get("id"): row for row in literature}
    if set(literature_map) != set(expected_literature):
        issues.append(issue("literature", "bad-id-set", sorted(literature_map)))
    for row_id, classification in expected_literature.items():
        if literature_map.get(row_id, {}).get("classification") != classification:
            issues.append(issue("literature", "bad-classification", row_id))
        if not str(literature_map.get(row_id, {}).get("url", "")).startswith("https://"):
            issues.append(issue("literature", "bad-url", row_id))

    issues.extend(direct_coordinate_audit())
    issues.extend(direct_endpoint_pf_audit())

    try:
        canonical = build_artifact()
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if canonical != artifact:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[CoordinateIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: exact normalized deep-Schur coordinate",
        "ordinary Schur theory has `h_0=1`",
        "Q_(m,n)(-100)=A_0(-100)^m*s_((n+m-1)^m)(h)",
        "lambda_(q+1)=n+j_(m-1-q)+q",
        "lambda_m-m+1>=0 exactly on D_m",
        "Within this candidate hierarchy, only the rectangles are independent",
        "s_(1,1,1)(h)=h_1^3-2*h_1*h_2+h_3",
        "genuinely additional inequality",
        "## Deep Rectangle Counterexample",
        "s_((N^10))(h)<0 for N=9,10,11,12",
        "s_((N^m))(h)>0 for every m>=10 and N>=m-1 is false",
        "antecedent is false",
        "No direct closing theorem was found in this audited primary-source",
        "not a claim about all literature",
        "separate Jensen-window obligation remains",
    )
    issues: list[CoordinateIssue] = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "rh is proved",
        "the endpoint hierarchy is proved",
        "pf-infinity is proved",
        "all relevant literature",
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
                f"ENDPOINT-DEEP-SCHUR {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"endpoint deep-Schur coordinate: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated endpoint deep-Schur coordinate: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['rectangle_index_checks']} rectangle checks, "
        f"{summary['arbitrary_column_checks']} arbitrary-column checks, "
        f"{summary['inverse_partition_checks']} inverse checks, "
        "1 rigorous endpoint PF_3 counterexample, "
        "4 negative deep rectangles, 1 rejected m>=10 rectangle hierarchy, "
        "1 separate Jensen/PF bridge"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
