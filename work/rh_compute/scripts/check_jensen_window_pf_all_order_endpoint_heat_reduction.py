#!/usr/bin/env python3
"""Validate the all-order signed-Hankel endpoint-to-heat reduction."""

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

from jensen_window_pf_all_order_endpoint_heat_reduction import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    PARITY_MAX_ORDER,
    SOURCE_PATHS,
    SYMBOLIC_ORDERS,
    build_artifact,
    sha256,
)


@dataclass(frozen=True)
class ReductionIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> ReductionIssue:
    return ReductionIssue(section, code, str(detail))


def hankel(values: tuple[sp.Symbol, ...], order: int, shift: int) -> sp.Expr:
    if order == 0:
        return sp.Integer(1)
    return sp.Matrix(
        [[values[shift + i + j] for j in range(order)] for i in range(order)]
    ).det(method="domain-ge")


def delta(expression: sp.Expr, values: tuple[sp.Symbol, ...]) -> sp.Expr:
    return sp.expand(
        sum(
            sp.diff(expression, values[index]) * values[index + 1]
            for index in range(len(values) - 1)
        )
    )


def direct_symbolic_audit() -> list[ReductionIssue]:
    issues: list[ReductionIssue] = []
    n = sp.symbols("n", integer=True, nonnegative=True)
    for order in SYMBOLIC_ORDERS:
        values = sp.symbols(f"x0:{2 * order + 2}")
        hm = hankel(values, order, 0)
        hm1 = hankel(values, order, 1)
        low0 = hankel(values, order - 1, 0)
        low1 = hankel(values, order - 1, 1)
        low2 = hankel(values, order - 1, 2)
        lowlow2 = hankel(values, order - 2, 2)
        condensation = sp.factor(hm * lowlow2 - low0 * low2 + low1**2)
        plucker = sp.factor(low1 * delta(hm, values) - low0 * hm1 - delta(low1, values) * hm)
        heat = sp.expand(
            sum(
                sp.diff(hm, values[index])
                * (4 * (n + index) + 2)
                * values[index + 1]
                for index in range(len(values) - 1)
            )
        )
        affine = sp.factor(heat - (4 * n + 8 * order - 6) * delta(hm, values))
        residuals = {
            "condensation": condensation,
            "flag-plucker": plucker,
            "affine-heat": affine,
        }
        for name, residual in residuals.items():
            if residual != 0:
                issues.append(
                    issue("symbolic", f"bad-{name}-order-{order}", residual)
                )
    return issues


def validate_artifact(path: Path) -> tuple[dict, list[ReductionIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[ReductionIssue] = []

    if artifact.get("kind") != "jensen_window_pf_all_order_endpoint_heat_reduction":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    status = str(artifact.get("status", ""))
    if "endpoint-to-heat equivalence" not in status or "order-ten" not in status:
        issues.append(issue("artifact", "bad-status", status))
    boundary = str(artifact.get("proof_boundary", ""))
    for marker in (
        "does not prove",
        "all-order positivity antecedent is false",
        "Jensen-window PF bridge",
        "RH",
    ):
        if marker not in boundary:
            issues.append(issue("artifact", "bad-proof-boundary", marker))

    expected_summary = {
        "rows": 13,
        "ready_to_apply_rows": 11,
        "open_endpoint_rows": 0,
        "rejected_endpoint_rows": 1,
        "separate_bridge_rows": 1,
        "symbolic_specialization_orders": 4,
        "orientation_parity_checks": 255,
        "all_fixed_order_tail_theorems": 1,
        "all_order_affine_identities": 1,
        "all_order_flag_plucker_identities": 1,
        "all_order_cooperative_recursions": 1,
        "single_layer_propagation_lemmas": 1,
        "endpoint_interval_equivalences": 1,
        "arbitrary_column_consequences": 1,
        "completed_base_order": 9,
        "first_failed_order": 10,
        "negative_order10_endpoint_rows": 4,
        "open_endpoint_hierarchies": 0,
        "rejected_endpoint_hierarchies": 1,
        "separate_jensen_pf_bridges": 1,
    }
    summary = artifact.get("summary", {})
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    expected_ids = {
        "jwpfaoehr_01_coordinate",
        "jwpfaoehr_02_fixed_order_tail",
        "jwpfaoehr_03_signed_condensation",
        "jwpfaoehr_04_affine_heat",
        "jwpfaoehr_05_flag_plucker",
        "jwpfaoehr_06_cooperative_flow",
        "jwpfaoehr_07_single_layer_propagation",
        "jwpfaoehr_08_completed_base",
        "jwpfaoehr_09_endpoint_equivalence",
        "jwpfaoehr_10_arbitrary_columns",
        "jwpfaoehr_11_static_hierarchy",
        "jwpfaoehr_12_rejected_endpoint",
        "jwpfaoehr_13_bridge_boundary",
    }
    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows}
    if len(rows) != len(expected_ids) or ids != expected_ids:
        issues.append(issue("rows", "bad-id-set", sorted(str(value) for value in ids)))
    row_map = {row.get("id"): row for row in rows}
    if row_map.get("jwpfaoehr_12_rejected_endpoint", {}).get("readiness") != "rejected_by_counterexample":
        issues.append(issue("rows", "endpoint-not-rejected", row_map.get("jwpfaoehr_12_rejected_endpoint")))
    if row_map.get("jwpfaoehr_13_bridge_boundary", {}).get("readiness") != "separate_open_obligation":
        issues.append(issue("rows", "bridge-not-separated", row_map.get("jwpfaoehr_13_bridge_boundary")))

    exact = artifact.get("exact", {})
    expected_exact = {
        "affine_heat": (
            "Q_(m,n)'=c_(m,n)*delta(Q_(m,n)), c_(m,n)=4*n+8*m-6"
        ),
        "signed_condensation": (
            "Q_(m,n)*Q_(m-2,n+2)=Q_(m-1,n+1)^2-"
            "Q_(m-1,n)*Q_(m-1,n+2)"
        ),
        "flag_plucker": (
            "Q_(m-1,n+1)*delta(Q_(m,n))=Q_(m-1,n)*Q_(m,n+1)+"
            "delta(Q_(m-1,n+1))*Q_(m,n)"
        ),
        "candidate_endpoint": (
            "Q_(m,n)(-100)>0 for every integer m>=10 and n>=0"
        ),
        "order10_counterexample": "Q_(10,n)(-100)<0 for n=0,1,2,3",
        "endpoint_interval_equivalence": (
            "[Q_(m,n)(-100)>0 for every m>=10,n>=0] iff "
            "[Q_(m,n)(lambda)>0 for every m>=10,n>=0,-100<=lambda<=0]"
        ),
        "bridge_boundary": (
            "all-order shifted signed-Hankel positivity does not by itself prove "
            "PF-infinity of every binomially weighted Jensen window"
        ),
    }
    for key, expected in expected_exact.items():
        if exact.get(key) != expected:
            issues.append(issue("exact", f"bad-{key}", exact.get(key)))

    source_rows = artifact.get("sources", [])
    if len(source_rows) != len(SOURCE_PATHS):
        issues.append(issue("sources", "bad-count", len(source_rows)))
    source_by_path = {row.get("path"): row for row in source_rows}
    for source_path in SOURCE_PATHS:
        relative = source_path.relative_to(source_path.parents[3]).as_posix()
        row = source_by_path.get(relative)
        if row is None:
            issues.append(issue("sources", "missing-source", relative))
        elif row.get("sha256") != sha256(source_path):
            issues.append(issue("sources", "hash-mismatch", relative))

    arithmetic = artifact.get("arithmetic_audit", {})
    if arithmetic.get("parity_checks") != PARITY_MAX_ORDER - 1:
        issues.append(issue("arithmetic", "bad-parity-count", arithmetic))
    if arithmetic.get("parity_residual_values") != [1]:
        issues.append(issue("arithmetic", "bad-parity-residual", arithmetic))
    if arithmetic.get("coefficient_residual") != "0":
        issues.append(issue("arithmetic", "bad-coefficient-residual", arithmetic))
    for order in range(2, PARITY_MAX_ORDER + 1):
        residual = (
            order * (order - 1) // 2
            + (order - 2) * (order - 3) // 2
            - 2 * ((order - 1) * (order - 2) // 2)
        )
        if residual != 1:
            issues.append(issue("arithmetic", "direct-parity-failure", order))
            break

    symbolic = artifact.get("symbolic_identity_audit", {})
    if symbolic.get("orders") != list(SYMBOLIC_ORDERS):
        issues.append(issue("symbolic", "bad-orders", symbolic.get("orders")))
    if symbolic.get("all_residuals_zero") is not True:
        issues.append(issue("symbolic", "stored-residual-failure", symbolic))
    issues.extend(direct_symbolic_audit())

    try:
        canonical = build_artifact()
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[ReductionIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: exact all-order endpoint-to-heat equivalence",
        "for every fixed m exists N_m",
        "There is no exchange of `forall m` and `exists N_m`",
        "Q_(m,n)'=c_(m,n)*delta(Q_(m,n))",
        "flag-Plucker relation",
        "Q_(m,n)(-100)>0 for every integer m>=10 and n>=0",
        "Q_(10,n)(-100)<0 for n=0,1,2,3",
        "antecedent is false",
        "binomially weighted Jensen-window problem remains open",
        "cannot assume the rejected",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "rh is proved",
        "pf-infinity follows automatically",
        "one threshold uniform in m",
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
                f"ALL-ORDER-ENDPOINT-HEAT {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"all-order endpoint-to-heat reduction: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated all-order endpoint-to-heat reduction: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['symbolic_specialization_orders']} symbolic orders, "
        f"{summary['orientation_parity_checks']} parity checks, "
        f"{summary['endpoint_interval_equivalences']} endpoint/interval equivalence, "
        f"{summary['arbitrary_column_consequences']} arbitrary-column consequence, "
        f"{summary['rejected_endpoint_hierarchies']} rejected m>=10 endpoint hierarchy, "
        f"{summary['separate_jensen_pf_bridges']} separate Jensen/PF bridge"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
