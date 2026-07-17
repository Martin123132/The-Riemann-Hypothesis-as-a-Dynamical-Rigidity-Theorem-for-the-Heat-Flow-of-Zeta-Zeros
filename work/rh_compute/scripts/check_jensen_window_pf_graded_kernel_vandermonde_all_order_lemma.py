#!/usr/bin/env python3
"""Validate the all-fixed-order graded-kernel Vandermonde lemma."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_graded_kernel_vandermonde_all_order_lemma import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    MAX_COEFFICIENT_DEGREE,
    MAX_RECORDED_ORDER,
    MAX_STRESS_ORDER,
    build_artifact,
    signed_leading_constant,
)


@dataclass(frozen=True)
class LemmaIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> LemmaIssue:
    return LemmaIssue(section, code, str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_artifact(path: Path) -> tuple[dict, list[LemmaIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = load_json(path)
    issues: list[LemmaIssue] = []
    if artifact.get("kind") != "jensen_window_pf_graded_kernel_vandermonde_all_order_lemma":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    status = str(artifact.get("status", ""))
    if "no all-shift theorem" not in status:
        issues.append(issue("artifact", "bad-status", status))

    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 12,
        "ready_to_apply_rows": 12,
        "recorded_orders": MAX_RECORDED_ORDER,
        "coefficient_valuation_cells": (MAX_COEFFICIENT_DEGREE + 1) ** 2,
        "permutation_stress_cases": sum(
            __import__("math").factorial(order)
            for order in range(1, MAX_STRESS_ORDER + 1)
        ),
        "all_fixed_order_leading_theorems": 1,
        "all_fixed_order_eventual_tail_theorems": 1,
        "effective_thresholds": 0,
        "all_shift_theorems": 0,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))

    order_rows = artifact.get("order_rows", [])
    if [row.get("order") for row in order_rows] != list(
        range(1, MAX_RECORDED_ORDER + 1)
    ):
        issues.append(issue("orders", "bad-order-range", order_rows))
    for row in order_rows:
        order = row.get("order")
        if not isinstance(order, int):
            continue
        degree = order * (order - 1) // 2
        raw = signed_leading_constant(order)
        if row.get("leading_degree") != degree:
            issues.append(issue("orders", "bad-degree", row))
        if row.get("epsilon") != (-1) ** degree:
            issues.append(issue("orders", "bad-epsilon", row))
        if row.get("raw_leading_constant") != raw:
            issues.append(issue("orders", "bad-raw-constant", row))
        if row.get("positive_constant") != abs(raw):
            issues.append(issue("orders", "bad-positive-constant", row))

    if len(order_rows) >= 7:
        order7 = order_rows[6]
        if order7.get("leading_degree") != 21:
            issues.append(issue("order7", "bad-degree", order7))
        if order7.get("raw_leading_constant") != -52183852646400:
            issues.append(issue("order7", "bad-leading-constant", order7))

    valuation = artifact.get("valuation_lemma", {})
    table = valuation.get("coefficient_table", [])
    if len(table) != (MAX_COEFFICIENT_DEGREE + 1) ** 2:
        issues.append(issue("valuation", "bad-table-count", len(table)))
    for row in table:
        p = row.get("p")
        q = row.get("q")
        if not isinstance(p, int) or not isinstance(q, int):
            issues.append(issue("valuation", "bad-index", row))
            continue
        if p == q == 0:
            expected = 0
        elif p == 0 or q == 0:
            expected = "infinity"
        else:
            expected = max(p, q)
        if row.get("h_valuation_lower_bound") != expected:
            issues.append(issue("valuation", "bad-bound", row))
        if p == q and p > 0:
            expected_diagonal = f"(2*G2)**{p}/{__import__('math').factorial(p)}"
            if row.get("diagonal_leading_coefficient") != expected_diagonal:
                issues.append(issue("valuation", "bad-diagonal", row))

    stress = artifact.get("permutation_stress", {})
    stress_rows = stress.get("rows", [])
    if [row.get("order") for row in stress_rows] != list(
        range(1, MAX_STRESS_ORDER + 1)
    ):
        issues.append(issue("stress", "bad-order-range", stress_rows))
    if any(row.get("zero_penalties") != 1 for row in stress_rows):
        issues.append(issue("stress", "nonunique-zero", stress_rows))

    theorem = artifact.get("vandermonde_theorem", {})
    required_theorem = {
        "degree": "D=binom(m,2)",
        "raw_first_coefficient": (
            "epsilon_m*2^D*prod_(j=1)^(m-1)j!*G_2^D"
        ),
        "signed_first_coefficient": (
            "2^D*prod_(j=1)^(m-1)j!*G_2^D>0"
        ),
    }
    for key, expected in required_theorem.items():
        if theorem.get(key) != expected:
            issues.append(issue("theorem", f"bad-{key}", theorem.get(key)))

    heat = artifact.get("heat_tilt_all_fixed_order_extension", {})
    required_heat = {
        "lambert_recurrence": (
            "F_0(w)=w^2; F_(s+1)(w)=-s*F_s(w)+"
            "w/(1+w)*F_s'(w)"
        ),
        "lambert_derivative": (
            "d^s(w^2)/dk^s=k^(-s)*F_s(w), F_s(w)=O_s(w)"
        ),
        "remainder_choice": (
            "for fixed s choose R>s; 2^s O_R(w^(3R)/k^R)="
            "o(w/k^s)"
        ),
        "theorem": (
            "Delta^s log R_T^(1)(k)=O_s(log(k)/k^s) for every "
            "fixed s>=2, uniformly for 0<=T<=100"
        ),
    }
    for key, expected in required_heat.items():
        if heat.get(key) != expected:
            issues.append(issue("heat", f"bad-{key}", heat.get(key)))

    rows = artifact.get("rows", [])
    if len(rows) != 12 or len({row.get("id") for row in rows}) != 12:
        issues.append(issue("rows", "bad-rows", rows))
    if any(row.get("readiness") != "ready_to_apply" for row in rows):
        issues.append(issue("rows", "bad-readiness", rows))

    try:
        canonical = build_artifact()
    except Exception as exc:  # pragma: no cover - diagnostic path
        issues.append(issue("rebuild", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[LemmaIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: exact all-fixed-order leading-determinant theorem",
        "This is not a proof of",
        "p+q-ell>=max(p,q)",
        "[h^k]c_(k,k)=(2*G_2)^k/k!",
        "sum_i max(i,pi(i))=D+(1/2)*sum_i|i-pi(i)|",
        "epsilon_m*2^D*prod_(j=1)^(m-1)j!*G_2^D",
        "F_(s+1)(w)=-s*F_s(w)+w/(1+w)*F_s'(w)",
        "2^s O_R(w^(3R)/k^R)=o(w/k^s)",
        "for every fixed s>=2, uniformly for 0<=T<=100",
        "for every fixed m, exists N_m",
        "-52183852646400*G_2^21",
        "threshold may depend on",
        "sign-regularity bridge remain separate obligations",
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
        "all-shift sign regularity follows",
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
                f"GRADED-KERNEL-VANDERMONDE {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"graded-kernel all-order Vandermonde lemma: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated graded-kernel all-order Vandermonde lemma: "
        f"{summary['recorded_orders']} order rows, "
        f"{summary['permutation_stress_cases']} permutation stress cases, "
        f"{summary['coefficient_valuation_cells']} coefficient-valuation cells, "
        f"{summary['all_fixed_order_eventual_tail_theorems']} all-fixed-order "
        "eventual signed-tail theorem, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
