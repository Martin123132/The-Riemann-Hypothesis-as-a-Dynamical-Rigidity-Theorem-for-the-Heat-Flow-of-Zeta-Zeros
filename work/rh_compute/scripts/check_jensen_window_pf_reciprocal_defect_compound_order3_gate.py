#!/usr/bin/env python3
"""Validate the reciprocal-defect compound order-three gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_reciprocal_defect_compound_order3_gate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
    build_countermodel,
)


REQUIRED_IDS = {
    "rdco3_01_determinant_coordinate",
    "rdco3_02_defect_coordinate",
    "rdco3_03_reciprocal_coordinate",
    "rdco3_04_increment_budget",
    "rdco3_05_sufficient_increment",
    "rdco3_06_boundary_family",
    "rdco3_07_strict_countermodel",
    "rdco3_08_m100_entry",
    "rdco3_09_forward_invariance",
    "rdco3_10_noncontiguous_transfer",
    "rdco3_11_live_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Reciprocal-Defect Compound Order-Three Gate",
    "D_(3,n)=A_n^3*rho_n^6*x_(n+1)^3*F_n",
    "F_n=d_(n+1)*d_(n+3)*x_(n+2)^2-d_(n+2)^2",
    "C_n:=q_(n+1)*q_(n+3)-q_(n+2)^2+1>0",
    "C_n=1-a_n*b_n+q_(n+2)*(b_n-a_n)",
    "q_1,q_2,q_3=['10', '109/10', '58/5']",
    "C_0=-181/100",
    "strictly hyperbolic neighboring cubic Jensen windows",
    "C_n(-100)>0 and D_(3,n)(-100)<0 for every n>=0",
    "uniform analytic tail lower bound=57613471/66107054971",
    "C_n(lambda)>0 for every n>=0 and finite lambda>=-100",
    "D_(3,n)(0)<0 for every n>=0",
    "R_(3,n)(j_1,j_2,j_3)<0 for every n>=0",
    "All-column reshaped-Hankel order three is closed",
)

REQUIRED_CORE_STRINGS = (
    "### Lemma 11.26: Reciprocal-Defect Compound Order Three",
    "C_n=q_(n+1)*q_(n+3)-q_(n+2)^2+1",
    "q_1=10, q_2=109/10, q_3=58/5",
    "det[A_(i+j)]_(i,j=0..2)>0",
)


@dataclass(frozen=True)
class GateIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> GateIssue:
    return GateIssue(section=section, issue=name, detail=str(detail))


def symbolic_issues() -> list[GateIssue]:
    findings: list[GateIssue] = []
    a0, ratio, x1, x2, x3 = sp.symbols("a0 ratio x1 x2 x3", positive=True)
    ratios = (ratio, ratio * x1, ratio * x1 * x2, ratio * x1 * x2 * x3)
    coefficients = [a0]
    for item in ratios:
        coefficients.append(sp.expand(coefficients[-1] * item))
    determinant = sp.factor(
        sp.Matrix(
            [[coefficients[i + j] for j in range(3)] for i in range(3)]
        ).det()
    )
    frontier = x1 * x2**2 * x3 - x1 * x2**2 - x2**2 * x3 + 2 * x2 - 1
    expected = a0**3 * ratio**6 * x1**3 * frontier
    if sp.factor(determinant - expected) != 0:
        findings.append(issue("symbolic", "bad-determinant-coordinate", determinant))

    d1, d2, d3 = sp.symbols("d1 d2 d3", positive=True)
    defect = sp.factor(frontier.subs({x1: 1 - d1, x2: 1 - d2, x3: 1 - d3}))
    if sp.expand(defect - (d1 * d3 * (1 - d2) ** 2 - d2**2)) != 0:
        findings.append(issue("symbolic", "bad-defect-coordinate", defect))

    q1, q2, q3 = sp.symbols("q1 q2 q3", positive=True)
    margin = q1 * q3 - q2**2 + 1
    compound_gap = sp.factor(
        (d2**2 - (1 - d2) ** 2 * d1 * d3).subs(
            {d1: q1**-2, d2: q2**-2, d3: q3**-2}
        )
    )
    positive_factorization = (
        margin * (q1 * q3 + q2**2 - 1) / (q1**2 * q2**4 * q3**2)
    )
    if sp.factor(compound_gap - positive_factorization) != 0:
        findings.append(issue("symbolic", "bad-reciprocal-factorization", compound_gap))
    a, b = sp.symbols("a b", nonnegative=True)
    increment = sp.expand(
        margin.subs({q1: q2 - a, q3: q2 + b}, simultaneous=True)
    )
    if sp.expand(increment - (1 - a * b + q2 * (b - a))) != 0:
        findings.append(issue("symbolic", "bad-increment-coordinate", increment))

    countermodel = build_countermodel()
    if countermodel["compound_margin"] != "-181/100":
        findings.append(issue("countermodel", "bad-margin", countermodel["compound_margin"]))
    if countermodel["hankel_determinant"] != "4106267526339/1899424214416000000":
        findings.append(issue("countermodel", "bad-determinant", countermodel["hankel_determinant"]))
    return findings


def validate_artifact(artifact: dict) -> list[GateIssue]:
    findings: list[GateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_reciprocal_defect_compound_order3_gate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("date") != "2026-07-12":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    if artifact.get("status") != (
        "complete all-column order-three theorem with strict abstract "
        "countermodel gate"
    ):
        findings.append(issue("artifact", "bad-status", artifact.get("status")))

    boundary = str(artifact.get("proof_boundary", "")).lower()
    for phrase in (
        "factors",
        "proves",
        "countermodel",
        "entry",
        "does not prove",
        "order four",
        "order four",
        "all-order",
        "jensen bridge",
        "rh",
        "lambda<=0",
    ):
        if phrase not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", phrase))

    for ref in artifact.get("sources", []):
        if not (REPO_ROOT / ref).exists():
            findings.append(issue("artifact", "missing-source", ref))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows}
    if len(rows) != 11:
        findings.append(issue("rows", "bad-count", len(rows)))
    if ids != REQUIRED_IDS:
        findings.append(issue("rows", "bad-ids", sorted(str(value) for value in ids)))
    expected_roles = {
        "exact_identity": 2,
        "exact_equivalence": 2,
        "exact_theorem": 1,
        "exact_benchmark": 1,
        "exact_countermodel": 1,
        "interval_analytic_theorem": 1,
        "theorem_composition": 2,
        "open_handoff": 1,
    }
    actual_roles = {
        role: sum(row.get("role") == role for row in rows)
        for role in expected_roles
    }
    if actual_roles != expected_roles:
        findings.append(issue("rows", "bad-role-counts", actual_roles))
    if artifact != build_artifact():
        findings.append(issue("recompute", "mismatch", "stored artifact differs from exact rebuild"))
    return findings


def required_string_issues(path: Path, strings: tuple[str, ...], section: str) -> list[GateIssue]:
    text = path.read_text(encoding="utf-8")
    return [issue(section, "missing-string", value) for value in strings if value not in text]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    findings = symbolic_issues()
    findings.extend(validate_artifact(artifact))
    findings.extend(required_string_issues(args.note, REQUIRED_NOTE_STRINGS, "note"))
    findings.extend(
        required_string_issues(
            REPO_ROOT / "outputs/formal_core.md", REQUIRED_CORE_STRINGS, "formal-core"
        )
    )
    ok = not findings
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in findings]}, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF reciprocal-defect compound order-three gate: "
            f"11 rows, {len(findings)} issues, 2 exact coordinate identities, "
            "2 exact sign equivalences, 1 sufficient increment theorem, "
            "1 exact boundary benchmark, 1 strict cone countermodel, "
            "1 all-shift lambda=-100 entry theorem, 1 full forward propagation "
            "theorem, 1 arbitrary-column order-three theorem, 1 open order-four "
            "handoff"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
