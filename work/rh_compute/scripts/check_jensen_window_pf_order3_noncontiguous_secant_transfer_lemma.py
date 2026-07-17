#!/usr/bin/env python3
"""Validate the noncontiguous order-three secant transfer."""

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

from jensen_window_pf_order3_noncontiguous_secant_transfer_lemma import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    FORWARD_SOURCE,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_IDS = {
    "o3nst_01_column_normalization",
    "o3nst_02_strict_abscissas",
    "o3nst_03_orientation_identity",
    "o3nst_04_local_slope_order",
    "o3nst_05_secant_averaging",
    "o3nst_06_arbitrary_order_three",
    "o3nst_07_arbitrary_order_two",
    "o3nst_08_higher_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Order-Three Noncontiguous Secant-Transfer Lemma",
    "P_j=(u_j,v_j)=(r_(n+j),r_(n+j)*r_(n+j+1))",
    "u_(j+1)<u_j",
    "sigma_(j+1)<sigma_j for every j>=0",
    "S_(a,b)>S_(b,c)",
    "R_(3,n)(j_1,j_2,j_3)<0",
    "R_(2,n)(j_1,j_2)<0 for every j_1<j_2",
    "Order four has no automatic planar secant",
)

REQUIRED_CORE_STRINGS = (
    "### Lemma 11.29: Noncontiguous Order-Three Secant Transfer",
    "sigma_(j+1)<sigma_j",
    "R_(3,n)(j_1,j_2,j_3)<0",
    "R_(2,n)(j_1,j_2)<0",
)


@dataclass(frozen=True)
class LemmaIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> LemmaIssue:
    return LemmaIssue(section=section, issue=name, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def exact_issues() -> list[LemmaIssue]:
    findings: list[LemmaIssue] = []
    u0, u1, u2, v0, v1, v2 = sp.symbols("u0 u1 u2 v0 v1 v2")
    determinant = sp.Matrix([[1, 1, 1], [u0, u1, u2], [v0, v1, v2]]).det()
    s0 = (v1 - v0) / (u1 - u0)
    s1 = (v2 - v1) / (u2 - u1)
    expected = (u1 - u0) * (u2 - u1) * (s1 - s0)
    if sp.factor(determinant - expected) != 0:
        findings.append(issue("exact", "bad-orientation", determinant))

    points = ((5, 4), (3, -4), (2, -7), (0, -11))
    slopes = [
        sp.Rational(points[j + 1][1] - points[j][1])
        / (points[j + 1][0] - points[j][0])
        for j in range(3)
    ]
    if not slopes[0] > slopes[1] > slopes[2]:
        findings.append(issue("exact", "bad-benchmark-slopes", slopes))
    for indices in ((0, 1, 2), (1, 2, 3), (0, 1, 3), (0, 2, 3)):
        determinant = sp.Matrix(
            [[1, 1, 1], [points[j][0] for j in indices], [points[j][1] for j in indices]]
        ).det()
        if determinant >= 0:
            findings.append(issue("exact", "bad-benchmark-triple", (indices, determinant)))
    return findings


def artifact_issues(artifact: dict) -> list[LemmaIssue]:
    findings: list[LemmaIssue] = []
    if artifact.get("kind") != "jensen_window_pf_order3_noncontiguous_secant_transfer_lemma":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != (
        "exact all-column reshaped-Hankel order-two and order-three theorem "
        "at lambda=0 with order four open"
    ):
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for marker in (
        "proves",
        "does not prove",
        "orders two and three",
        "order four",
        "all-order",
        "rh",
        "lambda<=0",
    ):
        if marker not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", marker))
    for ref in artifact.get("sources", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(issue("artifact", "missing-source", ref))
    for key in ("generator", "checker"):
        ref = artifact.get(key)
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(issue("artifact", "missing-ref", (key, ref)))
    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_IDS:
        findings.append(issue("rows", "bad-ids", sorted(str(value) for value in ids)))
    expected_summary = {
        "lemma_rows": 8,
        "exact_identity_rows": 2,
        "secant_averaging_rows": 1,
        "arbitrary_order_two_rows": 1,
        "arbitrary_order_three_rows": 1,
        "open_handoffs": 1,
        "ready_to_apply_rows": 7,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(issue("summary", "mismatch", artifact.get("summary")))
    try:
        rebuilt = build_artifact()
    except Exception as exc:
        findings.append(issue("recompute", "failed", exc))
    else:
        if artifact != rebuilt:
            findings.append(issue("recompute", "mismatch", "stored artifact differs"))
    try:
        forward = load_json(FORWARD_SOURCE)
    except Exception as exc:
        findings.append(issue("source", "load-failed", exc))
    else:
        if forward.get("summary", {}).get("lambda_zero_theorem_rows") != 1:
            findings.append(issue("source", "forward-not-closed", forward.get("summary")))
    return findings


def required_string_issues(path: Path, required: tuple[str, ...], section: str) -> list[LemmaIssue]:
    if not path.exists():
        return [issue(section, "missing", path)]
    text = path.read_text(encoding="utf-8")
    return [issue(section, "missing-string", value) for value in required if value not in text]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    findings = exact_issues()
    try:
        artifact = load_json(args.artifact)
    except Exception as exc:
        artifact = {}
        findings.append(issue("artifact", "load-failed", exc))
    if artifact:
        findings.extend(artifact_issues(artifact))
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
            "validated Jensen-window PF order-three noncontiguous secant-transfer "
            f"lemma: 8 rows, {len(findings)} issues, 2 exact identities, "
            "1 secant-averaging lemma, 1 arbitrary-column order-two theorem, "
            "1 arbitrary-column order-three theorem, 1 open order-four handoff"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
