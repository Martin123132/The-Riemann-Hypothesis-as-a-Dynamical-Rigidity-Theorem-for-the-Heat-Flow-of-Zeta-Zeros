#!/usr/bin/env python3
"""Validate all-shift contiguous order-three entry at lambda=-100."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate import (  # noqa: E402
    ADAPTIVE_SOURCE,
    ANCHOR,
    COLLAR_321,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    LOCAL_REPAIR,
    REPO_ROOT,
    build_artifact,
    exact_tail_data,
)
import jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate as m100  # noqa: E402


REQUIRED_IDS = {
    "m100co3_01_repaired_source_merge",
    "m100co3_02_compound_coordinate",
    "m100co3_03_prefix_compound",
    "m100co3_04_prefix_defect_gap",
    "m100co3_05_prefix_shape",
    "m100co3_06_scaled_defect_anchor",
    "m100co3_07_adaptive_defect_input",
    "m100co3_08_tail_increment_bound",
    "m100co3_09_tail_compound",
    "m100co3_10_full_entry",
    "m100co3_11_forward_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda -100 Compound Order-Three Entry Certificate",
    "Status: all-shift contiguous order-three entry theorem at lambda=-100",
    "C_n=q_(n+1)*q_(n+3)-q_(n+2)^2+1>0, 0<=n<=317",
    "s_319-251/500>",
    "1-sqrt(1-2/m)<1/m+1/m^2",
    "m^2-2*m-1>0",
    "C_n>57613471/66107054971>0, n>=318",
    "D_(3,n)(-100)=det[A_(n+i+j)]_(i,j=0..2)<0",
    "Forward propagation to `lambda=0`",
)

REQUIRED_CORE_STRINGS = (
    "### Lemma 11.27: Contiguous Order-Three Entry at Lambda=-100",
    "s_319>251/500",
    "57613471/66107054971",
    "D_(3,n)(-100)<0",
)


@dataclass(frozen=True)
class CertificateIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> CertificateIssue:
    return CertificateIssue(section=section, issue=name, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def exact_issues() -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    q, u, v = sp.symbols("q u v", positive=True)
    q_prev = q - u
    q_next = q + v
    compound = sp.expand(q_prev * q_next - q**2 + 1)
    if sp.expand(compound - (1 - u * v + q * (v - u))) != 0:
        findings.append(issue("exact", "bad-increment-identity", compound))

    m = sp.symbols("m", positive=True)
    rhs_square_gap = sp.factor(
        (1 - 2 / m) - (1 - 1 / m - 1 / m**2) ** 2
    )
    expected_gap = sp.factor((m**2 - 2 * m - 1) / m**4)
    if sp.factor(rhs_square_gap - expected_gap) != 0:
        findings.append(issue("exact", "bad-square-root-reduction", rhs_square_gap))

    tail = exact_tail_data()
    if tail.get("uniform_tail_lower") != "57613471/66107054971":
        findings.append(issue("exact", "bad-tail-margin", tail))
    bracket = Fraction(1) + Fraction(2, 641) + Fraction(2, 641**2) + Fraction(1, 641**3)
    lower = Fraction(1) - Fraction(250, 251) * bracket
    if lower != Fraction(57_613_471, 66_107_054_971) or lower <= 0:
        findings.append(issue("exact", "bad-endpoint-comparison", lower))
    if ANCHOR != Fraction(251, 500):
        findings.append(issue("exact", "bad-anchor", ANCHOR))
    return findings


def source_issues() -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    for path, expected in ((LOCAL_REPAIR, (208, 217, 10)), (COLLAR_321, (321, 321, 1))):
        try:
            values = m100.load_source(path)
        except Exception as exc:
            findings.append(issue("sources", "load-failed", (path, exc)))
            continue
        actual = (min(values), max(values), len(values))
        if actual != expected:
            findings.append(issue("sources", "bad-range", (path, actual, expected)))
    try:
        adaptive = load_json(ADAPTIVE_SOURCE)
    except Exception as exc:
        findings.append(issue("sources", "adaptive-load-failed", exc))
    else:
        summary = adaptive.get("summary", {})
        if summary.get("open_requirements") != 0 or summary.get("defect_conclusion_rows") != 4:
            findings.append(issue("sources", "adaptive-not-closed", summary))
    return findings


def artifact_issues(artifact: dict) -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("date") != "2026-07-12":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    if artifact.get("status") != (
        "all-shift contiguous order-three entry theorem at lambda=-100 with open forward propagation"
    ):
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for marker in (
        "proves",
        "does not prove",
        "contiguous",
        "forward",
        "noncontiguous",
        "higher orders",
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
    if len(rows) != 11:
        findings.append(issue("rows", "bad-count", len(rows)))
    expected_summary = {
        "certificate_rows": 11,
        "positive_coefficients": 322,
        "prefix_compound_rows": 318,
        "prefix_defect_gap_rows": 318,
        "prefix_decreasing_increment_rows": 318,
        "prefix_increasing_compound_rows": 317,
        "tail_start_n": 318,
        "tail_theorem_rows": 1,
        "full_entry_rows": 1,
        "open_forward_handoffs": 1,
        "ready_to_apply_rows": 5,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(issue("summary", "mismatch", artifact.get("summary")))

    finite = artifact.get("finite", {})
    expected_finite = {
        "coefficient_range": [0, 321],
        "positive_coefficients": 322,
        "contraction_range": [1, 320],
        "positive_defects": 320,
        "increment_range": [1, 319],
        "positive_subunit_increments": 319,
        "decreasing_increment_range": [1, 318],
        "positive_increment_decreases": 318,
        "compound_shift_range": [0, 317],
        "positive_compound_margins": 318,
        "positive_defect_gaps": 318,
        "increasing_compound_range": [0, 316],
        "positive_compound_increases": 317,
        "minimum_compound_at_n": 0,
        "minimum_defect_gap_at_n": 317,
    }
    for key, value in expected_finite.items():
        if finite.get(key) != value:
            findings.append(issue("finite", f"bad-{key}", finite.get(key)))
    try:
        if not Decimal(finite["minimum_compound_lower"]) > Decimal("0.7707"):
            findings.append(issue("finite", "weak-compound-margin", finite["minimum_compound_lower"]))
        if not Decimal(finite["minimum_defect_gap_lower"]) > Decimal("7.7e-9"):
            findings.append(issue("finite", "weak-defect-gap", finite["minimum_defect_gap_lower"]))
        anchor_margin = Decimal(finite["scaled_defect_anchor"]["margin_lower"])
        if not Decimal("0.00075") < anchor_margin < Decimal("0.00077"):
            findings.append(issue("finite", "bad-anchor-margin", anchor_margin))
    except Exception as exc:
        findings.append(issue("finite", "bad-decimal", exc))

    try:
        rebuilt = build_artifact()
    except Exception as exc:
        findings.append(issue("recompute", "failed", exc))
    else:
        if artifact != rebuilt:
            findings.append(issue("recompute", "mismatch", "stored artifact differs from rebuild"))
    return findings


def required_string_issues(path: Path, required: tuple[str, ...], section: str) -> list[CertificateIssue]:
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
    findings.extend(source_issues())
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
            "validated Jensen-window PF negative-lambda -100 compound order-three "
            f"entry certificate: 11 rows, {len(findings)} issues, 322 positive "
            "coefficients, 318 prefix compound margins, 318 prefix defect gaps, "
            "1 exact tail theorem, 1 all-shift entry theorem, 1 open forward handoff"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
