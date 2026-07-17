#!/usr/bin/env python3
"""Validate the contiguous order-four condensation gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order4_condensation_gate import (  # noqa: E402
    COLLAR_321,
    COLLAR_322,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    LOCAL_REPAIR,
    REPO_ROOT,
    build_artifact,
    exact_countermodel,
)
import jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate as m100  # noqa: E402


REQUIRED_IDS = {
    "co4cg_01_condensation",
    "co4cg_02_log_concavity_equivalence",
    "co4cg_03_gap_coordinate",
    "co4cg_04_m100_prefix",
    "co4cg_05_lower_cone_countermodel",
    "co4cg_06_wrong_h4",
    "co4cg_07_tail_target",
    "co4cg_08_flow_and_columns",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Compound Order-Four Condensation Gate",
    "H_(4,n)*H_(2,n+2)=H_(3,n)*H_(3,n+2)-H_(3,n+1)^2",
    "H_(4,n)>0 iff T_(n+1)^2>T_n*T_(n+2)",
    "G_(n+1)^2>x_(n+3)^3*G_n*G_(n+2)",
    "all `317` available",
    "P_n=log(G_n*G_(n+2)/G_(n+1)^2)<2/(n+3)^2",
    "P_n<=4/(n+3)^2, n>=317",
    "76466200>0",
    "s_1,...,s_5=['3/10', '2/5', '41/100', '49/100', '11/20']",
    "G_1^2-x_3^3*G_0*G_2=",
    "H_(4,0)=",
    "complete lower compound layers do not imply order four",
)

REQUIRED_CORE_STRINGS = (
    "### Lemma 11.30: Contiguous Order-Four Condensation Frontier",
    "H_(4,n)*H_(2,n+2)",
    "G_(n+1)^2>x_(n+3)^3*G_n*G_(n+2)",
    "s_1=3/10, s_2=2/5, s_3=41/100, s_4=49/100, s_5=11/20",
)


@dataclass(frozen=True)
class GateIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> GateIssue:
    return GateIssue(section=section, issue=name, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def exact_issues() -> list[GateIssue]:
    findings: list[GateIssue] = []
    symbols = sp.symbols("a0:7")
    matrix = sp.Matrix([[symbols[i + j] for j in range(4)] for i in range(4)])
    h4 = matrix.det()
    h2_center = sp.Matrix([[symbols[2 + i + j] for j in range(2)] for i in range(2)]).det()
    h3_left = sp.Matrix([[symbols[i + j] for j in range(3)] for i in range(3)]).det()
    h3_mid = sp.Matrix([[symbols[1 + i + j] for j in range(3)] for i in range(3)]).det()
    h3_right = sp.Matrix([[symbols[2 + i + j] for j in range(3)] for i in range(3)]).det()
    if sp.factor(h4 * h2_center - h3_left * h3_right + h3_mid**2) != 0:
        findings.append(issue("exact", "bad-condensation", "Desnanot-Jacobi mismatch"))
    if 753 * 320**2 - 1000 * (2 * 320 + 1) != 76_466_200:
        findings.append(issue("exact", "bad-rational-endpoint", "k=320"))

    countermodel = exact_countermodel()
    if countermodel.get("h4_determinant") != "-6608596712914764288/582076609134674072265625":
        findings.append(issue("countermodel", "bad-h4", countermodel.get("h4_determinant")))
    if countermodel.get("order_four_gap_frontier") != "-933340447356927/58618164062500000000":
        findings.append(issue("countermodel", "bad-frontier", countermodel.get("order_four_gap_frontier")))
    if not all(countermodel.get("checks", {}).values()):
        findings.append(issue("countermodel", "failed-check", countermodel.get("checks")))
    return findings


def source_issues() -> list[GateIssue]:
    findings: list[GateIssue] = []
    for path, expected in (
        (LOCAL_REPAIR, (208, 217, 10)),
        (COLLAR_321, (321, 321, 1)),
        (COLLAR_322, (322, 322, 1)),
    ):
        try:
            values = m100.load_source(path)
        except Exception as exc:
            findings.append(issue("source", "load-failed", exc))
            continue
        actual = (min(values), max(values), len(values))
        if actual != expected:
            findings.append(issue("source", "bad-range", (actual, expected)))
    return findings


def artifact_issues(artifact: dict) -> list[GateIssue]:
    findings: list[GateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order4_condensation_gate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != (
        "exact contiguous order-four condensation coordinate with 317-row "
        "lambda=-100 prefix and strict lower-order countermodel"
    ):
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for marker in (
        "derives",
        "certifies",
        "countermodel",
        "does not prove",
        "all-index",
        "forward",
        "arbitrary",
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
        "gate_rows": 8,
        "exact_identity_rows": 2,
        "exact_equivalence_rows": 1,
        "prefix_order_four_margins": 317,
        "prefix_scaled_penalty_caps": 317,
        "strict_lower_order_countermodels": 1,
        "forbidden_promotions": 1,
        "open_handoffs": 2,
        "ready_to_apply_rows": 3,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(issue("summary", "mismatch", artifact.get("summary")))
    finite = artifact.get("finite", {})
    for key, value in {
        "coefficient_range": [0, 322],
        "positive_coefficients": 323,
        "order_three_gap_range": [0, 318],
        "positive_order_three_gaps": 319,
        "order_four_margin_range": [0, 316],
        "positive_order_four_margins": 317,
        "positive_scaled_penalty_caps": 317,
        "minimum_margin_at_n": 316,
        "maximum_margin_at_n": 0,
        "minimum_scaled_penalty_cap_at_n": 316,
        "maximum_scaled_penalty_at_n": 316,
    }.items():
        if finite.get(key) != value:
            findings.append(issue("finite", f"bad-{key}", finite.get(key)))
    try:
        minimum = Decimal(finite["minimum_margin_lower"])
        if not Decimal("0.0047") < minimum < Decimal("0.0048"):
            findings.append(issue("finite", "bad-minimum", minimum))
    except Exception as exc:
        findings.append(issue("finite", "bad-decimal", exc))
    try:
        rebuilt = build_artifact()
    except Exception as exc:
        findings.append(issue("recompute", "failed", exc))
    else:
        if artifact != rebuilt:
            findings.append(issue("recompute", "mismatch", "stored artifact differs"))
    return findings


def required_string_issues(path: Path, required: tuple[str, ...], section: str) -> list[GateIssue]:
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
            "validated Jensen-window PF compound order-four condensation gate: "
            f"8 rows, {len(findings)} issues, 2 exact identities, 1 exact sign "
            "equivalence, 317 positive lambda=-100 prefix margins, 1 strict "
            "lower-order countermodel, 1 forbidden promotion, 2 open handoffs"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
