#!/usr/bin/env python3
"""Validate the order-moment transport theorem-fit gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from pathlib import Path

from jensen_window_pf_order_moment_transport_fit_gate import (
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    SOURCE_URL,
    build_artifact,
)


REQUIRED_IDS = {
    "omtf_01_primary_scope",
    "omtf_02_exact_reparametrization",
    "omtf_03_cm_prerequisite",
    "omtf_04_origin_obstruction",
    "omtf_05_orientation_mismatch",
    "omtf_06_forbidden_promotion",
    "omtf_07_live_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Order-Moment Transport Fit Gate",
    "Status: primary-source theorem-fit rejection gate",
    "https://arxiv.org/abs/2606.31647v2",
    "Gamma-normalized average of a completely monotone",
    "A_t^(1)=sqrt(pi)/(10*400^t*Gamma(t+1/2))",
    "d_u log(phi_1)(0)=5+8*pi/(2*pi-3)-4*pi",
    "not completely monotone",
    "ordinary positive Hankel moment",
    "signed reciprocal-Gamma transport",
    "rejects one theorem application",
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


def artifact_issues(artifact: dict) -> list[GateIssue]:
    findings: list[GateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_order_moment_transport_fit_gate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    expected_status = (
        "primary-source theorem-fit audit with exact first-summand "
        "complete-monotonicity obstruction"
    )
    if artifact.get("status") != expected_status:
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("primary_source", {}).get("url") != SOURCE_URL:
        findings.append(issue("source", "bad-url", artifact.get("primary_source")))
    ids = {row.get("id") for row in artifact.get("rows", [])}
    if ids != REQUIRED_IDS:
        findings.append(issue("rows", "bad-ids", sorted(str(value) for value in ids)))
    expected_summary = {
        "fit_rows": 7,
        "exact_identity_rows": 1,
        "interval_counterchecks": 1,
        "structural_mismatches": 1,
        "forbidden_promotions": 1,
        "open_handoffs": 1,
        "ready_to_apply_rows": 1,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(issue("summary", "mismatch", artifact.get("summary")))
    try:
        lower = Decimal(artifact["exact"]["origin_log_derivative_lower"])
        if not Decimal("0.08") < lower < Decimal("0.09"):
            findings.append(issue("exact", "bad-origin-lower", lower))
    except Exception as exc:
        findings.append(issue("exact", "bad-origin-decimal", exc))
    try:
        rebuilt = build_artifact()
    except Exception as exc:
        findings.append(issue("recompute", "failed", exc))
    else:
        if artifact != rebuilt:
            findings.append(issue("recompute", "mismatch", "stored artifact differs"))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for marker in (
        "rejects direct use",
        "does not reject",
        "signed total positivity",
        "does not",
        "order-four entry",
        "rh",
        "lambda<=0",
    ):
        if marker not in boundary:
            findings.append(issue("artifact", "weak-boundary", marker))
    for ref in artifact.get("sources", []):
        if isinstance(ref, str) and ref.startswith("https://"):
            continue
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(issue("artifact", "missing-source", ref))
    for key in ("generator", "checker"):
        ref = artifact.get(key)
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(issue("artifact", "missing-ref", (key, ref)))
    return findings


def note_issues(path: Path) -> list[GateIssue]:
    if not path.exists():
        return [issue("note", "missing", path)]
    text = path.read_text(encoding="utf-8")
    return [
        issue("note", "missing-string", value)
        for value in REQUIRED_NOTE_STRINGS
        if value not in text
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    findings: list[GateIssue] = []
    try:
        artifact = load_json(args.artifact)
    except Exception as exc:
        artifact = {}
        findings.append(issue("artifact", "load-failed", exc))
    if artifact:
        findings.extend(artifact_issues(artifact))
    findings.extend(note_issues(args.note))
    ok = not findings
    if args.json:
        print(
            json.dumps(
                {"ok": ok, "issues": [asdict(item) for item in findings]},
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF order-moment transport fit gate: "
            f"7 rows, {len(findings)} issues, 1 exact reparametrization, "
            "1 positive origin derivative obstruction, 1 orientation mismatch, "
            "1 forbidden promotion, 1 open signed-transport handoff"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
