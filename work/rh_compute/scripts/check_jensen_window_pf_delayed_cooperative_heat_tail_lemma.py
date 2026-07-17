#!/usr/bin/env python3
"""Validate the shifted cooperative heat-flow descent lemma."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS = REPO_ROOT / "work/rh_compute/results"
HEAT_SOURCE = RESULTS / "jensen_window_pf_all_order_endpoint_heat_reduction.json"
ORDER10_SOURCE = RESULTS / "jensen_window_pf_compound_order10_m100_tail_curvature_reduction.json"
DEFAULT_ARTIFACT = RESULTS / "jensen_window_pf_delayed_cooperative_heat_tail_lemma.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_delayed_cooperative_heat_tail_lemma.md"

FIXED_TAIL = "for every fixed m exists N_m such that for every n>=N_m and -100<=lambda<=0, Q_(m,n)(lambda)>0"
COOPERATIVE_FLOW = "Q_(m,n)'=a_(m,n)*Q_(m,n+1)+b_(m,n)*Q_(m,n), a_(m,n)=c_(m,n)*Q_(m-1,n)/Q_(m-1,n+1)>0, b_(m,n)=c_(m,n)/(c_(m,n)-4)*(log Q_(m-1,n+1))'"
VARIATION = "Q_(m,n)(lambda)=E_(m,n)(lambda)*(Q_(m,n)(-100)+integral_(-100)^lambda E_(m,n)(s)^(-1)*a_(m,n)(s)*Q_(m,n+1)(s)ds), E_(m,n)>0"
SHIFTED_LEMMA = "[Q_(m-1,n)(lambda)>0 for every n>=n0 on -100<=lambda<=0, the fixed-order m eventual tail holds, and Q_(m,n)(-100)>0 for every n>=n0] => [Q_(m,n)(lambda)>0 for every n>=n0 and -100<=lambda<=0]"
ORDER10_HANDOFF = "[Q_(10,n)(-100)>0 for every integer n>=4] implies Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"
ORDER11_HANDOFF = "[Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0 and Q_(11,n)(-100)>0 for every n>=4] implies Q_(11,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"


@dataclass(frozen=True)
class Finding:
    section: str
    issue: str
    detail: str


def finding(section: str, issue: str, detail: object) -> Finding:
    return Finding(section, issue, str(detail))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def validate(artifact_path: Path, note_path: Path) -> list[Finding]:
    findings: list[Finding] = []
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    heat = json.loads(HEAT_SOURCE.read_text(encoding="utf-8"))
    order10 = json.loads(ORDER10_SOURCE.read_text(encoding="utf-8"))

    expected_sources = (HEAT_SOURCE, ORDER10_SOURCE)
    stored_sources = {row.get("path"): row.get("sha256") for row in artifact.get("sources", [])}
    for path in expected_sources:
        if stored_sources.get(rel(path)) != sha256(path):
            findings.append(finding("sources", "hash-mismatch", rel(path)))
    if set(stored_sources) != {rel(path) for path in expected_sources}:
        findings.append(finding("sources", "path-set-mismatch", sorted(stored_sources)))

    exact_heat = heat.get("exact", {})
    expected_heat = {
        "fixed_order_tail": FIXED_TAIL,
        "cooperative_flow": COOPERATIVE_FLOW,
        "variation_of_constants": VARIATION,
        "known_base": "Q_(m,n)(lambda)>0 for 1<=m<=9, every n>=0, and -100<=lambda<=0",
    }
    for key, expected in expected_heat.items():
        if exact_heat.get(key) != expected:
            findings.append(finding("heat-source", f"bad-{key}", exact_heat.get(key)))
    heat_summary = heat.get("summary", {})
    for key, expected in {
        "all_fixed_order_tail_theorems": 1,
        "all_order_cooperative_recursions": 1,
        "completed_base_order": 9,
    }.items():
        if heat_summary.get(key) != expected:
            findings.append(finding("heat-summary", f"bad-{key}", heat_summary.get(key)))
    if order10.get("exact", {}).get("conditional_heat_handoff") != ORDER10_HANDOFF:
        findings.append(finding("order10-source", "handoff-changed", order10.get("exact", {})))

    expected_exact = {
        "fixed_order_tail": FIXED_TAIL,
        "cooperative_flow": COOPERATIVE_FLOW,
        "variation_of_constants": VARIATION,
        "finite_descent": "N=max(n0,N_m); prove n=N-1 down to n0",
        "shifted_single_layer_implication": SHIFTED_LEMMA,
        "order10_delayed_handoff": ORDER10_HANDOFF,
        "order11_delayed_handoff": ORDER11_HANDOFF,
        "unused_prefix": "the proof uses no Q_(m,n) with n<n0",
    }
    if artifact.get("exact") != expected_exact:
        findings.append(finding("artifact", "exact-contract-mismatch", artifact.get("exact")))
    if artifact.get("summary") != {
        "rows": 7,
        "ready_rows": 7,
        "open_rows": 0,
        "shifted_heat_lemmas": 1,
        "order10_n4_specializations": 1,
        "order11_n4_specializations": 1,
        "endpoint_premises_proved": 0,
    }:
        findings.append(finding("artifact", "summary-mismatch", artifact.get("summary")))

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != [f"dcht_0{i}_{name}" for i, name in enumerate(("eventual_tail", "cooperative_equation", "scalar_solution", "finite_descent", "shifted_lemma", "order10_specialization", "order11_specialization"), 1)]:
        findings.append(finding("rows", "id-order-mismatch", [row.get("id") for row in rows]))
    if any(row.get("readiness") != "ready_to_apply" for row in rows):
        findings.append(finding("rows", "non-ready-row", rows))

    note = note_path.read_text(encoding="utf-8")
    for marker in (
        "Choose `N=max(n0,N_m)`",
        "finite descending induction",
        "proof uses no `Q_(m,n)` with `n<n0`",
        ORDER10_HANDOFF,
        ORDER11_HANDOFF,
        "conditional on proving the delayed endpoint premise",
        "does not assert either premise",
        "not\nPF-infinity, RH, or `Lambda<=0`",
    ):
        if marker not in note:
            findings.append(finding("note", "missing-marker", marker))
    return findings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    findings = validate(args.artifact, args.note)
    if findings:
        for row in findings:
            print(f"{row.section}: {row.issue}: {row.detail}")
        print(f"delayed cooperative heat-tail lemma: {len(findings)} issues")
        return 1
    print(
        "validated delayed cooperative heat-tail lemma: 1 shifted theorem, "
        "order-ten and order-eleven specializations, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
