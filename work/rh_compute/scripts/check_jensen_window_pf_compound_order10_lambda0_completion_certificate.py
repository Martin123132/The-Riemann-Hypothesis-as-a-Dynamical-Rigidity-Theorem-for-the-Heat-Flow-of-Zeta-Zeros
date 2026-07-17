#!/usr/bin/env python3
"""Validate the order-ten lambda-zero completion certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS = REPO_ROOT / "work/rh_compute/results"
DELAYED_ENTRY_SOURCE = RESULTS / "jensen_window_pf_compound_order10_m100_delayed_entry_certificate.json"
HEAT_LEMMA_SOURCE = RESULTS / "jensen_window_pf_delayed_cooperative_heat_tail_lemma.json"
PREFIX_SOURCE = RESULTS / "jensen_window_pf_compound_order10_lambda0_prefix_certificate.json"
ALL_ORDER_HEAT_SOURCE = RESULTS / "jensen_window_pf_all_order_endpoint_heat_reduction.json"
TRANSFER_SOURCE = RESULTS / "jensen_window_pf_order4_noncontiguous_total_positivity_transfer.json"
DEFAULT_ARTIFACT = RESULTS / "jensen_window_pf_compound_order10_lambda0_completion_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order10_lambda0_completion_certificate.md"

DELAYED_ENDPOINT = "Q_(10,n)(-100)>0 for every integer n>=4"
HEAT_RAY = "Q_(10,n)(lambda)>0 for every n>=4 and -100<=lambda<=0"
PREFIX = "Q_(10,n)(0)>0 for every integer 0<=n<=3"
ALL_SHIFT_Q10 = "Q_(10,n)(0)>0 for every integer n>=0"
LOWER_BASE = "Q_(m,n)(lambda)>0 for 1<=m<=9, every n>=0, and -100<=lambda<=0"
CONTIGUOUS_THROUGH_10 = "Q_(m,n)(0)>0 for every integer 1<=m<=10 and n>=0"
TRANSFER_THEOREM = "[epsilon_k H_(k,s)>0 for 1<=k<=m, all s] => [epsilon_k R_(k,n)(j_1,...,j_k)>0 for 1<=k<=m]"
ARBITRARY_THROUGH_10 = "epsilon_k*R_(k,n)(j_1,...,j_k)(0)>0 for 1<=k<=10, n>=0, and 0<=j_1<...<j_k"


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
    delayed = json.loads(DELAYED_ENTRY_SOURCE.read_text(encoding="utf-8"))
    lemma = json.loads(HEAT_LEMMA_SOURCE.read_text(encoding="utf-8"))
    prefix = json.loads(PREFIX_SOURCE.read_text(encoding="utf-8"))
    heat = json.loads(ALL_ORDER_HEAT_SOURCE.read_text(encoding="utf-8"))
    transfer = json.loads(TRANSFER_SOURCE.read_text(encoding="utf-8"))

    source_paths = (DELAYED_ENTRY_SOURCE, HEAT_LEMMA_SOURCE, PREFIX_SOURCE, ALL_ORDER_HEAT_SOURCE, TRANSFER_SOURCE)
    stored = {row.get("path"): row.get("sha256") for row in artifact.get("sources", [])}
    for path in source_paths:
        if stored.get(rel(path)) != sha256(path):
            findings.append(finding("sources", "hash-mismatch", rel(path)))
    if set(stored) != {rel(path) for path in source_paths}:
        findings.append(finding("sources", "path-set-mismatch", sorted(stored)))

    contracts = (
        ("delayed", delayed.get("exact", {}).get("delayed_entry"), DELAYED_ENDPOINT),
        ("lemma", lemma.get("exact", {}).get("order10_delayed_handoff"), f"[{DELAYED_ENDPOINT}] implies {HEAT_RAY}"),
        ("prefix", prefix.get("finite", {}).get("theorem"), PREFIX),
        ("lower-base", heat.get("exact", {}).get("known_base"), LOWER_BASE),
    )
    for name, actual, expected in contracts:
        if actual != expected:
            findings.append(finding("source-contract", name, actual))
    if delayed.get("summary", {}).get("delayed_m100_entry_theorems") != 1:
        findings.append(finding("source-summary", "delayed-not-closed", delayed.get("summary")))
    if lemma.get("summary", {}).get("order10_n4_specializations") != 1:
        findings.append(finding("source-summary", "lemma-not-closed", lemma.get("summary")))
    if prefix.get("summary", {}).get("positive_Q10_rows") != 4:
        findings.append(finding("source-summary", "prefix-not-closed", prefix.get("summary")))
    if heat.get("summary", {}).get("completed_base_order") != 9:
        findings.append(finding("source-summary", "lower-base-not-closed", heat.get("summary")))

    transfer_rows = {row.get("id"): row for row in transfer.get("rows", [])}
    if transfer_rows.get("o4ntp_09_fixed_order_transfer", {}).get("formula") != TRANSFER_THEOREM:
        findings.append(finding("transfer", "formula-changed", transfer_rows.get("o4ntp_09_fixed_order_transfer")))
    if transfer.get("summary", {}).get("fixed_order_transfer_theorems") != 1:
        findings.append(finding("transfer", "not-closed", transfer.get("summary")))

    expected_exact = {
        "delayed_endpoint": DELAYED_ENDPOINT,
        "delayed_heat_ray": HEAT_RAY,
        "lambda0_prefix": PREFIX,
        "all_shift_order10_lambda0": ALL_SHIFT_Q10,
        "lower_base": LOWER_BASE,
        "contiguous_through_order10_lambda0": CONTIGUOUS_THROUGH_10,
        "arbitrary_column_definition": "R_(k,n)(j_1,...,j_k)(0)=det[A_(n+i+j_l)(0)]_(0<=i<k,1<=l<=k)",
        "fixed_order_transfer": TRANSFER_THEOREM,
        "arbitrary_columns_through_order10_lambda0": ARBITRARY_THROUGH_10,
    }
    if artifact.get("exact") != expected_exact:
        findings.append(finding("artifact", "exact-mismatch", artifact.get("exact")))
    expected_summary = {
        "rows": 6,
        "ready_rows": 6,
        "open_rows": 0,
        "all_shift_order10_lambda0_theorems": 1,
        "contiguous_through_order10_lambda0_theorems": 1,
        "arbitrary_column_through_order10_lambda0_theorems": 1,
        "orders_above_10": 0,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("artifact", "summary-mismatch", artifact.get("summary")))
    rows = artifact.get("rows", [])
    expected_ids = [
        "co10l0c_01_delayed_endpoint",
        "co10l0c_02_heat_ray",
        "co10l0c_03_lambda0_prefix",
        "co10l0c_04_all_shift_order10",
        "co10l0c_05_contiguous_through10",
        "co10l0c_06_fixed_order_transfer",
    ]
    if [row.get("id") for row in rows] != expected_ids:
        findings.append(finding("rows", "id-order-mismatch", [row.get("id") for row in rows]))
    if any(row.get("readiness") != "ready_to_apply" for row in rows):
        findings.append(finding("rows", "non-ready-row", rows))

    note = note_path.read_text(encoding="utf-8")
    for marker in (
        DELAYED_ENDPOINT,
        HEAT_RAY,
        ALL_SHIFT_Q10,
        CONTIGUOUS_THROUGH_10,
        ARBITRARY_THROUGH_10,
        "nothing here proves order eleven or any higher order",
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
        print(f"order-ten lambda-zero completion: {len(findings)} issues")
        return 1
    print("validated order-ten lambda-zero completion: all shifts, arbitrary columns through order ten, 0 issues")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
