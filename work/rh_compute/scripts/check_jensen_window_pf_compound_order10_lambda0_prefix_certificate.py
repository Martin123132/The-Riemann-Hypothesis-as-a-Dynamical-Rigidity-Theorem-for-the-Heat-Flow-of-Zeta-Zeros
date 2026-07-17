#!/usr/bin/env python3
"""Independently validate the order-ten lambda-zero prefix certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


ORDER = 10
PROBE_M = 9
MAX_K = 21
PRECISION_DPS = 520
RESULTS = REPO_ROOT / "work/rh_compute/results"
COEFFICIENT_SOURCES = (
    RESULTS / "acb_enclosures_repro_hankel_15c_lamgrid_k0_k9.jsonl",
    RESULTS / "acb_enclosures_repro_hankel_15c_lamgrid_k10_k20.jsonl",
    RESULTS / "acb_enclosures_repro_hankel_15c_lamgrid_k21_k32.jsonl",
)
HANKEL_ROWS_SOURCE = RESULTS / "arb_hankel_enclosure_lam0_m0_m19_s0_s24_dps520.jsonl"
HANKEL_SUMMARY_SOURCE = RESULTS / "arb_hankel_enclosure_lam0_m0_m19_s0_s24_dps520_summary.json"
DEFINITION_SOURCE = RESULTS / "jensen_window_pf_all_order_endpoint_heat_reduction.json"
DEFAULT_ARTIFACT = RESULTS / "jensen_window_pf_compound_order10_lambda0_prefix_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_compound_order10_lambda0_prefix_certificate.md"


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


def relative_path(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT)).replace("\\", "/")


def decimal_equal(left: object, right: object) -> bool:
    return Decimal(str(left)) == Decimal(str(right))


def read_coefficients(findings: list[Finding]) -> dict[int, flint.arb]:
    coefficients: dict[int, flint.arb] = {}
    for path in COEFFICIENT_SOURCES:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                row = json.loads(line)
                if row.get("kind") != "acb_coefficient_enclosure":
                    continue
                if not decimal_equal(row.get("lam"), "0"):
                    continue
                k = int(row["k"])
                if k > MAX_K:
                    continue
                if k in coefficients:
                    findings.append(finding("coefficients", "duplicate", k))
                    continue
                coefficient = flint.arb(row["A_ball"])
                normalized = (
                    flint.arb(row["full_mu_ball"])
                    * math.factorial(k)
                    / math.factorial(2 * k)
                )
                if not coefficient.overlaps(normalized):
                    findings.append(finding("coefficients", "bad-normalization", k))
                if not coefficient > 0:
                    findings.append(finding("coefficients", "not-positive", k))
                coefficients[k] = coefficient
    missing = [k for k in range(MAX_K + 1) if k not in coefficients]
    if missing:
        findings.append(finding("coefficients", "missing", missing))
    return coefficients


def read_source_rows(findings: list[Finding]) -> dict[int, dict]:
    rows: dict[int, dict] = {}
    with HANKEL_ROWS_SOURCE.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            if (
                row.get("kind") == "arb_hankel_enclosure_sign_probe"
                and decimal_equal(row.get("lam"), "0")
                and int(row.get("m", -1)) == PROBE_M
                and 0 <= int(row.get("shift", -1)) <= 3
            ):
                n = int(row["shift"])
                if n in rows:
                    findings.append(finding("source-rows", "duplicate", n))
                rows[n] = row
    missing = [n for n in range(4) if n not in rows]
    if missing:
        findings.append(finding("source-rows", "missing", missing))
    return rows


def validate(artifact_path: Path, note_path: Path) -> list[Finding]:
    findings: list[Finding] = []
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    if artifact.get("kind") != "jensen_window_pf_compound_order10_lambda0_prefix_certificate":
        findings.append(finding("artifact", "bad-kind", artifact.get("kind")))

    expected_paths = (*COEFFICIENT_SOURCES, HANKEL_ROWS_SOURCE, HANKEL_SUMMARY_SOURCE, DEFINITION_SOURCE)
    actual_sources = {row.get("path"): row.get("sha256") for row in artifact.get("sources", [])}
    for path in expected_paths:
        rel = relative_path(path)
        if actual_sources.get(rel) != sha256(path):
            findings.append(finding("sources", "hash-mismatch", rel))
    if set(actual_sources) != {relative_path(path) for path in expected_paths}:
        findings.append(finding("sources", "path-set-mismatch", sorted(actual_sources)))

    summary = json.loads(HANKEL_SUMMARY_SOURCE.read_text(encoding="utf-8"))
    summary_expected = {
        "kind": "arb_hankel_enclosure_sign_probe_summary",
        "lam": "0",
        "m_min": 0,
        "m_max": 19,
        "shift_min": 0,
        "shift_max": 24,
        "dps": 520,
        "needed_max_k": 62,
        "rows": 500,
        "ok": 500,
        "failed_or_inconclusive": 0,
        "all_ok": True,
        "counts": {"positive": 500},
    }
    for key, expected in summary_expected.items():
        if summary.get(key) != expected:
            findings.append(finding("source-summary", f"bad-{key}", summary.get(key)))

    expected_definition = (
        "H_(m,n)=det[A_(n+i+j)]_(0<=i,j<m), H_(0,n)=1, "
        "epsilon_m=(-1)^binom(m,2), Q_(m,n)=epsilon_m*H_(m,n)"
    )
    definitions = json.loads(DEFINITION_SOURCE.read_text(encoding="utf-8"))
    if definitions.get("exact", {}).get("definitions") != expected_definition:
        findings.append(finding("definition", "contract-changed", "signed Hankel definition"))

    normalization = artifact.get("normalization", {})
    normalization_expected = {
        "coefficient": "A_k(lambda)=mu_(2k)(lambda)*k!/(2k)!",
        "definition": expected_definition,
        "order": 10,
        "probe_parameter": 9,
        "probe_matrix_dimension": 10,
        "epsilon_10": -1,
        "probe_sigma": -1,
        "index_identification": "probe m=9 and shift=n equals Q_(10,n)",
    }
    if normalization != normalization_expected:
        findings.append(finding("normalization", "mismatch", normalization))

    flint.ctx.dps = PRECISION_DPS
    coefficients = read_coefficients(findings)
    source_rows = read_source_rows(findings)
    artifact_rows = {int(row["n"]): row for row in artifact.get("finite", {}).get("rows", [])}
    if set(artifact_rows) != set(range(4)):
        findings.append(finding("artifact-rows", "bad-index-set", sorted(artifact_rows)))

    if len(coefficients) == MAX_K + 1 and len(source_rows) == 4:
        for n in range(4):
            raw = flint.arb_mat(
                [[coefficients[n + i + j] for j in range(ORDER)] for i in range(ORDER)]
            ).det()
            signed = -raw
            source = source_rows[n]
            source_signed = flint.arb(source["signed_det"])
            if int(source.get("sigma", 0)) != -1:
                findings.append(finding("source-row", "bad-sigma", n))
            if source.get("classification") != "positive" or source.get("ok") is not True:
                findings.append(finding("source-row", "not-certified-positive", n))
            if source.get("contains_zero") or source.get("sign_separated") is not True:
                findings.append(finding("source-row", "not-sign-separated", n))
            if not signed > 0:
                findings.append(finding("rebuild", "Q10-not-positive", n))
            if not signed.overlaps(source_signed):
                findings.append(finding("rebuild", "source-miss", n))
            stored = artifact_rows.get(n, {})
            if not flint.arb(stored.get("Q10_ball", "0")).contains(signed):
                findings.append(finding("artifact-row", "Q10-ball-does-not-contain-rebuild", n))
            if not flint.arb(stored.get("raw_H10_ball", "0")).contains(raw):
                findings.append(finding("artifact-row", "H10-ball-does-not-contain-rebuild", n))
            expected_fields = {
                "coefficient_range": [n, n + 18],
                "probe_m": 9,
                "matrix_dimension": 10,
                "probe_sigma": -1,
                "epsilon_10": -1,
                "source_Q10_ball": source["signed_det"],
                "source_overlap": True,
                "classification": "positive",
            }
            for key, expected in expected_fields.items():
                if stored.get(key) != expected:
                    findings.append(finding("artifact-row", f"bad-{key}", n))

    finite = artifact.get("finite", {})
    if finite.get("lambda") != "0" or finite.get("n_range") != [0, 3]:
        findings.append(finding("finite", "bad-domain", finite))
    if finite.get("coefficient_range") != [0, 21] or finite.get("precision_dps") != 520:
        findings.append(finding("finite", "bad-computation-contract", finite))
    if finite.get("all_Q10_positive") is not True:
        findings.append(finding("finite", "positivity-flag-false", finite))
    if finite.get("theorem") != "Q_(10,n)(0)>0 for every integer 0<=n<=3":
        findings.append(finding("finite", "bad-theorem", finite.get("theorem")))

    expected_artifact_summary = {
        "coefficient_rows": 22,
        "prefix_rows": 4,
        "positive_Q10_rows": 4,
        "source_overlaps": 4,
        "inconclusive_rows": 0,
    }
    if artifact.get("summary") != expected_artifact_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))

    note = note_path.read_text(encoding="utf-8")
    for marker in (
        "A_k(lambda)=mu_(2k)(lambda) k!/(2k)!",
        "sigma=(-1)^(9*10/2)=-1=epsilon_10",
        "Q_(10,n)(0)>0 for every integer 0<=n<=3",
        "not an all-shift",
        "not a proof of PF-infinity, RH, or `Lambda<=0`",
    ):
        if marker not in note:
            findings.append(finding("note", "missing-marker", marker))
    return findings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings = validate(args.artifact, args.note)
    if findings:
        for row in findings:
            print(f"{row.section}: {row.issue}: {row.detail}")
        print(f"order-ten lambda-zero prefix certificate: {len(findings)} issues")
        return 1
    print(
        "validated order-ten lambda-zero prefix: "
        "22 coefficients, 4 positive Q10 rows, 4 source overlaps, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
