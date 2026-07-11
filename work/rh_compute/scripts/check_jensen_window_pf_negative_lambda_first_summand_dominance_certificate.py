#!/usr/bin/env python3
"""Validate the all-k first-summand dominance certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_first_summand_dominance_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "fsdc_01_exact_original_ratio",
    "fsdc_02_ratio_monotonicity",
    "fsdc_03_zero_endpoint_tail_sum",
    "fsdc_04_first_integrand_strict_concavity",
    "fsdc_05_adaptive_saddle_bracket",
    "fsdc_06_low_region_probability_bound",
    "fsdc_07_high_region_kernel_tail_bound",
    "fsdc_08_full_moment_dominance",
    "fsdc_09_adjacent_wall_stability",
    "fsdc_10_dominant_wall_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda First-Summand Dominance Certificate",
    "Status: all-k analytic first-summand dominance certificate",
    "Artifact kind: `jensen_window_pf_negative_lambda_first_summand_dominance_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_dominance_certificate.py",
    "validated Jensen-window PF negative-lambda first-summand dominance certificate: 10 rows, 0 issues, 4 exact rows, 5 interval rows, 15 positive analytic gates, 1 open dominant-wall row, 0 ready-to-apply rows",
    "d_q log(r_n)=-(n^2-1)*(1+6/((2*n^2*q-3)*(2*q-3)))<0",
    "a(k)=log(k)/8",
    "S_k''(u)=-2*k/u^2-200-16*q-96*q/(2*q-3)^2<0",
    "1/2+1/(L-alpha)-2/(L+beta)>0",
    "(6-(3*pi/2)*sqrt(k))/k<0",
    "0<=delta_k<=2/k^6",
    "|L_k-L_k^(1)|<=16/(k-1)^6",
    "L_k^(1)>=1/(4*k^2)",
    "outputs/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
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


def validate_ref(ref: object) -> list[CertificateIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def validate_artifact(artifact: dict) -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_first_summand_dominance_certificate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "all-k analytic first-summand dominance certificate":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("date") != "2026-07-10":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in (
        "source_shift_lemma",
        "source_k300_audit",
        "source_cone_entry_target",
        "source_raw_corridor_target",
        "generator",
        "checker",
    ):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for text in ("all-k", "does not prove", "dominant n=1", "rh", "lambda <= 0"):
        if text not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", text))

    summary = artifact.get("summary", {})
    expected = {
        "certificate_rows": 10,
        "available_exact_rows": 4,
        "interval_validated_rows": 5,
        "open_requirement_rows": 1,
        "positive_analytic_gates": 15,
        "full_tail_power": 6,
        "full_tail_constant": 2,
        "tail_start_k": 300,
        "wall_transfer_start_k": 301,
        "ready_to_apply_rows": 0,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))

    diagnostics = artifact.get("diagnostics", {})
    params = diagnostics.get("parameters", {})
    if params.get("T") != 100 or params.get("K") != 300:
        findings.append(issue("parameters", "bad-T-or-K", params))
    gates = diagnostics.get("positive_gates", {})
    if len(gates) != 15 or not all(gates.values()) or diagnostics.get("all_positive_gates") is not True:
        findings.append(issue("diagnostics", "bad-positive-gates", gates))
    propagation = diagnostics.get("half_line_propagation", {})
    if len(propagation) != 9:
        findings.append(issue("diagnostics", "bad-propagation-row-count", len(propagation)))
    for required in ("Q'(k)", "strengthen for k>=300", "(6-(3*pi/2)*sqrt(k))/k<0"):
        if required not in " ".join(str(value) for value in propagation.values()):
            findings.append(issue("diagnostics", "missing-propagation-formula", required))
    try:
        eps_upper = Decimal(diagnostics["epsilon_zero"]["upper"])
        if not (Decimal("0.00217") < eps_upper < Decimal("0.0022")):
            findings.append(issue("diagnostics", "bad-epsilon-zero", eps_upper))
        if not Decimal(diagnostics["high_region"]["endpoint_log_margin_upper"]) < Decimal(-100):
            findings.append(issue("diagnostics", "weak-high-endpoint", diagnostics["high_region"]))
        if not Decimal(diagnostics["low_region"]["endpoint_log_margin_upper"]) < Decimal(-3):
            findings.append(issue("diagnostics", "weak-low-endpoint", diagnostics["low_region"]))
        if not Decimal(diagnostics["saddle_geometry"]["G_a_minus_100_lower"]) > Decimal(300):
            findings.append(issue("diagnostics", "weak-Ga-margin", diagnostics["saddle_geometry"]))
        if not Decimal(diagnostics["saddle_geometry"]["G_c_lower_ball"].split("+/-", 1)[0].strip("[ ")) > Decimal(200):
            findings.append(issue("diagnostics", "weak-Gc-margin", diagnostics["saddle_geometry"]))
    except Exception as exc:
        findings.append(issue("diagnostics", "bad-decimal", exc))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(item) for item in ids)))
    if len(rows) != 10:
        findings.append(issue("rows", "bad-row-count", len(rows)))
    open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
    if len(open_rows) != 1 or open_rows[0].get("id") != "fsdc_10_dominant_wall_handoff":
        findings.append(issue("rows", "bad-open-row", open_rows))
    for row in rows:
        if not isinstance(row, dict) or not row.get("claim") or not row.get("proof_boundary"):
            findings.append(issue("rows", "incomplete-row", row))
    try:
        recomputed = build_artifact()
        for key in ("diagnostics", "summary", "rows"):
            if artifact.get(key) != recomputed.get(key):
                findings.append(issue("recompute", f"mismatch-{key}", "stored differs"))
    except Exception as exc:
        findings.append(issue("recompute", "exception", exc))
    return findings


def validate_note(path: Path) -> list[CertificateIssue]:
    if not path.exists():
        return [issue("note", "missing", path)]
    text = path.read_text(encoding="utf-8")
    return [issue("note", "missing-text", value) for value in REQUIRED_NOTE_STRINGS if value not in text]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings: list[CertificateIssue] = []
    try:
        artifact = load_json(args.artifact)
    except Exception as exc:
        artifact = {}
        findings.append(issue("artifact", "load-failed", exc))
    if artifact:
        findings.extend(validate_artifact(artifact))
    findings.extend(validate_note(args.note))
    ok = not findings
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in findings]}, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF negative-lambda first-summand dominance certificate: "
            f"10 rows, {len(findings)} issues, 4 exact rows, 5 interval rows, "
            "15 positive analytic gates, 1 open dominant-wall row, 0 ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
