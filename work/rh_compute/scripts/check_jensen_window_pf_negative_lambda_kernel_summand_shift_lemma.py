#!/usr/bin/env python3
"""Validate the exact Newman-kernel summand-shift lemma."""

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

from jensen_window_pf_negative_lambda_kernel_summand_shift_lemma import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "ksl_01_kernel_summand_definition",
    "ksl_02_exact_summand_shift",
    "ksl_03_exact_moment_shift",
    "ksl_04_relative_integrand_ratio",
    "ksl_05_compact_monotonicity",
    "ksl_06_t100_k300_compact_certificate",
    "ksl_07_far_tail_partition",
    "ksl_08_fixed_t_saddle_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Kernel Summand-Shift Lemma",
    "Status: exact kernel-shift lemma with compact interval certificate",
    "Artifact kind: `jensen_window_pf_negative_lambda_kernel_summand_shift_lemma`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.py",
    "validated Jensen-window PF negative-lambda kernel summand-shift lemma: 8 rows, 0 issues, 6 exact rows, 1 compact interval row, 1 open far-tail row, 0 ready-to-apply rows",
    "phi_n(u)=n^(-1/2)*phi_1(u+a_n)",
    "rho_(n,k,T)(v)=n^(-1/2)*(1-a_n/v)^(2k)",
    "d_v log rho=2*k*a_n/(v*(v-a_n))+2*T*a_n>0",
    "a_21>3/2",
    "2.122e-29",
    "dominant n=1 moment",
    "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
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


def validate_ref(ref: object) -> list[LemmaIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def validate_artifact(artifact: dict) -> list[LemmaIssue]:
    findings: list[LemmaIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_kernel_summand_shift_lemma":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "exact kernel-shift lemma with compact interval certificate":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("date") != "2026-07-10":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in (
        "source_kernel_log_concavity",
        "source_t1156_counterexample",
        "source_k300_audit",
        "source_cone_entry_target",
        "generator",
        "checker",
    ):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for text in ("exact all-n", "does not bound", "dominant n=1", "rh", "lambda <= 0"):
        if text not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", text))

    params = artifact.get("parameters", {})
    if params != {
        "T": 100,
        "K": 300,
        "V": "1.5",
        "compact_n_range": [2, 20],
        "precision_bits": 512,
    }:
        findings.append(issue("parameters", "unexpected", params))

    summary = artifact.get("summary", {})
    expected = {
        "lemma_rows": 8,
        "available_exact_rows": 6,
        "interval_validated_rows": 1,
        "open_requirement_rows": 1,
        "compact_summand_count": 19,
        "first_far_only_summand": 21,
        "ready_to_apply_rows": 0,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))
    try:
        compact_upper = Decimal(summary["compact_ratio_sum_upper"])
        if not (Decimal("2.121e-29") < compact_upper < Decimal("2.122e-29")):
            findings.append(issue("summary", "bad-compact-upper", compact_upper))
    except Exception as exc:
        findings.append(issue("summary", "bad-decimal", exc))

    diagnostics = artifact.get("diagnostics", {})
    if diagnostics.get("compact_ratio_below_cap") is not True:
        findings.append(issue("diagnostics", "cap-not-certified", diagnostics.get("compact_ratio_below_cap")))
    if diagnostics.get("a_21_above_V") is not True:
        findings.append(issue("diagnostics", "a21-not-above-split", diagnostics.get("a_21_above_V")))
    if len(diagnostics.get("compact_ratio_rows", [])) != 19:
        findings.append(issue("diagnostics", "bad-ratio-row-count", len(diagnostics.get("compact_ratio_rows", []))))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(item) for item in ids)))
    if len(rows) != 8:
        findings.append(issue("rows", "bad-row-count", len(rows)))
    for row in rows:
        if not isinstance(row, dict) or not row.get("claim") or not row.get("proof_boundary"):
            findings.append(issue("rows", "incomplete-row", row))
    try:
        recomputed = build_artifact()
        for key in ("parameters", "diagnostics", "summary", "rows"):
            if artifact.get(key) != recomputed.get(key):
                findings.append(issue("recompute", f"mismatch-{key}", "stored differs"))
    except Exception as exc:
        findings.append(issue("recompute", "exception", exc))
    return findings


def validate_note(path: Path) -> list[LemmaIssue]:
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
    findings: list[LemmaIssue] = []
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
            "validated Jensen-window PF negative-lambda kernel summand-shift lemma: "
            f"8 rows, {len(findings)} issues, 6 exact rows, 1 compact interval row, "
            "1 open far-tail row, 0 ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
