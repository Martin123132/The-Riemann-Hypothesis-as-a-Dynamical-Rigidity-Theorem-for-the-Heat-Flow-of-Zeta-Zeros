#!/usr/bin/env python3
"""Validate the first-summand cumulant bridge and finite stress scout."""

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

from jensen_window_pf_negative_lambda_first_summand_cumulant_bridge import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "fscb_01_continuous_moment_family",
    "fscb_02_cumulant_identity",
    "fscb_03_bspline_identity",
    "fscb_04_gamma_wall",
    "fscb_05_uniform_cumulant_target",
    "fscb_06_exact_rational_transfer",
    "fscb_07_high_precision_scout",
    "fscb_08_full_kernel_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda First-Summand Cumulant Bridge",
    "Status: exact cumulant reduction with analytic hypothesis discharged",
    "Artifact kind: `jensen_window_pf_negative_lambda_first_summand_cumulant_bridge`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.py",
    "validated Jensen-window PF negative-lambda first-summand cumulant bridge: 8 rows, 0 issues, 4 exact identities, 1 conditional bridge, 9 positive samples, 0 open requirements, 3 ready-to-apply rows",
    "F'''(t)=kappa_3,t(2*log(U))",
    "Delta^3 F(k-1)=integral_[0,1]^3",
    "-Delta^3 log Gamma(k-1/2)=-log(1-1/(k+1/2)^2)",
    "kappa_3,t(2*log(U))>=-37/(50*t^2), t>=318",
    "4*n^4+4108*n^3+1489493*n^2+215582064*n+9130082706",
    "Generic strict log-concavity is not enough",
    "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md",
    "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class BridgeIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> BridgeIssue:
    return BridgeIssue(section=section, issue=name, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(ref: object) -> list[BridgeIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def validate_artifact(artifact: dict) -> list[BridgeIssue]:
    findings: list[BridgeIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_first_summand_cumulant_bridge":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "exact cumulant reduction with analytic hypothesis discharged":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("date") != "2026-07-10":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in (
        "source_saddle_target",
        "source_dominance_certificate",
        "source_collar_extension",
        "source_paired_ray_certificate",
        "generator",
        "checker",
    ):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for text in ("exact", "finite", "does not prove", "cumulant", "rh", "lambda <= 0"):
        if text not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", text))

    summary = artifact.get("summary", {})
    expected = {
        "target_rows": 8,
        "exact_identity_rows": 4,
        "conditional_bridge_rows": 1,
        "sample_rows": 9,
        "above_target_sample_rows": 9,
        "scaled_profile_strictly_increasing": True,
        "minimum_sample_margin_at_t": 318,
        "open_requirement_rows": 0,
        "ready_to_apply_rows": 3,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))
    try:
        margin = Decimal(summary["minimum_sample_margin"])
        if not Decimal("0.116") < margin < Decimal("0.117"):
            findings.append(issue("summary", "bad-minimum-margin", margin))
    except Exception as exc:
        findings.append(issue("summary", "bad-margin-decimal", exc))

    diagnostics = artifact.get("diagnostics", {})
    samples = diagnostics.get("sample_rows", [])
    if len(samples) != 9:
        findings.append(issue("diagnostics", "bad-sample-count", len(samples)))
    if not all(row.get("above_target") for row in samples):
        findings.append(issue("diagnostics", "sample-below-target", samples))
    try:
        scaled = [Decimal(row["scaled_t2_kappa3"]) for row in samples]
        if not all(right > left for left, right in zip(scaled, scaled[1:])):
            findings.append(issue("diagnostics", "nonincreasing-scaled-profile", scaled))
    except Exception as exc:
        findings.append(issue("diagnostics", "bad-scaled-decimal", exc))

    transfer = diagnostics.get("exact_transfer", {})
    if transfer.get("cumulant_constant") != "37/50":
        findings.append(issue("transfer", "bad-cumulant-constant", transfer))
    if transfer.get("shifted_coefficients_ascending") != [
        9_130_082_706,
        215_582_064,
        1_489_493,
        4_108,
        4,
    ]:
        findings.append(issue("transfer", "bad-shifted-coefficients", transfer))
    try:
        if not Decimal(transfer["rational_transfer_margin_at_k319_decimal"]) > 0:
            findings.append(issue("transfer", "nonpositive-endpoint-margin", transfer))
    except Exception as exc:
        findings.append(issue("transfer", "bad-endpoint-margin", exc))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(item) for item in ids)))
    if len(rows) != 8:
        findings.append(issue("rows", "bad-row-count", len(rows)))
    if sum(row.get("role") == "open_requirement" for row in rows) != 0:
        findings.append(issue("rows", "bad-open-requirement-count", rows))

    try:
        recomputed = build_artifact()
        for key in ("diagnostics", "summary", "rows"):
            if artifact.get(key) != recomputed.get(key):
                findings.append(issue("recompute", f"mismatch-{key}", "stored differs"))
    except Exception as exc:
        findings.append(issue("recompute", "exception", exc))
    return findings


def validate_note(path: Path) -> list[BridgeIssue]:
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
    findings: list[BridgeIssue] = []
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
            "validated Jensen-window PF negative-lambda first-summand cumulant bridge: "
            f"8 rows, {len(findings)} issues, 4 exact identities, 1 conditional bridge, "
            "9 positive samples, 0 open requirements, 3 ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
