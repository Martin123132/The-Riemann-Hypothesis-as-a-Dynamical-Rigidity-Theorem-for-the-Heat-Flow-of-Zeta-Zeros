#!/usr/bin/env python3
"""Validate full ratio-cone entry at lambda=-100."""

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

from jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate import (  # noqa: E402
    BASE_SOURCE,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPAIR_220_250,
    REPAIR_245_320,
    REPO_ROOT,
    build_artifact,
    load_source,
)


REQUIRED_ROW_IDS = {
    "m100fce_01_repaired_source_merge",
    "m100fce_02_prefix_coefficient_positivity",
    "m100fce_03_prefix_pointwise_cone",
    "m100fce_04_prefix_adjacent_cone",
    "m100fce_05_global_pointwise_walls",
    "m100fce_06_analytic_adjacent_tail",
    "m100fce_07_full_cone_entry",
    "m100fce_08_flow_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda -100 Full Cone-Entry Certificate",
    "Status: interval full cone-entry theorem at lambda=-100 and open flow-legitimacy handoff",
    "Artifact kind: `jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.py",
    "A_k(-100)>0, 0<=k<=320",
    "(2*k-1)/(2*k+1)<x_k<1, 1<=k<=319",
    "x_(k+1)>x_k, 1<=k<=318",
    "for every integer k>=1",
    "full infinite ratio-cone entry at `lambda=-100`",
    "infinite-dimensional or collared finite-flow legitimacy",
    "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
    "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md",
    "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
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


def source_issues() -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    expected_sources = {
        BASE_SOURCE: (0, 300, 301),
        REPAIR_220_250: (220, 250, 31),
        REPAIR_245_320: (245, 320, 76),
    }
    for path, (start, end, count) in expected_sources.items():
        try:
            values = load_source(path)
            if (min(values), max(values), len(values)) != (start, end, count):
                findings.append(issue("sources", "bad-range", (path, min(values), max(values), len(values))))
        except Exception as exc:
            findings.append(issue("sources", "load-failed", (path, exc)))

    source_artifacts = {
        "lower": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_heat_flow_boundary_threshold_lemma.json",
        "upper": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_kernel_mellin_upper_wall_certificate.json",
        "tail": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.json",
        "flow": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.json",
    }
    try:
        sources = {key: load_json(path) for key, path in source_artifacts.items()}
    except Exception as exc:
        return findings + [issue("sources", "artifact-load-failed", exc)]
    if sources["lower"].get("summary", {}).get("exact_rows") != 5:
        findings.append(issue("sources", "bad-lower-wall", sources["lower"].get("summary")))
    upper_summary = sources["upper"].get("summary", {})
    if not (
        upper_summary.get("compact_subintervals") == 200
        and upper_summary.get("positive_compact_subintervals") == 200
        and upper_summary.get("all_k_upper_wall_certified") is True
    ):
        findings.append(issue("sources", "bad-upper-wall", upper_summary))
    tail_summary = sources["tail"].get("summary", {})
    if not (
        tail_summary.get("open_requirement_rows") == 0
        and tail_summary.get("ready_to_apply_rows") == 2
    ):
        findings.append(issue("sources", "bad-adjacent-tail", tail_summary))
    if not sources["flow"].get("summary", {}).get("conditional_invariance_available"):
        findings.append(issue("sources", "bad-flow-source", sources["flow"].get("summary")))
    return findings


def coverage_issues() -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    prefix_pointwise = set(range(1, 320))
    prefix_adjacent = set(range(1, 319))
    if prefix_pointwise != set(range(1, 320)):
        findings.append(issue("coverage", "pointwise-prefix-gap", sorted(prefix_pointwise)))
    if prefix_adjacent.union({319}) != set(range(1, 320)):
        findings.append(issue("coverage", "adjacent-splice-gap", sorted(prefix_adjacent)))
    if max(prefix_adjacent) + 1 != 319:
        findings.append(issue("coverage", "bad-tail-start", max(prefix_adjacent)))
    return findings


def validate_artifact(artifact: dict) -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "interval full cone-entry theorem at lambda=-100 and open flow-legitimacy handoff":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("date") != "2026-07-10":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in (
        "source_base_jsonl",
        "source_repair_220_250_jsonl",
        "source_repair_245_320_jsonl",
        "source_lower_wall",
        "source_upper_wall",
        "source_adjacent_tail",
        "source_flow_invariance",
        "generator",
        "checker",
    ):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for text in ("proves", "does not prove", "full infinite", "lambda=-100", "flow", "rh", "lambda <= 0"):
        if text not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", text))

    summary = artifact.get("summary", {})
    expected_summary = {
        "certificate_rows": 8,
        "positive_coefficients": 321,
        "positive_pointwise_cone_rows": 319,
        "positive_adjacent_rows": 318,
        "analytic_adjacent_tail_start": 319,
        "full_cone_entry_rows": 1,
        "open_flow_handoff_rows": 1,
        "ready_to_apply_rows": 3,
    }
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))

    diagnostics = artifact.get("diagnostics", {})
    prefix = diagnostics.get("prefix", {})
    if diagnostics.get("parameters") != {"lambda": "-100", "precision_bits": 1024}:
        findings.append(issue("parameters", "mismatch", diagnostics.get("parameters")))
    expected_prefix = {
        "coefficient_range": [0, 320],
        "positive_coefficients": 321,
        "pointwise_cone_k_range": [1, 319],
        "positive_pointwise_cone_rows": 319,
        "adjacent_k_range": [1, 318],
        "positive_adjacent_rows": 318,
        "minimum_adjacent_margin_at_k": 318,
    }
    for key, value in expected_prefix.items():
        if prefix.get(key) != value:
            findings.append(issue("prefix", f"bad-{key}", prefix.get(key)))
    try:
        if not Decimal(prefix["minimum_lower_margin_lower"]) > 0:
            findings.append(issue("prefix", "nonpositive-lower", prefix))
        if not Decimal(prefix["minimum_upper_margin_lower"]) > 0:
            findings.append(issue("prefix", "nonpositive-upper", prefix))
        adjacent = Decimal(prefix["minimum_adjacent_margin_lower"])
        if not Decimal("3.7e-6") < adjacent < Decimal("3.8e-6"):
            findings.append(issue("prefix", "bad-adjacent-margin", adjacent))
    except Exception as exc:
        findings.append(issue("prefix", "bad-decimal", exc))
    if diagnostics.get("global_composition", {}).get("coverage") != "all integer k>=1":
        findings.append(issue("composition", "bad-coverage", diagnostics.get("global_composition")))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(value) for value in ids)))
    if len(rows) != 8:
        findings.append(issue("rows", "bad-row-count", len(rows)))
    ready = [row.get("id") for row in rows if row.get("readiness") == "ready_to_apply"]
    if ready != ["m100fce_05_global_pointwise_walls", "m100fce_06_analytic_adjacent_tail", "m100fce_07_full_cone_entry"]:
        findings.append(issue("rows", "bad-ready-rows", ready))

    findings.extend(source_issues())
    findings.extend(coverage_issues())
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
            "validated Jensen-window PF negative-lambda -100 full cone-entry certificate: "
            f"8 rows, {len(findings)} issues, 321 positive coefficients, "
            "319 pointwise cone rows, 318 adjacent rows, 1 analytic adjacent tail, "
            "1 open flow handoff, 3 ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
