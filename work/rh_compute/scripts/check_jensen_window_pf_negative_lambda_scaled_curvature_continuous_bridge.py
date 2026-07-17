#!/usr/bin/env python3
"""Validate the lambda=-100 scaled-curvature continuous bridge."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    certify_scaled_curvature_mode_block,
)
from jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    finite_prefix_certificate,
    high_eighth_envelope,
    ray_certificate,
)


REQUIRED_ROW_IDS = {
    "nlscb_01_normalized_curvature",
    "nlscb_02_tent_derivative_bridge",
    "nlscb_03_buffered_sufficient_condition",
    "nlscb_04_low_eighth_envelope",
    "nlscb_05_high_eighth_envelope",
    "nlscb_06_compact_buffer_certificate",
    "nlscb_07_ray_buffer_certificate",
    "nlscb_08_full_kernel_transfer",
    "nlscb_09_finite_prefix_certificate",
    "nlscb_10_full_scaled_curvature_theorem",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Scaled-Curvature Continuous Bridge",
    "Status: interval and analytic all-k scaled-curvature theorem at lambda=-100.",
    "This is not a proof of PF-infinity",
    "Artifact kind: `jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge`.",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.py",
    "16074 compact blocks",
    "318 prefix gaps",
    "1 all-k scaled-curvature theorem",
    "0 open requirements",
    "Q(r)-2*abs(H'''(r))>64/(r-3)^5",
    "u^2*t*((2*t+1)*V'''/V''^3-2/V'') >= 1/5",
    "C_(k+1)>=C_k for every integer k>=1.",
    "the two-sided raw-ratio corridor at lambda=-100",
)


@dataclass(frozen=True)
class Issue:
    section: str
    name: str
    detail: str


def issue(section: str, name: str, detail: object) -> Issue:
    return Issue(section=section, name=name, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[Issue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[Issue]:
    issues: list[Issue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "interval and analytic all-k scaled-curvature theorem at lambda=-100":
        issues.append(issue("artifact", "bad-status", artifact.get("status")))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for marker in (
        "proves c_(k+1)>=c_k",
        "lambda=-100",
        "does not prove pf-infinity",
        "all-order jensen bridge",
        "rh",
        "lambda <= 0",
    ):
        if marker not in boundary:
            issues.append(issue("artifact", "weak-proof-boundary", marker))
    for key in (
        "source_paired_compact",
        "source_paired_ray",
        "source_leading",
        "source_dominance",
        "source_full_cone",
        "source_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("artifact", artifact.get(key)))
    return issues


def validate_rows(artifact: dict) -> list[Issue]:
    rows = artifact.get("rows")
    if not isinstance(rows, list):
        return [issue("rows", "bad-rows", type(rows))]
    issues: list[Issue] = []
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue("rows", "missing-row", missing))
    if len(rows) != 10:
        issues.append(issue("rows", "bad-row-count", len(rows)))
    allowed_readiness = {
        "available_exact",
        "interval_validated",
        "ready_to_apply",
    }
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing>"))
        for key in ("role", "claim", "formula", "readiness", "proof_boundary"):
            if not row.get(key):
                issues.append(issue(row_id, "missing-field", key))
        if row.get("readiness") not in allowed_readiness:
            issues.append(issue(row_id, "bad-readiness", row.get("readiness")))
    final = rows_by_id.get("nlscb_10_full_scaled_curvature_theorem", {})
    if final.get("readiness") != "ready_to_apply":
        issues.append(issue("nlscb_10", "final-row-not-ready", final.get("readiness")))
    return issues


def validate_summary(artifact: dict) -> list[Issue]:
    expected = {
        "certificate_rows": 10,
        "compact_blocks": 16074,
        "positive_compact_buffer_blocks": 16074,
        "positive_compact_transfer_blocks": 16074,
        "high_envelope_intervals": 3000,
        "positive_prefix_gaps": 318,
        "analytic_ray_rows": 1,
        "full_scaled_curvature_theorem_rows": 1,
        "ready_to_apply_rows": 2,
        "open_requirements": 0,
    }
    summary = artifact.get("summary", {})
    issues: list[Issue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    return issues


def validate_recomputed(artifact: dict) -> list[Issue]:
    issues: list[Issue] = []
    diagnostics = artifact.get("diagnostics", {})
    try:
        high = high_eighth_envelope()
        ray = ray_certificate()
        prefix = finite_prefix_certificate()
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    if diagnostics.get("high_eighth_envelope") != high:
        issues.append(issue("recompute", "bad-high-envelope", "recorded result differs"))
    if diagnostics.get("ray_certificate") != ray:
        issues.append(issue("recompute", "bad-ray", "recorded result differs"))
    if diagnostics.get("finite_prefix") != prefix:
        issues.append(issue("recompute", "bad-prefix", "recorded result differs"))

    compact = diagnostics.get("compact_certificate", {})
    selected = compact.get("selected_blocks", [])
    if len(selected) < 5:
        issues.append(issue("recompute", "too-few-selected-blocks", len(selected)))
        return issues
    for row in selected:
        try:
            result = certify_scaled_curvature_mode_block(
                Fraction(row["left"]),
                Fraction(row["right"]),
                int(row["panels"]),
                window_y=int(row["window_y"]),
                eighth_envelope_bound=Fraction(row["eighth_envelope_bound"]),
            )
        except Exception as exc:
            issues.append(
                issue(
                    "recompute",
                    "selected-block-failed",
                    f"{row.get('left')}..{row.get('right')}: {type(exc).__name__}: {exc}",
                )
            )
            continue
        if not result.get("passed"):
            issues.append(
                issue(
                    "recompute",
                    "selected-block-not-positive",
                    f"{row.get('left')}..{row.get('right')}",
                )
            )
    return issues


def validate_note(path: Path) -> list[Issue]:
    if not path.exists():
        return [issue("note", "missing-note", path)]
    text = path.read_text(encoding="utf-8")
    return [
        issue("note", "missing-text", required)
        for required in REQUIRED_NOTE_STRINGS
        if required not in text
    ]


def validate_artifact(artifact: dict, note_path: Path) -> list[Issue]:
    issues: list[Issue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_rows(artifact))
    issues.extend(validate_summary(artifact))
    issues.extend(validate_recomputed(artifact))
    issues.extend(validate_note(note_path))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = load_json(args.artifact)
    issues = validate_artifact(artifact, args.note)
    if issues:
        print(f"scaled-curvature continuous bridge issues: {len(issues)}")
        for item in issues:
            print(f"- [{item.section}] {item.name}: {item.detail}")
        return 1
    summary = artifact["summary"]
    print(
        "validated Jensen-window PF negative-lambda scaled-curvature continuous bridge: "
        f"{summary['certificate_rows']} rows, 0 issues, {summary['compact_blocks']} compact blocks, "
        f"{summary['positive_prefix_gaps']} prefix gaps, 1 analytic ray, "
        "1 all-k scaled-curvature theorem, 0 open requirements"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
