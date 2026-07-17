#!/usr/bin/env python3
"""Validate the lambda=-100 raw-corridor theorem composition."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_m100_raw_corridor_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_IDS = {
    "m100rcc_01_raw_curvature_coordinates",
    "m100rcc_02_full_cone_input",
    "m100rcc_03_scaled_curvature_input",
    "m100rcc_04_linear_to_nonlinear_barrier",
    "m100rcc_05_curvature_corridor",
    "m100rcc_06_raw_ratio_corridor",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda -100 Raw-Corridor Certificate",
    "Status: exact interval-analytic raw-corridor theorem at lambda=-100.",
    "Artifact kind: `jensen_window_pf_negative_lambda_m100_raw_corridor_certificate`.",
    "validated Jensen-window PF negative-lambda -100 raw-corridor certificate: 6 rows, 0 issues, 2 theorem inputs, 1 raw-corridor theorem, 0 open requirements",
    "B_k>=0, B_(k+1)<=B_k.",
    "B_(k+1)>=((2*k+1)/(2*k+3))*B_k.",
    "This closes the zeta-specific raw-corridor target at lambda=-100.",
)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_calculus() -> list[str]:
    k, b = sp.symbols("k b", positive=True)
    alpha = (2 * k + 1) / (2 * k + 3)
    y = (2 * k + 1) * sp.exp(-b)
    left_derivative = y / (2 + y)
    derivative_gap = sp.factor(alpha - left_derivative)
    expected = sp.factor(
        2 * (2 * k + 1) * (1 - sp.exp(-b))
        / ((2 * k + 3) * (2 + (2 * k + 1) * sp.exp(-b)))
    )
    if sp.simplify(derivative_gap - expected) != 0:
        return [f"bad calculus derivative gap: {sp.simplify(derivative_gap - expected)}"]
    return []


def validate(artifact: dict, note_path: Path) -> list[str]:
    issues: list[str] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_m100_raw_corridor_certificate":
        issues.append(f"bad kind: {artifact.get('kind')!r}")
    if artifact.get("status") != "exact interval-analytic raw-corridor theorem at lambda=-100":
        issues.append(f"bad status: {artifact.get('status')!r}")
    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_IDS:
        issues.append(f"bad row ids: missing={sorted(REQUIRED_IDS-ids)}, extra={sorted(ids-REQUIRED_IDS)}")
    expected_summary = {
        "certificate_rows": 6,
        "exact_reduction_rows": 2,
        "theorem_input_rows": 2,
        "raw_corridor_theorem_rows": 1,
        "compact_source_blocks": 16074,
        "finite_prefix_source_gaps": 318,
        "open_requirements": 0,
        "ready_to_apply_rows": 4,
    }
    if artifact.get("summary") != expected_summary:
        issues.append("bad summary")
    for key in (
        "source_scaled_curvature",
        "source_full_cone",
        "source_exact_bridge",
        "source_target",
        "generator",
        "checker",
    ):
        ref = artifact.get(key)
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(f"missing ref {key}: {ref!r}")
    try:
        recomputed = build_artifact()
    except Exception as exc:
        issues.append(f"recompute failed: {type(exc).__name__}: {exc}")
    else:
        for key in ("rows", "summary", "status", "proof_boundary"):
            if artifact.get(key) != recomputed.get(key):
                issues.append(f"recomputed {key} differs")
    issues.extend(validate_calculus())
    if not note_path.exists():
        issues.append(f"missing note: {note_path}")
    else:
        text = note_path.read_text(encoding="utf-8")
        issues.extend(
            f"note missing: {required}"
            for required in REQUIRED_NOTE_STRINGS
            if required not in text
        )
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for marker in ("proves the full raw-ratio", "does not prove pf-infinity", "rh", "lambda <= 0"):
        if marker not in boundary:
            issues.append(f"weak proof boundary: {marker}")
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = load_json(args.artifact)
    issues = validate(artifact, args.note)
    if issues:
        print(f"lambda=-100 raw-corridor certificate issues: {len(issues)}")
        for item in issues:
            print(f"- {item}")
        return 1
    print(
        "validated Jensen-window PF negative-lambda -100 raw-corridor certificate: "
        "6 rows, 0 issues, 2 theorem inputs, 1 raw-corridor theorem, 0 open requirements"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
