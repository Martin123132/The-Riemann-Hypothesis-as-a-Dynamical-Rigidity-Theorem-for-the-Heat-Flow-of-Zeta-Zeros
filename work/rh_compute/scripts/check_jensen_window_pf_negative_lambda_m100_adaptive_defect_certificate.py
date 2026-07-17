#!/usr/bin/env python3
"""Validate the lambda=-100 adaptive-defect theorem composition."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_IDS = {
    "m100adc_01_defect_coordinates",
    "m100adc_02_full_cone_input",
    "m100adc_03_raw_corridor_input",
    "m100adc_04_defect_cone",
    "m100adc_05_defect_monotonicity",
    "m100adc_06_scaled_defect_cone",
    "m100adc_07_scaled_defect_growth",
    "m100adc_08_parameter_specific_closure",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda -100 Adaptive-Defect Certificate",
    "Status: exact interval-analytic adaptive-defect theorem at lambda=-100.",
    "Artifact kind: `jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate`.",
    "validated Jensen-window PF negative-lambda -100 adaptive-defect certificate: 8 rows, 0 issues, 2 theorem inputs, 4 defect conclusions, 0 open requirements",
    "0<=d_k<=2/(2*k+1),",
    "d_(k+1)<=d_k,",
    "0<=s_k<=1.",
    "s_(k+1)-s_k",
    "The older targets asked for the stronger simultaneous statement",
)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_algebra() -> list[str]:
    k, raw, raw_next = sp.symbols("k raw raw_next", positive=True)
    x = ((2 * k - 1) / (2 * k + 1)) * raw
    x_next = ((2 * k + 1) / (2 * k + 3)) * raw_next
    defect = 1 - x
    defect_next = 1 - x_next
    scaled = ((2 * k + 1) / 2) * defect
    scaled_next = ((2 * k + 3) / 2) * defect_next
    lower_raw = ((2 * k - 1) * (2 * k + 3) / (2 * k + 1) ** 2) * raw
    upper_raw = (2 + (2 * k - 1) * raw) / (2 * k + 1)

    identities = {
        "scaled coordinate": scaled - ((2 * k + 1) - (2 * k - 1) * raw) / 2,
        "defect monotonicity": defect - defect_next - (x_next - x),
        "raw lower transfer": (
            x_next
            - x
            - ((2 * k + 1) / (2 * k + 3)) * (raw_next - lower_raw)
        ),
        "raw upper transfer": (
            scaled_next
            - scaled
            - ((2 * k + 1) / 2) * (upper_raw - raw_next)
        ),
    }
    return [name for name, value in identities.items() if sp.simplify(value) != 0]


def validate(artifact: dict, note_path: Path) -> list[str]:
    issues: list[str] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_m100_adaptive_defect_certificate":
        issues.append(f"bad kind: {artifact.get('kind')!r}")
    if artifact.get("status") != "exact interval-analytic adaptive-defect theorem at lambda=-100":
        issues.append(f"bad status: {artifact.get('status')!r}")
    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_IDS:
        issues.append(
            f"bad row ids: missing={sorted(REQUIRED_IDS-ids)}, extra={sorted(ids-REQUIRED_IDS)}"
        )
    expected_summary = {
        "certificate_rows": 8,
        "theorem_input_rows": 2,
        "defect_conclusion_rows": 4,
        "parameter_specific_closure_rows": 1,
        "open_requirements": 0,
        "ready_to_apply_rows": 7,
    }
    if artifact.get("summary") != expected_summary:
        issues.append("bad summary")
    for key in (
        "source_full_cone",
        "source_raw_corridor",
        "source_adaptive_target",
        "source_defect_tail_target",
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
    issues.extend(f"bad algebra: {name}" for name in validate_algebra())
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
    for marker in (
        "all-k defect cone",
        "does not prove the stronger simultaneous",
        "pf-infinity",
        "rh",
        "lambda <= 0",
    ):
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
        print(f"lambda=-100 adaptive-defect certificate issues: {len(issues)}")
        for item in issues:
            print(f"- {item}")
        return 1
    print(
        "validated Jensen-window PF negative-lambda -100 adaptive-defect certificate: "
        "8 rows, 0 issues, 2 theorem inputs, 4 defect conclusions, 0 open requirements"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
