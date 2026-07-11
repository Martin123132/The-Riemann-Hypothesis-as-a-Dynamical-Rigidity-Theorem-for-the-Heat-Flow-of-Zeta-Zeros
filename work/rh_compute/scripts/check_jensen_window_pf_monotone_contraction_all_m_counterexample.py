#!/usr/bin/env python3
"""Validate the monotone-contraction all-m column counterexample."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import jensen_window_pf_monotone_contraction_all_m_counterexample as counterexample


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_all_m_counterexample.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_monotone_contraction_all_m_counterexample.md"

EXPECTED_VALUE = "-2611049410212561054670871/163840000000000000000"
EXPECTED_RESULT = (
    "validated Jensen-window PF monotone-contraction all-m counterexample: "
    "degree 7, m=11, exact full-cone witness, 6 lower walls, negative normalized value, 0 issues"
)


def load_payload(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_payload(payload: dict, expected_payload: dict) -> list[str]:
    issues: list[str] = []
    if payload != expected_payload:
        issues.append("stored payload differs from exact regenerated all-m counterexample")
        return issues

    if payload.get("degree") != 7 or payload.get("minor_size") != 11:
        issues.append("counterexample row is not degree 7 m=11")
    witness = payload.get("witness", {})
    if witness.get("x_contractions") != ["19/20", "19/20", "1", "1", "1", "1"]:
        issues.append("witness contractions drifted")
    if witness.get("x_from_s") != witness.get("x_contractions"):
        issues.append("s-parameter witness does not map to x witness")
    if witness.get("q_value") != EXPECTED_VALUE:
        issues.append("negative exact value drifted")
    if witness.get("q_value_is_negative") is not True:
        issues.append("witness is not marked negative")
    checks = witness.get("monotone_checks", {})
    for key in ("all_in_unit_interval", "weakly_increasing", "strictly_positive_first", "upper_bound_ok"):
        if checks.get(key) is not True:
            issues.append(f"monotone check {key} failed")
    full_cone = witness.get("full_cone_checks", {})
    if full_cone.get("shift") != 0:
        issues.append("full-cone witness shift drifted")
    if full_cone.get("lower_walls") != ["1/3", "3/5", "5/7", "7/9", "9/11", "11/13"]:
        issues.append("pointwise lower walls drifted")
    if full_cone.get("minimum_lower_wall_margin") != "2/13":
        issues.append("minimum lower-wall margin drifted")
    for key in ("all_lower_walls_ok", "infinite_extension_obeys_lower_wall"):
        if full_cone.get(key) is not True:
            issues.append(f"full-cone check {key} failed")
    summary = payload.get("summary", {})
    if summary.get("ready_to_apply_rows") != 0 or summary.get("target_closing") is not False:
        issues.append("proof-boundary readiness flags drifted")
    return issues


def validate_note(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    required = [
        EXPECTED_RESULT,
        "Status: exact countermodel gate. This is not evidence against RH,",
        "the full static ratio cone does not imply all-`m` column recurrence",
        "x = (19/20, 19/20, 1, 1, 1, 1)",
        "lower_walls = (1/3, 3/5, 5/7, 7/9, 9/11, 11/13)",
        "tail = x_k=1 for every k>=7",
        f"Q_11 = {EXPECTED_VALUE}",
        "This value is strictly negative.",
        "outputs/signed_hankel_jensen_dependency_graph.md",
    ]
    return [f"missing note text: {needle}" for needle in required if needle not in text]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-path", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = load_payload(args.json_path)
    expected_payload = counterexample.build_payload()
    issues = validate_payload(payload, expected_payload)
    issues.extend(validate_note(args.note))
    for issue in issues:
        print(f"ISSUE {issue}")
    if not issues:
        print(EXPECTED_RESULT)
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
