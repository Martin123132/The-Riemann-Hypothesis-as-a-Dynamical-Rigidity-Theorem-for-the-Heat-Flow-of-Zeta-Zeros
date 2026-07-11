#!/usr/bin/env python3
"""Validate the degree-7 m=10 sparse subdivision certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout as scout


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_monotone_contraction_sparse_degree7_subdivision_scout.md"

EXPECTED_RESULT = (
    "validated Jensen-window PF monotone-contraction sparse degree-7 subdivision scout: "
    "3 dyadic slabs, 785400 slab Bernstein coefficients, 0 negative slab coefficients, "
    "0 zero slab coefficients, repaired m=10 obstruction, 0 issues"
)
EXPECTED_SLABS = {
    "s0 in [0,1/2]": {
        "min": "297689420329/16384",
        "min_index": [24, 0, 10, 0, 0, 0],
        "negative": 0,
        "zero": 0,
        "extra_power": 24,
    },
    "s0 in [1/2,3/4]": {
        "min": "345475759715905/268435456",
        "min_index": [24, 0, 10, 0, 0, 0],
        "negative": 0,
        "zero": 0,
        "extra_power": 48,
    },
    "s0 in [3/4,1]": {
        "min": "225767871/129536",
        "min_index": [20, 0, 10, 0, 0, 0],
        "negative": 0,
        "zero": 0,
        "extra_power": 48,
    },
}


def load_payload(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_payload(payload: dict, expected_payload: dict) -> list[str]:
    issues: list[str] = []
    if payload != expected_payload:
        issues.append("stored payload differs from exact regenerated degree-7 subdivision certificate")
        return issues

    summary = payload.get("summary", {})
    expected_summary = {
        "repaired_degree": 7,
        "repaired_minor_size": 10,
        "global_negative_bernstein_coefficients": 126,
        "global_min_coefficient": "-4928",
        "dyadic_slab_count": 3,
        "total_slab_bernstein_coefficients": 785400,
        "negative_slab_coefficients": 0,
        "zero_slab_coefficients": 0,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(f"summary {key} expected {expected!r}, got {summary.get(key)!r}")
    if summary.get("ready_to_apply_rows") != 0 or summary.get("target_closing") is not False:
        issues.append("proof-boundary readiness flags drifted")

    global_stats = payload.get("global_one_shot_stats", {})
    if global_stats.get("bernstein_min_coefficient") != "-4928":
        issues.append("global obstruction minimum drifted")
    if global_stats.get("bernstein_negative_coefficient_count") != 126:
        issues.append("global obstruction negative coefficient count drifted")
    if global_stats.get("bernstein_zero_coefficient_count") != 0:
        issues.append("global obstruction zero count drifted")

    subdivision = payload.get("subdivision", {})
    if subdivision.get("certificate_success") is not True:
        issues.append("subdivision certificate is not marked successful")
    if subdivision.get("negative_slab_count") != 0 or subdivision.get("zero_slab_count") != 0:
        issues.append("subdivision has bad slab counts")

    slabs = subdivision.get("slabs", [])
    if len(slabs) != 3:
        issues.append("expected exactly 3 slabs")
    for slab in slabs:
        label = slab.get("label")
        expected = EXPECTED_SLABS.get(label)
        if expected is None:
            issues.append(f"unexpected slab {label!r}")
            continue
        if slab.get("bernstein_coefficient_count") != 261800:
            issues.append(f"{label} has wrong coefficient count")
        if slab.get("bernstein_min_coefficient") != expected["min"]:
            issues.append(f"{label} minimum coefficient drifted")
        if slab.get("bernstein_min_index") != expected["min_index"]:
            issues.append(f"{label} minimum index drifted")
        if slab.get("bernstein_negative_coefficient_count") != expected["negative"]:
            issues.append(f"{label} negative count drifted")
        if slab.get("bernstein_zero_coefficient_count") != expected["zero"]:
            issues.append(f"{label} zero count drifted")
        if slab.get("extra_power_of_two_denominator") != expected["extra_power"]:
            issues.append(f"{label} denominator power drifted")
        if slab.get("bernstein_coefficients_strictly_positive") is not True:
            issues.append(f"{label} is not marked strictly positive")
    return issues


def validate_note(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    required = [
        EXPECTED_RESULT,
        "Status: exact sparse degree-7 subdivision certificate. This is not a",
        "s0 in [0,1/2]: count=261800",
        "s0 in [1/2,3/4]: count=261800",
        "s0 in [3/4,1]: count=261800",
        "negative=0, zero=0",
        "repairs the finite degree-7 `m=10`",
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
    expected_payload = scout.build_payload()
    issues = validate_payload(payload, expected_payload)
    issues.extend(validate_note(args.note))
    for issue in issues:
        print(f"ISSUE {issue}")
    if not issues:
        print(EXPECTED_RESULT)
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
