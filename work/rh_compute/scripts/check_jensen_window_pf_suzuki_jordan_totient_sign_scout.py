#!/usr/bin/env python3
"""Validate the finite Suzuki Jordan-totient sign scout."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_jordan_totient_sign_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_suzuki_jordan_totient_sign_scout.md"
)
GENERATOR = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_suzuki_jordan_totient_sign_scout.py"
)

REQUIRED_NOTE_TEXT = (
    "finite double-precision reconnaissance report",
    "c_omega(n)",
    "h_omega^<1>(x)",
    "All 6,000 sampled values are positive",
    "g_(1/2)^<1>(exp(-4))",
    "g_(1/2)^<1>(1/2)",
    "genuinely signed arithmetic convolution",
    "Finite positivity cannot be promoted to eventual sign",
    "no uniformity as `omega->0`",
    "not a proof",
)


def load_generator():
    spec = importlib.util.spec_from_file_location("suzuki_sign_scout", GENERATOR)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load scout generator")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def validate(artifact: Path, note: Path) -> list[str]:
    issues: list[str] = []
    if not artifact.exists():
        return [f"missing artifact: {artifact}"]
    if not note.exists():
        return [f"missing note: {note}"]

    payload = json.loads(artifact.read_text(encoding="utf-8"))
    if (
        payload.get("kind")
        != "jensen_window_pf_suzuki_jordan_totient_sign_scout"
    ):
        issues.append("unexpected artifact kind")

    grid_spec = payload.get("grid", {})
    expected_grid = {
        "linear_start": 1.001,
        "linear_stop": 20.0,
        "linear_count": 600,
        "geometric_start": 20.01,
        "geometric_stop": 5000.0,
        "geometric_count": 600,
        "unique_count": 1200,
    }
    for key, expected in expected_grid.items():
        if grid_spec.get(key) != expected:
            issues.append(
                f"grid mismatch for {key}: "
                f"{grid_spec.get(key)!r} != {expected!r}"
            )

    summaries = payload.get("summaries", [])
    if len(summaries) != 5:
        issues.append(f"expected 5 omega summaries, found {len(summaries)}")
    labels = [summary.get("omega") for summary in summaries]
    if labels != ["1/2", "1/4", "1/8", "1/16", "1/32"]:
        issues.append(f"unexpected omega labels: {labels}")

    generator = load_generator()
    grid = generator.x_grid()
    for omega, label, stored in zip(
        generator.OMEGAS,
        generator.OMEGA_LABELS,
        summaries,
        strict=True,
    ):
        coefficients = generator.jordan_coefficients(omega, 5000)
        if not np.all(coefficients[1:] > 0.0):
            issues.append(f"nonpositive coefficient at omega={label}")

        probe_x = (1.001, 10.0, 100.0, 1000.0, 5000.0)
        probe_values = {
            x: generator.scaled_h(omega, x, coefficients)
            for x in probe_x
        }
        if not all(math.isfinite(value) and value > 0.0 for value in probe_values.values()):
            issues.append(f"nonpositive or nonfinite anchor at omega={label}")

        stored_anchors = stored.get("anchors", {})
        for x in (10.0, 100.0, 1000.0, 5000.0):
            key = format(x, ".0f")
            if not math.isclose(
                float(stored_anchors.get(key, math.nan)),
                probe_values[x],
                rel_tol=2e-12,
                abs_tol=2e-12,
            ):
                issues.append(f"anchor mismatch omega={label}, x={key}")

        if stored.get("sample_count") != int(grid.size):
            issues.append(f"sample count mismatch at omega={label}")
        if stored.get("negative_count") != 0:
            issues.append(f"stored negative count is nonzero at omega={label}")
        if stored.get("nonfinite_count") != 0:
            issues.append(f"stored nonfinite count is nonzero at omega={label}")
        if float(stored.get("minimum_scaled_h", 0.0)) <= 0.08:
            issues.append(f"minimum sampled margin too small at omega={label}")

    guard = payload.get("signed_weight_guard", {})
    negative = float(guard.get("negative_probe_value", math.nan))
    positive = float(guard.get("positive_probe_value", math.nan))
    if not (
        math.isfinite(negative)
        and math.isfinite(positive)
        and negative < -1.0
        and positive > 1.0
        and guard.get("sign_change_certified_by_double_probes") is True
    ):
        issues.append("signed-weight guard failed")

    audit = payload.get("audit", {})
    expected_audit = {
        "omega_count": 5,
        "total_sample_count": 6000,
        "total_negative_count": 0,
        "all_coefficients_positive": True,
        "all_sampled_h_positive": True,
        "eventual_sign_proved": False,
        "rh_proved": False,
        "lambda_le_zero_proved": False,
    }
    for key, expected in expected_audit.items():
        if audit.get(key) != expected:
            issues.append(
                f"audit mismatch for {key}: {audit.get(key)!r} != {expected!r}"
            )

    text = note.read_text(encoding="utf-8")
    lower = text.lower()
    for marker in REQUIRED_NOTE_TEXT:
        if marker.lower() not in lower:
            issues.append(f"note missing marker: {marker}")
    for forbidden in (
        "therefore eventual positivity holds",
        "this finite grid proves",
        "we have proved rh",
        "this proves rh",
        "proves lambda<=0",
        "rigorous rounding certificate",
    ):
        if forbidden in lower:
            issues.append(f"forbidden promotion language: {forbidden}")
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues = validate(args.artifact, args.note)
    for issue in issues:
        print(f"SUZUKI-JORDAN-SIGN issue: {issue}")
    print(
        "validated Suzuki Jordan-totient sign scout: "
        f"6000 samples, 5 omega values, 0 negative rows, "
        f"{len(issues)} issues"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
