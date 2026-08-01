#!/usr/bin/env python3
"""Check the one-contact uniform quartic length-14 interval certificate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_branch_uniform_length14_interval_certificate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_quartic_outer_branch_uniform_length14_interval_certificate.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_quartic_outer_branch_uniform_length14_interval_certificate.py"
)


def load_builder():
    spec = importlib.util.spec_from_file_location("quartic_uniform14", BUILDER)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load uniform length-14 builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, default=RESULT)
    parser.add_argument("--note", type=Path, default=NOTE)
    args = parser.parse_args()
    issues: list[str] = []
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    expected_kind = (
        "jensen_window_pf_quartic_outer_branch_"
        "uniform_length14_interval_certificate"
    )
    if payload.get("kind") != expected_kind:
        issues.append("artifact kind changed")
    if payload.get("status") != (
        "rigorous one-contact uniform length-fourteen obstruction"
    ):
        issues.append("artifact status changed")
    if payload.get("precision_bits") != 256:
        issues.append("Arb precision changed")
    if len(payload.get("rows", [])) != 9:
        issues.append("expected nine gate rows")

    builder = load_builder()
    builder.flint.ctx.prec = builder.PRECISION_BITS
    builder.lower_process_priority()
    certificate = payload.get("interval_certificate", {})
    cache = REPO_ROOT / certificate.get("event_log", "")
    if not cache.exists():
        issues.append("event log missing")
    else:
        if certificate.get("event_log_bytes") != cache.stat().st_size:
            issues.append("event-log byte count changed")
        if certificate.get("event_log_sha256") != hashlib.sha256(
            cache.read_bytes()
        ).hexdigest():
            issues.append("event-log file hash changed")

    stack = [""]
    previous = "0" * 64
    leaves = 0
    splits = 0
    max_depth = 0
    worst_upper = float("-inf")
    nonfinite_splits = 0
    with cache.open("r", encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            if not line.endswith("\n"):
                issues.append("event log has an incomplete final line")
                break
            event = json.loads(line)
            if event.get("index") != index:
                issues.append(f"event index changed at {index}")
                break
            if event.get("prev_sha256") != previous:
                issues.append(f"event predecessor changed at {index}")
                break
            body = {
                key: value
                for key, value in event.items()
                if key not in {"prev_sha256", "event_sha256"}
            }
            observed_hash = builder.event_hash(previous, body)
            if event.get("event_sha256") != observed_hash:
                issues.append(f"event digest changed at {index}")
                break
            if not stack or event.get("path") != stack.pop():
                issues.append(f"DFS coverage path changed at {index}")
                break
            path = event["path"]
            max_depth = max(max_depth, len(path))
            if event.get("action") == "split":
                dimension = event.get("dimension")
                if not isinstance(dimension, int) or not 0 <= dimension < 8:
                    issues.append(f"bad split dimension at {index}")
                    break
                stack.append(builder.child_path(path, dimension, 1))
                stack.append(builder.child_path(path, dimension, 0))
                splits += 1
                if event.get("reason") == "nonfinite_interval":
                    nonfinite_splits += 1
            elif event.get("action") == "leaf":
                assessment = builder.assess_box(builder.decode_box(path))
                if assessment is None:
                    issues.append(f"leaf {index} no longer evaluates finitely")
                    break
                if not builder.strictly_negative(assessment["mean_upper"]):
                    issues.append(f"leaf {index} is not strictly negative")
                    break
                if builder.serialize_arb(assessment["mean_upper"]) != event.get(
                    "mean_upper_ball"
                ):
                    issues.append(f"leaf {index} upper ball changed")
                    break
                serialized = builder.flint.arb(event["mean_upper_ball"])
                if not serialized.contains(assessment["mean_upper"]):
                    issues.append(
                        f"leaf {index} serialized enclosure lost its Arb ball"
                    )
                    break
                expected_endpoint = builder.serialize_arb(
                    assessment["mean_upper"].upper()
                )
                if expected_endpoint != event.get(
                    "mean_upper_upper_endpoint"
                ):
                    issues.append(
                        f"leaf {index} upper endpoint changed"
                    )
                    break
                upper_endpoint = builder.flint.arb(expected_endpoint)
                if not builder.strictly_negative(upper_endpoint):
                    issues.append(
                        f"leaf {index} serialized upper endpoint is not negative"
                    )
                    break
                if not upper_endpoint.contains(
                    assessment["mean_upper"].upper()
                ):
                    issues.append(
                        f"leaf {index} serialized upper endpoint lost its value"
                    )
                    break
                if assessment["monotone_signs"] != event.get(
                    "monotone_signs"
                ):
                    issues.append(f"leaf {index} monotonicity signs changed")
                    break
                worst_upper = max(
                    worst_upper, float(assessment["mean_upper"].upper())
                )
                leaves += 1
            else:
                issues.append(f"unknown action at {index}")
                break
            previous = observed_hash

    if stack:
        issues.append(f"interval tree has {len(stack)} uncovered boxes")
    observed = {
        "events": splits + leaves,
        "splits": splits,
        "leaves": leaves,
        "max_depth": max_depth,
        "worst_leaf_upper": worst_upper,
        "nonfinite_splits": nonfinite_splits,
    }
    for key, value in observed.items():
        if certificate.get(key) != value:
            issues.append(f"certificate statistic changed: {key}")
    if certificate.get("final_event_sha256") != previous:
        issues.append("final event-chain digest changed")

    exact = builder.exact_corner()
    if payload.get("exact") != exact:
        issues.append("exact corner data changed")
    if exact.get("delta14_sign") != -1:
        issues.append("exact corner is not negative")

    fixed, gap_1, gap_2 = builder.base_data()
    defects = {index: 1 - value for index, value in fixed.items()}
    cap_6 = gap_2**2 / (fixed[4] ** 3 * gap_1)
    derived_y6 = (
        defects[5] ** 2
        - Fraction(11, 13)
        * defects[5]
        * fixed[5] ** 2
        * defects[4]
    ) / cap_6
    if derived_y6 != builder.ROOT_Y6_UPPER:
        issues.append("root y_6 wall failed its independent exact derivation")
    boundary_d6 = (
        defects[5] ** 2 - derived_y6 * cap_6
    ) / (
        fixed[5] ** 2 * defects[4]
    )
    if 13 * boundary_d6 != 11 * defects[5]:
        issues.append("root y_6 wall does not saturate the scaled-defect step")
    if 1 - boundary_d6 != Fraction(31983, 32500):
        issues.append("root y_6 boundary contraction changed")

    # Independent raw-gap recurrence at an interior rational fixture.
    x = dict(fixed)
    d = {index: 1 - value for index, value in x.items()}
    gaps = {1: gap_1, 2: gap_2}
    fixture = {
        6: Fraction(1, 3),
        7: Fraction(2, 3),
        8: Fraction(3, 4),
        9: Fraction(4, 5),
        10: Fraction(5, 6),
        11: Fraction(6, 7),
        12: Fraction(7, 8),
        13: Fraction(8, 9),
    }
    for index in range(6, 14):
        cap = (
            gaps[index - 4] ** 2
            / (x[index - 2] ** 3 * gaps[index - 5])
        )
        gaps[index - 3] = fixture[index] * cap
        d[index] = (
            d[index - 1] ** 2 - gaps[index - 3]
        ) / (x[index - 1] ** 2 * d[index - 2])
        x[index] = 1 - d[index]
    raw_delta = (
        gaps[10] ** 2 / (x[12] ** 3 * gaps[9])
        - (d[13] ** 2 - x[13] ** 2 * d[12] * d[13])
    )
    factored = builder.evaluate_point(
        [fixture[index] for index in range(6, 14)]
    )
    if not factored.contains(
        builder.arb_exact(raw_delta)
    ):
        issues.append("factored recurrence failed the raw rational fixture")

    summary = payload.get("summary", {})
    if summary != {
        "rows": 9,
        "fixed_outer_contacts": 1,
        "tail_parameters": 8,
        "scaled_walls_used": 1,
        "certified_leaves": leaves,
        "uniform_length14_obstructions_for_fixed_contact": 1,
        "uniform_all_contact_theorems": 0,
    }:
        issues.append("summary counts changed")

    note = args.note.read_text(encoding="utf-8")
    for marker in (
        "one-contact uniform length-fourteen obstruction",
        "13*d_6>11*d_5",
        "[0,Y_6] x [0,1]^7",
        "256-bit Arb",
        "pending boxes: 0",
        "Every upper endpoint is strictly negative",
        "not uniform over the complete `(a,p,u)`",
        "`Lambda <= 0`",
        "RH",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated one-contact uniform quartic length-14 interval "
        f"certificate: 9 rows, 0 issues, {leaves} certified leaves, "
        f"{splits} dyadic splits, maximum depth {max_depth}, "
        "8 tail parameters, 1 scaled wall, "
        "1 uniform fixed-contact obstruction, 0 uniform all-contact theorems"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
