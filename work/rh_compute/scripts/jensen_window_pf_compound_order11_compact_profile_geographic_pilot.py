#!/usr/bin/env python3
"""Test the tuned compact H0-H23 profile across representative anchors."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order10_compact_h2_h24_unit_cache as h_source  # noqa: E402
import jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache as source  # noqa: E402
import jensen_window_pf_compound_order11_compact_high_precision_anchor_reach_pilot as reach  # noqa: E402
import jensen_window_pf_compound_order11_compact_profile_ladder_pilot as ladder  # noqa: E402
import jensen_window_pf_compound_order11_compact_stencil_precision_diagnostic as precision_diagnostic  # noqa: E402
from jensen_window_pf_compound_order11_shifted_taylor_model_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    PRECISION_BITS,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_profile_geographic_pilot.json"
)
PROFILE_NAME = "order11_ladder_p896_n96"
PROFILE = ladder.PROFILES[PROFILE_NAME]
TARGETS = tuple(
    Fraction(value)
    for value in (6000, 10000, 15004, 20004, 25004, 30004, 35004, 37800)
)
GENERATED_TARGETS = tuple(target for target in TARGETS if target != Fraction(20004))
TEST_OFFSETS = (Fraction(-2), Fraction(0), Fraction(2))
H_COLLAR = 12


def profile_task(task: tuple[int, Fraction]) -> dict:
    index, target = task
    source.profiles.PROFILE_SPECS[PROFILE_NAME] = PROFILE
    row = source.exact_task((index, target, PROFILE_NAME))
    return {**row, "kind": "order11_compact_profile_geographic_h0_h23_jet"}


def tasks() -> list[tuple[int, Fraction]]:
    return list(enumerate(GENERATED_TARGETS))


def central_ladder_record() -> dict:
    artifact = json.loads(ladder.DEFAULT_OUT.read_text(encoding="utf-8"))
    matches = [
        evaluation["source_row"]
        for evaluation in artifact.get("evaluations", [])
        if evaluation.get("profile") == PROFILE_NAME
    ]
    if len(matches) != 1:
        raise RuntimeError("profile ladder has no unique p896/n96 central row")
    precision_diagnostic.validate_point_record(matches[0], Fraction(20004))
    return matches[0]


def load_h_rows() -> list[dict]:
    selected = {
        Fraction(value)
        for target in TARGETS
        for value in range(int(target) - H_COLLAR, int(target) + H_COLLAR)
    }
    derivative_keys = {str(order) for order in range(2, 25)}
    rows: dict[Fraction, dict] = {}
    with h_source.DEFAULT_CACHE.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            left = Fraction(record["target_t_left"])
            if left not in selected:
                continue
            right = Fraction(record["target_t_right"])
            derivatives = record.get("h_derivatives", {})
            if (
                record.get("kind") != "order10_compact_h2_h24_unit_tile"
                or record.get("contract_id") != h_source.ROW_CONTRACT
                or record.get("passed") is not True
                or right != left + 1
                or set(derivatives) != derivative_keys
            ):
                raise RuntimeError(f"invalid compact H row at t={left}")
            rows[left] = {
                "target_t_left": left,
                "target_t_right": right,
                "H": {
                    order: precision_diagnostic.compact.interval_from_text(
                        derivatives[str(order)]
                    )
                    for order in range(2, 25)
                },
            }
    if set(rows) != selected:
        raise RuntimeError("compact H source misses a geographic-pilot collar")
    return [rows[target] for target in sorted(rows)]


def evaluate_record(record: dict, h_rows: list[dict]) -> dict:
    target = Fraction(record["target_t"])
    pair = precision_diagnostic.point_pair(record)
    expansion_rows = [
        reach.expansion_row(
            target + offset,
            pair,
            h_rows,
            exact_anchor=target,
        )
        for offset in TEST_OFFSETS
    ]
    completed = [
        block
        for row in expansion_rows
        for block in row["blocks"]
        if "scaled_curvature_upper" in block
    ]
    largest = (
        max(completed, key=lambda block: float(block["scaled_curvature_upper"]))
        if completed
        else None
    )
    return {
        "target": str(target),
        "source_row": record,
        "expansion_rows": expansion_rows,
        "largest_scaled_curvature_upper": (
            largest["scaled_curvature_upper"] if largest else None
        ),
        "passed": all(row["passed"] for row in expansion_rows),
    }


def build_artifact(*, workers: int = 4) -> dict:
    with ProcessPoolExecutor(max_workers=workers) as pool:
        generated = list(pool.map(profile_task, tasks()))
    if any(record.get("passed") is not True for record in generated):
        raise RuntimeError("a geographic-pilot exact source row failed")
    records = generated + [central_ladder_record()]
    records.sort(key=lambda record: Fraction(record["target_t"]))
    if [Fraction(record["target_t"]) for record in records] != list(TARGETS):
        raise RuntimeError("geographic exact source targets changed")
    flint.ctx.prec = PRECISION_BITS
    h_rows = load_h_rows()
    evaluations = [evaluate_record(record, h_rows) for record in records]
    passing = [row["target"] for row in evaluations if row["passed"]]
    return {
        "kind": "jensen_window_pf_compound_order11_compact_profile_geographic_pilot",
        "date": "2026-07-18",
        "status": "rigorous representative-anchor compact profile pilot",
        "proof_boundary": (
            "This covers six quarter blocks near each listed anchor and does not "
            "certify the compact intervals between representative anchors."
        ),
        "parameters": {
            "targets": [str(target) for target in TARGETS],
            "test_offsets": [str(offset) for offset in TEST_OFFSETS],
            "profile_name": PROFILE_NAME,
            "profile": PROFILE,
            "curvature_cap": CURVATURE_CONSTANT,
            "model_working_precision_bits": PRECISION_BITS,
        },
        "evaluations": evaluations,
        "summary": {
            "anchors": len(evaluations),
            "passing_anchors": passing,
            "all_anchors_passed": len(passing) == len(evaluations),
            "scaled_curvature_uppers": {
                row["target"]: row["largest_scaled_curvature_upper"]
                for row in evaluations
            },
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order11_compact_profile_geographic_pilot.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order11_compact_profile_geographic_pilot.py"
        ),
    }


def main() -> int:
    artifact = build_artifact()
    DEFAULT_OUT.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUT.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    summary = artifact["summary"]
    print(
        "compact geographic profile pilot: "
        f"passing={len(summary['passing_anchors'])}/{summary['anchors']}, "
        f"uppers={summary['scaled_curvature_uppers']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
