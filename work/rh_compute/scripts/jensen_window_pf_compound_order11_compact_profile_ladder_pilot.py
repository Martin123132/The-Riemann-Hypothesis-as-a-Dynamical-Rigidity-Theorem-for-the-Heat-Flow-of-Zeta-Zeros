#!/usr/bin/env python3
"""Tune exact H0-H23 source profiles for a step-four compact anchor grid."""

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

import jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache as source  # noqa: E402
import jensen_window_pf_compound_order11_compact_high_precision_anchor_reach_pilot as reach  # noqa: E402
import jensen_window_pf_compound_order11_compact_stencil_precision_diagnostic as precision_diagnostic  # noqa: E402
from jensen_window_pf_compound_order11_shifted_taylor_model_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    PRECISION_BITS,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_profile_ladder_pilot.json"
)
TARGET = Fraction(20004)
TEST_OFFSETS = (Fraction(-2), Fraction(0), Fraction(2))
PROFILES = {
    "order11_ladder_p896_n96": {
        "precision_bits": 896,
        "mode_bisections": 200,
        "panels": 96,
        "row_contract": "order11-compact-ladder-h0-h23-p896-b200-n96-w15-t30-v1",
    },
    "order11_ladder_p960_n88": {
        "precision_bits": 960,
        "mode_bisections": 200,
        "panels": 88,
        "row_contract": "order11-compact-ladder-h0-h23-p960-b200-n88-w15-t30-v1",
    },
    "order11_ladder_p1024_n80": {
        "precision_bits": 1024,
        "mode_bisections": 200,
        "panels": 80,
        "row_contract": "order11-compact-ladder-h0-h23-p1024-b200-n80-w15-t30-v1",
    },
    "order11_ladder_p1024_n88": {
        "precision_bits": 1024,
        "mode_bisections": 200,
        "panels": 88,
        "row_contract": "order11-compact-ladder-h0-h23-p1024-b200-n88-w15-t30-v1",
    },
    "order11_ladder_p1024_n96": {
        "precision_bits": 1024,
        "mode_bisections": 220,
        "panels": 96,
        "row_contract": "order11-compact-ladder-h0-h23-p1024-b220-n96-w15-t30-v1",
    },
    "order11_ladder_p1152_n104": {
        "precision_bits": 1152,
        "mode_bisections": 244,
        "panels": 104,
        "row_contract": "order11-compact-ladder-h0-h23-p1152-b244-n104-w15-t30-v1",
    },
    "order11_ladder_p1280_n112": {
        "precision_bits": 1280,
        "mode_bisections": 272,
        "panels": 112,
        "row_contract": "order11-compact-ladder-h0-h23-p1280-b272-n112-w15-t30-v1",
    },
    "order11_ladder_p1408_n120": {
        "precision_bits": 1408,
        "mode_bisections": 296,
        "panels": 120,
        "row_contract": "order11-compact-ladder-h0-h23-p1408-b296-n120-w15-t30-v1",
    },
}


def profile_task(task: tuple[int, str, dict]) -> dict:
    index, name, profile = task
    source.profiles.PROFILE_SPECS[name] = profile
    row = source.exact_task((index, TARGET, name))
    return {**row, "kind": "order11_compact_profile_ladder_h0_h23_jet"}


def tasks() -> list[tuple[int, str, dict]]:
    return [
        (index, name, profile)
        for index, (name, profile) in enumerate(PROFILES.items())
    ]


def evaluate_profile(record: dict, h_rows: list[dict]) -> dict:
    pair = precision_diagnostic.point_pair(record)
    expansion_rows = [
        reach.expansion_row(TARGET + offset, pair, h_rows)
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
        "profile": record["profile"],
        "source_row": record,
        "expansion_rows": expansion_rows,
        "largest_scaled_curvature_upper": (
            largest["scaled_curvature_upper"] if largest else None
        ),
        "passed": all(row["passed"] for row in expansion_rows),
    }


def build_artifact(*, workers: int = 4) -> dict:
    with ProcessPoolExecutor(max_workers=workers) as pool:
        records = list(pool.map(profile_task, tasks()))
    if any(record.get("passed") is not True for record in records):
        raise RuntimeError("a profile-ladder exact source row failed")
    flint.ctx.prec = PRECISION_BITS
    h_rows = reach.load_h_rows()
    evaluations = [evaluate_profile(record, h_rows) for record in records]
    passing = [row["profile"] for row in evaluations if row["passed"]]
    return {
        "kind": "jensen_window_pf_compound_order11_compact_profile_ladder_pilot",
        "date": "2026-07-18",
        "status": "rigorous local compact exact-source profile ladder",
        "proof_boundary": (
            "This tunes one anchor at t=20004 and tests offsets -2, 0, 2; "
            "it is not a full compact-interval certificate."
        ),
        "parameters": {
            "target": str(TARGET),
            "test_offsets": [str(offset) for offset in TEST_OFFSETS],
            "quarter_blocks_per_offset": 2,
            "curvature_cap": CURVATURE_CONSTANT,
            "model_working_precision_bits": PRECISION_BITS,
            "profiles": PROFILES,
        },
        "evaluations": evaluations,
        "summary": {
            "profiles": len(evaluations),
            "passing_profiles": passing,
            "first_passing_profile": passing[0] if passing else None,
            "scaled_curvature_uppers": {
                row["profile"]: row["largest_scaled_curvature_upper"]
                for row in evaluations
            },
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order11_compact_profile_ladder_pilot.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order11_compact_profile_ladder_pilot.py"
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
        "compact profile ladder: "
        f"passing={summary['passing_profiles']}, "
        f"uppers={summary['scaled_curvature_uppers']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
