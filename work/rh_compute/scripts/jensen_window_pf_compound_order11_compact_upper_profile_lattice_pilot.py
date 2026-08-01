#!/usr/bin/env python3
"""Tune rigorous step-four H0-H23 lattices on the upper compact region."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
import json
import os
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
import jensen_window_pf_compound_order11_compact_profile_geographic_pilot as geographic  # noqa: E402
import jensen_window_pf_compound_order11_compact_step4_geographic_pilot as step4  # noqa: E402
import jensen_window_pf_compound_order11_compact_stencil_precision_diagnostic as precision_diagnostic  # noqa: E402
from jensen_window_pf_compound_order11_shifted_taylor_model_core import (  # noqa: E402
    CURVATURE_CONSTANT,
    PRECISION_BITS,
)


DEFAULT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_upper_profile_lattice_h0_h23.jsonl"
)
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_compact_upper_profile_lattice_pilot.json"
)
CENTERS = tuple(Fraction(value) for value in (30004, 35004, 37800))
ANCHOR_SHIFTS = step4.ANCHOR_SHIFTS
TEST_OFFSETS = step4.TEST_OFFSETS
ROWS_PER_PROFILE = len(CENTERS) * len(ANCHOR_SHIFTS)
PROFILES = {
    "order11_upper_p1536_n128": {
        "precision_bits": 1536,
        "mode_bisections": 320,
        "panels": 128,
        "row_contract": "order11-compact-upper-h0-h23-p1536-b320-n128-w15-t30-v1",
    },
    "order11_upper_p1280_n120": {
        "precision_bits": 1280,
        "mode_bisections": 272,
        "panels": 120,
        "row_contract": "order11-compact-upper-h0-h23-p1280-b272-n120-w15-t30-v1",
    },
    "order11_upper_p1152_n112": {
        "precision_bits": 1152,
        "mode_bisections": 244,
        "panels": 112,
        "row_contract": "order11-compact-upper-h0-h23-p1152-b244-n112-w15-t30-v1",
    },
    "order11_upper_p1024_n112": {
        "precision_bits": 1024,
        "mode_bisections": 220,
        "panels": 112,
        "row_contract": "order11-compact-upper-h0-h23-p1024-b220-n112-w15-t30-v1",
    },
    "order11_upper_p896_n108": {
        "precision_bits": 896,
        "mode_bisections": 200,
        "panels": 108,
        "row_contract": "order11-compact-upper-h0-h23-p896-b200-n108-w15-t30-v1",
    },
    "order11_upper_p896_n104": {
        "precision_bits": 896,
        "mode_bisections": 200,
        "panels": 104,
        "row_contract": "order11-compact-upper-h0-h23-p896-b200-n104-w15-t30-v1",
    },
}


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def selected_profile_names(max_profiles: int) -> tuple[str, ...]:
    if not 1 <= max_profiles <= len(PROFILES):
        raise ValueError(f"max_profiles must be between 1 and {len(PROFILES)}")
    return tuple(PROFILES)[:max_profiles]


def all_tasks() -> list[tuple[int, str, Fraction]]:
    tasks = []
    index = 0
    for profile_name in PROFILES:
        for center in CENTERS:
            for shift in ANCHOR_SHIFTS:
                tasks.append((index, profile_name, center + shift))
                index += 1
    return tasks


def exact_task(task: tuple[int, str, Fraction]) -> dict:
    index, profile_name, target = task
    source.profiles.PROFILE_SPECS[profile_name] = PROFILES[profile_name]
    row = source.exact_task((index, target, profile_name))
    return {**row, "kind": "order11_compact_upper_profile_lattice_h0_h23_jet"}


def validate_source_row(record: dict, task: tuple[int, str, Fraction]) -> None:
    index, profile_name, target = task
    source.profiles.PROFILE_SPECS[profile_name] = PROFILES[profile_name]
    precision_diagnostic.validate_point_record(record, target)
    if (
        record.get("kind")
        != "order11_compact_upper_profile_lattice_h0_h23_jet"
        or record.get("index") != index
        or record.get("profile") != profile_name
        or record.get("contract_id") != source.row_contract(profile_name)
    ):
        raise RuntimeError(f"invalid upper compact source row {index}")


def load_source_cache(path: Path = DEFAULT_CACHE) -> list[dict]:
    if not path.exists():
        return []
    records = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise RuntimeError(
                    f"invalid upper compact JSONL row {line_number}"
                ) from exc
    expected = all_tasks()
    if len(records) > len(expected):
        raise RuntimeError("upper compact source cache has too many rows")
    for record, task in zip(records, expected):
        validate_source_row(record, task)
    return records


def build_source_cache(
    profile_names: tuple[str, ...],
    *,
    workers: int,
    generate: bool,
    path: Path = DEFAULT_CACHE,
) -> list[dict]:
    expected = all_tasks()
    stop = ROWS_PER_PROFILE * len(profile_names)
    records = load_source_cache(path)
    if len(records) >= stop:
        return records[:stop]
    if not generate:
        raise RuntimeError(
            f"upper compact source cache has {len(records)}/{stop} required rows"
        )
    remaining = expected[len(records) : stop]
    path.parent.mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(max_workers=workers) as pool, path.open(
        "a", encoding="utf-8", newline="\n"
    ) as handle:
        for record, task in zip(pool.map(exact_task, remaining), remaining):
            validate_source_row(record, task)
            handle.write(json.dumps(record, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
            records.append(record)
    return records[:stop]


def profile_evaluation(
    profile_name: str,
    records: list[dict],
    h_rows: list[dict],
) -> dict:
    profile_records = {
        Fraction(record["target_t"]): record
        for record in records
        if record["profile"] == profile_name
    }
    expected_targets = {
        center + shift for center in CENTERS for shift in ANCHOR_SHIFTS
    }
    if set(profile_records) != expected_targets:
        raise RuntimeError(f"profile {profile_name} has the wrong exact grid")
    centers = [
        step4.center_evaluation(center, profile_records, h_rows)
        for center in CENTERS
    ]
    passing_centers = [row["center"] for row in centers if row["passed"]]
    return {
        "profile_name": profile_name,
        "profile": PROFILES[profile_name],
        "centers": centers,
        "passing_centers": passing_centers,
        "all_centers_passed": len(passing_centers) == len(CENTERS),
        "scaled_curvature_uppers": {
            row["center"]: row["largest_scaled_curvature_upper"]
            for row in centers
        },
    }


def build_artifact(
    profile_names: tuple[str, ...],
    *,
    workers: int = 4,
    generate: bool = True,
) -> dict:
    if profile_names != tuple(PROFILES)[: len(profile_names)]:
        raise ValueError("evaluated profiles must be a deterministic prefix")
    records = build_source_cache(
        profile_names,
        workers=workers,
        generate=generate,
    )
    flint.ctx.prec = PRECISION_BITS
    h_rows = geographic.load_h_rows()
    evaluations = [
        profile_evaluation(profile_name, records, h_rows)
        for profile_name in profile_names
    ]
    passing_profiles = [
        row["profile_name"] for row in evaluations if row["all_centers_passed"]
    ]
    return {
        "kind": "jensen_window_pf_compound_order11_compact_upper_profile_lattice_pilot",
        "date": "2026-07-18",
        "status": "rigorous upper-compact step-four source-profile pilot",
        "proof_boundary": (
            "This covers six quarter blocks near each listed upper center for "
            "each evaluated profile; it is not contiguous compact coverage."
        ),
        "parameters": {
            "centers": [str(center) for center in CENTERS],
            "anchor_shifts": [str(shift) for shift in ANCHOR_SHIFTS],
            "test_offsets": [str(offset) for offset in TEST_OFFSETS],
            "rows_per_profile": ROWS_PER_PROFILE,
            "evaluated_profiles": list(profile_names),
            "profile_specs": {name: PROFILES[name] for name in profile_names},
            "curvature_cap": CURVATURE_CONSTANT,
            "model_working_precision_bits": PRECISION_BITS,
        },
        "evaluations": evaluations,
        "source_cache": {
            "path": relative(DEFAULT_CACHE),
            "sha256": sha256(DEFAULT_CACHE),
            "rows": len(records),
        },
        "summary": {
            "evaluated_profiles": len(evaluations),
            "exact_source_rows": len(records),
            "passing_profiles": passing_profiles,
            "first_passing_profile": passing_profiles[0] if passing_profiles else None,
            "all_profile_uppers": {
                row["profile_name"]: row["scaled_curvature_uppers"]
                for row in evaluations
            },
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order11_compact_upper_profile_lattice_pilot.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order11_compact_upper_profile_lattice_pilot.py"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-profiles", type=int, default=1)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    profile_names = selected_profile_names(args.max_profiles)
    artifact = build_artifact(
        profile_names,
        workers=max(1, args.workers),
        generate=True,
    )
    DEFAULT_OUT.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUT.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    summary = artifact["summary"]
    print(
        "upper compact profile lattice: "
        f"profiles={summary['evaluated_profiles']}, "
        f"passing={summary['passing_profiles']}, "
        f"uppers={summary['all_profile_uppers']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
