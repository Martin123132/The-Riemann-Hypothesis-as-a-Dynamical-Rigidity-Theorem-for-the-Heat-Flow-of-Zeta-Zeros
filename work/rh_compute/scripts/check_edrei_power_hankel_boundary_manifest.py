#!/usr/bin/env python3
"""Validate Edrei power-Hankel frontier blocker and repair artifacts.

This is a history-and-provenance guard.  It records the former lambda=1e-6
frontier cells that stayed Arb-width inconclusive with the old dps=2000
coefficient input, and it also validates the tighter dps=2400/tol=1e-140
repair rows that now promote the completed 2m+s<=51 staircase.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path


@dataclass(frozen=True)
class BoundarySpec:
    lam: str
    safe_lam: str
    m: int
    shift: int
    dps: int
    radius_hint: str


SPECS: tuple[BoundarySpec, ...] = (
    BoundarySpec("1e-6", "1em6", 24, 3, 2000, "3.46e-2728"),
    BoundarySpec("1e-6", "1em6", 25, 1, 2000, "6.56e-2732"),
)


@dataclass(frozen=True)
class RepairSpec:
    lam: str
    safe_lam: str
    m: int
    shift_min: int
    shift_max: int
    rows: int
    dps: int


REPAIR_SPECS: tuple[RepairSpec, ...] = (
    RepairSpec("1e-6", "1em6", 24, 2, 3, 2, 2400),
    RepairSpec("1e-6", "1em6", 25, 1, 1, 1, 2400),
)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def validate_summary(path: Path, spec: BoundarySpec) -> None:
    summary = load_json(path)
    checks = {
        "kind": "arb_edrei_power_hankel_probe_summary",
        "rigorous": True,
        "lam": spec.lam,
        "m_min": spec.m,
        "m_max": spec.m,
        "shift_min": spec.shift,
        "shift_max": spec.shift,
        "needed_max_n": 2 * spec.m + spec.shift,
        "dps": spec.dps,
        "rows": 1,
        "ok": 0,
        "failed_or_inconclusive": 1,
        "all_ok": False,
    }
    for key, expected in checks.items():
        actual = summary.get(key)
        if actual != expected:
            raise AssertionError(f"{path.name}: {key} {actual!r} != {expected!r}")
    if summary.get("counts") != {"inconclusive_contains_zero": 1}:
        raise AssertionError(f"{path.name}: unexpected counts {summary.get('counts')!r}")
    if "not PF-infinity or RH" not in str(summary.get("proof_boundary", "")):
        raise AssertionError(f"{path.name}: missing proof-boundary warning")


def validate_rows(path: Path, spec: BoundarySpec) -> None:
    rows = load_rows(path)
    if len(rows) != 1:
        raise AssertionError(f"{path.name}: expected 1 row, saw {len(rows)}")
    row = rows[0]
    checks = {
        "kind": "arb_edrei_power_hankel_probe",
        "lam": spec.lam,
        "m": spec.m,
        "shift": spec.shift,
        "classification": "inconclusive_contains_zero",
        "contains_zero": True,
        "ok": False,
    }
    for key, expected in checks.items():
        actual = row.get(key)
        if actual != expected:
            raise AssertionError(f"{path.name}: {key} {actual!r} != {expected!r}")
    det = str(row.get("det", ""))
    if "[+/-" not in det or spec.radius_hint not in det:
        raise AssertionError(f"{path.name}: unexpected determinant enclosure {det!r}")


def validate_spec(results_dir: Path, spec: BoundarySpec) -> None:
    stem = (
        f"arb_edrei_power_hankel_lam{spec.safe_lam}_"
        f"m{spec.m}_m{spec.m}_s{spec.shift}_s{spec.shift}_"
        f"dps{spec.dps}_frontier_tol1e-120"
    )
    power_dir = (
        results_dir / "arb_edrei" / "power_hankel" / f"lam{spec.safe_lam}"
    )
    summary_path = power_dir / f"{stem}_summary.json"
    rows_path = power_dir / f"{stem}.jsonl"
    validate_summary(summary_path, spec)
    validate_rows(rows_path, spec)


def validate_repair_sources(results_dir: Path) -> None:
    coeff_path = (
        results_dir
        / "acb_enclosures_edrei_boundary_lam1em6_k0_k51_dps220_tol1e-140_summary.json"
    )
    coeff_summary = load_json(coeff_path)
    coeff_checks = {
        "kind": "acb_coefficient_enclosure_summary",
        "lambdas": ["0.000001"],
        "k_min": 0,
        "k_max": 51,
        "dps": 220,
        "abs_tol": "1e-140",
        "rows": 52,
    }
    for key, expected in coeff_checks.items():
        actual = coeff_summary.get(key)
        if actual != expected:
            raise AssertionError(f"{coeff_path.name}: {key} {actual!r} != {expected!r}")

    log_path = (
        results_dir
        / "arb_edrei"
        / "log_sign"
        / "arb_edrei_log_sign_lam1em6_n1_n51_dps2400_boundary_tol1e-140_summary.json"
    )
    log_summary = load_json(log_path)
    log_checks = {
        "kind": "arb_edrei_log_sign_probe_summary",
        "lam": "1e-6",
        "max_n": 51,
        "dps": 2400,
        "rows": 51,
        "failed_or_inconclusive": 0,
        "all_ok": True,
    }
    for key, expected in log_checks.items():
        actual = log_summary.get(key)
        if actual != expected:
            raise AssertionError(f"{log_path.name}: {key} {actual!r} != {expected!r}")


def validate_repair_spec(results_dir: Path, spec: RepairSpec) -> None:
    stem = (
        f"arb_edrei_power_hankel_lam{spec.safe_lam}_"
        f"m{spec.m}_m{spec.m}_s{spec.shift_min}_s{spec.shift_max}_"
        f"dps{spec.dps}_boundary_tol1e-140"
    )
    power_dir = (
        results_dir / "arb_edrei" / "power_hankel" / f"lam{spec.safe_lam}"
    )
    summary_path = power_dir / f"{stem}_summary.json"
    rows_path = power_dir / f"{stem}.jsonl"
    summary = load_json(summary_path)
    checks = {
        "kind": "arb_edrei_power_hankel_probe_summary",
        "rigorous": True,
        "lam": spec.lam,
        "m_min": spec.m,
        "m_max": spec.m,
        "shift_min": spec.shift_min,
        "shift_max": spec.shift_max,
        "needed_max_n": 51,
        "dps": spec.dps,
        "rows": spec.rows,
        "ok": spec.rows,
        "failed_or_inconclusive": 0,
        "all_ok": True,
    }
    for key, expected in checks.items():
        actual = summary.get(key)
        if actual != expected:
            raise AssertionError(f"{summary_path.name}: {key} {actual!r} != {expected!r}")
    if summary.get("counts") != {"positive": spec.rows}:
        raise AssertionError(f"{summary_path.name}: unexpected counts {summary.get('counts')!r}")

    rows = load_rows(rows_path)
    if len(rows) != spec.rows:
        raise AssertionError(f"{rows_path.name}: expected {spec.rows} rows, saw {len(rows)}")
    expected_shifts = list(range(spec.shift_min, spec.shift_max + 1))
    seen_shifts: list[int] = []
    for row in rows:
        if row.get("kind") != "arb_edrei_power_hankel_probe":
            raise AssertionError(f"{rows_path.name}: bad row kind")
        if row.get("lam") != spec.lam or row.get("m") != spec.m:
            raise AssertionError(f"{rows_path.name}: row target mismatch {row!r}")
        if row.get("classification") != "positive":
            raise AssertionError(f"{rows_path.name}: non-positive row {row!r}")
        if row.get("contains_zero") or not row.get("ok"):
            raise AssertionError(f"{rows_path.name}: non-separated row {row!r}")
        seen_shifts.append(int(row["shift"]))
    if sorted(seen_shifts) != expected_shifts:
        raise AssertionError(f"{rows_path.name}: shifts {seen_shifts!r} != {expected_shifts!r}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="work/rh_compute/results", type=Path)
    args = parser.parse_args()

    for spec in SPECS:
        validate_spec(args.results_dir, spec)
        print(
            "OK retired Edrei power-Hankel boundary blocker: "
            f"lambda={spec.lam}, m={spec.m}, shift={spec.shift}, dps={spec.dps}"
        )
    validate_repair_sources(args.results_dir)
    repaired_rows = 0
    for spec in REPAIR_SPECS:
        validate_repair_spec(args.results_dir, spec)
        repaired_rows += spec.rows
        print(
            "OK repaired Edrei power-Hankel boundary: "
            f"lambda={spec.lam}, m={spec.m}, shifts={spec.shift_min}..{spec.shift_max}"
        )
    print(
        f"validated {len(SPECS)} retired inconclusive blocker rows "
        f"and {repaired_rows} repaired positive boundary rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
