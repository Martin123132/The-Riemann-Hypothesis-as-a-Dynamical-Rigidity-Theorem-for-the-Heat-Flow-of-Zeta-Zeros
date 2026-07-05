#!/usr/bin/env python3
"""Validate the non-rigorous Edrei power-Hankel midpoint warning scout.

This manifest intentionally checks a warning artifact, not a certificate.  It
records that the original loose-center midpoint scout at `m = 20, s = 1` had
mixed signs across the lambda grid.  A later tighter Arb enclosure pass
certifies the determinant positive, so this file is a guardrail against
trusting midpoint-only evidence.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path


@dataclass(frozen=True)
class FrontierSpec:
    lam: str
    safe_lam: str
    m20_midpoint: str


SPECS: tuple[FrontierSpec, ...] = (
    FrontierSpec("0", "0", "negative"),
    FrontierSpec("1e-6", "1em6", "negative"),
    FrontierSpec("1e-4", "1em4", "negative"),
    FrontierSpec("1e-2", "1em2", "positive"),
    FrontierSpec("1e-1", "1em1", "negative"),
)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_spec(results_dir: Path, spec: FrontierSpec) -> None:
    stem = f"edrei_power_hankel_midpoint_lam{spec.safe_lam}_m19_m20_s1_s1_dps2000"
    summary_path = results_dir / f"{stem}_summary.json"
    rows_path = results_dir / f"{stem}.jsonl"
    summary = load_json(summary_path)
    checks = {
        "kind": "edrei_power_hankel_midpoint_scout_summary",
        "rigorous": False,
        "lam": spec.lam,
        "m_min": 19,
        "m_max": 20,
        "shift_min": 1,
        "shift_max": 1,
        "needed_max_n": 41,
        "dps": 2000,
        "rows": 2,
    }
    for key, expected in checks.items():
        actual = summary.get(key)
        if actual != expected:
            raise AssertionError(f"{summary_path.name}: {key} {actual!r} != {expected!r}")
    if not rows_path.exists():
        raise FileNotFoundError(f"missing row log {rows_path.name}")

    seen: dict[int, dict] = {}
    with rows_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("kind") != "edrei_power_hankel_midpoint_scout":
                raise AssertionError(f"{rows_path.name}: bad row kind")
            if row.get("lam") != spec.lam:
                raise AssertionError(f"{rows_path.name}: lambda mismatch")
            if bool(row.get("rigorous")):
                raise AssertionError(f"{rows_path.name}: midpoint row marked rigorous")
            if int(row.get("shift")) != 1:
                raise AssertionError(f"{rows_path.name}: unexpected shift")
            seen[int(row["m"])] = row
    if sorted(seen) != [19, 20]:
        raise AssertionError(f"{rows_path.name}: expected m=19,20 rows, saw {sorted(seen)}")
    if seen[19].get("midpoint_classification") != "positive":
        raise AssertionError(f"{rows_path.name}: m=19 midpoint not positive")
    if seen[20].get("midpoint_classification") != spec.m20_midpoint:
        raise AssertionError(
            f"{rows_path.name}: m=20 midpoint {seen[20].get('midpoint_classification')!r} "
            f"!= {spec.m20_midpoint!r}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="work/rh_compute/results", type=Path)
    args = parser.parse_args()

    for spec in SPECS:
        validate_spec(args.results_dir, spec)
        print(
            "OK Edrei power-Hankel midpoint frontier: "
            f"lambda={spec.lam}, m=20 midpoint={spec.m20_midpoint}"
        )
    print(f"validated {len(SPECS)} non-rigorous Edrei midpoint frontier scouts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
