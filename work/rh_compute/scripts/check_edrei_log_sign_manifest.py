#!/usr/bin/env python3
"""Validate promoted finite Edrei-log sign diagnostic summaries.

This checker covers the Arb enclosure-backed necessary-condition diagnostics
for the coefficient PF route.  These rows are not PF-infinity certificates;
they only certify that a finite logarithmic sign obstruction was not found.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class LogSignSpec:
    lam: str
    safe_lam: str
    max_n: int = 32
    rows: int = 32


SPECS: tuple[LogSignSpec, ...] = (
    LogSignSpec("0", "0", max_n=64, rows=64),
    LogSignSpec("1e-6", "1em6", max_n=64, rows=64),
    LogSignSpec("1e-4", "1em4", max_n=64, rows=64),
    LogSignSpec("1e-2", "1em2", max_n=64, rows=64),
    LogSignSpec("1e-1", "1em1", max_n=64, rows=64),
)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_spec(results_dir: Path, spec: LogSignSpec) -> None:
    stem = f"arb_edrei_log_sign_lam{spec.safe_lam}_n1_n{spec.max_n}_dps340"
    summary_path = results_dir / f"{stem}_summary.json"
    rows_path = results_dir / f"{stem}.jsonl"
    summary = load_json(summary_path)

    checks = {
        "kind": "arb_edrei_log_sign_probe_summary",
        "lam": spec.lam,
        "max_n": spec.max_n,
        "rows": spec.rows,
    }
    for key, expected in checks.items():
        actual = summary.get(key)
        if actual != expected:
            raise AssertionError(f"{summary_path.name}: {key} {actual!r} != {expected!r}")
    if not bool(summary.get("rigorous")):
        raise AssertionError(f"{summary_path.name}: rigorous flag is false")
    if not bool(summary.get("all_ok")):
        raise AssertionError(f"{summary_path.name}: all_ok is false")
    if int(summary.get("failed_or_inconclusive", -1)) != 0:
        raise AssertionError(f"{summary_path.name}: failed_or_inconclusive is nonzero")
    counts = summary.get("counts", {})
    if int(counts.get("positive", 0)) != spec.rows:
        raise AssertionError(f"{summary_path.name}: positive count mismatch")

    if not rows_path.exists():
        raise FileNotFoundError(f"missing row log {rows_path}")
    seen = 0
    with rows_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            seen += 1
            n = int(row["n"])
            expected_sign = 1 if n % 2 else -1
            if row.get("kind") != "arb_edrei_log_sign_probe":
                raise AssertionError(f"{rows_path.name}: bad row kind at n={n}")
            if row.get("lam") != spec.lam:
                raise AssertionError(f"{rows_path.name}: lambda mismatch at n={n}")
            if int(row.get("expected_sign")) != expected_sign:
                raise AssertionError(f"{rows_path.name}: expected sign mismatch at n={n}")
            if row.get("classification") != "positive":
                raise AssertionError(f"{rows_path.name}: classification not positive at n={n}")
            if bool(row.get("contains_zero")):
                raise AssertionError(f"{rows_path.name}: contains_zero at n={n}")
            if not bool(row.get("ok")):
                raise AssertionError(f"{rows_path.name}: ok false at n={n}")
    if seen != spec.rows:
        raise AssertionError(f"{rows_path.name}: saw {seen} rows != {spec.rows}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="work/rh_compute/results", type=Path)
    args = parser.parse_args()

    for spec in SPECS:
        validate_spec(args.results_dir, spec)
        print(f"OK Edrei-log sign diagnostic: lambda={spec.lam}, n<={spec.max_n}")
    print(f"validated {sum(spec.rows for spec in SPECS)} finite Edrei-log sign diagnostics")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
