#!/usr/bin/env python3
"""Validate the Edrei moment-recurrence scout artifacts."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path


LAMBDAS = ["0", "1e-6", "1e-4", "1e-2", "1e-1"]


@dataclass(frozen=True)
class ScoutSpec:
    name: str
    summary_name: str
    jsonl_name: str
    order_min: int
    order_max: int
    rows: int
    positive_rows: int
    negative_rows: int
    inconclusive_rows: int
    all_positive: bool


SPECS = (
    ScoutSpec(
        name="positive Arb recurrence scout, orders 2..12",
        summary_name="edrei_moment_recurrence_lamgrid_order2_order12_arb_summary.json",
        jsonl_name="edrei_moment_recurrence_lamgrid_order2_order12_arb.jsonl",
        order_min=2,
        order_max=12,
        rows=55,
        positive_rows=55,
        negative_rows=0,
        inconclusive_rows=0,
        all_positive=True,
    ),
    ScoutSpec(
        name="frontier Arb recurrence scout, orders 2..20",
        summary_name="edrei_moment_recurrence_lamgrid_order2_order20_arb_summary.json",
        jsonl_name="edrei_moment_recurrence_lamgrid_order2_order20_arb.jsonl",
        order_min=2,
        order_max=20,
        rows=95,
        positive_rows=55,
        negative_rows=0,
        inconclusive_rows=40,
        all_positive=False,
    ),
)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_spec(results_dir: Path, spec: ScoutSpec) -> None:
    summary_path = results_dir / spec.summary_name
    row_path = results_dir / spec.jsonl_name
    summary = load_json(summary_path)
    checks = {
        "kind": "edrei_moment_recurrence_scout_summary",
        "rigorous": True,
        "lambdas": LAMBDAS,
        "order_min": spec.order_min,
        "order_max": spec.order_max,
        "rows": spec.rows,
        "positive_rows": spec.positive_rows,
        "negative_rows": spec.negative_rows,
        "inconclusive_rows": spec.inconclusive_rows,
        "all_positive": spec.all_positive,
    }
    for key, expected in checks.items():
        actual = summary.get(key)
        if actual != expected:
            raise AssertionError(f"{summary_path.name}: {key} {actual!r} != {expected!r}")

    if "not a proof of PF-infinity or RH" not in summary.get("proof_boundary", ""):
        raise AssertionError(f"{summary_path.name}: missing proof-boundary warning")

    if not row_path.exists():
        raise FileNotFoundError(f"missing row log {row_path.name}")
    row_count = 0
    positive = negative = inconclusive = 0
    with row_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            row_count += 1
            if row.get("kind") != "edrei_moment_recurrence_scout":
                raise AssertionError(f"{row_path.name}: bad row kind {row!r}")
            status = row.get("status")
            if status == "positive":
                positive += 1
                if not row.get("all_positive"):
                    raise AssertionError(f"{row_path.name}: inconsistent positive row {row!r}")
            elif status == "negative":
                negative += 1
                if row.get("all_positive"):
                    raise AssertionError(f"{row_path.name}: inconsistent negative row {row!r}")
            elif status in {"inconclusive_contains_zero", "unknown"}:
                inconclusive += 1
                if row.get("all_positive"):
                    raise AssertionError(f"{row_path.name}: inconsistent inconclusive row {row!r}")
            else:
                raise AssertionError(f"{row_path.name}: unknown row status {status!r}")

    if (row_count, positive, negative, inconclusive) != (
        spec.rows,
        spec.positive_rows,
        spec.negative_rows,
        spec.inconclusive_rows,
    ):
        raise AssertionError(
            f"{row_path.name}: counts {(row_count, positive, negative, inconclusive)} "
            f"!= {(spec.rows, spec.positive_rows, spec.negative_rows, spec.inconclusive_rows)}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=Path("work/rh_compute/results"))
    args = parser.parse_args()

    for spec in SPECS:
        validate_spec(args.results_dir, spec)
        print(f"OK Edrei quadrature scout: {spec.name}")
    print("validated 1 positive Arb recurrence scout and 1 inconclusive frontier scout")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
