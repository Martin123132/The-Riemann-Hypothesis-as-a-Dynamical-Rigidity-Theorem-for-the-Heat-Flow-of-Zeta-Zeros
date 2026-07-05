#!/usr/bin/env python3
"""Validate promoted signed-Hankel finite certificate summaries."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path


LAMBDAS: tuple[str, ...] = ("0", "1e-6", "1e-4", "1e-2", "1e-1")


@dataclass(frozen=True)
class HankelSpec:
    lam: str
    summary_name: str
    jsonl_name: str
    rows: int
    m_min: int
    m_max: int
    shift_min: int
    shift_max: int
    needed_max_k: int


def label_for(lam: str) -> str:
    if lam == "0":
        return "lam0"
    return "lam" + lam.replace("-", "m").replace("+", "").replace(".", "p")


HANKEL_SPECS: tuple[HankelSpec, ...] = tuple(
    HankelSpec(
        lam=lam,
        summary_name=f"arb_hankel_enclosure_{label_for(lam)}_m0_m19_s0_s24_dps520_summary.json",
        jsonl_name=f"arb_hankel_enclosure_{label_for(lam)}_m0_m19_s0_s24_dps520.jsonl",
        rows=500,
        m_min=0,
        m_max=19,
        shift_min=0,
        shift_max=24,
        needed_max_k=62,
    )
    for lam in LAMBDAS
)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_spec(results_dir: Path, spec: HankelSpec) -> None:
    summary_path = results_dir / spec.summary_name
    row_path = results_dir / spec.jsonl_name
    summary = load_json(summary_path)
    checks = {
        "kind": "arb_hankel_enclosure_sign_probe_summary",
        "lam": spec.lam,
        "rows": spec.rows,
        "ok": spec.rows,
        "failed_or_inconclusive": 0,
        "m_min": spec.m_min,
        "m_max": spec.m_max,
        "shift_min": spec.shift_min,
        "shift_max": spec.shift_max,
        "needed_max_k": spec.needed_max_k,
        "all_ok": True,
    }
    for key, expected in checks.items():
        actual = summary.get(key)
        if actual != expected:
            raise AssertionError(f"{summary_path.name}: {key} {actual!r} != {expected!r}")
    counts = summary.get("counts", {})
    if counts.get("positive") != spec.rows or len(counts) != 1:
        raise AssertionError(f"{summary_path.name}: unexpected counts {counts!r}")

    if not row_path.exists():
        raise FileNotFoundError(f"missing row log {row_path.name}")
    with row_path.open("r", encoding="utf-8") as handle:
        row_count = 0
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            row_count += 1
            if not row.get("ok"):
                raise AssertionError(f"{row_path.name}: non-ok row {row!r}")
            if row.get("classification") != "positive":
                raise AssertionError(f"{row_path.name}: unexpected classification {row!r}")
            if row.get("contains_zero"):
                raise AssertionError(f"{row_path.name}: row contains zero {row!r}")
        if row_count != spec.rows:
            raise AssertionError(f"{row_path.name}: row count {row_count} != {spec.rows}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", default="work/rh_compute/results", type=Path)
    args = parser.parse_args()

    for spec in HANKEL_SPECS:
        validate_spec(args.results_dir, spec)
        print(f"OK signed-Hankel: lambda={spec.lam} :: {spec.summary_name}")
    print(f"validated {sum(spec.rows for spec in HANKEL_SPECS)} signed-Hankel finite certificates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
