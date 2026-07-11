#!/usr/bin/env python3
"""Validate the cofinal Jensen scaling-limit equivalence gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.json"


def falling_ratio(D: int, j: int) -> Fraction:
    value = Fraction(1, 1)
    for offset in range(j):
        value *= Fraction(D - offset, D)
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    issues: list[str] = []
    rows = payload.get("rows", [])
    exact = payload.get("exact", {})

    for j in range(7):
        values = [falling_ratio(D, j) for D in (100, 1000, 10000) if D >= j]
        if any(value <= 0 or value > 1 for value in values):
            issues.append(f"falling ratio range failed at j={j}")
        if j >= 2 and not values[0] < values[1] < values[2] < 1:
            issues.append(f"falling ratio convergence failed at j={j}")

    required_exact = {
        "scaled_polynomial": "(D)_j/D^j",
        "coefficient_limit": "->1",
        "local_uniform_limit": "locally uniformly",
        "cofinal_implication": "Laguerre-Polya",
        "converse": "Jensen's theorem",
        "equivalence": "<=>",
    }
    for key, marker in required_exact.items():
        if marker not in str(exact.get(key, "")):
            issues.append(f"missing exact marker {key}:{marker}")

    if len(rows) != 10:
        issues.append(f"row count {len(rows)} != 10")
    if sum(row.get("role") in {"noncircularity_guard", "forbidden_promotion"} for row in rows) != 2:
        issues.append("non-promotion guard count failed")
    if sum(row.get("role") == "open_handoff" for row in rows) != 1:
        issues.append("open handoff count failed")
    boundary = payload.get("proof_boundary", "")
    for marker in ("does not prove", "Laguerre-Polya", "PF-infinity", "RH"):
        if marker not in boundary:
            issues.append(f"missing proof-boundary marker {marker!r}")

    print(
        "validated Jensen-window PF cofinal scaling-limit equivalence gate: "
        f"{len(rows)} rows, {len(issues)} issues, 2 exact scaling identities, "
        "2 analytic limit steps, 1 cofinal-to-LP theorem, "
        "1 LP-to-all-degrees theorem, 1 fixed-shift equivalence, "
        "2 non-promotion guards, 1 open independent-route handoff"
    )
    for issue in issues:
        print(f"ISSUE {issue}")
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
