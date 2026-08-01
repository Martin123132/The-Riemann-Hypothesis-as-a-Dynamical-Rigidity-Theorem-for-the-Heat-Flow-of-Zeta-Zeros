#!/usr/bin/env python3
"""Validate fixed order-eleven signed-Hankel completion at lambda zero."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_lambda0_completion_certificate as certificate  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=certificate.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=certificate.DEFAULT_NOTE)
    args = parser.parse_args()
    issues: list[str] = []
    if not args.artifact.exists():
        print(f"missing artifact: {args.artifact}")
        return 1
    try:
        actual = json.loads(args.artifact.read_text(encoding="utf-8"))
        expected = certificate.build_artifact()
    except (OSError, json.JSONDecodeError, RuntimeError, KeyError, ValueError) as exc:
        print(f"order-eleven lambda-zero completion: source failure: {exc}")
        return 1
    if actual != expected:
        issues.append("artifact differs from deterministic live-source rebuild")
    if (
        actual.get("kind")
        != "jensen_window_pf_compound_order11_lambda0_completion_certificate"
        or actual.get("status")
        != (
            "rigorous all-shift signed Hankel order eleven and fixed-order "
            "arbitrary-column completion at lambda zero"
        )
        or actual.get("generator") != certificate.GENERATOR_PATH
        or actual.get("checker") != certificate.CHECKER_PATH
    ):
        issues.append("artifact identity or status changed")
    exact = actual.get("exact", {})
    expected_exact = {
        "endpoint": certificate.ENDPOINT,
        "delayed_order11_heat_ray": certificate.ORDER11_HEAT,
        "lambda0_prefix": certificate.PREFIX,
        "all_shift_order11_lambda0": certificate.ALL_SHIFT_Q11,
        "contiguous_through_order11_lambda0": certificate.CONTIGUOUS_THROUGH_11,
        "nonpromotion_guard": certificate.NONPROMOTION_GUARD,
        "arbitrary_columns_through_order11_lambda0": certificate.ARBITRARY_THROUGH_11,
    }
    if any(exact.get(key) != value for key, value in expected_exact.items()):
        issues.append("completion theorem chain changed")
    rows = actual.get("rows", [])
    if (
        len(rows) != 8
        or any(row.get("readiness") != "ready_to_apply" for row in rows)
        or rows[-1].get("formula") != certificate.ARBITRARY_THROUGH_11
    ):
        issues.append("completion ledger changed")
    summary = actual.get("summary", {})
    if (
        summary.get("open_rows") != 0
        or summary.get("delayed_order11_heat_theorems") != 1
        or summary.get("all_shift_order11_lambda0_theorems") != 1
        or summary.get("contiguous_through_order11_lambda0_theorems") != 1
        or summary.get("arbitrary_column_through_order11_lambda0_theorems") != 1
        or summary.get("countermodel_nonpromotion_guards") != 1
        or summary.get("orders_above_11") != 0
        or summary.get("pf_infinity_theorems") != 0
        or summary.get("rh_claims") != 0
    ):
        issues.append("claim-boundary summary changed")
    sources = actual.get("sources", [])
    if len(sources) != 5 or any(len(row.get("sha256", "")) != 64 for row in sources):
        issues.append("five-source hash contract changed")
    if not args.note.exists():
        issues.append("missing note")
    else:
        note = args.note.read_text(encoding="utf-8")
        for marker in (certificate.ALL_SHIFT_Q11, "not a proof of RH"):
            if marker not in note:
                issues.append(f"note missing marker: {marker}")
    if issues:
        print(f"order-eleven lambda-zero completion: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated fixed order-eleven lambda-zero completion: "
        "1 all-shift Q11 theorem, 1 arbitrary-column-through-11 theorem, "
        "0 PF-infinity claims, 0 RH claims"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
