#!/usr/bin/env python3
"""Validate the global order-eleven lambda=-100 entry certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_m100_entry_certificate as certificate  # noqa: E402


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
        print(f"order-eleven m100 entry certificate: source failure: {exc}")
        return 1
    if actual != expected:
        issues.append("artifact differs from deterministic live-source rebuild")
    if (
        actual.get("kind")
        != "jensen_window_pf_compound_order11_m100_entry_certificate"
        or actual.get("status")
        != "rigorous all-shift signed order-eleven entry at lambda=-100"
        or actual.get("generator") != certificate.GENERATOR_PATH
        or actual.get("checker") != certificate.CHECKER_PATH
    ):
        issues.append("artifact identity or status changed")
    exact = actual.get("exact", {})
    if (
        exact.get("order11_continuous") != certificate.ORDER11_CONTINUOUS
        or exact.get("analytic_tail") != certificate.TAIL_THEOREM
        or exact.get("finite_prefix") != certificate.FINITE_THEOREM
        or exact.get("global_endpoint") != certificate.GLOBAL_ENDPOINT
    ):
        issues.append("endpoint theorem chain changed")
    rows = actual.get("rows", [])
    if (
        len(rows) != 7
        or any(row.get("readiness") != "ready_to_apply" for row in rows)
        or rows[-1].get("formula") != certificate.GLOBAL_ENDPOINT
    ):
        issues.append("entry ledger changed")
    summary = actual.get("summary", {})
    if (
        summary.get("open_rows") != 0
        or summary.get("discharged_continuum_targets") != 1
        or summary.get("circularity_guards") != 1
        or summary.get("full_kernel_transfer_applications") != 1
        or summary.get("global_m100_order11_entry_theorems") != 1
        or summary.get("heat_interval_theorems") != 0
        or summary.get("orders_above_11") != 0
        or summary.get("rh_claims") != 0
    ):
        issues.append("claim-boundary summary changed")
    sources = actual.get("source_contract", {}).get("sources", [])
    if len(sources) != 4 or any(len(row.get("sha256", "")) != 64 for row in sources):
        issues.append("four-source hash contract changed")
    if not args.note.exists():
        issues.append("missing note")
    else:
        note = args.note.read_text(encoding="utf-8")
        for marker in (certificate.GLOBAL_ENDPOINT, "not a proof of RH"):
            if marker not in note:
                issues.append(f"note missing marker: {marker}")
    if issues:
        print(f"order-eleven m100 entry certificate: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated global order-eleven lambda=-100 entry: "
        "1 all-shift fixed-order endpoint theorem, 0 heat claims, 0 RH claims"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
