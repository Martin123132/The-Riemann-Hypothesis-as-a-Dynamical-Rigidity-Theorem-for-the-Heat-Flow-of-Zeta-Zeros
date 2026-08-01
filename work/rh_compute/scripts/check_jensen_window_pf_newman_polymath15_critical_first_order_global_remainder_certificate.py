#!/usr/bin/env python3
"""Validate the critical first-order global remainder certificate."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_critical_"
    "first_order_global_remainder_certificate"
)
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_builder():
    spec = importlib.util.spec_from_file_location("first_order_global", BUILDER)
    if spec is None or spec.loader is None:
        raise ImportError("cannot load first-order global builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    if not BUILDER.is_file() or not RESULT.is_file() or not NOTE.is_file():
        raise FileNotFoundError("builder, result, or note missing")
    builder = load_builder()
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    if artifact.get("kind") != STEM:
        raise RuntimeError("artifact kind mismatch")
    if artifact.get("date") != "2026-07-26":
        raise RuntimeError("artifact date mismatch")
    if artifact.get("builder_sha256") != file_hash(BUILDER):
        raise RuntimeError("builder hash mismatch")
    if artifact.get("source_sha256") != builder.source_hashes():
        raise RuntimeError("source hash drift")
    if artifact.get("exact") != builder.exact_identities():
        raise RuntimeError("exact dictionary drift")

    expected_summary = {
        "rows": 15,
        "exact_or_published_inputs": 4,
        "analytic_or_interval_bounds": 8,
        "proved_compositions": 4,
        "open_targets": 1,
        "eta_0": 100000,
        "eta_1": 200000,
    }
    if artifact.get("summary") != expected_summary:
        raise RuntimeError(f"summary mismatch: {artifact.get('summary')}")

    rows = artifact.get("rows", [])
    if len(rows) != 15:
        raise RuntimeError(f"expected 15 rows, found {len(rows)}")
    for index, row in enumerate(rows, start=1):
        if not row.get("id", "").startswith(f"np15f1grc_{index:02d}_"):
            raise RuntimeError(f"row order mismatch at {index}: {row.get('id')}")
    if sum(row.get("readiness") == "open" for row in rows) != 1:
        raise RuntimeError("open target accounting drifted")

    exact_text = "\n".join(str(value) for value in artifact["exact"].values())
    for marker in (
        "(C_1(p,sigma)-C_1(p,1/2))/a",
        "K_u=max(floor(-u)+3,floor(T_0/pi))",
        "2*10^-30/a^14",
        "30000*exp(-5L/4)",
        "100000*exp(-5L/4)",
        "200000*L*exp(-5L/4)",
        "strict Xi arithmetic reversal",
    ):
        if marker not in exact_text:
            raise RuntimeError(f"exact marker missing: {marker}")

    interval = artifact.get("interval", {})
    expected_numeric_fields = {
        "expanded_C0_ball": 5,
        "C0_prefactor_scaled_ball": 4,
        "C1_cross_moment_ball": 500,
        "C1_prefactor_scaled_ball": 1000,
        "positive_a2_constant_ball": 1,
        "negative_low_series_ball": 0.02,
        "tail_scaled_constant_ball": 10,
        "pair_lift_scaled_ball": 1000,
        "endpoint_pair_bracket_ball": 4000,
        "endpoint_normalized_ball": 30000,
        "fixed_raw_ball": 40000,
        "adjacent_total_ball": 20000,
        "global_raw_ball": 100000,
    }
    for key, upper in expected_numeric_fields.items():
        text = interval.get(key)
        if not isinstance(text, str):
            raise RuntimeError(f"interval field missing: {key}")
        midpoint = float(text.split("+/-")[0].strip(" []"))
        if midpoint >= upper:
            raise RuntimeError(f"interval guard failed: {key}={text}")
    if interval.get("derivative_bounds") != {
        "F": 5,
        "F1": 50,
        "F2": 1000,
        "F3": 30000,
        "F4": 1200000,
    }:
        raise RuntimeError("C0 derivative bounds drifted")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "strict Xi contact-normal inequality",
        "bounded-L shoulder",
        "contact exclusion",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            raise RuntimeError(f"proof-boundary marker missing: {marker}")

    note = NOTE.read_text(encoding="utf-8")
    for marker in (
        "# Newman Polymath-15 Critical First-Order Global Remainder Certificate",
        "Exact First-Order Cancellation",
        "Fourth-Derivative Strip",
        "Positive Shift",
        "Negative Shift",
        "Holomorphic Lift",
        "Adjacent Cutoffs",
        "Global First Jet",
        "Remaining Theorem",
    ):
        if marker not in note:
            raise RuntimeError(f"note marker missing: {marker}")

    print(
        "validated Newman critical first-order global remainder certificate: "
        "15 rows, eta_0=100000, eta_1=200000, "
        "1 open signed-contact target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
