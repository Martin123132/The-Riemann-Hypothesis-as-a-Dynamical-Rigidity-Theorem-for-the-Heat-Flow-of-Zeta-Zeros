#!/usr/bin/env python3
"""Validate the critical Dirichlet second-order remainder certificate."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_critical_"
    "dirichlet_second_order_remainder_certificate"
)
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_builder():
    spec = importlib.util.spec_from_file_location("dirichlet_second_order", BUILDER)
    if spec is None or spec.loader is None:
        raise ImportError("cannot load Dirichlet second-order builder")
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
    if artifact.get("exact") != builder.build_exact():
        raise RuntimeError("exact dictionary drift")

    expected_summary = {
        "rows": 12,
        "exact_or_published_reductions": 9,
        "proved_theorems": 2,
        "open_route_guards": 1,
        "per_term_constant": 400000,
        "fixed_cell_constant": 8000000,
    }
    if artifact.get("summary") != expected_summary:
        raise RuntimeError(f"summary mismatch: {artifact.get('summary')}")

    rows = artifact.get("rows", [])
    if len(rows) != 12:
        raise RuntimeError(f"expected 12 rows, found {len(rows)}")
    for index, row in enumerate(rows, start=1):
        if not row.get("id", "").startswith(f"np15cd2rc_{index:02d}_"):
            raise RuntimeError(f"row order mismatch at {index}: {row.get('id')}")
    if sum(row.get("readiness") == "ready_to_apply" for row in rows) != 2:
        raise RuntimeError("proved theorem accounting drifted")
    if sum(row.get("readiness") == "open" for row in rows) != 1:
        raise RuntimeError("open guard accounting drifted")

    audit = artifact.get("interval_audit", {})
    if audit.get("role") != (
        "high_precision_margin_audit_after_exact_rational_guards"
    ):
        raise RuntimeError("numeric audit was promoted above exact guards")
    rational_guards = audit.get("exact_rational_guards", {})
    for key in (
        "log_T_lt_14",
        "e_lt_2_72",
        "pi_gt_3_14",
        "exact_tail_log_upper",
        "mass_tail_log_upper",
        "second_tail_log_upper",
        "central_exponent_upper",
        "fixed_cell_constant_upper",
        "critical_T_min_guard",
    ):
        if key not in rational_guards:
            raise RuntimeError(f"exact rational guard missing: {key}")
    for key in (
        "exact_tail_log_margin_at_T_audit_min",
        "mass_tail_log_margin_at_T_audit_min",
        "second_tail_log_margin_at_T_audit_min",
        "exact_tail_log_derivative_at_T_audit_min",
        "mass_tail_log_derivative_at_T_audit_min",
        "second_tail_log_derivative_at_T_audit_min",
    ):
        if float(audit[key]) >= 0:
            raise RuntimeError(f"tail sign failed: {key}")
    if float(audit["central_exponent_upper_at_T_audit_min"]) >= 0.006:
        raise RuntimeError("central exponent audit failed")
    if float(audit["fixed_cell_constant_ball"]) >= 8_000_000:
        raise RuntimeError("fixed-cell constant audit failed")

    exact_text = "\n".join(str(value) for value in artifact["exact"].values())
    for marker in (
        "t|alpha_n|<=27",
        "|g+rho|<0.006",
        "312830/T^2",
        "400000/T^2",
        "8000000*exp(-7L/4)",
        "heat-integrated endpoint a^-2 bound",
    ):
        if marker not in exact_text:
            raise RuntimeError(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "Riemann-Siegel a^-2 endpoint remainder",
        "adjacent-cutoff signed lift",
        "contact exclusion",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            raise RuntimeError(f"proof-boundary marker missing: {marker}")

    note = NOTE.read_text(encoding="utf-8")
    for marker in (
        "# Newman Critical Dirichlet Second-Order Remainder Certificate",
        "Uniform Shift",
        "Central Range",
        "Gaussian Tail",
        "Per-Term Theorem",
        "Fixed-Cell Transfer",
        "Remaining Boundary",
    ):
        if marker not in note:
            raise RuntimeError(f"note marker missing: {marker}")

    print(
        "validated Newman critical Dirichlet second-order remainder "
        "certificate: 12 rows, C_D=400000, fixed-cell C1 constant "
        "8000000, 1 open global splice guard"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
