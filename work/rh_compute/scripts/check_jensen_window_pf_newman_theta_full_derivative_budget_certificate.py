#!/usr/bin/env python3
"""Independently validate the full finite-N derivative budgets."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

import flint
from flint import arb


STEM = "jensen_window_pf_newman_theta_full_derivative_budget_certificate"
DEFAULT_ARTIFACT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
COMPACT_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_stable_remainder_arb_quadratic_matrix.json"
)
OUTER_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate.json"
)
CONTRACT_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_adaptive_modular_c1_remainder_contract.json"
)
PRECISION_BITS = 256
N_VALUES = tuple(range(4, 11))
EXPECTED_IDS = [
    "ntfdbc_01_exact_domain_split",
    "ntfdbc_02_compact_full_arithmetic",
    "ntfdbc_03_outer_full_arithmetic",
    "ntfdbc_04_full_m9_budgets",
    "ntfdbc_05_direct_c1_errors",
    "ntfdbc_06_finite_monotone_calibration",
    "ntfdbc_07_cofinal_separation_handoff",
]


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def serialize(value: arb) -> dict:
    log10_value = value.log() / arb(10).log()
    return {
        "enclosure": value.str(65, more=True),
        "upper": value.upper().str(65),
        "log10_enclosure": log10_value.str(50, more=True),
        "relative_accuracy_bits": value.rel_accuracy_bits(),
    }


def rebuild_budgets() -> list[dict]:
    compact = json.loads(COMPACT_SOURCE.read_text(encoding="utf-8"))
    outer = json.loads(OUTER_SOURCE.read_text(encoding="utf-8"))
    compact_by_n = {
        int(row["N"]): row for row in compact["compact_rows"]
    }
    outer_by_n = {
        int(row["N"]): row for row in outer["outer_aggregates"]
    }
    x = arb(38)
    rows: list[dict] = []
    previous_d0: arb | None = None
    previous_d1: arb | None = None
    for retained in N_VALUES:
        compact_row = compact_by_n[retained]
        outer_row = outer_by_n[retained]
        compact_d0 = arb(compact_row["compact_full_d0"]["upper"])
        compact_d1 = arb(compact_row["compact_full_d1"]["upper"])
        outer_d0 = arb(outer_row["outer_d0"]["upper"])
        outer_d1 = arb(outer_row["outer_d1"]["upper"])
        d0 = compact_d0 + outer_d0
        d1 = compact_d1 + outer_d1
        if previous_d0 is not None and d0.upper() >= previous_d0.lower():
            raise RuntimeError("independent d0 monotonicity check failed")
        if previous_d1 is not None and d1.upper() >= previous_d1.lower():
            raise RuntimeError("independent d1 monotonicity check failed")
        previous_d0, previous_d1 = d0, d1
        rows.append(
            {
                "N": retained,
                "compact_d0": serialize(compact_d0),
                "outer_d0": serialize(outer_d0),
                "full_d0_m9_t1_5": serialize(d0),
                "compact_d1": serialize(compact_d1),
                "outer_d1": serialize(outer_d1),
                "full_d1_m9_t1_5": serialize(d1),
                "x38_J_error_upper": serialize(16 * d0 / x**5),
                "x38_J_prime_error_upper": serialize(
                    64 * d0 / x**6 + 16 * d1 / x**5
                ),
                "uniform_frequency_statement": (
                    "|J-J_N|<=16*d0/x^5; "
                    "|J'-J_N'|<=64*d0/x^6+16*d1/x^5, x>=38"
                ),
            }
        )
    return rows


def validate(path: Path) -> list[str]:
    flint.ctx.prec = PRECISION_BITS
    issues: list[str] = []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("kind") != STEM:
        issues.append("artifact kind mismatch")
    expected_parameters = {
        "precision_bits": PRECISION_BITS,
        "m_order": 9,
        "t_cap_exact": "1/5",
        "u_split_exact": "11/5",
        "n_values": list(N_VALUES),
        "frequency_calibration": 38,
    }
    if artifact.get("parameters") != expected_parameters:
        issues.append("parameter mismatch")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if [row.get("readiness") for row in rows] != (
        ["proved"] * 6 + ["not_ready_to_apply"]
    ):
        issues.append("row readiness mismatch")

    compact = json.loads(COMPACT_SOURCE.read_text(encoding="utf-8"))
    outer = json.loads(OUTER_SOURCE.read_text(encoding="utf-8"))
    contract = json.loads(CONTRACT_SOURCE.read_text(encoding="utf-8"))
    if compact.get("cache", {}).get("records") != 223:
        issues.append("compact cache row count mismatch")
    if not compact.get("cache", {}).get("complete"):
        issues.append("compact cache is not complete")
    if len(outer.get("outer_entries", [])) != 203:
        issues.append("outer entry count mismatch")

    expected_contract = (
        "d_(N,j,m)(T)=(1/2)*sup_(0<=t<=T)"
        "||partial_u^m[(iu)^j*exp(tu^2)r_N(u)]||_1"
    )
    stored_contract = (
        contract.get("exact", {})
        .get("tail_budgets", {})
        .get("derivative")
    )
    if stored_contract != expected_contract:
        issues.append("upstream derivative contract mismatch")
    audit = artifact.get("source_audit", {})
    expected_audit = {
        "compact_source_sha256": file_hash(COMPACT_SOURCE),
        "outer_source_sha256": file_hash(OUTER_SOURCE),
        "contract_source_sha256": file_hash(CONTRACT_SOURCE),
        "contract_formula": expected_contract,
    }
    if audit != expected_audit:
        issues.append("source audit mismatch")

    try:
        rebuilt = rebuild_budgets()
    except Exception as exc:
        issues.append(f"budget reconstruction failed: {exc}")
        rebuilt = []
    if artifact.get("budgets") != rebuilt:
        issues.append("independent full-budget reconstruction mismatch")

    for row in artifact.get("budgets", []):
        d0 = arb(row["full_d0_m9_t1_5"]["enclosure"])
        d1 = arb(row["full_d1_m9_t1_5"]["enclosure"])
        x = arb(38)
        expected_j = 16 * d0 / x**5
        expected_j_prime = 64 * d0 / x**6 + 16 * d1 / x**5
        stored_j = arb(row["x38_J_error_upper"]["enclosure"])
        stored_j_prime = arb(
            row["x38_J_prime_error_upper"]["enclosure"]
        )
        if not stored_j.overlaps(expected_j):
            issues.append(f"J error algebra mismatch at N={row['N']}")
        if not stored_j_prime.overlaps(expected_j_prime):
            issues.append(f"J' error algebra mismatch at N={row['N']}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "N=4..10",
        "first-jet",
        "unbounded adaptive N",
        "transition",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary guard missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman theta full derivative-budget certificate: "
        f"{len(artifact['rows'])} rows, "
        f"{len(artifact['budgets'])} full d0/d1 budgets, 0 issues, "
        "1 open cofinal first-jet handoff"
    )


if __name__ == "__main__":
    main()
