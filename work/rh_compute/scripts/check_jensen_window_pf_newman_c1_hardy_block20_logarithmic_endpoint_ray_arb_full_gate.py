#!/usr/bin/env python3
"""Validate the all-call logarithmic endpoint-ray Arb gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate.md"
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate.py"
)
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing log-ray full artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate",
        "log-ray full kind drift",
    )
    require(
        artifact["status"] == "rigorous_all_374_logarithmic_endpoint_ray_nonsaddle_enclosures_validated",
        "log-ray full status drift",
    )
    scope = artifact["scope"]
    aggregate = artifact["aggregate"]
    rows = artifact["rows"]
    require(scope["recursive_call_count"] == aggregate["call_count"] == len(rows) == 374, "log-ray full row count drift")
    require(scope["precision_ladder_decimal_digits"] == [70, 110], "log-ray full precision drift")
    require(scope["precision_overlap_count"] == 374, "log-ray full overlap count drift")
    require(scope["worker_count"] == 1, "log-ray full worker drift")
    require(sum(aggregate["cutoff_distribution"].values()) == 374, "log-ray full cutoff count drift")
    require(aggregate["maximum_ray_cutoff"] <= 65536, "log-ray full cutoff ceiling drift")
    require(aggregate["formula_target_identity_count"] == 374, "log-ray full target identity count drift")
    require(aggregate["correction_residual_identity_count"] == 374, "log-ray full residual identity count drift")
    require(aggregate["exact_correction_excluding_zero_count"] == 374, "log-ray full correction exclusion count drift")
    require(Decimal(aggregate["maximum_origin_radius_upper"]) < Decimal("1e-10"), "log-ray full origin bound drift")
    require(Decimal(aggregate["maximum_tail_radius_upper"]) < Decimal("1e-30"), "log-ray full tail bound drift")
    require(Decimal(aggregate["maximum_formula_target_gap_magnitude_upper"]) < Decimal("1e-10"), "log-ray full target gap drift")
    require(Decimal(aggregate["maximum_correction_residual_gap_magnitude_upper"]) < Decimal("1e-10"), "log-ray full residual gap drift")
    require(Decimal(aggregate["minimum_exact_correction_magnitude_lower"]) > Decimal("1e-7"), "log-ray full minimum correction drift")
    require(Decimal(aggregate["maximum_exact_correction_magnitude_upper"]) < Decimal("1e-2"), "log-ray full maximum correction drift")
    chain_roster = [row["chain"] for row in rows]
    require(chain_roster == sorted(chain_roster) and len(set(chain_roster)) == 374, "log-ray full row order drift")
    for row in rows:
        require(row["precision_overlap"], "log-ray full row precision drift")
        require(row["formula_target_identity_contains_zero"], "log-ray full target identity lost")
        require(row["correction_residual_identity_contains_zero"], "log-ray full residual identity lost")
        require(row["exact_minus_paper_excludes_zero"], "log-ray full exact correction lost")
        require(Decimal(row["tail_radius_upper"]) < Decimal("1e-30"), "log-ray full row tail drift")
        require(row["ray_cutoff"] in (64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536), "log-ray full row cutoff drift")
    for group in artifact["sources"].values():
        path = REPO_ROOT / group["path"]
        require(path.is_file(), f"log-ray full dependency missing: {group['path']}")
        require(file_hash(path) == group["sha256"], f"log-ray full dependency hash drift: {group['path']}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "log-ray full builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "log-ray full checker hash drift")
    for token in (
        "Status: rigorous all-374 finite-input logarithmic endpoint-ray enclosure",
        "Q_exact^- - (P^- - I69^-)",
        "374 / 374",
        "not a proof",
    ):
        require(token in note, f"log-ray full note token missing: {token}")
    require("one low-height block only" in artifact["proof_boundary"], "log-ray full boundary drift")
    print("validated logarithmic endpoint-ray full Arb gate: 374/374 calls, both exact identities enclose zero")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
