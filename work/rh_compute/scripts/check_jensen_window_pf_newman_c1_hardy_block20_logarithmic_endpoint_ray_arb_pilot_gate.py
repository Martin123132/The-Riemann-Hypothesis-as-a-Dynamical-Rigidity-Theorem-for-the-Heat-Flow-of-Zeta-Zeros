#!/usr/bin/env python3
"""Validate the four-branch logarithmic endpoint-ray Arb pilot."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate.md"
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate.py"
)
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing log-ray Arb pilot artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_pilot_gate",
        "log-ray Arb pilot kind drift",
    )
    require(
        artifact["status"] == "rigorous_four_branch_logarithmic_endpoint_ray_arb_pilot_validated",
        "log-ray Arb pilot status drift",
    )
    scope = artifact["scope"]
    aggregate = artifact["aggregate"]
    rows = artifact["rows"]
    require(scope["witness_chains"] == [1, 2, 3, 36], "log-ray Arb pilot witness drift")
    require(scope["precision_ladder_decimal_digits"] == [70, 110], "log-ray Arb pilot precision drift")
    require(scope["worker_count"] == 1, "log-ray Arb pilot worker drift")
    require(scope["precision_overlap_count"] == aggregate["witness_count"] == len(rows) == 4, "log-ray Arb pilot count drift")
    require(len({(row["phi1_sign"], row["phi3_sign"]) for row in rows}) == 4, "log-ray Arb pilot branch coverage drift")
    require(aggregate["formula_target_identity_count"] == 4, "log-ray Arb pilot target identity count drift")
    require(aggregate["correction_residual_identity_count"] == 4, "log-ray Arb pilot correction identity count drift")
    require(Decimal(aggregate["maximum_origin_radius_upper"]) < Decimal("1e-10"), "log-ray Arb pilot origin bound drift")
    require(Decimal(aggregate["maximum_tail_radius_upper"]) < Decimal("1e-30"), "log-ray Arb pilot tail bound drift")
    require(Decimal(aggregate["maximum_formula_target_gap_magnitude_upper"]) < Decimal("1e-10"), "log-ray Arb pilot target gap drift")
    require(Decimal(aggregate["maximum_correction_residual_gap_magnitude_upper"]) < Decimal("1e-10"), "log-ray Arb pilot correction gap drift")
    for row in rows:
        require(row["precision_overlap"], "log-ray Arb pilot row precision drift")
        require(row["formula_target_identity_contains_zero"], "log-ray Arb pilot target identity lost")
        require(row["correction_residual_identity_contains_zero"], "log-ray Arb pilot correction identity lost")
        require(row["b_angle"] in ("pi/2", "5*pi/6"), "log-ray Arb pilot b angle drift")
        require(row["c_angle"] in ("pi/2", "pi/6"), "log-ray Arb pilot c angle drift")
    for group in artifact["sources"].values():
        path = REPO_ROOT / group["path"]
        require(path.is_file(), f"log-ray Arb pilot dependency missing: {group['path']}")
        require(file_hash(path) == group["sha256"], f"log-ray Arb pilot dependency hash drift: {group['path']}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "log-ray Arb pilot builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "log-ray Arb pilot checker hash drift")
    for token in (
        "Status: rigorous four-branch finite-input Arb pilot",
        "K_m(w)=w^m*2F1",
        "Q_exact^- - (P^- - I69^-)",
        "not a proof",
    ):
        require(token in note, f"log-ray Arb pilot note token missing: {token}")
    require("four fixed low-height calls only" in artifact["proof_boundary"], "log-ray Arb pilot boundary drift")
    print("validated logarithmic endpoint-ray Arb pilot: 4/4 sign branches, both exact identities enclose zero")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
