#!/usr/bin/env python3
"""Independently check the equation-(10) global theta-current reduction."""

from __future__ import annotations

import hashlib
import json
import os
from decimal import Decimal
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import sympy as sp


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation10_global_theta_current_reduction_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation10_global_theta_current_reduction_gate.md"
SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/telemetry/zeta14cubicmult_telemetry.f90"
CHECKPOINT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/accumulation_weights/t1e10_crossblock/run/checkpoint.jsonl"
EXPECTED_BLOCKS = 36
M2 = 159577
FINAL_ALPHA = 5122421


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stage_zero() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with CHECKPOINT.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                row = json.loads(line, parse_float=str)
                if int(row["stage"]) == 0:
                    rows.append(row)
    return rows


def main() -> None:
    artifact = load(RESULT)
    require(artifact["passed"] is True, "artifact is not passed")
    require(artifact["scope"]["diagnostic_midpoint_used"] is False, "midpoint contamination")
    require(artifact["route_decision"]["continuous_derivative_alone_is_discrete_telescope"] is False, "telescope guard drift")
    require(artifact["route_decision"]["height_uniform_error_bound_proved"] is False, "uniform theorem overpromotion")

    alpha, j, x, y = sp.symbols("alpha j x y", real=True, nonzero=True)
    u = sp.pi * alpha * j * x / 2
    common = sp.exp(sp.I * sp.pi * (alpha**2 + j**2) * x / 4)
    direct_pair = (alpha + j) * sp.exp(sp.I * sp.pi * (alpha + j) ** 2 * x / 4) + (alpha - j) * sp.exp(sp.I * sp.pi * (alpha - j) ** 2 * x / 4)
    completed_pair = 2 * common * (alpha * sp.cos(u) + sp.I * j * sp.sin(u))
    require(sp.simplify(sp.expand_complex(direct_pair - completed_pair)) == 0, "independent direct-pair identity failed")

    potential = common * sp.cos(u)
    require(sp.simplify(4 * sp.diff(potential, alpha) / (sp.I * sp.pi * x) - completed_pair) == 0, "independent current identity failed")
    theta_term = sp.exp(sp.I * sp.pi * alpha**2 * x / 4 + sp.I * sp.pi * alpha * y / 2)
    require(sp.simplify(2 * sp.diff(theta_term, y) / (sp.I * sp.pi) - alpha * theta_term) == 0, "independent theta identity failed")

    rows = stage_zero()
    require(len(rows) == EXPECTED_BLOCKS, "checkpoint block count drift")
    total_alpha = 101
    pair_collections = 0
    gaussian_sums = 0
    for block in range(1, EXPECTED_BLOCKS):
        previous = int(Decimal(rows[block - 1]["rn1"]))
        endpoint = int(Decimal(rows[block]["rn1"]))
        mt = int(rows[block]["mt"])
        count = int(rows[block]["aenums"])
        require(mt % 2 == 1 and count > 0, f"block metadata drift at {block}")
        require(endpoint == previous + 2 * (mt + 1) * count, f"endpoint recurrence failed at {block}")
        first_pivot = previous + 2 + mt
        last_pivot = first_pivot + 2 * (mt + 1) * (count - 1)
        require(first_pivot - mt == previous + 2, f"first paired endpoint failed at {block}")
        require(last_pivot + mt == endpoint, f"last paired endpoint failed at {block}")
        total_alpha += (mt + 1) * count
        pair_collections += count
        gaussian_sums += 2 * count
        saved = artifact["block_rows"][block]
        require(saved["alpha_start"] == previous + 2 and saved["alpha_end"] == endpoint, f"saved block endpoints drift at {block}")
        require(saved["alpha_count"] == (mt + 1) * count, f"saved alpha count drift at {block}")

    require(total_alpha == (FINAL_ALPHA - M2) // 2 + 1, "global roster formula failed")
    require(artifact["aggregate"]["total_alpha_count"] == total_alpha, "saved total alpha drift")
    require(artifact["aggregate"]["pair_collection_count"] == pair_collections, "saved pair count drift")
    require(artifact["aggregate"]["gaussian_phase_sum_count"] == gaussian_sums, "saved Gaussian count drift")

    source_lines = SOURCE.read_text(encoding="utf-8", errors="replace").splitlines()
    for audit in artifact["source_audit"]:
        excerpt = "\n".join(source_lines[int(audit["line_start"]) - 1 : int(audit["line_end"])])
        require(all(token in excerpt for token in audit["required_tokens"]), f"source token drift: {audit['id']}")
        require(hashlib.sha256(excerpt.encode("utf-8")).hexdigest() == audit["excerpt_sha256"], f"source hash drift: {audit['id']}")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {dependency['path']}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {dependency['path']}")
    for source in artifact["sources"].values():
        path = REPO_ROOT / source["path"]
        require(file_hash(path) == source["sha256"], f"source hash drift: {source['path']}")

    note = NOTE.read_text(encoding="utf-8")
    for token in (
        "omitted sine term",
        "incomplete theta",
        "not automatically an endpoint telescope",
        "No finite-Poisson remainder",
    ):
        require(token in note, f"note token missing: {token}")

    print(
        "validated equation-(10) global theta-current reduction independently: "
        f"blocks={EXPECTED_BLOCKS}, alpha={total_alpha}, pairs={pair_collections}"
    )


if __name__ == "__main__":
    main()
