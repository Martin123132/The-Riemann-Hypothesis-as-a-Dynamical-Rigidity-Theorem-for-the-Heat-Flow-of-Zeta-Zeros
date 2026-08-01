#!/usr/bin/env python3
"""Validate the first-order cofinal boundary reduction."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_cofinal_boundary_reduction"
)
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_builder():
    spec = importlib.util.spec_from_file_location("cofinal_boundary", BUILDER)
    if spec is None or spec.loader is None:
        raise ImportError("cannot load cofinal-boundary builder")
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
    if artifact.get("audit") != builder.verify_exact_guards():
        raise RuntimeError("exact audit drift")
    if artifact.get("exact") != builder.exact_identities():
        raise RuntimeError("exact dictionary drift")

    expected_summary = {
        "rows": 13,
        "certified_inputs_or_compositions": 6,
        "exact_identities_or_bounds": 5,
        "open_boundary_targets": 3,
        "full_squared_budget": 50000000000,
        "half_squared_budget": 12500000000,
        "adjacent_squared_budget": 1300000000,
    }
    if artifact.get("summary") != expected_summary:
        raise RuntimeError(f"summary mismatch: {artifact.get('summary')}")
    rows = artifact.get("rows", [])
    if len(rows) != 13:
        raise RuntimeError(f"expected 13 rows, found {len(rows)}")
    for index, row in enumerate(rows, start=1):
        if not row.get("id", "").startswith(f"np15f1cbr_{index:02d}_"):
            raise RuntimeError(f"row order mismatch at {index}: {row.get('id')}")
    if sum(row.get("readiness") == "open" for row in rows) != 3:
        raise RuntimeError("open boundary target accounting drifted")

    exact_text = "\n".join(artifact["exact"].values())
    for marker in (
        "50000000000*exp(-5L/2)",
        "12500000000*exp(-5L/2)",
        "(3125/2)*exp(-L)",
        "t_j=25/(100+j)",
        "q<1 iff L<sqrt((100+j)/50)",
        "1300000000*exp(-5L/2)",
        "(13/500)*50000000000",
        "join continuously without a separate arithmetic theorem",
        "wind(V_J)<1",
    ):
        if marker not in exact_text:
            raise RuntimeError(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "q>=1 boundary margin",
        "q<1",
        "one-sided proxy phase bound",
        "bounded-L shoulders",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            raise RuntimeError(f"proof-boundary marker missing: {marker}")

    note = NOTE.read_text(encoding="utf-8")
    for marker in (
        "# Newman Polymath-15 First-Order Cofinal Boundary Reduction",
        "First-Order Boundary Budget",
        "Ray-Aligned Boundary",
        "Bottom Split",
        "Continuity And Winding",
        "Proof Boundary",
    ):
        if marker not in note:
            raise RuntimeError(f"note marker missing: {marker}")

    print(
        "validated Newman first-order cofinal boundary reduction: "
        "13 rows, full budget 50000000000, "
        "3 open boundary targets"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
