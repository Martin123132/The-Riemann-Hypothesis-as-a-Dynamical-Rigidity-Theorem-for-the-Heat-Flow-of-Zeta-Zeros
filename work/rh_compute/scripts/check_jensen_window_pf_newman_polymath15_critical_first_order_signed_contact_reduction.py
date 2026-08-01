#!/usr/bin/env python3
"""Validate the critical first-order signed-contact reduction."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_critical_"
    "first_order_signed_contact_reduction"
)
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_builder():
    spec = importlib.util.spec_from_file_location("signed_contact", BUILDER)
    if spec is None or spec.loader is None:
        raise ImportError("cannot load signed-contact builder")
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
    builder.verify_symbolic_identities()
    if artifact.get("exact") != builder.exact_identities():
        raise RuntimeError("exact dictionary drift")

    expected_summary = {
        "rows": 15,
        "certified_inputs_or_compositions": 3,
        "exact_definitions_or_identities": 10,
        "route_guards": 1,
        "open_frequency_targets": 1,
        "eta_0": 100000,
        "eta_1": 200000,
    }
    if artifact.get("summary") != expected_summary:
        raise RuntimeError(f"summary mismatch: {artifact.get('summary')}")
    rows = artifact.get("rows", [])
    if len(rows) != 15:
        raise RuntimeError(f"expected 15 rows, found {len(rows)}")
    for index, row in enumerate(rows, start=1):
        if not row.get("id", "").startswith(f"np15f1scr_{index:02d}_"):
            raise RuntimeError(f"row order mismatch at {index}: {row.get('id')}")
    if sum(row.get("readiness") == "open" for row in rows) != 1:
        raise RuntimeError("open frequency target accounting drifted")

    exact_text = "\n".join(artifact["exact"].values())
    for marker in (
        "f_n=e_n*(1+d_n)",
        "F'''(p)/(12*pi^2*a)",
        "J_[1]=2*Re(E_[1])",
        "|X_[1]|<50000*exp(-5L/4)",
        "|U_[1]|<100000*L*exp(-5L/4)",
        "sum_(n=1)^N log(n/a)f_n",
        "q>=1",
        "endpoint simplicity",
    ):
        if marker not in exact_text:
            raise RuntimeError(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "strict Xi crossing-band inequality",
        "q<1",
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
        "# Newman Polymath-15 Critical First-Order Signed-Contact Reduction",
        "Corrected Complex Main",
        "Rate-Free Contact Identity",
        "Certified Contact Box",
        "Saddle-Centered Observable",
        "Domain Split",
        "Proof Boundary",
    ):
        if marker not in note:
            raise RuntimeError(f"note marker missing: {marker}")

    print(
        "validated Newman critical first-order signed-contact reduction: "
        "15 rows, eta_0=100000, eta_1=200000, "
        "1 open frequency target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
