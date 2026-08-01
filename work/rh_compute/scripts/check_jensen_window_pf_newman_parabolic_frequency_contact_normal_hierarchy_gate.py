#!/usr/bin/env python3
"""Validate the parabolic-frequency contact-normal hierarchy gate."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_parabolic_frequency_contact_normal_hierarchy_gate"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
HISTORICAL_SCALED_SUCCESSOR_SHA256 = (
    "7a165fd517ef6442a9a8f4fcde8397bfb443c3ad04001d5d7a0d06220e4842a5"
)
HISTORICAL_RAY_SHA256 = (
    "2f2c7e69957446d461354e673ff27af66c346863e9264cf2d573be0b470c10e5"
)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_builder():
    spec = importlib.util.spec_from_file_location("contact_normal_builder", BUILDER)
    if spec is None or spec.loader is None:
        raise ImportError("cannot load contact-normal builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    if not BUILDER.is_file() or not RESULT.is_file() or not NOTE.is_file():
        raise FileNotFoundError("builder, result, or note missing")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    builder = load_builder()

    expected_status = (
        "exact parabolic-frequency contact normal form, hierarchy route "
        "guards, and one open Xi C1 target"
    )
    if artifact.get("kind") != STEM:
        raise RuntimeError("artifact kind mismatch")
    if artifact.get("date") != "2026-07-26":
        raise RuntimeError("artifact date mismatch")
    if artifact.get("status") != expected_status:
        raise RuntimeError("artifact status mismatch")
    if artifact.get("builder_sha256") != file_hash(BUILDER):
        raise RuntimeError("builder hash mismatch")
    stored_sources = artifact.get("source_sha256", {})
    current_sources = builder.source_hashes()
    normalized_sources = dict(current_sources)
    for key, historical in (
        ("scaled_successor", HISTORICAL_SCALED_SUCCESSOR_SHA256),
        ("ray_aligned_reduction", HISTORICAL_RAY_SHA256),
    ):
        if stored_sources.get(key) not in {
            current_sources.get(key),
            historical,
        }:
            raise RuntimeError(f"unrecognized source snapshot: {key}")
        normalized_sources[key] = stored_sources.get(key)
    if stored_sources != normalized_sources:
        raise RuntimeError("source hash drift")
    scaled = json.loads(
        builder.SOURCE_FILES["scaled_successor"].read_text(encoding="utf-8")
    )
    ray = json.loads(
        builder.SOURCE_FILES["ray_aligned_reduction"].read_text(
            encoding="utf-8"
        )
    )
    if scaled.get("exact_summary", {}).get("all_j_xi_theorem") is not False:
        raise RuntimeError("current scaled-successor semantic guard drifted")
    if ray.get("summary", {}).get("new_strip_open_antecedents") != 0:
        raise RuntimeError("current ray new-strip semantic source drifted")
    if ray.get("summary", {}).get("old_collar_open_antecedents") != 1:
        raise RuntimeError("current ray old-collar semantic source drifted")
    if artifact.get("exact") != builder.build_exact():
        raise RuntimeError("stored exact identities drifted")

    rows = artifact.get("rows", [])
    if len(rows) != 17:
        raise RuntimeError(f"expected 17 rows, found {len(rows)}")
    expected_prefixes = [f"pfcn_{index:02d}_" for index in range(1, 18)]
    for prefix, row in zip(expected_prefixes, rows, strict=True):
        if not row.get("id", "").startswith(prefix):
            raise RuntimeError(f"row order mismatch: {row.get('id')}")
    readiness = [row.get("readiness") for row in rows]
    if readiness.count("available_exact") != 11:
        raise RuntimeError("exact-row accounting drifted")
    if readiness.count("closed_route") != 3:
        raise RuntimeError("route-guard accounting drifted")
    if readiness.count("closed_counterexample") != 1:
        raise RuntimeError("countermodel accounting drifted")
    if readiness.count("literature_nonpromotion") != 1:
        raise RuntimeError("literature-guard accounting drifted")
    if readiness.count("open") != 1:
        raise RuntimeError("open-target accounting drifted")

    expected_summary = {
        "rows": 17,
        "exact_identities": 12,
        "exact_route_guards": 3,
        "literature_guards": 1,
        "open_xi_targets": 1,
    }
    if artifact.get("summary") != expected_summary:
        raise RuntimeError(f"summary drifted: {artifact.get('summary')}")

    exact = artifact["exact"]
    required_exact = (
        "epsilon_0^2+beta^2*epsilon_1^2",
        "epsilon_0=2500*exp(-3L/4)",
        "|J_x|<=5000*L*exp(-3L/4)",
        "a contact requires |J|<=epsilon_0",
        "on |2X-Q|<=epsilon_0 prove",
        "|M(h_+-h_-)+D_0|>L*epsilon_1/2",
        "625*exp(-3L/4)*(4L+|h_++h_-|)",
        "M(h_+-h_-)+D_0=-r_x/2",
        "r=r_0+delta_0",
        "L*eta_1/2+(eta_0/4)|h_++h_-|",
        "full-wedge pointwise C1-box exclusion",
        "q<=1 iff t<=1/(2L^2)",
        "min(A0+A1/2,A0+A1)",
        "the s_pf,t H_x term vanishes",
        "F_1*F_3-F_2^2=-9*H_xx^4/16<0",
        "F_1F_3-F_2^2=-9",
    )
    exact_text = "\n".join(str(value) for value in exact.values())
    missing_exact = [text for text in required_exact if text not in exact_text]
    if missing_exact:
        raise RuntimeError(f"exact contract missing: {missing_exact}")

    note = NOTE.read_text(encoding="utf-8")
    required_note = (
        "# Newman Parabolic-Frequency Contact-Normal Hierarchy Gate",
        "Corrected Contact Normal Form",
        "Exact Chart Split",
        "Why Scale Alone Cannot Close Theorem",
        "Correlation-Hierarchy Translation",
        "Route Decision",
        "A strip-wide cone is no longer needed",
        "phase-matched equation",
        "peeling a signed leading remainder",
        "endpoint quantifier guard",
        "do not prove (2a)",
    )
    missing_note = [text for text in required_note if text not in note]
    if missing_note:
        raise RuntimeError(f"note contract missing: {missing_note}")

    print(
        "validated Newman parabolic-frequency contact-normal hierarchy gate: "
        "17 rows, 12 exact identities, 3 exact route guards, "
        "1 literature guard, 1 open Xi target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
