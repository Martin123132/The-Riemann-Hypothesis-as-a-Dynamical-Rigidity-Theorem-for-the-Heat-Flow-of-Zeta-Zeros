#!/usr/bin/env python3
"""Check the ray-aligned exhaustion and parabolic-frequency reduction."""

from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_ray_aligned_parabolic_frequency_reduction"
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts" / f"{STEM}.py"


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_builder():
    spec = importlib.util.spec_from_file_location("ray_scale_builder", BUILDER)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def check_symbolic_identities() -> None:
    t, ell = sp.symbols("t ell", positive=True)
    q = 2 * t * ell**2
    scale = sp.sqrt(2 * t / (1 + q))
    alpha = sp.simplify(scale / sp.sqrt(2 * t))
    beta = sp.simplify(ell * scale)
    checks = (
        sp.simplify(scale ** -2 - ell**2 - 1 / (2 * t)),
        sp.simplify(alpha**2 + beta**2 - 1),
        sp.simplify(
            sp.diff(sp.log(scale), t) - 1 / (2 * t * (1 + q))
        ),
        sp.simplify(
            sp.diff(sp.log(alpha), t) + ell**2 / (1 + q)
        ),
        sp.simplify(
            sp.diff(sp.log(beta), t) - 1 / (2 * t * (1 + q))
        ),
    )
    if any(value != 0 for value in checks):
        raise RuntimeError(f"symbolic scale failure: {checks}")


def check_schedule() -> None:
    c = Fraction(25)
    a = 100
    previous_time: Fraction | None = None
    for stage in range(10_001):
        time = c / (a + stage)
        next_time = c / (a + stage + 1)
        entry_log = Fraction(a + stage + 1)
        if next_time * entry_log != c:
            raise RuntimeError(f"new-strip ray failure at stage {stage}")
        expected_delta = c / ((a + stage) * (a + stage + 1))
        if time - next_time != expected_delta:
            raise RuntimeError(f"time-step identity failure at stage {stage}")
        if previous_time is not None and time >= previous_time:
            raise RuntimeError(f"time schedule is not decreasing at {stage}")
        if entry_log < 101:
            raise RuntimeError(f"log radius below outer threshold at {stage}")
        entry_q = 2 * next_time * entry_log**2
        if entry_q < 5050:
            raise RuntimeError(f"entry chart is not frequency-dominant at {stage}")
        previous_time = time


def main() -> int:
    if not RESULT.is_file() or not NOTE.is_file() or not BUILDER.is_file():
        raise FileNotFoundError("builder output, note, or builder is missing")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    builder = load_builder()

    expected_status = (
        "exact ray-aligned exhaustion and parabolic-frequency scale "
        "reduction with one open Xi descendant antecedent"
    )
    if artifact.get("kind") != STEM:
        raise RuntimeError("artifact kind mismatch")
    if artifact.get("date") != "2026-07-26":
        raise RuntimeError("artifact date mismatch")
    if artifact.get("status") != expected_status:
        raise RuntimeError("artifact status mismatch")
    if artifact.get("builder_sha256") != file_hash(BUILDER):
        raise RuntimeError("builder hash mismatch")
    if artifact.get("exact") != builder.build_exact():
        raise RuntimeError("stored exact identities drifted")
    if artifact.get("source_sha256") != builder.source_hashes():
        raise RuntimeError("source hashes drifted")

    constants = artifact.get("constants", {})
    if constants != {
        "c_star": "4911678521/1933561194",
        "explicit_ray": "25",
        "schedule_offset": 100,
        "top_time": "1/2",
        "base_time": "1/4",
        "minimum_log_radius": 101,
    }:
        raise RuntimeError("constant contract drifted")

    rows = artifact.get("rows", [])
    expected_ids = [f"rapf_{index:02d}_" for index in range(1, 15)]
    if len(rows) != 14:
        raise RuntimeError(f"expected 14 rows, found {len(rows)}")
    for prefix, row in zip(expected_ids, rows, strict=True):
        if not row.get("id", "").startswith(prefix):
            raise RuntimeError(f"row ordering mismatch: {row.get('id')}")
    readiness = [row["readiness"] for row in rows]
    if readiness.count("open") != 1 or readiness.count("open_antecedent") != 1:
        raise RuntimeError("open antecedent accounting drifted")
    wedge_formula = rows[12].get("formula", "")
    if "tau(x)=min(1/4,25/L(x))" not in wedge_formula:
        raise RuntimeError("direct wedge reduction missing")
    if "must not be defined as" not in rows[13].get("formula", ""):
        raise RuntimeError("relative-majorant nonpromotion guard missing")

    summary = artifact.get("summary", {})
    if summary != {
        "rows": 14,
        "exact_reductions": 11,
        "asymptotic_compositions": 1,
        "open_antecedents": 2,
        "new_strip_open_antecedents": 0,
        "old_collar_open_antecedents": 1,
    }:
        raise RuntimeError(f"summary drifted: {summary}")

    samples = artifact.get("schedule_samples", [])
    if [row["stage"] for row in samples] != [0, 1, 10, 100, 1000, 10_000]:
        raise RuntimeError("schedule samples drifted")
    if any(row["new_strip_min_tL"] != "25" for row in samples):
        raise RuntimeError("stored schedule leaves the explicit ray")

    check_symbolic_identities()
    check_schedule()

    note = NOTE.read_text(encoding="utf-8")
    required = (
        "# Newman Ray-Aligned Parabolic-Frequency Reduction",
        "Every New Strip Is Already Closed",
        "Parabolic-Frequency Scale",
        "Relative Two-Chart Handoff",
        "Conditional Cofinal Theorem",
        "equivalent direct wedge formulation",
        "tau(x)=min(1/4,25/L(x))",
        "essential nonpromotion guard",
        "without dividing by the unknown jet norm",
        "The result does not prove the relative Xi bound",
    )
    missing = [text for text in required if text not in note]
    if missing:
        raise RuntimeError(f"note contract missing: {missing}")
    if "this is not a proof of `Lambda<=0` or RH" not in note:
        raise RuntimeError("note nonpromotion guard missing")

    print(
        "validated Newman ray-aligned parabolic-frequency reduction: "
        "14 rows, 11 exact reductions, 1 asymptotic composition, "
        "0 open new-strip antecedents, 1 open Xi descendant theorem"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
