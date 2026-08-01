#!/usr/bin/env python3
"""Independently validate the resumable modular-retained Q31 certificate."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR_DIR = REPO_ROOT / "work" / "rh_compute" / "vendor"
if str(VENDOR_DIR) not in sys.path:
    sys.path.insert(0, str(VENDOR_DIR))

from flint import arb


STEM = "jensen_window_pf_newman_theta_modular_retained_q31_interval_certificate"
SCHEMA = "newman_theta_modular_retained_q31_interval_v1"
DEFAULT_ARTIFACT = (
    REPO_ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
)
BUILDER = (
    REPO_ROOT / "work" / "rh_compute" / "scripts" / f"{STEM}.py"
)
FULL_BUDGET_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_full_derivative_budget_certificate.json"
)
MODULAR_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_modular_blend_gate.json"
)
OUTER_SOURCE = (
    REPO_ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate.json"
)
PRECISION_BITS = 192
RETAINED_N = 7
DIAGONAL_J = 31
TIME_LOWER = Fraction(1, 155)
TIME_UPPER = Fraction(1, 5)
X_LOWER = Fraction(38)
X_UPPER = Fraction(69)
TIME_STEP = Fraction(1, 100)
X_STEP = Fraction(1, 2)
MAX_DEPTH = 7
EXPECTED_IDS = [
    "ntmrq31_01_even_modular_retained_sum",
    "ntmrq31_02_full_m9_tail_bars",
    "ntmrq31_03_retained_cutoff_tail",
    "ntmrq31_04_interval_slab",
    "ntmrq31_05_q31_no_contact",
    "ntmrq31_06_q31_winding",
    "ntmrq31_07_cofinal_handoff",
]


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def fraction_range(
    lower: Fraction, upper: Fraction, step: Fraction
) -> list[tuple[Fraction, Fraction]]:
    rows: list[tuple[Fraction, Fraction]] = []
    left = lower
    while left < upper:
        right = min(left + step, upper)
        rows.append((left, right))
        left = right
    return rows


def initial_tasks() -> list[dict]:
    return [
        {
            "t_low": str(t_low),
            "t_high": str(t_high),
            "x_low": str(x_low),
            "x_high": str(x_high),
        }
        for t_low, t_high in fraction_range(
            TIME_LOWER, TIME_UPPER, TIME_STEP
        )
        for x_low, x_high in fraction_range(
            X_LOWER, X_UPPER, X_STEP
        )
    ]


def task_id(task: dict) -> str:
    return (
        f"t{task['t_low']}..{task['t_high']}:"
        f"x{task['x_low']}..{task['x_high']}"
    )


def expected_config() -> dict:
    return {
        "schema": SCHEMA,
        "implementation_sha256": file_hash(BUILDER),
        "precision_bits": PRECISION_BITS,
        "retained_n": RETAINED_N,
        "diagonal_j": DIAGONAL_J,
        "domain": {
            "t": [str(TIME_LOWER), str(TIME_UPPER)],
            "x": [str(X_LOWER), str(X_UPPER)],
        },
        "initial_time_step": str(TIME_STEP),
        "initial_x_step": str(X_STEP),
        "max_depth": MAX_DEPTH,
        "integration_cutoff": "11/5",
        "integration_options": {
            "deg_limit": 80,
            "eval_limit": 250_000,
            "depth_limit": 50,
            "use_heap": True,
            "abs_tol": "2^-150",
            "rel_tol": "2^-145",
        },
        "taylor_time_order": 2,
        "taylor_x_order": 3,
        "full_budget_sha256": file_hash(FULL_BUDGET_SOURCE),
        "modular_source_sha256": file_hash(MODULAR_SOURCE),
        "outer_source_sha256": file_hash(OUTER_SOURCE),
        "task_ids": [task_id(task) for task in initial_tasks()],
    }


def read_cache(
    path: Path,
    config_hash: str,
    issues: list[str],
) -> tuple[list[dict], dict[str, dict]]:
    if not path.exists():
        return [], {}
    raw = path.read_bytes()
    if raw and not raw.endswith(b"\n"):
        issues.append("cache does not end with a line terminator")
    records: list[dict] = []
    by_id: dict[str, dict] = {}
    previous = "0" * 64
    for line_number, line in enumerate(
        raw.decode("utf-8").splitlines(), start=1
    ):
        try:
            record = json.loads(line)
        except Exception as exc:
            issues.append(f"cache JSON failed at line {line_number}: {exc}")
            continue
        stored_hash = record.get("row_sha256")
        payload = dict(record)
        payload.pop("row_sha256", None)
        actual_hash = sha256(
            canonical_json(payload).encode("utf-8")
        ).hexdigest()
        if stored_hash != actual_hash:
            issues.append(f"cache hash mismatch at line {line_number}")
        if record.get("sequence") != line_number:
            issues.append(f"cache sequence mismatch at line {line_number}")
        if record.get("schema") != SCHEMA:
            issues.append(f"cache schema mismatch at line {line_number}")
        if record.get("config_sha256") != config_hash:
            issues.append(f"cache config mismatch at line {line_number}")
        if record.get("previous_row_sha256") != previous:
            issues.append(f"cache chain mismatch at line {line_number}")
        previous = stored_hash or ""
        try:
            identifier = task_id(record["task"])
        except Exception as exc:
            issues.append(f"bad task at line {line_number}: {exc}")
            continue
        if identifier in by_id:
            issues.append(f"duplicate cache task: {identifier}")
        by_id[identifier] = record
        records.append(record)
    return records, by_id


def rectangle(record: dict) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    return (
        Fraction(record["t_low"]),
        Fraction(record["t_high"]),
        Fraction(record["x_low"]),
        Fraction(record["x_high"]),
    )


def area(bounds: tuple[Fraction, Fraction, Fraction, Fraction]) -> Fraction:
    t_low, t_high, x_low, x_high = bounds
    return (t_high - t_low) * (x_high - x_low)


def interiors_overlap(
    left: tuple[Fraction, Fraction, Fraction, Fraction],
    right: tuple[Fraction, Fraction, Fraction, Fraction],
) -> bool:
    return (
        max(left[0], right[0]) < min(left[1], right[1])
        and max(left[2], right[2]) < min(left[3], right[3])
    )


def endpoint_abs_lower(lower: arb, upper: arb) -> arb:
    if lower.lower() > 0:
        return lower
    if upper.upper() < 0:
        return -upper
    return arb(0)


def validate_evaluation(row: dict, label: str, issues: list[str]) -> None:
    try:
        j_lower = arb(row["j_retained_lower"])
        j_upper = arb(row["j_retained_upper"])
        jp_lower = arb(row["j_retained_prime_lower"])
        jp_upper = arb(row["j_retained_prime_upper"])
        tail_value = arb(row["tail_value_upper"])
        tail_derivative = arb(row["tail_derivative_upper"])
        if j_lower > j_upper:
            issues.append(f"reversed J endpoints: {label}")
        if jp_lower > jp_upper:
            issues.append(f"reversed J' endpoints: {label}")
        if tail_value.lower() <= 0 or tail_derivative.lower() <= 0:
            issues.append(f"nonpositive tail bar: {label}")
        value_lower = endpoint_abs_lower(j_lower, j_upper)
        derivative_lower = endpoint_abs_lower(jp_lower, jp_upper)
        value_pass = value_lower.lower() > tail_value.upper()
        derivative_pass = (
            derivative_lower.lower() > tail_derivative.upper()
        )
        if bool(row.get("certified")) != (value_pass or derivative_pass):
            issues.append(f"certification predicate mismatch: {label}")
        if row.get("certified"):
            branch = row.get("branch")
            if branch == "value" and not value_pass:
                issues.append(f"invalid value branch: {label}")
            if branch == "derivative" and not derivative_pass:
                issues.append(f"invalid derivative branch: {label}")
            if branch not in ("value", "derivative"):
                issues.append(f"unknown certified branch: {label}")
            if row.get("full_value_negative") and not (
                j_upper.upper() < -tail_value.upper()
            ):
                issues.append(f"invalid negative-value flag: {label}")
            if row.get("full_value_positive") and not (
                j_lower.lower() > tail_value.upper()
            ):
                issues.append(f"invalid positive-value flag: {label}")
    except Exception as exc:
        issues.append(f"evaluation parse failed for {label}: {exc}")


def validate_result(
    task: dict, result: dict, label: str, issues: list[str]
) -> None:
    if result.get("initial_bounds") != task:
        issues.append(f"initial bounds mismatch: {label}")
    evaluations = result.get("evaluations", [])
    if result.get("evaluated_boxes") != len(evaluations):
        issues.append(f"evaluation count mismatch: {label}")
    for index, evaluation in enumerate(evaluations):
        validate_evaluation(evaluation, f"{label}/eval{index}", issues)
        if int(evaluation.get("depth", -1)) > MAX_DEPTH:
            issues.append(f"depth exceeds maximum: {label}/eval{index}")

    certified = [row for row in evaluations if row.get("certified")]
    initial = (
        Fraction(task["t_low"]),
        Fraction(task["t_high"]),
        Fraction(task["x_low"]),
        Fraction(task["x_high"]),
    )
    leaves = [rectangle(row) for row in certified]
    for bounds in leaves:
        if not (
            initial[0] <= bounds[0] < bounds[1] <= initial[1]
            and initial[2] <= bounds[2] < bounds[3] <= initial[3]
        ):
            issues.append(f"certified leaf escapes initial box: {label}")
    for index, left in enumerate(leaves):
        for right in leaves[index + 1 :]:
            if interiors_overlap(left, right):
                issues.append(f"certified leaves overlap: {label}")
                break
    covered_area = sum((area(bounds) for bounds in leaves), Fraction(0))
    status = result.get("status")
    if status == "certified":
        if covered_area != area(initial):
            issues.append(f"certified leaves do not cover box: {label}")
        if result.get("unresolved_boxes") != 0:
            issues.append(f"certified result reports unresolved: {label}")
    elif status == "unresolved":
        if result.get("unresolved_boxes", 0) <= 0:
            issues.append(f"unresolved result lacks failures: {label}")
    else:
        issues.append(f"unknown initial-box status: {label}")
    if result.get("certified_leaf_boxes") != len(certified):
        issues.append(f"certified leaf count mismatch: {label}")


def validate(path: Path) -> list[str]:
    issues: list[str] = []
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("kind") != STEM:
        issues.append("artifact kind mismatch")
    config = artifact.get("config", {})
    expected = expected_config()
    if config != expected:
        issues.append("configuration mismatch")
    expected_hash = sha256(
        canonical_json(expected).encode("utf-8")
    ).hexdigest()
    if artifact.get("config_sha256") != expected_hash:
        issues.append("configuration hash mismatch")

    cache_info = artifact.get("cache", {})
    cache_path = REPO_ROOT / cache_info.get("path", "")
    try:
        cache_path.resolve().relative_to(REPO_ROOT.resolve())
    except Exception:
        issues.append("cache path escapes repository")
        return issues
    records, by_id = read_cache(cache_path, expected_hash, issues)
    if cache_info.get("records") != len(records):
        issues.append("artifact cache count mismatch")
    tasks = initial_tasks()
    expected_ids = [task_id(task) for task in tasks]
    if cache_info.get("expected_records") != len(tasks):
        issues.append("expected cache count mismatch")
    if cache_info.get("last_row_sha256") != (
        records[-1]["row_sha256"] if records else None
    ):
        issues.append("terminal cache hash mismatch")
    missing = [
        identifier for identifier in expected_ids if identifier not in by_id
    ]
    if cache_info.get("missing_count") != len(missing):
        issues.append("missing count mismatch")
    if cache_info.get("missing_task_ids") != missing:
        issues.append("missing task list mismatch")
    if cache_info.get("complete") != (not missing):
        issues.append("completion flag mismatch")

    allowed = set(expected_ids)
    for record in records:
        identifier = task_id(record["task"])
        if identifier not in allowed:
            issues.append(f"unexpected task: {identifier}")
        validate_result(
            record["task"], record["result"], identifier, issues
        )

    results = [record["result"] for record in records]
    unresolved = [
        result for result in results if result.get("status") != "certified"
    ]
    theorem_ready = not missing and not unresolved
    summary = artifact.get("summary", {})
    if summary.get("theorem_ready") != theorem_ready:
        issues.append("theorem-ready flag mismatch")
    if summary.get("certified_initial_boxes") != sum(
        result.get("status") == "certified" for result in results
    ):
        issues.append("certified initial count mismatch")
    if summary.get("unresolved_initial_boxes") != len(unresolved):
        issues.append("unresolved initial count mismatch")
    if summary.get("evaluated_boxes") != sum(
        result.get("evaluated_boxes", 0) for result in results
    ):
        issues.append("evaluated total mismatch")
    if summary.get("subdivisions") != sum(
        result.get("subdivisions", 0) for result in results
    ):
        issues.append("subdivision total mismatch")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    expected_readiness = [
        "proved",
        "proved",
        "proved",
        "proved" if theorem_ready else "not_ready_to_apply",
        "proved" if theorem_ready else "not_ready_to_apply",
        "proved" if theorem_ready else "not_ready_to_apply",
        "not_ready_to_apply",
    ]
    if [row.get("readiness") for row in rows] != expected_readiness:
        issues.append("row readiness mismatch")

    source_audit = artifact.get("source_audit", {})
    for entry in source_audit.values():
        source_path = REPO_ROOT / entry.get("path", "")
        if not source_path.exists():
            issues.append(f"audited source missing: {source_path}")
        elif entry.get("sha256") != file_hash(source_path):
            issues.append(f"audited source hash mismatch: {source_path}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "Q_31",
        "unbounded adaptive-N",
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
        "validated Newman theta modular-retained Q31 interval certificate: "
        f"{artifact['cache']['records']}/"
        f"{artifact['cache']['expected_records']} initial boxes, "
        f"{artifact['summary']['evaluated_boxes']} Taylor boxes, "
        f"{artifact['summary']['unresolved_initial_boxes']} unresolved, "
        "0 issues"
    )


if __name__ == "__main__":
    main()
