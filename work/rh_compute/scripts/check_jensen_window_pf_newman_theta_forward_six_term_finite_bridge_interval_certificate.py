#!/usr/bin/env python3
"""Validate the resumable six-term finite-bridge interval certificate."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_compact_transversality_interval_certificate as compact  # noqa: E402
import check_jensen_window_pf_newman_theta_forward_six_term_bridge_tail_gate as tail_check  # noqa: E402


STEM = (
    "jensen_window_pf_newman_theta_forward_six_term_"
    "finite_bridge_interval_certificate"
)
RESULT_DIR = REPO_ROOT / "work" / "rh_compute" / "results"
DEFAULT_ARTIFACT = RESULT_DIR / f"{STEM}.json"
DEFAULT_CACHE = RESULT_DIR / f"{STEM}.jsonl"
BUILDER = SCRIPT_DIR / f"{STEM}.py"
TAIL_STEM = "jensen_window_pf_newman_theta_forward_six_term_bridge_tail_gate"
TAIL_SOURCE = RESULT_DIR / f"{TAIL_STEM}.json"
TAIL_BUILDER = SCRIPT_DIR / f"{TAIL_STEM}.py"
TAIL_CHECKER = SCRIPT_DIR / f"check_{TAIL_STEM}.py"
MODULAR_SOURCE = (
    RESULT_DIR / "jensen_window_pf_newman_theta_modular_blend_gate.json"
)
OUTER_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_theta_stable_remainder_outer_tail_gate.json"
)
COMPACT_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_theta_compact_transversality_interval_certificate.json"
)
Q31_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_theta_modular_retained_q31_interval_certificate.json"
)
DIAGONAL_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.json"
)
WINDING_SOURCE = (
    RESULT_DIR / "jensen_window_pf_newman_first_jet_winding_gate.json"
)
BOUNDARY_SOURCE = (
    RESULT_DIR
    / "jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
)

SCHEMA = "newman_theta_forward_six_term_finite_bridge_interval_v1"
PRECISION_BITS = 352
TIME_GLOBAL_LOWER = Fraction(1, 1035)
TIME_GLOBAL_UPPER = Fraction(1, 5)
LOW_STRIP_TIME_UPPER = Fraction(1, 155)
LOW_STRIP_X_LOWER = Fraction(38)
LOW_STRIP_X_UPPER = Fraction(69)
RIGHT_STRIP_X_LOWER = Fraction(69)
RIGHT_STRIP_X_UPPER = Fraction(245)
INITIAL_TIME_STEP = Fraction(1, 100)
INITIAL_X_STEP = Fraction(1, 2)
TAYLOR_TIME_ORDER = 30
TAYLOR_X_ORDER = 24
MAX_SUBDIVISION_DEPTH = 8
MAX_OSCILLATORY_MOMENT = 85
MAX_ABSOLUTE_MOMENT = 86
EXPECTED_IDS = [
    "ntfsfbic_01_six_term_forward_partition",
    "ntfsfbic_02_direct_six_term_tail",
    "ntfsfbic_03_shared_taylor_panel",
    "ntfsfbic_04_finite_bridge_cover",
    "ntfsfbic_05_q207_no_contact",
    "ntfsfbic_06_q207_winding",
    "ntfsfbic_07_cofinal_handoff",
]


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def fraction_range(
    lower: Fraction,
    upper: Fraction,
    step: Fraction,
) -> list[tuple[Fraction, Fraction]]:
    rows: list[tuple[Fraction, Fraction]] = []
    current = lower
    while current < upper:
        following = min(current + step, upper)
        rows.append((current, following))
        current = following
    return rows


def expected_tasks() -> list[dict]:
    tasks: list[dict] = []
    for lower, upper in fraction_range(
        LOW_STRIP_X_LOWER,
        LOW_STRIP_X_UPPER,
        INITIAL_X_STEP,
    ):
        tasks.append(
            {
                "region": "low_time_left_strip",
                "x_low": str(lower),
                "x_high": str(upper),
            }
        )
    for lower, upper in fraction_range(
        RIGHT_STRIP_X_LOWER,
        RIGHT_STRIP_X_UPPER,
        INITIAL_X_STEP,
    ):
        tasks.append(
            {
                "region": "full_time_right_strip",
                "x_low": str(lower),
                "x_high": str(upper),
            }
        )
    return tasks


def task_id(task: dict) -> str:
    return f"{task['region']}:x{task['x_low']}..{task['x_high']}"


def domain_for_task(
    task: dict,
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    time_upper = (
        LOW_STRIP_TIME_UPPER
        if task["region"] == "low_time_left_strip"
        else TIME_GLOBAL_UPPER
    )
    return (
        TIME_GLOBAL_LOWER,
        time_upper,
        Fraction(task["x_low"]),
        Fraction(task["x_high"]),
    )


def positive_overlap(
    left: tuple[Fraction, Fraction, Fraction, Fraction],
    right: tuple[Fraction, Fraction, Fraction, Fraction],
) -> bool:
    return (
        max(left[0], right[0]) < min(left[1], right[1])
        and max(left[2], right[2]) < min(left[3], right[3])
    )


def parse_leaf_bounds(
    leaf: dict,
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    return (
        Fraction(leaf["t_low"]),
        Fraction(leaf["t_high"]),
        Fraction(leaf["x_low"]),
        Fraction(leaf["x_high"]),
    )


def interval_abs_lower(lower: compact.arb, upper: compact.arb) -> compact.arb:
    if lower.lower() <= 0 <= upper.upper():
        return compact.arb(0)
    return min(abs(lower), abs(upper))


def validate_leaf(
    issues: list[str],
    label: str,
    leaf: dict,
    panel_domain: tuple[Fraction, Fraction, Fraction, Fraction],
) -> tuple[Fraction, Fraction, Fraction, Fraction] | None:
    try:
        bounds = parse_leaf_bounds(leaf)
    except Exception as exc:
        issues.append(f"{label}: invalid rational bounds: {exc}")
        return None
    if (
        bounds[0] < panel_domain[0]
        or bounds[1] > panel_domain[1]
        or bounds[2] < panel_domain[2]
        or bounds[3] > panel_domain[3]
        or bounds[0] >= bounds[1]
        or bounds[2] >= bounds[3]
    ):
        issues.append(f"{label}: leaf lies outside its panel")
    if not leaf.get("certified"):
        issues.append(f"{label}: stored certified leaf is not certified")
    if not 0 <= int(leaf.get("depth", -1)) <= MAX_SUBDIVISION_DEPTH:
        issues.append(f"{label}: invalid subdivision depth")

    try:
        j_lower = compact.arb(leaf["j_retained_lower"])
        j_upper = compact.arb(leaf["j_retained_upper"])
        d_lower = compact.arb(leaf["j_retained_prime_lower"])
        d_upper = compact.arb(leaf["j_retained_prime_upper"])
        tail_value = compact.arb(leaf["tail_value_upper"])
        tail_derivative = compact.arb(leaf["tail_derivative_upper"])
        value_ratio = compact.arb(leaf["value_ratio_lower"])
        derivative_ratio = compact.arb(leaf["derivative_ratio_lower"])
        certified_ratio = compact.arb(
            leaf["certified_ratio_lower"]
        )
    except Exception as exc:
        issues.append(f"{label}: directed endpoint parse failed: {exc}")
        return bounds
    if (
        tail_value.lower() <= 0
        or tail_derivative.lower() <= 0
    ):
        issues.append(f"{label}: nonpositive tail bar")
    if j_lower.upper() > j_upper.lower():
        issues.append(f"{label}: reversed J interval")
    if d_lower.upper() > d_upper.lower():
        issues.append(f"{label}: reversed J' interval")

    value_lower = interval_abs_lower(j_lower, j_upper)
    derivative_lower = interval_abs_lower(d_lower, d_upper)
    rebuilt_value_ratio = value_lower / tail_value
    rebuilt_derivative_ratio = derivative_lower / tail_derivative
    if value_ratio.lower() > rebuilt_value_ratio.upper():
        issues.append(f"{label}: value ratio is overstated")
    if derivative_ratio.lower() > rebuilt_derivative_ratio.upper():
        issues.append(f"{label}: derivative ratio is overstated")

    value_pass = value_lower.lower() > tail_value.upper()
    derivative_pass = derivative_lower.lower() > tail_derivative.upper()
    branch = leaf.get("branch")
    if branch == "value" and not value_pass:
        issues.append(f"{label}: invalid value branch")
    elif branch == "derivative" and not derivative_pass:
        issues.append(f"{label}: invalid derivative branch")
    elif branch not in ("value", "derivative"):
        issues.append(f"{label}: invalid certification branch")
    if not (value_pass or derivative_pass):
        issues.append(f"{label}: neither strict disjunction is certified")
    branch_ratio = (
        rebuilt_value_ratio
        if branch == "value"
        else rebuilt_derivative_ratio
    )
    if certified_ratio.lower() > branch_ratio.upper():
        issues.append(f"{label}: certified ratio is overstated")

    expected_negative = bool(j_upper.upper() < -tail_value.lower())
    expected_positive = bool(j_lower.lower() > tail_value.upper())
    if bool(leaf.get("full_value_negative")) != expected_negative:
        issues.append(f"{label}: full-value-negative flag mismatch")
    if bool(leaf.get("full_value_positive")) != expected_positive:
        issues.append(f"{label}: full-value-positive flag mismatch")
    return bounds


def validate_panel(
    issues: list[str],
    sequence: int,
    task: dict,
    result: dict,
) -> None:
    label = f"panel {sequence}"
    domain = domain_for_task(task)
    expected_time_cells = len(
        fraction_range(
            domain[0],
            domain[1],
            INITIAL_TIME_STEP,
        )
    )
    if result.get("initial_time_cells") != expected_time_cells:
        issues.append(f"{label}: initial time-cell count mismatch")
    leaves = result.get("certified_leaves", [])
    unresolved = result.get("unresolved_leaves", [])
    if result.get("certified_leaf_boxes") != len(leaves):
        issues.append(f"{label}: certified leaf count mismatch")
    if result.get("unresolved_boxes") != len(unresolved):
        issues.append(f"{label}: unresolved leaf count mismatch")
    if result.get("status") == "certified" and unresolved:
        issues.append(f"{label}: certified status has unresolved leaves")
    if result.get("status") == "unresolved" and not unresolved:
        issues.append(f"{label}: unresolved status has no unresolved leaves")
    if result.get("evaluated_boxes", 0) < len(leaves) + len(unresolved):
        issues.append(f"{label}: evaluated-box count is too small")
    if result.get("subdivisions", 0) < 0:
        issues.append(f"{label}: negative subdivision count")

    model = result.get("model", {})
    x_low = Fraction(task["x_low"])
    x_high = Fraction(task["x_high"])
    if model.get("x_center") != str((x_low + x_high) / 2):
        issues.append(f"{label}: model center mismatch")
    if model.get("x_radius") != str((x_high - x_low) / 2):
        issues.append(f"{label}: model radius mismatch")
    if model.get("oscillatory_moments") != MAX_OSCILLATORY_MOMENT + 1:
        issues.append(f"{label}: oscillatory moment count mismatch")
    if int(model.get("minimum_integral_accuracy_bits", -1)) < 128:
        issues.append(f"{label}: weak oscillatory integral accuracy")
    for remainder_name in ("h_remainder", "h_prime_remainder"):
        remainder = model.get(remainder_name, {})
        parts: dict[str, compact.arb] = {}
        for key in ("time", "frequency", "total"):
            try:
                parts[key] = compact.arb(
                    remainder[key]["enclosure"]
                )
            except Exception as exc:
                issues.append(
                    f"{label}: {remainder_name} {key} parse failed: {exc}"
                )
                continue
            if parts[key].lower() <= 0:
                issues.append(
                    f"{label}: {remainder_name} {key} is not positive"
                )
        if set(parts) == {"time", "frequency", "total"}:
            if not parts["total"].overlaps(
                parts["time"] + parts["frequency"]
            ):
                issues.append(
                    f"{label}: {remainder_name} total mismatch"
                )

    leaf_bounds: list[
        tuple[Fraction, Fraction, Fraction, Fraction]
    ] = []
    for index, leaf in enumerate(leaves):
        bounds = validate_leaf(
            issues,
            f"{label} leaf {index}",
            leaf,
            domain,
        )
        if bounds is not None:
            leaf_bounds.append(bounds)
    for index, leaf in enumerate(unresolved):
        try:
            bounds = parse_leaf_bounds(leaf)
            if (
                bounds[0] < domain[0]
                or bounds[1] > domain[1]
                or bounds[2] < domain[2]
                or bounds[3] > domain[3]
            ):
                issues.append(
                    f"{label} unresolved {index}: outside panel"
                )
        except Exception as exc:
            issues.append(
                f"{label} unresolved {index}: invalid bounds: {exc}"
            )

    if result.get("status") == "certified":
        target_area = (
            (domain[1] - domain[0])
            * (domain[3] - domain[2])
        )
        leaf_area = sum(
            (bounds[1] - bounds[0])
            * (bounds[3] - bounds[2])
            for bounds in leaf_bounds
        )
        if leaf_area != target_area:
            issues.append(f"{label}: certified leaf area mismatch")
        for left_index, left in enumerate(leaf_bounds):
            for right in leaf_bounds[left_index + 1 :]:
                if positive_overlap(left, right):
                    issues.append(
                        f"{label}: certified leaves overlap"
                    )
                    return


def expected_source_hashes() -> dict[str, str]:
    return {
        "six_term_tail": digest(TAIL_SOURCE),
        "modular_partition": digest(MODULAR_SOURCE),
        "outer_tail": digest(OUTER_SOURCE),
        "compact_core": digest(COMPACT_SOURCE),
        "q31": digest(Q31_SOURCE),
        "diagonal_exhaustion": digest(DIAGONAL_SOURCE),
        "winding": digest(WINDING_SOURCE),
        "boundary": digest(BOUNDARY_SOURCE),
        "tail_builder_sha256": digest(TAIL_BUILDER),
        "tail_checker_sha256": digest(TAIL_CHECKER),
    }


def validate(
    artifact_path: Path,
    cache_path: Path,
) -> list[str]:
    compact.flint.ctx.prec = PRECISION_BITS
    issues: list[str] = []
    artifact = json.loads(
        artifact_path.read_text(encoding="utf-8")
    )
    if artifact.get("kind") != STEM:
        issues.append("artifact kind mismatch")
    config = artifact.get("config", {})
    if config.get("schema") != SCHEMA:
        issues.append("config schema mismatch")
    if config.get("implementation_sha256") != digest(BUILDER):
        issues.append("builder hash mismatch")
    if config.get("precision_bits") != PRECISION_BITS:
        issues.append("precision mismatch")
    if config.get("retained_terms") != 6:
        issues.append("retained-term count mismatch")
    if config.get("first_omitted") != 7:
        issues.append("first omitted index mismatch")
    if config.get("taylor_time_order") != TAYLOR_TIME_ORDER:
        issues.append("time Taylor order mismatch")
    if config.get("taylor_x_order") != TAYLOR_X_ORDER:
        issues.append("frequency Taylor order mismatch")
    if config.get("max_subdivision_depth") != MAX_SUBDIVISION_DEPTH:
        issues.append("subdivision-depth mismatch")
    if config.get("max_oscillatory_moment") != MAX_OSCILLATORY_MOMENT:
        issues.append("maximum oscillatory moment mismatch")
    if config.get("max_absolute_moment") != MAX_ABSOLUTE_MOMENT:
        issues.append("maximum absolute moment mismatch")
    tasks = expected_tasks()
    expected_ids = [task_id(task) for task in tasks]
    if config.get("task_ids") != expected_ids:
        issues.append("task sequence mismatch")
    rebuilt_config_hash = sha256(
        canonical_json(config).encode("utf-8")
    ).hexdigest()
    if artifact.get("config_sha256") != rebuilt_config_hash:
        issues.append("artifact config hash mismatch")

    source_audit = artifact.get("source_audit", {})
    for key, expected in expected_source_hashes().items():
        if key.endswith("_sha256"):
            observed = source_audit.get(key)
        else:
            observed = source_audit.get(key, {}).get("sha256")
        if observed != expected:
            issues.append(f"source hash mismatch: {key}")
    issues.extend(
        f"six-term tail: {issue}"
        for issue in tail_check.validate(TAIL_SOURCE)
    )

    records: list[dict] = []
    previous = "0" * 64
    if cache_path.exists():
        with cache_path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.endswith("\n"):
                    issues.append(
                        f"cache line {line_number} lacks terminator"
                    )
                    break
                try:
                    record = json.loads(line)
                except Exception as exc:
                    issues.append(
                        f"cache line {line_number} parse failed: {exc}"
                    )
                    break
                stored_hash = record.pop("row_sha256", None)
                actual_hash = sha256(
                    canonical_json(record).encode("utf-8")
                ).hexdigest()
                record["row_sha256"] = stored_hash
                if stored_hash != actual_hash:
                    issues.append(
                        f"cache line {line_number} hash mismatch"
                    )
                if record.get("sequence") != line_number:
                    issues.append(
                        f"cache line {line_number} sequence mismatch"
                    )
                if record.get("config_sha256") != rebuilt_config_hash:
                    issues.append(
                        f"cache line {line_number} config mismatch"
                    )
                if record.get("previous_row_sha256") != previous:
                    issues.append(
                        f"cache line {line_number} chain mismatch"
                    )
                if line_number <= len(tasks):
                    if record.get("task") != tasks[line_number - 1]:
                        issues.append(
                            f"cache line {line_number} task-order mismatch"
                        )
                else:
                    issues.append("cache has excess records")
                if record.get("schema") != SCHEMA:
                    issues.append(
                        f"cache line {line_number} schema mismatch"
                    )
                validate_panel(
                    issues,
                    line_number,
                    record.get("task", {}),
                    record.get("result", {}),
                )
                records.append(record)
                previous = stored_hash

    cache_meta = artifact.get("cache", {})
    if cache_meta.get("records") != len(records):
        issues.append("artifact/cache record count mismatch")
    if cache_meta.get("expected_records") != len(tasks):
        issues.append("expected record count mismatch")
    if cache_meta.get("last_row_sha256") != (
        records[-1]["row_sha256"] if records else None
    ):
        issues.append("last cache hash mismatch")
    complete = len(records) == len(tasks)
    unresolved_panels = sum(
        record["result"].get("status") != "certified"
        for record in records
    )
    theorem_ready = complete and unresolved_panels == 0
    summary = artifact.get("summary", {})
    if summary.get("theorem_ready") != theorem_ready:
        issues.append("theorem-ready status mismatch")
    if summary.get("certified_panels") != (
        len(records) - unresolved_panels
    ):
        issues.append("certified panel summary mismatch")
    if summary.get("unresolved_panels") != unresolved_panels:
        issues.append("unresolved panel summary mismatch")
    if cache_meta.get("complete") != complete:
        issues.append("cache complete flag mismatch")
    if cache_meta.get("missing_count") != len(tasks) - len(records):
        issues.append("missing panel count mismatch")
    if cache_meta.get("missing_task_ids") != expected_ids[len(records) :]:
        issues.append("missing task suffix mismatch")

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

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "Q207",
        "x=245",
        "cofinal transition",
        "strict Laguerre",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary guard missing: {marker}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--artifact",
        type=Path,
        default=DEFAULT_ARTIFACT,
    )
    parser.add_argument(
        "--cache",
        type=Path,
        default=DEFAULT_CACHE,
    )
    args = parser.parse_args()
    issues = validate(args.artifact, args.cache)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        return 1
    artifact = json.loads(
        args.artifact.read_text(encoding="utf-8")
    )
    print(
        "validated Newman theta forward six-term finite-bridge "
        "interval certificate: "
        f"{len(artifact['rows'])} rows, "
        f"{artifact['cache']['records']}/"
        f"{artifact['cache']['expected_records']} panels, "
        f"{artifact['summary']['certified_leaf_boxes']} certified leaves, "
        f"{artifact['summary']['unresolved_panels']} unresolved panels, "
        "0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
