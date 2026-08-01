#!/usr/bin/env python3
"""Validate the stored two-block Newman shell frontier scout."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_two_block_shell_frontier_scout as target  # noqa: E402


def arb(text: str):
    return target.two_block.compact.arb(text)


def parse_bounds(row: dict) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    return (
        Fraction(row["t_low"]),
        Fraction(row["t_high"]),
        Fraction(row["x_low"]),
        Fraction(row["x_high"]),
    )


def separation(row: dict) -> tuple[bool, bool, bool]:
    j_lower = arb(row["j_two_block_lower"])
    j_upper = arb(row["j_two_block_upper"])
    jp_lower = arb(row["j_two_block_prime_lower"])
    jp_upper = arb(row["j_two_block_prime_upper"])
    value_tail = arb(row["tail_value_upper"]).upper()
    derivative_tail = arb(row["tail_derivative_upper"]).upper()
    if not j_lower.upper() <= j_upper.lower():
        raise ValueError("reversed value endpoints")
    if not jp_lower.upper() <= jp_upper.lower():
        raise ValueError("reversed derivative endpoints")
    value = (
        j_lower.lower() > value_tail
        or j_upper.upper() < -value_tail
    )
    derivative = (
        jp_lower.lower() > derivative_tail
        or jp_upper.upper() < -derivative_tail
    )
    negative = j_upper.upper() < -value_tail
    return value, derivative, negative


def validate_shell(row: dict, header: dict) -> list[str]:
    issues: list[str] = []
    j = row.get("j")
    result_path = REPO_ROOT / row.get("result_path", "")
    try:
        resolved = result_path.resolve()
        root = (
            REPO_ROOT / "work/rh_compute/results" / target.STEM
        ).resolve()
        if root not in resolved.parents:
            issues.append(f"j={j}: result path escapes frontier directory")
            return issues
    except OSError as exc:
        issues.append(f"j={j}: invalid result path: {exc}")
        return issues
    if not resolved.exists():
        issues.append(f"j={j}: result file is missing")
        return issues
    if target.sha256_path(resolved) != row.get("result_sha256"):
        issues.append(f"j={j}: result hash mismatch")
        return issues
    try:
        result = json.loads(resolved.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"j={j}: invalid result JSON: {exc}")
        return issues

    if result.get("kind") != f"{target.STEM}_shell":
        issues.append(f"j={j}: result kind mismatch")
    if result.get("j") != j:
        issues.append(f"j={j}: result shell index mismatch")
    if result.get("config") != header["config"]:
        issues.append(f"j={j}: result config mismatch")
    if result.get("source_sha256") != header["source_sha256"]:
        issues.append(f"j={j}: result source hashes mismatch")
    if result.get("summary") != row.get("summary"):
        issues.append(f"j={j}: row/result summary mismatch")

    expected_domain = target.shell_domain(j)
    domain = result.get("domain", {})
    stored_domain = tuple(
        Fraction(domain[key])
        for key in (
            "time_lower", "time_upper", "x_lower", "x_upper"
        )
    )
    if stored_domain != expected_domain:
        issues.append(f"j={j}: shell domain mismatch")

    config = header["config"]
    time_step = Fraction(config["time_step"])
    x_step = Fraction(config["x_step"])
    queue = [
        (bounds, 0)
        for bounds in target.initial_boxes(j, time_step, x_step)
    ]
    cursor = 0
    evaluations = result.get("evaluations", [])
    subdivisions = 0
    certified: list[dict] = []
    uncertified: list[dict] = []
    value_count = 0
    derivative_count = 0
    negative_count = 0

    for index, evaluation in enumerate(evaluations):
        if cursor >= len(queue):
            issues.append(f"j={j}: evaluation {index} exceeds queue")
            break
        expected_bounds, expected_depth = queue[cursor]
        cursor += 1
        try:
            bounds = parse_bounds(evaluation)
            depth = int(evaluation["depth"])
        except (KeyError, ValueError, ZeroDivisionError) as exc:
            issues.append(f"j={j}: invalid evaluation {index}: {exc}")
            continue
        if bounds != expected_bounds or depth != expected_depth:
            issues.append(
                f"j={j}: evaluation {index} breaks task order"
            )
        try:
            value, derivative, negative = separation(evaluation)
        except (KeyError, ValueError, ZeroDivisionError) as exc:
            issues.append(
                f"j={j}: invalid endpoints at evaluation {index}: {exc}"
            )
            continue
        stored_certified = evaluation.get("certified") is True
        if stored_certified:
            certified.append(evaluation)
            value_count += int(value)
            derivative_count += int(derivative)
            negative_count += int(negative)
            if not (value or derivative):
                issues.append(
                    f"j={j}: certified evaluation {index} is not separated"
                )
            branch = evaluation.get("branch")
            if branch == "value" and not value:
                issues.append(
                    f"j={j}: invalid value branch at evaluation {index}"
                )
            elif branch == "derivative" and not derivative:
                issues.append(
                    f"j={j}: invalid derivative branch at evaluation {index}"
                )
            elif branch not in ("value", "derivative"):
                issues.append(
                    f"j={j}: invalid branch at evaluation {index}"
                )
            if (
                evaluation.get("two_block_value_negative") is True
            ) != negative:
                issues.append(
                    f"j={j}: invalid sign flag at evaluation {index}"
                )
            if not arb(
                evaluation["certified_ratio_lower"]
            ).lower() > 1:
                issues.append(
                    f"j={j}: ratio lost at evaluation {index}"
                )
        else:
            uncertified.append(evaluation)
            if value or derivative:
                issues.append(
                    f"j={j}: rejected evaluation {index} is separated"
                )
            if depth < config["max_depth"]:
                queue.extend(
                    (child, depth + 1)
                    for child in target.two_block.split_box(
                        bounds, time_step, x_step
                    )
                )
                subdivisions += 1

    summary = result.get("summary", {})
    expected_counts = {
        "initial_boxes": len(
            target.initial_boxes(j, time_step, x_step)
        ),
        "evaluated_boxes": len(evaluations),
        "certified_leaf_boxes": len(certified),
        "subdivisions": subdivisions,
        "uncertified_evaluations": len(uncertified),
        "terminal_unresolved_boxes": sum(
            row.get("certified") is False
            and row.get("depth") == config["max_depth"]
            for row in evaluations
        ),
        "queue_remaining": len(queue) - cursor,
        "negative_value_boxes": negative_count,
        "value_separated_boxes": value_count,
        "derivative_separated_boxes": derivative_count,
        "derivative_only_boxes": len(certified) - value_count,
    }
    for key, expected in expected_counts.items():
        if summary.get(key) != expected:
            issues.append(
                f"j={j}: summary {key} mismatch "
                f"{summary.get(key)!r}!={expected!r}"
            )

    if certified:
        minimum = min(
            certified,
            key=lambda item: float(
                arb(item["certified_ratio_lower"]).lower()
            ),
        )
        if (
            summary.get("minimum_certified_ratio_lower")
            != minimum["certified_ratio_lower"]
        ):
            issues.append(f"j={j}: minimum ratio mismatch")
    status = summary.get("status")
    if status != row.get("status"):
        issues.append(f"j={j}: row status mismatch")
    if status == "certified":
        if cursor != len(queue):
            issues.append(f"j={j}: certified shell leaves queued boxes")
        if any(
            item.get("certified") is False
            and item.get("depth") == config["max_depth"]
            for item in evaluations
        ):
            issues.append(
                f"j={j}: certified shell has terminal unresolved boxes"
            )
        if summary.get("complete_cover") is not True:
            issues.append(f"j={j}: certified shell lost cover flag")
    elif status == "unresolved":
        if not evaluations:
            issues.append(f"j={j}: unresolved shell has no evaluations")
        else:
            last = evaluations[-1]
            if (
                last.get("certified") is not False
                or last.get("depth") != config["max_depth"]
            ):
                issues.append(
                    f"j={j}: unresolved status lacks depth frontier"
                )
    elif status == "evaluation_limit":
        if len(evaluations) != config["max_evaluations_per_shell"]:
            issues.append(f"j={j}: evaluation limit mismatch")
    elif status == "time_limit":
        if result.get("elapsed_seconds", 0) < config["max_shell_seconds"]:
            issues.append(f"j={j}: time limit mismatch")
    else:
        issues.append(f"j={j}: invalid status {status!r}")

    boundary = result.get("proof_boundary", "")
    for marker in ("Finite diagnostic", "cofinal", "Lambda<=0", "RH"):
        if marker not in boundary:
            issues.append(
                f"j={j}: proof-boundary marker missing: {marker}"
            )
    return issues


def validate(rows_path: Path, note_path: Path) -> list[str]:
    issues: list[str] = []
    try:
        rows = target.load_rows(rows_path)
    except (OSError, RuntimeError) as exc:
        return [f"invalid frontier JSONL: {exc}"]
    if not rows:
        return ["frontier JSONL is empty"]
    header = rows[0]
    if header.get("kind") != f"{target.STEM}_header":
        issues.append("header kind mismatch")
    if header.get("source_sha256") != target.source_hashes():
        issues.append("header source hashes drifted")
    config = header.get("config", {})
    if config.get("schema_version") != target.SCHEMA_VERSION:
        issues.append("schema version mismatch")
    if config.get("endpoint_serialization_digits") != 50:
        issues.append("endpoint serialization contract drifted")
    shell_rows = rows[1:]
    for offset, row in enumerate(shell_rows):
        expected_j = config.get("start_j", 5) + offset
        if row.get("j") != expected_j:
            issues.append("shell rows are not contiguous")
        issues.extend(validate_shell(row, header))
    failed = [
        row for row in shell_rows if row.get("status") != "certified"
    ]
    if len(failed) > 1:
        issues.append("more than one finite frontier row is stored")
    if failed and shell_rows[-1] is not failed[0]:
        issues.append("rows continue beyond the first finite frontier")

    try:
        note = note_path.read_text(encoding="utf-8")
    except OSError as exc:
        issues.append(f"missing frontier note: {exc}")
        note = ""
    for marker in (
        "finite interval diagnostic",
        "Passing finitely many",
        "cofinal theorem",
        "not a proof of RH",
        str(rows_path.relative_to(REPO_ROOT)).replace("\\", "/"),
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    expected_success = target.success_line(rows)
    if expected_success not in note:
        issues.append("success line missing from note")
    return issues


def main() -> int:
    issues = validate(target.DEFAULT_ROWS, target.DEFAULT_NOTE)
    if issues:
        print(
            f"Newman two-block shell frontier scout: {len(issues)} issues"
        )
        for issue in issues:
            print(f"- {issue}")
        return 1
    rows = target.load_rows(target.DEFAULT_ROWS)
    print(target.success_line(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
