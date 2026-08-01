#!/usr/bin/env python3
"""Scout the finite frontier of the two-block Newman shell method."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from time import perf_counter

import psutil


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_second_diagonal_shell_two_block_interval_certificate as two_block  # noqa: E402


STEM = "jensen_window_pf_newman_theta_two_block_shell_frontier_scout"
DEFAULT_RESULT_DIR = (
    REPO_ROOT / "work/rh_compute/results" / STEM
)
DEFAULT_ROWS = DEFAULT_RESULT_DIR / "frontier.jsonl"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SCHEMA_VERSION = 1


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_bytes(value: dict) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")


def row_sha256(row: dict) -> str:
    payload = {
        key: value for key, value in row.items()
        if key != "row_sha256"
    }
    return hashlib.sha256(canonical_bytes(payload)).hexdigest()


def source_hashes() -> dict[str, str]:
    paths = (
        Path(two_block.__file__).resolve(),
        Path(two_block.compact.__file__).resolve(),
        two_block.DEFAULT_OUT.resolve(),
    )
    return {
        str(path.relative_to(REPO_ROOT)).replace("\\", "/"):
        sha256_path(path)
        for path in paths
    }


def build_config(args: argparse.Namespace) -> dict:
    return {
        "schema_version": SCHEMA_VERSION,
        "start_j": args.start_j,
        "time_step": str(args.time_step),
        "x_step": str(args.x_step),
        "max_depth": args.max_depth,
        "max_evaluations_per_shell": args.max_evaluations_per_shell,
        "max_shell_seconds": args.max_shell_seconds,
        "cpu_sample_every_evaluations": (
            args.cpu_sample_every_evaluations
        ),
        "cpu_sample_interval_seconds": args.cpu_sample_interval_seconds,
        "cpu_park_threshold_percent": args.cpu_park_threshold_percent,
        "cpu_consecutive_park_samples": (
            args.cpu_consecutive_park_samples
        ),
        "precision_bits": two_block.compact.PRECISION_BITS,
        "endpoint_serialization_digits": two_block.ENDPOINT_DIGITS,
        "retained_theta_blocks": [1, 2],
        "tail_moment_zero_upper": str(
            two_block.TAIL_MOMENT_ZERO
        ),
        "tail_moment_one_upper": str(
            two_block.TAIL_MOMENT_ONE
        ),
    }


def header_row(config: dict) -> dict:
    row = {
        "kind": f"{STEM}_header",
        "date": "2026-07-24",
        "proof_boundary": (
            "Finite shell diagnostics only. Passing any finite number of "
            "shells does not prove a cofinal theorem, Lambda<=0, or RH."
        ),
        "config": config,
        "source_sha256": source_hashes(),
        "previous_row_sha256": None,
    }
    row["row_sha256"] = row_sha256(row)
    return row


def load_rows(path: Path) -> list[dict]:
    if not path.exists():
        return []
    rows: list[dict] = []
    previous: str | None = None
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as exc:
                raise RuntimeError(
                    f"invalid JSONL row {line_number}: {exc}"
                ) from exc
            if row.get("previous_row_sha256") != previous:
                raise RuntimeError(
                    f"row {line_number} breaks the hash chain"
                )
            if row.get("row_sha256") != row_sha256(row):
                raise RuntimeError(
                    f"row {line_number} has an invalid row hash"
                )
            previous = row["row_sha256"]
            rows.append(row)
    return rows


def append_row(path: Path, row: dict, previous: str | None) -> dict:
    row = {**row, "previous_row_sha256": previous}
    row["row_sha256"] = row_sha256(row)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    return row


def write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def shell_domain(j: int) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    return (
        Fraction(1, 5 * j),
        Fraction(1, 5),
        Fraction(38),
        Fraction(38 + j),
    )


def initial_boxes(
    j: int, time_step: Fraction, x_step: Fraction
) -> list[tuple[Fraction, Fraction, Fraction, Fraction]]:
    time_lower, time_upper, x_lower, x_upper = shell_domain(j)
    return [
        (time_low, time_high, x_low, x_high)
        for time_low, time_high in two_block.compact.fraction_range(
            time_lower, time_upper, time_step
        )
        for x_low, x_high in two_block.compact.fraction_range(
            x_lower, x_upper, x_step
        )
    ]


def lower_float(row: dict, key: str) -> float:
    return float(two_block.compact.arb(row[key]).lower())


def minimum_record(records: list[dict], key: str) -> dict | None:
    if not records:
        return None
    return min(records, key=lambda row: lower_float(row, key))


def record_location(row: dict | None) -> dict | None:
    if row is None:
        return None
    return {
        key: row[key]
        for key in ("depth", "t_low", "t_high", "x_low", "x_high")
    }


def summarize(
    status: str,
    initial_count: int,
    evaluations: list[dict],
    subdivisions: int,
    queue_remaining: int,
) -> dict:
    certified = [
        row for row in evaluations if row.get("certified") is True
    ]
    uncertified = [
        row for row in evaluations
        if row.get("certified") is False
    ]
    value_separated = [
        row for row in certified
        if lower_float(row, "value_ratio_lower") > 1
    ]
    derivative_separated = [
        row for row in certified
        if lower_float(row, "derivative_ratio_lower") > 1
    ]
    minimum = minimum_record(certified, "certified_ratio_lower")
    minimum_value = minimum_record(certified, "value_ratio_lower")
    minimum_derivative = minimum_record(
        certified, "derivative_ratio_lower"
    )
    return {
        "status": status,
        "initial_boxes": initial_count,
        "evaluated_boxes": len(evaluations),
        "certified_leaf_boxes": len(certified),
        "subdivisions": subdivisions,
        "uncertified_evaluations": len(uncertified),
        "terminal_unresolved_boxes": (
            1 if status == "unresolved" else 0
        ),
        "queue_remaining": queue_remaining,
        "maximum_certified_depth": max(
            (row["depth"] for row in certified), default=-1
        ),
        "branch_counts": {
            "value": sum(
                row.get("branch") == "value" for row in certified
            ),
            "derivative": sum(
                row.get("branch") == "derivative"
                for row in certified
            ),
        },
        "negative_value_boxes": sum(
            row.get("two_block_value_negative") is True
            for row in certified
        ),
        "value_separated_boxes": len(value_separated),
        "derivative_separated_boxes": len(derivative_separated),
        "derivative_only_boxes": (
            len(certified) - len(value_separated)
        ),
        "minimum_certified_ratio_lower": (
            minimum["certified_ratio_lower"]
            if minimum is not None else None
        ),
        "minimum_certified_location": record_location(minimum),
        "minimum_value_ratio_lower": (
            minimum_value["value_ratio_lower"]
            if minimum_value is not None else None
        ),
        "minimum_value_location": record_location(minimum_value),
        "minimum_derivative_ratio_lower": (
            minimum_derivative["derivative_ratio_lower"]
            if minimum_derivative is not None else None
        ),
        "minimum_derivative_location": record_location(
            minimum_derivative
        ),
        "complete_cover": (
            status == "certified" and queue_remaining == 0
        ),
    }


def run_shell(
    j: int,
    config: dict,
    priority_lowered: bool,
) -> dict:
    time_step = Fraction(config["time_step"])
    x_step = Fraction(config["x_step"])
    boxes = initial_boxes(j, time_step, x_step)
    queue = [(bounds, 0) for bounds in boxes]
    certifier = two_block.TwoBlockCertifier()
    evaluations: list[dict] = []
    cpu_samples: list[float] = []
    high_cpu_samples = 0
    cpu_park_requested = False
    subdivisions = 0
    cursor = 0
    status = "running"
    started = perf_counter()

    while cursor < len(queue):
        elapsed = perf_counter() - started
        if len(evaluations) >= config["max_evaluations_per_shell"]:
            status = "evaluation_limit"
            break
        if elapsed >= config["max_shell_seconds"]:
            status = "time_limit"
            break

        bounds, depth = queue[cursor]
        cursor += 1
        record = certifier.certify_two_block_box(*bounds, depth)
        evaluations.append(record)
        if not record["certified"]:
            if depth < config["max_depth"]:
                queue.extend(
                    (child, depth + 1)
                    for child in two_block.split_box(
                        bounds, time_step, x_step
                    )
                )
                subdivisions += 1
            else:
                status = "unresolved"
                break

        sample_every = config["cpu_sample_every_evaluations"]
        if len(evaluations) % sample_every == 0:
            sample = float(
                psutil.cpu_percent(
                    interval=config["cpu_sample_interval_seconds"]
                )
            )
            cpu_samples.append(sample)
            if sample > config["cpu_park_threshold_percent"]:
                high_cpu_samples += 1
            else:
                high_cpu_samples = 0
            if (
                high_cpu_samples
                >= config["cpu_consecutive_park_samples"]
            ):
                cpu_park_requested = True

    if status == "running":
        status = "certified"
    elapsed = perf_counter() - started
    summary = summarize(
        status,
        len(boxes),
        evaluations,
        subdivisions,
        len(queue) - cursor,
    )
    time_lower, time_upper, x_lower, x_upper = shell_domain(j)
    return {
        "kind": f"{STEM}_shell",
        "date": "2026-07-24",
        "j": j,
        "domain": {
            "time_lower": str(time_lower),
            "time_upper": str(time_upper),
            "x_lower": str(x_lower),
            "x_upper": str(x_upper),
        },
        "config": config,
        "source_sha256": source_hashes(),
        "below_normal_priority_applied": priority_lowered,
        "elapsed_seconds": round(elapsed, 6),
        "cpu_samples_percent": cpu_samples,
        "cpu_park_requested": cpu_park_requested,
        "summary": summary,
        "evaluations": evaluations,
        "proof_boundary": (
            f"Finite diagnostic for Q_{j} only. Even a certified shell "
            "does not prove the cofinal family, Lambda<=0, or RH."
        ),
    }


def shell_row(
    result: dict, result_path: Path, result_sha256: str
) -> dict:
    summary = result["summary"]
    return {
        "kind": f"{STEM}_row",
        "date": "2026-07-24",
        "j": result["j"],
        "status": summary["status"],
        "result_path": str(
            result_path.relative_to(REPO_ROOT)
        ).replace("\\", "/"),
        "result_sha256": result_sha256,
        "elapsed_seconds": result["elapsed_seconds"],
        "cpu_samples_percent": result["cpu_samples_percent"],
        "cpu_park_requested": result["cpu_park_requested"],
        "summary": summary,
        "proof_boundary": result["proof_boundary"],
    }


def success_line(rows: list[dict]) -> str:
    shell_rows = rows[1:]
    certified = sum(row["status"] == "certified" for row in shell_rows)
    if shell_rows:
        span = f"j={shell_rows[0]['j']}..{shell_rows[-1]['j']}"
        frontier = (
            "none"
            if shell_rows[-1]["status"] == "certified"
            else f"j={shell_rows[-1]['j']}:{shell_rows[-1]['status']}"
        )
    else:
        span = "empty"
        frontier = "none"
    return (
        "validated Newman theta two-block shell frontier scout: "
        f"{len(shell_rows)} shell rows, 0 issues, {span}, "
        f"{certified} certified shells, first finite frontier={frontier}"
    )


def render_note(rows: list[dict], rows_path: Path) -> str:
    header = rows[0]
    shell_rows = rows[1:]
    lines = [
        "# Newman Theta Two-Block Shell Frontier Scout",
        "",
        "Date: 2026-07-24",
        "",
        "Status: finite interval diagnostic, not a proof of RH or "
        "`Lambda <= 0`.",
        "",
        "This scout measures where the fixed two-block first-jet method "
        "continues to certify finite diagonal shells. Passing finitely many "
        "rows cannot establish a cofinal theorem.",
        "",
        "```text",
        str(rows_path.relative_to(REPO_ROOT)).replace("\\", "/"),
        "python work/rh_compute/scripts/"
        "jensen_window_pf_newman_theta_two_block_shell_frontier_scout.py",
        "python work/rh_compute/scripts/"
        "check_jensen_window_pf_newman_theta_two_block_shell_frontier_scout.py",
        "```",
        "",
        "## Fixed Configuration",
        "",
        "```text",
        json.dumps(header["config"], sort_keys=True),
        "```",
        "",
        "## Shell Rows",
        "",
        "| j | status | initial | evaluated | subdiv | negative | "
        "value-separated | derivative-only | minimum ratio | seconds |",
        "|---:|:---|---:|---:|---:|---:|---:|---:|:---|---:|",
    ]
    for row in shell_rows:
        summary = row["summary"]
        lines.append(
            f"| {row['j']} | {row['status']} | "
            f"{summary['initial_boxes']} | "
            f"{summary['evaluated_boxes']} | "
            f"{summary['subdivisions']} | "
            f"{summary['negative_value_boxes']} | "
            f"{summary['value_separated_boxes']} | "
            f"{summary['derivative_only_boxes']} | "
            f"{summary['minimum_certified_ratio_lower']} | "
            f"{row['elapsed_seconds']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Proof Boundary",
            "",
            "The rows are rigorous finite interval diagnostics for the "
            "stored boxes, but this artifact is not promoted as a theorem "
            "for any new `Q_j`. A cofinal result needs an analytic "
            "`j`-dependent separation theorem or another uniform mechanism.",
            "",
            success_line(rows),
            "",
        ]
    )
    return "\n".join(lines)


def write_note(path: Path, rows: list[dict], rows_path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(render_note(rows, rows_path))
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def sample_baseline(seconds: int) -> list[float]:
    return [
        float(psutil.cpu_percent(interval=1.0))
        for _ in range(seconds)
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-j", type=int, default=5)
    parser.add_argument("--max-j", type=int, default=20)
    parser.add_argument(
        "--time-step", type=Fraction, default=Fraction(1, 100)
    )
    parser.add_argument(
        "--x-step", type=Fraction, default=Fraction(1, 2)
    )
    parser.add_argument("--max-depth", type=int, default=6)
    parser.add_argument(
        "--max-evaluations-per-shell", type=int, default=20000
    )
    parser.add_argument(
        "--max-shell-seconds", type=float, default=180.0
    )
    parser.add_argument(
        "--max-run-seconds", type=float, default=9000.0
    )
    parser.add_argument(
        "--baseline-sample-seconds", type=int, default=5
    )
    parser.add_argument(
        "--baseline-max-percent", type=float, default=62.5
    )
    parser.add_argument(
        "--cpu-sample-every-evaluations", type=int, default=25
    )
    parser.add_argument(
        "--cpu-sample-interval-seconds", type=float, default=0.2
    )
    parser.add_argument(
        "--cpu-park-threshold-percent", type=float, default=85.0
    )
    parser.add_argument(
        "--cpu-consecutive-park-samples", type=int, default=2
    )
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.start_j < 2 or args.max_j < args.start_j:
        raise SystemExit("invalid shell range")
    if args.cpu_sample_every_evaluations < 1:
        raise SystemExit("CPU sample cadence must be positive")

    priority_lowered = (
        two_block.compact.request_below_normal_priority()
    )
    config = build_config(args)
    rows = load_rows(args.rows)
    if not rows:
        header = append_row(args.rows, header_row(config), None)
        rows = [header]
    else:
        expected = header_row(config)
        if rows[0].get("config") != expected["config"]:
            raise RuntimeError(
                "existing frontier config differs from this invocation"
            )
        if rows[0].get("source_sha256") != expected["source_sha256"]:
            raise RuntimeError(
                "frontier source hashes changed; start a new cache"
            )

    shell_rows = rows[1:]
    for offset, row in enumerate(shell_rows):
        expected_j = args.start_j + offset
        if row.get("j") != expected_j:
            raise RuntimeError("frontier shell rows are not contiguous")
    if shell_rows and shell_rows[-1]["status"] != "certified":
        write_note(args.note, rows, args.rows)
        print(success_line(rows))
        return 0

    next_j = (
        args.start_j if not shell_rows else shell_rows[-1]["j"] + 1
    )
    if next_j > args.max_j:
        write_note(args.note, rows, args.rows)
        print(success_line(rows))
        return 0

    baseline = sample_baseline(args.baseline_sample_seconds)
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "frontier baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, "
        f"mean={baseline_mean:.1f}"
    )
    if baseline_mean > args.baseline_max_percent:
        print(
            "frontier compute deferred: baseline exceeds "
            f"{args.baseline_max_percent:.1f}%"
        )
        return 2

    run_started = perf_counter()
    for j in range(next_j, args.max_j + 1):
        if perf_counter() - run_started >= args.max_run_seconds:
            print("frontier parked before a new shell: run-time limit")
            break
        result = run_shell(j, config, priority_lowered)
        result_path = args.rows.parent / f"shell_{j:04d}.json"
        if result_path.exists():
            raise RuntimeError(
                f"unindexed shell result already exists: {result_path}"
            )
        write_json_atomic(result_path, result)
        row = shell_row(
            result, result_path, sha256_path(result_path)
        )
        row = append_row(
            args.rows, row, rows[-1]["row_sha256"]
        )
        rows.append(row)
        write_note(args.note, rows, args.rows)
        summary = result["summary"]
        print(
            f"frontier j={j}: status={summary['status']}, "
            f"evaluated={summary['evaluated_boxes']}, "
            f"subdivisions={summary['subdivisions']}, "
            "minimum_ratio="
            f"{summary['minimum_certified_ratio_lower']}, "
            f"seconds={result['elapsed_seconds']:.3f}"
        )
        if summary["status"] != "certified":
            break
        if result["cpu_park_requested"]:
            print(
                "frontier parked after certified shell: sustained "
                "total-CPU threshold"
            )
            break

    print(success_line(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
