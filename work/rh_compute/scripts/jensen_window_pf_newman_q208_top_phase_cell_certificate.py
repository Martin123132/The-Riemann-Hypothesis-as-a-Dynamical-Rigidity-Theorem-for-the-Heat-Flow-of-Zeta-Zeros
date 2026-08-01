#!/usr/bin/env python3
"""Certify the complete Q208 top edge by fixed-time phase cells."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import time

import psutil

import jensen_window_pf_newman_q208_bottom_phase_cell_certificate as core


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_q208_top_phase_cell_certificate"
RESULT_DIR = REPO_ROOT / "work/rh_compute/results"
DEFAULT_CACHE = RESULT_DIR / f"{STEM}.jsonl"
DEFAULT_OUT = RESULT_DIR / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"

SCHEMA = "newman_q208_top_phase_cell_v1"
DATE = "2026-07-25"
TIME_VALUE = Fraction(1, 5)
DEFAULT_RUNTIME_LIMIT_SECONDS = 10_800.0


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def source_hashes() -> dict:
    return {
        **core.source_hashes(),
        "fixed_time_core": core.file_hash(Path(core.__file__).resolve()),
    }


def contract_payload() -> dict:
    return {
        "schema": SCHEMA,
        "time": str(TIME_VALUE),
        "x_domain": [str(core.X_LOWER), str(core.X_UPPER)],
        "initial_x_step": str(core.INITIAL_X_STEP),
        "taylor_x_order": core.TAYLOR_X_ORDER,
        "max_subdivision_depth": core.MAX_SUBDIVISION_DEPTH,
        "precision_bits": core.bridge.PRECISION_BITS,
        "retained_terms": core.bridge.RETAINED_TERMS,
        "first_omitted": core.bridge.FIRST_OMITTED,
        "integration_cutoff": str(core.bridge.INTEGRATION_CUTOFF),
        "proxy": "F_t(x)=16*(1+x^4)*H_t(x)",
        "proxy_derivative": (
            "F_t'(x)=64*x^3*H_t(x)+16*(1+x^4)*H_t'(x)"
        ),
        "source_sha256": source_hashes(),
    }


def contract_hash() -> str:
    return sha256(canonical_json(contract_payload()).encode()).hexdigest()


def cache_record(
    task: dict,
    result: dict,
    previous_hash: str,
) -> dict:
    payload = {
        "schema": SCHEMA,
        "contract_sha256": contract_hash(),
        "sequence": task["sequence"],
        "task": task,
        "previous_hash": previous_hash,
        "result": result,
    }
    payload["record_sha256"] = sha256(
        canonical_json(payload).encode()
    ).hexdigest()
    return payload


def validate_cache_record(
    record: dict,
    task: dict,
    previous_hash: str,
) -> None:
    if record.get("schema") != SCHEMA:
        raise RuntimeError("top cache schema mismatch")
    if record.get("contract_sha256") != contract_hash():
        raise RuntimeError("top cache contract hash mismatch")
    if record.get("sequence") != task["sequence"]:
        raise RuntimeError("top cache sequence mismatch")
    if record.get("task") != task:
        raise RuntimeError("top cache task mismatch")
    if record.get("previous_hash") != previous_hash:
        raise RuntimeError("top cache hash chain is broken")
    claimed = record.get("record_sha256")
    replay = dict(record)
    replay.pop("record_sha256", None)
    actual = sha256(canonical_json(replay).encode()).hexdigest()
    if claimed != actual:
        raise RuntimeError("top cache record hash mismatch")


def load_cache(path: Path, tasks: list[dict]) -> list[dict]:
    if not path.exists():
        return []
    records: list[dict] = []
    previous_hash = "GENESIS"
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) > len(tasks):
        raise RuntimeError("top cache has too many records")
    for index, line in enumerate(lines):
        if not line.strip():
            raise RuntimeError("top cache contains a blank row")
        record = json.loads(line)
        validate_cache_record(record, tasks[index], previous_hash)
        records.append(record)
        previous_hash = record["record_sha256"]
    return records


def append_cache(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(canonical_json(record) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def summarize(
    records: list[dict],
    tasks: list[dict],
    stop_reason: str,
    resource: dict,
) -> dict:
    leaves = [
        leaf
        for record in records
        for leaf in record["result"].get("certified_leaves", [])
    ]
    unresolved = sum(
        record["result"].get("unresolved_boxes", 0)
        for record in records
    )
    complete = len(records) == len(tasks) and unresolved == 0
    margins = [
        core.compact.arb(leaf["closest_point_norm_lower"])
        for leaf in leaves
    ]
    return {
        "tasks_total": len(tasks),
        "tasks_completed": len(records),
        "certified_leaf_cells": len(leaves),
        "unresolved_boxes": unresolved,
        "branch_counts": dict(
            Counter(leaf["branch"] for leaf in leaves)
        ),
        "maximum_subdivision_depth": max(
            (
                record["result"].get("maximum_depth", 0)
                for record in records
            ),
            default=0,
        ),
        "minimum_closest_point_norm_lower": (
            core.serialize_interval(min(margins)) if margins else None
        ),
        "complete_top_edge_certificate": complete,
        "stop_reason": stop_reason,
        "open_phase_chain": (
            core.build_open_chain(records) if complete else None
        ),
        "resource": resource,
    }


def render_note(artifact: dict) -> str:
    summary = artifact["summary"]
    chain = summary["open_phase_chain"]
    lines = [
        "# Newman Q208 Top Phase-Cell Certificate",
        "",
        f"Date: {DATE}",
        "",
        (
            "Status: rigorous complete fixed-time top certificate."
            if summary["complete_top_edge_certificate"]
            else "Status: resumable partial fixed-time computation."
        ),
        "This is not a proof of Q208, `Lambda<=0`, RH, or the Clay prize.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"work/rh_compute/results/{STEM}.jsonl",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Contract",
        "",
        "The six-term, 352-bit fixed-time engine certifies",
        "`F_t+iF_t'` at `t=1/5` on `0<=x<=246`, using the same",
        "axis-safe proxy and raw omitted-tail theorem as the bottom chain.",
        "",
        "## Progress",
        "",
        "```text",
        f"panels={summary['tasks_completed']}/{summary['tasks_total']}",
        f"certified cells={summary['certified_leaf_cells']}",
        f"unresolved={summary['unresolved_boxes']}",
        f"branches={summary['branch_counts']}",
        (
            "minimum closest-point norm lower="
            f"{summary['minimum_closest_point_norm_lower']}"
        ),
        f"maximum subdivision depth={summary['maximum_subdivision_depth']}",
        f"stop reason={summary['stop_reason']}",
        "```",
        "",
    ]
    if chain is not None:
        lines.extend(
            [
                "## Exact Open Phase Chain",
                "",
                "```text",
                f"cells={chain['cell_count']}",
                f"witnesses={chain['witness_count']}",
                (
                    "forward positive-ray crossing count="
                    f"{chain['open_positive_ray_crossing_count']}"
                ),
                "```",
                "",
                "The Q208 boundary traverses this top chain in reverse.",
                "Its endpoints must still be matched to the transformed",
                "right edge and positive-real symmetry axis before the",
                "closed winding is evaluated.",
                "",
            ]
        )
    lines.extend(
        [
            "## Scope",
            "",
            "The bottom chain, right/top corner witnesses, and cyclic exact",
            "winding are separate obligations. No incomplete chain is",
            "promoted to Q208 or a cofinal theorem.",
            "",
        ]
    )
    return "\n".join(lines)


def write_artifact(
    out: Path,
    note: Path,
    cache: Path,
    records: list[dict],
    tasks: list[dict],
    stop_reason: str,
    resource: dict,
) -> dict:
    summary = summarize(records, tasks, stop_reason, resource)
    complete = summary["complete_top_edge_certificate"]
    artifact = {
        "kind": STEM,
        "date": DATE,
        "status": (
            "rigorous complete Q208 top-edge phase-cell certificate"
            if complete
            else "resumable partial Q208 top-edge phase-cell computation"
        ),
        "proof_boundary": (
            "The complete status certifies only the fixed-time top path "
            "t=1/5 on 0<=x<=246 and its exact rational open phase chain. "
            "The bottom chain, transformed right-edge joins, cyclic winding, "
            "Q208, every later shell, Lambda<=0, RH, and the Clay prize "
            "remain open."
        ),
        "contract": contract_payload(),
        "contract_sha256": contract_hash(),
        "builder_sha256": core.file_hash(Path(__file__).resolve()),
        "cache": str(cache.relative_to(REPO_ROOT)).replace("\\", "/"),
        "cache_sha256": core.file_hash(cache) if cache.exists() else None,
        "cache_last_record_sha256": (
            records[-1]["record_sha256"] if records else "GENESIS"
        ),
        "summary": summary,
        "records": records,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    note.write_text(render_note(artifact), encoding="utf-8")
    return artifact


def baseline_samples() -> list[float]:
    return [
        float(psutil.cpu_percent(interval=1.0))
        for _ in range(core.BASELINE_SECONDS)
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument(
        "--max-panels",
        type=int,
        default=-1,
        help="Maximum new panels; -1 means unlimited and 0 initializes only.",
    )
    parser.add_argument(
        "--runtime-limit-seconds",
        type=float,
        default=DEFAULT_RUNTIME_LIMIT_SECONDS,
    )
    args = parser.parse_args()

    args.cache.parent.mkdir(parents=True, exist_ok=True)
    args.cache.touch(exist_ok=True)
    tasks = core.panel_tasks()
    records = load_cache(args.cache, tasks)
    priority_lowered = core.compact.request_below_normal_priority()
    baseline = baseline_samples()
    baseline_mean = sum(baseline) / len(baseline)
    resource = {
        "mode": "day",
        "worker_count": 1,
        "below_normal_priority_applied": priority_lowered,
        "baseline_samples": baseline,
        "baseline_mean": baseline_mean,
        "runtime_cpu_samples": [],
        "elapsed_seconds": 0.0,
    }
    if baseline_mean > core.BASELINE_LIMIT:
        write_artifact(
            args.out,
            args.note,
            args.cache,
            records,
            tasks,
            "deferred_baseline_cpu",
            resource,
        )
        print(
            "deferred Q208 top phase-cell certificate: "
            f"baseline mean {baseline_mean:.2f}% exceeds "
            f"{core.BASELINE_LIMIT:.2f}%"
        )
        return 0

    started = time.monotonic()
    certifier = core.FixedTimePanelCertifier(TIME_VALUE)
    previous_hash = (
        records[-1]["record_sha256"] if records else "GENESIS"
    )
    high_cpu_streak = 0
    stop_reason = "complete"
    newly_completed = 0
    for task in tasks[len(records) :]:
        if (
            args.max_panels >= 0
            and newly_completed >= args.max_panels
        ):
            stop_reason = "max_panels"
            break
        if time.monotonic() - started >= args.runtime_limit_seconds:
            stop_reason = "runtime_limit"
            break
        result = certifier.certify_panel(task)
        record = cache_record(task, result, previous_hash)
        append_cache(args.cache, record)
        records.append(record)
        previous_hash = record["record_sha256"]
        newly_completed += 1
        cpu_sample = float(psutil.cpu_percent(interval=0.25))
        resource["runtime_cpu_samples"].append(cpu_sample)
        print(
            f"Q208 top {task['sequence']}/{len(tasks)} "
            f"x={task['x_low']}..{task['x_high']} "
            f"status={result['status']} "
            f"leaves={result['certified_leaf_boxes']} "
            f"depth={result['maximum_depth']} "
            f"in {result['elapsed_seconds']:.3f}s "
            f"cpu={cpu_sample:.1f}%",
            flush=True,
        )
        if result["status"] != "certified":
            stop_reason = "unresolved_panel"
            break
        if cpu_sample > core.CPU_PARK_THRESHOLD:
            high_cpu_streak += 1
        else:
            high_cpu_streak = 0
        if high_cpu_streak >= core.CPU_CONSECUTIVE_PARK_SAMPLES:
            stop_reason = "sustained_cpu_limit"
            break
    if len(records) < len(tasks) and stop_reason == "complete":
        stop_reason = "partial_existing_cache"
    resource["elapsed_seconds"] = time.monotonic() - started
    artifact = write_artifact(
        args.out,
        args.note,
        args.cache,
        records,
        tasks,
        stop_reason,
        resource,
    )
    summary = artifact["summary"]
    print(
        "stored Q208 top phase-cell certificate: "
        f"{summary['tasks_completed']}/{summary['tasks_total']} panels, "
        f"{summary['certified_leaf_cells']} cells, "
        f"{summary['unresolved_boxes']} unresolved, "
        f"complete={summary['complete_top_edge_certificate']}, "
        f"stop={summary['stop_reason']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
