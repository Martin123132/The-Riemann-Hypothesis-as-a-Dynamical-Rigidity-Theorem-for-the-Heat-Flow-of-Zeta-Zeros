#!/usr/bin/env python3
"""Wait for the point cache, then complete and validate the order-nine chain."""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from time import monotonic, sleep


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
POINT_CACHE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_shifted_point_h0_h8_half_grid.jsonl"
)
POINT_MANIFEST = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_shifted_point_h0_h8_cache.json"
)
DEFAULT_STATUS = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_completion_pipeline_status.json"
)
DEFAULT_LOCK = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_completion_pipeline.lock"
)
EXPECTED_POINT_ROWS = 8929


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def timestamp() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def write_status(path: Path, status: dict) -> None:
    status["updated_at"] = timestamp()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(status, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def point_manifest_ready() -> tuple[bool, str]:
    if not POINT_MANIFEST.exists():
        return False, "manifest absent"
    try:
        manifest = json.loads(POINT_MANIFEST.read_text(encoding="utf-8"))
        cache = manifest["cache"]
        if cache.get("row_count") != EXPECTED_POINT_ROWS:
            return False, f"manifest rows={cache.get('row_count')}"
        if cache.get("all_rows_passed") is not True:
            return False, "manifest does not mark every row passed"
        if not POINT_CACHE.exists():
            return False, "cache absent"
        if sha256(POINT_CACHE) != cache.get("sha256"):
            return False, "cache hash does not match manifest"
    except Exception as exc:
        return False, f"manifest validation failed: {exc}"
    return True, "complete hash-bound 8929-row cache"


def wait_for_point_cache(
    status_path: Path,
    status: dict,
    *,
    poll_seconds: int,
    timeout_hours: float,
    stale_minutes: float,
) -> None:
    started = monotonic()
    last_size = -1
    last_change = started
    last_report = 0.0
    while True:
        ready, detail = point_manifest_ready()
        if ready:
            print(f"point cache ready: {detail}", flush=True)
            status["point_cache"] = {
                "state": "complete",
                "detail": detail,
                "bytes": POINT_CACHE.stat().st_size,
                "sha256": sha256(POINT_CACHE),
            }
            write_status(status_path, status)
            return
        now = monotonic()
        if now - started > timeout_hours * 3600:
            raise RuntimeError(
                f"point-cache wait exceeded {timeout_hours} hours: {detail}"
            )
        size = POINT_CACHE.stat().st_size if POINT_CACHE.exists() else 0
        if size != last_size:
            last_size = size
            last_change = now
        elif now - last_change > stale_minutes * 60:
            raise RuntimeError(
                f"point cache has not grown for {stale_minutes} minutes: {detail}"
            )
        if now - last_report >= max(60, poll_seconds):
            elapsed = (now - started) / 60
            print(
                f"waiting for point cache: {size} bytes, {elapsed:.1f} min, "
                f"{detail}",
                flush=True,
            )
            status["point_cache"] = {
                "state": "waiting",
                "detail": detail,
                "bytes": size,
                "elapsed_minutes": elapsed,
            }
            write_status(status_path, status)
            last_report = now
        sleep(poll_seconds)


def command_steps(workers: int, *, rebuild_lower: bool) -> list[tuple[str, list[str]]]:
    lower_checker = [
        "check_jensen_window_pf_compound_order9_localized_lower_bridge_certificate.py"
    ]
    steps = [
        (
            "check exact-point H0-H8 cache",
            [
                "check_jensen_window_pf_compound_order9_shifted_point_h0_h8_cache.py"
            ],
        ),
        (
            "build localized lower bridge",
            [
                "jensen_window_pf_compound_order9_localized_lower_bridge_certificate.py",
                "--workers",
                str(workers),
            ],
        ),
        ("check localized lower bridge", lower_checker),
    ]
    if rebuild_lower:
        steps.append(
            (
                "rebuild localized lower bridge",
                lower_checker + ["--rebuild", "--workers", str(workers)],
            )
        )
    steps.extend(
        [
            (
                "check compact upper bridge",
                [
                    "check_jensen_window_pf_compound_order9_nested_curvature_compact_certificate.py"
                ],
            ),
            (
                "check finite saddle ray",
                [
                    "check_jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate.py"
                ],
            ),
            (
                "check asymptotic saddle ray",
                [
                    "check_jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate.py"
                ],
            ),
            (
                "build global first-summand theorem",
                [
                    "jensen_window_pf_compound_order9_first_summand_curvature_certificate.py"
                ],
            ),
            (
                "check global first-summand theorem",
                [
                    "check_jensen_window_pf_compound_order9_first_summand_curvature_certificate.py"
                ],
            ),
            (
                "check full-kernel transfer",
                [
                    "check_jensen_window_pf_compound_order9_first_summand_curvature_bridge.py"
                ],
            ),
            (
                "check finite endpoint splice",
                [
                    "check_jensen_window_pf_compound_order9_m100_finite_splice_certificate.py"
                ],
            ),
            (
                "build lambda=-100 endpoint entry",
                ["jensen_window_pf_compound_order9_m100_entry_certificate.py"],
            ),
            (
                "check lambda=-100 endpoint entry",
                ["check_jensen_window_pf_compound_order9_m100_entry_certificate.py"],
            ),
            (
                "check uniform-tail flow reduction",
                [
                    "check_jensen_window_pf_compound_order9_uniform_tail_flow_reduction.py"
                ],
            ),
            (
                "check lower order-eight heat theorem",
                [
                    "check_jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate.py"
                ],
            ),
            (
                "build order-nine heat-forward theorem",
                [
                    "jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate.py"
                ],
            ),
            (
                "check order-nine heat-forward theorem",
                [
                    "check_jensen_window_pf_compound_order9_uniform_heat_forward_invariance_certificate.py"
                ],
            ),
        ]
    )
    return steps


def run_step(label: str, arguments: list[str]) -> dict:
    command = [sys.executable, str(SCRIPT_DIR / arguments[0]), *arguments[1:]]
    print(f"\n[{timestamp()}] START {label}", flush=True)
    started = monotonic()
    process = subprocess.Popen(
        command,
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        bufsize=1,
    )
    output = []
    assert process.stdout is not None
    for line in process.stdout:
        print(line, end="", flush=True)
        output.append(line.rstrip("\n"))
        if len(output) > 100:
            output.pop(0)
    return_code = process.wait()
    elapsed = monotonic() - started
    result = {
        "label": label,
        "command": command,
        "return_code": return_code,
        "elapsed_seconds": elapsed,
        "tail": output,
    }
    if return_code:
        raise RuntimeError(
            f"pipeline step failed ({return_code}): {label}"
        )
    print(f"[{timestamp()}] PASS {label} ({elapsed:.1f}s)", flush=True)
    return result


def acquire_lock(path: Path, *, force: bool) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    if force and path.exists():
        path.unlink()
    try:
        descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as exc:
        raise RuntimeError(f"completion pipeline lock already exists: {path}") from exc
    os.write(descriptor, f"pid={os.getpid()} started={timestamp()}\n".encode("ascii"))
    return descriptor


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--poll-seconds", type=int, default=60)
    parser.add_argument("--wait-timeout-hours", type=float, default=8.0)
    parser.add_argument("--stale-minutes", type=float, default=30.0)
    parser.add_argument("--skip-lower-rebuild", action="store_true")
    parser.add_argument("--status", type=Path, default=DEFAULT_STATUS)
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--force-lock", action="store_true")
    args = parser.parse_args()

    lock_descriptor = acquire_lock(args.lock, force=args.force_lock)
    status = {
        "kind": "jensen_window_pf_compound_order9_completion_pipeline_status",
        "started_at": timestamp(),
        "state": "running",
        "pid": os.getpid(),
        "steps": [],
    }
    write_status(args.status, status)
    try:
        wait_for_point_cache(
            args.status,
            status,
            poll_seconds=max(5, args.poll_seconds),
            timeout_hours=args.wait_timeout_hours,
            stale_minutes=args.stale_minutes,
        )
        for label, command in command_steps(
            max(1, args.workers),
            rebuild_lower=not args.skip_lower_rebuild,
        ):
            status["active_step"] = label
            write_status(args.status, status)
            result = run_step(label, command)
            status["steps"].append(result)
            write_status(args.status, status)
        status["state"] = "complete"
        status["completed_at"] = timestamp()
        status.pop("active_step", None)
        write_status(args.status, status)
        print("order-nine completion pipeline finished successfully", flush=True)
        return 0
    except Exception as exc:
        status["state"] = "failed"
        status["failure"] = str(exc)
        status["failed_at"] = timestamp()
        write_status(args.status, status)
        raise
    finally:
        os.close(lock_descriptor)
        if args.lock.exists():
            args.lock.unlink()


if __name__ == "__main__":
    raise SystemExit(main())
