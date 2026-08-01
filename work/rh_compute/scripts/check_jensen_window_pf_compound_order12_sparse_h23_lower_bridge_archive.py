#!/usr/bin/env python3
"""Validate the published gzip archive with the canonical segment checker."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS = REPO_ROOT / "work/rh_compute/results"
DEFAULT_ARCHIVE = (
    RESULTS
    / "jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments.jsonl.gz"
)
DEFAULT_RUN_CONTRACT = (
    RESULTS
    / "jensen_window_pf_compound_order12_sparse_h23_lower_bridge_run_contract.json"
)
CHECKER = (
    REPO_ROOT
    / "work/rh_compute/scripts/check_jensen_window_pf_compound_order12_sparse_h23_lower_bridge_segments.py"
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--run-contract", type=Path, default=DEFAULT_RUN_CONTRACT)
    args = parser.parse_args()

    if not args.archive.is_file():
        parser.error(f"archive does not exist: {args.archive}")

    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".jsonl", delete=False) as target:
            temporary_path = Path(target.name)
            with gzip.open(args.archive, "rb") as source:
                shutil.copyfileobj(source, target, length=8 * 1024 * 1024)

        completed = subprocess.run(
            [
                sys.executable,
                str(CHECKER),
                "--cache",
                str(temporary_path),
                "--run-contract",
                str(args.run_contract),
            ],
            cwd=REPO_ROOT,
            check=False,
            creationflags=getattr(subprocess, "BELOW_NORMAL_PRIORITY_CLASS", 0),
        )
        return completed.returncode
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
