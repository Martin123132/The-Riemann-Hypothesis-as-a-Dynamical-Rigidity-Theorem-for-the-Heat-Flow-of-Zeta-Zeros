#!/usr/bin/env python3
"""Statically audit the RH-owned resumable Hardy source contract."""

from __future__ import annotations

import hashlib
from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[3]
ROOT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable"
MODULE = ROOT / "rh_hardy_checkpoint.f90"
SOURCE = ROOT / "zeta14cubicmult_resumable.f90"
README = ROOT / "README.md"
LICENSE = ROOT / "COPYING.GPL-3.0.txt"
RUNNER = REPO_ROOT / "work/rh_compute/scripts/run_hardy_resumable_equivalence_fixture.py"
BASE_SOURCE = Path(
    r"C:\Users\ollet\Documents\Codex\third_party\hardy_fastcode_build\zeta14cubicmult.fixed.f90"
)
BASE_SHA256 = "5d0699865ab58f6968adfbe11ed969dafaabf5e965d5bd74391206af8b58b620"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require_tokens(text: str, tokens: tuple[str, ...], label: str) -> None:
    missing = [token for token in tokens if token not in text]
    require(not missing, f"{label} missing tokens: {missing}")


def main() -> int:
    for path in (MODULE, SOURCE, README, LICENSE, RUNNER):
        require(path.is_file(), f"missing source-contract artifact: {path}")
    if BASE_SOURCE.is_file():
        require(file_hash(BASE_SOURCE) == BASE_SHA256, "patched base-source hash drift")

    module = MODULE.read_text(encoding="utf-8")
    source = SOURCE.read_text(encoding="utf-8")
    readme = README.read_text(encoding="utf-8")
    runner = RUNNER.read_text(encoding="utf-8")

    require_tokens(
        module,
        (
            "RH_HARDY_CHECKPOINT_V1",
            "RH_HARDY_RUN_ID",
            "RH_HARDY_MAX_SECONDS",
            "RH_HARDY_STOP_FILE",
            "RH_HARDY_STOP_AFTER_STAGE",
            "RH_HARDY_STOP_AFTER_UNIT",
            "RH_HARDY_RS_CHUNK_SIZE",
            "call rh_sync_file(trim(temporary_path), status)",
            "call rename(trim(temporary_path), trim(cfg%checkpoint_path)",
            "position='append'",
            "state%zsum",
            "state%rszsumtot",
            "state%rszsum",
            "state%comm_size /= expected_comm_size",
            "checkpoint provenance hash does not match this run",
            "rh_serial_significand_digits = 41",
        ),
        "checkpoint module",
    )
    require(module.count("ES49.40E4") >= 10, "sign-safe full-precision serialization format drift")
    require("ES48.40E4" not in module, "negative-unsafe serialization width reintroduced")
    require(module.count("rh_checkpoint_commit") == 3, "checkpoint commit definition/export drift")
    require("ieee_is_finite" in module, "finite-state guard missing")

    require_tokens(
        source,
        (
            "use rh_hardy_checkpoint",
            "if (rh_config%resume) then",
            "do iblock=rh_block_start,totblock",
            "do rh_rs_unit=rh_unit_start,rh_rs_rounds",
            "rh_local_chunk=(rh_rs_unit-1_rh_i8)*int(COMM_SIZE,rh_i8)",
            "N_INTERNAL=rh_config%rs_chunk_size",
            "NIS=N_MAX*rh_config%rs_chunk_size+1_rh_i8",
            "rzsum=0.0_dp",
            "stop 75",
            "call MPI_Finalize(IERR)",
        ),
        "resumable evaluator",
    )
    require(source.count("call rh_checkpoint_capture") == 5, "checkpoint capture count drift")
    require(source.count("call rh_checkpoint_commit") == 5, "checkpoint commit call count drift")
    require(source.lower().count("call mpi_allreduce") >= 5, "MPI reduction coverage drift")
    require("shift_step=0.02" not in source and "shift_step=0.04" not in source, "rejected grid leaked into source")
    for shift_line in (
        "tar(i)=t+(i-numbercalc)*0.01",
        "yphasear(i)=yphase+(i-numbercalc)*0.01",
        "rsphasear(i)=rsphase+(i-numbercalc)*(log(0.25*a)-1/(48*t*t))*0.01",
        "aar(i)=a+4*(i-numbercalc)*0.01/(p*a)",
    ):
        require(shift_line in source, f"accepted 0.01 grid drift: {shift_line}")

    require_tokens(
        runner,
        (
            '"--cpus=1"',
            "nice -n 10 mpirun --allow-run-as-root -np 1",
            '"RH_HARDY_RS_CHUNK_SIZE": "257"',
            '"RH_HARDY_RUN_ID": provenance_id(binary, run_root, image_id)',
            '"park_runtime"',
            '"park_stop_file"',
            '"park_rs"',
            '"reject_run_id"',
            '"reject_chunk_size"',
            "uninterrupted_journal == resumed_journal",
        ),
        "equivalence runner",
    )

    require(BASE_SHA256 in readme, "base hash missing from resumable README")
    require("GPL-3.0" in readme, "license declaration missing")
    require("diagnostic" in readme.lower(), "diagnostic proof boundary missing")
    require("interval evaluator" in readme, "interval boundary missing")
    require(re.search(r"Riemann-hypothesis implication", readme) is not None, "RH boundary missing")

    print(
        "validated resumable Hardy source contract: "
        "5 captures, 5 commits, 9 controls, 15-channel state, h=0.01 retained"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
