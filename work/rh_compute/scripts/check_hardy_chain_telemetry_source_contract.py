#!/usr/bin/env python3
"""Audit the non-invasive source boundary for Hardy chain telemetry."""

from __future__ import annotations

import hashlib
from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[3]
ACCEPTED = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/resumable"
    / "zeta14cubicmult_resumable.f90"
)
TELEMETRY_ROOT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/telemetry"
DERIVATIVE = TELEMETRY_ROOT / "zeta14cubicmult_telemetry.f90"
MODULE = TELEMETRY_ROOT / "rh_hardy_chain_telemetry.f90"
README = TELEMETRY_ROOT / "README.md"
ACCEPTED_SHA256 = "0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d"

BEGIN = re.compile(r"^\s*! RH_TELEMETRY_BEGIN ([a-z0-9_]+)\s*$")
END = re.compile(r"^\s*! RH_TELEMETRY_END ([a-z0-9_]+)\s*$")

EXPECTED_BLOCKS = {
    "program_use": (
        "  use rh_hardy_chain_telemetry, only : rh_tel_init, rh_tel_close, rh_tel_set_position",
    ),
    "program_init": ("  call rh_tel_init(RANK)",),
    "program_position": ("        call rh_tel_set_position(iblock,jsum)",),
    "program_close": ("  call rh_tel_close()",),
    "alphasum_use": (
        "  use rh_hardy_chain_telemetry, only : rh_tel_set_branch, rh_tel_set_initial_coefficients",
    ),
    "alphasum_branch_one": (
        "     call rh_tel_set_branch(1)",
        "     call rh_tel_set_initial_coefficients(mc,initialcoeff)",
    ),
    "alphasum_branch_two": (
        "     call rh_tel_set_branch(2)",
        "     call rh_tel_set_initial_coefficients(mc,initialcoeffcc)",
    ),
    "mgausssum_use": (
        "  use rh_hardy_chain_telemetry, only : rh_tel_begin_chain, rh_tel_step, rh_tel_end_chain",
    ),
    "mgausssum_locals": (
        "  complex(kind=16)   :: rh_tel_before,rh_tel_pre",
    ),
    "mgausssum_begin": (
        "  call rh_tel_begin_chain(m,ip,MIT,mmax,L,m1,phicoeff,fracL,xr,icj,tpm)",
    ),
    "mgausssum_step_before": ("     rh_tel_before=csum",),
    "mgausssum_step_pre": ("     rh_tel_pre=csum",),
    "mgausssum_step_after": (
        "     call rh_tel_step(k,rh_tel_before,c1,qq,rh_tel_pre,csum,L,m1,phicoeff,icj,tpm)",
    ),
    "mgausssum_end": ("  call rh_tel_end_chain(csum)",),
    "q_use": ("  use rh_hardy_chain_telemetry, only : rh_tel_q",),
    "q_terms": ("  call rh_tel_q(nit,t1,t2,t3,t4,t5,endpoint,qq)",),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def split_instrumented_source(text: str) -> tuple[list[str], dict[str, tuple[str, ...]]]:
    outside: list[str] = []
    blocks: dict[str, tuple[str, ...]] = {}
    active_name: str | None = None
    active_lines: list[str] = []

    for line_number, line in enumerate(text.splitlines(), start=1):
        begin = BEGIN.match(line)
        end = END.match(line)
        if begin:
            require(active_name is None, f"nested telemetry block at line {line_number}")
            active_name = begin.group(1)
            require(active_name not in blocks, f"duplicate telemetry block: {active_name}")
            active_lines = []
            continue
        if end:
            require(active_name is not None, f"orphan telemetry end at line {line_number}")
            require(end.group(1) == active_name, f"mismatched telemetry end at line {line_number}")
            blocks[active_name] = tuple(active_lines)
            active_name = None
            active_lines = []
            continue
        if active_name is None:
            outside.append(line)
        else:
            active_lines.append(line)

    require(active_name is None, f"unterminated telemetry block: {active_name}")
    return outside, blocks


def nonblank(lines: list[str]) -> list[str]:
    return [line for line in lines if line.strip()]


def require_tokens(text: str, tokens: tuple[str, ...], label: str) -> None:
    missing = [token for token in tokens if token not in text]
    require(not missing, f"{label} missing tokens: {missing}")


def main() -> int:
    for path in (ACCEPTED, DERIVATIVE, MODULE, README):
        require(path.is_file(), f"missing telemetry source-contract artifact: {path}")

    require(sha256(ACCEPTED) == ACCEPTED_SHA256, "accepted evaluator hash drift")

    accepted = ACCEPTED.read_text(encoding="utf-8")
    derivative = DERIVATIVE.read_text(encoding="utf-8")
    module = MODULE.read_text(encoding="utf-8")
    readme = README.read_text(encoding="utf-8")

    outside, blocks = split_instrumented_source(derivative)
    require(set(blocks) == set(EXPECTED_BLOCKS), "telemetry marker inventory drift")
    for name, expected in EXPECTED_BLOCKS.items():
        require(blocks[name] == expected, f"telemetry block body drift: {name}")

    require(
        nonblank(outside) == nonblank(accepted.splitlines()),
        "telemetry derivative changes accepted nonblank source outside marked observer blocks",
    )
    require("rh_hardy_chain_telemetry" not in accepted, "telemetry leaked into accepted source")

    require_tokens(
        module,
        (
            "RH_HARDY_CHAIN_TELEMETRY",
            "RH_HARDY_TELEMETRY_MAX_CHAINS",
            "RH_HARDY_TELEMETRY_MAX_DIRECT_TERMS",
            "',\"tpm\":\"' // trim(r_text(tpm))",
            "selected_int_kind(32)",
            "transfer(value,bits)",
            "write(text,'(Z32.32)') bits",
            "',\"coefficients_hex\":'",
            "',\"initial_coefficients_hex\":'",
            "',\"local_defect_hex\":'",
            "position='append'",
            "flush(telemetry_unit)",
            "schema\":\"rh_hardy_chain_telemetry_v1",
            "call direct_level_sum(nit",
            "call direct_level_sum(nit+1",
            "phase = (phase+phicoeff(j,level))*y",
            "exp(cmplx(0.0_dp,tpm*phase,kind=16))",
            "local_defect = adapted_parent - model_after",
            "child_state_defect = state_before - adapted_child",
            "parent_state_defect = state_after - adapted_parent",
        ),
        "telemetry module",
    )
    require(module.count("call write_record(line)") >= 6, "telemetry record coverage drift")
    require("stop " not in module.lower(), "telemetry module must not stop the evaluator")
    require("error stop" not in module.lower(), "telemetry module must fail inertly")
    require(ACCEPTED_SHA256 in readme, "accepted source hash missing from telemetry README")
    require("p=4*ATAN(C)" in readme, "phase-constant provenance missing from README")
    require("exact 32-hex-digit" in readme, "exact binary payload contract missing")
    require("not establish interval" in readme, "telemetry proof boundary missing")

    print(
        "validated Hardy telemetry source boundary: accepted hash pinned, "
        f"{len(blocks)} observer blocks exact, direct parent/child defects instrumented"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
