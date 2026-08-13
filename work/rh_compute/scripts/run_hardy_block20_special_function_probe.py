#!/usr/bin/env python3
"""Probe the source PSI/ERF routines at all recursive block-20 calls."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
TELEMETRY = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full/enabled/chain_telemetry.jsonl"
)
FIXTURE_ROOT = (
    REPO_ROOT
    / "work/rh_compute/external/hardy_fastcode/fixtures/special_function_probe/t1e10_block20"
)
PROBE_SOURCE = FIXTURE_ROOT / "probe_generated.f90"
PROBE_INPUT = FIXTURE_ROOT / "probe_input.txt"
PROBE_OUTPUT = FIXTURE_ROOT / "probe_output.txt"
PROBE_STDERR = FIXTURE_ROOT / "probe.stderr.txt"
RESULT = FIXTURE_ROOT / "probe_result.json"
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_hardy_block20_special_function_probe.py"
DEFAULT_BUILD_ROOT = Path(
    r"C:\Users\ollet\Documents\Codex\third_party\hardy_fastcode_build\special_function_probe_20260806"
)
IMAGE = "rh-hardy-fastcode:bookworm"
SOURCE_SHA256 = "0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d"
HEX128 = re.compile(r"[0-9A-F]{32}\Z")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def docker_mount(path: Path) -> str:
    return str(path.resolve())


def run_command(command: list[str]) -> subprocess.CompletedProcess[str]:
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    completed.elapsed_seconds = time.perf_counter() - started  # type: ignore[attr-defined]
    require(
        completed.returncode == 0,
        f"command exited {completed.returncode}: {' '.join(command)}\n"
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
    )
    return completed


def extract_inclusive(source: str, start_token: str, end_token: str) -> str:
    start = source.index(start_token)
    end = source.index(end_token, start) + len(end_token)
    while end < len(source) and source[end] in " \t":
        end += 1
    if end < len(source) and source[end] == "\r":
        end += 1
    if end < len(source) and source[end] == "\n":
        end += 1
    return source[start:end]


def source_derivative() -> tuple[str, dict[str, str]]:
    source = SOURCE.read_text(encoding="utf-8")
    init_start = "! Set p=pi, tpp=2*pi and tpm=-2*pi"
    init_end = "!     ** All of the above code starting from ** can be excluded"
    init_begin = source.index(init_start)
    init_stop = source.index(init_end, init_begin)
    initialization = source[init_begin:init_stop]
    erf_routine = extract_inclusive(source, "SUBROUTINE ERF(Z,N,ER)", "END SUBROUTINE ERF")
    psi_routine = extract_inclusive(source, "SUBROUTINE PSI(X,PS)", "END SUBROUTINE PSI")

    main = f"""program rh_special_function_probe
  implicit none

  integer, parameter :: dp = selected_real_kind(33)
  integer, parameter :: i128 = selected_int_kind(38)
  integer :: chain, child_length, ios, i, ifact, jbot, mmax
  real(dp) :: p, sp, tpp, tpm, c, gam, zet(3:13)
  real(dp) :: a1, xr, fracL, sx, con1, con2, con3, sav1, sav2, z
  real(dp) :: psi_args(6), psi_values(6)
  complex(kind=16) :: coeff(7,0:6), zlarcoeff(-1:6), zsmalcoeff(-1:6)
  complex(kind=16) :: c1, c2, c3, c4, c5, c6, c7, epi4
  complex(kind=16) :: endpoint, erf_values(2), erf_weights(2)
  complex(kind=16) :: fn, cr2, cr3, source_t1, source_t2, source_t4
  character(len=32) :: h_a1, h_xr, h_frac, h_endpoint_re, h_endpoint_im

  common /PARMS/ coeff,zlarcoeff,zsmalcoeff
  common /PARS2/ gam,zet
  common /PARS4/ p,sp,tpp,tpm,epi4

{initialization}
  write(*,'(A)',advance='no') '#constants'
  write(*,'(1X,A)',advance='no') real_hex(p)
  write(*,'(1X,A)',advance='no') real_hex(sp)
  write(*,'(1X,A)',advance='no') real_hex(tpp)
  write(*,'(1X,A)',advance='no') real_hex(tpm)
  write(*,'(1X,A)',advance='no') real_hex(real(epi4,kind=dp))
  write(*,'(1X,A)') real_hex(aimag(epi4))

  do
     read(*,*,iostat=ios) chain, child_length, h_a1, h_xr, h_frac, h_endpoint_re, h_endpoint_im
     if (ios.lt.0) exit
     if (ios.ne.0) error stop 'special-function probe input failure'
     call real_from_hex(h_a1,a1)
     call real_from_hex(h_xr,xr)
     call real_from_hex(h_frac,fracL)
     call complex_from_hex(h_endpoint_re,h_endpoint_im,endpoint)

     jbot=ceiling(a1)
     psi_args(1)=child_length+fracL+1.0
     psi_args(2)=1.0-fracL
     psi_args(3)=a1+1.0
     psi_args(4)=child_length+1.0-a1
     psi_args(5)=jbot-a1
     psi_args(6)=jbot-(child_length+fracL)
     do i=1,6
        call psi(psi_args(i),psi_values(i))
     enddo

     sav1=psi_values(2)
     con2=psi_values(1)-psi_values(2)
     sav2=psi_values(4)
     con3=psi_values(3)-psi_values(4)
     source_t1=(0.0,1.0)*(conjg(endpoint)*con2-con3)/tpp

     sx=sqrt(xr)
     con1=1.0-fracL
     z=sp*con1/sx
     call erf(z,6,erf_values(1))
     cr2=1.0-erf_values(1)
     fn=-(0.0,1.0)*p*(con1**2)/xr
     cr3=p*exp(fn)*cr2/(sx*epi4)-1/con1-1/(child_length*1.0+1.0)
     source_t2=(0.0,1.0)*conjg(endpoint)*cr3/tpp
     erf_weights(1)=conjg(-(0.0,1.0)*conjg(endpoint)*p*exp(fn)/(sx*epi4*tpp))

     if (a1.gt.0.0) then
        con1=a1
     else
        con1=a1+1.0
     endif
     z=sp*con1/sx
     call erf(z,6,erf_values(2))
     cr2=1.0-erf_values(2)
     fn=-(0.0,1.0)*p*(con1**2)/xr
     cr3=(0.0,1.0)*exp(fn)*cr2/(2*sx*epi4)
     if (a1.gt.0.0) then
        source_t4=-(0.0,1.0)*conjg(endpoint)/(tpp*(child_length+fracL))+cr3
     else
        source_t4=(0.0,1.0)*(a1/con1-1.0)/tpp+cr3
     endif
     erf_weights(2)=conjg(-(0.0,1.0)*exp(fn)/(2*sx*epi4))

     write(*,'(I0)',advance='no') chain
     do i=1,6
        write(*,'(1X,A,1X,A)',advance='no') real_hex(psi_args(i)),real_hex(psi_values(i))
     enddo
     write(*,'(1X,A,1X,A,1X,A)',advance='no') real_hex(sp*(1.0-fracL)/sx), &
          real_hex(real(erf_values(1),kind=dp)),real_hex(aimag(erf_values(1)))
     write(*,'(1X,A,1X,A,1X,A)',advance='no') real_hex(z), &
          real_hex(real(erf_values(2),kind=dp)),real_hex(aimag(erf_values(2)))
     do i=1,2
        write(*,'(1X,A,1X,A)',advance='no') real_hex(real(erf_weights(i),kind=dp)), &
             real_hex(aimag(erf_weights(i)))
     enddo
     write(*,'(1X,A,1X,A)',advance='no') real_hex(real(source_t1,kind=dp)),real_hex(aimag(source_t1))
     write(*,'(1X,A,1X,A)',advance='no') real_hex(real(source_t2,kind=dp)),real_hex(aimag(source_t2))
     write(*,'(1X,A,1X,A)') real_hex(real(source_t4,kind=dp)),real_hex(aimag(source_t4))
  enddo

contains

  function real_hex(value) result(text)
    real(dp), intent(in) :: value
    character(len=32) :: text
    integer(i128) :: bits
    bits=transfer(value,bits)
    write(text,'(Z32.32)') bits
  end function real_hex

  subroutine real_from_hex(text,value)
    character(len=32), intent(in) :: text
    real(dp), intent(out) :: value
    integer(i128) :: bits
    read(text,'(Z32)') bits
    value=transfer(bits,value)
  end subroutine real_from_hex

  subroutine complex_from_hex(real_text,imag_text,value)
    character(len=32), intent(in) :: real_text,imag_text
    complex(kind=16), intent(out) :: value
    real(dp) :: real_part,imag_part
    call real_from_hex(real_text,real_part)
    call real_from_hex(imag_text,imag_part)
    value=cmplx(real_part,imag_part,kind=16)
  end subroutine complex_from_hex

end program rh_special_function_probe

"""
    derivative = main + erf_routine + "\n" + psi_routine + "\n"
    hashes = {
        "initialization_sha256": hashlib.sha256(initialization.encode()).hexdigest(),
        "erf_routine_sha256": hashlib.sha256(erf_routine.encode()).hexdigest(),
        "psi_routine_sha256": hashlib.sha256(psi_routine.encode()).hexdigest(),
    }
    return derivative, hashes


def load_recursive_rows() -> tuple[list[str], dict[int, dict[str, Any]]]:
    chains: dict[int, dict[str, Any]] = {}
    for line in TELEMETRY.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        chain = int(record.get("chain", 0))
        if record["type"] == "chain":
            chains[chain] = {"header": record, "levels": {}, "q_terms": {}}
        elif record["type"] == "level":
            chains[chain]["levels"][int(record["level"])] = record
        elif record["type"] == "q_terms":
            chains[chain]["q_terms"][int(record["nit"])] = record

    recursive = {chain: data for chain, data in chains.items() if int(data["header"]["mit"]) == 2}
    require(len(recursive) == 374, "recursive probe roster drift")
    input_lines: list[str] = []
    for chain, data in sorted(recursive.items()):
        level_one = data["levels"][1]
        child_length = int(data["levels"][2]["length"])
        endpoint = data["q_terms"][1]["endpoint_hex"]
        fields = [
            str(chain),
            str(child_length),
            level_one["coefficients_hex"][0],
            level_one["xr_hex"],
            level_one["frac_length_hex"],
            endpoint[0],
            endpoint[1],
        ]
        require(all(HEX128.fullmatch(value) for value in fields[2:]), f"chain {chain} input hex drift")
        input_lines.append(" ".join(fields))
    return input_lines, recursive


def parse_output(text: str, recursive: dict[int, dict[str, Any]]) -> dict[str, Any]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    require(lines and lines[0].startswith("#constants "), "probe constants header missing")
    constants = lines[0].split()
    require(len(constants) == 7 and all(HEX128.fullmatch(value) for value in constants[1:]), "probe constants drift")
    rows = lines[1:]
    require(len(rows) == 374, "probe output row count drift")
    bit_matches = 0
    for line, chain in zip(rows, sorted(recursive)):
        fields = line.split()
        require(len(fields) == 29, f"chain {chain} probe field count drift")
        require(int(fields[0]) == chain, f"chain {chain} probe order drift")
        require(all(HEX128.fullmatch(value) for value in fields[1:]), f"chain {chain} probe output hex drift")
        q_terms = recursive[chain]["q_terms"][1]
        for offset, name in ((23, "t1_hex"), (25, "t2_hex"), (27, "t4_hex")):
            require(fields[offset : offset + 2] == q_terms[name], f"chain {chain} {name} replay drift")
            bit_matches += 1
    return {
        "row_count": len(rows),
        "psi_call_count": len(rows) * 6,
        "erf_call_count": len(rows) * 2,
        "saved_component_bit_match_count": bit_matches,
        "constants_hex": {
            "pi": constants[1],
            "sqrt_pi": constants[2],
            "two_pi": constants[3],
            "minus_two_pi": constants[4],
            "epi4_real": constants[5],
            "epi4_imag": constants[6],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-root", type=Path, default=DEFAULT_BUILD_ROOT)
    args = parser.parse_args()
    for path in (SOURCE, TELEMETRY, CHECKER):
        require(path.is_file(), f"missing probe dependency: {path}")
    require(file_hash(SOURCE) == SOURCE_SHA256, "accepted source hash drift")

    FIXTURE_ROOT.mkdir(parents=True, exist_ok=True)
    args.build_root.mkdir(parents=True, exist_ok=True)
    derivative, extraction_hashes = source_derivative()
    PROBE_SOURCE.write_text(derivative, encoding="utf-8")
    input_lines, recursive = load_recursive_rows()
    PROBE_INPUT.write_text("\n".join(input_lines) + "\n", encoding="ascii")

    image = run_command(["docker", "image", "inspect", "--format={{.Id}}", IMAGE]).stdout.strip()
    require(re.fullmatch(r"sha256:[0-9a-f]{64}", image) is not None, "probe image-id drift")
    binary = args.build_root / "rh_special_function_probe"
    build = run_command(
        [
            "docker", "run", "--rm", "--cpus=1", "--network=none",
            "-v", f"{docker_mount(FIXTURE_ROOT)}:/fixture:ro",
            "-v", f"{docker_mount(args.build_root)}:/build",
            "-w", "/build", IMAGE, "bash", "-lc",
            "nice -n 10 gfortran -O3 -ffree-line-length-none /fixture/probe_generated.f90 -o /build/rh_special_function_probe",
        ]
    )
    require(binary.is_file(), "probe binary was not produced")
    run = run_command(
        [
            "docker", "run", "--rm", "--cpus=1", "--network=none",
            "-v", f"{docker_mount(FIXTURE_ROOT)}:/fixture:ro",
            "-v", f"{docker_mount(args.build_root)}:/build:ro",
            "-w", "/build", IMAGE, "bash", "-lc",
            "nice -n 10 ./rh_special_function_probe < /fixture/probe_input.txt",
        ]
    )
    PROBE_OUTPUT.write_text(run.stdout, encoding="ascii")
    PROBE_STDERR.write_text(run.stderr, encoding="utf-8")
    validation = parse_output(run.stdout, recursive)

    artifact = {
        "kind": "hardy_block20_special_function_probe",
        "status": "source_derived_374_call_probe_with_exact_t1_t2_t4_replay",
        "block": 20,
        "recursive_chain_count": 374,
        "validation": validation,
        "build": {
            "image": IMAGE,
            "image_id": image,
            "cpu_cap": 1,
            "nice": 10,
            "compiler_flags": ["-O3", "-ffree-line-length-none"],
            "elapsed_seconds": build.elapsed_seconds,  # type: ignore[attr-defined]
            "binary_sha256": file_hash(binary),
        },
        "run": {
            "cpu_cap": 1,
            "nice": 10,
            "elapsed_seconds": run.elapsed_seconds,  # type: ignore[attr-defined]
        },
        "source_extraction": extraction_hashes,
        "artifacts": {
            "generated_probe": {"path": relative(PROBE_SOURCE), "sha256": file_hash(PROBE_SOURCE)},
            "input": {"path": relative(PROBE_INPUT), "sha256": file_hash(PROBE_INPUT)},
            "output": {"path": relative(PROBE_OUTPUT), "sha256": file_hash(PROBE_OUTPUT)},
            "stderr": {"path": relative(PROBE_STDERR), "sha256": file_hash(PROBE_STDERR)},
        },
        "sources": {
            "accepted_source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE)},
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "runner": {"path": relative(Path(__file__).resolve()), "sha256": file_hash(Path(__file__).resolve())},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "The standalone one-CPU probe copies the accepted initialization and PSI/ERF routines and "
            "replays saved t1, t2, and t4 bits at 374 finite recursive calls. It does not itself prove "
            "special-function accuracy, uniform cell bounds, recurrence accuracy, RH, or a prize-level theorem."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        "built Hardy block-20 special-function probe: "
        "374 calls, 2992 source values, 1122 exact t1/t2/t4 bit replays"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
