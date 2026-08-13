#!/usr/bin/env python3
"""Probe every source-level t5 component at recursive block-20 calls."""

from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_CEILING
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
    / "work/rh_compute/external/hardy_fastcode/fixtures/t5_component_probe/t1e10_block20"
)
PROBE_SOURCE = FIXTURE_ROOT / "probe_generated.f90"
PROBE_INPUT = FIXTURE_ROOT / "probe_input.txt"
PROBE_OUTPUT = FIXTURE_ROOT / "probe_output.txt"
PROBE_STDERR = FIXTURE_ROOT / "probe.stderr.txt"
RESULT = FIXTURE_ROOT / "probe_result.json"
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_hardy_block20_t5_component_probe.py"
DEFAULT_BUILD_ROOT = Path(
    r"C:\Users\ollet\Documents\Codex\third_party\hardy_fastcode_build\t5_component_probe_20260806"
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


def extract_inclusive(source: str, start_token: str, end_token: str, start_at: int = 0) -> str:
    start = source.index(start_token, start_at)
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
    psi_routine = extract_inclusive(source, "SUBROUTINE PSI(X,PS)", "END SUBROUTINE PSI")
    q_start = source.index("subroutine q(m,nit,ip,qq)")
    t5_source_block = extract_inclusive(
        source,
        "  sum0=(0.0,0.0)",
        "  t5=t5+(SUM0+SUM1)",
        q_start,
    )

    main = f"""program rh_t5_component_probe
  implicit none

  integer, parameter :: dp = selected_real_kind(33)
  integer, parameter :: i128 = selected_int_kind(38)
  integer :: chain, child_length, parent_length, ip, degree, ios
  integer :: i, j, ifact, jbot, mmax, ip1, slot_count
  integer :: slot_active(3), slot_region(3), slot_index(3)
  real(dp) :: p, sp, tpp, tpm, c, gam, zet(3:13)
  real(dp) :: a1, a2, a3, xr, fracL, sx, con1, con2, con3, wm
  real(dp) :: beta, xbeta, c1q, ecor, gc, sav1, sav2, ps1, ps2, tol
  real(dp) :: erfc_args(3), erfc_values(3), phi(3)
  complex(kind=16) :: coeff(7,0:6), zlarcoeff(-1:6), zsmalcoeff(-1:6)
  complex(kind=16) :: c1, c2, c3, c4, c5, c6, c7, epi4
  complex(kind=16) :: endpoint, cr1, cr2, cr3, sum0, sum1, t5, t5psi
  complex(kind=16) :: erfc_weights(3), erfc_terms(3)
  character(len=32) :: h_a1, h_a2, h_a3, h_xr, h_frac
  character(len=32) :: h_endpoint_re, h_endpoint_im

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
     read(*,*,iostat=ios) chain, child_length, parent_length, ip, degree, &
          h_a1, h_a2, h_a3, h_xr, h_frac, h_endpoint_re, h_endpoint_im
     if (ios.lt.0) exit
     if (ios.ne.0) error stop 't5 component probe input failure'
     if (degree.ne.3) error stop 't5 component probe degree drift'
     call real_from_hex(h_a1,a1)
     call real_from_hex(h_a2,a2)
     call real_from_hex(h_a3,a3)
     call real_from_hex(h_xr,xr)
     call real_from_hex(h_frac,fracL)
     call complex_from_hex(h_endpoint_re,h_endpoint_im,endpoint)
     phi(1)=a1
     phi(2)=a2
     phi(3)=a3

     ! Replay the source state inherited by the t5 block.
     sx=sqrt(xr)
     con1=1.0-fracL
     con2=child_length+fracL+1.0
     call psi(con2,ps2)
     call psi(con1,ps1)
     sav1=ps1
     con2=ps2-ps1
     call psi(a1+1.0,ps2)
     call psi(child_length+1.0-a1,ps1)
     sav2=ps1
     con3=ps2-ps1

     slot_count=0
     slot_active=0
     slot_region=0
     slot_index=0
     erfc_args=0.0
     erfc_values=0.0
     erfc_weights=(0.0,0.0)
     erfc_terms=(0.0,0.0)

     sum0=(0.0,0.0)
     jbot=ceiling(a1)

     con1=jbot-a1
     call psi(con1,ps1)
     sav2=sav2-ps1
     con2=jbot-(child_length+fracL)
     call psi(con2,ps2)
     sav1=ps2-sav1
     t5=(0.0,1.0)*(sav2+endpoint*sav1)/tpp
     t5psi=t5

     ip1=ip
     if (ip1.gt.(0.5*child_length)) then
        ip1=max1(0.2*child_length,1.0)
     endif
     do i=jbot,ip1
        con1=i*1.0-a1
        wm=3*a3*(con1**2)/(xr**3)
        c=con1/xr-wm
        sav1=sp*sx*c*sqrt(2.0)
        slot_count=slot_count+1
        if (slot_count.gt.3) error stop 't5 component probe slot overflow'
        slot_active(slot_count)=1
        slot_region(slot_count)=0
        slot_index(slot_count)=i
        erfc_args(slot_count)=sav1
        sav1=erfc(sav1)
        erfc_values(slot_count)=sav1
        gc=i*c
        do j=1,degree
           gc=gc-phi(j)*(c**j)
        enddo
        cr1=(0.0,1.0)*tpp*gc
        cr1=exp(cr1)
        sav2=exp(tpm*con1*c)
        cr2=(0.0,1.0)*a2/p
        cr3=(1/con1+cr2*(1+tpp*c*con1*(1.0+p*con1*c))/(con1**3))
        erfc_weights(slot_count)=-cr1/(2*sx*epi4)
        erfc_terms(slot_count)=-cr1*sav1/(2*sx*epi4)
        sum0=sum0+(-cr1*sav1/(2*sx*EPI4)-(0.0,1.0)*(sav2*cr3-cr2/(con1**3))/tpp)
     enddo

     sum1=(0.0,0.0)
     c1q=child_length+fracL
     ecor=0.0
     do i=child_length,child_length-ip1+1,-1
        con1=i*1.0-a1
        wm=3*a3*(con1**2)/(xr**3)
        c=con1/xr-wm
        beta=c/parent_length
        xbeta=parent_length*(1-beta)
        sav1=sp*sx*xbeta*sqrt(2.0)
        slot_count=slot_count+1
        if (slot_count.gt.3) error stop 't5 component probe slot overflow'
        slot_active(slot_count)=1
        slot_region(slot_count)=1
        slot_index(slot_count)=i
        erfc_args(slot_count)=sav1
        sav1=erfc(sav1)
        erfc_values(slot_count)=sav1
        gc=i*c
        do j=1,degree
           gc=gc-phi(j)*(c**j)
        enddo
        cr1=(0.0,1.0)*tpp*gc
        cr1=exp(cr1)
        cr2=(0.0,1.0)*a2/p
        sav2=exp(tpm*xbeta*(c1q-i))
        cr3=(1/(c1q-i)+cr2*(1+tpp*xbeta*(c1q-i)*(1.0+p*(c1q-i)*xbeta))/((c1q-i)**3))
        erfc_weights(slot_count)=-cr1/(2*sx*epi4)
        erfc_terms(slot_count)=(-sav1)*cr1/(2*sx*epi4)
        sum1=sum1+(-sav1)*cr1/(2*sx*EPI4)
        sum1=sum1-((0.0,1.0)*endpoint*(sav2*cr3-cr2/((c1q-i)**3))/tpp)
     enddo

     t5=t5+(sum0+sum1)

     write(*,'(I0)',advance='no') chain
     write(*,'(1X,A,1X,A)',advance='no') real_hex(real(t5psi,kind=dp)),real_hex(aimag(t5psi))
     write(*,'(1X,A,1X,A)',advance='no') real_hex(real(sum0,kind=dp)),real_hex(aimag(sum0))
     write(*,'(1X,A,1X,A)',advance='no') real_hex(real(sum1,kind=dp)),real_hex(aimag(sum1))
     write(*,'(1X,A,1X,A)',advance='no') real_hex(real(t5,kind=dp)),real_hex(aimag(t5))
     do i=1,3
        if (i.lt.3) then
           write(*,'(1X,I0,1X,I0,1X,I0,6(1X,A))',advance='no') &
                slot_active(i),slot_region(i),slot_index(i),real_hex(erfc_args(i)), &
                real_hex(erfc_values(i)),real_hex(real(erfc_weights(i),kind=dp)), &
                real_hex(aimag(erfc_weights(i))),real_hex(real(erfc_terms(i),kind=dp)), &
                real_hex(aimag(erfc_terms(i)))
        else
           write(*,'(1X,I0,1X,I0,1X,I0,6(1X,A))') &
                slot_active(i),slot_region(i),slot_index(i),real_hex(erfc_args(i)), &
                real_hex(erfc_values(i)),real_hex(real(erfc_weights(i),kind=dp)), &
                real_hex(aimag(erfc_weights(i))),real_hex(real(erfc_terms(i),kind=dp)), &
                real_hex(aimag(erfc_terms(i)))
        endif
     enddo
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

end program rh_t5_component_probe

"""
    derivative = main + psi_routine + "\n"
    hashes = {
        "initialization_sha256": hashlib.sha256(initialization.encode()).hexdigest(),
        "psi_routine_sha256": hashlib.sha256(psi_routine.encode()).hexdigest(),
        "t5_source_block_sha256": hashlib.sha256(t5_source_block.encode()).hexdigest(),
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
    require(len(recursive) == 374, "recursive t5-probe roster drift")
    input_lines: list[str] = []
    for chain, data in sorted(recursive.items()):
        header = data["header"]
        level_one = data["levels"][1]
        child_length = int(data["levels"][2]["length"])
        endpoint = data["q_terms"][1]["endpoint_hex"]
        coefficients = level_one["coefficients_hex"]
        require(len(coefficients) == 3, f"chain {chain} degree drift")
        fields = [
            str(chain),
            str(child_length),
            str(int(header["initial_length"])),
            str(int(header["ip"])),
            str(int(header["base_degree"])),
            *coefficients,
            level_one["xr_hex"],
            level_one["frac_length_hex"],
            endpoint[0],
            endpoint[1],
        ]
        require(all(HEX128.fullmatch(value) for value in fields[5:]), f"chain {chain} input hex drift")
        input_lines.append(" ".join(fields))
    return input_lines, recursive


def expected_slots(data: dict[str, Any]) -> list[tuple[int, int]]:
    header = data["header"]
    level_one = data["levels"][1]
    child_length = int(data["levels"][2]["length"])
    jbot = int(Decimal(level_one["coefficients"][0]).to_integral_value(rounding=ROUND_CEILING))
    ip1 = int(header["ip"])
    if ip1 > 0.5 * child_length:
        ip1 = int(max(Decimal("0.2") * child_length, Decimal(1)))
    return [*( (0, i) for i in range(jbot, ip1 + 1) ), *( (1, i) for i in range(child_length, child_length - ip1, -1) )]


def parse_output(text: str, recursive: dict[int, dict[str, Any]]) -> dict[str, Any]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    require(lines and lines[0].startswith("#constants "), "t5 probe constants header missing")
    constants = lines[0].split()
    require(len(constants) == 7 and all(HEX128.fullmatch(value) for value in constants[1:]), "t5 probe constants drift")
    rows = lines[1:]
    require(len(rows) == 374, "t5 probe output row count drift")
    bit_matches = 0
    lower_calls = 0
    upper_calls = 0
    two_call_rows = 0
    three_call_rows = 0
    for line, chain in zip(rows, sorted(recursive)):
        fields = line.split()
        require(len(fields) == 36, f"chain {chain} t5-probe field count drift")
        require(int(fields[0]) == chain, f"chain {chain} t5-probe order drift")
        require(all(HEX128.fullmatch(value) for value in fields[1:9]), f"chain {chain} component hex drift")
        require(fields[7:9] == recursive[chain]["q_terms"][1]["t5_hex"], f"chain {chain} t5 replay drift")
        bit_matches += 1

        wanted = expected_slots(recursive[chain])
        require(len(wanted) in (2, 3), f"chain {chain} unexpected source erfc roster")
        two_call_rows += len(wanted) == 2
        three_call_rows += len(wanted) == 3
        for slot in range(3):
            offset = 9 + 9 * slot
            active = int(fields[offset])
            region = int(fields[offset + 1])
            index = int(fields[offset + 2])
            require(all(HEX128.fullmatch(value) for value in fields[offset + 3 : offset + 9]), f"chain {chain} slot {slot + 1} hex drift")
            if slot < len(wanted):
                require(active == 1, f"chain {chain} active slot {slot + 1} drift")
                require((region, index) == wanted[slot], f"chain {chain} slot {slot + 1} identity drift")
                lower_calls += region == 0
                upper_calls += region == 1
            else:
                require((active, region, index) == (0, 0, 0), f"chain {chain} inactive slot drift")
    return {
        "row_count": len(rows),
        "saved_t5_bit_match_count": bit_matches,
        "intrinsic_erfc_call_count": lower_calls + upper_calls,
        "lower_erfc_call_count": lower_calls,
        "upper_erfc_call_count": upper_calls,
        "two_call_row_count": two_call_rows,
        "three_call_row_count": three_call_rows,
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
        require(path.is_file(), f"missing t5-probe dependency: {path}")
    require(file_hash(SOURCE) == SOURCE_SHA256, "accepted source hash drift")

    FIXTURE_ROOT.mkdir(parents=True, exist_ok=True)
    args.build_root.mkdir(parents=True, exist_ok=True)
    derivative, extraction_hashes = source_derivative()
    PROBE_SOURCE.write_text(derivative, encoding="utf-8")
    input_lines, recursive = load_recursive_rows()
    PROBE_INPUT.write_text("\n".join(input_lines) + "\n", encoding="ascii")

    image = run_command(["docker", "image", "inspect", "--format={{.Id}}", IMAGE]).stdout.strip()
    require(re.fullmatch(r"sha256:[0-9a-f]{64}", image) is not None, "t5-probe image-id drift")
    binary = args.build_root / "rh_t5_component_probe"
    build = run_command(
        [
            "docker", "run", "--rm", "--cpus=1", "--network=none",
            "-v", f"{docker_mount(FIXTURE_ROOT)}:/fixture:ro",
            "-v", f"{docker_mount(args.build_root)}:/build",
            "-w", "/build", IMAGE, "bash", "-lc",
            "nice -n 10 gfortran -O3 -ffree-line-length-none /fixture/probe_generated.f90 -o /build/rh_t5_component_probe",
        ]
    )
    require(binary.is_file(), "t5-probe binary was not produced")
    run = run_command(
        [
            "docker", "run", "--rm", "--cpus=1", "--network=none",
            "-v", f"{docker_mount(FIXTURE_ROOT)}:/fixture:ro",
            "-v", f"{docker_mount(args.build_root)}:/build:ro",
            "-w", "/build", IMAGE, "bash", "-lc",
            "nice -n 10 ./rh_t5_component_probe < /fixture/probe_input.txt",
        ]
    )
    PROBE_OUTPUT.write_text(run.stdout, encoding="ascii")
    PROBE_STDERR.write_text(run.stderr, encoding="utf-8")
    validation = parse_output(run.stdout, recursive)

    artifact = {
        "kind": "hardy_block20_t5_component_probe",
        "status": "source_derived_374_call_probe_with_exact_t5_replay",
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
            "The standalone one-CPU probe copies accepted initialization and PSI code, reproduces the source t5 "
            "statement order, and requires exact saved t5 replay at all 374 finite recursive calls. It exposes "
            "finite component values and intrinsic real-erfc arguments, outputs, and weights; it does not itself "
            "prove erfc accuracy, endpoint-saddle accuracy, recurrence accuracy, RH, or a prize-level theorem."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        "built Hardy block-20 t5 component probe: "
        f"374 calls, {validation['intrinsic_erfc_call_count']} intrinsic erfc values, "
        "374 exact t5 bit replays"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
