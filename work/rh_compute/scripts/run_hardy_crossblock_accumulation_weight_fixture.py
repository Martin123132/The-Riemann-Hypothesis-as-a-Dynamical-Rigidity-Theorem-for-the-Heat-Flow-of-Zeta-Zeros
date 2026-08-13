#!/usr/bin/env python3
"""Emit exact binary128 outer accumulation weights for Hardy blocks 20--35."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import run_hardy_block20_accumulation_weight_fixture as block20
import run_hardy_chain_telemetry_fixture as chain


EXTERNAL_ROOT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode"
ACCEPTED_SOURCE = EXTERNAL_ROOT / "resumable/zeta14cubicmult_resumable.f90"
REFERENCE_ROOT = EXTERNAL_ROOT / "fixtures/chain_telemetry/t1e10_crossblock_probe/enabled"
FIXTURE_ROOT = EXTERNAL_ROOT / "fixtures/accumulation_weights/t1e10_crossblock"
RUN_ROOT = FIXTURE_ROOT / "run"
GENERATED_SOURCE = FIXTURE_ROOT / "generated/zeta14cubicmult_crossblock_accumulation_weights.f90"
WEIGHTS = RUN_ROOT / "crossblock_accumulation_weights.txt"
RESULT = FIXTURE_ROOT / "fixture_result.json"
CHECKER = SCRIPT_ROOT / "check_hardy_crossblock_accumulation_weight_fixture.py"
DEFAULT_BUILD_ROOT = Path(
    r"C:\Users\ollet\Documents\Codex\third_party\hardy_fastcode_build\crossblock_accumulation_weights_20260809"
)
BINARY_NAME = "zeta14cubicmult.crossblock_accumulation_weights"
FIRST_BLOCK = 20
LAST_BLOCK = 35
SUM_COUNT = 212
OUTPUT_COUNT = 15
EXPECTED_ROWS = (LAST_BLOCK - FIRST_BLOCK + 1) * SUM_COUNT * OUTPUT_COUNT


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def replace_once(source: str, old: str, new: str, label: str) -> str:
    require(source.count(old) == 1, f"{label} source anchor count drift")
    return source.replace(old, new, 1)


def generate_source() -> dict[str, str]:
    original = ACCEPTED_SOURCE.read_text(encoding="utf-8")
    generated = replace_once(
        original,
        "call alphasum(MTM,rae,pcsum)",
        "call alphasum(iblock,jsum,MTM,rae,pcsum)",
        "alphasum call",
    )
    generated = replace_once(
        generated,
        "subroutine alphasum(MTM,rae,pcsum)",
        "subroutine alphasum(iblock,jsum,MTM,rae,pcsum)",
        "alphasum signature",
    )
    generated = replace_once(
        generated,
        "  integer, parameter :: dp = selected_real_kind(33)  \n"
        "  integer, parameter :: dp1 = selected_int_kind(16) \n\n"
        "  integer       :: K,ip,mc,i\n"
        "  integer(dp1)  :: MTM,numbercalc,nchalf,jj",
        "  integer, parameter :: dp = selected_real_kind(33)  \n"
        "  integer, parameter :: dp1 = selected_int_kind(16) \n"
        "  integer, parameter :: hex_kind = selected_int_kind(32)\n\n"
        "  integer       :: K,ip,mc,i,iblock,rh_weight_unit\n"
        "  integer(dp1)  :: jsum,MTM,numbercalc,nchalf,jj",
        "alphasum declarations",
    )
    anchor = (
        "  do i=1,nchalf\n"
        "   c1=(0.0,1.0)*wp(i)\n"
        "   c1=exp(c1)\n"
        "   c2=(0.0,1.0)*wm(i)\n"
        "   c2=exp(c2)\n"
        "   pcsum(i)=amp(i)*(real(c1*pcgsum(1)+c2*pcgsum(2)))\n"
        "  enddo \n"
        "  return"
    )
    observer = (
        "  do i=1,nchalf\n"
        "   c1=(0.0,1.0)*wp(i)\n"
        "   c1=exp(c1)\n"
        "   c2=(0.0,1.0)*wm(i)\n"
        "   c2=exp(c2)\n"
        "   pcsum(i)=amp(i)*(real(c1*pcgsum(1)+c2*pcgsum(2)))\n"
        "  enddo \n\n"
        "  ! RH_CROSSBLOCK_ACCUMULATION_WEIGHT_OBSERVER_BEGIN\n"
        "  if ((iblock.ge.20).and.(iblock.le.35)) then\n"
        "     open(newunit=rh_weight_unit,file='crossblock_accumulation_weights.txt', &\n"
        "          status='unknown',position='append',action='write')\n"
        "     do i=1,nchalf\n"
        "        c1=exp((0.0,1.0)*wp(i))\n"
        "        c2=exp((0.0,1.0)*wm(i))\n"
        "        write(rh_weight_unit,'(I0,1X,I0,1X,I0,1X,8(Z32.32,1X))') iblock,jsum,i, &\n"
        "             transfer(rae,0_hex_kind),transfer(amp(i),0_hex_kind), &\n"
        "             transfer(wp(i),0_hex_kind),transfer(wm(i),0_hex_kind), &\n"
        "             transfer(real(c1,kind=dp),0_hex_kind), &\n"
        "             transfer(aimag(c1),0_hex_kind), &\n"
        "             transfer(real(c2,kind=dp),0_hex_kind), &\n"
        "             transfer(aimag(c2),0_hex_kind)\n"
        "     enddo\n"
        "     flush(rh_weight_unit)\n"
        "     close(rh_weight_unit)\n"
        "  endif\n"
        "  ! RH_CROSSBLOCK_ACCUMULATION_WEIGHT_OBSERVER_END\n"
        "  return"
    )
    generated = replace_once(generated, anchor, observer, "accumulation observer")
    GENERATED_SOURCE.parent.mkdir(parents=True, exist_ok=True)
    GENERATED_SOURCE.write_text(generated, encoding="utf-8")
    require(ACCEPTED_SOURCE.read_text(encoding="utf-8") == original, "accepted source changed")
    return {
        "accepted_source_sha256": hashlib.sha256(original.encode()).hexdigest(),
        "generated_source_sha256": hashlib.sha256(generated.encode()).hexdigest(),
    }


def configure_shared_runner() -> None:
    block20.REFERENCE_ROOT = REFERENCE_ROOT
    block20.FIXTURE_ROOT = FIXTURE_ROOT
    block20.RUN_ROOT = RUN_ROOT
    block20.GENERATED_SOURCE = GENERATED_SOURCE
    block20.WEIGHTS = WEIGHTS
    block20.RESULT = RESULT
    block20.CHECKER = CHECKER
    block20.BINARY_NAME = BINARY_NAME


def validate_weights() -> dict[str, object]:
    lines = WEIGHTS.read_text(encoding="ascii").splitlines()
    require(len(lines) == EXPECTED_ROWS, "cross-block weight row-count drift")
    pattern = re.compile(r"([0-9]+) ([0-9]+) ([0-9]+)(?: ([0-9A-F]{32})){8} ?")
    roster: list[tuple[int, int, int]] = []
    for line_number, line in enumerate(lines, 1):
        require(pattern.fullmatch(line) is not None, f"invalid weight row {line_number}")
        fields = line.split()
        roster.append((int(fields[0]), int(fields[1]), int(fields[2])))
    expected = [
        (block, sum_index, output_index)
        for block in range(FIRST_BLOCK, LAST_BLOCK + 1)
        for sum_index in range(1, SUM_COUNT + 1)
        for output_index in range(1, OUTPUT_COUNT + 1)
    ]
    require(roster == expected, "cross-block accumulation-weight roster drift")
    return {
        "row_count": len(lines),
        "block_count": LAST_BLOCK - FIRST_BLOCK + 1,
        "first_block": FIRST_BLOCK,
        "last_block": LAST_BLOCK,
        "sum_count_per_block": SUM_COUNT,
        "output_count": OUTPUT_COUNT,
        "roster_exact": True,
        "binary128_payload_count": len(lines) * 8,
        "sha256": file_hash(WEIGHTS),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-root", type=Path, default=DEFAULT_BUILD_ROOT)
    args = parser.parse_args()

    configure_shared_runner()
    source = generate_source()
    image_id = chain.inspect_image_id()
    binary, build_record = block20.build(args.build_root)
    block20.reset_run_root()
    run_record = block20.invoke(binary)
    weights = validate_weights()

    new_journal = chain.load_jsonl(RUN_ROOT / "checkpoint.jsonl", decimal=True)
    reference_journal = chain.load_jsonl(REFERENCE_ROOT / "checkpoint.jsonl", decimal=True)
    require(
        chain.without_run_id(new_journal) == chain.without_run_id(reference_journal),
        "cross-block weight observer changed evaluator checkpoint state",
    )
    new_values = chain.parse_grand_totals(RUN_ROOT / "run.output.txt")
    reference_values = chain.parse_grand_totals(REFERENCE_ROOT / "run.output.txt")
    require(new_values == reference_values, "cross-block weight observer changed displayed Hardy values")

    artifact = {
        "kind": "hardy_crossblock_accumulation_weight_fixture",
        "status": "exact_binary128_blocks_20_35_weight_fixture_with_evaluator_equivalence_validated",
        "source": source,
        "image_id": image_id,
        "build": build_record,
        "run": run_record,
        "weights": weights,
        "equivalence": {
            "checkpoint_records": len(new_journal),
            "checkpoint_state_equal_ignoring_run_id": True,
            "displayed_hardy_values_exact": True,
            "displayed_value_count": len(new_values),
            "reference_checkpoint_sha256": file_hash(REFERENCE_ROOT / "checkpoint.jsonl"),
            "reference_output_sha256": file_hash(REFERENCE_ROOT / "run.output.txt"),
        },
        "artifacts": {
            "generated_source": {"path": relative(GENERATED_SOURCE), "sha256": file_hash(GENERATED_SOURCE)},
            "weights": {"path": relative(WEIGHTS), "sha256": file_hash(WEIGHTS)},
            "checkpoint": {"path": relative(RUN_ROOT / "checkpoint.jsonl"), "sha256": file_hash(RUN_ROOT / "checkpoint.jsonl")},
            "output": {"path": relative(RUN_ROOT / "run.output.txt"), "sha256": file_hash(RUN_ROOT / "run.output.txt")},
            "runner": {"path": relative(Path(__file__)), "sha256": file_hash(Path(__file__))},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "proof_boundary": (
            "Finite source-equivalent observer fixture for exact binary128 outer amplitude and phase "
            "factors in blocks 20--35 at t=10^10. It proves no local kernel error, accumulated "
            "accuracy estimate, height-uniform theorem, outer Hardy remainder, RH, or prize result."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(f"built cross-block accumulation-weight fixture: {EXPECTED_ROWS} rows, evaluator state unchanged")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
