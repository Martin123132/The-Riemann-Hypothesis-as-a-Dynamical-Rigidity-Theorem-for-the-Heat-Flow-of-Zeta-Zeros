#!/usr/bin/env python3
"""Build the exact binary structural atlas for tractable t=1e10 blocks 20-35."""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
import run_hardy_chain_telemetry_fixture as telemetry_base
from flint import acb, ctx


FIXTURE_ROOT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_crossblock_probe"
RUN_ROOT = FIXTURE_ROOT / "enabled"
TELEMETRY = RUN_ROOT / "chain_telemetry.jsonl"
CHECKPOINT = RUN_ROOT / "checkpoint.jsonl"
OUTPUT = RUN_ROOT / "run.output.txt"
INPUT = RUN_ROOT / "inputs3.nml"
REFERENCE_ROOT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10/enabled"
BLOCK20_FIXTURE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/chain_telemetry/t1e10_block20_full/fixture_result.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate.py"
EXPECTED_BLOCKS = tuple(range(20, 36))
EXPECTED_CALLS_PER_BLOCK = 424


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def magnitude_upper(payload: list[str]) -> Fraction:
    return cells.bound_fraction(abs(point_balls.binary128_complex(payload)).upper())


def build() -> dict[str, Any]:
    for path in (TELEMETRY, CHECKPOINT, OUTPUT, INPUT, BLOCK20_FIXTURE, CHECKER):
        require(path.is_file(), f"missing cross-block atlas dependency: {path}")
    ctx.dps = 120
    ctx.threads = 1
    validation = telemetry_base.validate_telemetry(TELEMETRY, 20000)
    require(validation["chain_count"] == len(EXPECTED_BLOCKS) * EXPECTED_CALLS_PER_BLOCK, "cross-block call count drift")
    records = telemetry_base.load_jsonl(TELEMETRY)
    chains = [row for row in records if row.get("type") == "chain"]
    levels: dict[int, dict[int, dict[str, Any]]] = defaultdict(dict)
    recurrences: dict[int, dict[str, Any]] = {}
    for row in records:
        chain = int(row.get("chain", 0))
        if row.get("type") == "level":
            levels[chain][int(row["level"])] = row
        elif row.get("type") == "recurrence":
            recurrences[chain] = row

    reference_journal = telemetry_base.load_jsonl(REFERENCE_ROOT / "checkpoint.jsonl", decimal=True)
    scout_journal = telemetry_base.load_jsonl(CHECKPOINT, decimal=True)
    require(
        telemetry_base.without_run_id(reference_journal) == telemetry_base.without_run_id(scout_journal),
        "cross-block telemetry changed checkpoint state",
    )
    require(
        telemetry_base.parse_grand_totals(REFERENCE_ROOT / "run.output.txt")
        == telemetry_base.parse_grand_totals(OUTPUT),
        "cross-block telemetry changed displayed Hardy values",
    )

    rows: list[dict[str, Any]] = []
    block_accumulator: dict[int, dict[str, Any]] = {
        block: {
            "mit": Counter(),
            "child": Counter(),
            "parent_lengths": set(),
            "min_a1_margin": None,
            "min_frac_margin": None,
            "min_curvature": None,
            "min_one_plus_z": None,
            "max_abs_z": Fraction(0),
            "max_frac_rounding": Fraction(0),
            "max_local_defect": Fraction(0),
            "max_local_defect_chain": None,
        }
        for block in EXPECTED_BLOCKS
    }

    for record in chains:
        chain = int(record["chain"])
        block = int(record["block"])
        require(block in block_accumulator, f"unexpected telemetry block {block}")
        mit = int(record["mit"])
        require(mit in (1, 2), f"chain {chain} unexpected MIT {mit}")
        parent = levels[chain][1]
        require(int(parent["degree"]) == 3, f"chain {chain} noncubic parent")
        a1, a2, a3 = (cells.binary128_fraction(value) for value in parent["coefficients_hex"])
        y = cells.binary128_fraction(parent["xr_hex"])
        require(y == 2 * a2 and y > 0, f"chain {chain} xr coefficient drift")
        require(Fraction(-1, 2) < a1 < Fraction(1, 2) and a1 != 0, f"chain {chain} a1 cell drift")
        require(a3 != 0, f"chain {chain} zero cubic coefficient")
        length = int(parent["length"])
        require(length == int(record["initial_length"]), f"chain {chain} parent length drift")
        curvature_endpoint = y + 6 * a3 * length
        curvature_floor = min(y, curvature_endpoint)
        require(curvature_floor > 0, f"chain {chain} curvature sign drift")
        a1_margin = min(a1 + Fraction(1, 2), Fraction(1, 2) - a1, abs(a1))
        child_length: int | None = None
        frac_margin: Fraction | None = None
        max_abs_z = Fraction(0)
        min_one_plus_z = Fraction(1)
        frac_rounding = Fraction(0)
        local_defect_upper: Fraction | None = None
        if mit == 2:
            child = levels[chain][2]
            child_length = int(child["length"])
            xi = a1 + y * length + 3 * a3 * length**2
            require(xi.numerator // xi.denominator == child_length, f"chain {chain} reconstructed dual selector drift")
            frac = xi - child_length
            require(0 < frac < 1, f"chain {chain} dual fractional boundary drift")
            source_frac = cells.binary128_fraction(parent["frac_length_hex"])
            require(0 < source_frac < 1, f"chain {chain} stored frac boundary drift")
            frac_margin = min(frac, 1 - frac)
            frac_rounding = abs(frac - source_frac)
            z_values = [12 * a3 * (Fraction(index) - a1) / y**2 for index in (0, child_length)]
            max_abs_z = max(abs(value) for value in z_values)
            min_one_plus_z = min(1 + value for value in z_values)
            require(min_one_plus_z > 0, f"chain {chain} cubic radical domain drift")
            recurrence = recurrences[chain]
            require(bool(recurrence["direct_ok"]), f"chain {chain} direct audit missing")
            local_defect_upper = magnitude_upper(recurrence["local_defect_hex"])

        acc = block_accumulator[block]
        acc["mit"][mit] += 1
        if child_length is not None:
            acc["child"][child_length] += 1
        acc["parent_lengths"].add(length)
        acc["min_a1_margin"] = a1_margin if acc["min_a1_margin"] is None else min(acc["min_a1_margin"], a1_margin)
        if frac_margin is not None:
            acc["min_frac_margin"] = frac_margin if acc["min_frac_margin"] is None else min(acc["min_frac_margin"], frac_margin)
            acc["min_one_plus_z"] = min_one_plus_z if acc["min_one_plus_z"] is None else min(acc["min_one_plus_z"], min_one_plus_z)
            acc["max_abs_z"] = max(acc["max_abs_z"], max_abs_z)
            acc["max_frac_rounding"] = max(acc["max_frac_rounding"], frac_rounding)
        acc["min_curvature"] = curvature_floor if acc["min_curvature"] is None else min(acc["min_curvature"], curvature_floor)
        if local_defect_upper is not None and local_defect_upper > acc["max_local_defect"]:
            acc["max_local_defect"] = local_defect_upper
            acc["max_local_defect_chain"] = chain

        rows.append(
            {
                "chain": chain,
                "block": block,
                "sum_index": int(record["sum_index"]),
                "branch": int(record["branch"]),
                "mit": mit,
                "parent_length": length,
                "child_length": child_length,
                "a1_selector_margin": decimal(a1_margin),
                "curvature_floor": decimal(curvature_floor),
                "dual_fractional_margin": None if frac_margin is None else decimal(frac_margin),
                "frac_source_rounding_gap": None if frac_margin is None else decimal(frac_rounding),
                "maximum_child_radical_abs_z": None if frac_margin is None else decimal(max_abs_z),
                "minimum_child_radical_one_plus_z": None if frac_margin is None else decimal(min_one_plus_z),
                "source_local_defect_magnitude_upper": None if local_defect_upper is None else decimal(local_defect_upper),
            }
        )

    block_rows: list[dict[str, Any]] = []
    for block in EXPECTED_BLOCKS:
        acc = block_accumulator[block]
        require(sum(acc["mit"].values()) == EXPECTED_CALLS_PER_BLOCK, f"block {block} call roster drift")
        require(len(acc["parent_lengths"]) == 1, f"block {block} parent length not constant")
        recursive_count = acc["mit"][2]
        block_rows.append(
            {
                "block": block,
                "call_count": EXPECTED_CALLS_PER_BLOCK,
                "recursive_call_count": recursive_count,
                "direct_call_count": acc["mit"][1],
                "parent_length": next(iter(acc["parent_lengths"])),
                "child_length_histogram": {str(key): value for key, value in sorted(acc["child"].items())},
                "minimum_a1_selector_margin": decimal(acc["min_a1_margin"]),
                "minimum_curvature_floor": decimal(acc["min_curvature"]),
                "minimum_dual_fractional_margin": None if recursive_count == 0 else decimal(acc["min_frac_margin"]),
                "maximum_frac_source_rounding_gap": None if recursive_count == 0 else decimal(acc["max_frac_rounding"]),
                "maximum_child_radical_abs_z": None if recursive_count == 0 else decimal(acc["max_abs_z"]),
                "minimum_child_radical_one_plus_z": None if recursive_count == 0 else decimal(acc["min_one_plus_z"]),
                "maximum_source_local_defect_magnitude_upper": None if recursive_count == 0 else decimal(acc["max_local_defect"]),
                "maximum_source_local_defect_witness": acc["max_local_defect_chain"],
            }
        )

    recursive_rows = [row for row in rows if row["mit"] == 2]
    aggregate = {
        "block_count": len(EXPECTED_BLOCKS),
        "first_block": min(EXPECTED_BLOCKS),
        "last_block": max(EXPECTED_BLOCKS),
        "call_count": len(rows),
        "calls_per_block": EXPECTED_CALLS_PER_BLOCK,
        "recursive_call_count": len(recursive_rows),
        "direct_call_count": len(rows) - len(recursive_rows),
        "last_recursive_block": max(row["block"] for row in recursive_rows),
        "first_all_direct_block": min(row["block"] for row in block_rows if row["recursive_call_count"] == 0),
        "maximum_parent_length": max(row["parent_length"] for row in rows),
        "minimum_recursive_a1_selector_margin": min(
            (row["a1_selector_margin"] for row in recursive_rows), key=Fraction
        ),
        "minimum_recursive_dual_fractional_margin": min(
            (row["dual_fractional_margin"] for row in recursive_rows), key=Fraction
        ),
        "minimum_recursive_curvature_floor": min(
            (row["curvature_floor"] for row in recursive_rows), key=Fraction
        ),
        "maximum_recursive_child_radical_abs_z": max(
            (row["maximum_child_radical_abs_z"] for row in recursive_rows), key=Fraction
        ),
        "minimum_recursive_child_radical_one_plus_z": min(
            (row["minimum_child_radical_one_plus_z"] for row in recursive_rows), key=Fraction
        ),
        "maximum_recursive_frac_source_rounding_gap": max(
            (row["frac_source_rounding_gap"] for row in recursive_rows), key=Fraction
        ),
        "maximum_source_local_defect_magnitude_upper": max(
            (row["source_local_defect_magnitude_upper"] for row in recursive_rows), key=Fraction
        ),
        "checkpoint_state_equal_ignoring_run_id": True,
        "displayed_hardy_values_exact": True,
    }
    return {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate",
        "status": "rigorous_t1e10_blocks_20_35_binary_structural_atlas_with_route_transition",
        "scope": "All 6,784 tractable cubic calls in blocks 20 through 35 of the accepted t=1e10 evaluator run.",
        "route_contract": {
            "recursive": "MIT=2 has one q recurrence and a length-one-or-two transformed child.",
            "direct": "MIT=1 evaluates the finite cubic parent kernel directly and has no q recurrence.",
            "structural_checks": "Exact binary coefficients preserve a1, dual-floor, curvature, cubic-sign, and radical-domain selectors.",
        },
        "blocks": block_rows,
        "rows": rows,
        "aggregate": aggregate,
        "telemetry_validation": validation,
        "artifacts": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "checkpoint": {"path": relative(CHECKPOINT), "sha256": file_hash(CHECKPOINT)},
            "output": {"path": relative(OUTPUT), "sha256": file_hash(OUTPUT)},
            "input": {"path": relative(INPUT), "sha256": file_hash(INPUT)},
        },
        "dependencies": {
            "block20_fixture": {"path": relative(BLOCK20_FIXTURE), "sha256": file_hash(BLOCK20_FIXTURE)},
            "reference_checkpoint": {"path": relative(REFERENCE_ROOT / "checkpoint.jsonl"), "sha256": file_hash(REFERENCE_ROOT / "checkpoint.jsonl")},
            "reference_output": {"path": relative(REFERENCE_ROOT / "run.output.txt"), "sha256": file_hash(REFERENCE_ROOT / "run.output.txt")},
            "telemetry_runner": {"path": relative(Path(telemetry_base.__file__).resolve()), "sha256": file_hash(Path(telemetry_base.__file__).resolve())},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Replay the block-20 corrected endpoint theorem on the 1,040 recursive calls in blocks 21-28 and "
            "certify the 5,370 direct kernels separately, then derive a height-uniform route-transition theorem."
        ),
        "proof_boundary": (
            "Rigorous finite exact-binary structural atlas and observer-equivalence result at t=1e10 only. Source local "
            "defects are diagnostics, not corrected bounds. No corrected endpoint theorem is yet claimed for blocks "
            "21-28, no direct-kernel error bound for blocks 21-35, and no height-uniform, outer Hardy, Lambda<=0, "
            "PF-infinity, RH, or prize-level conclusion follows."
        ),
    }


def render_note(artifact: dict[str, Any]) -> str:
    lines = [
        "# Hardy t=1e10 cross-block structural atlas",
        "",
        "Date: 2026-08-09",
        f"Status: {artifact['status']}; finite structural atlas only, not a proof of RH",
        "",
        "## Route transition",
        "",
        "The one-CPU observer replay leaves the checkpoint state and all fifteen",
        "displayed Hardy values exactly unchanged. It captures 424 cubic calls in",
        "each block from 20 through 35.",
        "",
        "```text",
        "block  parent N  recursive  direct  child lengths",
    ]
    for row in artifact["blocks"]:
        lines.append(
            f"{row['block']:>5} {row['parent_length']:>9} {row['recursive_call_count']:>10} "
            f"{row['direct_call_count']:>7}  {json.dumps(row['child_length_histogram'], sort_keys=True)}"
        )
    aggregate = artifact["aggregate"]
    lines.extend(
        [
            "```",
            "",
            f"There are `{aggregate['recursive_call_count']}` recursive and `{aggregate['direct_call_count']}` direct calls. "
            f"The final recursive calls occur in block `{aggregate['last_recursive_block']}`; blocks "
            f"`{aggregate['first_all_direct_block']}` through `{aggregate['last_block']}` are entirely direct.",
            "",
            "## Exact structural margins",
            "",
            "Every recursive exact-binary parent preserves its fundamental a1 cell,",
            "dual floor, positive curvature, and real cubic stationary radical.",
            "",
            "```text",
            f"minimum a1 selector margin       {aggregate['minimum_recursive_a1_selector_margin']}",
            f"minimum dual fractional margin  {aggregate['minimum_recursive_dual_fractional_margin']}",
            f"minimum curvature floor         {aggregate['minimum_recursive_curvature_floor']}",
            f"maximum child |z|               {aggregate['maximum_recursive_child_radical_abs_z']}",
            f"minimum child (1+z)             {aggregate['minimum_recursive_child_radical_one_plus_z']}",
            f"maximum stored-frac rounding    {aggregate['maximum_recursive_frac_source_rounding_gap']}",
            "```",
            "",
            "## Boundary",
            "",
            "This identifies the finite scaling domain and its route transition. It",
            "does not promote the block-20 endpoint estimate to blocks 21-28, certify",
            "the later direct-kernel errors, or prove a height-uniform route theorem,",
            "outer Hardy estimate, `Lambda<=0`, PF-infinity, RH, or a prize result.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    artifact = build()
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built t=1e10 cross-block structural atlas: "
        f"{artifact['aggregate']['call_count']} calls, "
        f"{artifact['aggregate']['recursive_call_count']} recursive, "
        f"blocks {artifact['aggregate']['first_all_direct_block']}-35 all direct"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
