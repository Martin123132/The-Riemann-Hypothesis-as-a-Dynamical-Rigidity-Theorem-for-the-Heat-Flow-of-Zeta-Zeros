#!/usr/bin/env python3
"""Transport all t=1e10 recursive endpoint majorants through outer weights."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, ctx


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate.py"
WEIGHT_FIXTURE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/accumulation_weights/t1e10_crossblock/fixture_result.json"
WEIGHTS = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/accumulation_weights/t1e10_crossblock/run/crossblock_accumulation_weights.txt"
ENDPOINTS = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_retained_endpoint_gate.json"
BLOCK20 = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_retained_decay_endpoint_majorant_gate.json"
STRUCTURE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_structural_atlas_gate.json"
REQUESTED_ERROR_SCALE = Fraction(5, 1000)
PRECISION_DIGITS = 90


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def upper(value: Any) -> Fraction:
    return cells.bound_fraction(value.upper())


def exact_real(payload: str):
    return cells.fraction_ball(cells.binary128_fraction(payload))


def output_label(index: int) -> str:
    return f"t{(index - 8) / 100:+.2f}"


def load_weight_factors() -> dict[tuple[int, int, int, int], Fraction]:
    factors: dict[tuple[int, int, int, int], Fraction] = {}
    lines = WEIGHTS.read_text(encoding="ascii").splitlines()
    require(len(lines) == 50880, "cross-block weight row-count drift")
    for line_number, line in enumerate(lines, 1):
        fields = line.split()
        require(len(fields) == 11, f"weight field-count drift at line {line_number}")
        block, sum_index, output_index = map(int, fields[:3])
        amplitude = cells.binary128_fraction(fields[4])
        require(amplitude > 0, f"nonpositive amplitude at line {line_number}")
        for branch, offset in ((1, 7), (2, 9)):
            phase = acb(exact_real(fields[offset]), exact_real(fields[offset + 1]))
            key = (block, sum_index, output_index, branch)
            require(key not in factors, f"duplicate weight key {key}")
            factors[key] = amplitude * upper(abs(phase))
    require(len(factors) == 16 * 212 * 15 * 2, "expanded cross-block weight roster drift")
    return factors


def compute() -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    ctx.dps = PRECISION_DIGITS
    ctx.threads = 1
    factors = load_weight_factors()
    endpoint_artifact = json.loads(ENDPOINTS.read_text(encoding="utf-8"))
    block20_artifact = json.loads(BLOCK20.read_text(encoding="utf-8"))
    endpoint_rows = endpoint_artifact["rows"]
    require(len(endpoint_rows) == 1414, "recursive endpoint roster drift")

    by_block_output: dict[tuple[int, int], Fraction] = {
        (block, output): Fraction(0)
        for block in range(20, 29)
        for output in range(1, 16)
    }
    call_counts = {block: 0 for block in range(20, 29)}
    for row in endpoint_rows:
        block = int(row["block"])
        sum_index = int(row["sum_index"])
        branch = int(row["branch"])
        local_bound = Fraction(row["complete_majorant_upper"])
        call_counts[block] += 1
        for output in range(1, 16):
            by_block_output[(block, output)] += (
                factors[(block, sum_index, output, branch)] * local_bound
            )

    block20_expected = {
        int(row["output_index"]): row["retained_decay_transported_majorant_upper"]
        for row in block20_artifact["transported_outputs"]
    }
    require(len(block20_expected) == 15, "block-20 output anchor drift")
    for output in range(1, 16):
        require(
            decimal(by_block_output[(20, output)]) == block20_expected[output],
            "block-20 transported output mismatch at "
            f"{output}: got {decimal(by_block_output[(20, output)])}, "
            f"expected {block20_expected[output]}",
        )

    block_rows: list[dict[str, Any]] = []
    for block in range(20, 29):
        values = [(by_block_output[(block, output)], output) for output in range(1, 16)]
        maximum = max(values)
        minimum = min(values)
        block_rows.append(
            {
                "block": block,
                "recursive_call_count": call_counts[block],
                "minimum_transported_majorant_upper": decimal(minimum[0]),
                "minimum_witness_output": minimum[1],
                "maximum_transported_majorant_upper": decimal(maximum[0]),
                "maximum_witness_output": maximum[1],
            }
        )

    output_rows: list[dict[str, Any]] = []
    for output in range(1, 16):
        total = sum((by_block_output[(block, output)] for block in range(20, 29)), Fraction(0))
        output_rows.append(
            {
                "output_index": output,
                "output_label": output_label(output),
                "block20_majorant_upper": decimal(by_block_output[(20, output)]),
                "later_blocks_majorant_upper": decimal(total - by_block_output[(20, output)]),
                "all_recursive_blocks_majorant_upper": decimal(total),
                "to_requested_scale_ratio": decimal(total / REQUESTED_ERROR_SCALE),
                "within_requested_scale": total < REQUESTED_ERROR_SCALE,
            }
        )

    max_total = max(
        (Fraction(row["all_recursive_blocks_majorant_upper"]), int(row["output_index"]))
        for row in output_rows
    )
    max_later = max(
        (Fraction(row["later_blocks_majorant_upper"]), int(row["output_index"]))
        for row in output_rows
    )
    max_ratio = max(
        (Fraction(row["to_requested_scale_ratio"]), int(row["output_index"]))
        for row in output_rows
    )
    aggregate = {
        "recursive_call_count": len(endpoint_rows),
        "block_count": 9,
        "output_count": 15,
        "requested_error_scale": decimal(REQUESTED_ERROR_SCALE),
        "block20_output_anchor_exact": True,
        "maximum_later_blocks_majorant_upper": decimal(max_later[0]),
        "maximum_later_blocks_witness_output": max_later[1],
        "maximum_all_recursive_blocks_majorant_upper": decimal(max_total[0]),
        "maximum_all_recursive_blocks_witness_output": max_total[1],
        "maximum_to_requested_scale_ratio": decimal(max_ratio[0]),
        "maximum_to_requested_scale_ratio_witness_output": max_ratio[1],
        "outputs_within_requested_scale": sum(bool(row["within_requested_scale"]) for row in output_rows),
    }
    return block_rows, output_rows, aggregate


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    block_table = "\n".join(
        f"{row['block']:>5} {row['recursive_call_count']:>6}  {row['maximum_transported_majorant_upper']}"
        for row in artifact["block_rows"]
    )
    output_table = "\n".join(
        f"{row['output_label']:>7}  {row['block20_majorant_upper']}  "
        f"{row['later_blocks_majorant_upper']}  {row['all_recursive_blocks_majorant_upper']}"
        for row in artifact["output_rows"]
    )
    return f"""# Hardy t=1e10 cross-block recursive output budget gate

Date: 2026-08-09
Status: rigorous finite recursive-call outer-weight transport; not a proof of the complete evaluator or RH

The exact binary128 outer weights are composed without signed cancellation
with every retained endpoint majorant in blocks 20--28.  The block-20 subtotal
reproduces all fifteen previously admitted transported values exactly.

```text
block  calls  maximum over outputs
{block_table}
```

```text
 output  block 20  blocks 21--28  all recursive blocks
{output_table}
```

The largest blocks-21--28 contribution is
`{aggregate['maximum_later_blocks_majorant_upper']}`.  The largest complete
recursive contribution is `{aggregate['maximum_all_recursive_blocks_majorant_upper']}`,
or `{aggregate['maximum_to_requested_scale_ratio']}` times the nominal 0.005
scale.  Exactly `{aggregate['outputs_within_requested_scale']}/15` outputs
close at that scale.

This is a finite exact-point recursive endpoint budget.  It does not include
coefficient-cell inflation outside block 20, any of the 5370 direct-kernel
errors, source arithmetic outside the admitted local model, height uniformity,
the outer Hardy remainder, Lambda<=0, PF-infinity, RH, or a prize theorem.
"""


def main() -> int:
    dependencies = (WEIGHT_FIXTURE, WEIGHTS, ENDPOINTS, BLOCK20, STRUCTURE, CHECKER)
    for path in dependencies:
        require(path.is_file(), f"missing recursive-output dependency: {path}")
    fixture = json.loads(WEIGHT_FIXTURE.read_text(encoding="utf-8"))
    require(fixture["weights"]["row_count"] == 50880, "weight fixture row-count drift")
    require(fixture["equivalence"]["displayed_hardy_values_exact"], "weight fixture lost equivalence")
    block_rows, output_rows, aggregate = compute()
    closes = aggregate["outputs_within_requested_scale"] == 15
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate",
        "status": (
            "rigorous_t1e10_all_recursive_endpoint_output_budget_closes_nominal_scale"
            if closes
            else "rigorous_t1e10_all_recursive_endpoint_output_budget_exceeds_nominal_scale"
        ),
        "scope": {
            "first_block": 20,
            "last_recursive_block": 28,
            "recursive_calls": 1414,
            "outer_weight_rows": 50880,
            "output_count": 15,
            "arb_precision_decimal_digits": PRECISION_DIGITS,
        },
        "block_rows": block_rows,
        "output_rows": output_rows,
        "aggregate": aggregate,
        "dependencies": {
            path.stem: {"path": relative(path), "sha256": file_hash(path)}
            for path in dependencies[:-1]
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": (
            "Certify and transport all 5370 direct-kernel errors, then combine them with this recursive "
            "column. Separately extend coefficient-cell and source-arithmetic transport beyond block 20."
        ),
        "proof_boundary": (
            "Rigorous finite exact-point absolute transport of recursive endpoint majorants at t=10^10 only. "
            "No direct-kernel, later coefficient-cell, height-uniform, outer Hardy, Lambda <= 0, "
            "PF-infinity, RH, or prize-level conclusion follows."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    print(
        "built t=1e10 cross-block recursive output budget: "
        f"max {aggregate['maximum_all_recursive_blocks_majorant_upper']}, "
        f"{aggregate['outputs_within_requested_scale']}/15 below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
