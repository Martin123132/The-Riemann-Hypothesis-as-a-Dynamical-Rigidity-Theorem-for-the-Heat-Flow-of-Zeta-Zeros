#!/usr/bin/env python3
"""Close the t=1e10 endpoint/direct column with exact block-20 corrections."""

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

import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
import jensen_window_pf_newman_c1_hardy_t1e10_crossblock_recursive_output_budget_gate as transport
from flint import ctx


ENDPOINTS = transport.ENDPOINTS
WEIGHTS = transport.WEIGHTS
WEIGHT_FIXTURE = transport.WEIGHT_FIXTURE
LOG_FULL = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_logarithmic_endpoint_ray_arb_full_gate.json"
SOURCE_BRIDGE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_source_to_corrected_model_bridge_gate.json"
RECURSIVE_OUTPUT = transport.RESULT
DIRECT_OUTPUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_direct_kernel_output_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_hybrid_exact_output_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_crossblock_hybrid_exact_output_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = SCRIPT_ROOT / "check_jensen_window_pf_newman_c1_hardy_t1e10_crossblock_hybrid_exact_output_gate.py"
REQUESTED_ERROR_SCALE = Fraction(5, 1000)
SERIALIZATION_PAD = Fraction(1, 10**50)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def upper_abs_record(record: dict[str, Any]) -> Fraction:
    value = native_q.complex_from_record(record)
    return cells.bound_fraction(abs(value).upper())


def compute() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    ctx.dps = transport.PRECISION_DIGITS
    ctx.threads = 1
    factors = transport.load_weight_factors()
    endpoints = json.loads(ENDPOINTS.read_text(encoding="utf-8"))
    log_full = json.loads(LOG_FULL.read_text(encoding="utf-8"))
    bridge = json.loads(SOURCE_BRIDGE.read_text(encoding="utf-8"))
    direct = json.loads(DIRECT_OUTPUT.read_text(encoding="utf-8"))
    prior = json.loads(RECURSIVE_OUTPUT.read_text(encoding="utf-8"))

    require(log_full["aggregate"]["call_count"] == 374, "full-contour call count drift")
    require(log_full["aggregate"]["formula_target_identity_count"] == 374, "formula identity drift")
    require(log_full["aggregate"]["correction_residual_identity_count"] == 374, "correction identity drift")
    require(bridge["aggregate"]["identity_handoff_count"] == 374, "source bridge identity drift")
    require(bridge["aggregate"]["endpoint_decomposition_overlap_count"] == 374, "source bridge overlap drift")
    require(len(endpoints["rows"]) == 1414, "endpoint roster drift")

    exact_block20 = {
        int(row["chain"]): upper_abs_record(row["exact_minus_paper"])
        for row in log_full["rows"]
    }
    require(len(exact_block20) == 374, "exact block-20 correction roster drift")
    direct_by_output = {
        int(row["output_index"]): Fraction(row["direct_kernel_majorant_upper"]) + SERIALIZATION_PAD
        for row in direct["output_rows"]
    }
    prior_by_output = {
        int(row["output_index"]): Fraction(row["all_recursive_blocks_majorant_upper"]) + SERIALIZATION_PAD
        for row in prior["output_rows"]
    }

    output_rows: list[dict[str, Any]] = []
    for output in range(1, 16):
        block20_total = Fraction(0)
        later_total = Fraction(0)
        for row in endpoints["rows"]:
            block = int(row["block"])
            factor = factors[(block, int(row["sum_index"]), output, int(row["branch"]))]
            if block == 20:
                local = exact_block20[int(row["chain"])]
                block20_total += factor * local
            else:
                local = Fraction(row["complete_majorant_upper"]) + SERIALIZATION_PAD
                later_total += factor * local
        direct_total = direct_by_output[output]
        combined = block20_total + later_total + direct_total
        prior_total = prior_by_output[output] + direct_total
        output_rows.append(
            {
                "output_index": output,
                "output_label": transport.output_label(output),
                "block20_exact_correction_majorant_upper": decimal(block20_total),
                "blocks21_28_retained_majorant_upper": decimal(later_total),
                "direct_kernel_majorant_upper": decimal(direct_total),
                "hybrid_complete_majorant_upper": decimal(combined),
                "prior_complete_majorant_upper": decimal(prior_total),
                "absolute_improvement": decimal(prior_total - combined),
                "improvement_factor": decimal(prior_total / combined),
                "margin_below_requested_scale": decimal(REQUESTED_ERROR_SCALE - combined),
                "to_requested_scale_ratio": decimal(combined / REQUESTED_ERROR_SCALE),
                "within_requested_scale": combined < REQUESTED_ERROR_SCALE,
            }
        )

    def maximum(field: str) -> tuple[Fraction, int]:
        return max((Fraction(row[field]), int(row["output_index"])) for row in output_rows)

    def minimum(field: str) -> tuple[Fraction, int]:
        return min((Fraction(row[field]), int(row["output_index"])) for row in output_rows)

    max_block20 = maximum("block20_exact_correction_majorant_upper")
    max_later = maximum("blocks21_28_retained_majorant_upper")
    max_direct = maximum("direct_kernel_majorant_upper")
    max_hybrid = maximum("hybrid_complete_majorant_upper")
    max_ratio = maximum("to_requested_scale_ratio")
    min_margin = minimum("margin_below_requested_scale")
    min_improvement = minimum("improvement_factor")
    aggregate = {
        "recursive_call_count": 1414,
        "direct_call_count": 5370,
        "total_call_count": 6784,
        "exact_block20_call_count": 374,
        "retained_later_call_count": 1040,
        "output_count": 15,
        "requested_error_scale": decimal(REQUESTED_ERROR_SCALE),
        "maximum_block20_exact_correction_majorant_upper": decimal(max_block20[0]),
        "maximum_block20_exact_correction_witness_output": max_block20[1],
        "maximum_blocks21_28_retained_majorant_upper": decimal(max_later[0]),
        "maximum_blocks21_28_retained_witness_output": max_later[1],
        "maximum_direct_kernel_majorant_upper": decimal(max_direct[0]),
        "maximum_direct_kernel_witness_output": max_direct[1],
        "maximum_hybrid_complete_majorant_upper": decimal(max_hybrid[0]),
        "maximum_hybrid_complete_witness_output": max_hybrid[1],
        "maximum_to_requested_scale_ratio": decimal(max_ratio[0]),
        "maximum_to_requested_scale_ratio_witness_output": max_ratio[1],
        "minimum_margin_below_requested_scale": decimal(min_margin[0]),
        "minimum_margin_witness_output": min_margin[1],
        "minimum_improvement_factor": decimal(min_improvement[0]),
        "minimum_improvement_factor_witness_output": min_improvement[1],
        "outputs_within_requested_scale": sum(bool(row["within_requested_scale"]) for row in output_rows),
        "uses_signed_cross_call_cancellation": False,
        "block20_complex_intervals_decoded_from_dyadic_endpoints": True,
    }
    return output_rows, aggregate


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    rows = "\n".join(
        f"{row['output_label']:>7}  {row['block20_exact_correction_majorant_upper']}  "
        f"{row['blocks21_28_retained_majorant_upper']}  {row['hybrid_complete_majorant_upper']}"
        for row in artifact["output_rows"]
    )
    return f"""# Hardy t=1e10 cross-block hybrid exact output gate

Date: 2026-08-09
Status: rigorous finite corrected-model endpoint/direct column below 0.005; not a proof of the complete evaluator, a height-uniform theorem, or RH

The 374 block-20 retained bounds are replaced by the already certified
full-contour intervals for `Q_exact-Q_paper`.  Their dyadic real and imaginary
endpoints are decoded directly before taking magnitudes.  Blocks 21--28 retain
their selector-gap-uniform local majorants, and all 5370 direct kernels retain
their independent certified column.  No signed cancellation between calls is
used.

```text
 output  exact block 20  retained blocks 21--28  complete 6784-call bound
{rows}
```

The maximum exact block-20 column is
`{aggregate['maximum_block20_exact_correction_majorant_upper']}` and the
maximum later recursive column is
`{aggregate['maximum_blocks21_28_retained_majorant_upper']}`.  Including every
direct kernel gives `{aggregate['maximum_hybrid_complete_majorant_upper']}`,
only `{aggregate['maximum_to_requested_scale_ratio']}` of the nominal 0.005
scale.  All `{aggregate['outputs_within_requested_scale']}/15` outputs close,
with minimum remaining margin
`{aggregate['minimum_margin_below_requested_scale']}`.

This is a finite exact-input corrected-model endpoint/direct theorem at
t=10^10.  It does not extend the full block-20 correction to coefficient
cells or later blocks, certify later source arithmetic, control the outer
Hardy remainder, establish height uniformity, Lambda<=0, PF-infinity, RH, or
a prize theorem.
"""


def main() -> int:
    dependencies = (
        ENDPOINTS,
        WEIGHTS,
        WEIGHT_FIXTURE,
        LOG_FULL,
        SOURCE_BRIDGE,
        RECURSIVE_OUTPUT,
        DIRECT_OUTPUT,
        CHECKER,
    )
    for path in dependencies:
        require(path.is_file(), f"missing hybrid output dependency: {path}")
    output_rows, aggregate = compute()
    require(aggregate["outputs_within_requested_scale"] == 15, "hybrid exact-point budget did not close")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_crossblock_hybrid_exact_output_gate",
        "status": "rigorous_t1e10_complete_6784_call_endpoint_direct_exact_point_column_below_0_005",
        "scope": {
            "first_block": 20,
            "last_block": 35,
            "exact_full_contour_block": 20,
            "retained_majorant_blocks": [21, 22, 23, 24, 25, 26, 27, 28],
            "all_direct_blocks": list(range(20, 36)),
            "serialization_pad_per_imported_decimal": [SERIALIZATION_PAD.numerator, SERIALIZATION_PAD.denominator],
        },
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
            "Extend the exact-correction or coefficient-cell theorem through blocks 21--28 and certify "
            "the remaining source arithmetic before attacking a height-uniform block transition and outer Hardy remainder."
        ),
        "proof_boundary": (
            "Rigorous finite exact-point corrected-model endpoint/direct column at t=10^10 only. It does not prove "
            "coefficient-cell or height-uniform control, the outer Hardy representation, Lambda <= 0, "
            "PF-infinity, RH, or a prize-level conclusion."
        ),
    }
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built t=1e10 cross-block hybrid exact output gate: "
        f"6784 calls, max {aggregate['maximum_hybrid_complete_majorant_upper']}, 15/15 below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
