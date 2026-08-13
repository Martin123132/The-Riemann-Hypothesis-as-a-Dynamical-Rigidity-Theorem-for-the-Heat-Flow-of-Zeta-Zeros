#!/usr/bin/env python3
"""Enclose the t5 correction from default-real sqrt(2.0) argument scaling."""

from __future__ import annotations

from collections import Counter
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

import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
import jensen_window_pf_newman_c1_hardy_block20_t5_real_erfc_residual_gate as t5_erfc
import jensen_window_pf_newman_c1_hardy_block20_w1_finite_ip_completion_gate as finite_ip
import jensen_window_pf_newman_c1_hardy_point_ball_defect_gate as point_balls
import jensen_window_pf_newman_c1_hardy_q_selector_cell_gate as cells
from flint import acb, arb, ctx


SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
TELEMETRY = finite_ip.TELEMETRY
PROBE_RESULT = finite_ip.PROBE_RESULT
PROBE_OUTPUT = finite_ip.PROBE_OUTPUT
NATIVE_Q = finite_ip.NATIVE_Q
FINITE_IP = finite_ip.RESULT
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_t5_sqrt2_literal_kind_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_t5_sqrt2_literal_kind_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = REPO_ROOT / "work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_block20_t5_sqrt2_literal_kind_gate.py"
PRECISIONS = (180, 260)
SOURCE_SQRT_TWO = finite_ip.binary32_fraction(finite_ip.SOURCE_SQRT_TWO_BINARY32_BITS)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def fraction_decimal(value: Fraction) -> str:
    return cells.decimal_text(value)


def magnitude_upper(value: arb | acb) -> Fraction:
    return cells.bound_fraction(abs(value).upper())


def magnitude_lower(value: acb) -> Fraction:
    return cells.bound_fraction(abs(value).lower())


def source_locations() -> dict[str, list[int]]:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    matches = [number for number, line in enumerate(lines, 1) if "sqrt(2.0)" in line]
    require(matches == [2439, 2484], f"sqrt(2.0) source locations drift: {matches}")
    return {"default_real_sqrt_two_calls": matches}


def evaluate(
    data: dict[str, Any],
    probe_row: dict[str, Any],
    constants: dict[str, str],
    native_record: dict[str, Any],
    finite_ip_record: dict[str, Any],
    dps: int,
) -> dict[str, Any]:
    ctx.dps = dps
    ctx.threads = 1
    parent = data["levels"][1]
    child = data["levels"][2]
    require(int(parent["length"]) == 104 and int(child["length"]) in (1, 2), "sqrt2 chain shape drift")
    a1, _, a3 = (finite_ip.exact_arb(value) for value in parent["coefficients_hex"])
    xr = finite_ip.exact_arb(parent["xr_hex"])
    sqrt_p = finite_ip.exact_arb(constants["sqrt_pi"])
    source_sqrt_two = cells.fraction_ball(SOURCE_SQRT_TWO)
    exact_sqrt_two = arb(2).sqrt()
    sx = xr.sqrt()

    total_delta = acb(0)
    call_rows: list[dict[str, Any]] = []
    for call in probe_row["calls"]:
        index = int(call["index"])
        con1 = arb(index) - a1
        c = con1 / xr - 3 * a3 * con1**2 / xr**3
        if call["region"] == "lower":
            base_argument = sqrt_p * sx * c
        else:
            base_argument = sqrt_p * sx * (arb(104) - c)
        source_model_argument = base_argument * source_sqrt_two
        source_argument = finite_ip.exact_arb(call["argument_hex"])
        argument_generation_gap = source_model_argument - source_argument
        argument_delta = base_argument * (exact_sqrt_two - source_sqrt_two)
        corrected_argument = source_argument + argument_delta
        erfc_delta = corrected_argument.erfc() - source_argument.erfc()
        weight = finite_ip.exact_acb(call["weight_hex"])
        q_delta = weight * erfc_delta
        total_delta += q_delta
        call_rows.append(
            {
                "slot": call["slot"],
                "region": call["region"],
                "index": index,
                "source_argument": t5_erfc.interval_record(source_argument),
                "source_model_argument": t5_erfc.interval_record(source_model_argument),
                "argument_generation_gap": t5_erfc.interval_record(argument_generation_gap),
                "literal_kind_argument_delta": t5_erfc.interval_record(argument_delta),
                "corrected_argument": t5_erfc.interval_record(corrected_argument),
                "erfc_delta": t5_erfc.interval_record(erfc_delta),
                "weight": point_balls.acb_record(weight, 55),
                "q_delta": point_balls.acb_record(q_delta, 55),
            }
        )

    required = native_q.complex_from_record(native_record["required_native_q_compensation"])
    ip_completion = native_q.complex_from_record(finite_ip_record["omitted_displayed_w1_completion"])
    after_sqrt_two = required - total_delta
    after_stack = after_sqrt_two - ip_completion
    return {
        "total_delta": total_delta,
        "required": required,
        "ip_completion": ip_completion,
        "after_sqrt_two": after_sqrt_two,
        "after_stack": after_stack,
        "calls": call_rows,
    }


def build_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    source_root = artifact["normalization"]["source_sqrt_two"]
    return f"""# Hardy block-20 t5 `sqrt(2.0)` literal-kind gate

Date: 2026-08-06

Status: rigorous finite exact-argument correction for the default-real square-root factor; not a contour theorem

## Source issue

Both intrinsic real-erfc arguments in `t5` contain the unsuffixed expression
`sqrt(2.0)`.  Under the admitted gfortran build, `2.0` is default real and the
square root is therefore first rounded to binary32:

```text
bits       = 0x{source_root['bits']}
exact      = {source_root['numerator']} / {source_root['denominator']}
decimal    = {source_root['decimal']}
```

It is then promoted into the binary128 expression.  This is distinct from the
principal positive mathematical root `sqrt(2)`.  No circle or polygon supplies
this constant: it is the algebraic root of `x^2=2`.  The source's `sqrt(pi)`
factor remains fixed, so this gate isolates only the literal-kind error.

For each of the 913 admitted retained calls, write the source argument as
`x_src` and reconstruct its pre-`sqrt(2)` factor `B` from exact binary128
inputs.  Define

```text
delta x = B * (sqrt(2) - sqrt2_binary32),
x_corr  = x_src + delta x,
delta q = W * (erfc(x_corr) - erfc(x_src)).             (1)
```

Here `W` is the independently emitted exact binary128 source weight.  The
previous intrinsic-erfc gate corrected function evaluation at `x_src`; (1)
corrects the generation of the argument and does not count that residual twice.

## Result

```text
maximum source-model/emitted argument gap       <= {aggregate['maximum_argument_generation_gap_abs_upper']}
maximum literal-kind argument shift             <= {aggregate['maximum_argument_delta_abs_upper']}
maximum erfc-value shift                        <= {aggregate['maximum_erfc_delta_abs_upper']}
maximum one-call q shift                        <= {aggregate['maximum_call_q_delta_abs_upper']}
maximum correlated per-chain q shift            <= {aggregate['maximum_total_q_delta_abs_upper']}
minimum residual after sqrt(2) correction       >= {aggregate['minimum_after_sqrt_two_magnitude_lower']}
minimum residual after sqrt(2)+finite-ip stack  >= {aggregate['minimum_after_stack_magnitude_lower']}
stacked residuals excluding zero                   {aggregate['after_stack_excluding_zero_count']} / 374
```

The stacked residual is `R_q-delta_q_sqrt2-T_ip`, where `T_ip` is the separately
admitted displayed-W1 finite-index completion.  Every operation is repeated at
180 and 260 decimal digits with overlap required.

## Boundary

This closes one source literal-kind normalization on one finite block.  It does
not prove that the displayed W1 approximation equals the exact equation-(69)
integrals, bound higher saddle phase or vertical-leg remainders, reconstruct
the unavailable Maple code, control W2--W4 or outer Hardy errors, establish a
height-uniform recurrence, or prove `Lambda<=0`, PF-infinity, RH, or a
prize-level theorem.
"""


def main() -> int:
    for path in (SOURCE, TELEMETRY, PROBE_RESULT, PROBE_OUTPUT, NATIVE_Q, FINITE_IP, CHECKER):
        require(path.is_file(), f"missing sqrt2 dependency: {path}")
    probe_artifact = json.loads(PROBE_RESULT.read_text(encoding="utf-8"))
    require(probe_artifact["validation"]["saved_t5_bit_match_count"] == 374, "sqrt2 probe not admitted")
    require(probe_artifact["build"]["compiler_flags"] == ["-O3", "-ffree-line-length-none"], "sqrt2 compiler flags drift")
    constants, probe = t5_erfc.load_probe()
    recursive = finite_ip.load_recursive_chains()
    native_artifact = json.loads(NATIVE_Q.read_text(encoding="utf-8"))
    native = {row["chain"]: row for row in native_artifact["rows"]}
    ip_artifact = json.loads(FINITE_IP.read_text(encoding="utf-8"))
    ip_rows = {row["chain"]: row for row in ip_artifact["rows"]}
    require(set(recursive) == set(probe) == set(native) == set(ip_rows), "sqrt2 roster mismatch")

    rows: list[dict[str, Any]] = []
    maxima: dict[str, tuple[Fraction, dict[str, Any]]] = {}
    minima: dict[str, tuple[Fraction, int]] = {}
    region_histogram: Counter[str] = Counter()
    after_sqrt_nonzero = 0
    after_stack_nonzero = 0
    sqrt_improved = 0
    sqrt_worsened = 0
    sqrt_indeterminate = 0
    stack_improved = 0
    stack_worsened = 0
    stack_indeterminate = 0

    def observe_max(name: str, value: Fraction, witness: dict[str, Any]) -> None:
        if name not in maxima or value > maxima[name][0]:
            maxima[name] = (value, witness)

    def observe_min(name: str, value: Fraction, chain: int) -> None:
        if name not in minima or value < minima[name][0]:
            minima[name] = (value, chain)

    for chain in sorted(recursive):
        low = evaluate(recursive[chain], probe[chain], constants, native[chain], ip_rows[chain], PRECISIONS[0])
        high = evaluate(recursive[chain], probe[chain], constants, native[chain], ip_rows[chain], PRECISIONS[1])
        for name in ("total_delta", "required", "ip_completion", "after_sqrt_two", "after_stack"):
            require(low[name].overlaps(high[name]), f"chain {chain} sqrt2 {name} precision nonoverlap")
        require(len(low["calls"]) == len(high["calls"]), f"chain {chain} sqrt2 call roster drift")

        for call in high["calls"]:
            region_histogram[call["region"]] += 1
            generation_gap = native_q.interval_from_record(call["argument_generation_gap"])
            argument_delta = native_q.interval_from_record(call["literal_kind_argument_delta"])
            erfc_delta = native_q.interval_from_record(call["erfc_delta"])
            q_delta = native_q.complex_from_record(call["q_delta"])
            witness = {"chain": chain, "slot": call["slot"], "region": call["region"], "index": call["index"]}
            observe_max("generation_gap", magnitude_upper(generation_gap), witness)
            observe_max("argument_delta", magnitude_upper(argument_delta), witness)
            observe_max("erfc_delta", magnitude_upper(erfc_delta), witness)
            observe_max("call_q_delta", magnitude_upper(q_delta), witness)

        total_upper = magnitude_upper(high["total_delta"])
        observe_max("total_delta", total_upper, {"chain": chain})
        before_abs = abs(high["required"])
        sqrt_abs = abs(high["after_sqrt_two"])
        stack_abs = abs(high["after_stack"])
        before_lower = cells.bound_fraction(before_abs.lower())
        before_upper = cells.bound_fraction(before_abs.upper())
        sqrt_lower = cells.bound_fraction(sqrt_abs.lower())
        sqrt_upper = cells.bound_fraction(sqrt_abs.upper())
        stack_lower = cells.bound_fraction(stack_abs.lower())
        stack_upper = cells.bound_fraction(stack_abs.upper())
        observe_min("sqrt", sqrt_lower, chain)
        observe_min("stack", stack_lower, chain)
        observe_max("sqrt_residual", sqrt_upper, {"chain": chain})
        observe_max("stack_residual", stack_upper, {"chain": chain})
        sqrt_excludes = not high["after_sqrt_two"].contains(0)
        stack_excludes = not high["after_stack"].contains(0)
        after_sqrt_nonzero += sqrt_excludes
        after_stack_nonzero += stack_excludes

        if sqrt_upper < before_lower:
            sqrt_class = "improved"
            sqrt_improved += 1
        elif sqrt_lower > before_upper:
            sqrt_class = "worsened"
            sqrt_worsened += 1
        else:
            sqrt_class = "indeterminate"
            sqrt_indeterminate += 1
        if stack_upper < before_lower:
            stack_class = "improved"
            stack_improved += 1
        elif stack_lower > before_upper:
            stack_class = "worsened"
            stack_worsened += 1
        else:
            stack_class = "indeterminate"
            stack_indeterminate += 1

        rows.append(
            {
                "chain": chain,
                "sum_index": int(recursive[chain]["header"]["sum_index"]),
                "branch": int(recursive[chain]["header"]["branch"]),
                "call_count": len(high["calls"]),
                "total_q_delta_sqrt_two": point_balls.acb_record(high["total_delta"], 55),
                "total_q_delta_sqrt_two_magnitude_upper": fraction_decimal(total_upper),
                "required_native_q_compensation": point_balls.acb_record(high["required"], 55),
                "finite_ip_completion": point_balls.acb_record(high["ip_completion"], 55),
                "residual_after_sqrt_two": point_balls.acb_record(high["after_sqrt_two"], 55),
                "residual_after_sqrt_two_magnitude_lower": fraction_decimal(sqrt_lower),
                "residual_after_sqrt_two_magnitude_upper": fraction_decimal(sqrt_upper),
                "residual_after_sqrt_two_excludes_zero": sqrt_excludes,
                "sqrt_two_magnitude_classification": sqrt_class,
                "residual_after_sqrt_two_and_finite_ip": point_balls.acb_record(high["after_stack"], 55),
                "residual_after_stack_magnitude_lower": fraction_decimal(stack_lower),
                "residual_after_stack_magnitude_upper": fraction_decimal(stack_upper),
                "residual_after_stack_excludes_zero": stack_excludes,
                "stack_magnitude_classification": stack_class,
                "precision_overlap": True,
                "calls": high["calls"],
            }
        )

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_block20_t5_sqrt2_literal_kind_gate",
        "status": "rigorous_913_call_default_real_sqrt_two_argument_correction_enclosed",
        "scope": {
            "block": 20,
            "recursive_call_count": 374,
            "retained_erfc_call_count": sum(len(row["calls"]) for row in rows),
            "precision_ladder_decimal_digits": list(PRECISIONS),
            "precision_overlap_count": 374,
        },
        "normalization": {
            "source_sqrt_two": {
                "bits": f"{finite_ip.SOURCE_SQRT_TWO_BINARY32_BITS:08X}",
                "numerator": SOURCE_SQRT_TWO.numerator,
                "denominator": SOURCE_SQRT_TWO.denominator,
                "decimal": fraction_decimal(SOURCE_SQRT_TWO),
            },
            "corrected_sqrt_two": "Arb principal positive sqrt(2)",
            "pi_boundary": "The exact emitted source sqrt(pi) factor is held fixed; no independent pi is introduced.",
        },
        "identity": {
            "argument_delta": "delta_x=B*(sqrt(2)-sqrt2_binary32)",
            "corrected_argument": "x_corr=x_source+delta_x",
            "q_delta": "delta_q=sum W*(erfc(x_corr)-erfc(x_source))",
            "stacked_residual": "R_stack=R_q-delta_q-T_ip",
        },
        "aggregate": {
            "region_histogram": dict(region_histogram),
            "maximum_argument_generation_gap_abs_upper": fraction_decimal(maxima["generation_gap"][0]),
            "maximum_argument_generation_gap_witness": maxima["generation_gap"][1],
            "maximum_argument_delta_abs_upper": fraction_decimal(maxima["argument_delta"][0]),
            "maximum_argument_delta_witness": maxima["argument_delta"][1],
            "maximum_erfc_delta_abs_upper": fraction_decimal(maxima["erfc_delta"][0]),
            "maximum_erfc_delta_witness": maxima["erfc_delta"][1],
            "maximum_call_q_delta_abs_upper": fraction_decimal(maxima["call_q_delta"][0]),
            "maximum_call_q_delta_witness": maxima["call_q_delta"][1],
            "maximum_total_q_delta_abs_upper": fraction_decimal(maxima["total_delta"][0]),
            "maximum_total_q_delta_witness": maxima["total_delta"][1],
            "minimum_after_sqrt_two_magnitude_lower": fraction_decimal(minima["sqrt"][0]),
            "minimum_after_sqrt_two_witness": minima["sqrt"][1],
            "maximum_after_sqrt_two_magnitude_upper": fraction_decimal(maxima["sqrt_residual"][0]),
            "after_sqrt_two_excluding_zero_count": after_sqrt_nonzero,
            "sqrt_two_improved_call_count": sqrt_improved,
            "sqrt_two_worsened_call_count": sqrt_worsened,
            "sqrt_two_indeterminate_call_count": sqrt_indeterminate,
            "minimum_after_stack_magnitude_lower": fraction_decimal(minima["stack"][0]),
            "minimum_after_stack_witness": minima["stack"][1],
            "maximum_after_stack_magnitude_upper": fraction_decimal(maxima["stack_residual"][0]),
            "after_stack_excluding_zero_count": after_stack_nonzero,
            "stack_improved_call_count": stack_improved,
            "stack_worsened_call_count": stack_worsened,
            "stack_indeterminate_call_count": stack_indeterminate,
        },
        "rows": rows,
        "source": {"path": relative(SOURCE), "sha256": file_hash(SOURCE), "locations": source_locations()},
        "sources": {
            "telemetry": {"path": relative(TELEMETRY), "sha256": file_hash(TELEMETRY)},
            "t5_probe_result": {"path": relative(PROBE_RESULT), "sha256": file_hash(PROBE_RESULT)},
            "t5_probe_output": {"path": relative(PROBE_OUTPUT), "sha256": file_hash(PROBE_OUTPUT)},
            "native_q_target": {"path": relative(NATIVE_Q), "sha256": file_hash(NATIVE_Q)},
            "finite_ip_completion": {"path": relative(FINITE_IP), "sha256": file_hash(FINITE_IP)},
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "python_flint": {"version": point_balls.flint.__version__, "flint_version": point_balls.flint.__FLINT_VERSION__},
        },
        "next_handoff": {
            "analytic_target": "With function, argument-kind, and displayed finite-index corrections separated, derive the exact-contour minus displayed-W1 saddle/vertical-leg remainder.",
            "source_fix": "A future evaluator repair should suffix the literal at source level and revalidate equivalence; this proof artifact does not modify the accepted evaluator.",
        },
        "proof_boundary": (
            "Rigorous finite correction for replacing the two source default-real sqrt(2.0) factors by the principal "
            "mathematical root at 913 admitted calls, plus a finite stack with the displayed-W1 index completion. "
            "It does not bound the exact contour, saddle, vertical-leg, W2-W4, hierarchy, outer Hardy, or uniform-height "
            "remainders and does not establish Lambda <= 0, PF-infinity, RH, or a prize-level theorem."
        ),
    }
    NOTE.write_text(build_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        "built Hardy block-20 t5 sqrt2 literal-kind gate: "
        f"{artifact['scope']['retained_erfc_call_count']} calls, "
        f"{after_stack_nonzero} stacked residuals exclude zero"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
