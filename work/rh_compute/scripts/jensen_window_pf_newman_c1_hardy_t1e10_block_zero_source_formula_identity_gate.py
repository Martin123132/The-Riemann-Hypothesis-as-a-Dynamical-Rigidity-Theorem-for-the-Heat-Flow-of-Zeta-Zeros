#!/usr/bin/env python3
"""Localize the saved block-zero source/formula residual at its actual cutoff."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp

import jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate as partition


CHECKPOINT = partition.CHECKPOINT
PARTITION = partition.RESULT
SOURCE = partition.SOURCE
PAPER = partition.PAPER
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_block_zero_source_formula_identity_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_block_zero_source_formula_identity_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISIONS = (70, 110)
EXPECTED_OUTPUTS = 15
M2 = 159577
N1 = 159777
ALPHA_COUNT = (N1 - M2) // 2 + 1
REQUESTED_TOLERANCE = arb("0.005")
SOURCE_MODEL_GUARD = arb("1e-20")
SHIFT_LITERAL_GUARD = arb("2e-12")
SOURCE_MARKERS = (
    "call start(a,t6,yphase,transit,M2)",
    "do ial=M2,N1,2",
    "call ter(ial,1/aar(i),yphasear(i),tar(i),term)",
    "zsum(i)=zsum(i)+term",
    "zsum(i)=zsum(i)*sqrt(8/aar(i))",
    "v=t*(b*(b-s)+log(s+b))+yphase",
    "term=cos(v)/sqrt(b)",
    "tar(i)=t+(i-numbercalc)*0.01",
    "yphasear(i)=yphase+(i-numbercalc)*0.01",
    "aar(i)=a+4*(i-numbercalc)*0.01/(p*a)",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def audit_source() -> None:
    source = SOURCE.read_text(encoding="utf-8")
    for marker in SOURCE_MARKERS:
        require(marker in source, f"block-zero source marker drift: {marker}")


def exact_algebra_audit() -> dict[str, Any]:
    s, b = sp.symbols("s b")
    relation = b**2 - s**2 + 1
    pc = (s + b) ** 2
    pc_expanded = 2 * s**2 + 2 * s * b - 1
    reciprocal_product = (s + b) * (s - b) - 1
    phase_numerator = sp.expand(2 * pc * (b * (b - s) + 1) - pc - 1)
    return {
        "domain": "a>0, alpha>a, s=alpha/a>1, b=sqrt(s^2-1)>0",
        "relation": "b^2-s^2+1=0",
        "pc_identity": "pc=2s^2+2sb-1=(s+b)^2",
        "pc_polynomial_remainder": str(sp.rem(sp.expand(pc - pc_expanded), relation, b)),
        "reciprocal_identity": "(s+b)(s-b)=1",
        "reciprocal_polynomial_remainder": str(sp.rem(sp.expand(reciprocal_product), relation, b)),
        "phase_identity": "b(b-s)+log(s+b)+1=(log(pc)+1/pc+1)/2",
        "phase_polynomial_remainder_after_log_identity": str(sp.rem(phase_numerator, relation, b)),
        "log_identity": "log(pc)=2log(s+b), since s+b>0 and pc=(s+b)^2",
        "amplitude_identity": "sqrt(8/a)/sqrt(b)=2sqrt(2)/(alpha^2-a^2)^(1/4)",
        "amplitude_fourth_power_relation": "b^2=s^2-1 and alpha=a*s",
        "all_polynomial_remainders_zero": all(
            sp.rem(expr, relation, b) == 0
            for expr in (sp.expand(pc - pc_expanded), sp.expand(reciprocal_product), phase_numerator)
        ),
    }


def central_parameters(precision: int) -> dict[str, arb | int]:
    ctx.dps = precision
    t = arb("1e10")
    pi = arb.pi()
    a = (8 * t / pi).sqrt()
    floor_a = 159576
    require(arb(floor_a) < a and a < arb(floor_a + 1), "central a floor drift")
    quotient = math.floor((1e10 + math.pi / 8) / (2 * math.pi))
    yphase = t + pi / 8 - 2 * pi * quotient
    require(yphase > 0 and yphase < 2 * pi, "central yphase reduction drift")
    t6 = ((pi / 8).sqrt() * a) ** (arb(1) / 3)
    g = (M2 - a) * t6
    gm = (a - (M2 - 2)) * t6
    require(g > arb("3.2"), "upper transition branch activated")
    require(gm > arb("3.2"), "lower transition branch activated")
    return {"t": t, "a": a, "yphase": yphase, "quotient": quotient, "t6": t6, "g": g, "gm": gm}


def default_real_shift(offset: Decimal) -> Decimal:
    shift_index = int(offset * 100)
    require(Decimal(shift_index) / 100 == offset, "shift is not an integer hundredth")
    literal = struct.unpack("f", struct.pack("f", 0.01))[0]
    product = struct.unpack("f", struct.pack("f", shift_index * literal))[0]
    return Decimal.from_float(product)


def ideal_block_zero(offset: Decimal, precision: int, source_literal: bool = True) -> dict[str, arb]:
    ctx.dps = precision
    params = central_parameters(precision)
    pi = arb.pi()
    effective_offset = default_real_shift(offset) if source_literal else offset
    delta = arb(format(effective_offset, "f"))
    t = params["t"] + delta
    a = params["a"] + 4 * delta / (pi * params["a"])
    yphase = params["yphase"] + delta
    require(yphase > 0 and yphase < 2 * pi, "shifted yphase wrap drift")
    raw = arb(0)
    for alpha in range(M2, N1 + 1, 2):
        s = arb(alpha) / a
        b = (s * s - 1).sqrt()
        phase = t * (b * (b - s) + (s + b).log()) + yphase
        raw += phase.cos() / b.sqrt()
    scale = (8 / a).sqrt()
    return {
        "t": t,
        "a": a,
        "yphase": yphase,
        "effective_offset": delta,
        "raw_sum": raw,
        "normalization": scale,
        "normalized_sum": raw * scale,
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    transition = artifact["transition_certificate"]["high_precision"]
    return f"""# Block-zero source-formula identity gate

Date: 2026-08-09

Status: finite source localization validated; not a proof of RH

Source-attribution notice: the `0.0161` quantity below is the residual of the
published/source hybrid with its equation-(124)/(125) endpoint convention.
The later midpoint comparison is an introduced diagnostic alternative, not a
cutoff rule stated in the paper.  It therefore cannot be subtracted here and
called a repair of equations (126)--(127) without an independent theorem.

The source starts block zero at `M2={M2}` and tests two transition collars.
Arb certifies

```text
g  = {transition['g_ball']}
gm = {transition['gm_ball']}
```

Both are strictly above `3.2`, so neither branch executes and the transition
contribution is exactly zero in this run.

Put `s=alpha/a`, `b=sqrt(s^2-1)`, and `pc=(s+b)^2`.  Exact polynomial
reduction modulo `b^2=s^2-1`, together with positivity, proves

```text
sqrt(8/a)/sqrt(b) = 2sqrt(2)/(alpha^2-a^2)^(1/4),
b(b-s)+log(s+b)+1 = (log(pc)+1/pc+1)/2.
```

Thus each normalized source `ter` summand is exactly the left side of paper
equation (62) in ideal arithmetic.  The source normalization and phase are not
independent approximations or mismatched conventions.

An independent {PRECISIONS[0]}/{PRECISIONS[1]}-digit Arb reconstruction of all
{ALPHA_COUNT} odd-alpha summands at each of the fifteen outputs differs from
the saved binary128 block-zero checkpoint by at most

```text
{aggregate['maximum_source_model_residual_absolute_ball']}
```

which is below `1e-20` on every output.  Yet the reconstructed equation-(62)
sum also exposes the source's unsuffixed default-real `0.01` shift literal.
Replacing that rounded literal by the intended exact hundredths changes the
normalized block by at most

```text
{aggregate['maximum_default_real_shift_effect_absolute_ball']}
```

which is below `2e-12` on every output.  The reconstructed equation-(62)
sum is below its exact equation-(124) classical target by more than `0.005` on
all fifteen, with absolute residual between

```text
{aggregate['minimum_model_target_residual_absolute_ball']}
{aggregate['maximum_model_target_residual_absolute_ball']}
```

The actual source discrepancy is therefore not created by the transition,
the displayed amplitude normalization, phase translation, or observable
binary128 departure from their ideal formula.  It is the source-aligned
equation-(62)/(124)--(127) hybrid residual at these saved heights.  This gate
does not decompose that residual into a proved approximation error and a
proved cutoff correction.  It supplies no height-uniform remainder bound,
does not repair the hybrid formula, and has no RH implication.

The paper's equation (57) supplies only the unsigned global estimate
`E(t)<6.15*t^(-1/12)`.  At `t=10^10` its right side is
`{artifact['published_remainder_context']['equation57_bound_at_1e10_ball']}`,
or `{artifact['published_remainder_context']['bound_to_0p005_ratio_ball']}` times
the saved tolerance.  It falls below `0.005` only for
`t>1230^12={artifact['published_remainder_context']['strict_threshold_integer']}`.
That global estimate is not a signed or cellwise block-zero remainder.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    audit_source()
    algebra = exact_algebra_audit()
    require(algebra["all_polynomial_remainders_zero"], "equation-(62) algebra audit failed")
    source_blocks = partition.load_stage_zero()
    source_block_zero = source_blocks[0]
    parent = json.loads(PARTITION.read_text(encoding="utf-8"))
    require(len(parent["output_rows"]) == EXPECTED_OUTPUTS, "partition output roster drift")

    low_transition = central_parameters(PRECISIONS[0])
    high_transition = central_parameters(PRECISIONS[1])
    for key in ("a", "yphase", "t6", "g", "gm"):
        require(low_transition[key].overlaps(high_transition[key]), f"transition precision drift: {key}")

    rows: list[dict[str, Any]] = []
    source_model_abs: list[arb] = []
    shift_effect_abs: list[arb] = []
    model_target_abs: list[arb] = []
    source_target_abs: list[arb] = []
    for output_index in range(1, EXPECTED_OUTPUTS + 1):
        offset = Decimal(output_index - 8) / Decimal(100)
        target_t = format(Decimal("1e10") + offset, "f")
        low_model = ideal_block_zero(offset, PRECISIONS[0], source_literal=True)
        high_model = ideal_block_zero(offset, PRECISIONS[1], source_literal=True)
        low_intended = ideal_block_zero(offset, PRECISIONS[0], source_literal=False)
        high_intended = ideal_block_zero(offset, PRECISIONS[1], source_literal=False)
        for key in ("a", "yphase", "raw_sum", "normalization", "normalized_sum"):
            require(low_model[key].overlaps(high_model[key]), f"model precision drift at output {output_index}: {key}")
            require(low_intended[key].overlaps(high_intended[key]), f"intended-shift model precision drift at output {output_index}: {key}")

        parent_row = parent["output_rows"][output_index - 1]
        cutoff = int(parent_row["block_rows"][0]["cutoff_n"])
        low_target = partition.cumulative_targets(target_t, [cutoff], PRECISIONS[0])[cutoff]
        high_target = partition.cumulative_targets(target_t, [cutoff], PRECISIONS[1])[cutoff]
        require(low_target.overlaps(high_target), f"target precision drift at output {output_index}")

        ctx.dps = PRECISIONS[1]
        source_value = arb(source_block_zero["zsum"][output_index - 1])
        source_minus_model = source_value - high_model["normalized_sum"]
        intended_minus_source_model = high_intended["normalized_sum"] - high_model["normalized_sum"]
        model_minus_target = high_model["normalized_sum"] - high_target
        source_minus_target = source_value - high_target
        decomposition = source_minus_model + model_minus_target - source_minus_target
        require(decomposition.contains(0), f"residual decomposition failed at output {output_index}")
        require(abs(source_minus_model) < SOURCE_MODEL_GUARD, f"source/model guard failed at output {output_index}")
        require(abs(intended_minus_source_model) < SHIFT_LITERAL_GUARD, f"default-real shift guard failed at output {output_index}")
        require(model_minus_target < -REQUESTED_TOLERANCE, f"model target does not reject tolerance at output {output_index}")
        parent_residual = arb(parent_row["block_rows"][0]["cumulative_residual_ball"])
        require(source_minus_target.overlaps(parent_residual), f"parent residual drift at output {output_index}")
        source_model_abs.append(abs(source_minus_model))
        shift_effect_abs.append(abs(intended_minus_source_model))
        model_target_abs.append(abs(model_minus_target))
        source_target_abs.append(abs(source_minus_target))
        rows.append(
            {
                "output_index": output_index,
                "output_label": f"t{offset:+.2f}",
                "target_t": target_t,
                "alpha_count": ALPHA_COUNT,
                "classical_cutoff_n": cutoff,
                "low_precision": {
                    "precision_decimal_digits": PRECISIONS[0],
                    "raw_ter_sum_ball": low_model["raw_sum"].str(PRECISIONS[0], more=True),
                    "normalization_ball": low_model["normalization"].str(PRECISIONS[0], more=True),
                    "normalized_equation62_sum_ball": low_model["normalized_sum"].str(PRECISIONS[0], more=True),
                    "intended_exact_shift_normalized_sum_ball": low_intended["normalized_sum"].str(PRECISIONS[0], more=True),
                    "exact_classical_target_ball": low_target.str(PRECISIONS[0], more=True),
                },
                "high_precision": {
                    "precision_decimal_digits": PRECISIONS[1],
                    "shifted_a_ball": high_model["a"].str(PRECISIONS[1], more=True),
                    "shifted_yphase_ball": high_model["yphase"].str(PRECISIONS[1], more=True),
                    "source_default_real_effective_offset_ball": high_model["effective_offset"].str(PRECISIONS[1], more=True),
                    "raw_ter_sum_ball": high_model["raw_sum"].str(PRECISIONS[1], more=True),
                    "normalization_ball": high_model["normalization"].str(PRECISIONS[1], more=True),
                    "normalized_equation62_sum_ball": high_model["normalized_sum"].str(PRECISIONS[1], more=True),
                    "intended_exact_shift_normalized_sum_ball": high_intended["normalized_sum"].str(PRECISIONS[1], more=True),
                    "exact_classical_target_ball": high_target.str(PRECISIONS[1], more=True),
                },
                "saved_source_block_zero_binary128_decimal": source_block_zero["zsum"][output_index - 1],
                "source_minus_model_ball": source_minus_model.str(PRECISIONS[1], more=True),
                "intended_minus_source_default_real_model_ball": intended_minus_source_model.str(PRECISIONS[1], more=True),
                "model_minus_target_ball": model_minus_target.str(PRECISIONS[1], more=True),
                "source_minus_target_ball": source_minus_target.str(PRECISIONS[1], more=True),
                "residual_decomposition_ball": decomposition.str(PRECISIONS[1], more=True),
                "source_model_below_1e_minus_20": True,
                "default_real_shift_effect_below_2e_minus_12": True,
                "model_target_rejects_0p005": True,
            }
        )

    aggregate = {
        "output_count": EXPECTED_OUTPUTS,
        "alpha_terms_per_output": ALPHA_COUNT,
        "total_reconstructed_alpha_terms": EXPECTED_OUTPUTS * ALPHA_COUNT,
        "transition_zero_count": EXPECTED_OUTPUTS,
        "source_model_below_1e_minus_20_count": sum(row["source_model_below_1e_minus_20"] for row in rows),
        "default_real_shift_effect_below_2e_minus_12_count": sum(row["default_real_shift_effect_below_2e_minus_12"] for row in rows),
        "model_target_rejects_0p005_count": sum(row["model_target_rejects_0p005"] for row in rows),
        "minimum_source_model_residual_absolute_ball": min(source_model_abs, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
        "maximum_source_model_residual_absolute_ball": max(source_model_abs, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
        "maximum_default_real_shift_effect_absolute_ball": max(shift_effect_abs, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
        "minimum_model_target_residual_absolute_ball": min(model_target_abs, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
        "maximum_model_target_residual_absolute_ball": max(model_target_abs, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
        "minimum_source_target_residual_absolute_ball": min(source_target_abs, key=lambda x: x.lower()).str(PRECISIONS[1], more=True),
        "maximum_source_target_residual_absolute_ball": max(source_target_abs, key=lambda x: x.upper()).str(PRECISIONS[1], more=True),
        "maximum_source_model_to_source_target_ratio_ball": max(
            (source_model_abs[i] / source_target_abs[i] for i in range(EXPECTED_OUTPUTS)), key=lambda x: x.upper()
        ).str(PRECISIONS[1], more=True),
    }
    require(aggregate["source_model_below_1e_minus_20_count"] == EXPECTED_OUTPUTS, "source/model aggregate drift")
    require(aggregate["default_real_shift_effect_below_2e_minus_12_count"] == EXPECTED_OUTPUTS, "default-real shift aggregate drift")
    require(aggregate["model_target_rejects_0p005_count"] == EXPECTED_OUTPUTS, "model/target aggregate drift")

    ctx.dps = PRECISIONS[1]
    equation57_bound = arb("6.15") * arb("1e10") ** (-arb(1) / 12)
    bound_ratio = equation57_bound / REQUESTED_TOLERANCE
    strict_threshold = 1230**12
    require(equation57_bound > REQUESTED_TOLERANCE, "published equation-(57) bound unexpectedly closes tolerance")

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_block_zero_source_formula_identity_gate",
        "status": "finite_block_zero_equation62_approximation_error_localized",
        "passed": True,
        "scope": {
            "height_center": "1e10",
            "outputs": EXPECTED_OUTPUTS,
            "alpha_range": f"{M2}-{N1} odd",
            "alpha_terms_per_output": ALPHA_COUNT,
            "precisions_decimal_digits": list(PRECISIONS),
        },
        "transition_certificate": {
            "source_threshold": "3.2",
            "decision": "both branches false; transition=0",
            "low_precision": {
                "precision_decimal_digits": PRECISIONS[0],
                "g_ball": low_transition["g"].str(PRECISIONS[0], more=True),
                "gm_ball": low_transition["gm"].str(PRECISIONS[0], more=True),
            },
            "high_precision": {
                "precision_decimal_digits": PRECISIONS[1],
                "a_ball": high_transition["a"].str(PRECISIONS[1], more=True),
                "t6_ball": high_transition["t6"].str(PRECISIONS[1], more=True),
                "g_ball": high_transition["g"].str(PRECISIONS[1], more=True),
                "gm_ball": high_transition["gm"].str(PRECISIONS[1], more=True),
                "yphase_ball": high_transition["yphase"].str(PRECISIONS[1], more=True),
                "yphase_reduction_quotient": high_transition["quotient"],
            },
        },
        "equation62_source_identity": algebra,
        "output_rows": rows,
        "aggregate": aggregate,
        "published_remainder_context": {
            "paper_equation": "57",
            "bound": "E(t)<6.15*t^(-1/12)",
            "equation57_bound_at_1e10_ball": equation57_bound.str(PRECISIONS[1], more=True),
            "bound_to_0p005_ratio_ball": bound_ratio.str(PRECISIONS[1], more=True),
            "strict_threshold_condition": "t>(6.15/0.005)^12=1230^12",
            "strict_threshold_integer": str(strict_threshold),
            "is_signed_block_zero_remainder": False,
            "closes_0p005_at_1e10": False,
        },
        "decision": {
            "transition_creates_block_zero_defect": False,
            "amplitude_normalization_mismatch_creates_block_zero_defect": False,
            "phase_translation_mismatch_creates_block_zero_defect": False,
            "observable_source_ideal_departure_creates_block_zero_defect": False,
            "default_real_shift_literal_mismatch_creates_block_zero_defect": False,
            "source_aligned_equation62_hybrid_residual_rejects_0p005_at_saved_outputs": True,
            "diagnostic_midpoint_is_a_published_cutoff_rule": False,
            "height_uniform_equation62_remainder_bound_proved": False,
            "finite_gate_has_rh_implication": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "dependencies": {
            "source_checkpoint": {"path": relative(CHECKPOINT), "sha256": file_hash(CHECKPOINT)},
            "equation124_partition": {"path": relative(PARTITION), "sha256": file_hash(PARTITION)},
            "source_file": {"path": relative(SOURCE), "sha256": file_hash(SOURCE)},
            "paper": {"path": relative(PAPER), "sha256": file_hash(PAPER)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": "Derive a signed, height-uniform error bound for the source-aligned equations-(126)--(127) hybrid under its equation-(124)/(125) endpoint convention. Any alternative midpoint cutoff must first be derived as a separate valid hybrid theorem.",
        "proof_boundary": "Rigorous finite source/formula localization at fifteen saved heights only. It certifies the source-aligned hybrid residual but does not assign it to a proved cutoff defect, prove a height-uniform remainder, repair the hybrid representation, prove Lambda<=0 or PF-infinity, prove RH, or reach a prize-level conclusion.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "validated block-zero source identity: "
        f"{aggregate['total_reconstructed_alpha_terms']} terms, transition-zero={aggregate['transition_zero_count']}, "
        f"source-model<1e-20={aggregate['source_model_below_1e_minus_20_count']}, "
        f"eq62-rejects-0.005={aggregate['model_target_rejects_0p005_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
