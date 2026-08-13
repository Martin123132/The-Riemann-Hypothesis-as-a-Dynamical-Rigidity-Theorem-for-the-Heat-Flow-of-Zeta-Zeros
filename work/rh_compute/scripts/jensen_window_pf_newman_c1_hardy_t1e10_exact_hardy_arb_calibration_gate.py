#!/usr/bin/env python3
"""Rigorous Arb calibration of the fifteen saved source outputs against Hardy Z."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


CHECKPOINT = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/fixtures/accumulation_weights/t1e10_crossblock/run/checkpoint.jsonl"
OUTER_LEDGER = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_outer_hardy_representation_dependency_ledger_gate.json"
INNER_BOUND = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_rigorous_arb_source_q_dependency_cut_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_exact_hardy_arb_calibration_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_exact_hardy_arb_calibration_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISIONS = (70, 110)
REQUESTED_TOLERANCE = arb("0.005")
EXPECTED_OUTPUTS = 15


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


def load_last_checkpoint() -> dict[str, Any]:
    last = ""
    with CHECKPOINT.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                last = line
    require(bool(last), "checkpoint is empty")
    return json.loads(last, parse_float=str)


def hardy_ball(target_t: str, precision: int) -> dict[str, Any]:
    ctx.dps = precision
    ctx.threads = 1
    t = arb(target_t)
    s = acb(arb("0.5"), t)
    theta = acb(arb("0.25"), t / 2).lgamma().imag - t * arb.pi().log() / 2
    value = acb(0, theta).exp() * s.zeta()
    require(value.imag.contains(0), f"Hardy imaginary enclosure excludes zero at t={target_t}")
    return {
        "precision_decimal_digits": precision,
        "theta_ball": theta.str(precision, more=True),
        "hardy_real_ball": value.real.str(precision, more=True),
        "hardy_imag_ball": value.imag.str(precision, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Exact Hardy-Z Arb calibration at t=10^10

Date: 2026-08-09

Status: finite interval certificate validated; not a proof of RH

The fifteen preserved binary128 source outputs at
`t=10^10-0.07,...,10^10+0.07` were compared directly with

```text
Z(t)=exp(i theta(t)) zeta(1/2+i t),
theta(t)=Im log Gamma(1/4+i t/2)-t log(pi)/2.
```

Arb evaluated each target independently at {PRECISIONS[0]} and
{PRECISIONS[1]} decimal digits.  Every pair of Hardy enclosures overlaps,
every imaginary enclosure contains zero, and every signed
`source-Hardy` interval is strictly negative.

All {EXPECTED_OUTPUTS} source outputs are rigorously farther than `0.005`
from the exact Hardy target.  The smallest certified absolute gap occurs at
output {aggregate['minimum_absolute_error_witness_output']} and is enclosed by

```text
{aggregate['minimum_absolute_error_ball']}
```

The largest occurs at output
{aggregate['maximum_absolute_error_witness_output']} and is enclosed by

```text
{aggregate['maximum_absolute_error_ball']}
```

This falsifies the absolute `0.005` accuracy target for this saved legacy
run.  It does not falsify RH, and it does not prove that every implementation
or every parameter choice in the paper has the same error.

The existing blocks-20--35 internal corrected-model column is too small to
explain the discrepancy by itself: even an adjustment with its most helpful
possible sign leaves a reverse-triangle residual greater than `0.005` on all
fifteen outputs.  At least one still-open outer column from Section 11.295 is
therefore essential.  The next proof-facing task is to replace the non-exact
hybrid outer formula by an exact identity with explicit remainders; fitting or
retuning the source output would not close that theorem obligation.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    checkpoint = load_last_checkpoint()
    outer = json.loads(OUTER_LEDGER.read_text(encoding="utf-8"))
    inner = json.loads(INNER_BOUND.read_text(encoding="utf-8"))
    require(outer["aggregate"]["outputs_hardy_z_certified"] == 0, "outer ledger was overpromoted")
    require(int(checkpoint["stage"]) == 2 and int(checkpoint["nchalf"]) == EXPECTED_OUTPUTS, "checkpoint roster drift")
    require(len(checkpoint["zsum"]) == EXPECTED_OUTPUTS, "source output roster drift")
    require(len(inner["output_rows"]) == EXPECTED_OUTPUTS, "internal-bound roster drift")

    base_t = Decimal(checkpoint["t"])
    rows: list[dict[str, Any]] = []
    for index, (source_text, inner_row) in enumerate(zip(checkpoint["zsum"], inner["output_rows"]), start=1):
        offset = Decimal(index - 8) / Decimal(100)
        target_t = format(base_t + offset, "f")
        low = hardy_ball(target_t, PRECISIONS[0])
        high = hardy_ball(target_t, PRECISIONS[1])
        ctx.dps = PRECISIONS[1]
        low_real = arb(low["hardy_real_ball"])
        high_real = arb(high["hardy_real_ball"])
        require(low_real.overlaps(high_real), f"precision ladder missed at output {index}")
        source = arb(source_text)
        signed = source - high_real
        absolute = abs(signed)
        internal = arb(inner_row["source_q_free_complete_majorant_upper"])
        residual = absolute - internal
        rejected = absolute.lower() > REQUESTED_TOLERANCE
        internal_insufficient = residual.lower() > REQUESTED_TOLERANCE
        require(signed.upper() < 0, f"source-minus-Hardy sign is not strict at output {index}")
        require(rejected, f"requested tolerance is not rejected at output {index}")
        require(internal_insufficient, f"internal reverse-triangle residual does not reject tolerance at output {index}")
        expected_label = f"t{offset:+.2f}"
        require(inner_row["output_label"] == expected_label, f"output label drift at {index}")
        rows.append(
            {
                "output_index": index,
                "output_label": expected_label,
                "target_t": target_t,
                "source_binary128_decimal": source_text,
                "low_precision": low,
                "high_precision": high,
                "precision_enclosures_overlap": True,
                "signed_source_minus_hardy_ball": signed.str(PRECISIONS[1], more=True),
                "absolute_error_ball": absolute.str(PRECISIONS[1], more=True),
                "requested_absolute_tolerance": "0.005",
                "requested_tolerance_rejected": True,
                "current_internal_column_upper": inner_row["source_q_free_complete_majorant_upper"],
                "residual_after_any_internal_bounded_adjustment_ball": residual.str(PRECISIONS[1], more=True),
                "internal_column_alone_cannot_restore_tolerance": True,
            }
        )

    absolute_balls = [arb(row["absolute_error_ball"]) for row in rows]
    minimum_index = min(range(EXPECTED_OUTPUTS), key=lambda i: absolute_balls[i].lower())
    maximum_index = max(range(EXPECTED_OUTPUTS), key=lambda i: absolute_balls[i].upper())
    aggregate = {
        "output_count": len(rows),
        "precision_overlap_count": sum(row["precision_enclosures_overlap"] for row in rows),
        "strict_negative_source_minus_hardy_count": sum(arb(row["signed_source_minus_hardy_ball"]).upper() < 0 for row in rows),
        "requested_tolerance_rejected_count": sum(row["requested_tolerance_rejected"] for row in rows),
        "internal_column_insufficient_count": sum(row["internal_column_alone_cannot_restore_tolerance"] for row in rows),
        "minimum_absolute_error_witness_output": minimum_index + 1,
        "minimum_absolute_error_ball": rows[minimum_index]["absolute_error_ball"],
        "maximum_absolute_error_witness_output": maximum_index + 1,
        "maximum_absolute_error_ball": rows[maximum_index]["absolute_error_ball"],
        "absolute_0p005_target_holds_for_saved_source_run": False,
        "rh_implication": False,
    }
    require(aggregate["requested_tolerance_rejected_count"] == EXPECTED_OUTPUTS, "not all source outputs reject tolerance")
    require(aggregate["internal_column_insufficient_count"] == EXPECTED_OUTPUTS, "internal-column guard drift")

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_exact_hardy_arb_calibration_gate",
        "status": "rigorous_fifteen_point_arb_calibration_rejects_absolute_0p005_for_saved_source_run",
        "passed": True,
        "scope": {
            "height_center": "1e10",
            "offsets": "-0.07 through +0.07 in steps of 0.01",
            "outputs": EXPECTED_OUTPUTS,
            "precisions_decimal_digits": list(PRECISIONS),
        },
        "target_definition": {
            "hardy": "Z(t)=exp(i*theta(t))*zeta(1/2+i*t)",
            "theta": "theta(t)=Im(log Gamma(1/4+i*t/2))-t*log(pi)/2",
            "pi": "Arb mathematical pi",
        },
        "output_rows": rows,
        "aggregate": aggregate,
        "decision": {
            "saved_source_outputs_meet_absolute_0p005": False,
            "current_internal_column_alone_can_explain_gap": False,
            "outer_remainder_still_required": True,
            "finite_calibration_is_evidence_for_rh": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "dependencies": {
            "source_checkpoint": {"path": relative(CHECKPOINT), "sha256": file_hash(CHECKPOINT)},
            "outer_dependency_ledger": {"path": relative(OUTER_LEDGER), "sha256": file_hash(OUTER_LEDGER)},
            "internal_bound": {"path": relative(INNER_BOUND), "sha256": file_hash(INNER_BOUND)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": "Identify the outer source-to-Hardy discrepancy by starting from an exact classical Riemann-Siegel or Hardy identity with explicit constants. Use the fifteen signed residuals as a regression oracle, not as a fitted proof.",
        "proof_boundary": "Rigorous finite calibration of one saved fifteen-point source run only. It disproves that run's absolute 0.005 target, but proves no height-uniform error theorem, no Lambda<=0 statement, no PF-infinity statement, no RH, and no prize-level conclusion.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "validated exact Hardy Arb calibration: "
        f"15/15 reject 0.005, min={aggregate['minimum_absolute_error_ball']}, "
        f"max={aggregate['maximum_absolute_error_ball']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
