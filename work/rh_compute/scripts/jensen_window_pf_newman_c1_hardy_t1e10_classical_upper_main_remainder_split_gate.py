#!/usr/bin/env python3
"""Compare source ZP with the classical upper Riemann-Siegel main sum."""

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
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_lower_rs_hybrid_upper_component_split_gate as component_split


CHECKPOINT = component_split.CHECKPOINT
COMPONENT_SPLIT = component_split.RESULT
HARDY_CALIBRATION = component_split.HARDY_CALIBRATION
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_classical_upper_main_remainder_split_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_classical_upper_main_remainder_split_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISIONS = (70, 110)
LOWER_START = 622
UPPER_END = 39894
EXPECTED_OUTPUTS = 15
REQUESTED_TOLERANCE = arb("0.005")


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


def classical_upper_main(target_t: str, precision: int) -> arb:
    ctx.dps = precision
    ctx.threads = 1
    t = arb(target_t)
    theta = acb(arb("0.25"), t / 2).lgamma().imag - t * arb.pi().log() / 2
    total = arb(0)
    for n in range(LOWER_START, UPPER_END + 1):
        nn = arb(n)
        total += (theta - t * nn.log()).cos() / nn.sqrt()
    return 2 * total


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Classical upper-main versus hybrid-ZP split

Date: 2026-08-09

Status: finite interval diagnostic validated; not a proof of RH

The exact complementary upper target from Section 11.297 is split once more:

```text
U(t)=M_upper(t)+rho_RS(t),
M_upper(t)=2 sum_(n=622)^39894 cos(theta(t)-t log n)/sqrt(n).
```

Here `rho_RS=U-M_upper` is defined exactly.  It contains all classical
correction and remainder content left after the source-style leading endpoint
term already included in the lower component.  The source hybrid residual
therefore satisfies the exact identity

```text
P_source-U=(P_source-M_upper)-rho_RS.
```

The 70- and 110-digit Arb enclosures overlap for all fifteen upper main sums
and remainder balls; the independent checker repeats them at higher precision.
All {aggregate['hybrid_minus_classical_main_rejects_0p005_count']} source
hybrid-minus-classical-main discrepancies exceed `0.005`.

Maximum exact remaining correction/remainder:

```text
{aggregate['maximum_exact_remaining_correction_absolute_ball']}
```

Minimum source-hybrid versus classical-main discrepancy:

```text
{aggregate['minimum_hybrid_minus_classical_main_absolute_ball']}
```

The minimum discrepancy-to-remainder ratio is

```text
{aggregate['minimum_hybrid_gap_to_maximum_remainder_ratio_ball']}
```

Thus the tiny exact Riemann--Siegel correction/remainder left outside the
classical upper main sum cannot explain the saved `ZP` failure.  The error is
inside the source's approximation to the upper main sum itself.  This finite
localization does not identify which hybrid block, transition, normalization,
or asymptotic replacement creates it and has no RH implication.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    split = json.loads(COMPONENT_SPLIT.read_text(encoding="utf-8"))
    calibration = json.loads(HARDY_CALIBRATION.read_text(encoding="utf-8"))
    require(split["aggregate"]["hybrid_upper_component_dominant_count"] == EXPECTED_OUTPUTS, "component-split dependency drift")
    require(calibration["aggregate"]["requested_tolerance_rejected_count"] == EXPECTED_OUTPUTS, "Hardy dependency drift")
    require(len(split["output_rows"]) == EXPECTED_OUTPUTS, "split output roster drift")

    rows: list[dict[str, Any]] = []
    for source_row in split["output_rows"]:
        index = int(source_row["output_index"])
        target_t = source_row["target_t"]
        low_main = classical_upper_main(target_t, PRECISIONS[0])
        low_components = component_split.exact_components(target_t, PRECISIONS[0])
        high_main = classical_upper_main(target_t, PRECISIONS[1])
        high_components = component_split.exact_components(target_t, PRECISIONS[1])
        ctx.dps = PRECISIONS[1]
        low_upper = arb(low_components["exact_complementary_upper_ball"])
        high_upper = arb(high_components["exact_complementary_upper_ball"])
        require(low_main.overlaps(high_main), f"upper-main precision ladder missed at output {index}")
        require(low_upper.overlaps(high_upper), f"upper-target precision ladder missed at output {index}")

        exact_remainder = high_upper - high_main
        source_hybrid = arb(source_row["source_upper_zp_binary128_decimal"])
        hybrid_minus_main = source_hybrid - high_main
        hybrid_minus_target = source_hybrid - high_upper
        identity_residual = hybrid_minus_main - exact_remainder - hybrid_minus_target
        require(identity_residual.contains(0), f"upper identity excludes zero at output {index}")
        main_gap_dominates = abs(hybrid_minus_main).lower() > abs(exact_remainder).upper()
        rejects = abs(hybrid_minus_main).lower() > REQUESTED_TOLERANCE
        require(main_gap_dominates, f"classical remainder is not rigorously subordinate at output {index}")
        require(rejects, f"hybrid-minus-main does not reject tolerance at output {index}")
        rows.append(
            {
                "output_index": index,
                "output_label": source_row["output_label"],
                "target_t": target_t,
                "source_upper_zp_binary128_decimal": source_row["source_upper_zp_binary128_decimal"],
                "low_precision": {
                    "precision_decimal_digits": PRECISIONS[0],
                    "classical_upper_main_ball": low_main.str(PRECISIONS[0], more=True),
                    "exact_complementary_upper_ball": low_components["exact_complementary_upper_ball"],
                },
                "high_precision": {
                    "precision_decimal_digits": PRECISIONS[1],
                    "classical_upper_main_ball": high_main.str(PRECISIONS[1], more=True),
                    "exact_complementary_upper_ball": high_components["exact_complementary_upper_ball"],
                },
                "precision_enclosures_overlap": True,
                "exact_remaining_correction_ball": exact_remainder.str(PRECISIONS[1], more=True),
                "hybrid_minus_classical_main_ball": hybrid_minus_main.str(PRECISIONS[1], more=True),
                "hybrid_minus_exact_upper_ball": hybrid_minus_target.str(PRECISIONS[1], more=True),
                "identity_residual_ball": identity_residual.str(PRECISIONS[1], more=True),
                "hybrid_main_gap_dominates_exact_remainder": True,
                "hybrid_minus_classical_main_rejects_0p005": True,
            }
        )

    remainder_balls = [abs(arb(row["exact_remaining_correction_ball"])) for row in rows]
    hybrid_balls = [abs(arb(row["hybrid_minus_classical_main_ball"])) for row in rows]
    remainder_max_index = max(range(EXPECTED_OUTPUTS), key=lambda i: remainder_balls[i].upper())
    hybrid_min_index = min(range(EXPECTED_OUTPUTS), key=lambda i: hybrid_balls[i].lower())
    hybrid_max_index = max(range(EXPECTED_OUTPUTS), key=lambda i: hybrid_balls[i].upper())
    remainder_max = remainder_balls[remainder_max_index]
    hybrid_min = hybrid_balls[hybrid_min_index]
    aggregate = {
        "output_count": len(rows),
        "precision_overlap_count": sum(row["precision_enclosures_overlap"] for row in rows),
        "identity_contains_zero_count": sum(arb(row["identity_residual_ball"]).contains(0) for row in rows),
        "remainder_subordinate_count": sum(row["hybrid_main_gap_dominates_exact_remainder"] for row in rows),
        "hybrid_minus_classical_main_rejects_0p005_count": sum(row["hybrid_minus_classical_main_rejects_0p005"] for row in rows),
        "maximum_exact_remaining_correction_absolute_ball": remainder_max.str(PRECISIONS[1], more=True),
        "maximum_exact_remaining_correction_witness_output": remainder_max_index + 1,
        "minimum_hybrid_minus_classical_main_absolute_ball": hybrid_min.str(PRECISIONS[1], more=True),
        "minimum_hybrid_minus_classical_main_witness_output": hybrid_min_index + 1,
        "maximum_hybrid_minus_classical_main_absolute_ball": hybrid_balls[hybrid_max_index].str(PRECISIONS[1], more=True),
        "maximum_hybrid_minus_classical_main_witness_output": hybrid_max_index + 1,
        "minimum_hybrid_gap_to_maximum_remainder_ratio_ball": (hybrid_min / remainder_max).str(PRECISIONS[1], more=True),
    }
    require(aggregate["precision_overlap_count"] == EXPECTED_OUTPUTS, "precision aggregate drift")
    require(aggregate["identity_contains_zero_count"] == EXPECTED_OUTPUTS, "identity aggregate drift")
    require(aggregate["remainder_subordinate_count"] == EXPECTED_OUTPUTS, "remainder dominance aggregate drift")
    require(aggregate["hybrid_minus_classical_main_rejects_0p005_count"] == EXPECTED_OUTPUTS, "tolerance aggregate drift")

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_classical_upper_main_remainder_split_gate",
        "status": "rigorous_fifteen_point_source_zp_vs_classical_upper_main_and_exact_remaining_correction_split",
        "passed": True,
        "scope": {
            "height_center": "1e10",
            "outputs": EXPECTED_OUTPUTS,
            "classical_upper_index_start": LOWER_START,
            "classical_upper_index_end": UPPER_END,
            "classical_upper_summands": UPPER_END - LOWER_START + 1,
            "precisions_decimal_digits": list(PRECISIONS),
        },
        "identity": {
            "exact_upper": "U=M_upper+rho_RS",
            "classical_upper_main": "M_upper=2*sum_(n=622)^39894 cos(theta-t*log(n))/sqrt(n)",
            "exact_remaining_correction": "rho_RS=U-M_upper",
            "hybrid_residual": "P_source-U=(P_source-M_upper)-rho_RS",
        },
        "output_rows": rows,
        "aggregate": aggregate,
        "decision": {
            "exact_remaining_correction_explains_saved_zp_error": False,
            "saved_error_is_inside_hybrid_approximation_to_upper_main": True,
            "finite_split_is_height_uniform_theorem": False,
            "finite_split_has_rh_implication": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "dependencies": {
            "component_split": {"path": relative(COMPONENT_SPLIT), "sha256": file_hash(COMPONENT_SPLIT)},
            "exact_hardy_calibration": {"path": relative(HARDY_CALIBRATION), "sha256": file_hash(HARDY_CALIBRATION)},
            "source_checkpoint": {"path": relative(CHECKPOINT), "sha256": file_hash(CHECKPOINT)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "component_evaluator": {"path": relative(Path(component_split.__file__).resolve()), "sha256": file_hash(Path(component_split.__file__).resolve())},
        },
        "next_obligation": "Split source ZP into block-zero and per-block increments and derive an exact classical-upper target partition matching those blocks, or replace the hybrid route by a new exact transformation with explicit per-block remainders.",
        "proof_boundary": "Rigorous finite localization for one saved fifteen-point run only. It proves no exact global hybrid transformation, no per-block analytic cause, no height-uniform error theorem, no Lambda<=0 statement, no PF-infinity statement, no RH, and no prize-level conclusion.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "validated classical upper-main split: "
        f"15 identities, remainder-subordinate={aggregate['remainder_subordinate_count']}, "
        f"hybrid-main-rejects-0.005={aggregate['hybrid_minus_classical_main_rejects_0p005_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
