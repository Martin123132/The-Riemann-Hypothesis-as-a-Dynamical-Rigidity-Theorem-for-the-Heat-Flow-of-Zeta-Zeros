#!/usr/bin/env python3
"""Split the saved Hardy residual into lower-RS, hybrid-upper, and final-addition parts."""

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
HARDY_CALIBRATION = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_exact_hardy_arb_calibration_gate.json"
OUTER_LEDGER = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_outer_hardy_representation_dependency_ledger_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_lower_rs_hybrid_upper_component_split_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_lower_rs_hybrid_upper_component_split_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISIONS = (70, 110)
EXPECTED_OUTPUTS = 15
EXPECTED_NC = 621
EXPECTED_RSN = 39894
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


def load_source_stages() -> tuple[dict[str, Any], dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with CHECKPOINT.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                records.append(json.loads(line, parse_float=str))
    stage_one = [row for row in records if int(row["stage"]) == 1]
    stage_two = [row for row in records if int(row["stage"]) == 2]
    require(bool(stage_one) and len(stage_two) == 1, "source stage roster drift")
    return stage_one[-1], stage_two[0]


def exact_components(target_t: str, precision: int) -> dict[str, Any]:
    ctx.dps = precision
    ctx.threads = 1
    t = arb(target_t)
    pi = arb.pi()
    theta = acb(arb("0.25"), t / 2).lgamma().imag - t * pi.log() / 2
    hardy = acb(0, theta).exp() * acb(arb("0.5"), t).zeta()
    require(hardy.imag.contains(0), f"Hardy imaginary enclosure excludes zero at t={target_t}")

    lower_main = arb(0)
    for n in range(1, EXPECTED_NC + 1):
        nn = arb(n)
        lower_main += (theta - t * nn.log()).cos() / nn.sqrt()
    lower_main *= 2

    root = (t / (2 * pi)).sqrt()
    require(root > EXPECTED_RSN and root < EXPECTED_RSN + 1, f"Riemann-Siegel root bracket drift at t={target_t}")
    frac = root - EXPECTED_RSN
    endpoint_c = (2 * pi * (frac * (frac - 1) - arb(1) / 16)).cos() / (2 * pi * frac).cos()
    leading_endpoint = ((-1) ** (EXPECTED_RSN - 1)) * root.rsqrt() * endpoint_c
    lower = lower_main + leading_endpoint
    upper = hardy.real - lower
    return {
        "precision_decimal_digits": precision,
        "theta_ball": theta.str(precision, more=True),
        "hardy_real_ball": hardy.real.str(precision, more=True),
        "hardy_imag_ball": hardy.imag.str(precision, more=True),
        "rs_root_ball": root.str(precision, more=True),
        "rs_fraction_ball": frac.str(precision, more=True),
        "lower_main_ball": lower_main.str(precision, more=True),
        "leading_endpoint_ball": leading_endpoint.str(precision, more=True),
        "exact_lower_component_ball": lower.str(precision, more=True),
        "exact_complementary_upper_ball": upper.str(precision, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    aggregate = artifact["aggregate"]
    return f"""# Lower-RS versus hybrid-upper residual split

Date: 2026-08-09

Status: finite interval diagnostic validated; not a proof of RH

For each of the fifteen saved heights, Arb reconstructs the exact mathematical
version of the source's visible lower component,

```text
L(t)=2 sum_(n=1)^621 cos(theta(t)-t log n)/sqrt(n)+C_0(t),
```

where `C_0(t)` is the leading Riemann--Siegel endpoint formula used by the
source.  The exact complementary upper target is defined by `U(t)=Z(t)-L(t)`.
The source values are read separately from the last stage-1 `ZP` checkpoint,
the stage-2 `rszsum`, and the final rounded sum.

The residual therefore splits identically as

```text
S-Z = (R_source-L) + (P_source-U) + (S-P_source-R_source).
```

This is an exact identity; the interval calculation only encloses its
components.

The 70- and 110-digit component enclosures overlap on all
{aggregate['output_count']} outputs, and every reconstructed identity interval
contains zero.  The independent checker repeats the calculation at higher
precision.

Lower-component dominance count: {aggregate['lower_component_dominant_count']}.
Hybrid-upper dominance count: {aggregate['hybrid_upper_component_dominant_count']}.
Final-addition dominance count: {aggregate['final_addition_component_dominant_count']}.
Outputs whose hybrid-upper discrepancy alone rigorously exceeds `0.005`:
{aggregate['hybrid_upper_alone_rejects_0p005_count']}.

The largest lower, upper, and final-addition discrepancy balls are recorded in
the machine artifact.  This finite split localizes the saved error; it does not
turn `U(t)` into a proved hybrid formula and does not identify which analytic
term inside block zero, blocks 1--19, `H(t)`, or the outer truncations is
responsible.  The next step must refine the dominant component from an exact
identity with explicit remainders rather than fit the observed residual.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    stage_one, stage_two = load_source_stages()
    calibration = json.loads(HARDY_CALIBRATION.read_text(encoding="utf-8"))
    outer = json.loads(OUTER_LEDGER.read_text(encoding="utf-8"))
    require(int(stage_one["nc"]) == EXPECTED_NC and int(stage_two["nc"]) == EXPECTED_NC, "NC drift")
    require(int(stage_one["nchalf"]) == EXPECTED_OUTPUTS and int(stage_two["nchalf"]) == EXPECTED_OUTPUTS, "output roster drift")
    require(len(stage_one["zsum"]) == EXPECTED_OUTPUTS, "stage-one ZP roster drift")
    require(len(stage_two["rszsum"]) == EXPECTED_OUTPUTS and len(stage_two["zsum"]) == EXPECTED_OUTPUTS, "stage-two roster drift")
    require(calibration["aggregate"]["requested_tolerance_rejected_count"] == EXPECTED_OUTPUTS, "Hardy calibration dependency drift")
    require(outer["aggregate"]["outputs_hardy_z_certified"] == 0, "outer ledger was overpromoted")

    base_t = Decimal(stage_two["t"])
    rows: list[dict[str, Any]] = []
    for index in range(1, EXPECTED_OUTPUTS + 1):
        offset = Decimal(index - 8) / Decimal(100)
        target_t = format(base_t + offset, "f")
        low = exact_components(target_t, PRECISIONS[0])
        high = exact_components(target_t, PRECISIONS[1])
        ctx.dps = PRECISIONS[1]
        overlap_fields = (
            "hardy_real_ball",
            "exact_lower_component_ball",
            "exact_complementary_upper_ball",
        )
        require(all(arb(low[field]).overlaps(arb(high[field])) for field in overlap_fields), f"precision ladder missed at output {index}")

        source_upper = arb(stage_one["zsum"][index - 1])
        source_lower = arb(stage_two["rszsum"][index - 1])
        source_final = arb(stage_two["zsum"][index - 1])
        exact_hardy = arb(high["hardy_real_ball"])
        exact_lower = arb(high["exact_lower_component_ball"])
        exact_upper = arb(high["exact_complementary_upper_ball"])

        lower_gap = source_lower - exact_lower
        upper_gap = source_upper - exact_upper
        addition_gap = source_final - source_upper - source_lower
        direct_total = source_final - exact_hardy
        component_total = lower_gap + upper_gap + addition_gap
        identity_residual = component_total - direct_total
        require(identity_residual.contains(0), f"component identity excludes zero at output {index}")

        magnitudes = {
            "lower": abs(lower_gap),
            "upper": abs(upper_gap),
            "addition": abs(addition_gap),
        }
        dominant = max(magnitudes, key=lambda key: magnitudes[key].upper())
        require(
            magnitudes[dominant].lower() > max(magnitudes[key].upper() for key in magnitudes if key != dominant),
            f"component dominance is not rigorous at output {index}",
        )
        rows.append(
            {
                "output_index": index,
                "output_label": f"t{offset:+.2f}",
                "target_t": target_t,
                "source_final_binary128_decimal": stage_two["zsum"][index - 1],
                "source_lower_rszsum_binary128_decimal": stage_two["rszsum"][index - 1],
                "source_upper_zp_binary128_decimal": stage_one["zsum"][index - 1],
                "low_precision": low,
                "high_precision": high,
                "precision_enclosures_overlap": True,
                "lower_source_minus_exact_ball": lower_gap.str(PRECISIONS[1], more=True),
                "hybrid_upper_source_minus_exact_ball": upper_gap.str(PRECISIONS[1], more=True),
                "final_addition_rounding_ball": addition_gap.str(PRECISIONS[1], more=True),
                "direct_source_minus_hardy_ball": direct_total.str(PRECISIONS[1], more=True),
                "component_sum_ball": component_total.str(PRECISIONS[1], more=True),
                "identity_residual_ball": identity_residual.str(PRECISIONS[1], more=True),
                "dominant_component": dominant,
                "hybrid_upper_alone_rejects_0p005": magnitudes["upper"].lower() > REQUESTED_TOLERANCE,
            }
        )

    def maximum_witness(field: str) -> tuple[int, arb]:
        balls = [abs(arb(row[field])) for row in rows]
        index = max(range(len(balls)), key=lambda i: balls[i].upper())
        return index + 1, balls[index]

    def minimum_witness(field: str) -> tuple[int, arb]:
        balls = [abs(arb(row[field])) for row in rows]
        index = min(range(len(balls)), key=lambda i: balls[i].lower())
        return index + 1, balls[index]

    lower_max_index, lower_max = maximum_witness("lower_source_minus_exact_ball")
    upper_min_index, upper_min = minimum_witness("hybrid_upper_source_minus_exact_ball")
    upper_max_index, upper_max = maximum_witness("hybrid_upper_source_minus_exact_ball")
    addition_max_index, addition_max = maximum_witness("final_addition_rounding_ball")
    aggregate = {
        "output_count": len(rows),
        "precision_overlap_count": sum(row["precision_enclosures_overlap"] for row in rows),
        "identity_contains_zero_count": sum(arb(row["identity_residual_ball"]).contains(0) for row in rows),
        "lower_component_dominant_count": sum(row["dominant_component"] == "lower" for row in rows),
        "hybrid_upper_component_dominant_count": sum(row["dominant_component"] == "upper" for row in rows),
        "final_addition_component_dominant_count": sum(row["dominant_component"] == "addition" for row in rows),
        "hybrid_upper_alone_rejects_0p005_count": sum(row["hybrid_upper_alone_rejects_0p005"] for row in rows),
        "maximum_lower_component_absolute_ball": lower_max.str(PRECISIONS[1], more=True),
        "maximum_lower_component_witness_output": lower_max_index,
        "minimum_hybrid_upper_absolute_ball": upper_min.str(PRECISIONS[1], more=True),
        "minimum_hybrid_upper_witness_output": upper_min_index,
        "maximum_hybrid_upper_absolute_ball": upper_max.str(PRECISIONS[1], more=True),
        "maximum_hybrid_upper_witness_output": upper_max_index,
        "minimum_hybrid_upper_to_maximum_lower_ratio_ball": (upper_min / lower_max).str(PRECISIONS[1], more=True),
        "maximum_final_addition_absolute_ball": addition_max.str(PRECISIONS[1], more=True),
        "maximum_final_addition_witness_output": addition_max_index,
    }
    require(aggregate["precision_overlap_count"] == EXPECTED_OUTPUTS, "precision overlap aggregate drift")
    require(aggregate["identity_contains_zero_count"] == EXPECTED_OUTPUTS, "identity aggregate drift")
    require(
        aggregate["lower_component_dominant_count"]
        + aggregate["hybrid_upper_component_dominant_count"]
        + aggregate["final_addition_component_dominant_count"]
        == EXPECTED_OUTPUTS,
        "dominance aggregate drift",
    )

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_lower_rs_hybrid_upper_component_split_gate",
        "status": "rigorous_fifteen_point_lower_rs_hybrid_upper_and_final_addition_residual_split",
        "passed": True,
        "scope": {
            "height_center": "1e10",
            "outputs": EXPECTED_OUTPUTS,
            "lower_sum_cutoff": EXPECTED_NC,
            "riemann_siegel_index": EXPECTED_RSN,
            "precisions_decimal_digits": list(PRECISIONS),
        },
        "component_definition": {
            "exact_lower": "L=2*sum_(n=1)^621 cos(theta-t*log(n))/sqrt(n)+leading source-style RS endpoint term",
            "exact_upper": "U=Z-L",
            "source_upper": "last stage-1 checkpoint zsum before RS addition",
            "source_lower": "stage-2 rszsum after the leading endpoint term",
            "source_final": "stage-2 final zsum after binary128 addition",
            "identity": "S-Z=(R_source-L)+(P_source-U)+(S-P_source-R_source)",
        },
        "output_rows": rows,
        "aggregate": aggregate,
        "decision": {
            "component_split_is_exact_identity": True,
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
            "source_checkpoint": {"path": relative(CHECKPOINT), "sha256": file_hash(CHECKPOINT)},
            "exact_hardy_calibration": {"path": relative(HARDY_CALIBRATION), "sha256": file_hash(HARDY_CALIBRATION)},
            "outer_dependency_ledger": {"path": relative(OUTER_LEDGER), "sha256": file_hash(OUTER_LEDGER)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "next_obligation": "Refine the rigorously dominant finite component into exact analytic subcolumns, beginning from a primary exact representation with explicit constants. Retain the signed residual vector only as a regression oracle.",
        "proof_boundary": "Rigorous finite component localization for one saved fifteen-point run only. It proves no exact global hybrid formula, no height-uniform error theorem, no Lambda<=0 statement, no PF-infinity statement, no RH, and no prize-level conclusion.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "validated lower-RS/hybrid-upper split: "
        f"15 identities, lower-dominant={aggregate['lower_component_dominant_count']}, "
        f"upper-dominant={aggregate['hybrid_upper_component_dominant_count']}, "
        f"upper-rejects-0.005={aggregate['hybrid_upper_alone_rejects_0p005_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
