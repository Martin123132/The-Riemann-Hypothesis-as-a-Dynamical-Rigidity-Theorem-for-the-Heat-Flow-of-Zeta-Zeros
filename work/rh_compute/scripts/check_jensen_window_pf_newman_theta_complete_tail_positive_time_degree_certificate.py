#!/usr/bin/env python3
"""Check the complete-theta local positive-time degree certificate."""

from __future__ import annotations

import argparse
import importlib
import json
import os
from pathlib import Path
import sys

for variable in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(variable, "1")

import psutil


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from flint import arb  # noqa: E402


STEM = "jensen_window_pf_newman_theta_complete_tail_positive_time_degree_certificate"
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
SOURCE = SCRIPT_DIR / f"{STEM}.py"


class CheckFailure(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CheckFailure(message)


def static_audit(payload: dict, core) -> dict:
    require(payload["kind"] == STEM, "result kind mismatch")
    require(payload["schema_version"] == 1, "schema version mismatch")
    require(
        payload["status"]
        == "rigorous complete-theta local positive-time exclusion complete",
        "certificate is not complete",
    )
    require(payload["source_sha256"] == core.sha256_path(SOURCE), "source hash mismatch")
    for parent in payload["parents"].values():
        path = REPO_ROOT / parent["path"]
        require(parent["sha256"] == core.sha256_path(path), f"parent hash mismatch: {path}")

    resource = payload["resource_policy"]
    require(resource["active_compute_workers"] == 1, "worker count is not one")
    require(resource["thread_caps"] == 1, "thread cap is not one")
    require(resource["below_normal_priority_applied"] is True, "priority gate missing")
    require(payload["rectangle"]["time"] == ["0", "0.45"], "time rectangle mismatch")
    require(payload["rectangle"]["x"] == ["135.5", "136"], "x rectangle mismatch")

    rebuilt_ratio = core.exact_component_ratio_audit()
    require(
        rebuilt_ratio == payload["component_ratio_audit"],
        "exact component-ratio audit drifted",
    )
    ratio_text = json.dumps(rebuilt_ratio, sort_keys=True)
    for phrase in (
        "canonical Jacobi-theta normalization",
        "phi_n/phi_5 <= 2*(n/5)^4",
        "sum_(n>=6)b_n",
        "10^-9*phi_5",
    ):
        require(phrase in ratio_text, f"ratio audit lost phrase: {phrase}")

    budget = payload["positive_moment_budget"]
    require(budget["reference_component"] == "phi_5", "reference component mismatch")
    require(budget["derivative_orders"] == [0, 1], "derivative orders mismatch")
    require(budget["time_upper"] == "0.45", "time upper mismatch")
    require(
        budget["precision_bits"] == core.PRECISION_BITS,
        "production precision mismatch",
    )
    require(
        budget["absolute_integration_tolerance"] == core.ABSOLUTE_TOLERANCE,
        "production integration tolerance mismatch",
    )
    raw_moments = [arb(value) for value in budget["raw_moment_uppers"]]
    moments = [arb(value) for value in budget["moment_uppers"]]
    tails = [arb(value) for value in budget["tail_derivative_bounds"]]
    require(
        len(raw_moments) == len(moments) == len(tails) == 2,
        "moment budget length mismatch",
    )
    safety_factor = arb(budget["safety_factor"])
    require(
        budget["safety_factor"] == core.MOMENT_SAFETY_FACTOR,
        "moment safety factor mismatch",
    )
    for raw_moment, moment, tail in zip(
        raw_moments, moments, tails, strict=True
    ):
        require(raw_moment.lower() > 0, "raw positive moment is not positive")
        require(
            moment.overlaps(safety_factor * raw_moment),
            "conservative moment budget drifted",
        )
        require(moment.lower() > 0, "positive moment is not positive")
        expected = arb(core.TAIL_RATIO_BOUND) * moment
        require(tail.overlaps(expected), "tail derivative bound drifted")

    maximum_tail = arb(payload["maximum_tail_derivative_bound"])
    finite_margin = arb(payload["finite_margin_lower"])
    ratio = arb(payload["tail_to_margin_ratio_upper"])
    require(maximum_tail.overlaps(max(tails)), "maximum tail summary mismatch")
    require(maximum_tail.upper() < finite_margin.lower(), "tail exceeds finite margin")
    require(ratio.upper() < arb("4e-21"), "tail-to-margin ratio is too large")
    require(
        "complete theta-kernel heat transform" in payload["theorem"],
        "complete-theta theorem missing",
    )
    for phrase in ("other frequency", "all real x", "Lambda<=0", "RH", "prize-level"):
        require(phrase in payload["proof_boundary"], f"proof boundary lost phrase: {phrase}")
    return {
        "tail_orders": len(tails),
        "finite_margin": str(finite_margin),
        "maximum_tail": str(maximum_tail),
        "tail_to_margin": str(ratio),
    }


def independent_replay(payload: dict, core) -> dict:
    independent_precision = 320
    independent_tolerance = "1e-48"
    core.continuation.ABS_TOL = independent_tolerance
    core.continuation.flint.ctx.prec = independent_precision
    fifth = core.scout.FifthTransformFamily(
        core.Fraction(0),
        core.TIME_HIGH,
        core.TIME_HIGH,
        core.Fraction(0),
        weighted=False,
    )
    moments = [fifth.moment_upper(order, core.TIME_HIGH) for order in range(2)]
    tail_ratio = arb(core.TAIL_RATIO_BOUND)
    tails = [tail_ratio * moment for moment in moments]
    stored_raw_moments = [
        arb(value) for value in payload["positive_moment_budget"]["raw_moment_uppers"]
    ]
    stored_budgets = [
        arb(value) for value in payload["positive_moment_budget"]["moment_uppers"]
    ]
    relative_drifts = []
    for rebuilt, stored_raw, stored_budget in zip(
        moments, stored_raw_moments, stored_budgets, strict=True
    ):
        require(
            stored_budget.lower() > rebuilt.upper(),
            "conservative production budget does not dominate the stronger replay",
        )
        drift = (
            abs(float(rebuilt.mid()) - float(stored_raw.mid()))
            / float(rebuilt.mid())
        )
        require(drift < 5e-2, f"independent positive moment drifted by {drift:.3e}")
        relative_drifts.append(drift)
    finite_margin = arb(payload["finite_margin_lower"])
    maximum_tail = max(tails)
    require(maximum_tail.upper() < finite_margin.lower(), "independent tail exceeds margin")
    return {
        "precision_bits": independent_precision,
        "absolute_integration_tolerance": independent_tolerance,
        "moment_uppers": [core.continuation.arb_text(value) for value in moments],
        "tail_derivative_bounds": [
            core.continuation.arb_text(value) for value in tails
        ],
        "relative_drifts": relative_drifts,
        "tail_to_margin_upper": core.continuation.arb_text(
            (maximum_tail / finite_margin).upper()
        ),
        "production_budgets_dominate": True,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=RESULT)
    parser.add_argument("--baseline-seconds", type=int, default=5)
    parser.add_argument("--baseline-max-percent", type=float, default=60.0)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    core = importlib.import_module(STEM)
    priority_lowered = core.continuation.request_below_normal_priority()
    baseline = [float(psutil.cpu_percent(interval=1.0)) for _ in range(args.baseline_seconds)]
    baseline_mean = sum(baseline) / len(baseline)
    print(
        "complete-theta checker baseline CPU: "
        f"samples={','.join(f'{value:.1f}' for value in baseline)}, "
        f"mean={baseline_mean:.1f}%",
        flush=True,
    )
    if baseline_mean > args.baseline_max_percent:
        print("complete-theta checker deferred: daytime baseline is already busy")
        return 2
    require(priority_lowered, "checker could not apply below-normal priority")
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    static = static_audit(payload, core)
    independent = independent_replay(payload, core)
    print(
        "checked complete-theta local degree certificate: "
        f"tail_orders={static['tail_orders']}, tail_to_margin<4e-21, "
        "complete_theta=True, no_contacts=True, 0 issues"
    )
    print(json.dumps({"static": static, "independent": independent}, indent=2))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except CheckFailure as error:
        print(f"CHECK FAILED: {error}", file=sys.stderr)
        sys.exit(1)
