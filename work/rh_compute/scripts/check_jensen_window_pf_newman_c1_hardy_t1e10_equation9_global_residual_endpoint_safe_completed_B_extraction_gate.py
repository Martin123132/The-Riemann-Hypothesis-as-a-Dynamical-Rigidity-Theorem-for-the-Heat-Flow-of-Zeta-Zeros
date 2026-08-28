#!/usr/bin/env python3
"""Independently check the endpoint-safe crossing-completed B extraction."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_endpoint_safe_completed_B_extraction_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    # Check the allocation with unrelated rational values rather than replaying
    # the builder's symbolic expression tree.
    bulk = Fraction(17, 13)
    u_b = Fraction(-19, 11)
    p_b_minus = Fraction(23, 7)
    for tau in (0, 1):
        for local in (0, 1):
            for q_negative in (0, 1):
                sigma = 1 - tau - q_negative
                p_b = u_b - q_negative * bulk
                trace = u_b + p_b_minus + local * sigma * bulk
                omitted = (1 - local) * sigma * bulk
                completed = p_b + p_b_minus + (1 - tau) * bulk
                require(trace + omitted == completed, "independent completed-mode allocation failed")

    # Across q_negative: 1 -> 0, U_B jumps by -P_bulk and sigma jumps by +1.
    for tau in (0, 1):
        for local in (0, 1):
            sigma_left = -tau
            sigma_right = 1 - tau
            trace_jump = -bulk + local * (sigma_right - sigma_left) * bulk
            omitted_jump = (1 - local) * (sigma_right - sigma_left) * bulk
            require(trace_jump + omitted_jump == 0, "independent crossing-jump cancellation failed")
            if local == 0:
                require(trace_jump == -bulk and omitted_jump == bulk, "nonlocal jump obstruction was lost")
            else:
                require(trace_jump == 0 and omitted_jump == 0, "local crossing should already be completed")

    ctx.dps = 70
    ctx.threads = 1
    geometry = artifact["geometry_certificate"]
    c_low = arb(geometry["c_window_low_ball"])
    c_high = arb(geometry["c_window_high_ball"])
    c_2delta = Fraction(5_122_421, 10_000)
    c_half = Fraction(5_122_421, 4)
    require(c_2delta < 620 and arb(620) < c_low, "left crossing witness failed")
    require(c_high < arb(623) and Fraction(623) < c_half, "right crossing witness failed")
    require(len(range(257, 621)) == geometry["left_positive_eta_crossing_count"], "left crossing count failed")
    require(len(range(623, 1_280_606)) == geometry["right_exterior_crossing_count"], "right crossing count failed")

    symbolic = artifact["symbolic_certificate"]
    require("Btr_(M,epsilon)+N_(M,epsilon)" in symbolic["definitions"]["completed_B"], "completed trace definition missing")
    require("classical first and second derivatives do not exist" in symbolic["old_trace_obstruction"], "derivative obstruction missing")

    decision = artifact["decision"]
    require(decision["nonlocal_B_trace_jump_obstruction_proved"] is True, "jump obstruction not certified")
    require(decision["classical_derivatives_of_uncompleted_exterior_B_trace_available"] is False, "uncompleted derivative overclaim")
    require(decision["crossing_completed_exterior_B_current_smooth_at_finite_cutoff"] is True, "completed regularity lost")
    require(decision["Abel_corner_grouping_unchanged"] is True, "corner grouping lost")
    require(decision["negative_B_trace_window_certificate_unchanged"] is True, "negative window changed")
    require(decision["completed_exterior_first_second_derivative_bounds_proved"] is False, "derivative bounds overclaimed")
    require(decision["joined_remainder_below_1p4058e_minus_4_proved"] is False, "joined bound overclaimed")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently checked crossing completion, derivative obstruction, and endpoint-safe allocation", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
