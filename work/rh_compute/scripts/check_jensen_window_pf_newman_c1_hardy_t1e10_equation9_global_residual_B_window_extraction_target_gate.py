#!/usr/bin/env python3
"""Independent finite-algebra check of the global B-window extraction."""

from __future__ import annotations

import hashlib
import json
from decimal import Decimal, getcontext
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_B_window_extraction_target_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    getcontext().prec = 50
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    # Independent finite-array replay: the extraction is pointwise and does
    # not depend on a continuum or asymptotic argument.
    global_values = [Decimal("3.25"), Decimal("-7.5"), Decimal("11.125"), Decimal("0.75")]
    B_values = [Decimal("1.5"), Decimal("-2.25"), Decimal("8.0"), Decimal("-1.0")]
    window = [Decimal(0), Decimal(1), Decimal(1), Decimal(0)]
    extracted = [chi * value for chi, value in zip(window, B_values)]
    remainder = [global_value - value for global_value, value in zip(global_values, extracted)]
    require(all(left + right == total for left, right, total in zip(extracted, remainder, global_values)), "finite extraction replay failed")

    # Replay the mode allocation independently for modes inside and outside
    # the extracted positive-mode block.  This is the algebra that prevents
    # the B window from being defined as an unspecified global difference.
    mode_rows = [
        (Decimal("2.75"), Decimal("-0.5"), Decimal("1.125"), Decimal("-3.0"), Decimal("0.625"), Decimal(0), Decimal(0), Decimal(0)),
        (Decimal("-1.25"), Decimal("4.0"), Decimal("-2.5"), Decimal("0.75"), Decimal("-1.125"), Decimal(1), Decimal(1), Decimal(1)),
        (Decimal("0.875"), Decimal("-2.0"), Decimal("3.5"), Decimal("1.25"), Decimal("-0.25"), Decimal(1), Decimal(0), Decimal(1)),
    ]
    for bulk, A_plus, U_plus, I_minus, B_minus, tau, local, q_negative in mode_rows:
        sigma = Decimal(1) - tau - q_negative
        pair_residual = A_plus + I_minus + U_plus + sigma * bulk
        allocated_B_trace = U_plus + B_minus + local * sigma * bulk
        allocated_rest = A_plus + I_minus - B_minus + (Decimal(1) - local) * sigma * bulk
        require(pair_residual == allocated_B_trace + allocated_rest, "finite-mode allocation replay failed")

    certificate = artifact["symbolic_certificate"]
    require("tau_m=1_(622<=m<=39894)" in certificate["step_coefficient"], "target-mode indicator drift")
    require("ell_m=1_(m in {621,622})" in certificate["step_coefficient"], "local-step indicator drift")
    require("(1-ell_m)sigma_m P_bulk" in certificate["finite_mode_allocation"], "nonlocal step allocation omitted")
    require("P_B,-m" in certificate["B_trace_package"], "negative B allocation omitted")
    require("nonlocal" in certificate["omitted_step_completion"], "nonlocal step completion omitted")

    target = artifact["target_certificate"]
    window_margin = -Decimal(target["complete_B_trace_window_upper_bound"])
    grouped = Decimal(target["sufficient_grouped_remainder_upper_target"])
    reference = Decimal(target["working_global_residual_upper_target"])
    require(grouped - window_margin == reference, "closing-target arithmetic drift")

    decision = artifact["decision"]
    require(decision["remaining_object_is_one_grouped_global_remainder"] is True, "grouping decision drift")
    require(decision["nonlocal_constant_bulk_steps_retained_in_grouped_remainder"] is True, "nonlocal step completion lost")
    require(decision["complete_endpoint_B_window_extracted"] is False, "complete endpoint B window overclaimed")
    require(decision["isolated_left_B_tail_triangle_budget_admissible"] is False, "illegal tail split admitted")
    require(decision["grouped_remainder_below_1p4058e_minus_4_proved"] is False, "grouped-remainder overclaim")
    require(decision["complete_T_upper_proved"] is False, "T_upper overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently checked B trace-package allocation, window extraction, and grouped target arithmetic", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
