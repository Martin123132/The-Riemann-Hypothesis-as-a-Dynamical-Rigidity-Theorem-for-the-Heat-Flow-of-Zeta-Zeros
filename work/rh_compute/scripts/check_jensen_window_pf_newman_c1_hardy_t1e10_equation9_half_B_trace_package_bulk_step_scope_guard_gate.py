#!/usr/bin/env python3
"""Independently check B trace-package and bulk-step ownership."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_trace_package_bulk_step_scope_guard_gate"
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

    bulk, A_plus, U_plus, I_minus, B_minus = sp.symbols("bulk A U Iminus Bminus")
    for tau, local, q_negative in ((0, 0, 0), (0, 1, 1), (0, 1, 0), (1, 1, 1), (1, 1, 0), (1, 0, 1), (0, 0, 1)):
        sigma = 1 - tau - q_negative
        pair = A_plus + I_minus + U_plus + sigma * bulk
        trace = U_plus + B_minus + local * sigma * bulk
        grouped = A_plus + I_minus - B_minus + (1 - local) * sigma * bulk
        require(sp.expand(pair - trace - grouped) == 0, "independent mode allocation failed")

    cases = {
        "m_le_620": (0, 0, 0, 1),
        "m_621_left": (0, 1, 1, 0),
        "m_621_right": (0, 1, 0, 1),
        "m_622_left": (1, 1, 1, -1),
        "m_622_right": (1, 1, 0, 0),
        "m_623_to_39894": (1, 0, 1, -1),
        "m_ge_39895": (0, 0, 1, 0),
    }
    for name, (tau, _local, q_negative, expected_sigma) in cases.items():
        require(1 - tau - q_negative == expected_sigma, f"sigma roster failed: {name}")

    ctx.dps = 70
    ctx.threads = 1
    certificate = artifact["interval_certificate"]
    c_low = arb(certificate["c_window_low_ball"])
    c_high = arb(certificate["c_window_high_ball"])
    require(arb(620) < c_low < arb(621), "independent c_low enclosure failed")
    require(arb(622) < c_high < arb(623), "independent c_high enclosure failed")

    symbolic = artifact["symbolic_certificate"]
    require("sum_(m=1)^620" in symbolic["nonlocal_window_step_block"], "positive omitted block missing")
    require("sum_(m=623)^39894" in symbolic["nonlocal_window_step_block"], "negative omitted block missing")

    decision = artifact["decision"]
    require(decision["negative_trace_package_bound_preserved"] is True, "negative package lost")
    require(decision["nonlocal_bulk_steps_owned_by_grouped_remainder"] is True, "nonlocal ownership lost")
    require(decision["independent_nonlocal_bulk_step_bound_proved"] is False, "step bound overclaimed")
    require(decision["complete_endpoint_B_window_bound_proved"] is False, "endpoint B window overclaimed")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently checked B trace-package scope and nonlocal bulk-step ownership", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
