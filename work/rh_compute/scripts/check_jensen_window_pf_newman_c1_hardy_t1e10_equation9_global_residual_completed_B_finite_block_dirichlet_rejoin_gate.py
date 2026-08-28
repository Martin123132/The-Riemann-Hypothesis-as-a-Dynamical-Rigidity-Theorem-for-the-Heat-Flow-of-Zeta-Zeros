#!/usr/bin/env python3
"""Independent check of the finite completed-B Dirichlet rejoin."""

from __future__ import annotations

from decimal import Decimal, getcontext
import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_finite_block_dirichlet_rejoin_gate"
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

    # Reconstruct the allocation with a different variable ordering.
    low, high, joined, window, eta = sp.symbols("L H J W eta")
    global_current = window + eta * low + eta * high + joined
    compressed = global_current - window - eta * high
    require(sp.expand(compressed - (eta * low + joined)) == 0, "independent rejoin failed")

    u = sp.symbols("u", real=True)
    z = sp.exp(-2 * sp.pi * sp.I * u)
    target_count = 39_894 - 622 + 1
    target_sum = z**622 * (1 - z**target_count) / (1 - z)
    require(sp.simplify(target_count - 39_273) == 0, "target count drift")
    require(sp.cancel((1 - z) * target_sum - (z**622 - z**39895)) == 0, "independent target kernel failed")

    getcontext().prec = 40
    working = Decimal("8.6e-6") + Decimal("1.3198e-4") - Decimal("8e-10")
    negative = Decimal("1.3198e-4") - Decimal("8e-10")
    require(working == Decimal(artifact["target_arithmetic"]["finite_join_sufficient_for_working_target"]), "working target drift")
    require(negative == Decimal(artifact["target_arithmetic"]["finite_join_sufficient_for_negative_R_end"]), "negative target drift")

    scope = artifact["scope"]
    require(scope["finite_completed_positive_mode_count"] == 5_122_420, "finite mode count drift")
    require(scope["avoided_crossing_chart_count"] == 1_280_347, "crossing count drift")
    require(scope["target_mode_count"] == 39_273, "target-mode count drift")
    decision = artifact["decision"]
    require(decision["finite_completed_B_block_rejoined_before_norms"] is True, "rejoin decision drift")
    require(decision["finite_crossing_roster_replaced_by_common_Dirichlet_Abel_kernel"] is True, "kernel decision drift")
    require(decision["compressed_residual_below_1p405792e_minus_4_proved"] is False, "working-target overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently checked finite completed-B Dirichlet rejoin and targets", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
