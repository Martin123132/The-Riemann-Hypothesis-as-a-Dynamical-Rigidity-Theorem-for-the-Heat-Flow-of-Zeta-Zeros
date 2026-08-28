#!/usr/bin/env python3
"""Test pole-free centred-cell recombination of the target-notch dual kernel."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import time


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import numpy as np
from scipy.special import erfc, erfcinv, erfcx


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_pilot"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
TARGET_START = 622
TARGET_END = 39_894
LEFT_EDGE = TARGET_START - 0.5
RIGHT_EDGE = TARGET_END + 0.5
EPSILONS = (1.0e-6, 1.0e-8, 1.0e-10)
CELL_SAMPLES = (0.0, 1.0e-8, 1.0e-5, 1.0e-3, 0.13, 0.49, -0.49)
DUAL_CUTOFF = 64
DIRECT_TAIL_TARGET = 2.0e-11


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def direct_cutoff(epsilon: float) -> tuple[int, float]:
    cutoff = math.ceil(float(erfcinv(DIRECT_TAIL_TARGET * math.sqrt(epsilon))) / math.sqrt(math.pi * epsilon))
    tail = float(erfc(math.sqrt(math.pi * epsilon) * cutoff) / math.sqrt(epsilon))
    require(tail <= DIRECT_TAIL_TARGET, "direct tail cutoff failed")
    return cutoff, tail


def direct_complement(epsilon: float, s: float, cutoff: int) -> complex:
    modes = np.arange(1, cutoff + 1, dtype=np.float64)
    weights = np.exp(-math.pi * epsilon * modes * modes)
    full = 1.0 + 2.0 * float(np.dot(weights, np.cos(2.0 * math.pi * s * modes)))
    target = np.arange(TARGET_START, TARGET_END + 1, dtype=np.float64)
    target_sum = np.dot(
        np.exp(-math.pi * epsilon * target * target),
        np.exp(-2j * math.pi * s * target),
    )
    return complex(full - target_sum)


def edge_data(epsilon: float, edge: float) -> tuple[float, float, float]:
    q = math.exp(-math.pi * epsilon * edge * edge)
    q1 = -2.0 * math.pi * epsilon * edge * q
    q2 = (4.0 * math.pi**2 * epsilon**2 * edge**2 - 2.0 * math.pi * epsilon) * q
    return q, q1, q2


def stable_cell_coefficients(s: float) -> tuple[complex, float]:
    z = math.pi * s
    if abs(z) < 1.0e-3:
        z2 = z * z
        c0 = z * (1.0 / 6.0 + z2 * (7.0 / 360.0 + z2 * (31.0 / 15120.0 + z2 * 127.0 / 604800.0)))
        c1 = 1.0 / 6.0 + z2 * (7.0 / 120.0 + z2 * (31.0 / 3024.0 + z2 * 127.0 / 86400.0))
    else:
        c0 = 1.0 / math.sin(z) - 1.0 / z
        c1 = 1.0 / (z * z) - math.cos(z) / math.sin(z) ** 2
    return c0 / (2j), c1 / 4.0


def h_hat(epsilon: float, xi: float, s: float, k: int) -> complex:
    root_pi_epsilon = math.sqrt(math.pi * epsilon)
    imag = math.sqrt(math.pi / epsilon) * xi
    qc = math.exp(-math.pi * epsilon * LEFT_EDGE**2)
    qd = math.exp(-math.pi * epsilon * RIGHT_EDGE**2)
    sign = -1.0 if k % 2 else 1.0
    phase_c = sign * np.exp(-2j * math.pi * s * LEFT_EDGE)
    phase_d = sign * np.exp(-2j * math.pi * s * RIGHT_EDGE)
    left = qc * phase_c * erfcx(-root_pi_epsilon * LEFT_EDGE - 1j * imag)
    right = qd * phase_d * erfcx(root_pi_epsilon * RIGHT_EDGE + 1j * imag)
    return complex((left + right) / (2.0 * math.sqrt(epsilon)))


def centred_dual(epsilon: float, s: float, cutoff: int) -> tuple[complex, float]:
    qc, q1c, q2c = edge_data(epsilon, LEFT_EDGE)
    qd, q1d, q2d = edge_data(epsilon, RIGHT_EDGE)
    phase_c_s = np.exp(-2j * math.pi * s * LEFT_EDGE)
    phase_d_s = np.exp(-2j * math.pi * s * RIGHT_EDGE)
    e0 = qd * phase_d_s - qc * phase_c_s
    e1 = q1d * phase_d_s - q1c * phase_c_s
    g0, g1 = stable_cell_coefficients(s)
    value = h_hat(epsilon, s, s, 0) + e0 * g0 + e1 * g1

    for k in range(-cutoff, cutoff + 1):
        if k == 0:
            continue
        xi = k + s
        sign = -1.0 if k % 2 else 1.0
        phase_c = sign * phase_c_s
        phase_d = sign * phase_d_s
        denominator = 2j * math.pi * xi
        b0 = (qd * phase_d - qc * phase_c) / denominator
        b1 = (q1d * phase_d - q1c * phase_c) / (denominator * denominator)
        value += h_hat(epsilon, xi, s, k) - b0 - b1

    k3 = abs(q2c) + abs(q2d) + 20.0 * math.pi * epsilon
    tail = k3 / (2.0 * math.pi) ** 3 / (cutoff - 0.5) ** 2
    return complex(value), tail


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    rows = []
    maximum_absolute_discrepancy = 0.0
    maximum_relative_discrepancy = 0.0

    for epsilon in EPSILONS:
        cutoff, direct_tail = direct_cutoff(epsilon)
        for s in CELL_SAMPLES:
            direct = direct_complement(epsilon, s, cutoff)
            dual, dual_tail = centred_dual(epsilon, s, DUAL_CUTOFF)
            discrepancy = abs(direct - dual)
            allowance = direct_tail + dual_tail + 2.0e-8 * max(1.0, abs(direct))
            require(discrepancy <= allowance, f"centred-cell mismatch at epsilon={epsilon}, s={s}")
            maximum_absolute_discrepancy = max(maximum_absolute_discrepancy, discrepancy)
            maximum_relative_discrepancy = max(maximum_relative_discrepancy, discrepancy / max(1.0, abs(direct)))
            rows.append(
                {
                    "epsilon": epsilon,
                    "cell_coordinate": s,
                    "direct_cutoff": cutoff,
                    "dual_cutoff": DUAL_CUTOFF,
                    "direct_tail_upper": direct_tail,
                    "dual_tail_upper": dual_tail,
                    "direct_value": {"real": direct.real, "imag": direct.imag},
                    "dual_value": {"real": dual.real, "imag": dual.imag},
                    "absolute_discrepancy": discrepancy,
                    "relative_discrepancy": discrepancy / max(1.0, abs(direct)),
                }
            )

    artifact = {
        "kind": STEM,
        "status": "exploratory_pole_free_punctured_cell_dual_recombination_passed_not_interval_certified",
        "passed": True,
        "scope": {
            "target_positive_modes": [TARGET_START, TARGET_END],
            "half_integer_edges": [LEFT_EDGE, RIGHT_EDGE],
            "centred_cell": "-1/2<=s<=1/2",
            "epsilons": list(EPSILONS),
            "cell_samples": list(CELL_SAMPLES),
        },
        "formula": {
            "punctured_recombination": "C_epsilon(s)=hat h_epsilon(s)+E_0(s)g_0(s)+E_1(s)g_1(s)+sum_(k!=0)rho_2(k+s)",
            "g0": "g_0(s)=[csc(pi*s)-1/(pi*s)]/(2i), g_0(0)=0",
            "g1": "g_1(s)={1/(pi*s)^2-cos(pi*s)/sin(pi*s)^2}/4, g_1(0)=1/24",
            "uniform_tail": "sum_(|k|>K)|rho_2(k+s)|<=K3/[(2pi)^3(K-1/2)^2] on |s|<=1/2",
        },
        "rows": rows,
        "summary": {
            "row_count": len(rows),
            "integer_rows": sum(row["cell_coordinate"] == 0.0 for row in rows),
            "near_integer_rows": sum(abs(row["cell_coordinate"]) <= 1.0e-5 for row in rows),
            "largest_direct_cutoff": max(row["direct_cutoff"] for row in rows),
            "fixed_dual_cutoff": DUAL_CUTOFF,
            "maximum_absolute_discrepancy": maximum_absolute_discrepancy,
            "maximum_relative_discrepancy": maximum_relative_discrepancy,
            "route_decision": "The k=0 recombination removes the apparent poles and remains stable at and near the integer lattice. Promote the exact centred-cell identity, then fold the u integral to one cell before constructing Kummer-weighted quadrature.",
        },
        "decision": {
            "integer_and_near_integer_rows_passed": True,
            "single_pole_free_dual_cell_route_worth_promoting": True,
            "interval_certified": False,
            "Kummer_source_integration_completed": False,
            "R_after_A_evaluated": False,
            "R_after_A_bound_proved": False,
            "rh_implication": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "blas_threads": 1,
            "process_priority": priority,
        },
        "source": {"path": str(Path(__file__).resolve().relative_to(REPO_ROOT.resolve())).replace("\\", "/"), "sha256": file_hash(Path(__file__).resolve())},
        "proof_boundary": "Floating route-selection pilot for pole-free cell recombination only. It is not interval arithmetic, does not perform the Kummer-source integral or evaluate R_after_A, and proves no R_after_A, R_Dir, Q_K-T, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"passed {len(rows)} punctured-cell rows; max relative discrepancy "
        f"{maximum_relative_discrepancy:.3e}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
