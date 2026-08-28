#!/usr/bin/env python3
"""Scout a two-edge modular evaluator for the target-notched Abel kernel."""

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
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_dual_pilot"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
TARGET_START = 622
TARGET_END = 39_894
LEFT_EDGE = TARGET_START - 0.5
RIGHT_EDGE = TARGET_END + 0.5
EPSILONS = (1.0e-6, 1.0e-8, 1.0e-10)
FRACTIONAL_U = (0.13, 0.29, 0.47)
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
    argument = DIRECT_TAIL_TARGET * math.sqrt(epsilon)
    cutoff = math.ceil(float(erfcinv(argument)) / math.sqrt(math.pi * epsilon))
    tail = float(erfc(math.sqrt(math.pi * epsilon) * cutoff) / math.sqrt(epsilon))
    require(tail <= DIRECT_TAIL_TARGET, "direct Gaussian tail cutoff failed")
    return cutoff, tail


def direct_complement(epsilon: float, u: float, cutoff: int) -> complex:
    modes = np.arange(1, cutoff + 1, dtype=np.float64)
    weights = np.exp(-math.pi * epsilon * modes * modes)
    full_theta = 1.0 + 2.0 * float(np.dot(weights, np.cos(2.0 * math.pi * u * modes)))

    target = np.arange(TARGET_START, TARGET_END + 1, dtype=np.float64)
    target_weights = np.exp(-math.pi * epsilon * target * target)
    target_sum = np.dot(target_weights, np.exp(-2j * math.pi * u * target))
    return complex(full_theta - target_sum)


def edge_data(epsilon: float, edge: float) -> tuple[float, float, float]:
    q = math.exp(-math.pi * epsilon * edge * edge)
    q1 = -2.0 * math.pi * epsilon * edge * q
    q2 = (4.0 * math.pi**2 * epsilon**2 * edge**2 - 2.0 * math.pi * epsilon) * q
    return q, q1, q2


def dual_complement(epsilon: float, u: float, cutoff: int) -> tuple[complex, float]:
    sqrt_epsilon = math.sqrt(epsilon)
    root_pi_epsilon = math.sqrt(math.pi * epsilon)
    root_pi_over_epsilon = math.sqrt(math.pi / epsilon)
    qc, q1c, q2c = edge_data(epsilon, LEFT_EDGE)
    qd, q1d, q2d = edge_data(epsilon, RIGHT_EDGE)
    phase_c_u = np.exp(-2j * math.pi * u * LEFT_EDGE)
    phase_d_u = np.exp(-2j * math.pi * u * RIGHT_EDGE)

    sine = math.sin(math.pi * u)
    cosine = math.cos(math.pi * u)
    require(abs(sine) > 0.05, "pilot point too close to an integer")
    edge0 = (qd * phase_d_u - qc * phase_c_u) / (2j * sine)
    edge1 = -(q1d * phase_d_u - q1c * phase_c_u) * cosine / (4.0 * sine * sine)

    residual = 0j
    for k in range(-cutoff, cutoff + 1):
        xi = k + u
        sign = -1.0 if k % 2 else 1.0
        phase_c = sign * phase_c_u
        phase_d = sign * phase_d_u
        imag = root_pi_over_epsilon * xi
        left = qc * phase_c * erfcx(-root_pi_epsilon * LEFT_EDGE - 1j * imag)
        right = qd * phase_d * erfcx(root_pi_epsilon * RIGHT_EDGE + 1j * imag)
        h_hat = (left + right) / (2.0 * sqrt_epsilon)
        denominator = 2j * math.pi * xi
        boundary0 = (qd * phase_d - qc * phase_c) / denominator
        boundary1 = (q1d * phase_d - q1c * phase_c) / (denominator * denominator)
        residual += h_hat - boundary0 - boundary1

    q3_l1_upper = 20.0 * math.pi * epsilon
    k3 = abs(q2c) + abs(q2d) + q3_l1_upper
    reciprocal_tail = 0.5 / (cutoff + u) ** 2 + 0.5 / (cutoff - u) ** 2
    tail_bound = k3 * reciprocal_tail / (2.0 * math.pi) ** 3
    return complex(edge0 + edge1 + residual), tail_bound


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    rows = []
    maximum_discrepancy = 0.0
    maximum_normalized_discrepancy = 0.0

    for epsilon in EPSILONS:
        cutoff, direct_tail = direct_cutoff(epsilon)
        for u in FRACTIONAL_U:
            direct = direct_complement(epsilon, u, cutoff)
            dual, dual_tail = dual_complement(epsilon, u, DUAL_CUTOFF)
            discrepancy = abs(direct - dual)
            floating_allowance = 3.0e-8
            comparison_allowance = direct_tail + dual_tail + floating_allowance
            require(discrepancy <= comparison_allowance, f"direct/dual mismatch at epsilon={epsilon}, u={u}")
            maximum_discrepancy = max(maximum_discrepancy, discrepancy)
            maximum_normalized_discrepancy = max(maximum_normalized_discrepancy, discrepancy / comparison_allowance)
            rows.append(
                {
                    "epsilon": epsilon,
                    "u_fraction": u,
                    "direct_cutoff": cutoff,
                    "direct_tail_upper": direct_tail,
                    "dual_cutoff": DUAL_CUTOFF,
                    "dual_tail_upper": dual_tail,
                    "direct_value": {"real": direct.real, "imag": direct.imag},
                    "dual_value": {"real": dual.real, "imag": dual.imag},
                    "absolute_discrepancy": discrepancy,
                    "comparison_allowance": comparison_allowance,
                }
            )

    artifact = {
        "kind": STEM,
        "status": "exploratory_two_edge_modular_target_notch_pilot_passed_not_interval_certified",
        "passed": True,
        "scope": {
            "target_positive_modes": [TARGET_START, TARGET_END],
            "half_integer_edges": [LEFT_EDGE, RIGHT_EDGE],
            "epsilons": list(EPSILONS),
            "fractional_u": list(FRACTIONAL_U),
        },
        "formula": {
            "sampled_complement": "C_epsilon(u)=sum_(m not in T)exp(-pi*epsilon*m^2)exp(-2pi*i*m*u)",
            "continuous_complement": "h_epsilon(y)=exp(-pi*epsilon*y^2)[1_(y<c)+1_(y>d)]",
            "poisson_dual": "C_epsilon(u)=sum_(k in Z)hat h_epsilon(k+u)",
            "boundary_terms": "hat h=b_0+b_1+r_2, with b_j=[q^(j)(d)e^(-2pi*i*xi*d)-q^(j)(c)e^(-2pi*i*xi*c)]/(2pi*i*xi)^(j+1)",
            "closed_boundary_sums": "sum b_0=E_0/(2i sin(pi*u)); sum b_1=-E_1 cos(pi*u)/(4 sin(pi*u)^2)",
            "dual_remainder_bound": "|r_2(xi)|<=[|q''(c)|+|q''(d)|+20pi*epsilon]/(2pi|xi|)^3",
        },
        "rows": rows,
        "summary": {
            "row_count": len(rows),
            "maximum_absolute_discrepancy": maximum_discrepancy,
            "maximum_allowance_fraction": maximum_normalized_discrepancy,
            "largest_direct_cutoff": max(row["direct_cutoff"] for row in rows),
            "fixed_dual_cutoff": DUAL_CUTOFF,
            "route_decision": "The two-edge boundary-subtracted modular form is numerically consistent with the direct target-notched kernel and has a cubic dual tail away from integer u. Promote the exact identity and a hybrid near-integer/direct versus away-from-integer/dual evaluator; do not infer an R_after_A bound from these floating samples.",
        },
        "decision": {
            "direct_and_dual_samples_agree_within_declared_allowance": True,
            "two_edge_modular_route_worth_promoting": True,
            "interval_certified": False,
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
        "proof_boundary": "Floating route-selection pilot for the regulated target-notched kernel only. It is not interval arithmetic, does not evaluate the joined R_after_A object, and proves no R_after_A, R_Dir, Q_K-T, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"passed {len(rows)} target-notch direct/dual rows; "
        f"largest direct cutoff {artifact['summary']['largest_direct_cutoff']}, "
        f"max discrepancy {maximum_discrepancy:.3e}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
