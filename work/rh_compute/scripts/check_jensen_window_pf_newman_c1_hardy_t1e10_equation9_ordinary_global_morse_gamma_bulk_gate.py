#!/usr/bin/env python3
"""Independent replay of the exact global-Morse Gamma bulk gate."""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_gamma_bulk_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
CHECKER = Path(__file__).resolve()
T = 10_000_000_000
LOWER_START = 622
UPPER_END = 39_894


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stored_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def overlap_complex(left: acb, right: acb) -> bool:
    return left.real.overlaps(right.real) and left.imag.overlaps(right.imag)


def main() -> int:
    ctx.dps = 140
    ctx.threads = 1
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "Gamma bulk artifact is not passed")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(CHECKER), "checker hash drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(dependency["sha256"] == file_hash(path), f"dependency hash drift: {path}")

    t = arb(T)
    y = t / 2
    pi = arb.pi()
    imaginary_unit = acb(0, 1)
    exponent = acb(arb("0.75"), y)
    log_integral = exponent.lgamma() + imaginary_unit * y - exponent * (y.log() + imaginary_unit * pi / 2)
    log_gaussian = arb(2).log() + (pi.log() - t.log()) / 2 - imaginary_unit * pi / 4
    log_factor = log_integral - log_gaussian
    factor = log_factor.exp()
    theta_exact = acb(arb("0.25"), y).lgamma().imag - y * pi.log()
    theta_zero = y * ((t / (2 * pi)).log() - 1) - pi / 8
    phase_factor = (imaginary_unit * (theta_exact - theta_zero)).exp()
    defect = factor - phase_factor

    reciprocal_sqrt_sum = arb(0)
    for mode in range(UPPER_END, LOWER_START - 1, -1):
        reciprocal_sqrt_sum += arb(mode).rsqrt()
    complex_bound = abs(defect) * reciprocal_sqrt_sum
    two_real_bound = 2 * complex_bound
    physical_bound = arb(2).sqrt() / pi * complex_bound

    certificate = artifact["certificate"]
    require(overlap_complex(log_factor, stored_complex(certificate["log_Gamma_bulk_factor_ball"])), "log-factor replay mismatch")
    require(overlap_complex(factor, stored_complex(certificate["Gamma_bulk_factor_ball"])), "factor replay mismatch")
    require(overlap_complex(defect, stored_complex(certificate["bulk_coefficient_defect_ball"])), "defect replay mismatch")
    require(reciprocal_sqrt_sum.overlaps(arb(certificate["reciprocal_sqrt_mode_sum_ball"])), "mode-sum replay mismatch")
    require(complex_bound.overlaps(arb(certificate["aggregate_complex_absolute_bound_ball"])), "complex-bound replay mismatch")
    require(two_real_bound.overlaps(arb(certificate["aggregate_two_real_absolute_bound_ball"])), "two-real replay mismatch")
    require(physical_bound.overlaps(arb(certificate["aggregate_physical_absolute_bound_ball"])), "physical replay mismatch")
    require(abs(defect) < arb("3.2e-22"), "coefficient replay exceeds 3.2e-22")
    require(two_real_bound < arb("3e-19"), "aggregate replay exceeds 3e-19")
    require(artifact["decision"]["finite_lower_endpoint_current_closed_here"] is False, "lower-endpoint overclaim")
    require(artifact["decision"]["finite_upper_endpoint_current_closed_here"] is False, "upper-endpoint overclaim")
    require(artifact["decision"]["complete_T_upper_proved"] is False, "T_upper overclaim")

    print("checked exact global-Morse Gamma bulk over 39273 modes at 140 digits", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
