#!/usr/bin/env python3
"""Independent check of the uniform Gamma-target insertion gate."""

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

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_uniform_gamma_target_insertion_gate"
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

    # Rebuild the finite weighted identity with an independently ordered set
    # of symbols.
    source_kernel, target_integral, bulk_sum, tail_sum, half = sp.symbols(
        "D I B Q H"
    )
    lhs = half + source_kernel - target_integral + tail_sum
    rhs = half + source_kernel - bulk_sum
    require(sp.expand(lhs.subs(target_integral, bulk_sum + tail_sum) - rhs) == 0, "independent insertion identity failed")

    qk, gamma, target = sp.symbols("QK Gamma Target")
    require(sp.expand((qk - target) - ((qk - gamma) + (gamma - target))) == 0, "independent target correction failed")

    require(39_694 - 622 + 1 == 39_073, "ordinary count drift")
    require(39_894 - 39_695 + 1 == 200, "fold count drift")
    require(39_073 + 200 == 39_273, "total target count drift")

    # Recompute the Gamma-to-classical correction instead of trusting the
    # production artifact's saved interval.
    ctx.dps = 90
    ctx.threads = 1
    t = arb(10_000_000_000)
    pi = arb.pi()
    imaginary_unit = acb(0, 1)
    y = t / 2
    exponent = acb(arb("0.75"), y)
    log_integral = exponent.lgamma() + imaginary_unit * y - exponent * (
        y.log() + imaginary_unit * pi / 2
    )
    log_gaussian = arb(2).log() + (pi.log() - t.log()) / 2 - imaginary_unit * pi / 4
    factor = (log_integral - log_gaussian).exp()
    theta_exact = acb(arb("0.25"), y).lgamma().imag - y * pi.log()
    theta_zero = y * ((t / (2 * pi)).log() - 1) - pi / 8
    target_factor = (imaginary_unit * (theta_exact - theta_zero)).exp()
    reciprocal_sqrt_sum = arb(0)
    for mode in range(622, 39_895):
        reciprocal_sqrt_sum += arb(mode).rsqrt()
    physical = arb(2).sqrt() / pi * abs(factor - target_factor) * reciprocal_sqrt_sum
    require(physical.upper() < arb("5e-20"), "independent physical correction exceeds 5e-20")
    stored = arb(artifact["correction_certificate"]["Gamma_to_classical_physical_absolute_bound_ball"])
    require(physical.overlaps(stored), "independent correction misses stored interval")

    getcontext().prec = 50
    require(
        Decimal("1.405792e-4") - Decimal("5e-20") == Decimal("0.00014057919999999995"),
        "independent working target drift",
    )
    require(
        Decimal("1.319792e-4") - Decimal("5e-20") == Decimal("0.00013197919999999995"),
        "independent negative target drift",
    )

    decision = artifact["decision"]
    require(decision["Dirichlet_target_kernel_is_exact_Gamma_bulk_insertion"] is True, "carrier insertion decision drift")
    require(decision["A_face_exception_remains_in_endpoint_defect_not_bulk_carrier"] is True, "A-face ownership drift")
    require(decision["compressed_R_Dir_target_proved"] is False, "compressed residual overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently checked uniform Gamma target insertion and correction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
