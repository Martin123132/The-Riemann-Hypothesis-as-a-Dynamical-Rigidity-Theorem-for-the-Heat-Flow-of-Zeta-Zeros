#!/usr/bin/env python3
"""Independently check the grouped characteristic Airy-boundary gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))

from flint import arb, ctx


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate.md"

PRECISION = 120
T = 10_000_000_000
C = 159577
B = 5122421
MODE_LO = 39853
MODE_HI = 39936
Z_CORE = 4
Y_CORE = 64


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    require(RESULT.is_file() and NOTE.is_file(), "missing result or note")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "gate is not marked passed")
    require(artifact.get("status") == "exact_grouped_characteristic_Airy_boundary_normal_form_and_compact_core_complete", "status drift")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"source hash drift: {path}")

    ctx.dps = PRECISION
    t = arb(T)
    c = arb(C)
    pi = arb.pi()
    eta = pi * c**2 / (8 * t)
    beta = (t * eta) ** (arb(1) / 3)
    sigma = 4 * beta / (pi * c)
    h = 4 * beta / c
    airy_lambda = (eta - 1) * t / beta
    d_min = h * (arb(MODE_LO) - c / 4)
    d_max = h * (arb(MODE_HI) - c / 4)

    endpoint = arb(64) / 15 * arb(Z_CORE) ** 5 / beta**2
    coupling = arb(Y_CORE) * arb(Z_CORE) ** 3 / (3 * beta**2)
    quadratic = arb(Z_CORE) * arb(Y_CORE) ** 2 / (4 * beta**2)
    total = endpoint + coupling + quadratic
    amp_upper = sigma * arb(Y_CORE) / c
    amp_lower = 1 - (arb(Z_CORE) / beta).cosh() ** (-arb(3) / 2)
    amp_error = max(amp_upper, amp_lower)
    first_zero = 2 * pi / (arb(MODE_HI - MODE_LO + 1) * h)
    period = 2 * pi / h
    y_max = arb(B - C) / sigma
    fresnel_scale = beta.sqrt()

    core = artifact["certified_core"]
    checks = {
        "eta_ball": eta,
        "beta_ball": beta,
        "sigma_ball": sigma,
        "Airy_lambda_ball": airy_lambda,
        "mode_detuning_spacing_ball": h,
        "mode_detuning_min_ball": d_min,
        "mode_detuning_max_ball": d_max,
        "Dirichlet_first_zero_ball": first_zero,
        "Dirichlet_period_ball": period,
        "full_y_span_ball": y_max,
        "Fresnel_y_scale_sqrt_beta_ball": fresnel_scale,
        "y_core_over_Fresnel_scale_ball": arb(Y_CORE) / fresnel_scale,
        "endpoint_phase_remainder_bound_ball": endpoint,
        "normal_coupling_remainder_bound_ball": coupling,
        "quadratic_coefficient_remainder_bound_ball": quadratic,
        "total_joint_phase_remainder_bound_ball": total,
        "relative_amplitude_upper_error_bound_ball": amp_upper,
        "relative_amplitude_lower_error_bound_ball": amp_lower,
        "relative_amplitude_error_bound_ball": amp_error,
    }
    for key, value in checks.items():
        require(arb(core[key]).contains(value), f"independent value escaped recorded ball: {key}")

    require(core["mode_count"] == 84, "mode count drift")
    require(total < arb("0.0022"), "phase threshold failed")
    require(amp_error < arb("0.00001"), "amplitude threshold failed")
    require(arb("-2.3") < d_min < arb("-2.2") < arb(0) < arb("2.2") < d_max < arb("2.3"), "detuning atlas failed")
    require(artifact["decision"]["scalar_slow_Airy_amplitude_model_admissible"] is False, "scalar-amplitude guard drift")
    require(artifact["decision"]["canonical_integral_bound_proved"] is False, "proof boundary drift")

    note = NOTE.read_text(encoding="utf-8")
    for token in ("D_84(y)", "Theta_m^0(z,y)", "< 0.0022", "< 1e-5", "not admissible"):
        require(token in note, f"note token missing: {token}")

    print(
        "validated grouped characteristic Airy-boundary form independently: "
        "modes=84, phase<0.0022, amplitude<1e-5, scalar model rejected"
    )


if __name__ == "__main__":
    main()
