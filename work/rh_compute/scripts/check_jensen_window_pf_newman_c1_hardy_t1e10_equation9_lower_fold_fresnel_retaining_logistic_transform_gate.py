#!/usr/bin/env python3
"""Independently validate the exact Fresnel-retaining logistic transform."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, REPO_ROOT / "work/rh_compute/scripts"):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_fresnel_retaining_logistic_transform_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_fresnel_retaining_logistic_transform_gate.md"
PRECISION = 120
C = 159_577
B = 5_122_423
M = 39_894


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


def main() -> None:
    priority = set_low_priority()
    require(RESULT.is_file() and NOTE.is_file(), "missing logistic-transform artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "logistic-transform artifact is not passed")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file() and file_hash(path) == record["sha256"], f"hash drift: {path}")

    x, m, endpoint = sp.symbols("x m D", positive=True)
    q = sp.sqrt(x / 2) * (endpoint - 2 * m / x)
    require(
        sp.simplify(-sp.pi * m**2 / x + sp.pi * q**2 / 2 - (sp.pi * x * endpoint**2 / 4 - sp.pi * m * endpoint)) == 0,
        "independent quadratic completion failed",
    )
    t, v = sp.symbols("t v", positive=True)
    r = t / (2 * sp.pi * m**2)
    psi_difference = -sp.pi * m**2 * r * (v - 1) + t * sp.log(v) / 2
    require(sp.simplify(psi_difference + t * (v - 1 - sp.log(v)) / 2) == 0, "independent logistic phase failed")
    tau = sp.pi * m * (endpoint - 2 * m)
    x_event = sp.simplify(2 * sp.pi * m**2 / (tau + 2 * sp.pi * m**2))
    require(sp.simplify(x_event - 2 * m / endpoint) == 0, "independent event x failed")
    require(sp.simplify(q.subs(x, x_event)) == 0, "independent event q failed")

    ctx.dps = PRECISION
    pi = arb.pi()
    tau_a = pi * M * (C - 2 * M)
    r_a = tau_a / (2 * pi * M**2)
    x_a = 1 / (1 + r_a)
    q_c = (x_a / 2).sqrt() * (C - 2 * M / x_a)
    q_b = (x_a / 2).sqrt() * (B - 2 * M / x_a)
    jacobian = r_a * x_a * (x_a * (1 - x_a)) ** (-arb(1) / 4)
    Delta = (2 * pi * M**2 - tau_a) / (2 * pi * M**2 + tau_a)
    require(q_c.contains(0), "independent q_C event regularity failed")
    require(q_b > arb("2480000"), "independent q_B margin failed")
    require((Delta + arb(1) / C).contains(0), "independent Delta identity failed")
    require(arb("0.70") < jacobian < arb("0.72"), "independent Jacobian regularity failed")

    saved = artifact["prototype_certificate"]
    require(x_a.overlaps(arb(saved["saddle_x_ball"])), "saved saddle x drift")
    require(q_c.overlaps(arb(saved["lower_endpoint_Fresnel_coordinate_ball"])), "saved q_C drift")
    require(q_b.overlaps(arb(saved["upper_endpoint_Fresnel_coordinate_ball"])), "saved q_B drift")
    require(jacobian.overlaps(arb(saved["finite_Jacobian_weight_ball"])), "saved Jacobian drift")

    decisions = artifact["decision"]
    require(decisions["endpoint_current_and_Fresnel_term_remain_paired"] is True, "paired-current decision missing")
    require(decisions["global_logistic_Morse_phase_is_exact"] is True, "exact-Morse decision missing")
    require(decisions["transform_is_regular_at_turning_event"] is True, "event regularity missing")
    require(decisions["division_by_characteristic_defect_used"] is False, "characteristic division introduced")
    require(decisions["Airy_to_logistic_remainder_proved"] is False, "Airy/logistic remainder overpromoted")
    require(decisions["propagation_across_399_events_proved"] is False, "atlas propagation overpromoted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")

    note = NOTE.read_text(encoding="utf-8")
    for token in (
        "P_m(x)=",
        "psi_m(x)-psi_m(x_m)=-t*s^2/4",
        "q_C(m,x_m)=0",
        "does not bound the",
    ):
        require(token in note, f"note token missing: {token}")
    print(
        "validated exact Fresnel-retaining logistic transform independently: "
        f"mode={M}, q_C(event)=0, no characteristic division; priority={priority}"
    )


if __name__ == "__main__":
    main()
