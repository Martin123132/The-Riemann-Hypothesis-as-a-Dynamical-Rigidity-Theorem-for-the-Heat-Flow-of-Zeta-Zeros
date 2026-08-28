#!/usr/bin/env python3
"""Certify the natural endpoint-integral Hardy lift on transition cells."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import time
from pathlib import Path
from typing import Any

import mpmath as mp


ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_cell_natural_A_carrier_lift_gate"
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
NOTE = ROOT / "outputs" / f"{STEM}.md"

HEIGHT = 10_000_000_000
ENDPOINT = 159_577
EPSILON_A = -1
MODES = (39_853, 39_894, 39_936)

COMMON_CELL = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_cell_"
    "common_carrier_subtraction_reduction_gate.json"
)
A_ASSEMBLY = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_"
    "hyperbolic_wedge_A_complete_endpoint_projector_assembly_gate.json"
)
LOCAL_AFFINE = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_"
    "hyperbolic_wedge_affine_A_local_exact_domain_gate.json"
)
LOCAL_AFFINE_ROWS = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_"
    "hyperbolic_wedge_affine_A_local_exact_domain_gate_rows.jsonl"
)
EXACT_H = ROOT / "work" / "rh_compute" / "results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_exact_H_"
    "continuation_defect_target_gate.json"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


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


def ball_midpoint(text: str) -> mp.mpf:
    cleaned = text.strip()
    if cleaned.startswith("[") and cleaned.endswith("]"):
        cleaned = cleaned[1:-1]
    return mp.mpf(cleaned.split(" +/- ", 1)[0])


def complex_midpoint(record: dict[str, str]) -> mp.mpc:
    return mp.mpc(ball_midpoint(record["real_ball"]), ball_midpoint(record["imag_ball"]))


def decimal(value: Any, digits: int = 36) -> str:
    return mp.nstr(value, digits)


def load_selected_rows() -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {}
    for line in LOCAL_AFFINE_ROWS.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        record = json.loads(line)
        row = record.get("row", {})
        mode = row.get("mode")
        if mode in MODES:
            rows[int(mode)] = row
    require(tuple(sorted(rows)) == MODES, "selected local A rows are incomplete")
    return rows


def endpoint_factors(t: mp.mpf, endpoint: int, mode: int) -> dict[str, Any]:
    pi = mp.pi
    m = mp.mpf(mode)
    r = t / (2 * pi * m**2)
    c_t = (pi / (32 * t)) ** mp.mpf("0.25")
    theta_zero = t * (mp.log(t / (2 * pi)) - 1) / 2 - pi / 8
    reduced_phase = theta_zero + pi / 8 - t * mp.log(m)
    full_phase = pi * m * endpoint - pi * m**2 - t / 2 + t * mp.log(t / (2 * pi * m**2)) / 2
    W0 = 4 * (pi * endpoint * m - 1j) * r ** mp.mpf("0.75") / (
        endpoint * mp.sqrt(pi * t)
    )
    return {
        "m": m,
        "r": r,
        "c_t": c_t,
        "theta_zero": theta_zero,
        "reduced_phase": reduced_phase,
        "full_phase": full_phase,
        "W0": W0,
    }


def lift_witness(
    t: mp.mpf,
    endpoint: int,
    mode: int,
    canonical_integral: mp.mpc,
    H: mp.mpf,
    theta: mp.mpf,
) -> dict[str, Any]:
    f = endpoint_factors(t, endpoint, mode)
    m = f["m"]
    s = mp.mpf("0.5") + 1j * t
    raw = EPSILON_A * mp.exp(1j * f["full_phase"]) * f["W0"] * canonical_integral
    X = f["c_t"] * mp.exp(-1j * mp.pi / 8) * raw
    eta = (
        EPSILON_A
        * f["c_t"]
        * mp.sqrt(m)
        * mp.exp(1j * f["theta_zero"])
        * f["W0"]
        * canonical_integral
    )
    m_minus_s = mp.exp(-s * mp.log(m))
    coefficient = H * mp.exp(-1j * theta) * eta
    natural_lift = m_minus_s * coefficient
    physical = 2 * mp.re(X)
    hardy_natural = 2 * mp.re(mp.exp(1j * theta) * natural_lift)
    canonical_lift = H * mp.exp(-1j * theta) * physical / 2
    hardy_null = 2 * mp.re(mp.exp(1j * theta) * (natural_lift - canonical_lift))
    return {
        "mode": mode,
        "phase_parity_discrepancy_absolute": decimal(
            abs(mp.exp(1j * f["full_phase"]) - mp.exp(1j * f["reduced_phase"]))
        ),
        "m_minus_s_factorization_discrepancy_absolute": decimal(abs(X - m_minus_s * eta)),
        "natural_Hardy_lift_discrepancy_absolute": decimal(abs(hardy_natural - H * physical)),
        "natural_minus_canonical_Hardy_null_discrepancy_absolute": decimal(abs(hardy_null)),
        "natural_minus_canonical_lift_absolute": decimal(abs(natural_lift - canonical_lift)),
        "endpoint_coefficient_absolute": decimal(abs(coefficient)),
        "physical_projection": decimal(physical),
        "raw_endpoint_real": decimal(mp.re(raw)),
        "raw_endpoint_imag": decimal(mp.im(raw)),
    }


def production_row_witnesses(H: mp.mpf) -> list[dict[str, Any]]:
    rows = load_selected_rows()
    t = mp.mpf(HEIGHT)
    theta = mp.siegeltheta(t)
    witnesses: list[dict[str, Any]] = []
    for mode in MODES:
        row = rows[mode]
        canonical = complex_midpoint(row["canonical_local_exact_affine_ball"])
        saved_raw = complex_midpoint(row["raw_A_endpoint_local_exact_affine_ball"])
        saved_physical = ball_midpoint(row["physical_local_exact_affine_ball"])
        f = endpoint_factors(t, ENDPOINT, mode)
        reconstructed_raw = EPSILON_A * mp.exp(1j * f["full_phase"]) * f["W0"] * canonical
        reconstructed_physical = 2 * f["c_t"] * mp.re(mp.exp(-1j * mp.pi / 8) * reconstructed_raw)
        witness = lift_witness(t, ENDPOINT, mode, canonical, H, theta)
        witness["saved_raw_reconstruction_discrepancy_absolute"] = decimal(
            abs(reconstructed_raw - saved_raw)
        )
        witness["saved_physical_reconstruction_discrepancy_absolute"] = decimal(
            abs(reconstructed_physical - saved_physical)
        )
        require(abs(reconstructed_raw - saved_raw) < mp.mpf("5e-10"), "saved raw A row drift")
        require(
            abs(reconstructed_physical - saved_physical) < mp.mpf("5e-12"),
            "saved physical A row drift",
        )
        witnesses.append(witness)
    return witnesses


def surrogate_canonical_integral(t: mp.mpf, endpoint: int, mode: int) -> mp.mpc:
    m = mp.mpf(mode)
    lower_s = mp.mpf("-0.42") + m / 100
    upper_s = lower_s + mp.mpf("0.63")
    rho = mp.mpf("-0.73") + m / 200
    face_zero = mp.mpf("0.31") + endpoint / 1000
    lambda_p = (3 * mp.pi * endpoint * m - 1j) / (
        2 * (mp.pi * endpoint * m - 1j) * mp.sqrt(mp.pi * endpoint * m)
    )
    mu_s = mp.mpf(5) / 12 * mp.sqrt(2 / t)

    def outer(S: mp.mpf) -> mp.mpc:
        face = face_zero + rho * S

        def inner(P: mp.mpf) -> mp.mpc:
            amplitude = 1 + lambda_p * P + mu_s * S
            return amplitude * mp.exp(1j * (P**2 - S**2) / 2)

        return mp.quad(inner, [-mp.mpf("0.88"), face])

    return mp.quad(outer, [lower_s, upper_s]) / (2 * mp.pi)


def altered_surrogate_witnesses() -> dict[str, Any]:
    t = mp.mpf("23.75")
    endpoint = 19
    H = mp.mpf("1.03125")
    theta = mp.siegeltheta(t)
    modes = (2, 4, 7)
    rows = [
        lift_witness(t, endpoint, mode, surrogate_canonical_integral(t, endpoint, mode), H, theta)
        for mode in modes
    ]
    tolerance = mp.mpf("1e-65")
    for row in rows:
        require(mp.mpf(row["phase_parity_discrepancy_absolute"]) < tolerance, "altered phase parity failed")
        require(
            mp.mpf(row["m_minus_s_factorization_discrepancy_absolute"]) < tolerance,
            "altered carrier factorization failed",
        )
        require(
            mp.mpf(row["natural_Hardy_lift_discrepancy_absolute"]) < tolerance,
            "altered natural Hardy lift failed",
        )
    return {
        "height": decimal(t),
        "odd_endpoint": endpoint,
        "modes": list(modes),
        "canonical_integral": "finite altered affine wedge integral",
        "rows": rows,
    }


def render_note(artifact: dict[str, Any]) -> str:
    return f"""# Natural A endpoint-integral carrier lift

Date: 2026-08-27

Status: exact transition-cell carrier alignment; no transition coefficient bound.

Let `s=1/2+it`, `c_t=(pi/(32t))^(1/4)`, and let
`C_A,m^exact` denote the convergent exact endpoint-wedge integral from
Formal Core (11.428.1).  For the odd endpoint `A`, endpoint orientation
`epsilon_A=-1`, phase `Phi_A,m`, and prefactor `W0_A,m`, define

```text
X_A,m=c_t exp(-i*pi/8) epsilon_A exp(i Phi_A,m)
      W0_A,m C_A,m^exact,
a_m=2 Re X_A,m.
```

Odd-endpoint parity gives

```text
exp(-i*pi/8)exp(i Phi_A,m)=exp(i theta_0)m^(-it),
theta_0=t[log(t/(2pi))-1]/2-pi/8.
```

Consequently

```text
X_A,m=m^(-s) eta_A,m,
eta_A,m=epsilon_A c_t sqrt(m) exp(i theta_0)
        W0_A,m C_A,m^exact.
```

The natural Hardy lift and its same-cell coefficient are

```text
mathcal A_A,m^nat=H(t)exp(-i theta(t))X_A,m
                 =m^(-s)C_A,m,
C_A,m=H(t)exp(-i theta(t))eta_A,m,
Hardy_t[mathcal A_A,m^nat]=H(t)a_m.
```

It differs from the earlier real canonical lift by a Hardy-null complex
component.  Keeping that component is essential because transition cells now
reduce exactly to the single pre-norm coefficient

```text
m^(-s)[beta_m-C_G-C_A,m], 39853<=m<=39936.
```

The production replay reconstructs saved complex A rows at modes
`{', '.join(str(m) for m in MODES)}`.  A separate altered witness uses height
`{artifact['altered_surrogate']['height']}`, odd endpoint
`{artifact['altered_surrogate']['odd_endpoint']}`, and modes
`{artifact['altered_surrogate']['modes']}`.

Pi provenance: every `pi` above descends from the equation-(9) endpoint phase,
odd-endpoint Fourier parity, the exact bi-Morse normalization, the paper
rotation, or the Riemann-Siegel phase.  No fitted constant is introduced.

Proof boundary: exact natural complex lift and common `m^(-s)` coefficient
alignment only.  This gate does not provide 84 modewise full-endpoint complex
enclosures, bound `beta_m-C_G-C_A,m`, join `P_W` to unowned cells, enclose
`J_Z` or `D_K`, or prove a non-A, all-height, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    mp.mp.dps = 95
    dependencies = {
        "common_cell": COMMON_CELL,
        "A_assembly": A_ASSEMBLY,
        "local_affine": LOCAL_AFFINE,
        "exact_H": EXACT_H,
    }
    loaded = {name: load_json(path) for name, path in dependencies.items()}
    require(all(item.get("passed") is True for item in loaded.values()), "a dependency is not passed")
    require(
        loaded["common_cell"]["decision"]["natural_A_same_cell_integral_lift_derived"] is False,
        "common-cell dependency boundary drift",
    )
    require(
        loaded["A_assembly"]["decision"]["complete_exact_A_endpoint_carrier_certified"] is True,
        "exact A endpoint dependency drift",
    )
    H = ball_midpoint(loaded["exact_H"]["exact_H_certificate"]["high_precision_H_ball"])
    production = production_row_witnesses(H)
    altered = altered_surrogate_witnesses()

    source_builder = Path(__file__).resolve()
    source_checker = ROOT / "work" / "rh_compute" / "scripts" / f"check_{STEM}.py"
    artifact = {
        "kind": STEM,
        "status": "exact_natural_A_endpoint_integral_Hardy_lift_and_transition_coefficient_alignment_certified",
        "passed": True,
        "exact_reduction": {
            "canonical_endpoint_integral": "C_A,m^exact=(2pi)^(-1) integral_(S>h_m) integral_(P<P_A,m(S)) A_A,m(P,S) exp(i(P^2-S^2)/2)dP dS",
            "complex_preprojection": "X_A,m=c_t exp(-i*pi/8) epsilon_A exp(i Phi_A,m) W0_A,m C_A,m^exact",
            "physical_mode": "a_m=2Re[X_A,m]=P_t[A_m+A_-m]",
            "carrier_factorization": "X_A,m=m^(-s)eta_A,m",
            "endpoint_amplitude": "eta_A,m=epsilon_A c_t sqrt(m) exp(i theta_0) W0_A,m C_A,m^exact",
            "natural_Hardy_lift": "mathcal_A_A,m^nat=H exp(-i theta)X_A,m=m^(-s)C_A,m",
            "same_cell_coefficient": "C_A,m=H exp(-i theta)eta_A,m",
            "Hardy_projection": "Hardy_t[mathcal_A_A,m^nat]=H a_m",
            "Hardy_null_difference": "Hardy_t[mathcal_A_A,m^nat-(H/2)exp(-i theta)a_m]=0",
            "transition_cell": "m^(-s)[beta_m-C_G-C_A,m] for 39853<=m<=39936",
        },
        "production_local_integral_witnesses": production,
        "altered_surrogate": altered,
        "decision": {
            "natural_A_endpoint_integral_Hardy_lift_exact": True,
            "A_endpoint_aligned_on_common_m_minus_s_cell_carrier": True,
            "arbitrary_real_canonical_lift_eliminated_from_exact_reduction": True,
            "transition_cells_reduce_to_one_joined_complex_coefficient": True,
            "modewise_full_endpoint_complex_enclosures_certified": False,
            "transition_joined_coefficients_bounded": False,
            "P_W_unowned_cell_Mordell_join_proved": False,
            "actual_height_J_Z_enclosed": False,
            "actual_height_D_K_enclosed": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in dependencies.items()
        },
        "witness_rows": {"path": relative(LOCAL_AFFINE_ROWS), "sha256": file_hash(LOCAL_AFFINE_ROWS)},
        "sources": {
            "builder": {"path": relative(source_builder), "sha256": file_hash(source_builder)},
            "checker": {"path": relative(source_checker), "sha256": file_hash(source_checker)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.time() - started, 3),
        },
        "next_obligation": "Produce rigorous modewise complex enclosures for the complete 84-mode endpoint coefficients C_A,m, then bound beta_m-C_G-C_A,m as joined packets. In parallel, derive the exact Mordell join of P_W with cells 1..621 and 39937..infinity.",
        "proof_boundary": "Exact natural endpoint-integral Hardy lift, common m^(-s) carrier alignment, and transition coefficient reduction only. No modewise full-endpoint complex enclosure, joined coefficient bound, P_W/unowned-cell Mordell join, actual-height J_Z or D_K enclosure, non-A, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified natural A endpoint-integral carrier lift", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
