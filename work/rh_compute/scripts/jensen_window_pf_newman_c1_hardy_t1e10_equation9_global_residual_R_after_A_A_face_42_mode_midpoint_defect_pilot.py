#!/usr/bin/env python3
"""Scout the phase-stripped midpoint defect for the 42 shifted A modes."""

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
from scipy.integrate import quad
from scipy.special import erfcx


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_A_face_42_mode_midpoint_defect_pilot"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"

A = 159_577.0
T = 10_000_000_000.0
LEFT_EDGE = 39_894.5
RIGHT_EDGE = 39_936.5
MODES = np.arange(39_895.0, 39_937.0)
X_STAR = (1.0 - math.sqrt(1.0 - 8.0 * T / (math.pi * A * A))) / 2.0
X_SAMPLES = (0.45, 0.49, 0.499, X_STAR, 0.4999, 0.5)


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


def phase_stripped_A_tail(mode: np.ndarray | float, x: float) -> np.ndarray | complex:
    mode_array = np.asarray(mode, dtype=np.float64)
    delta = mode_array - A * x / 2.0
    q = -delta * math.sqrt(2.0 / x)
    require(bool(np.all(q < 0.0)), "A-tail scout left the negative-q sheet")
    a = mode_array * math.sqrt(2.0 / x)
    w = np.exp(-1j * math.pi / 4.0) * math.sqrt(math.pi / 2.0) * (-q)
    fresnel_tail_over_carrier = np.exp(1j * math.pi / 4.0) / math.sqrt(2.0) * erfcx(w)
    value = -1.0 / (1j * math.pi) - a * fresnel_tail_over_carrier
    if np.ndim(mode) == 0:
        return complex(value)
    return value


def two_current_edge(mode: float, x: float) -> tuple[complex, float]:
    delta = mode - A * x / 2.0
    c0 = A * x / (2j * math.pi * delta)
    c1 = -mode * x / (2.0 * math.pi**2 * delta**3)
    remainder = 3.0 * mode * x * x / (2.0 * math.pi**3 * delta**5)
    return complex(c0 + c1), remainder


def integrate_strip(x: float) -> tuple[complex, float, float]:
    def real_part(mode: float) -> float:
        return phase_stripped_A_tail(mode, x).real

    def imag_part(mode: float) -> float:
        return phase_stripped_A_tail(mode, x).imag

    real_value, real_error = quad(real_part, LEFT_EDGE, RIGHT_EDGE, epsabs=2.0e-9, epsrel=2.0e-12, limit=300)
    imag_value, imag_error = quad(imag_part, LEFT_EDGE, RIGHT_EDGE, epsabs=2.0e-9, epsrel=2.0e-12, limit=300)
    return complex(real_value, imag_value), real_error, imag_error


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    rows = []
    maximum_two_current_ratio = 0.0
    roundoff_limited_two_current_rows = 0
    for x in X_SAMPLES:
        exact_modes = phase_stripped_A_tail(MODES, x)
        midpoint_sum = complex(np.sum(exact_modes))
        strip_integral, real_error, imag_error = integrate_strip(x)
        defect = midpoint_sum - strip_integral

        edge_exact = phase_stripped_A_tail(RIGHT_EDGE, x)
        edge_two, edge_remainder = two_current_edge(RIGHT_EDGE, x)
        edge_error = abs(edge_exact - edge_two)
        require(edge_error <= edge_remainder * (1.0 + 2.0e-8) + 2.0e-11, f"two-current edge mismatch at x={x}")
        maximum_two_current_ratio = max(maximum_two_current_ratio, edge_error / edge_remainder)
        if edge_remainder < 2.0e-11:
            roundoff_limited_two_current_rows += 1

        rows.append(
            {
                "x": x,
                "is_tangential_saddle": abs(x - X_STAR) < 1.0e-14,
                "midpoint_sum": {"real": midpoint_sum.real, "imag": midpoint_sum.imag, "modulus": abs(midpoint_sum)},
                "strip_integral": {"real": strip_integral.real, "imag": strip_integral.imag, "modulus": abs(strip_integral)},
                "quadrature_error_estimate": {"real": real_error, "imag": imag_error},
                "midpoint_defect": {"real": defect.real, "imag": defect.imag, "modulus": abs(defect)},
                "relative_defect": abs(defect) / max(1.0, abs(midpoint_sum)),
                "translated_edge_exact": {"real": edge_exact.real, "imag": edge_exact.imag, "modulus": abs(edge_exact)},
                "translated_edge_two_current": {"real": edge_two.real, "imag": edge_two.imag, "modulus": abs(edge_two)},
                "translated_edge_two_current_error": edge_error,
                "translated_edge_remainder_bound": edge_remainder,
            }
        )

    saddle_row = next(row for row in rows if row["is_tangential_saddle"])
    artifact = {
        "kind": STEM,
        "status": "exploratory_phase_stripped_42_mode_midpoint_defect_and_translated_edge_two_current_crosscheck_not_interval_certified",
        "passed": True,
        "scope": {
            "height": int(T),
            "A": int(A),
            "shifted_modes": [39895, 39936],
            "continuous_strip": [LEFT_EDGE, RIGHT_EDGE],
            "x_samples": list(X_SAMPLES),
            "x_star": X_STAR,
            "Abel_regulator": 0.0,
        },
        "formula": {
            "phase_stripped_tail": "G_A(y,x)=exp(-i*pi*q_A(y,x)^2/2)L_A(y,x)=-1/(i*pi)-y*sqrt(2/x)*exp(i*pi/4)erfcx[e^(-i*pi/4)sqrt(pi/2)(-q_A)]/sqrt(2)",
            "parity_cancellation": "The natural continuous factor exp(i*pi*A*y) cancels exp(-i*pi*A*y) from quadratic completion; at integer y this is the inherited odd-A parity.",
            "midpoint_defect": "D_42(x)=sum_(m=39895)^39936 G_A(m,x)-integral_(39894.5)^39936.5 G_A(y,x)dy",
            "interpretation": "D_42 is the k!=0 midpoint-Poisson content of the shifted strip. It is diagnostic until the full fixed-regulator edge/block identity and Kummer amplitude are reassembled.",
        },
        "rows": rows,
        "summary": {
            "row_count": len(rows),
            "maximum_midpoint_defect_modulus": max(row["midpoint_defect"]["modulus"] for row in rows),
            "minimum_midpoint_defect_modulus": min(row["midpoint_defect"]["modulus"] for row in rows),
            "saddle_midpoint_defect": saddle_row["midpoint_defect"],
            "maximum_two_current_error_to_bound_ratio": maximum_two_current_ratio,
            "roundoff_limited_two_current_rows": roundoff_limited_two_current_rows,
            "two_current_ratio_warning": "Ratios on rows whose rigorous remainder is below 2e-11 are dominated by double-precision subtraction and are telemetry, not a validation of the interval bound.",
            "route_decision": "The phase-stripped block and strip are individually large. Use their signed midpoint defect, not either modulus, as the next A-face amplitude. Derive its exact Peano/Poisson representation and derivative envelope before the tangential Morse integral.",
        },
        "decision": {
            "two_current_crosschecks_consistent_with_certified_formula_and_float_roundoff_guard": True,
            "midpoint_defect_route_worth_deriving_exactly": True,
            "floating_values_used_as_proof": False,
            "interval_certified": False,
            "fixed_regulator_reassembly_proved": False,
            "tangential_Morse_integral_bounded": False,
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
        "proof_boundary": "Floating phase-stripped midpoint-defect route selection only. No interval enclosure, fixed-regulator edge/block theorem, Peano derivative bound, tangential Morse estimate, R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"scouted {len(rows)} A-face midpoint rows; saddle defect "
        f"{saddle_row['midpoint_defect']['modulus']:.6e}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
