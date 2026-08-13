#!/usr/bin/env python3
"""Float64 FFT scout for the two exact Fourier packets in D_fold."""

from __future__ import annotations

import json
import os
from pathlib import Path
import time

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import numpy as np
from scipy.special import airy


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_increment_fft_scout"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
C = 159_577
Q = 39_894
GRIDS = (262_144, 524_288)
HEIGHT_NODES = 5


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


def beta4_profile(s: np.ndarray, lam: float, beta: float, y_max: float) -> np.ndarray:
    y = y_max * s
    x = lam + y
    ai, aip, _, _ = airy(-x)
    a_x = -aip
    epsilon = beta**-2
    u1 = -(13 * lam + 3 * y) / 60
    v1 = (8 * lam**2 - 4 * lam * y + 3 * y**2) / 60
    n_u2 = 448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2 - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2
    n_v2 = 40 * lam**3 - 20 * lam**2 * y + lam * y**2 - 9 * y**3 + 27
    u = 1 + epsilon * u1 - epsilon**2 * n_u2 / 50400
    v = epsilon * v1 + epsilon**2 * n_v2 / 1680
    phase = -1.5 * np.pi * s + y_max**2 * s**2 / (4 * beta)
    return y_max * np.exp(1j * phase) * (u * ai + v * a_x)


def packet_sum(profile: np.ndarray, frequencies: np.ndarray) -> complex:
    panels = profile.size
    transform = np.fft.fft(profile) / panels
    indices = np.mod(frequencies, panels)
    midpoint_phase = np.exp(-1j * np.pi * frequencies / panels)
    return complex(np.sum(transform[indices] * midpoint_phase))


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    tstar = np.pi * C**2 / 8
    beta = tstar ** (1 / 3)
    mu_radius = np.pi / (16 * tstar)
    lambda_max = beta**2 * mu_radius
    y_max = np.pi * C / (2 * beta)
    low_frequencies = np.arange(-200, 0, dtype=np.int64)
    high_frequencies = np.arange(-79_789, -79_589, dtype=np.int64)
    rows = []
    previous: dict[int, complex] = {}
    for panels in GRIDS:
        for index in range(HEIGHT_NODES):
            lam = lambda_max * index / (HEIGHT_NODES - 1)
            s = (np.arange(panels, dtype=np.float64) + 0.5) / panels
            profile = beta4_profile(s, lam, beta, y_max)
            low = packet_sum(profile, low_frequencies)
            high = packet_sum(profile, high_frequencies)
            total = low + high
            physical = 4 * np.sqrt(2) * total.real
            key = index
            prior = previous.get(key)
            rows.append({
                "panels": panels,
                "height_fraction_index_over_4": index,
                "lambda": lam,
                "low_packet": {"real": low.real, "imag": low.imag, "abs": abs(low)},
                "high_packet": {"real": high.real, "imag": high.imag, "abs": abs(high)},
                "total_packet": {"real": total.real, "imag": total.imag, "abs": abs(total)},
                "physical_fold_increment_beta4": physical,
                "difference_from_previous_grid": None if prior is None else abs(total - prior),
            })
            previous[key] = total

    fine_rows = [row for row in rows if row["panels"] == GRIDS[-1]]
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "nonrigorous_float64_two_packet_fft_scout",
        "passed": True,
        "geometry": {
            "profile_frequency_index": "n=m-39895",
            "positive_fold_packet": [-200, -1],
            "negative_fold_packet": [-79_789, -79_590],
            "packet_size_each": 200,
            "identity": "D_fold=kappa_C*int_0^1 g_ex(s)[sum_low exp(-2pi*i*n*s)+sum_high exp(-2pi*i*n*s)]ds",
            "physical_projection": "Delta Q_fold=4sqrt(2)Re(packet sum)",
        },
        "rows": rows,
        "fine_grid_summary": {
            "physical_min": min(row["physical_fold_increment_beta4"] for row in fine_rows),
            "physical_max": max(row["physical_fold_increment_beta4"] for row in fine_rows),
            "maximum_grid_difference": max(row["difference_from_previous_grid"] for row in fine_rows),
        },
        "decision": {
            "exact_two_packet_frequency_identity_recorded": True,
            "float64_fft_is_rigorous_certificate": False,
            "exact_profile_remainder_applied": False,
            "signed_fold_increment_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_action": "Use the observed scale and grid convergence only to decide whether to build a ball-arithmetic two-packet quadrature with an analytic midpoint/alias remainder and the certified exact-profile sup norm.",
        "proof_boundary": "Float64 FFT scout and exact frequency bookkeeping only. No interval quadrature, exact-profile remainder transfer, signed fold increment, T_upper, Lambda<=0, RH, or prize-level conclusion.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(
        "# Paired-increment two-packet FFT scout\n\n"
        "Date: 2026-08-13\n\n"
        "Status: nonrigorous float64 diagnostic only\n\n"
        "The exact common-profile frequency blocks are `-200..-1` and "
        "`-79789..-79590`. Two midpoint FFT grids were compared.\n\n"
        f"Fine-grid physical range: `{artifact['fine_grid_summary']['physical_min']}` to "
        f"`{artifact['fine_grid_summary']['physical_max']}`.\n\n"
        f"Maximum complex grid difference: `{artifact['fine_grid_summary']['maximum_grid_difference']}`.\n\n"
        "These values are not certified bounds and are not promoted into the formal core.\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(artifact["fine_grid_summary"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
