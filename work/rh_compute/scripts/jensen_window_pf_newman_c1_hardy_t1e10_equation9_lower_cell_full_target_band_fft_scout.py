#!/usr/bin/env python3
"""Float64 scout for the completed one-cell source minus full target band."""

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
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_cell_full_target_band_fft_scout"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
C = 159_577
PROFILE_ORIGIN = 39_895
TARGET_FIRST = 622
ORDINARY_LAST = 39_694
TARGET_LAST = 39_894
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


def beta4_profile(s: np.ndarray | float, lam: float, beta: float, y_max: float) -> np.ndarray:
    y = y_max * s
    x = lam + y
    ai, aip, _, _ = airy(-x)
    a_x = -aip
    epsilon = beta**-2
    u1 = -(13 * lam + 3 * y) / 60
    v1 = (8 * lam**2 - 4 * lam * y + 3 * y**2) / 60
    n_u2 = (
        448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2
        - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2
    )
    n_v2 = 40 * lam**3 - 20 * lam**2 * y + lam * y**2 - 9 * y**3 + 27
    u = 1 + epsilon * u1 - epsilon**2 * n_u2 / 50400
    v = epsilon * v1 + epsilon**2 * n_v2 / 1680
    phase = -1.5 * np.pi * s + y_max**2 * s**2 / (4 * beta)
    return y_max * np.exp(1j * phase) * (u * ai + v * a_x)


def packet_sum(transform: np.ndarray, frequencies: np.ndarray) -> complex:
    panels = transform.size
    indices = np.mod(frequencies, panels)
    midpoint_phase = np.exp(-1j * np.pi * frequencies / panels)
    return complex(np.sum(transform[indices] * midpoint_phase))


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    tstar = np.pi * C**2 / 8
    beta = tstar ** (1 / 3)
    lambda_max = np.pi / (16 * beta)
    y_max = np.pi * C / (2 * beta)
    ordinary_frequencies = np.arange(TARGET_FIRST - PROFILE_ORIGIN, ORDINARY_LAST - PROFILE_ORIGIN + 1)
    fold_frequencies = np.arange(ORDINARY_LAST + 1 - PROFILE_ORIGIN, TARGET_LAST - PROFILE_ORIGIN + 1)
    target_frequencies = np.arange(TARGET_FIRST - PROFILE_ORIGIN, TARGET_LAST - PROFILE_ORIGIN + 1)

    rows = []
    previous: dict[int, complex] = {}
    for panels in GRIDS:
        s = (np.arange(panels, dtype=np.float64) + 0.5) / panels
        for index in range(HEIGHT_NODES):
            lam = lambda_max * index / (HEIGHT_NODES - 1)
            profile = beta4_profile(s, lam, beta, y_max)
            transform = np.fft.fft(profile) / panels
            ordinary = packet_sum(transform, ordinary_frequencies)
            fold = packet_sum(transform, fold_frequencies)
            target = packet_sum(transform, target_frequencies)
            source = complex(beta4_profile(0.0, lam, beta, y_max))
            residual = source - target
            physical_source = 4 * np.sqrt(2) * source.real
            physical_target = 4 * np.sqrt(2) * target.real
            physical_residual = 4 * np.sqrt(2) * residual.real
            prior = previous.get(index)
            rows.append({
                "panels": panels,
                "height_fraction_index_over_4": index,
                "lambda": lam,
                "source_completed_value": {"real": source.real, "imag": source.imag},
                "ordinary_packet": {"real": ordinary.real, "imag": ordinary.imag},
                "fold_packet": {"real": fold.real, "imag": fold.imag},
                "full_target_packet": {"real": target.real, "imag": target.imag},
                "completed_source_minus_target_packet": {
                    "real": residual.real,
                    "imag": residual.imag,
                    "abs": abs(residual),
                },
                "physical_source": physical_source,
                "physical_target": physical_target,
                "physical_residual": physical_residual,
                "residual_difference_from_previous_grid": None if prior is None else abs(residual - prior),
                "ordinary_plus_fold_consistency": abs((ordinary + fold) - target),
            })
            previous[index] = residual

    fine_rows = [row for row in rows if row["panels"] == GRIDS[-1]]
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "nonrigorous_float64_completed_one_cell_full_target_band_scout",
        "passed": True,
        "geometry": {
            "profile_frequency_index": "n=m-39895",
            "ordinary_target_modes": [TARGET_FIRST, ORDINARY_LAST],
            "ordinary_frequencies": [TARGET_FIRST - PROFILE_ORIGIN, ORDINARY_LAST - PROFILE_ORIGIN],
            "fold_target_modes": [ORDINARY_LAST + 1, TARGET_LAST],
            "fold_frequencies": [ORDINARY_LAST + 1 - PROFILE_ORIGIN, TARGET_LAST - PROFILE_ORIGIN],
            "full_target_modes": [TARGET_FIRST, TARGET_LAST],
            "full_target_frequencies": [TARGET_FIRST - PROFILE_ORIGIN, TARGET_LAST - PROFILE_ORIGIN],
            "target_mode_count": int(target_frequencies.size),
            "canonical_identity": "R_cell=g_ex(0)-sum_(m=622)^39894 G_ex,m",
            "physical_projection": "R_cell_phys=4sqrt(2)Re(R_cell)",
        },
        "rows": rows,
        "fine_grid_summary": {
            "physical_source_min": min(row["physical_source"] for row in fine_rows),
            "physical_source_max": max(row["physical_source"] for row in fine_rows),
            "physical_target_min": min(row["physical_target"] for row in fine_rows),
            "physical_target_max": max(row["physical_target"] for row in fine_rows),
            "physical_residual_min": min(row["physical_residual"] for row in fine_rows),
            "physical_residual_max": max(row["physical_residual"] for row in fine_rows),
            "maximum_residual_grid_difference": max(row["residual_difference_from_previous_grid"] for row in fine_rows),
            "maximum_packet_split_error": max(row["ordinary_plus_fold_consistency"] for row in fine_rows),
        },
        "decision": {
            "full_target_band_frequency_identity_recorded": True,
            "ordinary_and_fold_target_packets_kept_in_one_completed_cell": True,
            "float64_fft_is_rigorous_certificate": False,
            "exact_profile_remainder_applied": False,
            "completed_one_cell_residual_bounded": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_action": (
            "Use the observed full-band residual scale to choose between a rigorous endpoint-corrected Fourier "
            "packet enclosure and a proof that the one-cell localization still omits a mandatory global current."
        ),
        "proof_boundary": (
            "Float64 midpoint-FFT scout and exact frequency bookkeeping only. No interval quadrature, exact-profile "
            "transfer, fixed-state residual, T_upper, Lambda<=0, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(
        "# Completed one-cell full-target-band FFT scout\n\n"
        "Date: 2026-08-13\n\n"
        "Status: nonrigorous float64 diagnostic only\n\n"
        "The exact target band `622..39894` is kept intact and split only for reporting. "
        "The local canonical residual is `g_ex(0)-sum_target G_ex,m`.\n\n"
        f"Fine-grid physical source range: `{artifact['fine_grid_summary']['physical_source_min']}` to "
        f"`{artifact['fine_grid_summary']['physical_source_max']}`.\n\n"
        f"Fine-grid physical target-packet range: `{artifact['fine_grid_summary']['physical_target_min']}` to "
        f"`{artifact['fine_grid_summary']['physical_target_max']}`.\n\n"
        f"Fine-grid physical residual range: `{artifact['fine_grid_summary']['physical_residual_min']}` to "
        f"`{artifact['fine_grid_summary']['physical_residual_max']}`.\n\n"
        f"Maximum complex residual grid difference: `{artifact['fine_grid_summary']['maximum_residual_grid_difference']}`.\n\n"
        "These values are not certified bounds and are not promoted into the formal core.\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(artifact["fine_grid_summary"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
