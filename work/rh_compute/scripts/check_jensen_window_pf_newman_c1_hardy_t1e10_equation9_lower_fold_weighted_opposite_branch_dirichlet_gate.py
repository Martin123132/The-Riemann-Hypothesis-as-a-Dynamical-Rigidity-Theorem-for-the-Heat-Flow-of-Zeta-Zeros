#!/usr/bin/env python3
"""Independently check the weighted opposite-branch Dirichlet bound."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_opposite_branch_dirichlet_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
L = 39_696
Q = 39_894
U = 40_094
PANELS = 32_768


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_ball(text: str) -> tuple[float, float]:
    body = text.strip("[] ")
    if "+/-" not in body:
        return float(body), 0.0
    midpoint, radius = body.split("+/-", 1)
    return float(midpoint), float(radius)


def parse_complex(record: dict[str, str]) -> tuple[complex, float]:
    real, real_radius = parse_ball(record["real_ball"])
    imag, imag_radius = parse_ball(record["imag_ball"])
    return complex(real, imag), math.hypot(real_radius, imag_radius)


def add(left: tuple[complex, float], right: tuple[complex, float]) -> tuple[complex, float]:
    return left[0] + right[0], left[1] + right[1]


def subtract(left: tuple[complex, float], right: tuple[complex, float]) -> tuple[complex, float]:
    return left[0] - right[0], left[1] + right[1]


def conjugate(value: tuple[complex, float]) -> tuple[complex, float]:
    return value[0].conjugate(), value[1]


def reconstruct(pairing: dict) -> dict[int, tuple[complex, float]]:
    weights = {}
    c = pairing["leading_multiplier_certificate"]
    for row in c["pair_rows"]:
        coherent = parse_complex(row["coherent_coefficient_ball"])
        leakage = parse_complex(row["conjugacy_leakage_coefficient_ball"])
        weights[int(row["minus_mode"])] = add(coherent, leakage)
        weights[int(row["plus_mode"])] = conjugate(subtract(coherent, leakage))
    for row in c["edge_rows"]:
        weights[int(row["mode"])] = parse_complex(row["leading_defect_ball"])
    require(set(weights) == set(range(L, Q)) | set(range(Q + 1, U + 1)), "independent roster drift")
    return weights


def l1_bound(coefficients: list[tuple[complex, float]]) -> float:
    midpoints = np.array([value[0] for value in coefficients], dtype=np.complex128)
    radii = np.array([value[1] for value in coefficients], dtype=np.float64)
    exponents = np.arange(len(coefficients), dtype=np.float64)
    derivative = 2 * np.pi * np.sum(exponents * (np.abs(midpoints) + radii))
    uncertainty = float(np.sum(radii))
    total = 0.0
    chunk = 1024
    for start in range(0, PANELS, chunk):
        stop = min(start + chunk, PANELS)
        s = (2 * np.arange(start, stop, dtype=np.float64) + 1) / (2 * PANELS)
        values = np.exp(-2j * np.pi * np.outer(s, exponents)) @ midpoints
        total += float(np.sum(np.abs(values) + uncertainty + 1e-15))
    return total / PANELS + float(derivative) / (2 * PANELS)


def ball_upper(text: str) -> float:
    midpoint, radius = parse_ball(text)
    return midpoint + radius


def main() -> None:
    require(RESULT.is_file(), "missing result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "stored gate is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")
    for source in artifact["sources"].values():
        path = REPO_ROOT / source["path"]
        require(path.is_file(), f"missing source: {path}")
        require(file_hash(path) == source["sha256"], f"source hash drift: {path}")

    pairing_path = REPO_ROOT / artifact["dependencies"]["event_pairing"]["path"]
    weights = reconstruct(json.loads(pairing_path.read_text(encoding="utf-8")))
    minus = l1_bound([weights[mode] for mode in range(L, Q)])
    plus = l1_bound([weights[mode] for mode in range(Q + 1, U + 1)])
    require(minus < 1.6e-5, "independent negative kernel L1 failed")
    require(plus < 2.4e-5, "independent positive kernel L1 failed")
    require(minus + plus < 4.0e-5, "independent two-kernel L1 failed")

    c = artifact["certificate"]
    require(ball_upper(c["two_branch_kernel_L1_bound_ball"]) < 4.0e-5, "stored L1 bound failed")
    require(ball_upper(c["weighted_opposite_branch_bound_ball"]) < 0.00166, "stored opposite bound failed")
    require(ball_upper(c["selected_weighted_branch_bound_ball"]) < 0.02324, "stored selected bound failed")
    decision = artifact["decision"]
    require(decision["weighted_opposite_branch_bound_proved"] is True, "opposite decision drift")
    require(decision["exact_finite_integral_amplitude_remainder_proved"] is False, "proof boundary drift")
    print("independently checked weighted opposite-branch Dirichlet bound", flush=True)


if __name__ == "__main__":
    main()
