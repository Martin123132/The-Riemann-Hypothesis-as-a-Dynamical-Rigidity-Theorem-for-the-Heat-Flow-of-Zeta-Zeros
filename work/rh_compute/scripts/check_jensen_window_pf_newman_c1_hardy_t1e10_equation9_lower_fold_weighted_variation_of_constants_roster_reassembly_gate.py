#!/usr/bin/env python3
"""Independently check weighted variation-of-constants roster reassembly."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_variation_of_constants_roster_reassembly_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ball_upper(text: str) -> mp.mpf:
    if "+/-" not in text:
        return mp.mpf(text.strip("[] "))
    midpoint, radius = text.strip("[] ").split("+/-")
    return mp.mpf(midpoint.strip()) + mp.mpf(radius.strip())


def finite_fubini_check(count: int, negative: bool) -> None:
    weights = {j: Fraction(3 * j + 2, 5 * j + 7) for j in range(1, count + 1)}
    kernels = {(j, k): Fraction(7 * j - 2 * k + 1, 11 * j + 13 * k + 17) for j in range(1, count + 1) for k in range(1, j + 1)}
    direct = sum(weights[j] * sum((kernels[(j, k)] for k in range(1, j + 1)), Fraction()) for j in range(1, count + 1))
    swapped = sum(sum((weights[j] * kernels[(j, k)] for j in range(k, count + 1)), Fraction()) for k in range(1, count + 1))
    if negative:
        direct = -direct
        swapped = -swapped
    require(direct == swapped, "finite Fubini orientation failed")


def main() -> None:
    require(RESULT.is_file(), "missing result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "stored gate is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")

    finite_fubini_check(200, False)
    finite_fubini_check(198, True)
    algebra = artifact["exact_algebra"]
    require(algebra["positive_source_cell_incidences"] == 20100, "positive incidence count failed")
    require(algebra["negative_source_cell_incidences"] == 19701, "negative incidence count failed")

    beta, c, m = sp.symbols("beta C m", nonzero=True)
    d = beta * (4 * m / c - 1)
    chi = sp.expand(2 * beta**2 * d + beta * d**2)
    beta3 = sp.pi * c**2 / 8
    require(sp.simplify(chi.subs(beta**3, beta3) - (2 * sp.pi * m**2 - beta3)) == 0, "gauge check failed")

    mp.mp.dps = 60
    operator = artifact["operator_certificate"]
    minus_k = ball_upper(operator["minus_branch_weighted_K_coefficient_ball"])
    plus_k = ball_upper(operator["plus_branch_weighted_K_coefficient_ball"])
    minus_d = ball_upper(operator["minus_branch_weighted_partial_s_K_coefficient_ball"])
    plus_d = ball_upper(operator["plus_branch_weighted_partial_s_K_coefficient_ball"])
    require(minus_k + plus_k < mp.mpf("7.58e-7"), "combined K coefficient failed")
    require(minus_d + plus_d < mp.mpf("0.001634"), "combined derivative coefficient failed")
    require(artifact["decision"]["completed_initial_data_cancellation_proved"] is False, "proof boundary drift")
    print("independently checked weighted variation-of-constants roster reassembly", flush=True)


if __name__ == "__main__":
    main()
