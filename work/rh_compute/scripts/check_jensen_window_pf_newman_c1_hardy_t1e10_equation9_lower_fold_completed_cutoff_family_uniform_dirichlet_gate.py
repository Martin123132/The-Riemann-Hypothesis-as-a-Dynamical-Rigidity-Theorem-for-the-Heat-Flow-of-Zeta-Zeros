#!/usr/bin/env python3
"""Independently check the uniform completed cutoff-family Dirichlet bound."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_completed_cutoff_family_uniform_dirichlet_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
A = 159_577
N = 399


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ball_upper(text: str) -> mp.mpf:
    body = text.strip("[] ")
    if "+/-" not in body:
        return mp.mpf(body)
    midpoint, radius = body.split("+/-", 1)
    return mp.mpf(midpoint.strip()) + mp.mpf(radius.strip())


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

    mp.mp.dps = 80
    pi = mp.pi
    beta = (pi * A * A / 8) ** (mp.mpf(1) / 3)
    epsilon = beta**-2
    y = pi * A / (2 * beta)
    lam = pi / (16 * beta)
    u_bound = 1 + epsilon * (13 * lam + 3 * y) / 60
    u_bound += epsilon**2 * (
        448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2 + 105 * lam * y**4
        + 30 * lam * y + 63 * y**5 + 405 * y**2
    ) / 50400
    v_bound = epsilon * (8 * lam**2 + 4 * lam * y + 3 * y**2) / 60
    v_bound += epsilon**2 * (40 * lam**3 + 20 * lam**2 * y + lam * y**2 + 9 * y**3 + 27) / 1680

    c = artifact["certificate"]
    require(abs(ball_upper(c["U_absolute_bound_ball"]) - u_bound) < mp.mpf("1e-70"), "U bound mismatch")
    require(abs(ball_upper(c["V_absolute_bound_ball"]) - v_bound) < mp.mpf("1e-70"), "V bound mismatch")
    dirichlet = 1 + mp.log(pi * N / 2)
    require(abs(ball_upper(c["maximum_one_sided_Dirichlet_L1_bound_ball"]) - dirichlet) < mp.mpf("1e-70"), "Dirichlet bound mismatch")

    airy = ball_upper(c["Ai_absolute_bound_ball"])
    airy_prime = ball_upper(c["Ai_prime_absolute_bound_ball"])
    require(airy < mp.mpf("0.711"), "Airy value envelope failed")
    require(airy_prime < mp.mpf("2.066"), "Airy derivative envelope failed")
    independent_profile = y * (u_bound * airy + v_bound * airy_prime)
    independent_completed = independent_profile * (1 + dirichlet)
    variation = ball_upper(
        json.loads((REPO_ROOT / artifact["dependencies"]["weighted_completion"]["path"]).read_text(encoding="utf-8"))["coefficient_certificate"]["endpoint_augmented_total_variation_ball"]
    )
    require(independent_profile < mp.mpf("82.66"), "independent profile bound failed")
    require(independent_completed < mp.mpf("698"), "independent completed-family bound failed")
    require(independent_completed * variation < mp.mpf("0.02159"), "independent weighted roster bound failed")
    require(ball_upper(c["weighted_completed_roster_bound_ball"]) < mp.mpf("0.02159"), "stored weighted bound failed")

    decision = artifact["decision"]
    require(decision["uniform_completed_partial_remainder_bound_proved"] is True, "completed-family decision drift")
    require(decision["exact_finite_integral_amplitude_remainder_proved"] is False, "proof boundary drift")
    print("independently checked uniform completed cutoff-family Dirichlet bound", flush=True)


if __name__ == "__main__":
    main()
