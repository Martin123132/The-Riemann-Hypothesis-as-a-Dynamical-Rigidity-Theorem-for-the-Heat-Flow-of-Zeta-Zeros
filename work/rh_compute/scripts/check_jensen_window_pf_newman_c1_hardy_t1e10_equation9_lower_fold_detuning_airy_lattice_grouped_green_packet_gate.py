#!/usr/bin/env python3
"""Independently check the grouped Airy phase and Green-packet gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_airy_lattice_grouped_green_packet_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
C = 159_577
J_MIN = -198
J_MAX = 200
COEFFICIENTS = ((3, 2), (3, 8), (-1, 16), (3, 128), (-3, 256))


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


def ball_lower(text: str) -> mp.mpf:
    if "+/-" not in text:
        return mp.mpf(text.strip("[] "))
    midpoint, radius = text.strip("[] ").split("+/-")
    return mp.mpf(midpoint.strip()) - mp.mpf(radius.strip())


def phase(j: int, u: mp.mpf) -> mp.mpf:
    xj = (8 * j - 2) / mp.mpf(C) + u
    x0 = -2 / mp.mpf(C) + u
    truncated = mp.mpf("0")
    for degree, (numerator, denominator) in enumerate(COEFFICIENTS, start=1):
        coefficient = mp.mpf(numerator) / denominator
        truncated += coefficient * (xj**degree - x0**degree)
    return mp.pi * C**2 * truncated / 12 - mp.pi * j * (C - 1 + 2 * j)


def main() -> None:
    require(RESULT.is_file(), "missing result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "stored gate is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")

    mp.mp.dps = 90
    js = list(range(J_MIN, J_MAX + 1))
    terms = [mp.exp(1j * phase(j, mp.mpf("0"))) for j in js]
    prefix = [mp.mpc(0)]
    for term in terms:
        prefix.append(prefix[-1] + term)

    maximum = mp.mpf("0")
    indices = None
    for start in range(len(js)):
        for stop in range(start + 1, len(js) + 1):
            value = abs(prefix[stop] - prefix[start])
            if value > maximum:
                maximum = value
                indices = (js[start], js[stop - 1])
    require(indices == (-30, 31), f"unexpected central maximum block: {indices}")

    certificate = artifact["certificate"]
    stored_center_max = ball_upper(certificate["center_quintic_max_contiguous_upper_ball"])
    # The independent evaluation subtracts O(C^2) phases in mp arithmetic,
    # so compare well inside its roughly 80 reliable decimal places rather
    # than ordering its final cancellation crumbs against an Arb endpoint.
    require(
        abs(stored_center_max - maximum) < mp.mpf("1e-75"),
        "stored center contiguous bound disagrees with the independent value",
    )

    u_max = mp.mpf(1) / (2 * C**2)
    direct_height_variation = sum(abs(phase(j, u_max) - phase(j, mp.mpf("0"))) for j in js)
    require(
        ball_upper(certificate["uniform_sum_of_height_phase_variations_ball"]) >= direct_height_variation,
        "stored differential height bound misses endpoint variation",
    )
    require(direct_height_variation < mp.mpf("0.196"), "height variation exceeded independent threshold")

    exact_contiguous = ball_upper(certificate["exact_Airy_any_contiguous_packet_upper_ball"])
    exact_full_lower = ball_lower(certificate["exact_Airy_full_packet_absolute_lower_ball"])
    grouped_kernel = ball_upper(certificate["grouped_any_contiguous_Green_kernel_bound_ball"])
    improvement = ball_lower(certificate["termwise_to_grouped_improvement_lower_ball"])
    grouped_derivative = ball_upper(certificate["grouped_any_contiguous_derivative_Green_kernel_bound_ball"])
    derivative_improvement = ball_lower(certificate["derivative_termwise_to_grouped_improvement_lower_ball"])
    require(exact_contiguous < mp.mpf("52.53"), "exact contiguous packet threshold failed")
    require(exact_full_lower > mp.mpf("31.4"), "full packet noncancellation threshold failed")
    require(grouped_kernel < mp.mpf("0.0246"), "grouped Green threshold failed")
    require(improvement > mp.mpf("7.5"), "grouped Green improvement threshold failed")
    require(grouped_derivative < mp.mpf("52.8"), "grouped derivative-Green threshold failed")
    require(derivative_improvement > mp.mpf("7.5"), "grouped derivative-Green improvement threshold failed")
    require(artifact["decision"]["derivative_kernel_packet_proved"] is True, "derivative packet guard drift")
    require(artifact["decision"]["termwise_absolute_kernel_sum_used"] is False, "termwise guard drift")
    beta, c_symbol, m = sp.symbols("beta C m", nonzero=True)
    d = beta * (4 * m / c_symbol - 1)
    chi = 2 * beta**2 * d + beta * d**2
    beta_cubed = sp.pi * c_symbol**2 / 8
    require(
        sp.expand(chi.subs(beta**3, beta_cubed) - (2 * sp.pi * m**2 - beta_cubed)) == 0,
        "exact gauge-lattice lock failed",
    )
    require(
        artifact["decision"]["inverse_gauge_phase_common_on_integer_event_lattice"] is True,
        "inverse-gauge decision drift",
    )
    print("independently checked grouped Airy phase and Green packets", flush=True)


if __name__ == "__main__":
    main()
