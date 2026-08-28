#!/usr/bin/env python3
"""Independently check the finite Q_K geometric compression and saddle route guard."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_QK_geometric_half_line_and_saddle_route_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def direct_half_line(t: mp.mpf, alpha: int) -> mp.mpc:
    s = mp.mpf("0.5") + 1j * t
    lam = mp.pi * (1 + 1j) / mp.sqrt(2)
    f = lambda x: x ** (-s) * mp.e ** (-mp.pi * x * x - lam * alpha * x)
    return mp.e ** (3 * mp.pi * t / 4 + 3j * mp.pi / 8) * mp.quad(
        f, [0, mp.mpf("0.01"), mp.mpf("0.1"), 1, 4, mp.inf]
    )


def compressed_half_line(t: mp.mpf, alpha_min: int, alpha_max: int) -> mp.mpc:
    s = mp.mpf("0.5") + 1j * t
    lam = mp.pi * (1 + 1j) / mp.sqrt(2)
    count = (alpha_max - alpha_min) // 2 + 1

    def source(x: mp.mpf) -> mp.mpc:
        if x == 0:
            return mp.mpf(count)
        return (
            mp.e ** (-lam * alpha_min * x)
            * mp.expm1(-2 * lam * count * x)
            / mp.expm1(-2 * lam * x)
        )

    f = lambda x: x ** (-s) * mp.e ** (-mp.pi * x * x) * source(x)
    return mp.e ** (3 * mp.pi * t / 4 + 3j * mp.pi / 8) * mp.quad(
        f, [0, mp.mpf("0.01"), mp.mpf("0.1"), 1, 4, mp.inf]
    )


def main() -> int:
    require(gate.RESULT.is_file(), "missing geometric-compression result")
    require(gate.NOTE.is_file(), "missing geometric-compression note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == gate.STEM, "kind drift")
    require(
        artifact["status"]
        == "exact_finite_QK_geometric_half_line_compression_and_guarded_saddle_route_certified",
        "status drift",
    )
    require(artifact.get("passed") is True, "geometric-compression gate failed")

    mp.mp.dps = 70
    t = mp.mpf("3.75")
    alpha_min = 3
    alpha_max = 9
    direct = mp.fsum(direct_half_line(t, alpha) for alpha in range(alpha_min, alpha_max + 1, 2))
    compressed = compressed_half_line(t, alpha_min, alpha_max)
    require(abs(direct - compressed) < mp.mpf("1e-35"), "fresh geometric compression failed")

    threshold = mp.sqrt(8 * mp.mpf(gate.HEIGHT) / mp.pi)
    require(gate.SOURCE_ALPHA_MIN - 2 < threshold < gate.SOURCE_ALPHA_MIN, "fresh transition bracket failed")
    discriminant_before = -4 * mp.pi + 1j * mp.pi * (mp.pi * (gate.SOURCE_ALPHA_MIN - 2) ** 2 - 8 * gate.HEIGHT)
    discriminant_after = -4 * mp.pi + 1j * mp.pi * (mp.pi * gate.SOURCE_ALPHA_MIN**2 - 8 * gate.HEIGHT)
    require(mp.im(discriminant_before) < 0 < mp.im(discriminant_after), "fresh discriminant crossing failed")

    test_t = mp.mpf("2.75")
    alpha = 5
    y = mp.mpf("1.3")
    s = mp.mpf("0.5") + 1j * test_t
    lam = mp.pi * (1 + 1j) / mp.sqrt(2)
    K0 = mp.e ** (3 * mp.pi * test_t / 4 + 3j * mp.pi / 8)
    rotation = mp.e ** (-3j * mp.pi / 4)
    x = rotation * y
    left = K0 * mp.e ** (-s * (mp.log(y) - 3j * mp.pi / 4)) * mp.e ** (-mp.pi * x * x - lam * alpha * x) * rotation
    right = y ** (-s) * mp.e ** (1j * mp.pi * (alpha * y - y * y))
    require(abs(left - right) < mp.mpf("1e-65"), "fresh candidate-ray integrand identity failed")

    exact = artifact["exact_geometric_half_line"]
    require(exact["denominator_nonzero_for_real_x_positive"] is True, "denominator guard lost")
    require(exact["finite_sum_integral_interchange_requires_limit_theorem"] is False, "finite sum treated as infinite")
    candidate = artifact["candidate_saddle_ray"]
    require(candidate["integrand_normalization_cancellation_exact"] is True, "candidate normalization identity lost")
    require(candidate["positive_half_line_to_saddle_ray_contour_deformation_certified"] is False, "contour deformation overclaim")
    require(candidate["Stokes_multiplier_and_large_arc_contributions_owned"] is False, "Stokes ownership overclaim")
    decision = artifact["decision"]
    require(decision["finite_QK_roster_compressed_to_one_exact_integral"] is True, "finite compression lost")
    require(decision["direct_Kummer_roster_evaluation_required"] is False, "direct Kummer route reintroduced")
    require(decision["candidate_saddle_ray_permitted_as_integral_identity"] is False, "candidate ray promoted")
    require(decision["actual_height_QK_enclosed"] is False, "Q_K overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in (
        "exact scalar integral",
        "first odd label above",
        "not yet permission to rotate",
        "Stokes multiplier",
        "not inserted or fitted",
        "prize-level enclosure is proved",
    ):
        require(token in note, f"note boundary token missing: {token}")

    print("independently checked finite Q_K geometric compression and saddle-route guard", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
