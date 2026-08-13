#!/usr/bin/env python3
"""Independent replay of the B face-window Fresnel-remainder certificate."""

from __future__ import annotations

import hashlib
import json
import re
from decimal import Decimal
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

T = 10_000_000_000
B = 5_122_421


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ball_upper(text: str) -> Decimal:
    match = re.fullmatch(r"\[([^\]]+) \+/- ([^\]]+)\]", text)
    require(match is not None, f"unparsed Arb ball: {text}")
    midpoint, radius = (Decimal(value) for value in match.groups())
    return midpoint + radius


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")

    mp.mp.dps = 100
    t = mp.mpf(T)
    endpoint = mp.mpf(B)
    x0 = (1 - mp.sqrt(1 - 8 * t / (mp.pi * endpoint**2))) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_h = mp.sqrt(hessian)
    m0 = endpoint * x0 / 2
    spacing = 2 * root_h / endpoint
    xis = {mode: spacing * (mode - m0) for mode in (620, 621, 622, 623)}
    require(xis[620] < -70 < xis[621], "left crossing roster replay failed")
    require(xis[622] < 70 < xis[623], "right crossing roster replay failed")

    x_low = x0 - 70 / root_h
    x_high = x0 + 70 / root_h
    c_low = endpoint * x_low / 2
    c_high = endpoint * x_high / 2
    q620 = (endpoint * x_low - 2 * 620) / mp.sqrt(2 * x_low)
    q623 = abs((endpoint * x_high - 2 * 623) / mp.sqrt(2 * x_high))
    require(q620 > 85 and q623 > 75, "independent q-margin replay failed")

    weighted = mp.mpf(0)
    for mode in range(1, 39_895):
        if mode in (621, 622):
            continue
        distance = c_low - mode if mode <= 620 else mode - c_high
        weighted += mode / distance**7
    first_outer = mp.mpf(39_895)
    y0 = first_outer - c_high
    outer = first_outer / y0**7 + 1 / (5 * y0**5) + c_high / (6 * y0**6)
    bound = 15 * x_high**2 * (weighted + outer) / (4 * mp.pi**4)
    require(bound < mp.mpf("7.48e-6"), "independent infinite remainder bound failed")

    # Directly test the three-term expansion at the certified worst q scale.
    r = mp.mpf(75)
    e = mp.exp(mp.j * mp.pi * r**2 / 2)
    exact_tail = (1 + mp.j) * mp.erfc(mp.e ** (-mp.j * mp.pi / 4) * mp.sqrt(mp.pi / 2) * r) / 2
    approximation = mp.j * e / (mp.pi * r) + e / (mp.pi**2 * r**3) - 3 * mp.j * e / (mp.pi**3 * r**5)
    require(abs(exact_tail - approximation) <= 30 / (mp.pi**4 * r**7), "direct Fresnel remainder witness failed")

    certificate = artifact["interval_certificate"]
    require(ball_upper(certificate["infinite_nonlocal_completed_remainder_ball"]) < Decimal("7.48e-6"), "stored remainder ball too large")
    decision = artifact["decision"]
    require(decision["xi_70_window_contains_exactly_B_crossings_621_622"] is True, "window decision drift")
    require(decision["nonlocal_completed_remainder_absolutely_summable"] is True, "summability decision drift")
    require(decision["leading_rational_terms_summed_without_pole_subtraction"] is False, "pole guard lost")
    require(decision["complete_B_face_estimate_proved"] is False, "B-face overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked B face window and Fresnel remainder independently", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
