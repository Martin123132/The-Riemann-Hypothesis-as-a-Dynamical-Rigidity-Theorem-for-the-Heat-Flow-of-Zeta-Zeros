#!/usr/bin/env python3
"""Independent replay of the local B normal-safe Fresnel remainder gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_local_normal_safe_fresnel_remainder_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

T = 10_000_000_000
B = 5_122_421
SPLITS = {621: Fraction("0.00024238523659"), 622: Fraction("0.000242916437095")}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ball(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def three_halves_bound(mode: int) -> arb:
    endpoint = Fraction(B)
    split = SPLITS[mode]
    if mode == 621:
        current = Fraction(mode) - endpoint * split / 2
        last = Fraction(mode)
        sign = -1
    else:
        current = endpoint * split / 2 - Fraction(mode)
        last = endpoint / 4 - Fraction(mode)
        sign = 1

    total = arb(0)
    while current < last:
        following = min(Fraction(3, 2) * current, last)
        endpoint_distance = current if sign < 0 else following
        x_endpoint = 2 * (Fraction(mode) + sign * endpoint_distance) / endpoint
        x = ball(x_endpoint)
        h = x ** (arb(7) / 4) * (1 - x) ** (-arb(1) / 4)
        total += h * (ball(current) ** -6 - ball(following) ** -6) / 6
        current = following

    pi = arb.pi()
    raw = arb(15 * mode) * total / (2 * pi**4 * B)
    return 2 * (pi / (32 * arb(T))) ** (arb(1) / 4) * raw


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    # Replay the Jacobian and weight algebra from d=|Bx/2-m| rather than
    # importing production helpers.
    x, d, mode, endpoint = sp.symbols("x d m B", positive=True, real=True)
    density = (x * (1 - x)) ** (-sp.Rational(1, 4)) * 15 * mode * x**2 / (4 * sp.pi**4 * d**7)
    changed = sp.simplify(density * 2 / endpoint)
    expected = 15 * mode * x ** sp.Rational(7, 4) * (1 - x) ** (-sp.Rational(1, 4)) / (2 * sp.pi**4 * endpoint * d**7)
    require(sp.simplify(changed - expected) == 0, "independent distance-Jacobian algebra failed")

    ctx.dps = 120
    ctx.threads = 1
    b621 = three_halves_bound(621)
    b622 = three_halves_bound(622)
    combined = b621 + b622
    require(b621 < arb("4.02e-11"), "independent mode-621 bound failed")
    require(b622 < arb("1.56e-10"), "independent mode-622 bound failed")
    require(combined < arb("1.97e-10"), "independent combined bound failed")

    saved = arb(artifact["interval_certificate"]["combined_physically_normalized_remainder_ball"])
    # The meshes are intentionally different upper sums. They need not overlap,
    # but each must prove the same advertised strict threshold.
    require(saved < arb("1.97e-10"), "saved combined bound drift")
    decision = artifact["decision"]
    require(decision["exact_equation9_normalization_restored"] is True, "normalization decision drift")
    require(decision["complementary_outer_safe_local_current_bounded"] is False, "outer-safe overclaim")
    require(decision["complete_B_face_estimate_proved"] is False, "complete-B overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print(
        "independently checked local B normal-safe remainders with 3/2 distance panels: "
        f"621={b621}, 622={b622}, combined={combined}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
