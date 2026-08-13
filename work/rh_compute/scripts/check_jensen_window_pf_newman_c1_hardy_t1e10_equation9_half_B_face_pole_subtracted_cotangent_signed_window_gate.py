#!/usr/bin/env python3
"""Independent coarser-panel replay of the signed B cotangent window."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_pole_subtracted_cotangent_signed_window_gate"
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


def main() -> int:
    import psutil

    psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS if os.name == "nt" else 10)
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    # Independent algebraic derivation of the regular cotangent coefficients.
    delta = sp.symbols("delta")
    expansion = sp.series(sp.pi * sp.cot(sp.pi * delta) - 1 / delta, delta, 0, 12)
    require(expansion.removeO().coeff(delta, 1) == -sp.pi**2 / 3, "cotangent linear coefficient failed")
    require(expansion.removeO().coeff(delta, 3) == -sp.pi**4 / 45, "cotangent cubic coefficient failed")

    ctx.dps = 82
    ctx.threads = 1
    t, endpoint, pi, imaginary = arb(T), arb(B), arb.pi(), acb(0, 1)
    x0 = (1 - (1 - 8 * t / (pi * endpoint**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_hessian = hessian.sqrt()
    log0 = ((1 - x0) / x0).log()
    phase0 = pi * endpoint**2 * x0 / 4 + t * log0 / 2
    normalizer = (pi / (32 * t)) ** (arb(1) / 4)

    def regular(delta_ball: acb) -> acb:
        radius = abs(delta_ball).upper()
        if radius < arb("0.3"):
            value = acb(0)
            for index in range(1, 13):
                value -= 2 * arb(2 * index).zeta() * delta_ball ** (2 * index - 1)
            remainder = 4 * radius**25 / (1 - radius**2)
            return value + acb(arb(0, remainder), arb(0, remainder))
        return pi * (pi * delta_ball).cot() - 1 / delta_ball

    def function(local_mode: int):
        other = 1243 - local_mode

        def integrand(xi: acb, _: bool) -> acb:
            x = x0 + xi / root_hessian
            c = endpoint * x / 2
            reg = regular(c - local_mode)
            s_hat = reg / (8 * c) - 1 / (8 * c * (c + local_mode)) - 1 / (8 * c**2) - 1 / (4 * (c**2 - other**2))
            current = -2 * (2 + imaginary * pi * endpoint**2 * x) * s_hat / pi**2
            phase = pi * endpoint**2 * (x - x0) / 4 + t * (((1 - x) / x).log() - log0) / 2
            return (x * (1 - x)) ** (-arb(1) / 4) * current * (imaginary * phase).exp() / root_hessian

        return integrand

    cuts = [arb(-70 + 2 * index) for index in range(71)]
    cuts.extend(
        root_hessian * (2 * arb(value) / endpoint - x0)
        for value in ("621", "621.5", "622")
    )
    cuts = sorted(cuts, key=float)
    total = acb(0)
    for left, right in zip(cuts, cuts[1:]):
        if right <= left:
            continue
        center = (left + right) / 2
        c_center = endpoint * (x0 + center / root_hessian) / 2
        local_mode = 621 if c_center < arb("621.5") else 622
        total += acb.integral(
            function(local_mode),
            left,
            right,
            abs_tol=arb("2e-14"),
            rel_tol=arb("2e-14"),
            eval_limit=300_000,
            depth_limit=42,
        )

    rotation = (imaginary * phase0 - imaginary * pi / 8).exp()
    physical = 2 * normalizer * (rotation * total).real
    saved = arb(artifact["numerical_certificate"]["physically_normalized_signed_projection_ball"])
    require(physical.overlaps(saved), "independent signed projection misses production enclosure")
    require(arb("4.791e-6") < physical < arb("4.792e-6"), "independent physical bound failed")
    require(2 * normalizer * abs(total) > arb("2.07e-5"), "modulus guard disappeared")

    decision = artifact["decision"]
    require(decision["first_pole_subtracted_current_integrated_before_absolute_value"] is True, "signed-integration decision drift")
    require(decision["exact_local_621_622_replacements_added"] is False, "local-replacement overclaim")
    require(decision["complete_B_face_estimate_proved"] is False, "complete-B overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print(f"independently validated signed B cotangent window projection: {physical}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
