#!/usr/bin/env python3
"""Independently check the beta^-4 completed selector branch projection."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate as point_gate


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
PRECISION = 110


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def correction_second_derivative(beta: arb, left: arb, right: arb, pi: arb) -> arb:
    lam = arb((left + right) / 2, (right - left) / 2)
    lam2 = lam * lam
    lam3 = lam2 * lam
    lam4 = lam3 * lam
    lam5 = lam4 * lam
    epsilon = 1 / (beta * beta)
    du = -arb(13) * epsilon * lam / 60 - epsilon * epsilon * (arb(448) * lam5 + arb(4565) * lam2) / 50400
    v = arb(2) * epsilon * lam2 / 15 + epsilon * epsilon * (arb(40) * lam3 + 27) / 1680
    dup = -arb(13) * epsilon / 60 - epsilon * epsilon * (arb(2240) * lam4 + arb(9130) * lam) / 50400
    dupp = -epsilon * epsilon * (arb(8960) * lam3 + arb(9130)) / 50400
    vp = arb(4) * epsilon * lam / 15 + epsilon * epsilon * arb(120) * lam2 / 1680
    vpp = arb(4) * epsilon / 15 + epsilon * epsilon * arb(240) * lam / 1680
    p = dup - lam * v
    q = du + vp
    pa = dupp - v - lam * vp - lam * q
    qa = p + dup + vpp
    airy = (-lam).airy_ai()
    airy_x = -(-lam).airy_ai(derivative=1)
    return abs(2 * pi * (pa * airy + qa * airy_x))


def main() -> None:
    priority = set_low_priority()
    require(RESULT.is_file(), "missing beta4 completed selector artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "beta4 completed selector artifact is not passed")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file() and file_hash(path) == record["sha256"], f"hash drift: {path}")

    ctx.dps = PRECISION
    p = point_gate.parameters()
    beta = p["beta"]
    pi = p["pi"]
    epsilon = 1 / (beta * beta)
    radius = pi / (16 * beta)
    point = json.loads((REPO_ROOT / artifact["dependencies"]["point_completed_strip"]["path"]).read_text(encoding="utf-8"))
    cell = json.loads((REPO_ROOT / artifact["dependencies"]["leading_completed_strip"]["path"]).read_text(encoding="utf-8"))

    exact_minus_leading = arb(point["numerical_certificate"]["reduced_exact_minus_canonical_fold_ball"])
    beta4 = -arb(9) / 560 * 2 * pi * arb(0).airy_ai(derivative=1) / (beta * beta * beta * beta)
    residual0 = exact_minus_leading - beta4
    require(residual0.overlaps(arb(artifact["center_certificate"]["exact_minus_beta4_reduced_ball"])), "center residual does not overlap")
    require(abs(residual0) < arb("6.6e-22"), "independent center residual ceiling failed")

    derivative = arb(cell["center_derivative"]["rigorous_derivative_ball"])
    beta4_derivative = -arb(13) / 60 * 2 * pi * arb(0).airy_ai() / (beta * beta)
    residual1 = derivative - beta4_derivative
    require(residual1.overlaps(arb(artifact["height_cell_certificate"]["center_first_derivative_residual_ball"])), "derivative residual does not overlap")

    panel_bound = arb(0)
    panels = 16
    for index in range(panels):
        left = -radius + 2 * radius * index / panels
        right = -radius + 2 * radius * (index + 1) / panels
        candidate = correction_second_derivative(beta, left, right, pi)
        if candidate.upper() > panel_bound:
            panel_bound = candidate.upper()
    require(panel_bound < arb("5.85e-8"), "independent beta4 second-derivative bound failed")
    saved_second = arb(artifact["height_cell_certificate"]["beta4_correction_second_derivative_bound_ball"])
    require(panel_bound <= saved_second.upper(), "independent beta4 second-derivative bound exceeds saved enclosure")

    exact_second = arb(cell["second_derivative_majorant"]["two_ray_second_derivative_bound"])
    uniform = abs(residual0) + radius * abs(residual1) + radius * radius * (exact_second + panel_bound) / 2
    require(uniform < arb("5.3e-15"), "independent uniform reduced bound failed")
    strip_scale = arb(2).sqrt() / pi * p["Y"]
    max_normalizer = (pi / (32 * (p["height"] - pi / 16))) ** (arb(1) / 4)
    require(max_normalizer * strip_scale * uniform < arb("4.9e-16"), "independent normalized bound failed")

    q = sp.symbols("q", nonzero=True)
    cp, cm = sp.symbols("C_plus C_minus")
    negative = sum(q ** (4 * j + 1) for j in range(199))
    positive = sum(q ** (-(4 * j + 3)) for j in range(200))
    selected = cm * negative + cp * positive
    opposite = cp * negative + cm * positive
    require(sp.expand(selected + opposite - (cp + cm) * (negative + positive)) == 0, "independent branch split failed")
    require(artifact["decision"]["zero_negative_outer_positive_complement_retained_as_one_object"] is True, "complement retention lost")
    require(artifact["decision"]["ordinary_Morse_corridor_splice_proved"] is False, "ordinary splice overpromoted")
    print(f"independently validated beta^-4 completed selector projection; priority={priority}", flush=True)


if __name__ == "__main__":
    main()
