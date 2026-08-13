#!/usr/bin/env python3
"""Independently check the endpoint-safe exterior B trace extraction."""

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

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_tangential_exterior_cutoff_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    u = sp.symbols("u", real=True)
    smoothstep = 6 * u**5 - 15 * u**4 + 10 * u**3
    first = sp.diff(smoothstep, u)
    second = sp.diff(first, u)
    require([smoothstep.subs(u, k) for k in (0, 1)] == [0, 1], "smoothstep endpoint drift")
    require(all(first.subs(u, k) == 0 and second.subs(u, k) == 0 for k in (0, 1)), "C2 join drift")
    require(sp.simplify(first.subs(u, sp.Rational(1, 2)) - sp.Rational(15, 8)) == 0, "first derivative maximum drift")
    critical = (sp.Integer(3) - sp.sqrt(3)) / 6
    require(sp.simplify(second.subs(u, critical) - 10 * sp.sqrt(3) / 3) == 0, "second derivative maximum drift")

    # Independent finite-scalar replay of the three-part allocation.
    rows = [
        (sp.Rational(7, 3), sp.Rational(-5, 4), 0, 0),
        (sp.Rational(-11, 5), sp.Rational(9, 7), 1, 0),
        (sp.Rational(13, 8), sp.Rational(-3, 2), 0, sp.Rational(2, 5)),
    ]
    for global_value, B_value, window, cutoff in rows:
        B_window = window * B_value
        B_tangent = cutoff * (1 - window) * B_value
        joined = global_value - B_window - B_tangent
        require(B_window + B_tangent + joined == global_value, "finite allocation replay failed")

    ctx.dps = 90
    ctx.threads = 1
    t = arb(10_000_000_000)
    endpoint = arb(5_122_421)
    pi = arb.pi()
    delta = arb(1) / 10_000
    x_trace = (1 - (1 - 8 * t / (pi * endpoint**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x_trace) / (2 * x_trace**2 * (1 - x_trace) ** 2)
    root_h = hessian.sqrt()
    x_low = x_trace - arb(70) / root_h
    x_high = x_trace + arb(70) / root_h

    def G(x: arb) -> arb:
        return (pi * endpoint**2 / 4 - t / (2 * x * (1 - x))) / root_h

    require(2 * delta < x_low, "independent cutoff/window separation failed")
    require(-G(x_low) > arb("70.06"), "independent left gap failed")
    require(G(x_high) > arb("69.92"), "independent right gap failed")
    require((arb(15) / (8 * delta)) / root_h < arb("6.5e-5"), "independent first cutoff derivative failed")
    require((arb(10) * arb(3).sqrt() / (3 * delta**2)) / hessian < arb("7e-9"), "independent second cutoff derivative failed")

    decision = artifact["decision"]
    require(decision["Abel_endpoint_difference_retained_in_joined_remainder"] is True, "endpoint grouping lost")
    require(decision["separate_exterior_and_joined_Abel_limits_proved"] is False, "separate Abel limits overclaimed")
    require(decision["exterior_B_trace_tangential_bound_proved"] is False, "exterior trace bound overclaim")
    require(decision["joined_remainder_bound_proved"] is False, "joined-remainder overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently checked C2 corner cutoff, exterior tangential gaps, and global allocation", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
