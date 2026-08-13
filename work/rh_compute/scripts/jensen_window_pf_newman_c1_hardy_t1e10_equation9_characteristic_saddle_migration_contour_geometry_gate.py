#!/usr/bin/env python3
"""Certify the saddle-migration threshold and outgoing Airy contour rays."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


QUADRATURE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_saddle_migration_contour_geometry_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_saddle_migration_contour_geometry_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 110
T = 10_000_000_000
C = 159577
Y_MAX = 64
OLD_RADIUS = 4
NEW_RADIUS = 9


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


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


def symbolic_geometry() -> dict[str, Any]:
    radius, s, a = sp.symbols("R s a", real=True, positive=True)
    i = sp.I
    plus = radius + s * sp.exp(i * sp.pi / 6)
    minus = -radius + s * sp.exp(5 * i * sp.pi / 6)
    expected = (
        (radius**2 - a) * s / 2
        + sp.sqrt(3) * radius * s**2 / 2
        + s**3 / 3
    )
    plus_imaginary = sp.simplify(sp.im(plus**3 / 3 - a * plus))
    minus_imaginary = sp.simplify(sp.im(minus**3 / 3 - a * minus))
    require(sp.simplify(plus_imaginary - expected) == 0, "positive ray identity failed")
    require(sp.simplify(minus_imaginary - expected) == 0, "negative ray identity failed")
    return {
        "canonical_z_phase": "phi_a(z)=z^3/3-a*z, a=lambda+y",
        "stationary_points": "z_+-=+/-sqrt(a)",
        "positive_ray": "z=R+exp(i*pi/6)s, s>=0",
        "negative_ray": "z=-R+exp(5i*pi/6)s, s>=0",
        "ray_imaginary_phase": "Im phi_a=[R^2-a]s/2+[sqrt(3)R]s^2/2+s^3/3",
        "ray_modulus": "|exp(i phi_a)|=exp(-Im phi_a)",
        "deformation_guard": "For R^2>max(a), both real tails deform to these entire-function rays and decay immediately from s=0.",
    }


def certified_geometry() -> dict[str, Any]:
    ctx.dps = PRECISION
    t = arb(T)
    c = arb(C)
    pi = arb.pi()
    eta = pi * c**2 / (8 * t)
    beta = (t * eta) ** (arb(1) / 3)
    lam = (eta - 1) * t / beta
    a_max = lam + Y_MAX
    old_crossing = arb(OLD_RADIUS) ** 2 - lam
    new_crossing = arb(NEW_RADIUS) ** 2 - lam
    maximum_saddle = a_max.sqrt()
    saddle_clearance = arb(NEW_RADIUS) - maximum_saddle
    quadratic_margin = arb(NEW_RADIUS) ** 2 - a_max
    linear_decay = quadratic_margin / 2

    z = arb(NEW_RADIUS)
    y = arb(Y_MAX)
    endpoint_error = arb(64) / 15 * z**5 / beta**2
    coupling_error = y * z**3 / (3 * beta**2)
    quadratic_error = z * y**2 / (4 * beta**2)
    total_error = endpoint_error + coupling_error + quadratic_error
    relative_amplitude_upper = 4 * beta / (pi * c) * y / c
    relative_amplitude_lower = 1 - (z / beta).cosh() ** (-arb(3) / 2)
    relative_amplitude_error = max(relative_amplitude_upper, relative_amplitude_lower)

    require(arb(10) < old_crossing < arb(11), "old saddle crossing drift")
    require(old_crossing < arb(Y_MAX), "old radius unexpectedly contains all saddles")
    require(new_crossing > arb(Y_MAX), "new radius does not contain all saddles")
    require(maximum_saddle < arb(NEW_RADIUS), "maximum saddle escaped new radius")
    require(saddle_clearance > arb("0.68"), "saddle clearance below 0.68")
    require(quadratic_margin > arb("11.8"), "quadratic ray margin below 11.8")
    require(linear_decay > arb("5.9"), "ray decay coefficient below 5.9")
    require(total_error < arb("0.06"), "expanded-core phase error exceeds 0.06")
    require(relative_amplitude_error < arb("0.000014"), "expanded-core amplitude error exceeds 1.4e-5")

    return {
        "height": T,
        "y_range": [0, Y_MAX],
        "old_radius": OLD_RADIUS,
        "new_radius": NEW_RADIUS,
        "lambda_ball": lam.str(PRECISION, more=True),
        "old_radius_saddle_crossing_y_ball": old_crossing.str(PRECISION, more=True),
        "new_radius_saddle_crossing_y_ball": new_crossing.str(PRECISION, more=True),
        "maximum_saddle_radius_ball": maximum_saddle.str(PRECISION, more=True),
        "new_radius_saddle_clearance_ball": saddle_clearance.str(PRECISION, more=True),
        "outgoing_ray_quadratic_margin_ball": quadratic_margin.str(PRECISION, more=True),
        "outgoing_ray_linear_decay_coefficient_ball": linear_decay.str(PRECISION, more=True),
        "expanded_core_endpoint_phase_error_ball": endpoint_error.str(PRECISION, more=True),
        "expanded_core_coupling_phase_error_ball": coupling_error.str(PRECISION, more=True),
        "expanded_core_quadratic_phase_error_ball": quadratic_error.str(PRECISION, more=True),
        "expanded_core_total_phase_error_ball": total_error.str(PRECISION, more=True),
        "expanded_core_relative_amplitude_error_ball": relative_amplitude_error.str(PRECISION, more=True),
        "decision": "The measured |z|>4 contribution mixes a saddle-migration connector with a true tail. Radius 9 contains every Airy saddle for 0<=y<=64; only |z|>9 should be treated by the outgoing pi/6 and 5pi/6 decay rays.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_geometry"]
    return f"""# Characteristic saddle migration and contour geometry

Date: 2026-08-11

Status: exact contour geometry validated; not a proof of the contour-tail
integral bound or the outer saddle join

For fixed `y`, the canonical `z` phase is

```text
phi_a(z)=z^3/3-a*z,       a=lambda+y,
z_+-=+/-sqrt(a).                                           (CG1)
```

At the previous radius `R=4`, a saddle crosses the boundary when

```text
y=16-lambda
 ={c['old_radius_saddle_crossing_y_ball']}.              (CG2)
```

Because (CG2) lies inside `0<=y<=64`, the measured `|z|>4` contribution is
not purely a nonstationary tail.  It contains the stationary connector for
larger `y`.

Choose `R=9`.  The largest saddle on this normal range is

```text
sqrt(lambda+64)={c['maximum_saddle_radius_ball']}<9,      (CG3)
```

with clearance `{c['new_radius_saddle_clearance_ball']}>0.68`.

The outgoing rays are

```text
z=9+exp(i*pi/6)s,
z=-9+exp(5i*pi/6)s,       s>=0.                          (CG4)
```

and exact expansion on either ray gives

```text
Im phi_a(z)
 =[81-a]s/2+(9sqrt(3))s^2/2+s^3/3.                      (CG5)
```

Uniformly for `0<=y<=64`,

```text
81-a >= {c['outgoing_ray_quadratic_margin_ball']},
(81-a)/2 >= {c['outgoing_ray_linear_decay_coefficient_ball']}>5.9. (CG6)
```

Thus both rays decay immediately and cubically at infinity.  The integrand is
entire in `z`, so the real tails outside `9` may be deformed to (CG4); the
connector `4<|z|<=9` must instead remain part of the stationary core.

The exact finite-`t` comparison can still be continued on the expanded
rectangle.  The conservative phase and relative-amplitude errors obey

```text
|R_phase| <= {c['expanded_core_total_phase_error_ball']} <0.06,
|A/(sqrt(2)/pi)-1|
 <= {c['expanded_core_relative_amplitude_error_ball']} <1.4e-5. (CG7)
```

These are admissibility bounds, not a replacement for direct
cancellation-preserving quadrature.

Pi provenance: the ray angles are the cubic Airy steepest-descent angles,
while `lambda` and the finite-`t` amplitude retain the Kummer and integer
Fourier--Poisson normalization already traced.  No geometric fit is used.

Proof boundary: exact contour algebra and source-height interval geometry
only.  No connector quadrature, contour-tail integral constant, `y>64`
partition, complete `T_upper`, `Lambda<=0`, RH, or prize-level conclusion is
proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    require(QUADRATURE_GATE.is_file() and CHECKER.is_file(), "missing dependency or checker")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_saddle_migration_contour_geometry_gate",
        "status": "characteristic_saddle_migration_and_R9_outgoing_contour_geometry_complete",
        "passed": True,
        "symbolic_geometry": symbolic_geometry(),
        "certified_geometry": certified_geometry(),
        "decision": {
            "old_abs_z_gt_4_term_is_pure_nonstationary_tail": False,
            "saddle_migration_threshold_certified": True,
            "R9_contains_all_y_le_64_saddles": True,
            "R9_outgoing_rays_decay_immediately": True,
            "connector_quadrature_proved": False,
            "contour_tail_integral_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Integrate the canonical and exact finite-t connector 4<|z|<=9 with resumable panels, subtract from the full-line value to measure the true R=9 contour tail, and then certify that ray integral directly.",
        "proof_boundary": "Exact contour algebra and source-height interval geometry only. No connector quadrature, contour-tail integral constant, y>64 partition, complete T_upper, Lambda<=0, RH, or prize-level conclusion is proved.",
        "dependencies": {
            "quadrature_gate": {"path": relative(QUADRATURE_GATE), "sha256": file_hash(QUADRATURE_GATE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("built characteristic contour geometry: crossing=10.8901, R9 margin=11.8901")


if __name__ == "__main__":
    main()
