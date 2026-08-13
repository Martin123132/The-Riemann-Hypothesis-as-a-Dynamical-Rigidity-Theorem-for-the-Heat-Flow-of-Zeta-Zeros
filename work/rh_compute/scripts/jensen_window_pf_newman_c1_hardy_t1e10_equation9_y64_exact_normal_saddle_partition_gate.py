#!/usr/bin/env python3
"""Certify the exact y=64 endpoint/Fresnel/Morse saddle partition."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


NORMAL_FORM_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate.json"
RAY_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_R9_outgoing_ray_quadrature_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_exact_normal_saddle_partition_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_exact_normal_saddle_partition_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 100
T = 10_000_000_000
C = 159577
B = 5122421
MODE_LO = 39853
MODE_HI = 39936
Y0 = 64
RADIUS = 9
FRESNEL_WIDTHS = 4


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


def atanh_ball(value: arb) -> arb:
    return ((1 + value) / (1 - value)).log() / 2


def symbolic_partition() -> dict[str, str]:
    y, beta, d, tau, x, y0 = sp.symbols("y beta d tau x Y0", positive=True, real=True)
    phase = -(d + beta * tau) * y + x * y**2 / (2 * beta)
    saddle = beta * (d + beta * tau) / x
    completed = x * (y - saddle) ** 2 / (2 * beta) - beta * (d + beta * tau) ** 2 / (2 * x)
    require(sp.simplify(phase - completed) == 0, "normal square completion failed")

    g = d + beta * tau - x * y0 / beta
    boundary_q = sp.simplify((y0 - saddle) * sp.sqrt(x / beta))
    require(sp.simplify(boundary_q + g * sp.sqrt(beta / x)) == 0, "boundary coordinate failed")

    return {
        "exact_y_phase": "-(d_m+beta*tanh(z/beta))*y+x(z)*y^2/(2beta)",
        "positive_curvature": "partial_y^2 Theta=x(z)/beta>0",
        "exact_saddle": "y_m(z)=beta[d_m+beta*tanh(z/beta)]/x(z)",
        "boundary_function": "g_m(z)=d_m+beta*tanh(z/beta)-Y0*x(z)/beta",
        "boundary_derivative": "g_m'(z)=sech(z/beta)^2[1+Y0/(2beta^2)]>0",
        "Fresnel_coordinate": "q_m(y,z)=(y-y_m(z))*sqrt(x(z)/beta)",
        "boundary_coordinate": "q_m(Y0,z)=-g_m(z)*sqrt(beta/x(z))",
        "closed_root": "z_m(c)=beta*atanh([c-d_m+Y0/(2beta)]/[beta+Y0/(2beta)])",
        "partition": "g<=-delta: nonstationary; -delta<g<delta: endpoint/Fresnel; g>=delta: ordinary Morse",
    }


def certified_partition() -> dict[str, Any]:
    ctx.dps = PRECISION
    t = arb(T)
    c = arb(C)
    pi = arb.pi()
    eta = pi * c**2 / (8 * t)
    beta = (t * eta) ** (arb(1) / 3)
    sigma = 4 * beta / (pi * c)
    h = 4 * beta / c
    detunings = [h * (arb(mode) - c / 4) for mode in range(MODE_LO, MODE_HI + 1)]

    def tau(z: arb) -> arb:
        return (z / beta).tanh()

    def x(z: arb) -> arb:
        return (1 - tau(z)) / 2

    x_max = x(arb(-RADIUS))
    x_min = x(arb(RADIUS))
    delta = arb(FRESNEL_WIDTHS) * (x_max / beta).sqrt()
    denominator = beta + arb(Y0) / (2 * beta)

    def root(level: arb, detuning: arb) -> arb:
        argument = (level - detuning + arb(Y0) / (2 * beta)) / denominator
        require(abs(argument) < arb(1), "closed-root argument left (-1,1)")
        return beta * atanh_ball(argument)

    def saddle(z: arb, detuning: arb) -> arb:
        return beta * (detuning + beta * tau(z)) / x(z)

    left_roots = [root(-delta, detuning) for detuning in detunings]
    center_roots = [root(arb(0), detuning) for detuning in detunings]
    right_roots = [root(delta, detuning) for detuning in detunings]
    canonical_roots = [arb(Y0) / (2 * beta) - detuning for detuning in detunings]
    threshold_shifts = [abs(exact - canonical) for exact, canonical in zip(center_roots, canonical_roots)]
    band_widths = [right - left for left, right in zip(left_roots, right_roots)]

    normalized_left = [delta * (beta / x(z)).sqrt() for z in left_roots]
    normalized_right = [delta * (beta / x(z)).sqrt() for z in right_roots]
    min_normalized_cut = min(normalized_left + normalized_right)

    y_span = arb(B - C) / sigma
    max_saddle = saddle(arb(RADIUS), detunings[-1])
    upper_normalized_clearance = (y_span - max_saddle) * (x_min / beta).sqrt()

    require(max(left_roots) < max(center_roots) < max(right_roots), "upper mode partition ordering failed")
    require(min(left_roots) < min(center_roots) < min(right_roots), "lower mode partition ordering failed")
    require(min(left_roots) > arb(-RADIUS) and max(right_roots) < arb(RADIUS), "transition band escaped R=9")
    require(min_normalized_cut >= arb(FRESNEL_WIDTHS), "four-width cut failed")
    require(max(threshold_shifts) < arb("0.00002"), "canonical threshold shift exceeds 2e-5")
    require(max_saddle < arb("50000"), "normal saddle exceeds 50000")
    require(upper_normalized_clearance > arb("1000000"), "upper endpoint is not uniformly remote")
    require(y_span > arb("100000000"), "full y span drift")

    roster = []
    for offset, mode in enumerate(range(MODE_LO, MODE_HI + 1)):
        roster.append(
            {
                "mode": mode,
                "detuning_ball": detunings[offset].str(PRECISION, more=True),
                "nonstationary_to_Fresnel_z_ball": left_roots[offset].str(PRECISION, more=True),
                "saddle_crossing_z_ball": center_roots[offset].str(PRECISION, more=True),
                "Fresnel_to_Morse_z_ball": right_roots[offset].str(PRECISION, more=True),
                "canonical_saddle_crossing_z_ball": canonical_roots[offset].str(PRECISION, more=True),
            }
        )

    return {
        "height": T,
        "lower_endpoint": C,
        "upper_endpoint": B,
        "mode_range": [MODE_LO, MODE_HI],
        "mode_count": len(detunings),
        "z_working_range": [-RADIUS, RADIUS],
        "outer_y_range": [Y0, "(B-C)/sigma"],
        "beta_ball": beta.str(PRECISION, more=True),
        "sigma_ball": sigma.str(PRECISION, more=True),
        "full_y_span_ball": y_span.str(PRECISION, more=True),
        "x_minus_R_ball": x_max.str(PRECISION, more=True),
        "x_plus_R_ball": x_min.str(PRECISION, more=True),
        "uniform_g_buffer_delta_ball": delta.str(PRECISION, more=True),
        "minimum_normalized_boundary_cut_ball": min_normalized_cut.str(PRECISION, more=True),
        "left_cut_global_min_ball": min(left_roots).str(PRECISION, more=True),
        "left_cut_global_max_ball": max(left_roots).str(PRECISION, more=True),
        "crossing_global_min_ball": min(center_roots).str(PRECISION, more=True),
        "crossing_global_max_ball": max(center_roots).str(PRECISION, more=True),
        "right_cut_global_min_ball": min(right_roots).str(PRECISION, more=True),
        "right_cut_global_max_ball": max(right_roots).str(PRECISION, more=True),
        "minimum_transition_band_width_ball": min(band_widths).str(PRECISION, more=True),
        "maximum_transition_band_width_ball": max(band_widths).str(PRECISION, more=True),
        "maximum_exact_minus_canonical_crossing_ball": max(threshold_shifts).str(PRECISION, more=True),
        "maximum_normal_saddle_on_R9_ball": max_saddle.str(PRECISION, more=True),
        "minimum_upper_endpoint_normalized_clearance_ball": upper_normalized_clearance.str(PRECISION, more=True),
        "roster": roster,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_partition"]
    return f"""# Exact y=64 normal-saddle partition

Date: 2026-08-11

Status: exact saved-height endpoint/Fresnel/Morse partition validated; not a
proof of uniform estimates for the three resulting integrals

For a transition mode `m`, the complete finite-`t` dependence on the normal
coordinate `y` is

```text
Theta_m(z,y)=constant-(d_m+beta*tanh(z/beta))*y
             +x(z)y^2/(2beta),
x(z)=[1-tanh(z/beta)]/2.                                (Y64P1)
```

It has positive curvature `x(z)/beta` and the unique normal saddle

```text
y_m(z)=beta[d_m+beta*tanh(z/beta)]/x(z).                (Y64P2)
```

At `Y0=64`, define

```text
g_m(z)=d_m+beta*tanh(z/beta)-64x(z)/beta.               (Y64P3)
```

Then `y_m(z)>=64` exactly when `g_m(z)>=0`, and

```text
g_m'(z)=sech(z/beta)^2[1+32/beta^2]>0.                 (Y64P4)
```

Thus every mode has one crossing.  More generally `g_m(z)=c` has the closed
solution

```text
z_m(c)=beta*atanh([c-d_m+32/beta]/[beta+32/beta]).      (Y64P5)
```

No fitted root finder is used.

The normalized Fresnel coordinate at the boundary is

```text
q_m(64,z)=-g_m(z)sqrt(beta/x(z)).                       (Y64P6)
```

Choose the common buffer

```text
delta={c['uniform_g_buffer_delta_ball']}.               (Y64P7)
```

It guarantees `|q_m(64,z)|>=4` whenever `|g_m(z)|>=delta` throughout
`|z|<=9`.  The exact disjoint partition is therefore

```text
g_m(z)<=-delta:             nonstationary outer chart,
-delta<g_m(z)<delta:        endpoint/Fresnel chart,
g_m(z)>=delta:              ordinary Morse chart.       (Y64P8)
```

Assign equality to the two outer charts as displayed.  The 84 modewise bands
all lie inside

```text
{c['left_cut_global_min_ball']} <= z <=
{c['right_cut_global_max_ball']}                              (Y64P9)
```

and hence well inside `|z|<9`.  The largest exact normal saddle on the R=9
rectangle is

```text
{c['maximum_normal_saddle_on_R9_ball']} <50000,         (Y64P10)
```

whereas the upper endpoint is at

```text
y_B=(B-C)/sigma={c['full_y_span_ball']}>1e8.            (Y64P11)
```

Its minimum normalized clearance exceeds one million, so no transition-mode
saddle in `|z|<=9` can be assigned accidentally to the upper endpoint chart.
The exact crossing differs from the canonical threshold
`64/(2beta)-d_m` by less than `2e-5` for every mode.

Pi provenance: `beta`, `sigma`, and `d_m` retain the Kummer and integer
Fourier--Poisson normalization of Section 11.314.  The four-width buffer is a
dimensionless Fresnel-coordinate choice, not a geometric fit.

Proof boundary: exact finite-`t` normal square completion and a disjoint
saved-height chart assignment only.  No uniform integral estimate on any
chart, join to all ordinary interior modes, complete `T_upper`, `Lambda<=0`,
RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    require(NORMAL_FORM_GATE.is_file() and RAY_GATE.is_file() and CHECKER.is_file(), "missing dependency or checker")
    symbolic = symbolic_partition()
    certified = certified_partition()
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_exact_normal_saddle_partition_gate",
        "status": "exact_y64_nonstationary_Fresnel_Morse_partition_complete",
        "passed": True,
        "symbolic_partition": symbolic,
        "certified_partition": certified,
        "decision": {
            "exact_normal_square_completion_proved": True,
            "unique_modewise_y64_crossing_proved": True,
            "closed_crossing_formula_proved": True,
            "three_chart_partition_disjoint_and_exhaustive": True,
            "every_displaced_R9_normal_saddle_assigned": True,
            "upper_endpoint_excluded_for_transition_saddles": True,
            "three_chart_integral_estimates_proved": False,
            "ordinary_interior_mode_join_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Derive explicit nonstationary integration-by-parts, uniform incomplete-Fresnel, and ordinary Morse remainder constants on the three charts, then match the Morse main phase and amplitude to the 39231 lower interior modes without double counting.",
        "proof_boundary": "Exact finite-t normal square completion and a disjoint saved-height chart assignment only. No uniform integral estimate on any chart, join to all ordinary interior modes, complete T_upper, Lambda<=0, RH, or prize-level conclusion is proved.",
        "dependencies": {
            "normal_form_gate": {"path": relative(NORMAL_FORM_GATE), "sha256": file_hash(NORMAL_FORM_GATE)},
            "ray_gate": {"path": relative(RAY_GATE), "sha256": file_hash(RAY_GATE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {"workers": 1, "process_priority": priority, "elapsed_seconds": round(time.perf_counter() - started, 3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("built exact y=64 normal-saddle partition: 84 modes, three disjoint charts")


if __name__ == "__main__":
    main()
