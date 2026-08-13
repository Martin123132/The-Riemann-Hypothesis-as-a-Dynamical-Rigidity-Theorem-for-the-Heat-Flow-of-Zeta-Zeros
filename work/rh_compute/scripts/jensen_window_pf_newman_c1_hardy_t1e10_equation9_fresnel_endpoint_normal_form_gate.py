#!/usr/bin/env python3
"""Certify the exact Fresnel endpoint normal form for finite-Poisson modes."""

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


INTERCHANGE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json"
SADDLE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fresnel_endpoint_normal_form_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fresnel_endpoint_normal_form_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 110
T = 10_000_000_000
A = 159577
B = 5122421


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


def symbolic_reduction() -> dict[str, Any]:
    alpha, x, m, t, C, z = sp.symbols("alpha x m t C z", positive=True)
    pi = sp.pi
    z_alpha = sp.sqrt(x / 2) * (alpha - 2 * m / x)
    completed_phase = x * (alpha - 2 * m / x) ** 2 / 4
    require(sp.simplify(completed_phase - z_alpha**2 / 2) == 0, "Fresnel scaling failed")
    require(sp.simplify(sp.diff(z_alpha, alpha) - sp.sqrt(x / 2)) == 0, "Fresnel Jacobian failed")

    k = sp.exp(-sp.I * pi / 4) * sp.sqrt(pi / 2)
    fresnel = sp.exp(sp.I * pi / 4) * sp.erf(k * z) / sp.sqrt(2)
    require(sp.simplify(sp.diff(fresnel, z) - sp.exp(sp.I * pi * z**2 / 2)) == 0, "Fresnel primitive failed")

    x_star = 2 * pi * m**2 / (t + 2 * pi * m**2)
    alpha_star = 2 * m + t / (pi * m)
    z_c = sp.sqrt(x / 2) * (C - 2 * m / x)
    z_c_star = sp.factor(sp.simplify(z_c.subs(x, x_star)))
    expected_z = m * sp.sqrt(pi / (t + 2 * pi * m**2)) * (C - alpha_star)
    require(sp.simplify(z_c_star - expected_z) == 0, "stationary endpoint coordinate failed")

    psi = -pi * m**2 / x + t * sp.log((1 - x) / x) / 2
    psi1 = sp.simplify(sp.diff(psi, x).subs(x, x_star))
    psi2 = sp.factor(sp.simplify(sp.diff(psi, x, 2).subs(x, x_star)))
    expected_psi2 = -(t + 2 * pi * m**2) ** 4 / (8 * pi**2 * m**4 * t)
    require(psi1 == 0, "reduced x saddle failed")
    require(sp.simplify(psi2 - expected_psi2) == 0, "reduced x curvature failed")

    return {
        "normalized_endpoints": "z_C(m,x)=sqrt(x/2)(C-2m/x), C in {A,B}",
        "fresnel_primitive": "F(z)=exp(i*pi/4)/sqrt(2)*erf(exp(-i*pi/4)*sqrt(pi/2)*z)",
        "fresnel_derivative": "F'(z)=exp(i*pi*z^2/2)",
        "exact_mode_fresnel": "F_m(A,B;x)=exp(-i*pi*m^2/x)*sqrt(2/x)[F(z_B)-F(z_A)]",
        "mode_current": "J_m=(-1)^m/x{[E_m(B)-E_m(A)]/(i*pi)+m*F_m(A,B;x)}",
        "stationary_endpoint_coordinate": "z_C*=m*sqrt(pi/(t+2*pi*m^2))[C-2m-t/(pi*m)]",
        "reduced_x_phase": "psi_m(x)=-pi*m^2/x+(t/2)log((1-x)/x)",
        "reduced_x_saddle": "x_m=2*pi*m^2/(t+2*pi*m^2)",
        "reduced_x_curvature": str(psi2),
        "curvature_sign": "strictly_negative for m,t>0",
        "fresnel_infinity": "F(+infinity)=(1+i)/2 under Abel damping",
        "positive_tail_bound": "|F(+infinity)-F(z)|<=2/(pi*z) for z>0",
        "turning_interpretation": "The alpha integral is a smooth Fresnel endpoint transition and the reduced x saddle does not degenerate. This alone does not exclude a characteristic boundary tangent or composite fold in the coupled x/Fresnel integral.",
    }


def z_star(mode: int, endpoint: int) -> arb:
    m = arb(mode)
    t = arb(T)
    pi = arb.pi()
    return m * (pi / (t + 2 * pi * m**2)).sqrt() * (arb(endpoint) - 2 * m - t / (pi * m))


def certified_atlas() -> dict[str, Any]:
    ctx.dps = PRECISION
    gap = [(mode, z_star(mode, A), z_star(mode, B)) for mode in range(39853, 39937)]
    min_a = min(gap, key=lambda row: row[1])
    max_a = max(gap, key=lambda row: row[1])
    min_b = min(gap, key=lambda row: row[2])
    max_b = max(gap, key=lambda row: row[2])
    half_saddle_error = max_a[1] + 2 / (arb.pi() * min_b[2])
    require(min_a[1] > 0, "A-gap endpoint coordinate is not positive")
    require(max_a[1] < arb("0.044"), "A-gap endpoint coordinate exceeds compact chart")
    require(min_b[2] > arb("2480000"), "far endpoint is not asymptotic")
    require(half_saddle_error < arb("0.044"), "half-saddle error exceeds target chart")

    key_values = {
        "zA_mode_39852": z_star(39852, A),
        "zA_mode_39937": z_star(39937, A),
        "zB_mode_621": z_star(621, B),
        "zB_mode_622": z_star(622, B),
        "zB_mode_2560588": z_star(2560588, B),
        "zB_mode_2560589": z_star(2560589, B),
    }
    require(key_values["zA_mode_39852"] < 0 and key_values["zA_mode_39937"] < 0, "A interior neighbor signs failed")
    require(key_values["zB_mode_621"] < 0 < key_values["zB_mode_622"], "lower B transition signs failed")
    require(key_values["zB_mode_2560589"] < 0 < key_values["zB_mode_2560588"], "upper B transition signs failed")

    return {
        "precision_decimal_digits": PRECISION,
        "positive_mode_partition": {
            "lower_nonstationary_tail": "1..621: z_A*<z_B*<0",
            "lower_interior": "622..39852: z_A*<0<z_B*",
            "lower_A_endpoint_transition": "39853..39894: 0<z_A*<z_B*",
            "upper_A_endpoint_transition": "39895..39936: 0<z_A*<z_B*",
            "upper_interior": "39937..2560588: z_A*<0<z_B*",
            "upper_nonstationary_tail": "2560589..infinity: z_A*<z_B*<0",
            "zero_and_negative_modes": "no positive-alpha stationary point",
        },
        "A_transition_mode_count": len(gap),
        "zA_transition_min_mode": min_a[0],
        "zA_transition_min_ball": min_a[1].str(PRECISION, more=True),
        "zA_transition_max_mode": max_a[0],
        "zA_transition_max_ball": max_a[1].str(PRECISION, more=True),
        "zB_transition_min_mode": min_b[0],
        "zB_transition_min_ball": min_b[2].str(PRECISION, more=True),
        "zB_transition_max_mode": max_b[0],
        "zB_transition_max_ball": max_b[2].str(PRECISION, more=True),
        "half_saddle_factor_error_bound_ball": half_saddle_error.str(PRECISION, more=True),
        "key_neighbor_balls": {name: value.str(PRECISION, more=True) for name, value in key_values.items()},
        "decision": "All 84 A-endpoint modes lie in one compact Fresnel chart and their saddle-point Fresnel values are within 0.044 of (1+i)/2. This resolves the alpha integral itself but does not decide whether the coupled x/Fresnel integral requires a fold-uniform normal form.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    a = artifact["certified_atlas"]
    s = artifact["symbolic_reduction"]
    return f"""# Exact Fresnel endpoint normal form for the portcullis modes

Date: 2026-08-10

Status: exact endpoint-uniform alpha normal form validated; not a proof of the quantitative x-stationary remainder

For each positive Poisson mode, complete the square in the alpha integral and
introduce

```text
z_C(m,x)=sqrt(x/2)(C-2m/x), C in {{A,B}},
F(z)=exp(i*pi/4)/sqrt(2)
     erf(exp(-i*pi/4)*sqrt(pi/2)z).
```

Then `F'(z)=exp(i*pi*z^2/2)` and the mode is exactly

```text
F_m(A,B;x)=exp(-i*pi*m^2/x)*sqrt(2/x)
            [F(z_B)-F(z_A)],                            (FN1)

J_m=(-1)^m/x{{[E_m(B)-E_m(A)]/(i*pi)+m F_m(A,B;x)}}.    (FN2)
```

No asymptotic expansion has been used in (FN1)--(FN2).  At the reduced
Kummer saddle

```text
x_m=2*pi*m^2/(t+2*pi*m^2),
z_C*=m sqrt(pi/(t+2*pi*m^2))[C-2m-t/(pi*m)].             (FN3)
```

The remaining x phase is

```text
psi_m(x)=-pi*m^2/x+(t/2)log((1-x)/x),
psi_m''(x_m)={s['reduced_x_curvature']}<0.               (FN4)
```

Thus the alpha integral is a smooth Fresnel endpoint transition, while the
reduced x saddle remains nondegenerate.  These two facts do not by themselves
exclude a characteristic boundary tangent in the coupled x/Fresnel integral;
that composite question must be tested separately.

At `t=10^10`, all 84 modes `39853..39936` lie in one certified compact chart:

```text
{a['zA_transition_min_ball']} <= z_A* <=
{a['zA_transition_max_ball']},

z_B* >= {a['zB_transition_min_ball']}.
```

Using `F(+infinity)=(1+i)/2` and
`|F(+infinity)-F(z)|<=2/(pi*z)` for `z>0`, every one of these
84 factors satisfies

```text
|[F(z_B*)-F(z_A*)]-(1+i)/2|
 <= {a['half_saddle_factor_error_bound_ball']} < 0.044. (FN5)
```

This explains the 42 classical lower endpoint terms: they are nearly
half-saddle Fresnel contributions, not failed or divergent saddles.  The exact
positive-mode partition is

```text
1..621             lower nonstationary B-tail,
622..39852         lower interior,
39853..39894       lower A-endpoint transition,
39895..39936       reflected A-endpoint transition,
39937..2560588     upper interior,
2560589..infinity  upper nonstationary B-tail.
```

The zero and negative modes have no positive-alpha stationary point.  The
next theorem must estimate the x integral in these six positive ranges plus
the nonpositive modes, keeping the boundary and Fresnel pieces in (FN2)
grouped.

Proof boundary: exact endpoint-normal-form algebra and a certified finite
transition atlas at `t=10^10`.  No explicit x-stationary remainder,
height-uniform source error, `Lambda<=0`, RH, or prize-level conclusion is
proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    require(INTERCHANGE_GATE.is_file() and SADDLE_GATE.is_file() and CHECKER.is_file(), "missing dependency or checker")
    symbolic = symbolic_reduction()
    atlas = certified_atlas()
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_fresnel_endpoint_normal_form_gate",
        "status": "exact_fresnel_endpoint_normal_form_and_six_range_mode_partition_complete",
        "passed": True,
        "scope": {"height_center": "1e10", "source_roster": f"{A}..{B} odd", "diagnostic_midpoint_used": False},
        "symbolic_reduction": symbolic,
        "certified_atlas": atlas,
        "decision": {
            "exact_alpha_fresnel_normal_form_proved": True,
            "all_84_A_endpoint_modes_share_one_compact_chart": True,
            "alpha_integral_Airy_normal_form_required": False,
            "composite_x_endpoint_fold_uniformization_required": "open_near_characteristic",
            "reduced_x_saddle_nondegenerate": True,
            "positive_mode_partition_complete": True,
            "x_stationary_remainder_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Prove explicit x-stationary and nonstationary estimates across the six positive-mode ranges and nonpositive modes, retaining the exact grouped current J_m and summing before absolute values wherever required by convergence.",
        "proof_boundary": "Exact Fresnel alpha normal form and certified finite endpoint atlas at t=10^10 only. No explicit x-stationary remainder, height-uniform source error, Lambda<=0, RH, or prize-level conclusion is proved.",
        "dependencies": {
            "interchange_gate": {"path": relative(INTERCHANGE_GATE), "sha256": file_hash(INTERCHANGE_GATE)},
            "saddle_gate": {"path": relative(SADDLE_GATE), "sha256": file_hash(SADDLE_GATE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {"workers": 1, "sympy_threads": 1, "flint_threads": 1, "process_priority": priority, "elapsed_seconds": round(time.perf_counter() - started, 3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("built Fresnel endpoint normal form: modes=84, chart<0.044, ranges=6")


if __name__ == "__main__":
    main()
