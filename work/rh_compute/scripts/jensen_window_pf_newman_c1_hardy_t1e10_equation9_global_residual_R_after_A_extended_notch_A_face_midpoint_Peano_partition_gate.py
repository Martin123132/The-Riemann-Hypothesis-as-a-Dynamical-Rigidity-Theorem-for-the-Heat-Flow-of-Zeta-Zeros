#!/usr/bin/env python3
"""Certify the signed A-face midpoint-Peano reduction and x partition."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_midpoint_Peano_partition_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
PILOT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_A_face_42_mode_midpoint_defect_pilot.json"

DEPENDENCIES = {
    "finite_regulator_equivalence": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_finite_regulator_equivalence_gate.json",
    "normal_tangential_split": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_tangential_split_gate.json",
    "two_current_expansion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_two_current_expansion_gate.json",
}

A = 159_577
T = 10_000_000_000
FIRST_MODE = 39_895
LAST_MODE = 39_936
LEFT_NUMERATOR = 79_789
RIGHT_NUMERATOR = 79_873
Q_CUT = 16
PRECISION = 100


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


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


def symbolic_certificate() -> dict[str, Any]:
    left = sp.Rational(LEFT_NUMERATOR, 2)
    right = sp.Rational(RIGHT_NUMERATOR, 2)
    modes = tuple(range(FIRST_MODE, LAST_MODE + 1))
    require(len(modes) == 42, "shifted A-block cardinality drift")
    require(right - left == 42, "continuous strip length drift")
    require(all(sp.Rational(m, 1) == left + sp.Rational(2 * j + 1, 2) for j, m in enumerate(modes)), "midpoint roster drift")

    z = sp.symbols("z", real=True)
    half = sp.Rational(1, 2)
    kernel_left = sp.Rational(1, 2) * (half + z) ** 2
    kernel_right = sp.Rational(1, 2) * (half - z) ** 2
    kernel_mass = sp.integrate(kernel_left, (z, -half, 0)) + sp.integrate(kernel_right, (z, 0, half))
    require(kernel_mass == sp.Rational(1, 24), "one-cell Peano mass drift")
    require(42 * kernel_mass == sp.Rational(7, 4), "42-cell Peano mass drift")

    # The Peano formula is independently pinned on the complete cubic test
    # space; constants and linears vanish, while quadratics and cubics fix the
    # kernel sign and normalization.
    a0, a1, a2, a3 = sp.symbols("a0 a1 a2 a3")
    polynomial = a0 + a1 * z + a2 * z**2 + a3 * z**3
    midpoint_error = polynomial.subs(z, 0) - sp.integrate(polynomial, (z, -half, half))
    peano_error = -sp.integrate(kernel_left * sp.diff(polynomial, z, 2), (z, -half, 0)) - sp.integrate(kernel_right * sp.diff(polynomial, z, 2), (z, 0, half))
    require(sp.simplify(midpoint_error - peano_error) == 0, "Peano sign or normalization drift")

    q = sp.symbols("q", real=True)
    h = sp.Function("H")(q)
    h_q = 1 - sp.I * sp.pi * q * h
    h_qq = sp.diff(h_q, q).subs(sp.diff(h, q), h_q)
    expected_h_qq = -sp.I * sp.pi * q - (sp.I * sp.pi + sp.pi**2 * q**2) * h
    require(sp.simplify(h_qq - expected_h_qq) == 0, "phase-stripped Fresnel ODE drift")

    y, c = sp.symbols("y c", positive=True, real=True)
    h_formal = sp.Function("H")(q)

    def d_y(expression: sp.Expr) -> sp.Expr:
        return sp.diff(expression, y) - c * sp.diff(expression, q)

    g = -1 / (sp.I * sp.pi) - y * c * h_formal
    g_y = sp.simplify(d_y(g))
    g_yy = sp.simplify(d_y(g_y))
    expected_g_y = -c * h_formal + y * c**2 * sp.diff(h_formal, q)
    expected_g_yy = 2 * c**2 * sp.diff(h_formal, q) - y * c**3 * sp.diff(h_formal, q, 2)
    require(sp.simplify(g_y - expected_g_y) == 0, "first normal derivative drift")
    require(sp.simplify(g_yy - expected_g_yy) == 0, "second normal derivative drift")

    epsilon = sp.symbols("epsilon", nonnegative=True, real=True)
    g_function = sp.Function("G")(y)
    gaussian = sp.exp(-sp.pi * epsilon * y**2)
    weighted_second = sp.diff(gaussian * g_function, y, 2)
    expected_weighted_second = gaussian * (
        sp.diff(g_function, y, 2)
        - 4 * sp.pi * epsilon * y * sp.diff(g_function, y)
        + (4 * sp.pi**2 * epsilon**2 * y**2 - 2 * sp.pi * epsilon) * g_function
    )
    require(sp.simplify(weighted_second - expected_weighted_second) == 0, "Gaussian weighted derivative drift")

    x = sp.symbols("x", positive=True, real=True)
    normal_size = (y - sp.Rational(A, 2) * x) * sp.sqrt(2 / x)
    require(sp.simplify(sp.diff(normal_size, y) - sp.sqrt(2 / x)) == 0, "normal y monotonicity drift")
    require(sp.simplify(sp.diff(normal_size, x) + sp.sqrt(2) * y / (2 * x ** sp.Rational(3, 2)) + A / (2 * sp.sqrt(2 * x))) == 0, "normal x monotonicity drift")

    phase_prime = sp.pi * A**2 / 4 - T / (2 * x * (1 - x))
    phase_second = sp.diff(phase_prime, x)
    expected_phase_second = T * (1 - 2 * x) / (2 * x**2 * (1 - x) ** 2)
    require(sp.simplify(phase_second - expected_phase_second) == 0, "tangential curvature drift")

    q_left_endpoint = -normal_size.subs({y: left, x: sp.Rational(1, 2)})
    q_right_endpoint = -normal_size.subs({y: right, x: sp.Rational(1, 2)})
    require(q_left_endpoint == -sp.Rational(1, 2), "left endpoint q drift")
    require(q_right_endpoint == -sp.Rational(169, 2), "right endpoint q drift")

    phase_completion = sp.expand((y - sp.Rational(A, 2) * x) ** 2 / x - y**2 / x + A * y)
    require(sp.simplify(phase_completion - sp.Rational(A**2, 4) * x) == 0, "continuous phase cancellation drift")

    return {
        "strip": "[79789/2,79873/2]=[39894.5,39936.5]",
        "midpoint_roster": "m=39895,...,39936 are exactly the 42 unit-cell midpoints",
        "weighted_amplitude": "f_(epsilon,x)(y)=exp(-pi*epsilon*y^2)G_A(y,x)",
        "phase_stripped_tail": "G_A(y,x)=-1/(i*pi)-y*c(x)H(q_A(y,x)), c(x)=sqrt(2/x)",
        "Fresnel_ratio": "H(q)=exp(-i*pi*q^2/2)J_-(q)",
        "Fresnel_ODE": "H_q=1-i*pi*q*H; H_qq=-i*pi*q-(i*pi+pi^2*q^2)H",
        "normal_derivatives": "G_y=-cH+y*c^2*H_q; G_yy=2*c^2*H_q-y*c^3*H_qq",
        "weighted_second_derivative": "f_yy=exp(-pi*epsilon*y^2)[G_yy-4*pi*epsilon*y*G_y+(4*pi^2*epsilon^2*y^2-2*pi*epsilon)G]",
        "Peano_kernel": "K_m(y)=(1/2)(1/2-|y-m|)^2 on |y-m|<=1/2",
        "one_cell_kernel_mass": "1/24",
        "all_cell_kernel_mass": "7/4",
        "oriented_defect": "[-sum_m f(m)]-[-integral_strip f(y)dy]=sum_m integral_cell K_m(y)f_yy(y)dy",
        "zero_regulator_defect": "E_42(x)=sum_m integral_cell K_m(y)G_yy(y,x)dy",
        "continuous_phase_cancellation": "exp(i*pi*A*y)exp(-i*pi*y^2/x)exp(i*pi*q_A^2/2)=exp(i*pi*A^2*x/4)",
        "split_coordinate": "x_16 is the unique solution of (39894.5-A*x/2)sqrt(2/x)=16",
        "tail_region": "0<x<=x_16 implies -q_A(y,x)>=16 throughout the strip",
        "endpoint_region": "x_16<=x<=1/2 lies strictly after the unique half-domain tangential saddle",
    }


def numerical_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    a = arb(A)
    t = arb(T)
    left = arb(LEFT_NUMERATOR) / 2
    right = arb(RIGHT_NUMERATOR) / 2
    q_cut = arb(Q_CUT)

    sqrt_x_cut = ((2 * q_cut**2 + 8 * a * left).sqrt() - q_cut * arb(2).sqrt()) / (2 * a)
    x_cut = sqrt_x_cut**2
    x_star = (1 - (1 - 8 * t / (pi * a**2)).sqrt()) / 2
    phase_prime_cut = pi * a**2 / 4 - t / (2 * x_cut * (1 - x_cut))
    saddle_left_normal_size = (left - a * x_star / 2) * (2 / x_star).sqrt()
    endpoint_layer_width = arb(1) / 2 - x_cut
    endpoint_left_normal_size = (left - a / 4) * 2
    endpoint_right_normal_size = (right - a / 4) * 2

    require(x_cut.lower() > x_star.upper(), "x_16 does not lie after the saddle")
    require(phase_prime_cut.lower() > arb(21_000), "endpoint-layer phase floor drift")
    require(saddle_left_normal_size.lower() > arb(84), "saddle normal tail floor drift")
    require(endpoint_left_normal_size.contains(arb(1) / 2), "left endpoint transition drift")
    require(endpoint_right_normal_size.contains(arb(169) / 2), "right endpoint normal drift")
    require(endpoint_layer_width.upper() < arb("0.0001"), "endpoint layer width drift")

    return {
        "x_16_ball": x_cut.str(80, more=True),
        "x_star_ball": x_star.str(80, more=True),
        "x_16_minus_x_star_ball": (x_cut - x_star).str(80, more=True),
        "endpoint_layer_width_ball": endpoint_layer_width.str(80, more=True),
        "phase_prime_at_x_16_ball": phase_prime_cut.str(80, more=True),
        "phase_prime_floor_on_endpoint_layer": ">21000",
        "minimum_abs_q_on_strip_at_saddle_ball": saddle_left_normal_size.str(80, more=True),
        "endpoint_left_abs_q_ball": endpoint_left_normal_size.str(80, more=True),
        "endpoint_right_abs_q_ball": endpoint_right_normal_size.str(80, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    numerical = artifact["numerical_certificate"]
    telemetry = artifact["floating_route_telemetry"]
    return f"""# A-face midpoint-Peano reduction and saddle/endpoint partition

Date: 2026-08-23

Status: exact certificate for the fixed-regulator signed defect, Peano
reduction, and analytic tail/endpoint partition; quantitative derivative and
x-integral bounds remain open

After stripping the exact common A-face phase, put

```text
c(x)=sqrt(2/x),
q_A(y,x)=-(y-A*x/2)c(x),
H(q)=exp(-i*pi*q^2/2)J_-(q),
G_A(y,x)=-1/(i*pi)-y*c(x)H(q_A(y,x)).                (MP1)
```

For every positive Abel regulator define

```text
f_(epsilon,x)(y)=exp(-pi*epsilon*y^2)G_A(y,x).       (MP2)
```

The integers `39895,...,39936` are exactly the midpoints of the 42 unit
cells partitioning `[39894.5,39936.5]`.  Compact-strip Fubini therefore
reassembles the oriented primal-block/minus-translated-strip A-face channel as

```text
E_(42,epsilon)(x)
 =[-sum_(m=39895)^39936 f_(epsilon,x)(m)]
   -[-integral_(39894.5)^39936.5 f_(epsilon,x)(y)dy]. (MP3)
```

For `|y-m|<=1/2`, let

```text
K_m(y)=(1/2)(1/2-|y-m|)^2.                           (MP4)
```

Twice integrating the cell error, with the value and slope terms cancelling
by midpoint symmetry, gives the exact Peano identity

```text
E_(42,epsilon)(x)
 =sum_(m=39895)^39936
   integral_(m-1/2)^(m+1/2)K_m(y)
   partial_y^2 f_(epsilon,x)(y)dy.                   (MP5)
```

The kernel is positive and has the exact masses

```text
integral_cell K_m(y)dy=1/24,
sum_(42 cells) integral_cell K_m(y)dy=7/4.           (MP6)
```

Since the strip is finite and `G_A` is smooth there for `0<x<=1/2`, the
Abel-zero limit may be taken inside (MP5):

```text
E_42(x)=sum_m integral_cell K_m(y)G_(A,yy)(y,x)dy,
|E_42(x)|<=(1/24)sum_m sup_cell |G_(A,yy)(.,x)|.     (MP7)
```

This is the signed cancellation object.  Bounding the 42 modes or the strip
current separately discards (MP5).

The phase-stripped Fresnel ratio obeys

```text
H_q=1-i*pi*qH,
H_qq=-i*pi*q-(i*pi+pi^2*q^2)H,                       (MP8)

G_(A,y)=-cH+y*c^2 H_q,
G_(A,yy)=2*c^2 H_q-y*c^3 H_qq.                       (MP9)
```

At positive regulator the exact differentiated amplitude is

```text
partial_y^2 f_(epsilon,x)
 =exp(-pi*epsilon*y^2)
  [G_yy-4*pi*epsilon*y*G_y
   +(4*pi^2*epsilon^2*y^2-2*pi*epsilon)G].           (MP10)
```

No phase was fitted in this reduction.  Direct quadratic completion gives,
for real `y`,

```text
exp(i*pi*A*y)exp(-i*pi*y^2/x)exp(i*pi*q_A^2/2)
 =exp(i*pi*A^2*x/4).                                 (MP11)
```

The endpoint growth seen in floating telemetry requires a two-region
argument, not a false uniform-smallness claim.  Fix the rational normal
threshold `Q=16` and let

```text
x_16=
 [(sqrt(2*16^2+8*A*39894.5)-16sqrt(2))/(2A)]^2
 ={numerical['x_16_ball']}.                           (MP12)
```

Monotonicity of `(y-A*x/2)sqrt(2/x)` in both variables gives

```text
0<x<=x_16, 39894.5<=y<=39936.5
    ==> -q_A(y,x)>=16.                               (MP13)
```

The unique tangential saddle on the half-domain is

```text
x_*={numerical['x_star_ball']},
```

and `x_*<x_16`.  On the remaining endpoint layer, whose width is
`{numerical['endpoint_layer_width_ball']}`, the tangential phase is uniformly
nonstationary:

```text
Phi_A'(x)>=Phi_A'(x_16)
 ={numerical['phase_prime_at_x_16_ball']} >21000.     (MP14)
```

Thus the next rigorous estimate has a forced architecture:

1. On `0<x<=x_16`, use the uniform `|q_A|>=16` normal-tail hierarchy and
   retain the Morse treatment through `x_*`.
2. On `x_16<=x<=1/2`, keep the exact Peano amplitude and use tangential
   integration by parts, where (MP14) is valid.

The route pilot is intentionally nonrigorous.  It found
`|sum G_A(m)-integral G_A|` about
`{telemetry['saddle_midpoint_defect_modulus']:.12g}` at the saddle but
`{telemetry['endpoint_midpoint_defect_modulus']:.12g}` at `x=1/2`.  These
numbers select (MP12)--(MP14); neither is used in the proof.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (MP1)--(MP14) is inherited from the Gaussian
Abel regulator, Fourier character, canonical Fresnel primitive, and exact
Kummer phase.  The rational threshold `16` is an explicit analytic partition
choice, not a fitted occurrence of `pi` or a numerical theorem constant.

Proof boundary: exact fixed-regulator A-face block/strip orientation, Peano
identity, Abel-zero compact-strip passage, normal derivative ODE, continuous
phase cancellation, and tail/endpoint phase partition only.  No certified
`G_yy` envelope, endpoint-layer integration-by-parts remainder, tangential
Morse integral, quantitative `R_after_A`, `R_Dir`, or `Q_K-T` estimate,
all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "an exact dependency gate is not passed")
    require(dependencies["finite_regulator_equivalence"]["decision"]["finite_R_after_A_defined_on_one_common_regulator"] is True, "finite-regulator ownership drift")
    require(dependencies["normal_tangential_split"]["decision"]["edge_translation_equals_explicit_42_mode_block"] is True, "42-mode translation drift")
    require(dependencies["normal_tangential_split"]["decision"]["normal_then_tangential_analysis_licensed"] is True, "normal/tangential route drift")
    require(dependencies["two_current_expansion"]["decision"]["continuous_edge_A_tail_defined_by_exact_finite_endpoint_current"] is True, "continuous A-tail current drift")

    pilot = load_json(PILOT)
    require(pilot.get("passed") is True, "route pilot is unavailable")
    require(pilot["decision"]["floating_values_used_as_proof"] is False, "pilot proof-boundary drift")
    saddle = next(row for row in pilot["rows"] if row["is_tangential_saddle"])
    endpoint = next(row for row in pilot["rows"] if row["x"] == 0.5)

    artifact = {
        "kind": STEM,
        "status": "exact_A_face_fixed_regulator_midpoint_Peano_reduction_and_tail_endpoint_partition_certified_quantitative_bounds_open",
        "passed": True,
        "scope": {
            "height": T,
            "A": A,
            "shifted_modes": [FIRST_MODE, LAST_MODE],
            "continuous_strip": ["79789/2", "79873/2"],
            "Abel_regulator": "epsilon>0 followed by epsilon down to 0 on the finite compact strip",
            "normal_tail_threshold": Q_CUT,
        },
        "symbolic_certificate": symbolic_certificate(),
        "numerical_certificate": numerical_certificate(),
        "floating_route_telemetry": {
            "path": relative(PILOT),
            "sha256": file_hash(PILOT),
            "used_as_proof": False,
            "saddle_midpoint_defect_modulus": saddle["midpoint_defect"]["modulus"],
            "endpoint_midpoint_defect_modulus": endpoint["midpoint_defect"]["modulus"],
            "interpretation": "The endpoint sample rejects a uniform-smallness route and selects the exact tail/endpoint partition only.",
        },
        "decision": {
            "fixed_regulator_42_mode_and_translated_strip_channel_oriented_exactly": True,
            "midpoint_Peano_identity_proved": True,
            "zero_regulator_compact_strip_limit_proved": True,
            "phase_stripped_normal_derivative_ODE_proved": True,
            "continuous_y_phase_cancellation_proved": True,
            "uniform_abs_q_at_least_16_before_x_16_proved": True,
            "endpoint_layer_strictly_after_tangential_saddle_proved": True,
            "endpoint_layer_phase_derivative_above_21000_proved": True,
            "uniform_quantitative_G_yy_envelope_proved": False,
            "endpoint_layer_IBP_remainder_bounded": False,
            "tangential_Morse_integral_bounded": False,
            "signed_42_mode_cancellation_bound_proved": False,
            "R_after_A_bound_proved": False,
            "R_Dir_bound_proved": False,
            "QK_minus_T_bound_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "arb_threads": 1,
            "process_priority": priority,
            "precision_decimal_digits": PRECISION,
        },
        "next_obligation": "On 0<x<=x_16, differentiate a sufficiently deep exact Fresnel-tail remainder to certify the cellwise G_A,yy Peano envelope while retaining the tangential Morse core. On x_16<=x<=1/2, derive a complex interval envelope for E_42 and its x derivative and close one tangential integration-by-parts remainder using Phi_A'>21000. Only then reinsert the remaining B, zero, negative, half-current, and remote-positive channels.",
        "proof_boundary": "Exact fixed-regulator A-face block/strip orientation, Peano identity, Abel-zero compact-strip passage, normal derivative ODE, continuous phase cancellation, and tail/endpoint phase partition only. No certified G_yy envelope, endpoint-layer integration-by-parts remainder, tangential Morse integral, quantitative R_after_A, R_Dir, or Q_K-T estimate, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("certified A-face midpoint-Peano reduction and tail/endpoint partition", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
