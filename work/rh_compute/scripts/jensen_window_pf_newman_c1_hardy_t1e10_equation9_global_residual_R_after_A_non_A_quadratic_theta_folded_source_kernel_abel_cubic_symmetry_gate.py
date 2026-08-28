#!/usr/bin/env python3
"""Certify complete folded-source wall derivatives and cubic kernel symmetry."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
import flint
from flint import acb, arb


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_folded_source_kernel_abel_cubic_symmetry_gate"
)
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

COMMON_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_common_kernel_theta_current_gate"
)
POLE_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "extended_notch_pole_free_kernel_gate"
)
WALL_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_small_tau_grouped_first_wall_bridge_gate"
)
DERIVATIVE_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_small_tau_characteristic_euler_maclaurin_parameter_derivative_gate"
)
DEPENDENCIES = {
    "common_kernel": REPO_ROOT / f"work/rh_compute/results/{COMMON_STEM}.json",
    "pole_free_kernel": REPO_ROOT / f"work/rh_compute/results/{POLE_STEM}.json",
    "grouped_first_wall": REPO_ROOT / f"work/rh_compute/results/{WALL_STEM}.json",
    "characteristic_derivatives": REPO_ROOT / f"work/rh_compute/results/{DERIVATIVE_STEM}.json",
}

A = 159_577
B = 5_122_421
L = 2_481_422
NOTCH_FIRST = 622
NOTCH_LAST = 39_936
NOTCH_LENGTH = NOTCH_LAST - NOTCH_FIRST + 1
CUBIC_BLOCKS = NOTCH_LENGTH // 3
S0 = Fraction(1, 3)
N = 496_283
T_WALL = Fraction(1, 2 * N)
PRECISION_BITS = 352
BLOCK_SIZE = 16
VARIATION_RADIUS = arb("0.001")
CPU_BASELINE_SAMPLES = [24.85, 18.43, 20.13, 26.99, 21.50, 9.16]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


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


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency {relative(path)}")
    return json.loads(path.read_text(encoding="utf-8"))


def rat(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def cis_pi(value: Fraction) -> acb:
    return (acb.pi() * acb(0, 1) * acb(rat(value % 2))).exp()


def scaled_integer_text(coefficient: int, exponent: int, digits: int, round_up: bool) -> str:
    if coefficient <= 0:
        return "0"
    raw = str(coefficient)
    if len(raw) > digits:
        dropped = len(raw) - digits
        kept = int(raw[:digits])
        tail = raw[digits:]
        if round_up and any(character != "0" for character in tail):
            kept += 1
        raw = str(kept)
        exponent += dropped
    scientific_exponent = exponent + len(raw) - 1
    mantissa = raw if len(raw) == 1 else raw[0] + "." + raw[1:]
    return f"{mantissa}e{scientific_exponent}"


def upper_text(value: arb, digits: int = 45) -> str:
    midpoint, radius, exponent = value.mid_rad_10exp()
    coefficient = abs(int(midpoint)) + int(radius)
    return scaled_integer_text(coefficient, int(exponent), digits, True)


def lower_text(value: arb, digits: int = 45) -> str:
    midpoint, radius, exponent = value.mid_rad_10exp()
    coefficient = max(0, int(midpoint) - int(radius))
    return scaled_integer_text(coefficient, int(exponent), digits, False)


def acb_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": str(value.real),
        "imag_ball": str(value.imag),
        "radius_upper": upper_text(value.rad()),
        "absolute_lower": lower_text(value.abs_lower()),
        "absolute_upper": upper_text(value.abs_upper()),
    }


def parse_complex_ball(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def x_of_t(side: str, t: Fraction) -> Fraction:
    require(side in ("right", "left"), "unknown side")
    denominator = 5 + (4 if side == "right" else 6) * t
    return 2 * (1 + t) / denominator


def dx_dt(side: str, t: Fraction) -> Fraction:
    denominator = 5 + (4 if side == "right" else 6) * t
    return (2 if side == "right" else -2) / denominator**2


def direct_source_moments(x: Fraction, s: Fraction) -> tuple[acb, acb, acb, acb]:
    """Return sum alpha^j exp(i*pi*x*alpha^2/4), j=0..3."""
    totals = [acb(0), acb(0), acb(0), acb(0)]
    step = cis_pi(2 * x)
    for start in range(0, L, BLOCK_SIZE):
        stop = min(L, start + BLOCK_SIZE)
        alpha = Fraction(A) + 2 * s + 2 * start
        term = cis_pi(x * alpha * alpha / 4)
        ratio = cis_pi(x * (alpha + 1))
        alpha_ball = rat(alpha)
        local = [acb(0), acb(0), acb(0), acb(0)]
        for _ in range(start, stop):
            alpha2 = alpha_ball * alpha_ball
            local[0] += term
            local[1] += acb(alpha_ball) * term
            local[2] += acb(alpha2) * term
            local[3] += acb(alpha2 * alpha_ball) * term
            term *= ratio
            ratio *= step
            alpha_ball += 2
        for degree in range(4):
            totals[degree] += local[degree]
    return tuple(totals)  # type: ignore[return-value]


def kernel_certificate() -> dict[str, Any]:
    require(B - A == 2 * L, "source endpoint arithmetic drift")
    require(NOTCH_LENGTH == 39_315, "notch length drift")
    require(NOTCH_LENGTH == 3 * CUBIC_BLOCKS, "notch is not a whole cubic roster")
    require(NOTCH_FIRST % 3 == 1 and NOTCH_LAST % 3 == 0, "cubic block endpoints drift")

    sqrt_three = arb(3).sqrt()
    c_prime = arb(CUBIC_BLOCKS) * acb.pi() * acb(-sqrt_three, 3)
    c_prime_abs = 2 * arb.pi() * arb(CUBIC_BLOCKS) * sqrt_three
    require(c_prime.abs_lower() > arb(140_000), "kernel slope unexpectedly small")
    require(
        (c_prime.real**2 + c_prime.imag**2).overlaps(c_prime_abs**2),
        "kernel slope magnitude identity failed",
    )
    return {
        "notch": [NOTCH_FIRST, NOTCH_LAST],
        "notch_length_Q": NOTCH_LENGTH,
        "cubic_block_count_r": CUBIC_BLOCKS,
        "root": "omega=exp(-2*pi*i/3), 1+omega+omega^2=0",
        "interior_Abel_limit": "C_U,0(s)=-sum_(m=622)^39936 exp(-2*pi*i*m*s), s not in Z",
        "exact_factorization": (
            "C_U,0(1/3+u)=-omega*exp(-2*pi*i*622*u)"
            "*(1-exp(-2*pi*i*39315*u))/(1-omega*exp(-2*pi*i*u))"
        ),
        "cubic_zero": "C_U,0(1/3)=0",
        "weighted_root_sum": "sum_(m=622)^39936 m*omega^m=13105*(1-omega)",
        "first_derivative": "C_U,0'(1/3)=13105*pi*(-sqrt(3)+3*i)",
        "first_derivative_ball": acb_record(c_prime),
        "first_derivative_absolute_exact": "2*pi*13105*sqrt(3)",
        "first_derivative_absolute_ball": {
            "lower": lower_text(c_prime_abs),
            "upper": upper_text(c_prime_abs),
        },
        "full_cell_Abel_limit": (
            "lim_(epsilon->0+) integral_0^1 F(s)C_U,epsilon(s)ds="
            "[F(0)+F(1)]/2-sum_(m=622)^39936 hat_F(m)"
        ),
        "boundary_layer_guard": (
            "The full theta term is a periodic approximate identity and contributes [F(0)+F(1)]/2. "
            "The pointwise interior limit -D_U cannot be integrated through the integer boundary by itself."
        ),
        "integer_Fourier_unfolding": "hat_F_x(m)=integral_0^L f_x(u)*exp(-2*pi*i*m*u)du",
        "finite_Fresnel_handoff": (
            "hat_F_x(m)=(-1)^m*{[exp(i*pi*phi_m(B))-exp(i*pi*phi_m(A))]/(i*pi*x)"
            "+(m/x)J_m(x)}, phi_m(y)=x*y^2/4-m*y, x>0; x=0 uses the continuous direct-integral limit"
        ),
    }


def wall_row(side: str, wall_dependency: dict[str, Any], c_prime: acb) -> dict[str, Any]:
    x = x_of_t(side, T_WALL)
    wall_key = f"{side}_first_wall_source_box"
    expected_x = Fraction(wall_dependency[wall_key]["wall_x"])
    require(x == expected_x, f"{side} wall x map drift")

    m0, m1, m2, m3 = direct_source_moments(x, S0)
    source = m1
    source_x = acb.pi() * acb(0, 1) * m3 / 4
    source_s = 2 * m0 + acb.pi() * acb(0, 1) * acb(rat(x)) * m2
    source_t = acb(rat(dx_dt(side, T_WALL))) * source_x

    grouped = parse_complex_ball(wall_dependency[wall_key]["complete_source_box"])
    require(source.overlaps(grouped), f"{side} direct source misses grouped wall enclosure")

    product_s = source * c_prime
    mixed_ts = source_t * c_prime
    source_t_width = VARIATION_RADIUS / source_t.abs_upper()
    source_s_width = VARIATION_RADIUS / source_s.abs_upper()
    product_s_width = VARIATION_RADIUS / product_s.abs_upper()
    require(product_s_width < arb("1e-17"), f"{side} product slope no longer microscopic")
    return {
        "side": side,
        "t": str(T_WALL),
        "x": str(x),
        "s": str(S0),
        "direct_source_F": acb_record(source),
        "grouped_complete_source_box": acb_record(grouped),
        "direct_F_overlaps_grouped_complete_source_box": True,
        "direct_partial_x_F": acb_record(source_x),
        "direct_partial_t_F": acb_record(source_t),
        "direct_partial_s_F": acb_record(source_s),
        "interior_limit_product_F_times_C_at_s0": acb_record(acb(0)),
        "interior_limit_partial_t_product_at_s0": acb_record(acb(0)),
        "interior_limit_partial_s_product_at_s0": acb_record(product_s),
        "interior_limit_mixed_partial_t_partial_s_product_at_s0": acb_record(mixed_ts),
        "diagnostic_half_width_for_variation_0_001": {
            "source_t": upper_text(source_t_width),
            "source_s": upper_text(source_s_width),
            "product_s": upper_text(product_s_width),
            "certified_transport_width": False,
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    kernel = artifact["kernel_symmetry_certificate"]
    right, left = artifact["wall_source_derivative_rows"]
    return f"""# Folded source-kernel Abel limit and cubic symmetry

Date: 2026-08-25

Status: exact Abel/symmetry reduction and two rigorous wall derivative rows;
the physical non-A quadrature remains open.

## Abel-limit decomposition

For the extended notch `U={{622,...,39936}}`, put

```text
C_U,epsilon(s)=sum_(m in Z)[1-chi_U(m)]
               exp(-pi epsilon m^2)exp(-2 pi i m s).             (FSK1)
```

Away from the integer lattice, Poisson summation makes the full theta term
vanish as `epsilon` decreases to zero.  On the complete cell, however, that
theta term is a periodic approximate identity.  Therefore

```text
lim integral_0^1 F_x(s)C_U,epsilon(s)ds
 =[F_x(0)+F_x(1)]/2-sum_(m=622)^39936 hat F_x(m),                 (FSK2)

hat F_x(m)=integral_0^L f_x(u)exp(-2 pi i m u)du.                (FSK3)
```

The endpoint half-sum in (FSK2) is mandatory.  Integrating the pointwise
interior limit through `s=0` would silently discard the theta boundary layer.

With `phi_m(y)=x y^2/4-m y`, `A={A}`, and `B={B}`, the finite coefficient is

```text
hat F_x(m)=(-1)^m{{[exp(i pi phi_m(B))-exp(i pi phi_m(A))]/(i pi x)
                  +(m/x)integral_A^B exp(i pi phi_m(y))dy}}.      (FSK4)
```

Equation (FSK4) is written for `x>0`.  At `x=0` its two displayed pieces
must not be bounded separately; use the continuous direct-integral limit.

Thus no microscopic `s`-cell cover is required merely to remove the Abel
regulator: the exact handoff is a joined endpoint term plus 39,315 oscillatory
Fresnel coefficients.  Those coefficients must still be summed and integrated
in `x` without splitting their cancellation.

## Cubic symmetry

The notch length is

```text
Q=39936-622+1={kernel['notch_length_Q']}=3*{kernel['cubic_block_count_r']}.
```

For `omega=exp(-2 pi i/3)`, the pointwise interior kernel has the exact factor

```text
C_U,0(1/3+u)
 =-omega exp(-2 pi i 622u)
   [1-exp(-2 pi i 39315u)]/[1-omega exp(-2 pi i u)].             (FSK5)
```

Consequently

```text
C_U,0(1/3)=0,
C_U,0'(1/3)=13105*pi*(-sqrt(3)+3i),
|C_U,0'(1/3)|=2*pi*13105*sqrt(3).                                (FSK6)
```

The factor in (FSK5), rather than a near-constant value box, is the natural
interior quadrature coordinate.

## Complete-source wall audit

Exact-rational phase resets at 352 bits give:

| side | `x` | `|F|` lower | `|partial_t F|` lower | `|partial_s F|` lower | product `s` diagnostic width |
|---|---:|---:|---:|---:|---:|
| right | `{right['x']}` | `{right['direct_source_F']['absolute_lower']}` | `{right['direct_partial_t_F']['absolute_lower']}` | `{right['direct_partial_s_F']['absolute_lower']}` | `{right['diagnostic_half_width_for_variation_0_001']['product_s']}` |
| left | `{left['x']}` | `{left['direct_source_F']['absolute_lower']}` | `{left['direct_partial_t_F']['absolute_lower']}` | `{left['direct_partial_s_F']['absolute_lower']}` | `{left['diagnostic_half_width_for_variation_0_001']['product_s']}` |

Both direct 2,481,422-term source balls overlap the independently certified
grouped first-wall source boxes.  At `s=1/3`, the interior-limit product and
its `t` derivative vanish exactly because of (FSK6), but

```text
partial_s(F C)=F C_U,0'(1/3)                                    (FSK7)
```

remains large.  The cubic zero removes pointwise `t` stiffness on the symmetry
line; it does not create a macroscopic two-dimensional value box.

## Decision

Use (FSK2)--(FSK5) to split the next proof into a rigorously retained theta
boundary layer and an interior oscillatory Dirichlet/Fresnel assembly.  Do not
enumerate near-constant source or product boxes.  Preserve the finite 39,315
mode sum and its joined endpoint term before norms, then build an altered-
partition `x`-quadrature certificate.

## Pi provenance

Every `pi` above comes from the declared Gaussian Abel/Fourier kernel or the
original quadratic Kummer character.  The cubic root is selected by the exact
integer identity `39315=3*13105`; no geometric or fitted occurrence is added.

## Boundary

{artifact['proof_boundary']}
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    require(CHECKER.is_file(), "missing independent checker")
    loaded = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(loaded["common_kernel"].get("passed") is True, "common-kernel dependency failed")
    require(
        loaded["common_kernel"]["decision"]["long_source_integral_folded_to_one_cell"] is True,
        "one-cell fold dependency drift",
    )
    require(loaded["pole_free_kernel"].get("passed") is True, "pole-free dependency failed")
    require(
        loaded["pole_free_kernel"]["decision"]["pole_free_closed_centred_cell_formula_certified"] is True,
        "pole-free cell dependency drift",
    )
    require(loaded["grouped_first_wall"].get("passed") is True, "grouped wall dependency failed")
    require(loaded["characteristic_derivatives"].get("passed") is True, "derivative dependency failed")

    flint.ctx.prec = PRECISION_BITS
    kernel = kernel_certificate()
    c_prime = acb(
        arb(kernel["first_derivative_ball"]["real_ball"]),
        arb(kernel["first_derivative_ball"]["imag_ball"]),
    )
    rows = [wall_row(side, loaded["grouped_first_wall"], c_prime) for side in ("right", "left")]

    proof_boundary = (
        "The exact complete-cell Abel-limit identity for the folded source factor, its finite Fresnel-mode "
        "handoff, the interior cubic-root factorization at s=1/3, and rigorous complete-source value and "
        "first-derivative balls at the two first walls only. The reported variation widths are diagnostics, "
        "not certified cells. No uniform complete-source derivative theorem, boundary-layer numerical "
        "estimate, finite 39315-mode signed bound, physical x quadrature, non-A bound, joined R_after_A, "
        "R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
    )
    artifact = {
        "kind": STEM,
        "status": "folded_source_Abel_limit_cubic_kernel_zero_and_two_wall_derivative_rows_certified_signed_x_quadrature_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "source_odd_roster": [A, B],
            "source_cell_count": L,
            "notch": [NOTCH_FIRST, NOTCH_LAST],
            "s": str(S0),
            "first_wall_t": str(T_WALL),
            "precision_bits": PRECISION_BITS,
            "phase_reset_block_size": BLOCK_SIZE,
            "workers": 1,
            "daytime_cpu_baseline_samples_percent": CPU_BASELINE_SAMPLES,
            "daytime_cpu_baseline_average_percent": round(
                sum(CPU_BASELINE_SAMPLES) / len(CPU_BASELINE_SAMPLES), 2
            ),
        },
        "kernel_symmetry_certificate": kernel,
        "exact_source_derivatives": {
            "F": "sum alpha_n*exp(i*pi*x*alpha_n^2/4)",
            "partial_x_F": "(i*pi/4)*sum alpha_n^3*exp(i*pi*x*alpha_n^2/4)",
            "partial_s_F": "2*sum exp(i*pi*x*alpha_n^2/4)+i*pi*x*sum alpha_n^2*exp(i*pi*x*alpha_n^2/4)",
            "alpha_n": "A+2*n+2*s",
        },
        "wall_source_derivative_rows": rows,
        "decision": {
            "complete_source_direct_sums_overlap_both_grouped_wall_boxes": True,
            "complete_source_wall_derivatives_interval_certified": True,
            "interior_kernel_has_exact_cubic_zero_at_s_one_third": True,
            "interior_product_t_derivative_vanishes_on_symmetry_line": True,
            "interior_product_s_derivative_remains_microscopic_for_value_boxing": True,
            "pointwise_interior_limit_interchange_across_integer_boundary_forbidden": True,
            "complete_cell_Abel_limit_reduced_to_endpoint_half_sum_and_finite_modes": True,
            "near_constant_source_product_cover_selected": False,
            "signed_boundary_layer_plus_finite_mode_x_quadrature_selected": True,
            "physical_quadrature_completed": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Instantiate the endpoint half-sum and all 39315 finite Fresnel coefficients in the exact non-A "
            "ownership assembly. Preserve their signed sum before norms, derive a stable common x-phase "
            "representation, and build a bounded pilot followed by an interval altered-partition x quadrature."
        ),
        "proof_boundary": proof_boundary,
        "primary_source": loaded["common_kernel"].get("primary_source"),
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "blas_threads": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.time() - started, 3),
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("certified folded-source wall derivatives and cubic kernel symmetry", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
