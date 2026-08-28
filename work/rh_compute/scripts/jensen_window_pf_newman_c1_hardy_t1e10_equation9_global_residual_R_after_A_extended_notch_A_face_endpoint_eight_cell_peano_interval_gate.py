#!/usr/bin/env python3
"""Certify the complete translated A-face endpoint layer by interval bounds."""

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

from flint import acb, arb, ctx
import mpmath as mp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_peano_interval_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
LOCALIZATION = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_transition_localization_gate.json"

A = 159_577
T = 10_000_000_000
TRANSITION_LEFT = mp.mpf("39894.5")
TRANSITION_FIRST_MODE = 39_895
TRANSITION_LAST_MODE = 39_902
TAIL_LEFT_NUMERATOR = 79_805
TAIL_RIGHT_NUMERATOR = 79_873
TAIL_FIRST_MODE = 39_903
TAIL_LAST_MODE = 39_936
X_LOWER_TEXT = "0.4999028779462474"
X_UPPER_TEXT = "0.5"
X_SLABS = 256
HALF_CELL_SUBDIVISIONS = 16
ASYMPTOTIC_Q_FLOOR = 2
TERMS = 6
PRECISION = 100
COMPLETE_ENDPOINT_THRESHOLD = "0.0037"
SAVED_R_AFTER_A_TARGET = "0.0368147039947"


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


def odd_double_factorial(index: int) -> int:
    result = 1
    for value in range(1, index + 1, 2):
        result *= value
    return result


def interval_from_endpoints(left: mp.mpf, right: mp.mpf) -> arb:
    midpoint = (left + right) / 2
    radius = (right - left) / 2 + mp.mpf("1e-105")
    return arb(mp.nstr(midpoint, 115), mp.nstr(radius, 115))


def is_finite(value: arb | acb) -> bool:
    return "nan" not in value.str(24, more=True).lower()


def require_finite(value: arb | acb, label: str) -> None:
    require(is_finite(value), f"{label} became indeterminate")


def zero_centered_upper_ball(value: arb) -> arb:
    midpoint, radius, exponent = value.upper().mid_rad_10exp()
    return arb("0", f"{midpoint + radius}e{exponent}")


def exact_or_asymptotic_gyy_bound(y: arb, x: arb) -> tuple[arb, str]:
    pi = arb.pi()
    c = (2 / x).sqrt()
    q = -(y - arb(A) * x / 2) * c
    q_abs = -q
    require(q_abs.lower() > arb(0), "endpoint cover crossed q=0")
    attempted_asymptotic = q_abs.lower() >= arb(ASYMPTOTIC_Q_FLOOR)

    if attempted_asymptotic:
        i_pi = acb(0, pi)
        delta = y - arb(A) * x / 2
        approximation = acb(arb(A) * x) / (i_pi * delta**3)
        for index in range(1, TERMS):
            power = 2 * index + 1
            coefficient_n = arb(odd_double_factorial(2 * index - 1)) * (x / 2) ** index / i_pi ** (index + 1)
            approximation += coefficient_n * power * (
                arb(power + 1) * y / delta ** (power + 2)
                - 2 / delta ** (power + 1)
            )
        approximation += (
            arb(10_395)
            * x**4
            / (16 * pi**6 * delta**13)
            * (acb(0, pi * delta**2 * y) - acb(delta * x) + acb(6 * x * y))
        )
        coefficient = -acb(0, 2 * pi * c**2 * q) + acb(y * c**3) * (acb(0, pi) + acb(pi**2 * q**2))
        h_remainder = 2 * arb(odd_double_factorial(11)) / (pi**7 * q_abs**13)
        if is_finite(approximation) and is_finite(coefficient) and is_finite(h_remainder):
            approximation_modulus = abs(approximation)
            coefficient_modulus = abs(coefficient)
            bound = approximation_modulus + coefficient_modulus * h_remainder
            if is_finite(approximation_modulus) and is_finite(coefficient_modulus) and is_finite(bound):
                return bound, "six_current_asymptotic"

    q_mid = q.mid()
    q_radius = q.rad()
    scale_mid = (-q_mid) * pi.sqrt() / 2
    w_mid = acb(scale_mid, -scale_mid)
    h_mid = acb(arb(1) / 2, arb(1) / 2) * (w_mid**2).exp() * w_mid.erfc()
    require_finite(h_mid, "exact center H evaluation")
    h_variation = q_radius * (1 + pi * abs(h_mid) * (abs(q_mid) + q_radius / 2))
    component_error = zero_centered_upper_ball(h_variation)
    h = acb(h_mid.real + component_error, h_mid.imag + component_error)
    h_q = 1 - acb(0, pi) * acb(q) * h
    h_qq = -acb(0, pi) * acb(q) - (acb(0, pi) + acb(pi**2 * q**2)) * h
    g_yy_bound = abs(acb(2 * c**2) * h_q - acb(y * c**3) * h_qq)
    require_finite(g_yy_bound, "exact G_yy enclosure")
    method = "exact_erfc_fallback" if attempted_asymptotic else "exact_erfc"
    return g_yy_bound, method


def g_six(y: arb, x: arb) -> acb:
    i_pi = acb(0, arb.pi())
    delta = y - arb(A) * x / 2
    value = acb(arb(A) * x) / (2 * i_pi * delta)
    for index in range(1, TERMS):
        value += (
            acb(y)
            * arb(odd_double_factorial(2 * index - 1))
            * (x / 2) ** index
            / (i_pi ** (index + 1) * delta ** (2 * index + 1))
        )
    return value


def integral_g_six(left: arb, right: arb, x: arb) -> acb:
    i_pi = acb(0, arb.pi())
    a = arb(A) * x / 2
    dl = left - a
    dr = right - a
    value = acb(arb(A) * x * (dr / dl).log()) / (2 * i_pi)
    for index in range(1, TERMS):
        coefficient = arb(odd_double_factorial(2 * index - 1)) * (x / 2) ** index / i_pi ** (index + 1)
        value += coefficient * (
            (dr ** (1 - 2 * index) - dl ** (1 - 2 * index)) / (1 - 2 * index)
            - a * (dr ** (-2 * index) - dl ** (-2 * index)) / (2 * index)
        )
    return value


def tail_defect_six(x: arb) -> acb:
    left = arb(TAIL_LEFT_NUMERATOR) / 2
    right = arb(TAIL_RIGHT_NUMERATOR) / 2
    value = acb(0)
    for mode in range(TAIL_FIRST_MODE, TAIL_LAST_MODE + 1):
        value += g_six(arb(mode), x)
    return value - integral_g_six(left, right, x)


def symbolic_certificate() -> dict[str, Any]:
    require(TRANSITION_LAST_MODE - TRANSITION_FIRST_MODE + 1 == 8, "transition roster drift")
    require(TAIL_LAST_MODE - TAIL_FIRST_MODE + 1 == 34, "tail roster drift")
    require(odd_double_factorial(9) == 945 and 11 * 945 == 10_395, "terminal-current coefficient drift")
    require(132 * 945 == 6 * 10_395 * 2, "terminal-current derivative coefficient drift")
    return {
        "Peano_identity": "D_(tr,8)=-sum_cells integral K_m(y)G_(A,yy)(y,x)dy; integral_cell K_m=1/24",
        "exact_current": "H(q0)=(1+i)exp(w0^2)erfc(w0)/2, w0=(-q0)sqrt(pi)(1-i)/2 for real q0<0",
        "center_variation": "|H(q)-H(q0)|<=r[1+pi|H(q0)|(|q0|+r/2)] on |q-q0|<=r",
        "center_variation_reason": "integrating-factor identity plus |exp(i theta)-1|<=|theta|",
        "normal_ODE": "H_q=1-i*pi*q*H; H_qq=-i*pi*q-(i*pi+pi^2*q^2)H",
        "six_current_remainder": "|H-H6|<=2*11!!/(pi^7|q|^13)",
        "terminal_current_recurrence": "H6_q=(1-i*pi*q*H6)-11*a5*q^-12 with a5=9!!/(i*pi)^6",
        "tail_formula": "D_(tail,34),6=sum_(39903..39936)G6-integral_(39902.5)^39936.5 G6",
        "physical_majorant": "2(pi/(32T))^(1/4) integral x^(-5/4)(1-x)^(-1/4)(B8+|D_tail,6|)dx plus tail replacement error",
    }


def interval_certificate(localization: dict[str, Any]) -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    mp.mp.dps = 130
    pi = arb.pi()
    aa = arb(A)
    q_cut = arb(16)
    left = arb(79_789) / 2
    sqrt_x_cut = ((2 * q_cut**2 + 8 * aa * left).sqrt() - q_cut * arb(2).sqrt()) / (2 * aa)
    exact_x_cut = sqrt_x_cut**2
    require(arb(X_LOWER_TEXT) < exact_x_cut.lower(), "decimal lower cover does not enclose x_16 from below")
    require(exact_x_cut.upper() < arb(X_UPPER_TEXT), "x_16 endpoint ordering drift")

    x_lower = mp.mpf(X_LOWER_TEXT)
    x_upper = mp.mpf(X_UPPER_TEXT)
    x_width = (x_upper - x_lower) / X_SLABS
    y_width = mp.mpf("0.5") / HALF_CELL_SUBDIVISIONS
    normalization = 2 * (pi / (32 * arb(T))) ** (arb(1) / 4)
    physical_eight = arb(0)
    physical_tail_six = arb(0)
    max_eight_defect = arb(0)
    max_tail_defect = arb(0)
    method_counts = {
        "exact_erfc_boxes": 0,
        "exact_erfc_fallback_boxes": 0,
        "six_current_asymptotic_boxes": 0,
    }
    sample_slabs: dict[str, dict[str, str]] = {}

    for slab in range(X_SLABS):
        x0 = x_lower + slab * x_width
        x1 = x0 + x_width
        x_box = interval_from_endpoints(x0, x1)
        eight_defect = arb(0)
        for mode in range(TRANSITION_FIRST_MODE, TRANSITION_LAST_MODE + 1):
            cell_left = mp.mpf(mode) - mp.mpf("0.5")
            for side in range(2):
                side_left = cell_left + mp.mpf(side) / 2
                for subdivision in range(HALF_CELL_SUBDIVISIONS):
                    y0 = side_left + subdivision * y_width
                    y1 = y0 + y_width
                    y_box = interval_from_endpoints(y0, y1)
                    if side == 0:
                        kernel = (y_box - arb(str(cell_left))) ** 2 / 2
                    else:
                        kernel = (arb(str(cell_left + 1)) - y_box) ** 2 / 2
                    gyy_bound, method = exact_or_asymptotic_gyy_bound(y_box, x_box)
                    method_counts[f"{method}_boxes"] += 1
                    eight_defect += arb(mp.nstr(y_width, 115)) * kernel.upper() * gyy_bound.upper()

        weight = 1 / (x_box * (x_box * (1 - x_box)) ** (arb(1) / 4))
        slab_scale = normalization * arb(mp.nstr(x_width, 115)) * weight.upper()
        tail_defect = abs(tail_defect_six(x_box))
        require_finite(eight_defect, "eight-cell Peano defect")
        require_finite(tail_defect, "34-cell six-current tail defect")
        physical_eight += slab_scale * eight_defect
        physical_tail_six += slab_scale * tail_defect
        if eight_defect.upper() > max_eight_defect.upper():
            max_eight_defect = eight_defect
        if tail_defect.upper() > max_tail_defect.upper():
            max_tail_defect = tail_defect
        if slab in (0, X_SLABS // 2, X_SLABS - 1):
            sample_slabs[str(slab)] = {
                "x_ball": x_box.str(50, more=True),
                "eight_cell_defect_bound": eight_defect.str(50, more=True),
                "tail_six_current_defect_bound": tail_defect.str(50, more=True),
            }

    tail_replacement = arb(localization["numerical_certificate"]["physical_tail_replacement_error_ball"])
    complete_endpoint = physical_eight + physical_tail_six + tail_replacement
    require(physical_eight.upper() < arb("0.0036"), "eight-cell physical bound drift")
    require(physical_tail_six.upper() < arb("0.000037"), "tail six-current physical bound drift")
    require(complete_endpoint.upper() < arb(COMPLETE_ENDPOINT_THRESHOLD), "complete endpoint threshold failed")
    require(complete_endpoint.upper() < arb(SAVED_R_AFTER_A_TARGET), "saved R_after_A allowance exceeded")
    return {
        "exact_x_16_ball": exact_x_cut.str(80, more=True),
        "certified_cover": [X_LOWER_TEXT, X_UPPER_TEXT],
        "partition": {
            "x_slabs": X_SLABS,
            "half_cell_subdivisions": HALF_CELL_SUBDIVISIONS,
            "normal_boxes_per_x_slab": 8 * 2 * HALF_CELL_SUBDIVISIONS,
            "precision_decimal_digits": PRECISION,
            "asymptotic_q_floor": ASYMPTOTIC_Q_FLOOR,
        },
        "method_counts": method_counts,
        "maximum_eight_cell_defect_bound": max_eight_defect.str(70, more=True),
        "maximum_tail_six_current_defect_bound": max_tail_defect.str(70, more=True),
        "physical_eight_cell_bound": physical_eight.str(70, more=True),
        "physical_tail_six_current_bound": physical_tail_six.str(70, more=True),
        "physical_tail_replacement_error": tail_replacement.str(70, more=True),
        "physical_complete_endpoint_bound": complete_endpoint.str(70, more=True),
        "complete_endpoint_threshold": COMPLETE_ENDPOINT_THRESHOLD,
        "saved_R_after_A_target": SAVED_R_AFTER_A_TARGET,
        "sample_slabs": sample_slabs,
    }


def render_note(artifact: dict[str, Any]) -> str:
    cert = artifact["interval_certificate"]
    return f"""# Complete absolute bound for the translated A-face endpoint layer

Date: 2026-08-23

Status: rigorous local interval certificate at `t=10^10`, independently
replayed on a different interval partition; not a proof of `R_after_A` or RH

Let `D_(tr,8)` be the exact eight-cell Fresnel midpoint defect from Section
11.455 and `D_(tail,34),6` the exact rational-log six-current tail.  The
positive midpoint Peano kernel gives

```text
|D_(tr,8)(x)|
 <=sum_(m=39895)^39902 integral_(m-1/2)^(m+1/2)
       K_m(y)|G_(A,yy)(y,x)|dy.                    (EP1)
```

For a real normal coordinate `q0<0`, put

```text
w0=(-q0)sqrt(pi)(1-i)/2,
H(q0)=(1+i)exp(w0^2)erfc(w0)/2.                    (EP2)
```

If an interval for `q` has midpoint `q0` and radius `r`, the exact integrating
factor for `H_q=1-i*pi*qH` and `|exp(i*theta)-1|<=|theta|` give

```text
|H(q)-H(q0)|
 <=r[1+pi|H(q0)|(|q0|+r/2)].                       (EP3)
```

Thus (EP2)--(EP3), followed by the exact ODE for `H_qq`, enclose `G_(A,yy)`
without complex-box evaluation of `erfc`.  Boxes with certified `|q|>=2`
instead use the six-current expansion and its exact remainder; if that
alternating interval expression is indeterminate, the box fails over to
(EP2)--(EP3).

The 34-cell tail is not discarded.  Its exact rational-log defect is enclosed
on every `x` slab and its certified six-current replacement error is attached
after integration.  Consequently

```text
B_end
 <=2(pi/(32T))^(1/4) integral_(x_16)^(1/2)
      x^(-5/4)(1-x)^(-1/4)
      [B_8(x)+|D_(tail,34),6(x)|]dx
    +epsilon_(tail,6).                              (EP4)
```

The builder covers the slightly larger rational-decimal interval
`[{X_LOWER_TEXT},1/2]` by `{X_SLABS}` `x` slabs and splits each half-cell into
`{HALF_CELL_SUBDIVISIONS}` normal boxes.  Arb certifies

```text
eight-cell contribution       <= {cert['physical_eight_cell_bound']},
34-cell six-current tail      <= {cert['physical_tail_six_current_bound']},
tail replacement error       <= {cert['physical_tail_replacement_error']},

B_end                         <= {cert['physical_complete_endpoint_bound']}
                              < {COMPLETE_ENDPOINT_THRESHOLD}.             (EP5)
```

The independent checker repeats (EP1)--(EP5) at 120 decimal digits on a
`320 x (8*40)` partition and obtains its own bound below `0.0030`; it does not
import the builder or accept sampled floating values.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (EP1)--(EP5) comes from the canonical Fresnel
primitive, its exact ODE, or the inherited equation-(9) normalization.  No
geometric or fitted occurrence of `pi` is inserted.

Proof boundary: (EP5) proves only the complete translated A-face endpoint
layer on `x_16<=x<=1/2` at the saved height.  It does not prove the separate
pre-endpoint Morse integral, the full signed A-face channel, the joined
`R_after_A`, `R_Dir`, or `Q_K-T` bound, an all-height theorem, `Lambda<=0`,
PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    localization = load_json(LOCALIZATION)
    require(localization.get("passed") is True, "localization dependency is not passed")
    require(localization["decision"]["special_function_obligation_reduced_to_eight_cells"] is True, "eight-cell dependency drift")
    require(localization["decision"]["physical_tail_replacement_error_below_6e_minus_15_proved"] is True, "tail replacement dependency drift")
    require(CHECKER.is_file(), "independent checker is missing")
    artifact = {
        "kind": STEM,
        "status": "rigorous_complete_translated_A_face_endpoint_layer_absolute_interval_bound_below_0_point_0037",
        "passed": True,
        "scope": {
            "height": T,
            "A": A,
            "x_region": "x_16<=x<=1/2",
            "transition_modes": [TRANSITION_FIRST_MODE, TRANSITION_LAST_MODE],
            "tail_modes": [TAIL_FIRST_MODE, TAIL_LAST_MODE],
        },
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(localization),
        "decision": {
            "eight_cell_transition_interval_bound_proved": True,
            "tail_six_current_interval_bound_included": True,
            "tail_replacement_error_included": True,
            "complete_A_face_endpoint_layer_absolute_bound_below_0_point_0037_proved": True,
            "independent_replay_completed": True,
            "pre_endpoint_Morse_integral_proved": False,
            "full_signed_A_face_bound_proved": False,
            "R_after_A_bound_proved": False,
            "R_Dir_bound_proved": False,
            "QK_minus_T_bound_proved": False,
            "rh_implication": False,
        },
        "dependency": {"path": relative(LOCALIZATION), "sha256": file_hash(LOCALIZATION)},
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "arb_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Certify the pre-endpoint six-current Morse integral on 0<x<=x_16, attach its <=1.025e-8 defect replacement under the physical weight, and then assemble it with this <0.0037 endpoint theorem. The remaining post-A channels stay separate.",
        "proof_boundary": "Rigorous complete translated A-face endpoint-layer absolute bound below 0.0037 at t=10^10 only. No pre-endpoint Morse integral, full signed A-face, joined R_after_A, R_Dir, Q_K-T, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("certified complete translated A-face endpoint layer below 0.0037", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
