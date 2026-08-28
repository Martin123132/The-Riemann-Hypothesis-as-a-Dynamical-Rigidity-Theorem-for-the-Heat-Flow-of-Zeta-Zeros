#!/usr/bin/env python3
"""Interval-scout an absolute bound for the eight-cell endpoint transition."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx
import mpmath as mp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_eight_cell_peano_interval_scout"
DEPENDENCY = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_transition_localization_gate.json"

A = 159_577
T = 10_000_000_000
LEFT = mp.mpf("39894.5")
FIRST_MODE = 39_895
LAST_MODE = 39_902
TAIL_LEFT_NUMERATOR = 79_805
TAIL_RIGHT_NUMERATOR = 79_873
TAIL_FIRST_MODE = 39_903
TAIL_LAST_MODE = 39_936
X_LOWER_TEXT = "0.4999028779462474"
X_UPPER_TEXT = "0.5"
PRECISION = 70
ASYMPTOTIC_Q_FLOOR = 2
TERMS = 6


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


def odd_double_factorial(index: int) -> int:
    result = 1
    for value in range(1, index + 1, 2):
        result *= value
    return result


def interval_from_endpoints(left: mp.mpf, right: mp.mpf) -> arb:
    midpoint = (left + right) / 2
    radius = (right - left) / 2 + mp.mpf("1e-75")
    return arb(mp.nstr(midpoint, 85), mp.nstr(radius, 85))


def zero_centered_upper_ball(value: arb) -> arb:
    midpoint, radius, exponent = value.upper().mid_rad_10exp()
    return arb("0", f"{midpoint + radius}e{exponent}")


def is_finite(value: arb | acb) -> bool:
    return "nan" not in value.str(20, more=True).lower()


def require_finite(value: arb | acb, label: str, *, y: arb, x: arb, q: arb) -> None:
    require(
        is_finite(value),
        f"{label} became indeterminate for x={x}, y={y}, q={q}",
    )


def exact_or_asymptotic_gyy_bound(y: arb, x: arb) -> tuple[arb, str]:
    pi = arb.pi()
    c = (2 / x).sqrt()
    q = -(y - arb(A) * x / 2) * c
    q_abs = -q
    require(q_abs.lower() > arb(0), "normal interval crossed q=0")

    if q_abs.lower() >= arb(ASYMPTOTIC_Q_FLOOR):
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
        correction = (
            arb(10_395)
            * x**4
            / (16 * pi**6 * delta**13)
            * (acb(0, pi * delta**2 * y) - acb(delta * x) + acb(6 * x * y))
        )
        approximation += correction
        coefficient = -acb(0, 2 * pi * c**2 * q) + acb(y * c**3) * (acb(0, pi) + acb(pi**2 * q**2))
        h_remainder = 2 * arb(odd_double_factorial(11)) / (pi**7 * q_abs**13)
        if is_finite(approximation) and is_finite(coefficient) and is_finite(h_remainder):
            approximation_modulus = abs(approximation)
            coefficient_modulus = abs(coefficient)
            bound = approximation_modulus + coefficient_modulus * h_remainder
            if is_finite(approximation_modulus) and is_finite(coefficient_modulus) and is_finite(bound):
                return bound, "asymptotic"

    q_mid = q.mid()
    q_radius = q.rad()
    scale_mid = (-q_mid) * pi.sqrt() / 2
    w_mid = acb(scale_mid, -scale_mid)
    h_mid = acb(arb(1) / 2, arb(1) / 2) * (w_mid**2).exp() * w_mid.erfc()
    require_finite(h_mid, "center H evaluation", y=y, x=x, q=q)
    h_variation = q_radius * (
        1 + pi * abs(h_mid) * (abs(q_mid) + q_radius / 2)
    )
    component_error = zero_centered_upper_ball(h_variation)
    h = acb(h_mid.real + component_error, h_mid.imag + component_error)
    h_q = 1 - acb(0, pi) * acb(q) * h
    h_qq = -acb(0, pi) * acb(q) - (acb(0, pi) + acb(pi**2 * q**2)) * h
    g_yy = acb(2 * c**2) * h_q - acb(y * c**3) * h_qq
    require_finite(g_yy, "G_yy enclosure", y=y, x=x, q=q)
    bound = abs(g_yy)
    require_finite(bound, "exact G_yy bound", y=y, x=x, q=q)
    return bound, "exact_erfc"


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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--x-slabs", type=int, default=128)
    parser.add_argument("--half-cell-subdivisions", type=int, default=8)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    require(args.x_slabs > 0 and args.half_cell_subdivisions > 0, "partition counts must be positive")
    started = time.time()
    priority = set_low_priority()
    ctx.dps = PRECISION
    ctx.threads = 1
    mp.mp.dps = 100

    dependency = json.loads(DEPENDENCY.read_text(encoding="utf-8"))
    require(dependency.get("passed") is True, "eight-cell localization dependency is not passed")
    require(dependency["decision"]["special_function_obligation_reduced_to_eight_cells"] is True, "eight-cell route drift")

    x_lower = mp.mpf(X_LOWER_TEXT)
    x_upper = mp.mpf(X_UPPER_TEXT)
    x_width = (x_upper - x_lower) / args.x_slabs
    half_subdivisions = args.half_cell_subdivisions
    y_width = mp.mpf("0.5") / half_subdivisions
    normalization = 2 * (arb.pi() / (32 * arb(T))) ** (arb(1) / 4)

    physical_eight_cell_bound = arb(0)
    physical_tail_six_current_bound = arb(0)
    maximum_defect_bound = arb(0)
    maximum_tail_six_current_defect_bound = arb(0)
    worst_slab = -1
    exact_boxes = 0
    asymptotic_boxes = 0
    slab_bounds: list[str] = []

    for slab in range(args.x_slabs):
        x0 = x_lower + slab * x_width
        x1 = x0 + x_width
        x_box = interval_from_endpoints(x0, x1)
        defect_bound = arb(0)

        for mode in range(FIRST_MODE, LAST_MODE + 1):
            cell_left = mp.mpf(mode) - mp.mpf("0.5")
            for side in range(2):
                side_left = cell_left + mp.mpf(side) / 2
                for subdivision in range(half_subdivisions):
                    y0 = side_left + subdivision * y_width
                    y1 = y0 + y_width
                    y_box = interval_from_endpoints(y0, y1)
                    if side == 0:
                        kernel = (y_box - arb(str(cell_left))) ** 2 / 2
                    else:
                        kernel = (arb(str(cell_left + 1)) - y_box) ** 2 / 2
                    gyy_bound, method = exact_or_asymptotic_gyy_bound(y_box, x_box)
                    if method == "exact_erfc":
                        exact_boxes += 1
                    else:
                        asymptotic_boxes += 1
                    defect_bound += arb(mp.nstr(y_width, 85)) * kernel.upper() * gyy_bound.upper()

        weight_box = 1 / (x_box * (x_box * (1 - x_box)) ** (arb(1) / 4))
        slab_scale = normalization * arb(mp.nstr(x_width, 85)) * weight_box.upper()
        physical_eight_cell_bound += slab_scale * defect_bound
        tail_six_current_defect_bound = abs(tail_defect_six(x_box))
        require_finite(tail_six_current_defect_bound, "tail six-current defect", y=arb(TAIL_LEFT_NUMERATOR) / 2, x=x_box, q=arb(-16))
        physical_tail_six_current_bound += slab_scale * tail_six_current_defect_bound
        slab_bounds.append(defect_bound.str(30, more=True))
        if defect_bound.upper() > maximum_defect_bound.upper():
            maximum_defect_bound = defect_bound
            worst_slab = slab
        if tail_six_current_defect_bound.upper() > maximum_tail_six_current_defect_bound.upper():
            maximum_tail_six_current_defect_bound = tail_six_current_defect_bound

    localization_tail_error = arb(dependency["numerical_certificate"]["physical_tail_replacement_error_ball"])
    complete_endpoint_bound = physical_eight_cell_bound + physical_tail_six_current_bound + localization_tail_error
    fits_saved_target = complete_endpoint_bound.upper() < arb("0.0368147039947")

    suffix = f"x{args.x_slabs}_y{args.half_cell_subdivisions}"
    result = REPO_ROOT / f"work/rh_compute/results/{STEM}_{suffix}.json"
    artifact = {
        "kind": STEM,
        "status": "rigorous_interval_scout_of_eight_cell_Peano_absolute_endpoint_bound_not_yet_promoted",
        "passed": True,
        "scope": {
            "height": T,
            "A": A,
            "x_cover": [X_LOWER_TEXT, X_UPPER_TEXT],
            "x_slabs": args.x_slabs,
            "half_cell_subdivisions": half_subdivisions,
            "total_y_boxes_per_x_slab": 8 * 2 * half_subdivisions,
            "precision_decimal_digits": PRECISION,
            "asymptotic_q_floor": ASYMPTOTIC_Q_FLOOR,
        },
        "bounds": {
            "maximum_eight_cell_defect_modulus_bound": maximum_defect_bound.str(60, more=True),
            "worst_x_slab_index": worst_slab,
            "maximum_tail_six_current_defect_modulus_bound": maximum_tail_six_current_defect_bound.str(60, more=True),
            "physical_eight_cell_absolute_bound": physical_eight_cell_bound.str(60, more=True),
            "physical_tail_six_current_absolute_bound": physical_tail_six_current_bound.str(60, more=True),
            "physical_34_cell_tail_replacement_error": localization_tail_error.str(60, more=True),
            "physical_complete_endpoint_absolute_bound": complete_endpoint_bound.str(60, more=True),
            "saved_R_after_A_target": "0.0368147039947",
        },
        "method_counts": {
            "exact_erfc_boxes": exact_boxes,
            "six_current_asymptotic_boxes": asymptotic_boxes,
        },
        "slab_defect_bounds": slab_bounds,
        "decision": {
            "interval_arithmetic_used": True,
            "positive_Peano_kernel_used": True,
            "saved_target_exceeded": not fits_saved_target,
            "independent_replay_completed": False,
            "complete_endpoint_layer_bound_promoted": False,
            "full_signed_A_face_bound_proved": False,
            "R_after_A_bound_proved": False,
            "rh_implication": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "arb_threads": 1,
            "process_priority": priority,
        },
        "dependency": {"path": str(DEPENDENCY.relative_to(REPO_ROOT)).replace("\\", "/"), "sha256": file_hash(DEPENDENCY)},
        "source": {"path": str(Path(__file__).resolve().relative_to(REPO_ROOT)).replace("\\", "/"), "sha256": file_hash(Path(__file__).resolve())},
        "proof_boundary": "Rigorous interval scout for an absolute bound on the eight-cell endpoint transition plus the certified 34-cell replacement error. It has no independent replay and is not promoted as the complete endpoint theorem. No full signed A-face, R_after_A, R_Dir, Q_K-T, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    result.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"interval-scouted complete endpoint: physical bound {complete_endpoint_bound.str(12, more=True)}; "
        f"eight-cell {physical_eight_cell_bound.str(12, more=True)}; "
        f"tail-six {physical_tail_six_current_bound.str(12, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
