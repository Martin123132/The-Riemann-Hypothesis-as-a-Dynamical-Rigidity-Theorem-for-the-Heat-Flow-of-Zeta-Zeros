#!/usr/bin/env python3
"""Interval-scout a nonoscillatory absolute bound for the A-face Morse core."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_core_absolute_peano_interval_scout"
DEPENDENCY = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_midpoint_defect_reduction_gate.json"
A = 159_577
T = 10_000_000_000
LEFT = mp.mpf("39894.5")
FIRST_MODE = 39_895
LAST_MODE = 39_936
Q_LOWER_TEXT = "0.0003884882198961"
Q_TAIL = 40
Q_SEGMENTS = (
    (Q_LOWER_TEXT, "0.001", 32),
    ("0.001", "0.01", 64),
    ("0.01", "0.1", 64),
    ("0.1", "1", 64),
    ("1", "5", 64),
    ("5", "10", 32),
    ("10", "20", 32),
    ("20", "40", 32),
)
TERMS = 6
PRECISION = 90


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


def odd_df(index: int) -> int:
    value = 1
    for factor in range(1, index + 1, 2):
        value *= factor
    return value


def interval(left: mp.mpf, right: mp.mpf) -> arb:
    middle = (left + right) / 2
    radius = (right - left) / 2 + mp.mpf("1e-95")
    return arb(mp.nstr(middle, 105), mp.nstr(radius, 105))


def g_six_yy_over_x_bound(y: arb, x: arb) -> arb:
    pi = arb.pi()
    delta = y - arb(A) * x / 2
    require(delta.lower() > 0, "normal denominator crossed zero")
    delta_floor = delta.lower()
    x_ceiling = x.upper()
    y_ceiling = y.upper()
    pi_floor = pi.lower()
    value = arb(A) / (pi_floor * delta_floor**3)
    for index in range(1, TERMS):
        power = 2 * index + 1
        coefficient_modulus_over_x = (
            arb(odd_df(2 * index - 1))
            * x_ceiling ** (index - 1)
            / (arb(2) ** index * pi_floor ** (index + 1))
        )
        value += coefficient_modulus_over_x * power * (
            arb(power + 1) * y_ceiling / delta_floor ** (power + 2)
            + 2 / delta_floor ** (power + 1)
        )
    return value


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--panel-scale", type=int, default=1)
    result.add_argument("--half-cell-subdivisions", type=int, default=2)
    return result


def main() -> int:
    args = parser().parse_args()
    require(args.panel_scale > 0 and args.half_cell_subdivisions > 0, "partition counts must be positive")
    started = time.time()
    priority = set_low_priority()
    ctx.dps = PRECISION
    ctx.threads = 1
    mp.mp.dps = 115
    dependency = json.loads(DEPENDENCY.read_text(encoding="utf-8"))
    require(dependency.get("passed") is True, "six-current dependency is not passed")
    require(dependency["decision"]["uniform_signed_defect_error_below_1_point_025e_minus_8_proved"] is True, "six-current error dependency drift")

    dy = mp.mpf("0.5") / args.half_cell_subdivisions
    normalization = 2 * (arb.pi() / (32 * arb(T))) ** (arb(1) / 4)
    core_six_bound = arb(0)
    replacement_bound = arb(0)
    maximum_peano_over_x = arb(0)
    sample_panels: dict[str, dict[str, str]] = {}
    replacement_constant = (
        arb(84)
        * 2
        * arb(odd_df(11))
        * (arb(79_789) / 2)
        / (arb.pi() ** 7 * arb(2) ** 6)
    )

    panel_index = 0
    panel_total = sum(count for _, _, count in Q_SEGMENTS) * args.panel_scale
    for segment_left_text, segment_right_text, base_panels in Q_SEGMENTS:
        segment_left = mp.mpf(segment_left_text)
        segment_right = mp.mpf(segment_right_text)
        panels = base_panels * args.panel_scale
        dq = (segment_right - segment_left) / panels
        for local_panel in range(panels):
            qb = interval(segment_left + local_panel * dq, segment_left + (local_panel + 1) * dq)
            xb = 1 / (1 + qb.exp())
            peano_over_x = arb(0)
            for mode in range(FIRST_MODE, LAST_MODE + 1):
                cell_left = mp.mpf(mode) - mp.mpf("0.5")
                for side in range(2):
                    side_left = cell_left + mp.mpf(side) / 2
                    for part in range(args.half_cell_subdivisions):
                        y0 = side_left + part * dy
                        y1 = y0 + dy
                        yb = interval(y0, y1)
                        if side == 0:
                            kernel = (yb - arb(str(cell_left))) ** 2 / 2
                        else:
                            kernel = (arb(str(cell_left + 1)) - yb) ** 2 / 2
                        derivative_bound = g_six_yy_over_x_bound(yb, xb)
                        require(
                            "nan" not in derivative_bound.str(20, more=True).lower(),
                            f"G6_yy/x enclosure failed at q={qb}, mode={mode}",
                        )
                        peano_over_x += arb(mp.nstr(dy, 105)) * kernel.upper() * derivative_bound.upper()

            x_ceiling = xb.upper()
            logistic_weight = (x_ceiling * (1 - x_ceiling)) ** (arb(3) / 4)
            core_six_bound += normalization * arb(mp.nstr(dq, 105)) * logistic_weight.upper() * peano_over_x
            delta_left = arb(79_789) / 2 - arb(A) * xb / 2
            delta_floor = delta_left.lower()
            replacement_integrand = (
                replacement_constant
                * x_ceiling ** (arb(23) / 4)
                / delta_floor**13
            )
            replacement_bound += normalization * arb(mp.nstr(dq, 105)) * replacement_integrand.upper()
            if peano_over_x.upper() > maximum_peano_over_x.upper():
                maximum_peano_over_x = peano_over_x
            if panel_index in (0, panel_total // 2, panel_total - 1):
                sample_panels[str(panel_index)] = {
                    "q_ball": qb.str(45, more=True),
                    "x_ball": xb.str(45, more=True),
                    "Peano_D6_over_x_bound": peano_over_x.str(45, more=True),
                }
            panel_index += 1

    q_tail = arb(Q_TAIL)
    x_tail_upper = (-q_tail).exp()
    delta_floor = arb(79_789) / 2 - arb(A) * x_tail_upper / 2
    y_ceiling = arb(79_873) / 2
    derivative_tail = arb(A) / (arb.pi() * delta_floor**3)
    for index in range(1, TERMS):
        power = 2 * index + 1
        coefficient = (
            arb(odd_df(2 * index - 1))
            * x_tail_upper ** (index - 1)
            / (arb(2) ** index * arb.pi() ** (index + 1))
        )
        derivative_tail += coefficient * power * (
            arb(power + 1) * y_ceiling / delta_floor ** (power + 2)
            + 2 / delta_floor ** (power + 1)
        )
    peano_tail = arb(7) * derivative_tail / 4
    q_tail_core = normalization * peano_tail * arb(4) * (-arb(3) * q_tail / 4).exp() / 3
    q_tail_replacement = (
        normalization
        * replacement_constant
        / delta_floor**13
        * arb(4)
        * (-arb(23) * q_tail / 4).exp()
        / 23
    )
    core_six_bound += q_tail_core
    replacement_bound += q_tail_replacement

    total = core_six_bound + replacement_bound
    suffix = f"qscale{args.panel_scale}_y{args.half_cell_subdivisions}"
    result_path = REPO_ROOT / f"work/rh_compute/results/{STEM}_{suffix}.json"
    artifact = {
        "kind": STEM,
        "status": "rigorous_interval_scout_of_nonoscillatory_six_current_core_absolute_bound_not_yet_promoted",
        "passed": True,
        "scope": {
            "height": T,
            "A": A,
            "x_region": ["0", "x_16"],
            "q_cover": [Q_LOWER_TEXT, "infinity"],
            "finite_q_tail_cut": Q_TAIL,
            "panel_scale": args.panel_scale,
            "finite_q_panels": panel_total,
            "half_cell_subdivisions": args.half_cell_subdivisions,
            "normal_boxes_per_q_panel": 42 * 2 * args.half_cell_subdivisions,
        },
        "bounds": {
            "maximum_Peano_D6_over_x_bound": maximum_peano_over_x.str(65, more=True),
            "physical_six_current_core_absolute_bound": core_six_bound.str(65, more=True),
            "physical_D42_minus_D6_replacement_bound": replacement_bound.str(65, more=True),
            "physical_complete_pre_endpoint_core_absolute_bound": total.str(65, more=True),
            "saved_R_after_A_target": "0.0368147039947",
        },
        "sample_panels": sample_panels,
        "decision": {
            "interval_arithmetic_used": True,
            "oscillatory_cancellation_used": False,
            "logistic_q_regularization_used": True,
            "analytic_q_tail_from_40_used": True,
            "complete_pre_endpoint_core_bound_promoted": False,
            "full_signed_A_face_bound_proved": False,
            "R_after_A_bound_proved": False,
            "rh_implication": False,
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "arb_threads": 1,
            "process_priority": priority,
            "precision_decimal_digits": PRECISION,
        },
        "dependency": {"path": DEPENDENCY.relative_to(REPO_ROOT).as_posix(), "sha256": file_hash(DEPENDENCY)},
        "source": {"path": Path(__file__).resolve().relative_to(REPO_ROOT).as_posix(), "sha256": file_hash(Path(__file__).resolve())},
        "proof_boundary": "Rigorous interval scout for an absolute pre-endpoint translated A-face core bound only. It is not independently replayed or promoted. No complete A-face assembly, R_after_A, R_Dir, Q_K-T, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    result_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"interval-scouted absolute A-face core: total {total.str(12, more=True)}; "
        f"six-current {core_six_bound.str(12, more=True)}; replacement {replacement_bound.str(8, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
