#!/usr/bin/env python3
"""Independently replay the nonoscillatory translated A-face core bound."""

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
import mpmath as mp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_core_absolute_peano_interval_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
A = 159_577
T = 10_000_000_000
Q_LOWER_TEXT = "0.0003884882198961"
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
PANEL_SCALE = 3
HALF_SUBDIVISIONS = 3
TERMS = 6


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def set_low_priority() -> None:
    try:
        import psutil

        process = psutil.Process()
        process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS if os.name == "nt" else 10)
    except Exception:
        pass


def odd_df(index: int) -> int:
    answer = 1
    for value in range(1, index + 1, 2):
        answer *= value
    return answer


def make_box(left: mp.mpf, right: mp.mpf) -> arb:
    midpoint = (left + right) / 2
    radius = (right - left) / 2 + mp.mpf("1e-125")
    return arb(mp.nstr(midpoint, 135), mp.nstr(radius, 135))


def independent_derivative_majorant(y: arb, x: arb) -> arb:
    d = y - arb(A) * x / 2
    require(d.lower() > 0, "checker normal denominator crossed zero")
    d0 = d.lower()
    xx = x.upper()
    yy = y.upper()
    pp = arb.pi().lower()
    bound = arb(A) / (pp * d0**3)
    for n in range(1, TERMS):
        p = 2 * n + 1
        prefactor = arb(odd_df(2 * n - 1)) * xx ** (n - 1) / (arb(2) ** n * pp ** (n + 1))
        first = arb(p * (p + 1)) * yy / d0 ** (p + 2)
        second = arb(2 * p) / d0 ** (p + 1)
        bound += prefactor * (first + second)
    require("nan" not in bound.str(24, more=True).lower(), "checker derivative majorant failed")
    return bound


def replay(local_dependency: dict) -> arb:
    ctx.dps = 120
    ctx.threads = 1
    mp.mp.dps = 150
    pi = arb.pi()
    normalization = 2 * (pi / (32 * arb(T))) ** (arb(1) / 4)
    left = arb(79_789) / 2
    right = arb(79_873) / 2
    dy = mp.mpf("0.5") / HALF_SUBDIVISIONS
    core = arb(0)
    replacement = arb(0)
    replacement_constant = arb(84 * 2 * odd_df(11)) * left / (pi**7 * arb(2) ** 6)

    for left_text, right_text, base_panels in Q_SEGMENTS:
        qa = mp.mpf(left_text)
        qb = mp.mpf(right_text)
        panels = base_panels * PANEL_SCALE
        dq = (qb - qa) / panels
        for panel in range(panels):
            qbox = make_box(qa + panel * dq, qa + (panel + 1) * dq)
            xbox = 1 / (1 + qbox.exp())
            peano_over_x = arb(0)
            for mode in range(39_895, 39_937):
                cell_left = mp.mpf(mode) - mp.mpf("0.5")
                for side in range(2):
                    side_start = cell_left + mp.mpf(side) / 2
                    for part in range(HALF_SUBDIVISIONS):
                        ybox = make_box(side_start + part * dy, side_start + (part + 1) * dy)
                        if side == 0:
                            kernel_maximum = ((part + 1) * dy) ** 2 / 2
                        else:
                            kernel_maximum = ((HALF_SUBDIVISIONS - part) * dy) ** 2 / 2
                        peano_over_x += (
                            arb(mp.nstr(dy * kernel_maximum, 135))
                            * independent_derivative_majorant(ybox, xbox).upper()
                        )
            x_max = xbox.upper()
            logistic_weight = (x_max * (1 - x_max)) ** (arb(3) / 4)
            width = arb(mp.nstr(dq, 135))
            core += normalization * width * logistic_weight.upper() * peano_over_x
            delta_floor = (left - arb(A) * xbox / 2).lower()
            replacement += normalization * width * replacement_constant * x_max ** (arb(23) / 4) / delta_floor**13

    q0 = arb(40)
    x0 = (-q0).exp()
    delta0 = left - arb(A) * x0 / 2
    derivative = arb(A) / (pi.lower() * delta0**3)
    for n in range(1, TERMS):
        p = 2 * n + 1
        prefactor = arb(odd_df(2 * n - 1)) * x0 ** (n - 1) / (arb(2) ** n * pi.lower() ** (n + 1))
        derivative += prefactor * (
            arb(p * (p + 1)) * right / delta0 ** (p + 2)
            + arb(2 * p) / delta0 ** (p + 1)
        )
    core += normalization * arb(7) * derivative / 4 * arb(4) * (-arb(3) * q0 / 4).exp() / 3
    replacement += (
        normalization
        * replacement_constant
        / delta0**13
        * arb(4)
        * (-arb(23) * q0 / 4).exp()
        / 23
    )
    total = core + replacement
    require(core.upper() < arb("0.000012"), "checker six-current core threshold failed")
    require(replacement.upper() < arb("2e-15"), "checker replacement threshold failed")
    require(total.upper() < arb("0.000012"), "checker complete core threshold failed")
    require(local_dependency["decision"]["six_current_Fresnel_tail_with_finite_remainder_proved"] is True, "dependency theorem drift")
    return total


def main() -> int:
    set_low_priority()
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decision = artifact["decision"]
    for key in (
        "positive_Peano_kernel_used",
        "lower_endpoint_regularized_rigorously",
        "six_current_replacement_integrated_rigorously",
        "complete_pre_endpoint_A_face_core_absolute_bound_below_1_point_2e_minus_5_proved",
        "independent_replay_completed",
    ):
        require(decision.get(key) is True, f"missing proved decision: {key}")
    for key in (
        "oscillatory_Morse_cancellation_used",
        "endpoint_layer_included",
        "full_signed_A_face_bound_proved",
        "R_after_A_bound_proved",
        "R_Dir_bound_proved",
        "QK_minus_T_bound_proved",
        "rh_implication",
    ):
        require(decision.get(key) is False, f"proof-boundary drift: {key}")

    cert = artifact["interval_certificate"]
    require(arb(cert["physical_complete_pre_endpoint_core_absolute_bound"]).upper() < arb("0.000012"), "builder threshold drift")
    require(cert["partition"]["panel_scale"] == 2, "builder panel scale drift")
    require(cert["partition"]["half_cell_subdivisions"] == 4, "builder normal partition drift")

    dependency_record = artifact["dependency"]
    dependency_path = REPO_ROOT / dependency_record["path"]
    require(file_hash(dependency_path) == dependency_record["sha256"], "dependency hash drift")
    independent = replay(load_json(dependency_path))
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("No subtraction of the large point sum" in note, "stable Peano route missing")
    require("< 0.000012" in note, "core theorem threshold missing")
    require("does not by itself" in note and "RH" in note, "proof boundary missing")
    print(
        "independently checked nonoscillatory translated A-face core bound; "
        f"replay {independent.str(12, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
