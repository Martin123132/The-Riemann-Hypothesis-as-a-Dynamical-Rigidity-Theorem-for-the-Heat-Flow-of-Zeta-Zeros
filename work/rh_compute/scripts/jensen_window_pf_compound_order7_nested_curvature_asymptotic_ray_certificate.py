#!/usr/bin/env python3
"""Certify order-seven first-summand curvature on the ray u>=20."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate as order6  # noqa: E402
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_scale,
    series_sub,
)
from jensen_window_pf_compound_order7_high_cumulant_coarse_corridor import (  # noqa: E402
    DEFAULT_OUT as SOURCE_HIGH_CUMULANTS,
    RAY_EXACT_CAP as HIGH_CUMULANT_CAP,
)
from jensen_window_pf_compound_order7_nested_curvature_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
)
from jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate import (  # noqa: E402
    DEFAULT_OUT as SOURCE_FINITE_RAY,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_upper_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate.md"
)
SOURCE_ORDER4_RAY = order6.SOURCE_ORDER4_RAY
SOURCE_MID_CUMULANTS = order6.SOURCE_HIGH_CUMULANTS
PRECISION_BITS = 384
INVERSE_T_CAP = order6.INVERSE_T_CAP
COLLAR_RADIUS = 5
COLLAR_MODE_START = 19
NORMALIZED_CAP = order6.NORMALIZED_CAP
B0_FLOOR = order6.B0_FLOOR
B1_MAGNITUDE_FLOOR = order6.B1_MAGNITUDE_FLOOR
MID_CUMULANT_CAP = order6.HIGH_CUMULANT_CAP
MID_B_CAP = 2500
HIGH_B_CAP = 40000
FIRST_COMPOSITION_CAP = 100000
NESTED_COMPOSITION_CAP = 10**12
STABLE_LOG_ARGUMENT_CAP = Fraction(1, 1000)
SCALED_TARGET = 100


@dataclass(frozen=True)
class RayRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def source_contract() -> dict:
    order4 = load_json(SOURCE_ORDER4_RAY)
    middle = load_json(SOURCE_MID_CUMULANTS)
    high = load_json(SOURCE_HIGH_CUMULANTS)
    finite = load_json(SOURCE_FINITE_RAY)
    if order4.get("summary", {}).get("ray_corridor_to_curvature_closed") is not True:
        raise RuntimeError("normalized-H ray source is not closed")
    normalized = order4.get("normalized_h_boxes", {})
    if normalized.get("consequence") != (
        "0<x_r<=1 for 2<=r<=8, x_2>=97/100, x_3>=24/25"
    ):
        raise RuntimeError("normalized-H box contract changed")
    geometry = order4.get("geometry", {})
    if geometry.get("inverse_slope_cap") != str(INVERSE_T_CAP):
        raise RuntimeError("inverse-t ray cap changed")
    if middle.get("exact", {}).get("exact_corridor_cap") != MID_CUMULANT_CAP:
        raise RuntimeError("ninth/tenth cumulant source changed")
    if high.get("exact", {}).get("ray_exact_corridor_cap") != HIGH_CUMULANT_CAP:
        raise RuntimeError("eleventh/twelfth cumulant source changed")
    if finite.get("finite_ray", {}).get("mode_range") != ["2", "20"]:
        raise RuntimeError("order-seven finite-ray source changed")
    return {
        "normalized_H_boxes": normalized["consequence"],
        "inverse_t_cap": geometry["inverse_slope_cap"],
        "collar_t20_minus_t19_lower": geometry[
            "collar_t20_minus_t19_lower"
        ],
        "mid_cumulant_corridor": middle["exact"]["exact_corridor"],
        "high_cumulant_corridor": high["exact"]["ray_exact_corridor"],
        "finite_ray": "r_1''(t)<=600/t^2 for every saddle mode 2<=u<=20",
    }


def raw_high_normalized_cap(order: int, cumulant_cap: int) -> Fraction:
    geometry = (
        Fraction(201, 100) ** (order - 1)
        / Fraction(19, 10) ** order
    )
    return Fraction(1) + cumulant_cap * geometry / COLLAR_MODE_START


def exact_ratio_gates() -> dict:
    z = INVERSE_T_CAP
    low_order_upper_ratio = Fraction(1) / (1 - COLLAR_RADIUS * z) ** 7
    mid_order_upper_ratio = Fraction(1) / (1 - COLLAR_RADIUS * z) ** 9
    high_order_upper_ratio = Fraction(1) / (1 - COLLAR_RADIUS * z) ** 11
    lower_ratio = Fraction(1) / (1 + COLLAR_RADIUS * z) ** 2
    for name, ratio in (
        ("low", low_order_upper_ratio),
        ("mid", mid_order_upper_ratio),
        ("high", high_order_upper_ratio),
    ):
        if not ratio < NORMALIZED_CAP:
            raise RuntimeError(f"{name}-order collar ratio gate failed")
    if not lower_ratio > Fraction(999, 1000):
        raise RuntimeError("lower collar ratio gate failed")
    if not Fraction(97, 100) * lower_ratio > B0_FLOOR:
        raise RuntimeError("normalized B0 floor failed")
    if not Fraction(24, 25) * lower_ratio > B1_MAGNITUDE_FLOOR:
        raise RuntimeError("normalized B1 floor failed")

    high_caps = {}
    for order, cumulant_cap, ratio, target in (
        (9, MID_CUMULANT_CAP, mid_order_upper_ratio, MID_B_CAP),
        (10, MID_CUMULANT_CAP, mid_order_upper_ratio, MID_B_CAP),
        (11, HIGH_CUMULANT_CAP, high_order_upper_ratio, HIGH_B_CAP),
        (12, HIGH_CUMULANT_CAP, high_order_upper_ratio, HIGH_B_CAP),
    ):
        raw = raw_high_normalized_cap(order, cumulant_cap)
        shifted = raw * ratio
        if not shifted < target:
            raise RuntimeError(f"high normalized derivative cap failed at {order}")
        high_caps[str(order)] = {
            "raw_cap": str(raw),
            "shifted_cap": str(shifted),
            "target_cap": target,
            "cumulant_cap": cumulant_cap,
        }

    collar = Fraction(
        "684140403277229664153888231548284194636517090552223519946531E-23"
    )
    if collar <= COLLAR_RADIUS:
        raise RuntimeError("mode collar does not cover t+-5")
    return {
        "low_order_upper_ratio": str(low_order_upper_ratio),
        "mid_order_upper_ratio": str(mid_order_upper_ratio),
        "high_order_upper_ratio": str(high_order_upper_ratio),
        "upper_ratio_cap": str(NORMALIZED_CAP),
        "lower_ratio": str(lower_ratio),
        "lower_ratio_floor": "999/1000",
        "B0_floor": str(B0_FLOOR),
        "B1_magnitude_floor": str(B1_MAGNITUDE_FLOOR),
        "high_normalized_caps": high_caps,
        "collar_radius": COLLAR_RADIUS,
        "all_strict": True,
    }


def dimensionless_interval() -> dict:
    flint.ctx.prec = PRECISION_BITS
    zero = Fraction(0)
    b = [
        order6.arb_interval(B0_FLOOR, NORMALIZED_CAP),
        order6.arb_interval(-NORMALIZED_CAP, -B1_MAGNITUDE_FLOOR),
    ]
    for order in range(2, 7):
        b.append(
            order6.arb_interval(zero, NORMALIZED_CAP)
            if order % 2 == 0
            else order6.arb_interval(-NORMALIZED_CAP, zero)
        )
    b.extend(
        [
            order6.arb_interval(Fraction(-MID_B_CAP), Fraction(MID_B_CAP)),
            order6.arb_interval(Fraction(-MID_B_CAP), Fraction(MID_B_CAP)),
            order6.arb_interval(Fraction(-HIGH_B_CAP), Fraction(HIGH_B_CAP)),
            order6.arb_interval(Fraction(-HIGH_B_CAP), Fraction(HIGH_B_CAP)),
        ]
    )

    z_cap = flint.arb(INVERSE_T_CAP.numerator) / INVERSE_T_CAP.denominator
    z = z_cap / 2 + flint.arb(0, z_cap / 2)
    ell = order6.stable_normalized_log(b, 10, FIRST_COMPOSITION_CAP)

    j_hat = [
        2 * b[order]
        - z * (order + 1) * (order + 2) * ell[order + 2]
        for order in range(9)
    ]
    h = series_add(
        series_scale(ell[:9], 2),
        order6.stable_normalized_log(j_hat, 8, NESTED_COMPOSITION_CAP),
    )

    r_hat = [
        3 * b[order]
        - z * (order + 1) * (order + 2) * h[order + 2]
        for order in range(7)
    ]
    q_layer = series_add(
        series_sub(series_scale(h[:7], 2), ell[:7]),
        order6.stable_normalized_log(r_hat, 6, NESTED_COMPOSITION_CAP),
    )

    s_hat = [
        4 * b[order]
        - z * (order + 1) * (order + 2) * q_layer[order + 2]
        for order in range(5)
    ]
    p = series_add(
        series_sub(series_scale(q_layer[:5], 2), h[:5]),
        order6.stable_normalized_log(s_hat, 4, NESTED_COMPOSITION_CAP),
    )

    t_hat = [
        5 * b[order]
        - z * (order + 1) * (order + 2) * p[order + 2]
        for order in range(3)
    ]
    r_layer = series_add(
        series_sub(series_scale(p[:3], 2), q_layer[:3]),
        order6.stable_normalized_log(t_hat, 2, NESTED_COMPOSITION_CAP),
    )
    scaled_curvature = 2 * r_layer[2]

    if not bool(
        j_hat[0] > 0
        and r_hat[0] > 0
        and s_hat[0] > 0
        and t_hat[0] > 0
    ):
        raise RuntimeError("dimensionless stable-coordinate floor failed")
    if not bool(scaled_curvature < SCALED_TARGET):
        raise RuntimeError(
            "dimensionless scaled-curvature target failed: "
            f"{scaled_curvature}"
        )
    if not bool(
        z_cap * FIRST_COMPOSITION_CAP < flint.arb(1) / 1000
        and z_cap * NESTED_COMPOSITION_CAP < flint.arb(1) / 1000
    ):
        raise RuntimeError("stable-log argument cap failed")

    def coefficient_rows(values: list[flint.arb]) -> list[str]:
        return [value.str(40).replace("e", "E") for value in values]

    return {
        "scaled_variable": "z=1/t in (0,10^-30]",
        "normalized_B_definition": (
            "b_r=t^(r+1)*B^(r)/r!, so the y=s/t Taylor jet of B is z*b"
        ),
        "normalized_B_boxes": coefficient_rows(b),
        "normalized_ell_identity": (
            "ell=log(z)+log(b)+R(z*b), R(w)=log((1-exp(-w))/w)"
        ),
        "normalized_J_boxes": coefficient_rows(j_hat),
        "normalized_R_boxes": coefficient_rows(r_hat),
        "normalized_S_boxes": coefficient_rows(s_hat),
        "normalized_T_boxes": coefficient_rows(t_hat),
        "scaled_curvature_ball": scaled_curvature.str(60).replace("e", "E"),
        "scaled_curvature_upper": arb_upper_text(scaled_curvature),
        "scaled_target": SCALED_TARGET,
        "scaled_margin_lower": arb_lower_text(
            flint.arb(SCALED_TARGET) - scaled_curvature
        ),
        "j0_lower": arb_lower_text(j_hat[0]),
        "r0_lower": arb_lower_text(r_hat[0]),
        "s0_lower": arb_lower_text(s_hat[0]),
        "t0_lower": arb_lower_text(t_hat[0]),
        "first_composition_coefficient_cap": FIRST_COMPOSITION_CAP,
        "nested_composition_coefficient_cap": NESTED_COMPOSITION_CAP,
        "stable_log_argument_cap": str(STABLE_LOG_ARGUMENT_CAP),
        "bell_correction_bound": (
            "coefficient r error <=2^(r-1)*C*10^-30 for r>=1 "
            "and <=C*10^-30 for r=0"
        ),
    }


def build_artifact() -> dict:
    sources = source_contract()
    ratios = exact_ratio_gates()
    interval = dimensionless_interval()
    rows = [
        RayRow(
            "co7ncarc_01_shifted_low_H_boxes",
            "exact_theorem_composition",
            "ready_to_apply",
            "The signed normalized H-jet boxes through order eight remain valid on the full five-unit collar.",
            "0<x_r<=1 (2<=r<=8), x_2>=97/100, x_3>=24/25",
            "Existing exact ray theorem plus strict rational collar arithmetic.",
        ),
        RayRow(
            "co7ncarc_02_high_H_boxes",
            "exact_analytic_consequence",
            "ready_to_apply",
            "The coarse exact cumulant corridors give absolute normalized boxes for H^(9) through H^(12).",
            "|b_7|,|b_8|<2500 and |b_9|,|b_10|<40000",
            "Exact ray geometry and the 50000 and 700001 cumulant caps.",
            ratios["high_normalized_caps"],
        ),
        RayRow(
            "co7ncarc_03_normalized_B_boxes",
            "exact_interval_lemma",
            "ready_to_apply",
            "Every tent-averaged B derivative lies in one dimensionless box through order ten.",
            "alternating unit boxes through b_6 with four coarse high boxes",
            "Unit-mass tent averaging of shifted normalized H boxes.",
        ),
        RayRow(
            "co7ncarc_04_stable_log_majorant",
            "exact_analytic_lemma",
            "ready_to_apply",
            "The analytic stable-log correction is uniformly negligible at all five logarithmic evaluations.",
            interval["bell_correction_bound"],
            "Convergent defect series and the exact partial-Bell identity.",
        ),
        RayRow(
            "co7ncarc_05_dimensionless_interval",
            "interval_analytic_theorem",
            "ready_to_apply",
            "One outward-rounded dimensionless jet box proves the asymptotic order-seven curvature ceiling.",
            "t^2*r_1''(t)<100 for every mode u>=20",
            "Uniform interval theorem over z in [0,10^-30], not point sampling.",
            {
                "scaled_upper": interval["scaled_curvature_upper"],
                "scaled_margin_lower": interval["scaled_margin_lower"],
                "j0_lower": interval["j0_lower"],
                "r0_lower": interval["r0_lower"],
                "s0_lower": interval["s0_lower"],
                "t0_lower": interval["t0_lower"],
            },
        ),
        RayRow(
            "co7ncarc_06_asymptotic_consequence",
            "exact_theorem_composition",
            "ready_to_apply",
            "The asymptotic theorem is stronger than the continuous target needed by the order-seven bridge.",
            "t^2*r_1''(t)<100<600 for u>=20",
            "First Newman summand at lambda=-100 only.",
        ),
        RayRow(
            "co7ncarc_07_global_handoff",
            "exact_handoff",
            "ready_to_apply",
            "The shifted-jet, compact, finite-ray, and asymptotic certificates cover every real t>=320.",
            "[320<=t<=1000] union [1000<=t<=V'(2)] union [2<=u<=20] union [u>=20]",
            "The full-kernel composition is recorded separately.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate",
        "date": "2026-07-14",
        "status": "rigorous nested first-summand order-seven curvature theorem on the complete u>=20 ray",
        "proof_boundary": (
            "This artifact proves r_1''(t)<=600/t^2 on u>=20. It does not "
            "by itself compose the compact ranges, full-kernel transfer, order-seven "
            "entry, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.md",
            "outputs/jensen_window_pf_compound_order7_high_cumulant_coarse_corridor.md",
            "outputs/jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.md",
            "outputs/formal_core.md",
        ],
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate.py",
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "inverse_t_cap": str(INVERSE_T_CAP),
            "collar_radius": COLLAR_RADIUS,
            "normalized_cap": str(NORMALIZED_CAP),
            "B0_floor": str(B0_FLOOR),
            "B1_magnitude_floor": str(B1_MAGNITUDE_FLOOR),
            "mid_cumulant_cap": MID_CUMULANT_CAP,
            "high_cumulant_cap": HIGH_CUMULANT_CAP,
            "mid_B_cap": MID_B_CAP,
            "high_B_cap": HIGH_B_CAP,
            "first_composition_cap": FIRST_COMPOSITION_CAP,
            "nested_composition_cap": NESTED_COMPOSITION_CAP,
            "scaled_target": SCALED_TARGET,
        },
        "source_contract": sources,
        "ratio_gates": ratios,
        "dimensionless_interval": interval,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "normalized_H_box_theorems": 2,
            "stable_log_majorants": 1,
            "dimensionless_interval_theorems": 1,
            "asymptotic_curvature_theorems": 1,
            "open_rows": 0,
        },
    }


def write_note(path: Path, artifact: dict) -> None:
    ratios = artifact["ratio_gates"]
    interval = artifact["dimensionless_interval"]
    lines = [
        "# Jensen-Window PF Compound Order-Seven Nested Curvature Asymptotic-Ray Certificate",
        "",
        "Date: 2026-07-14",
        "",
        "Status: rigorous nested first-summand order-seven curvature theorem on",
        "the complete `u>=20` ray. This certificate is not a proof of full",
        "order-seven entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "## Shifted H Boxes",
        "",
        "The completed order-four ray theorem gives signed unit boxes through",
        "`H^(8)`. Exact coarse corridors add",
        "",
        "```text",
        "|kappa_r|*q^(r/2-1)/(r-2)!<50000, r=9,10,",
        "|kappa_r|*q^(r/2-1)/(r-2)!<700001, r=11,12.",
        "```",
        "",
        "The exact ray geometry and Hurwitz midpoint ceiling give, throughout",
        "the `t+-5` collar,",
        "",
        "```text",
        f"raw normalized H11 cap={ratios['high_normalized_caps']['11']['raw_cap']},",
        f"raw normalized H12 cap={ratios['high_normalized_caps']['12']['raw_cap']},",
        "|b_7|,|b_8|<2500, |b_9|,|b_10|<40000.",
        "```",
        "",
        "## Dimensionless Stable Recursion",
        "",
        "Put `z=1/t`. Every stable logarithm is evaluated as",
        "",
        "```text",
        "log(1-exp(-z*v))=log(z)+log(v)+R(z*v),",
        "R(w)=log((1-exp(-w))/w).",
        "```",
        "",
        "The convergent defect series and exact partial-Bell identity give",
        "",
        "```text",
        interval["bell_correction_bound"],
        "```",
        "",
        "for every layer. A single outward-rounded interval evaluation over",
        "`0<=z<=10^-30` gives",
        "",
        "```text",
        f"j_0 lower={interval['j0_lower']}",
        f"r_0 lower={interval['r0_lower']}",
        f"s_0 lower={interval['s0_lower']}",
        f"t_0 lower={interval['t0_lower']}",
        f"t^2*r_1'' ball={interval['scaled_curvature_ball']}",
        f"t^2*r_1'' upper={interval['scaled_curvature_upper']}<100<600,",
        f"margin below 100={interval['scaled_margin_lower']}.",
        "```",
        "",
        "This is a uniform interval theorem, not evaluation at a representative",
        "large value of `t`. Together with the preceding certificates it covers",
        "every real `t>=320` for the first Newman summand.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order7_nested_curvature_compact_certificate.md",
        "outputs/jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.md",
        "outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-seven nested curvature asymptotic-ray certificate: "
        f"{summary['rows']} rows, "
        f"{summary['dimensionless_interval_theorems']} dimensionless interval theorem, "
        f"{summary['asymptotic_curvature_theorems']} asymptotic curvature theorem, "
        f"{summary['open_rows']} open rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
