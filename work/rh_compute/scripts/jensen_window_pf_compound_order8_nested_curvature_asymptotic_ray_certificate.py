#!/usr/bin/env python3
"""Certify order-eight first-summand curvature on the ray u>=20."""

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
import jensen_window_pf_compound_order7_nested_curvature_asymptotic_ray_certificate as order7  # noqa: E402
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_scale,
    series_sub,
)
from jensen_window_pf_compound_order8_high_cumulant_coarse_corridor import (  # noqa: E402
    DEFAULT_OUT as SOURCE_HIGH_CUMULANTS,
    RAY_EXACT_CAP as ULTRA_CUMULANT_CAP,
)
from jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate import (  # noqa: E402
    DEFAULT_OUT as SOURCE_FINITE_RAY,
)
from jensen_window_pf_compound_order8_nested_curvature_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_upper_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate.md"
)
SOURCE_ORDER4_RAY = order6.SOURCE_ORDER4_RAY
SOURCE_MID_CUMULANTS = order6.SOURCE_HIGH_CUMULANTS
SOURCE_TOP_CUMULANTS = order7.SOURCE_HIGH_CUMULANTS
PRECISION_BITS = 384
INVERSE_T_CAP = order6.INVERSE_T_CAP
COLLAR_RADIUS = 6
COLLAR_MODE_START = 19
NORMALIZED_CAP = order6.NORMALIZED_CAP
B0_FLOOR = order6.B0_FLOOR
B1_MAGNITUDE_FLOOR = order6.B1_MAGNITUDE_FLOOR
MID_CUMULANT_CAP = order6.HIGH_CUMULANT_CAP
TOP_CUMULANT_CAP = order7.HIGH_CUMULANT_CAP
MID_B_CAP = order7.MID_B_CAP
TOP_B_CAP = order7.HIGH_B_CAP
ULTRA_B_CAP = 2
FIRST_COMPOSITION_CAP = 100000
NESTED_COMPOSITION_CAP = 10**6
STABLE_LOG_ARGUMENT_CAP = Fraction(1, 1000)
SCALED_TARGET = 200


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
    top = load_json(SOURCE_TOP_CUMULANTS)
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
    if top.get("exact", {}).get("ray_exact_corridor_cap") != TOP_CUMULANT_CAP:
        raise RuntimeError("eleventh/twelfth cumulant source changed")
    if high.get("exact", {}).get("ray_exact_corridor_cap") != ULTRA_CUMULANT_CAP:
        raise RuntimeError("thirteenth/fourteenth cumulant source changed")
    if finite.get("finite_ray", {}).get("mode_range") != ["2", "20"]:
        raise RuntimeError("order-eight finite-ray source changed")
    return {
        "normalized_H_boxes": normalized["consequence"],
        "inverse_t_cap": geometry["inverse_slope_cap"],
        "collar_t20_minus_t19_lower": geometry[
            "collar_t20_minus_t19_lower"
        ],
        "mid_cumulant_corridor": middle["exact"]["exact_corridor"],
        "top_cumulant_corridor": top["exact"]["ray_exact_corridor"],
        "ultra_cumulant_corridor": high["exact"]["ray_exact_corridor"],
        "finite_ray": (
            "s_1''(t)<=4000/t^2 for every saddle mode 2<=u<=20"
        ),
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
    top_order_upper_ratio = Fraction(1) / (1 - COLLAR_RADIUS * z) ** 11
    ultra_order_upper_ratio = Fraction(1) / (1 - COLLAR_RADIUS * z) ** 13
    lower_ratio = Fraction(1) / (1 + COLLAR_RADIUS * z) ** 2
    for name, ratio in (
        ("low", low_order_upper_ratio),
        ("mid", mid_order_upper_ratio),
        ("top", top_order_upper_ratio),
        ("ultra", ultra_order_upper_ratio),
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
        (11, TOP_CUMULANT_CAP, top_order_upper_ratio, TOP_B_CAP),
        (12, TOP_CUMULANT_CAP, top_order_upper_ratio, TOP_B_CAP),
        (13, ULTRA_CUMULANT_CAP, ultra_order_upper_ratio, ULTRA_B_CAP),
        (14, ULTRA_CUMULANT_CAP, ultra_order_upper_ratio, ULTRA_B_CAP),
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
        raise RuntimeError("mode collar does not cover t+-6")
    return {
        "low_order_upper_ratio": str(low_order_upper_ratio),
        "mid_order_upper_ratio": str(mid_order_upper_ratio),
        "top_order_upper_ratio": str(top_order_upper_ratio),
        "ultra_order_upper_ratio": str(ultra_order_upper_ratio),
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
    for cap in (
        MID_B_CAP,
        MID_B_CAP,
        TOP_B_CAP,
        TOP_B_CAP,
        ULTRA_B_CAP,
        ULTRA_B_CAP,
    ):
        b.append(order6.arb_interval(Fraction(-cap), Fraction(cap)))

    z_cap = flint.arb(INVERSE_T_CAP.numerator) / INVERSE_T_CAP.denominator
    z = z_cap / 2 + flint.arb(0, z_cap / 2)
    ell = order6.stable_normalized_log(b, 12, FIRST_COMPOSITION_CAP)

    j_hat = [
        2 * b[order]
        - z * (order + 1) * (order + 2) * ell[order + 2]
        for order in range(11)
    ]
    h = series_add(
        series_scale(ell[:11], 2),
        order6.stable_normalized_log(
            j_hat,
            10,
            NESTED_COMPOSITION_CAP,
        ),
    )

    r_hat = [
        3 * b[order]
        - z * (order + 1) * (order + 2) * h[order + 2]
        for order in range(9)
    ]
    q_layer = series_add(
        series_sub(series_scale(h[:9], 2), ell[:9]),
        order6.stable_normalized_log(
            r_hat,
            8,
            NESTED_COMPOSITION_CAP,
        ),
    )

    s_hat = [
        4 * b[order]
        - z * (order + 1) * (order + 2) * q_layer[order + 2]
        for order in range(7)
    ]
    p = series_add(
        series_sub(series_scale(q_layer[:7], 2), h[:7]),
        order6.stable_normalized_log(
            s_hat,
            6,
            NESTED_COMPOSITION_CAP,
        ),
    )

    t_hat = [
        5 * b[order]
        - z * (order + 1) * (order + 2) * p[order + 2]
        for order in range(5)
    ]
    r_layer = series_add(
        series_sub(series_scale(p[:5], 2), q_layer[:5]),
        order6.stable_normalized_log(
            t_hat,
            4,
            NESTED_COMPOSITION_CAP,
        ),
    )

    u_hat = [
        6 * b[order]
        - z * (order + 1) * (order + 2) * r_layer[order + 2]
        for order in range(3)
    ]
    s_layer = series_add(
        series_sub(series_scale(r_layer[:3], 2), p[:3]),
        order6.stable_normalized_log(
            u_hat,
            2,
            NESTED_COMPOSITION_CAP,
        ),
    )
    scaled_curvature = 2 * s_layer[2]

    coordinates = {
        "j": j_hat[0],
        "r": r_hat[0],
        "s": s_hat[0],
        "t": t_hat[0],
        "u": u_hat[0],
    }
    if not bool(all(value > 0 for value in coordinates.values())):
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
        "normalized_U_boxes": coefficient_rows(u_hat),
        "scaled_curvature_ball": scaled_curvature.str(60).replace("e", "E"),
        "scaled_curvature_upper": arb_upper_text(scaled_curvature),
        "scaled_target": SCALED_TARGET,
        "scaled_margin_lower": arb_lower_text(
            flint.arb(SCALED_TARGET) - scaled_curvature
        ),
        **{
            f"{name}0_lower": arb_lower_text(value)
            for name, value in coordinates.items()
        },
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
            "co8ncarc_01_shifted_low_H_boxes",
            "exact_theorem_composition",
            "ready_to_apply",
            "The signed normalized H-jet boxes through order eight remain valid on the full six-unit collar.",
            "0<x_r<=1 (2<=r<=8), x_2>=97/100, x_3>=24/25",
            "Existing exact ray theorem plus strict rational collar arithmetic.",
        ),
        RayRow(
            "co8ncarc_02_high_H_boxes",
            "exact_analytic_consequence",
            "ready_to_apply",
            "The coarse exact cumulant corridors give absolute normalized boxes for H^(9) through H^(14).",
            "|b_7|,|b_8|<2500; |b_9|,|b_10|<40000; |b_11|,|b_12|<2",
            "Exact ray geometry and the cumulant caps 50000, 700001, and 1.",
            ratios["high_normalized_caps"],
        ),
        RayRow(
            "co8ncarc_03_normalized_B_boxes",
            "exact_interval_lemma",
            "ready_to_apply",
            "Every tent-averaged B derivative lies in one dimensionless box through order twelve.",
            "alternating unit boxes through b_6 with six coarse high boxes",
            "Unit-mass tent averaging of shifted normalized H boxes.",
        ),
        RayRow(
            "co8ncarc_04_stable_log_majorant",
            "exact_analytic_lemma",
            "ready_to_apply",
            "The analytic stable-log correction is uniformly controlled at all six logarithmic evaluations.",
            interval["bell_correction_bound"],
            "Convergent defect series and the exact partial-Bell identity.",
        ),
        RayRow(
            "co8ncarc_05_dimensionless_interval",
            "interval_analytic_theorem",
            "ready_to_apply",
            "One outward-rounded dimensionless jet box proves the asymptotic order-eight curvature ceiling.",
            "t^2*s_1''(t)<200 for every mode u>=20",
            "Uniform interval theorem over z in [0,10^-30], not point sampling.",
            {
                "scaled_upper": interval["scaled_curvature_upper"],
                "scaled_margin_lower": interval["scaled_margin_lower"],
                **{
                    f"{name}0_lower": interval[f"{name}0_lower"]
                    for name in ("j", "r", "s", "t", "u")
                },
            },
        ),
        RayRow(
            "co8ncarc_06_asymptotic_consequence",
            "exact_theorem_composition",
            "ready_to_apply",
            "The asymptotic theorem is stronger than the continuous target needed by the order-eight bridge.",
            "t^2*s_1''(t)<200<4000 for u>=20",
            "First Newman summand at lambda=-100 only.",
        ),
        RayRow(
            "co8ncarc_07_global_handoff",
            "exact_handoff",
            "ready_to_apply",
            "The compact, finite-ray, and asymptotic certificates cover every real t>=999.",
            "[999<=t<=V'(2)] union [2<=u<=20] union [u>=20]",
            "The full-kernel composition is recorded separately.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate",
        "date": "2026-07-13",
        "status": "rigorous nested first-summand order-eight curvature theorem on the complete u>=20 ray",
        "proof_boundary": (
            "This artifact proves s_1''(t)<=4000/t^2 on u>=20 and records "
            "the exact composition with the preceding first-summand ranges. It "
            "does not by itself prove the full-kernel transfer, order-eight "
            "entry, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.md",
            "outputs/jensen_window_pf_compound_order7_high_cumulant_coarse_corridor.md",
            "outputs/jensen_window_pf_compound_order8_high_cumulant_coarse_corridor.md",
            "outputs/jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate.md",
            "outputs/formal_core.md",
        ],
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate.py",
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "inverse_t_cap": str(INVERSE_T_CAP),
            "collar_radius": COLLAR_RADIUS,
            "normalized_cap": str(NORMALIZED_CAP),
            "B0_floor": str(B0_FLOOR),
            "B1_magnitude_floor": str(B1_MAGNITUDE_FLOOR),
            "mid_cumulant_cap": MID_CUMULANT_CAP,
            "top_cumulant_cap": TOP_CUMULANT_CAP,
            "ultra_cumulant_cap": ULTRA_CUMULANT_CAP,
            "mid_B_cap": MID_B_CAP,
            "top_B_cap": TOP_B_CAP,
            "ultra_B_cap": ULTRA_B_CAP,
            "first_composition_cap": FIRST_COMPOSITION_CAP,
            "nested_composition_cap": NESTED_COMPOSITION_CAP,
            "scaled_target": SCALED_TARGET,
            "continuous_target": CURVATURE_CONSTANT,
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
            "global_first_summand_handoffs": 1,
            "open_rows": 0,
        },
    }


def write_note(path: Path, artifact: dict) -> None:
    ratios = artifact["ratio_gates"]
    interval = artifact["dimensionless_interval"]
    lines = [
        "# Jensen-Window PF Compound Order-Eight Nested Curvature Asymptotic-Ray Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous nested first-summand order-eight curvature theorem on",
        "the complete `u>=20` ray. This certificate is not a proof of full",
        "order-eight entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "## Shifted H Boxes",
        "",
        "The completed order-four ray theorem gives signed unit boxes through",
        "`H^(8)`. Exact coarse corridors through order fourteen give",
        "",
        "```text",
        "|kappa_r|*q^(r/2-1)/(r-2)!<50000, r=9,10,",
        "|kappa_r|*q^(r/2-1)/(r-2)!<700001, r=11,12,",
        "|kappa_r|*q^(r/2-1)/(r-2)!<1, r=13,14.",
        "```",
        "",
        "The exact ray geometry gives, throughout the `t+-6` collar,",
        "",
        "```text",
        f"raw normalized H13 cap={ratios['high_normalized_caps']['13']['raw_cap']},",
        f"raw normalized H14 cap={ratios['high_normalized_caps']['14']['raw_cap']},",
        "|b_7|,|b_8|<2500, |b_9|,|b_10|<40000, |b_11|,|b_12|<2.",
        "```",
        "",
        "## Dimensionless Stable Recursion",
        "",
        "Every stable logarithm uses the convergent defect series and the",
        "exact partial-Bell identity in the form",
        "",
        "```text",
        "log(1-exp(-z*v))=log(z)+log(v)+R(z*v),",
        "R(w)=log((1-exp(-w))/w).",
        interval["bell_correction_bound"],
        "```",
        "",
        "One outward-rounded interval evaluation over `0<=z<=10^-30` gives",
        "",
        "```text",
        f"j_0 lower={interval['j0_lower']}",
        f"r_0 lower={interval['r0_lower']}",
        f"s_0 lower={interval['s0_lower']}",
        f"t_0 lower={interval['t0_lower']}",
        f"u_0 lower={interval['u0_lower']}",
        f"t^2*s_1'' ball={interval['scaled_curvature_ball']}",
        f"t^2*s_1'' upper={interval['scaled_curvature_upper']}<200<4000,",
        f"margin below 200={interval['scaled_margin_lower']}.",
        "```",
        "",
        "This is a uniform interval theorem, not evaluation at a representative",
        "large value of `t`. Together with the compact and finite-ray",
        "certificates it covers every real `t>=999` for the first summand.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order8_nested_curvature_compact_certificate.md",
        "outputs/jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate.md",
        "outputs/jensen_window_pf_compound_order8_first_summand_curvature_bridge.md",
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
        "wrote order-eight nested curvature asymptotic-ray certificate: "
        f"{summary['rows']} rows, "
        f"{summary['dimensionless_interval_theorems']} dimensionless interval theorem, "
        f"{summary['asymptotic_curvature_theorems']} asymptotic curvature theorem, "
        f"{summary['open_rows']} open rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
