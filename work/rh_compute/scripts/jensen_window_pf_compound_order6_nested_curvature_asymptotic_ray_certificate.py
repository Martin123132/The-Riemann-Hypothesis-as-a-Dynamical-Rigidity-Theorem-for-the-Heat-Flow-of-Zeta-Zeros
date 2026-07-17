#!/usr/bin/env python3
"""Certify the nested first-summand order-six curvature on the u>=20 ray."""

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

from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_log,
    series_scale,
    series_sub,
)
from jensen_window_pf_compound_order6_nested_curvature_interval_core import (  # noqa: E402
    CURVATURE_CONSTANT,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_upper_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.md"
)
SOURCE_ORDER4_RAY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.json"
)
SOURCE_HIGH_CUMULANTS = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.json"
)
SOURCE_FINITE_RAY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.json"
)
PRECISION_BITS = 384
INVERSE_T_CAP = Fraction(1, 10**30)
COLLAR_RADIUS = 4
COLLAR_MODE_START = 19
NORMALIZED_CAP = Fraction(1001, 1000)
B0_FLOOR = Fraction(969, 1000)
B1_MAGNITUDE_FLOOR = Fraction(959, 1000)
HIGH_CUMULANT_CAP = 50000
HIGH_B_CAP = 2500
FIRST_COMPOSITION_CAP = 2501
NESTED_COMPOSITION_CAP = 100
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
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_contract() -> dict:
    order4 = load_json(SOURCE_ORDER4_RAY)
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
    if high.get("exact", {}).get("exact_corridor_cap") != HIGH_CUMULANT_CAP:
        raise RuntimeError("high-cumulant corridor source changed")
    if finite.get("finite_ray", {}).get("mode_range") != ["2", "20"]:
        raise RuntimeError("finite-ray source contract changed")
    return {
        "normalized_H_boxes": normalized["consequence"],
        "inverse_t_cap": geometry["inverse_slope_cap"],
        "collar_t20_minus_t19_lower": geometry[
            "collar_t20_minus_t19_lower"
        ],
        "high_cumulant_corridor": high["exact"]["exact_corridor"],
        "finite_ray": "p_1''(t)<=200/t^2 for every saddle mode 2<=u<=20",
    }


def arb_interval(lower: Fraction, upper: Fraction) -> flint.arb:
    low = flint.arb(lower.numerator) / lower.denominator
    high = flint.arb(upper.numerator) / upper.denominator
    return (low + high) / 2 + flint.arb(0, (high - low) / 2)


def raw_high_normalized_cap(order: int) -> Fraction:
    geometry = (
        Fraction(201, 100) ** (order - 1)
        / Fraction(19, 10) ** order
    )
    return Fraction(1) + HIGH_CUMULANT_CAP * geometry / COLLAR_MODE_START


def exact_ratio_gates() -> dict:
    z = INVERSE_T_CAP
    low_order_upper_ratio = Fraction(1) / (1 - COLLAR_RADIUS * z) ** 7
    high_order_upper_ratio = Fraction(1) / (1 - COLLAR_RADIUS * z) ** 9
    lower_ratio = Fraction(1) / (1 + COLLAR_RADIUS * z) ** 2
    if not low_order_upper_ratio < NORMALIZED_CAP:
        raise RuntimeError("low-order collar ratio gate failed")
    if not high_order_upper_ratio < NORMALIZED_CAP:
        raise RuntimeError("high-order collar ratio gate failed")
    if not lower_ratio > Fraction(999, 1000):
        raise RuntimeError("lower collar ratio gate failed")
    if not Fraction(97, 100) * lower_ratio > B0_FLOOR:
        raise RuntimeError("normalized B0 floor failed")
    if not Fraction(24, 25) * lower_ratio > B1_MAGNITUDE_FLOOR:
        raise RuntimeError("normalized B1 floor failed")

    high_caps = {}
    for order in (9, 10):
        raw = raw_high_normalized_cap(order)
        shifted = raw * high_order_upper_ratio
        if not shifted < HIGH_B_CAP:
            raise RuntimeError(f"high normalized derivative cap failed at {order}")
        high_caps[str(order)] = {
            "raw_cap": str(raw),
            "shifted_cap": str(shifted),
            "target_cap": HIGH_B_CAP,
        }

    collar = Fraction(
        "684140403277229664153888231548284194636517090552223519946531E-23"
    )
    if collar <= COLLAR_RADIUS:
        raise RuntimeError("mode collar does not cover t+-4")
    return {
        "low_order_upper_ratio": str(low_order_upper_ratio),
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


def correction_radius(order: int, coefficient_cap: int) -> flint.arb:
    factor = 1 if order == 0 else 2 ** (order - 1)
    value = Fraction(factor * coefficient_cap) * INVERSE_T_CAP
    return flint.arb(value.numerator) / value.denominator


def stable_normalized_log(
    values: list[flint.arb], order: int, coefficient_cap: int
) -> list[flint.arb]:
    """Enclose log(v)+R(z*v), omitting the irrelevant constant log(z)."""
    if len(values) <= order:
        raise ValueError("normalized stable series is too short")
    for index in range(order + 1):
        if not bool(abs(values[index]) < coefficient_cap):
            raise RuntimeError(
                f"normalized coefficient {index} exceeds cap {coefficient_cap}"
            )
    logarithm = series_log(values, order)
    return [
        logarithm[index]
        + flint.arb(0, correction_radius(index, coefficient_cap))
        for index in range(order + 1)
    ]


def dimensionless_interval() -> dict:
    flint.ctx.prec = PRECISION_BITS
    zero = Fraction(0)
    b = [
        arb_interval(B0_FLOOR, NORMALIZED_CAP),
        arb_interval(-NORMALIZED_CAP, -B1_MAGNITUDE_FLOOR),
    ]
    for order in range(2, 7):
        b.append(
            arb_interval(zero, NORMALIZED_CAP)
            if order % 2 == 0
            else arb_interval(-NORMALIZED_CAP, zero)
        )
    b.extend(
        [
            arb_interval(Fraction(-HIGH_B_CAP), Fraction(HIGH_B_CAP)),
            arb_interval(Fraction(-HIGH_B_CAP), Fraction(HIGH_B_CAP)),
        ]
    )

    z_cap = flint.arb(INVERSE_T_CAP.numerator) / INVERSE_T_CAP.denominator
    z = z_cap / 2 + flint.arb(0, z_cap / 2)
    ell = stable_normalized_log(b, 8, FIRST_COMPOSITION_CAP)

    j_hat = [
        2 * b[order]
        - z * (order + 1) * (order + 2) * ell[order + 2]
        for order in range(7)
    ]
    h = series_add(
        series_scale(ell[:7], 2),
        stable_normalized_log(j_hat, 6, NESTED_COMPOSITION_CAP),
    )

    r_hat = [
        3 * b[order]
        - z * (order + 1) * (order + 2) * h[order + 2]
        for order in range(5)
    ]
    q = series_add(
        series_sub(series_scale(h[:5], 2), ell[:5]),
        stable_normalized_log(r_hat, 4, NESTED_COMPOSITION_CAP),
    )

    s_hat = [
        4 * b[order]
        - z * (order + 1) * (order + 2) * q[order + 2]
        for order in range(3)
    ]
    p = series_add(
        series_sub(series_scale(q[:3], 2), h[:3]),
        stable_normalized_log(s_hat, 2, NESTED_COMPOSITION_CAP),
    )
    scaled_curvature = 2 * p[2]

    if not bool(j_hat[0] > 0 and r_hat[0] > 0 and s_hat[0] > 0):
        raise RuntimeError("dimensionless stable-coordinate floor failed")
    if not bool(scaled_curvature < SCALED_TARGET):
        raise RuntimeError("dimensionless scaled-curvature target failed")
    if not bool(z_cap * FIRST_COMPOSITION_CAP < flint.arb(1) / 1000):
        raise RuntimeError("first stable-log argument cap failed")
    if not bool(z_cap * NESTED_COMPOSITION_CAP < flint.arb(1) / 1000):
        raise RuntimeError("nested stable-log argument cap failed")

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
        "normalized_J_definition": (
            "J_y,r=z*j_r, j_r=2*b_r-z*(r+1)*(r+2)*ell_(r+2)"
        ),
        "normalized_J_boxes": coefficient_rows(j_hat),
        "normalized_h_identity": "h=2*ell+log(z)+log(j)+R(z*j)",
        "normalized_R_definition": (
            "R_y,r=z*r_r, r_r=3*b_r-z*(r+1)*(r+2)*h_(r+2)"
        ),
        "normalized_R_boxes": coefficient_rows(r_hat),
        "normalized_q_identity": "q=2*h-ell+log(z)+log(r)+R(z*r)",
        "normalized_S_definition": (
            "S_y,r=z*s_r, s_r=4*b_r-z*(r+1)*(r+2)*q_(r+2)"
        ),
        "normalized_S_boxes": coefficient_rows(s_hat),
        "normalized_p_identity": "p=2*q-h+log(z)+log(s)+R(z*s)",
        "scaled_curvature_ball": scaled_curvature.str(60).replace("e", "E"),
        "scaled_curvature_upper": arb_upper_text(scaled_curvature),
        "scaled_target": SCALED_TARGET,
        "scaled_margin_lower": arb_lower_text(
            flint.arb(SCALED_TARGET) - scaled_curvature
        ),
        "j0_lower": arb_lower_text(j_hat[0]),
        "r0_lower": arb_lower_text(r_hat[0]),
        "s0_lower": arb_lower_text(s_hat[0]),
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
            "co6ncarc_01_shifted_low_H_boxes",
            "exact_theorem_composition",
            "ready_to_apply",
            "The signed normalized H-jet boxes through order eight remain valid on the full four-unit collar.",
            "0<x_r<=1 (2<=r<=8), x_2>=97/100, x_3>=24/25",
            "Existing exact ray theorem plus strict rational collar arithmetic.",
        ),
        RayRow(
            "co6ncarc_02_high_H_boxes",
            "exact_analytic_consequence",
            "ready_to_apply",
            "The coarse exact cumulant corridors give absolute normalized boxes for H^(9) and H^(10).",
            "|b_7|<2500 and |b_8|<2500 on the four-unit collar",
            "Hurwitz midpoint ceiling, exact ray geometry, and the 50000 cumulant cap.",
            ratios["high_normalized_caps"],
        ),
        RayRow(
            "co6ncarc_03_normalized_B_boxes",
            "exact_interval_lemma",
            "ready_to_apply",
            "Every tent-averaged B derivative lies in one dimensionless box through order eight.",
            "alternating unit boxes through b_6; |b_7|,|b_8|<2500",
            "Unit-mass tent averaging of shifted normalized H boxes.",
        ),
        RayRow(
            "co6ncarc_04_stable_log_majorant",
            "exact_analytic_lemma",
            "ready_to_apply",
            "The analytic stable-log correction is uniformly negligible at all four layers.",
            interval["bell_correction_bound"],
            "Convergent defect series and the exact partial-Bell identity.",
        ),
        RayRow(
            "co6ncarc_05_dimensionless_interval",
            "interval_analytic_theorem",
            "ready_to_apply",
            "One outward-rounded dimensionless jet box proves the asymptotic order-six curvature ceiling.",
            "t^2*p_1''(t)<100 for every mode u>=20",
            "Uniform interval theorem over z in [0,10^-30], not point sampling.",
            {
                "scaled_upper": interval["scaled_curvature_upper"],
                "scaled_margin_lower": interval["scaled_margin_lower"],
                "j0_lower": interval["j0_lower"],
                "r0_lower": interval["r0_lower"],
                "s0_lower": interval["s0_lower"],
            },
        ),
        RayRow(
            "co6ncarc_06_asymptotic_consequence",
            "exact_theorem_composition",
            "ready_to_apply",
            "The asymptotic theorem is stronger than the continuous target needed by the order-six bridge.",
            "t^2*p_1''(t)<100<200 for u>=20",
            "First Newman summand at lambda=-100 only.",
        ),
        RayRow(
            "co6ncarc_07_global_handoff",
            "exact_handoff",
            "ready_to_apply",
            "The compact, finite-ray, and asymptotic certificates cover every real t>=321.",
            "[321<=t<=V'(2)] union [2<=u<=20] union [u>=20] covers t>=321",
            "The global composition is recorded separately.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate",
        "date": "2026-07-13",
        "status": "rigorous nested first-summand order-six curvature theorem on the complete u>=20 ray",
        "proof_boundary": (
            "This artifact proves p_1''(t)<=200/t^2 on u>=20. It does not "
            "by itself compose the compact range, full-kernel transfer, order-six "
            "entry, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order6_high_cumulant_coarse_corridor.md",
            "outputs/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order6_first_summand_curvature_bridge.md",
            "outputs/formal_core.md",
        ],
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.py",
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "inverse_t_cap": str(INVERSE_T_CAP),
            "collar_radius": COLLAR_RADIUS,
            "normalized_cap": str(NORMALIZED_CAP),
            "B0_floor": str(B0_FLOOR),
            "B1_magnitude_floor": str(B1_MAGNITUDE_FLOOR),
            "high_cumulant_cap": HIGH_CUMULANT_CAP,
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
        "# Jensen-Window PF Compound Order-Six Nested Curvature Asymptotic-Ray Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous nested first-summand order-six curvature theorem on",
        "the complete `u>=20` ray. This certificate is not a proof of full",
        "order-six entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.py",
        "```",
        "",
        "## Shifted H Boxes",
        "",
        "The completed order-four ray theorem gives signed unit boxes through",
        "`H^(8)`. The exact order-six coarse corridor adds",
        "",
        "```text",
        "|kappa_r|*q^(r/2-1)/(r-2)!<50000, r=9,10.",
        "```",
        "",
        "The proved geometry",
        "",
        "```text",
        "t^(r-1)q^(1-r/2)/V''^(r/2)<=D_r/u,",
        "D_r=(201/100)^(r-1)/(19/10)^r, u>=19,",
        "```",
        "",
        "and the midpoint Hurwitz ceiling give, throughout the `t+-4` collar,",
        "",
        "```text",
        f"raw normalized H9 cap={ratios['high_normalized_caps']['9']['raw_cap']},",
        f"raw normalized H10 cap={ratios['high_normalized_caps']['10']['raw_cap']},",
        "|b_7|<2500, |b_8|<2500.",
        "```",
        "",
        "The low derivatives retain",
        "",
        "```text",
        "b_0 in [969/1000,1001/1000],",
        "b_1 in [-1001/1000,-959/1000],",
        "(-1)^r*b_r in [0,1001/1000], 2<=r<=6.",
        "```",
        "",
        "## Dimensionless Stable Recursion",
        "",
        "Put `z=1/t` and use the scaled local variable `y=s/t`. Each stable",
        "logarithm is evaluated as",
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
        "`0<=z<=10^-30` then gives",
        "",
        "```text",
        f"j_0 lower={interval['j0_lower']}",
        f"r_0 lower={interval['r0_lower']}",
        f"s_0 lower={interval['s0_lower']}",
        f"t^2*p_1'' ball={interval['scaled_curvature_ball']}",
        f"t^2*p_1'' upper={interval['scaled_curvature_upper']}<100<200,",
        f"margin below 100={interval['scaled_margin_lower']}.",
        "```",
        "",
        "This is a uniform interval theorem, not evaluation at a representative",
        "large value of `t`. Together with the compact and finite-ray",
        "certificates it covers every real `t>=321`.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.md",
        "outputs/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.md",
        "outputs/jensen_window_pf_compound_order6_first_summand_curvature_bridge.md",
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
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-six nested curvature asymptotic-ray certificate: "
        f"{summary['rows']} rows, "
        f"{summary['dimensionless_interval_theorems']} dimensionless interval theorem, "
        f"{summary['asymptotic_curvature_theorems']} asymptotic curvature theorem, "
        f"{summary['open_rows']} open rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
