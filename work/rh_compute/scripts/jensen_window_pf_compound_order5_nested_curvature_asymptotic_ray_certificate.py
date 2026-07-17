#!/usr/bin/env python3
"""Certify the nested first-summand order-five curvature on the u>=20 ray."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
import math
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
    CURVATURE_CONSTANT,
    series_add,
    series_log,
    series_scale,
    series_sub,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_upper_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.md"
)
SOURCE_ORDER4_RAY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.json"
)
SOURCE_FINITE_RAY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.json"
)
PRECISION_BITS = 256
INVERSE_T_CAP = Fraction(1, 10**30)
COLLAR_RADIUS = 3
NORMALIZED_CAP = Fraction(1001, 1000)
B0_FLOOR = Fraction(969, 1000)
B1_MAGNITUDE_FLOOR = Fraction(959, 1000)
STABLE_LOG_ARGUMENT_CAP = Fraction(1, 1000)
FIRST_COMPOSITION_CAP = 2
NESTED_COMPOSITION_CAP = 100
SCALED_TARGET = 10


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
    if finite.get("finite_ray", {}).get("mode_range") != ["2001/1000", "20"]:
        raise RuntimeError("finite-ray source contract changed")
    return {
        "normalized_H_boxes": normalized["consequence"],
        "inverse_t_cap": geometry["inverse_slope_cap"],
        "collar_t20_minus_t19_lower": geometry[
            "collar_t20_minus_t19_lower"
        ],
        "defect_series": order4["logarithmic_defect"]["series"],
        "defect_derivative_majorant": order4["logarithmic_defect"][
            "series_majorant"
        ],
        "finite_ray": "q_1''(t)<=60/t^2 for every mode 2<=u<=20",
    }


def arb_interval(lower: Fraction, upper: Fraction) -> flint.arb:
    low = flint.arb(lower.numerator) / lower.denominator
    high = flint.arb(upper.numerator) / upper.denominator
    return (low + high) / 2 + flint.arb(0, (high - low) / 2)


def exact_ratio_gates() -> dict:
    z = INVERSE_T_CAP
    upper_ratio = Fraction(1, 1) / (1 - COLLAR_RADIUS * z) ** 7
    lower_ratio = Fraction(1, 1) / (1 + COLLAR_RADIUS * z) ** 2
    if not upper_ratio < NORMALIZED_CAP:
        raise RuntimeError("upper collar ratio gate failed")
    if not lower_ratio > Fraction(999, 1000):
        raise RuntimeError("lower collar ratio gate failed")
    if not Fraction(97, 100) * lower_ratio > B0_FLOOR:
        raise RuntimeError("normalized B0 floor failed")
    if not Fraction(24, 25) * lower_ratio > B1_MAGNITUDE_FLOOR:
        raise RuntimeError("normalized B1 floor failed")
    collar = Fraction(
        "684140403277229664153888231548284194636517090552223519946531E-23"
    )
    if collar <= COLLAR_RADIUS:
        raise RuntimeError("mode collar does not cover t+-3")
    return {
        "upper_ratio": str(upper_ratio),
        "upper_ratio_cap": str(NORMALIZED_CAP),
        "lower_ratio": str(lower_ratio),
        "lower_ratio_floor": "999/1000",
        "B0_floor": str(B0_FLOOR),
        "B1_magnitude_floor": str(B1_MAGNITUDE_FLOOR),
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
        scaled = abs(values[index])
        if not bool(scaled < coefficient_cap):
            raise RuntimeError(
                f"normalized coefficient {index} exceeds cap {coefficient_cap}"
            )
    logarithm = series_log(values, order)
    return [
        logarithm[index] + flint.arb(
            0, correction_radius(index, coefficient_cap)
        )
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

    z_cap = flint.arb(INVERSE_T_CAP.numerator) / INVERSE_T_CAP.denominator
    z = z_cap / 2 + flint.arb(0, z_cap / 2)
    ell = stable_normalized_log(b, 6, FIRST_COMPOSITION_CAP)

    j_hat = [
        2 * b[order]
        - z * (order + 1) * (order + 2) * ell[order + 2]
        for order in range(5)
    ]
    h = series_add(
        series_scale(ell[:5], 2),
        stable_normalized_log(j_hat, 4, NESTED_COMPOSITION_CAP),
    )

    r_hat = [
        3 * b[order]
        - z * (order + 1) * (order + 2) * h[order + 2]
        for order in range(3)
    ]
    q = series_add(
        series_sub(series_scale(h[:3], 2), ell[:3]),
        stable_normalized_log(r_hat, 2, NESTED_COMPOSITION_CAP),
    )
    scaled_curvature = 2 * q[2]

    if not bool(j_hat[0] > 0 and r_hat[0] > 0):
        raise RuntimeError("dimensionless stable-coordinate floor failed")
    if not bool(scaled_curvature < SCALED_TARGET):
        raise RuntimeError("dimensionless scaled-curvature target failed")
    if not bool(z_cap * NESTED_COMPOSITION_CAP < flint.arb(1) / 1000):
        raise RuntimeError("stable logarithm argument cap failed")

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
        "scaled_curvature_ball": scaled_curvature.str(60).replace("e", "E"),
        "scaled_curvature_upper": arb_upper_text(scaled_curvature),
        "scaled_target": SCALED_TARGET,
        "scaled_margin_lower": arb_lower_text(
            flint.arb(SCALED_TARGET) - scaled_curvature
        ),
        "j0_lower": arb_lower_text(j_hat[0]),
        "r0_lower": arb_lower_text(r_hat[0]),
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
            id="co5ncarc_01_shifted_H_boxes",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The normalized H-jet boxes remain valid on the complete three-unit collar around every central u>=20 point.",
            formula=(
                "0<x_r<=1 (2<=r<=8), x_2>=97/100, x_3>=24/25; "
                "t/(t-3)<(1001/1000)^(1/7)"
            ),
            proof_boundary="Existing exact ray theorem plus strict rational collar arithmetic.",
        ),
        RayRow(
            id="co5ncarc_02_normalized_B_boxes",
            role="exact_interval_lemma",
            readiness="ready_to_apply",
            claim="Every tent-averaged B derivative lies in one dimensionless alternating box.",
            formula=(
                "b_0 in [969/1000,1001/1000], "
                "b_1 in [-1001/1000,-959/1000], "
                "(-1)^r*b_r in [0,1001/1000] for 2<=r<=6"
            ),
            proof_boundary="Unit-mass tent averaging of the shifted normalized H boxes.",
        ),
        RayRow(
            id="co5ncarc_03_stable_log_majorant",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="The analytic stable-log correction is uniformly negligible in normalized variables.",
            formula=interval["bell_correction_bound"],
            proof_boundary=(
                "Uses |R^(k)(w)|<1, the exact partial-Bell identity, and "
                "0<=z*C<1/1000."
            ),
        ),
        RayRow(
            id="co5ncarc_04_dimensionless_interval",
            role="interval_analytic_theorem",
            readiness="ready_to_apply",
            claim="A single outward-rounded dimensionless jet box proves a ten-inverse-square curvature ceiling on the whole asymptotic ray.",
            formula="t^2*q_1''(t)<10 for every mode u>=20",
            proof_boundary="Uniform interval theorem over z in [0,10^-30], not a finite sample.",
            diagnostics={
                "scaled_upper": interval["scaled_curvature_upper"],
                "scaled_margin_lower": interval["scaled_margin_lower"],
                "j0_lower": interval["j0_lower"],
                "r0_lower": interval["r0_lower"],
            },
        ),
        RayRow(
            id="co5ncarc_05_asymptotic_consequence",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The asymptotic theorem is stronger than the continuous target needed by the order-five bridge.",
            formula=(
                "t^2*q_1''(t)<10<60 for u>=20, hence "
                "q_1''(t)<=60/t^2 on the asymptotic ray"
            ),
            proof_boundary="First Newman summand at lambda=-100 only.",
        ),
        RayRow(
            id="co5ncarc_06_global_handoff",
            role="exact_handoff",
            readiness="ready_to_apply",
            claim="The compact, finite-ray, and asymptotic certificates now cover every real t>=320.",
            formula=(
                "[320<=t<=V'(2)] union [2<=u<=20] union [u>=20] "
                "covers t>=320"
            ),
            proof_boundary="The global composition is recorded in a separate entry theorem.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate",
        "date": "2026-07-13",
        "status": (
            "rigorous nested first-summand order-five curvature theorem on "
            "the complete u>=20 ray"
        ),
        "proof_boundary": (
            "This artifact proves q_1''(t)<=60/t^2 on u>=20. It does not by "
            "itself compose the compact range, full-kernel transfer, order-five "
            "entry, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.py"
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "inverse_t_cap": str(INVERSE_T_CAP),
            "collar_radius": COLLAR_RADIUS,
            "normalized_cap": str(NORMALIZED_CAP),
            "B0_floor": str(B0_FLOOR),
            "B1_magnitude_floor": str(B1_MAGNITUDE_FLOOR),
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
            "normalized_H_box_theorems": 1,
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
        "# Jensen-Window PF Compound Order-Five Nested Curvature Asymptotic-Ray Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous nested first-summand order-five curvature theorem on",
        "the complete `u>=20` ray. This is not a proof of full order-five",
        "entry by itself, PF-infinity, RH, or `Lambda <= 0`.",
        "This is not by itself a proof of the full compact-to-kernel composition.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.py",
        "```",
        "",
        "## Shifted Normalized H Boxes",
        "",
        "The completed order-four ray theorem gives, throughout the enlarged",
        "`u>=19` collar,",
        "",
        "```text",
        "x_r=(-1)^r*t^(r-1)*H^(r)/(r-2)!",
        "0<x_r<=1 for 2<=r<=8, x_2>=97/100, x_3>=24/25,",
        "1/t<=10^-30.",
        "```",
        "",
        "Since `t(20)-t(19)>6.84e36`, every point within `t+-3` of a",
        "central `u>=20` point remains in that collar. Exact rational arithmetic",
        "also gives",
        "",
        "```text",
        f"(1-3*10^-30)^(-7)={ratios['upper_ratio']}<1001/1000,",
        f"(1+3*10^-30)^(-2)={ratios['lower_ratio']}>999/1000.",
        "```",
        "",
        "For `b_r=t^(r+1)B^(r)/r!`, unit-mass tent averaging therefore yields",
        "",
        "```text",
        "b_0 in [969/1000,1001/1000],",
        "b_1 in [-1001/1000,-959/1000],",
        "(-1)^r*b_r in [0,1001/1000], 2<=r<=6.",
        "```",
        "",
        "## Analytic Stable Logs",
        "",
        "Put `z=1/t` and scale the local variable by `y=s/t`. Then the Taylor",
        "jet of `B` is `z*b`. For",
        "",
        "```text",
        "R_0(w)=log((1-exp(-w))/w),",
        "```",
        "",
        "the existing convergent product series proves `|R_0^(k)(w)|<1`",
        "through order six whenever `0<=w<=1/1000`. If the normalized Taylor",
        "coefficients of `v` are bounded by `C`, the exact partial-Bell identity",
        "gives the coefficient correction",
        "",
        "```text",
        interval["bell_correction_bound"],
        "```",
        "",
        "The proof uses `C=2` at the defect layer and `C=100` at each nested",
        "stable layer. Arb verifies both coefficient caps and keeps `z*C<1/1000`.",
        "The constants `log(z)` disappear from every derivative used below.",
        "",
        "## Dimensionless Interval",
        "",
        "The normalized recursion is",
        "",
        "```text",
        interval["normalized_J_definition"],
        interval["normalized_R_definition"],
        interval["normalized_q_identity"],
        "```",
        "",
        "A single outward-rounded interval evaluation over the complete box",
        "`0<=z<=10^-30` gives",
        "",
        "```text",
        f"j_0 lower={interval['j0_lower']}",
        f"r_0 lower={interval['r0_lower']}",
        f"t^2*q_1'' ball={interval['scaled_curvature_ball']}",
        f"t^2*q_1'' upper={interval['scaled_curvature_upper']}<10,",
        f"margin below 10={interval['scaled_margin_lower']}.",
        "```",
        "",
        "This is a uniform interval theorem, not evaluation at a representative",
        "large value of `t`. Consequently",
        "",
        "```text",
        "t^2*q_1''(t)<10<60 for every mode u>=20.",
        "```",
        "",
        "Together with the compact and finite-ray certificates, this covers all",
        "real `t>=320`. The full composition is recorded separately.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md",
        "outputs/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.md",
        "outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-five nested curvature asymptotic-ray certificate: "
        f"{summary['rows']} rows, "
        f"{summary['dimensionless_interval_theorems']} dimensionless interval theorem, "
        f"{summary['asymptotic_curvature_theorems']} asymptotic curvature theorem, "
        f"{summary['open_rows']} open rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
