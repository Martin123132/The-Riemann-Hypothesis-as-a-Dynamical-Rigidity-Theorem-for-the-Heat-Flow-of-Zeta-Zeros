#!/usr/bin/env python3
"""Certify order-ten first-summand curvature on the asymptotic saddle ray."""

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
import jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate as order9  # noqa: E402
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_scale,
    series_sub,
)
from jensen_window_pf_compound_order10_high_cumulant_coarse_corridor import (  # noqa: E402
    DEFAULT_OUT as SOURCE_HIGH_CUMULANTS,
)
from jensen_window_pf_compound_order10_nested_curvature_finite_ray_certificate import (  # noqa: E402
    DEFAULT_OUT as SOURCE_FINITE_RAY,
    sha256,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_upper_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_nested_curvature_asymptotic_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order10_nested_curvature_asymptotic_ray_certificate.md"
)
PRECISION_BITS = 384
INVERSE_T_CAP = order9.INVERSE_T_CAP
COLLAR_RADIUS = order9.COLLAR_RADIUS
COLLAR_MODE_START = order9.COLLAR_MODE_START
NORMALIZED_CAP = order9.NORMALIZED_CAP
B0_FLOOR = order9.B0_FLOOR
B1_MAGNITUDE_FLOOR = order9.B1_MAGNITUDE_FLOOR
MID_CUMULANT_CAP = order9.MID_CUMULANT_CAP
TOP_CUMULANT_CAP = order9.TOP_CUMULANT_CAP
ULTRA_CUMULANT_CAP = order9.ULTRA_CUMULANT_CAP
MID_B_CAP = order9.MID_B_CAP
TOP_B_CAP = order9.TOP_B_CAP
ULTRA_B_CAP = order9.ULTRA_B_CAP
NEW_B_CAP = order9.NEW_B_CAP
SEVENTH_B_CAP = 2
FIRST_COMPOSITION_CAP = order9.FIRST_COMPOSITION_CAP
NESTED_COMPOSITION_CAP = order9.NESTED_COMPOSITION_CAP
SCALED_TARGET = 1000
ORDER_NINE_CURVATURE_CONSTANT = 4200


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
    high = load_json(SOURCE_HIGH_CUMULANTS)
    finite = load_json(SOURCE_FINITE_RAY)
    inherited = order9.source_contract()
    if high.get("exact", {}).get("ray_exact_corridor") != (
        "|kappa_r|*q^(r/2-1)/(r-2)!<1, r=17,18, u>=20"
    ):
        raise RuntimeError("order-ten high-cumulant ray source changed")
    if finite.get("finite_ray", {}).get("mode_range") != ["2001/1000", "20"]:
        raise RuntimeError("order-ten finite-ray source changed")
    return {
        "inherited_order9": inherited,
        "new_cumulants": high["exact"]["ray_exact_corridor"],
        "finite_ray": finite["finite_ray"]["theorem"],
        "order9_global_curvature": (
            "w_1''(t)<=4200/t^2 for every real t>=5700"
        ),
        "high_cumulant_sha256": sha256(SOURCE_HIGH_CUMULANTS),
        "finite_ray_sha256": sha256(SOURCE_FINITE_RAY),
    }


def raw_high_normalized_cap(order: int, cumulant_cap: int) -> Fraction:
    return order9.raw_high_normalized_cap(order, cumulant_cap)


def exact_ratio_gates() -> dict:
    inherited = order9.exact_ratio_gates()
    z = INVERSE_T_CAP
    seventh_ratio = Fraction(1) / (1 - COLLAR_RADIUS * z) ** 17
    if not seventh_ratio < NORMALIZED_CAP:
        raise RuntimeError("seventh-layer collar ratio gate failed")
    high_caps = {}
    for degree in (17, 18):
        raw = raw_high_normalized_cap(degree, 1)
        shifted = raw * seventh_ratio
        if not shifted < SEVENTH_B_CAP:
            raise RuntimeError(f"normalized derivative cap failed at {degree}")
        high_caps[str(degree)] = {
            "raw_cap": str(raw),
            "shifted_cap": str(shifted),
            "target_cap": SEVENTH_B_CAP,
            "cumulant_cap": 1,
        }
    return {
        "inherited": inherited,
        "seventh_order_upper_ratio": str(seventh_ratio),
        "upper_ratio_cap": str(NORMALIZED_CAP),
        "high_normalized_caps": high_caps,
        "all_strict": True,
    }


def dimensionless_w_floor() -> Fraction:
    z = INVERSE_T_CAP
    floor = (
        Fraction(8) / (2 + 3 * z)
        - ORDER_NINE_CURVATURE_CONSTANT * z / (1 - z**2)
    )
    if floor <= Fraction(399, 100):
        raise RuntimeError("asymptotic dimensionless W floor failed")
    return floor


def dimensionless_interval() -> dict:
    flint.ctx.prec = PRECISION_BITS
    zero = Fraction(0)
    b = [
        order6.arb_interval(B0_FLOOR, NORMALIZED_CAP),
        order6.arb_interval(-NORMALIZED_CAP, -B1_MAGNITUDE_FLOOR),
    ]
    for degree in range(2, 7):
        b.append(
            order6.arb_interval(zero, NORMALIZED_CAP)
            if degree % 2 == 0
            else order6.arb_interval(-NORMALIZED_CAP, zero)
        )
    for cap in (
        MID_B_CAP,
        MID_B_CAP,
        TOP_B_CAP,
        TOP_B_CAP,
        ULTRA_B_CAP,
        ULTRA_B_CAP,
        NEW_B_CAP,
        NEW_B_CAP,
        SEVENTH_B_CAP,
        SEVENTH_B_CAP,
    ):
        b.append(order6.arb_interval(Fraction(-cap), Fraction(cap)))

    z_cap = flint.arb(INVERSE_T_CAP.numerator) / INVERSE_T_CAP.denominator
    z = z_cap / 2 + flint.arb(0, z_cap / 2)
    ell = order6.stable_normalized_log(b, 16, FIRST_COMPOSITION_CAP)
    J = [
        2 * b[degree]
        - z * (degree + 1) * (degree + 2) * ell[degree + 2]
        for degree in range(15)
    ]
    h = series_add(
        series_scale(ell[:15], 2),
        order6.stable_normalized_log(J, 14, NESTED_COMPOSITION_CAP),
    )
    R = [
        3 * b[degree]
        - z * (degree + 1) * (degree + 2) * h[degree + 2]
        for degree in range(13)
    ]
    q = series_add(
        series_sub(series_scale(h[:13], 2), ell[:13]),
        order6.stable_normalized_log(R, 12, NESTED_COMPOSITION_CAP),
    )
    S = [
        4 * b[degree]
        - z * (degree + 1) * (degree + 2) * q[degree + 2]
        for degree in range(11)
    ]
    p = series_add(
        series_sub(series_scale(q[:11], 2), h[:11]),
        order6.stable_normalized_log(S, 10, NESTED_COMPOSITION_CAP),
    )
    T = [
        5 * b[degree]
        - z * (degree + 1) * (degree + 2) * p[degree + 2]
        for degree in range(9)
    ]
    r = series_add(
        series_sub(series_scale(p[:9], 2), q[:9]),
        order6.stable_normalized_log(T, 8, NESTED_COMPOSITION_CAP),
    )
    U = [
        6 * b[degree]
        - z * (degree + 1) * (degree + 2) * r[degree + 2]
        for degree in range(7)
    ]
    s = series_add(
        series_sub(series_scale(r[:7], 2), p[:7]),
        order6.stable_normalized_log(U, 6, NESTED_COMPOSITION_CAP),
    )
    V = [
        7 * b[degree]
        - z * (degree + 1) * (degree + 2) * s[degree + 2]
        for degree in range(5)
    ]
    w = series_add(
        series_sub(series_scale(s[:5], 2), r[:5]),
        order6.stable_normalized_log(V, 4, NESTED_COMPOSITION_CAP),
    )
    W = [
        8 * b[degree]
        - z * (degree + 1) * (degree + 2) * w[degree + 2]
        for degree in range(3)
    ]

    coordinates = {
        "J": J[0],
        "R": R[0],
        "S": S[0],
        "T": T[0],
        "U": U[0],
        "V": V[0],
        "W": W[0],
    }
    if not bool(all(value > 0 for value in coordinates.values())):
        raise RuntimeError("asymptotic stable-coordinate floor failed")
    floor = dimensionless_w_floor()
    floor_arb = flint.arb(floor.numerator) / floor.denominator
    stable_w_second = 2 * w[2]
    stable_s_second = 2 * s[2]
    positive_w_second = flint.arb(max(flint.arb(0), flint.arb(W[2].upper())))
    positive_term = 2 * positive_w_second / floor_arb
    formula_upper = (
        flint.arb((2 * stable_w_second - stable_s_second).upper())
        + positive_term
    )
    if not bool(formula_upper < SCALED_TARGET):
        raise RuntimeError(f"asymptotic scaled target failed: {formula_upper}")
    if not bool(
        z_cap * FIRST_COMPOSITION_CAP < flint.arb(1) / 1000
        and z_cap * NESTED_COMPOSITION_CAP < flint.arb(1) / 1000
    ):
        raise RuntimeError("stable-log argument cap failed")
    return {
        "scaled_variable": "z=1/t in (0,10^-30]",
        "normalized_B_boxes": [value.str(40).replace("e", "E") for value in b],
        "dimensionless_W_floor": str(floor),
        "stable_w_second_scaled_ball": stable_w_second.str(50).replace("e", "E"),
        "stable_s_second_scaled_ball": stable_s_second.str(50).replace("e", "E"),
        "coordinate_W_second_ball": W[2].str(50).replace("e", "E"),
        "positive_phi_W_second_upper": arb_upper_text(positive_term),
        "scaled_curvature_upper": arb_upper_text(formula_upper),
        "scaled_target": SCALED_TARGET,
        "scaled_margin_lower": arb_lower_text(
            flint.arb(SCALED_TARGET) - formula_upper
        ),
        **{
            f"{name}0_lower": arb_lower_text(value)
            for name, value in coordinates.items()
        },
        "first_composition_coefficient_cap": FIRST_COMPOSITION_CAP,
        "nested_composition_coefficient_cap": NESTED_COMPOSITION_CAP,
        "one_sided_formula": (
            "t^2*z''<=2*t^2*w''-t^2*s''+phi(W)*max(t^2*W'',0)"
        ),
        "singular_term_reduction": (
            "phi(W)<=1/W and t*W>=F imply "
            "phi(W)*max(t^2*W'',0)<=2*max(W_2,0)/F"
        ),
    }


def build_artifact() -> dict:
    ratios = exact_ratio_gates()
    interval = dimensionless_interval()
    rows = [
        RayRow(
            "co10ncarc_01_shifted_H_boxes",
            "exact_theorem_composition",
            "ready_to_apply",
            "Shifted normalized derivative boxes through H18 hold on every t+-7 collar for u>=20.",
            "|t^(r-1)H^(r)/(r-2)!|<C_r, r=2,...,18",
            "First Newman summand on the asymptotic saddle ray only.",
            ratios,
        ),
        RayRow(
            "co10ncarc_02_dimensionless_W_floor",
            "exact_analytic_bound",
            "ready_to_apply",
            "The global order-nine curvature theorem gives a strict dimensionless W floor.",
            "t*W>=8/(2+3/t)-4200/t/(1-1/t^2)>399/100",
            "Uses the certified order-nine first-summand theorem only.",
            {"dimensionless_W_floor": interval["dimensionless_W_floor"]},
        ),
        RayRow(
            "co10ncarc_03_dimensionless_box",
            "interval_analytic_theorem",
            "ready_to_apply",
            "One outward-rounded dimensionless box proves the asymptotic order-ten ceiling.",
            "t^2*z_1''(t)<1000 for every mode u>=20",
            "Uniform interval theorem over z in (0,10^-30].",
            interval,
        ),
        RayRow(
            "co10ncarc_04_range_composition",
            "exact_theorem_composition",
            "ready_to_apply",
            "The finite and asymptotic certificates cover every saddle mode u>=2001/1000.",
            "[2001/1000<=u<=20] union [u>=20]",
            "The compact range through mode 2001/1000 remains separate.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order10_nested_curvature_asymptotic_ray_certificate",
        "date": "2026-07-16",
        "status": "rigorous order-ten first-summand curvature theorem on u>=20",
        "proof_boundary": (
            "This artifact proves the displayed asymptotic first-summand ray and "
            "composes the ranges on u>=2001/1000. It does not prove the compact "
            "range, full-kernel transfer, entry, PF-infinity, or RH."
        ),
        "source_contract": source_contract(),
        "ratio_gates": ratios,
        "dimensionless_interval": interval,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": 4,
            "ready_rows": 4,
            "asymptotic_ray_theorems": 1,
            "global_above_mode_2001_1000_compositions": 1,
            "open_compact_ranges": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order10_nested_curvature_asymptotic_ray_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order10_nested_curvature_asymptotic_ray_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    interval = artifact["dimensionless_interval"]
    lines = [
        "# Jensen-Window PF Compound Order-Ten Nested Curvature Asymptotic-Ray Certificate",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous first-summand order-ten theorem on `u>=20`.",
        "This is not a proof of the compact range, full-kernel transfer, entry,",
        "PF-infinity, or RH.",
        "",
        "```text",
        "t^2*z_1''(t)<1000<5500 for every mode u>=20",
        f"dimensionless W floor={interval['dimensionless_W_floor']}",
        f"scaled upper={interval['scaled_curvature_upper']}",
        f"scaled margin to 1000={interval['scaled_margin_lower']}",
        "```",
        "",
        "Together with the finite-ray certificate, this covers every",
        "`u>=2001/1000`. The compact handoff remains open.",
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
    print(
        "wrote order-ten nested asymptotic certificate: scaled upper "
        f"{artifact['dimensionless_interval']['scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
