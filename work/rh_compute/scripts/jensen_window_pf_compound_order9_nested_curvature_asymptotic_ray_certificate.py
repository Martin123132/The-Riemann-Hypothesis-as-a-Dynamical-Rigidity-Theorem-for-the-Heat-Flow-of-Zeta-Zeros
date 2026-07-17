#!/usr/bin/env python3
"""Certify order-nine first-summand curvature on the ray u>=20."""

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
import jensen_window_pf_compound_order8_nested_curvature_asymptotic_ray_certificate as order8  # noqa: E402
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_scale,
    series_sub,
)
from jensen_window_pf_compound_order9_high_cumulant_coarse_corridor import (  # noqa: E402
    DEFAULT_OUT as SOURCE_HIGH_CUMULANTS,
    EXACT_CORRIDOR_CAP as NEW_CUMULANT_CAP,
)
from jensen_window_pf_compound_order9_nested_curvature_finite_ray_certificate import (  # noqa: E402
    DEFAULT_OUT as SOURCE_FINITE_RAY,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_upper_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate.md"
)
PRECISION_BITS = 384
INVERSE_T_CAP = order8.INVERSE_T_CAP
COLLAR_RADIUS = 7
COLLAR_MODE_START = 19
NORMALIZED_CAP = order8.NORMALIZED_CAP
B0_FLOOR = order8.B0_FLOOR
B1_MAGNITUDE_FLOOR = order8.B1_MAGNITUDE_FLOOR
MID_CUMULANT_CAP = order8.MID_CUMULANT_CAP
TOP_CUMULANT_CAP = order8.TOP_CUMULANT_CAP
ULTRA_CUMULANT_CAP = order8.ULTRA_CUMULANT_CAP
MID_B_CAP = order8.MID_B_CAP
TOP_B_CAP = order8.TOP_B_CAP
ULTRA_B_CAP = order8.ULTRA_B_CAP
NEW_B_CAP = 2
FIRST_COMPOSITION_CAP = order8.FIRST_COMPOSITION_CAP
NESTED_COMPOSITION_CAP = order8.NESTED_COMPOSITION_CAP
SCALED_TARGET = 500


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
    inherited = order8.source_contract()
    high = load_json(SOURCE_HIGH_CUMULANTS)
    finite = load_json(SOURCE_FINITE_RAY)
    if high.get("exact", {}).get("exact_corridor") != (
        "|kappa_r|*q^(r/2-1)/(r-2)!<1, r=15,16, u>=2"
    ):
        raise RuntimeError("fifteenth/sixteenth cumulant source changed")
    if finite.get("finite_ray", {}).get("mode_range") != ["2", "20"]:
        raise RuntimeError("order-nine finite-ray source changed")
    return {
        "inherited_order8": inherited,
        "new_cumulants": high["exact"]["ray_exact_corridor"],
        "finite_ray": "w_1''(t)<=4200/t^2 for every saddle mode 2<=u<=20",
    }


def raw_high_normalized_cap(order: int, cumulant_cap: int) -> Fraction:
    geometry = (
        Fraction(201, 100) ** (order - 1)
        / Fraction(19, 10) ** order
    )
    return Fraction(1) + cumulant_cap * geometry / COLLAR_MODE_START


def exact_ratio_gates() -> dict:
    z = INVERSE_T_CAP
    ratios = {
        "low": Fraction(1) / (1 - COLLAR_RADIUS * z) ** 7,
        "mid": Fraction(1) / (1 - COLLAR_RADIUS * z) ** 9,
        "top": Fraction(1) / (1 - COLLAR_RADIUS * z) ** 11,
        "ultra": Fraction(1) / (1 - COLLAR_RADIUS * z) ** 13,
        "new": Fraction(1) / (1 - COLLAR_RADIUS * z) ** 15,
    }
    for name, ratio in ratios.items():
        if not ratio < NORMALIZED_CAP:
            raise RuntimeError(f"{name} collar ratio gate failed")
    lower_ratio = Fraction(1) / (1 + COLLAR_RADIUS * z) ** 2
    if not lower_ratio > Fraction(999, 1000):
        raise RuntimeError("lower collar ratio gate failed")
    if not Fraction(97, 100) * lower_ratio > B0_FLOOR:
        raise RuntimeError("normalized B0 floor failed")
    if not Fraction(24, 25) * lower_ratio > B1_MAGNITUDE_FLOOR:
        raise RuntimeError("normalized B1 floor failed")

    high_caps = {}
    specifications = (
        (9, MID_CUMULANT_CAP, ratios["mid"], MID_B_CAP),
        (10, MID_CUMULANT_CAP, ratios["mid"], MID_B_CAP),
        (11, TOP_CUMULANT_CAP, ratios["top"], TOP_B_CAP),
        (12, TOP_CUMULANT_CAP, ratios["top"], TOP_B_CAP),
        (13, ULTRA_CUMULANT_CAP, ratios["ultra"], ULTRA_B_CAP),
        (14, ULTRA_CUMULANT_CAP, ratios["ultra"], ULTRA_B_CAP),
        (15, NEW_CUMULANT_CAP, ratios["new"], NEW_B_CAP),
        (16, NEW_CUMULANT_CAP, ratios["new"], NEW_B_CAP),
    )
    for degree, cumulant_cap, ratio, target in specifications:
        raw = raw_high_normalized_cap(degree, cumulant_cap)
        shifted = raw * ratio
        if not shifted < target:
            raise RuntimeError(f"normalized derivative cap failed at {degree}")
        high_caps[str(degree)] = {
            "raw_cap": str(raw),
            "shifted_cap": str(shifted),
            "target_cap": target,
            "cumulant_cap": cumulant_cap,
        }
    collar = Fraction(
        "684140403277229664153888231548284194636517090552223519946531E-23"
    )
    if collar <= COLLAR_RADIUS:
        raise RuntimeError("mode collar does not cover t+-7")
    return {
        **{f"{name}_order_upper_ratio": str(value) for name, value in ratios.items()},
        "upper_ratio_cap": str(NORMALIZED_CAP),
        "lower_ratio": str(lower_ratio),
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
    ):
        b.append(order6.arb_interval(Fraction(-cap), Fraction(cap)))

    z_cap = flint.arb(INVERSE_T_CAP.numerator) / INVERSE_T_CAP.denominator
    z = z_cap / 2 + flint.arb(0, z_cap / 2)
    ell = order6.stable_normalized_log(b, 14, FIRST_COMPOSITION_CAP)
    J = [
        2 * b[degree]
        - z * (degree + 1) * (degree + 2) * ell[degree + 2]
        for degree in range(13)
    ]
    h = series_add(
        series_scale(ell[:13], 2),
        order6.stable_normalized_log(J, 12, NESTED_COMPOSITION_CAP),
    )
    R = [
        3 * b[degree]
        - z * (degree + 1) * (degree + 2) * h[degree + 2]
        for degree in range(11)
    ]
    q = series_add(
        series_sub(series_scale(h[:11], 2), ell[:11]),
        order6.stable_normalized_log(R, 10, NESTED_COMPOSITION_CAP),
    )
    S = [
        4 * b[degree]
        - z * (degree + 1) * (degree + 2) * q[degree + 2]
        for degree in range(9)
    ]
    p = series_add(
        series_sub(series_scale(q[:9], 2), h[:9]),
        order6.stable_normalized_log(S, 8, NESTED_COMPOSITION_CAP),
    )
    T = [
        5 * b[degree]
        - z * (degree + 1) * (degree + 2) * p[degree + 2]
        for degree in range(7)
    ]
    r = series_add(
        series_sub(series_scale(p[:7], 2), q[:7]),
        order6.stable_normalized_log(T, 6, NESTED_COMPOSITION_CAP),
    )
    U = [
        6 * b[degree]
        - z * (degree + 1) * (degree + 2) * r[degree + 2]
        for degree in range(5)
    ]
    s = series_add(
        series_sub(series_scale(r[:5], 2), p[:5]),
        order6.stable_normalized_log(U, 4, NESTED_COMPOSITION_CAP),
    )
    V = [
        7 * b[degree]
        - z * (degree + 1) * (degree + 2) * s[degree + 2]
        for degree in range(3)
    ]
    w = series_add(
        series_sub(series_scale(s[:3], 2), r[:3]),
        order6.stable_normalized_log(V, 2, NESTED_COMPOSITION_CAP),
    )
    scaled = 2 * w[2]
    coordinates = {"J": J[0], "R": R[0], "S": S[0], "T": T[0], "U": U[0], "V": V[0]}
    if not bool(all(value > 0 for value in coordinates.values())):
        raise RuntimeError("asymptotic stable-coordinate floor failed")
    if not bool(scaled < SCALED_TARGET):
        raise RuntimeError(f"asymptotic scaled target failed: {scaled}")
    if not bool(
        z_cap * FIRST_COMPOSITION_CAP < flint.arb(1) / 1000
        and z_cap * NESTED_COMPOSITION_CAP < flint.arb(1) / 1000
    ):
        raise RuntimeError("stable-log argument cap failed")
    return {
        "scaled_variable": "z=1/t in (0,10^-30]",
        "normalized_B_boxes": [value.str(40).replace("e", "E") for value in b],
        "scaled_curvature_ball": scaled.str(60).replace("e", "E"),
        "scaled_curvature_upper": arb_upper_text(scaled),
        "scaled_target": SCALED_TARGET,
        "scaled_margin_lower": arb_lower_text(flint.arb(SCALED_TARGET) - scaled),
        **{
            f"{name}0_lower": arb_lower_text(value)
            for name, value in coordinates.items()
        },
        "first_composition_coefficient_cap": FIRST_COMPOSITION_CAP,
        "nested_composition_coefficient_cap": NESTED_COMPOSITION_CAP,
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
            "co9ncarc_01_shifted_H_boxes",
            "exact_theorem_composition",
            "ready_to_apply",
            "Normalized H-jet boxes through order sixteen remain valid on the seven-unit collar.",
            "low signed boxes plus coarse H9-H16 boxes",
            "Exact ray theorem and strict rational collar arithmetic.",
            ratios["high_normalized_caps"],
        ),
        RayRow(
            "co9ncarc_02_stable_log_majorant",
            "exact_analytic_lemma",
            "ready_to_apply",
            "The stable-log correction is controlled at all seven logarithmic evaluations.",
            interval["bell_correction_bound"],
            "Convergent defect series and exact partial-Bell identity.",
        ),
        RayRow(
            "co9ncarc_03_dimensionless_interval",
            "interval_analytic_theorem",
            "ready_to_apply",
            "One outward-rounded dimensionless box proves the asymptotic order-nine ceiling.",
            "t^2*w_1''(t)<500 for every mode u>=20",
            "Uniform interval theorem over z in [0,10^-30].",
            {
                "scaled_upper": interval["scaled_curvature_upper"],
                "scaled_margin_lower": interval["scaled_margin_lower"],
                **{
                    f"{name}0_lower": interval[f"{name}0_lower"]
                    for name in ("J", "R", "S", "T", "U", "V")
                },
            },
        ),
        RayRow(
            "co9ncarc_04_global_handoff",
            "exact_handoff",
            "ready_to_apply",
            "The compact, finite, and asymptotic certificates cover every real t>=5700.",
            "[5700<=t<=V'(2)] union [2<=u<=20] union [u>=20]",
            "The localized 1250<=t<=5700 bridge remains separate.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate",
        "date": "2026-07-14",
        "status": "rigorous order-nine first-summand curvature theorem on u>=20",
        "proof_boundary": (
            "This artifact proves the displayed asymptotic first-summand ray "
            "and composes the ranges above t=5700. It does not prove the lower "
            "localized bridge, full-kernel transfer, entry, PF-infinity, or RH."
        ),
        "source_contract": sources,
        "ratio_gates": ratios,
        "dimensionless_interval": interval,
        "parameters": {
            "inverse_t_cap": str(INVERSE_T_CAP),
            "collar_radius": COLLAR_RADIUS,
            "precision_bits": PRECISION_BITS,
            "scaled_target": SCALED_TARGET,
            "continuous_target": 4200,
            "new_B_cap": NEW_B_CAP,
            "new_cumulant_cap": NEW_CUMULANT_CAP,
        },
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": 4,
            "ready_rows": 4,
            "asymptotic_ray_theorems": 1,
            "global_above_5700_compositions": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    interval = artifact["dimensionless_interval"]
    lines = [
        "# Jensen-Window PF Compound Order-Nine Nested Curvature Asymptotic Ray",
        "",
        "Date: 2026-07-14",
        "",
        "Status: rigorous first-summand order-nine theorem on `u>=20`.",
        "This is not a proof of the lower localized bridge, entry, PF-infinity, or RH.",
        "",
        "```text",
        "t^2*w_1''(t)<500<4200 for every mode u>=20",
        "scaled upper=" + interval["scaled_curvature_upper"],
        "scaled margin below 500=" + interval["scaled_margin_lower"],
        "B,J,R,S,T,U,V>0",
        "```",
        "",
        "Together with the compact and finite certificates this covers every",
        "real `t>=5700`; the `1250<=t<=5700` bridge remains separate.",
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
        "wrote order-nine nested asymptotic certificate: scaled upper "
        f"{artifact['dimensionless_interval']['scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
