#!/usr/bin/env python3
"""Certify order-twelve first-summand curvature on the asymptotic saddle ray."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
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
import jensen_window_pf_compound_order11_nested_curvature_asymptotic_ray_certificate as order11  # noqa: E402
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_scale,
    series_sub,
)
from jensen_window_pf_compound_order12_high_cumulant_coarse_corridor import (  # noqa: E402
    DEFAULT_OUT as SOURCE_HIGH_CUMULANTS,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_upper_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order12_nested_curvature_asymptotic_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order12_nested_curvature_asymptotic_ray_certificate.md"
)
SOURCE_ORDER11_GLOBAL = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_first_summand_curvature_certificate.json"
)
PRECISION_BITS = 384
INVERSE_T_CAP = order11.INVERSE_T_CAP
COLLAR_RADIUS = order11.COLLAR_RADIUS
NORMALIZED_CAP = order11.NORMALIZED_CAP
B0_FLOOR = order11.B0_FLOOR
B1_MAGNITUDE_FLOOR = order11.B1_MAGNITUDE_FLOOR
MID_B_CAP = order11.MID_B_CAP
TOP_B_CAP = order11.TOP_B_CAP
ULTRA_B_CAP = order11.ULTRA_B_CAP
NEW_B_CAP = order11.NEW_B_CAP
SEVENTH_B_CAP = order11.SEVENTH_B_CAP
EIGHTH_B_CAP = order11.EIGHTH_B_CAP
NINTH_B_CAP = 2
FIRST_COMPOSITION_CAP = order11.FIRST_COMPOSITION_CAP
NESTED_COMPOSITION_CAP = order11.NESTED_COMPOSITION_CAP
SCALED_TARGET = 8000
ORDER_ELEVEN_CURVATURE_CONSTANT = 6000


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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def source_contract() -> dict:
    high = load_json(SOURCE_HIGH_CUMULANTS)
    global_order11 = load_json(SOURCE_ORDER11_GLOBAL)
    inherited = order11.source_contract()
    if high.get("exact", {}).get("ray_exact_corridor") != (
        "|kappa_r|*q^(r/2-1)/(r-2)!<1, r=21,22, u>=20"
    ):
        raise RuntimeError("order-twelve high-cumulant ray source changed")
    if (
        global_order11.get("status")
        != "rigorous global order-eleven first-summand curvature theorem on t>=1252"
        or global_order11.get("theorem")
        != "y_1''(t)<=6000/t^2 for every real t>=1252"
    ):
        raise RuntimeError("global order-eleven curvature source changed")
    return {
        "inherited_order11_asymptotic": inherited,
        "new_cumulants": high["exact"]["ray_exact_corridor"],
        "order11_global_curvature": global_order11["theorem"],
        "high_cumulant_sha256": sha256(SOURCE_HIGH_CUMULANTS),
        "order11_global_sha256": sha256(SOURCE_ORDER11_GLOBAL),
    }


def raw_high_normalized_cap(order: int, cumulant_cap: int) -> Fraction:
    return order11.raw_high_normalized_cap(order, cumulant_cap)


def exact_ratio_gates() -> dict:
    inherited = order11.exact_ratio_gates()
    z = INVERSE_T_CAP
    ninth_ratio = Fraction(1) / (1 - COLLAR_RADIUS * z) ** 21
    if not ninth_ratio < NORMALIZED_CAP:
        raise RuntimeError("ninth-layer collar ratio gate failed")
    high_caps = {}
    for degree in (21, 22):
        raw = raw_high_normalized_cap(degree, 1)
        shifted = raw * ninth_ratio
        if not shifted < NINTH_B_CAP:
            raise RuntimeError(f"normalized derivative cap failed at {degree}")
        high_caps[str(degree)] = {
            "raw_cap": str(raw),
            "shifted_cap": str(shifted),
            "target_cap": NINTH_B_CAP,
            "cumulant_cap": 1,
        }
    return {
        "inherited": inherited,
        "ninth_order_upper_ratio": str(ninth_ratio),
        "upper_ratio_cap": str(NORMALIZED_CAP),
        "high_normalized_caps": high_caps,
        "all_strict": True,
    }


def dimensionless_d9_floor() -> Fraction:
    z = INVERSE_T_CAP
    floor = (
        Fraction(10) / (2 + 3 * z)
        - ORDER_ELEVEN_CURVATURE_CONSTANT * z / (1 - z**2)
    )
    if floor <= Fraction(499, 100):
        raise RuntimeError("asymptotic dimensionless D9 floor failed")
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
        EIGHTH_B_CAP,
        EIGHTH_B_CAP,
        NINTH_B_CAP,
        NINTH_B_CAP,
    ):
        b.append(order6.arb_interval(Fraction(-cap), Fraction(cap)))

    z_cap = flint.arb(INVERSE_T_CAP.numerator) / INVERSE_T_CAP.denominator
    z = z_cap / 2 + flint.arb(0, z_cap / 2)
    ell = order6.stable_normalized_log(b, 20, FIRST_COMPOSITION_CAP)
    J = [
        2 * b[degree] - z * (degree + 1) * (degree + 2) * ell[degree + 2]
        for degree in range(19)
    ]
    h = series_add(
        series_scale(ell[:19], 2),
        order6.stable_normalized_log(J, 18, NESTED_COMPOSITION_CAP),
    )
    R = [
        3 * b[degree] - z * (degree + 1) * (degree + 2) * h[degree + 2]
        for degree in range(17)
    ]
    q = series_add(
        series_sub(series_scale(h[:17], 2), ell[:17]),
        order6.stable_normalized_log(R, 16, NESTED_COMPOSITION_CAP),
    )
    S = [
        4 * b[degree] - z * (degree + 1) * (degree + 2) * q[degree + 2]
        for degree in range(15)
    ]
    p = series_add(
        series_sub(series_scale(q[:15], 2), h[:15]),
        order6.stable_normalized_log(S, 14, NESTED_COMPOSITION_CAP),
    )
    T = [
        5 * b[degree] - z * (degree + 1) * (degree + 2) * p[degree + 2]
        for degree in range(13)
    ]
    r = series_add(
        series_sub(series_scale(p[:13], 2), q[:13]),
        order6.stable_normalized_log(T, 12, NESTED_COMPOSITION_CAP),
    )
    U = [
        6 * b[degree] - z * (degree + 1) * (degree + 2) * r[degree + 2]
        for degree in range(11)
    ]
    s = series_add(
        series_sub(series_scale(r[:11], 2), p[:11]),
        order6.stable_normalized_log(U, 10, NESTED_COMPOSITION_CAP),
    )
    V = [
        7 * b[degree] - z * (degree + 1) * (degree + 2) * s[degree + 2]
        for degree in range(9)
    ]
    w = series_add(
        series_sub(series_scale(s[:9], 2), r[:9]),
        order6.stable_normalized_log(V, 8, NESTED_COMPOSITION_CAP),
    )
    W = [
        8 * b[degree] - z * (degree + 1) * (degree + 2) * w[degree + 2]
        for degree in range(7)
    ]
    z_coordinate = series_add(
        series_sub(series_scale(w[:7], 2), s[:7]),
        order6.stable_normalized_log(W, 6, NESTED_COMPOSITION_CAP),
    )
    X = [
        9 * b[degree]
        - z * (degree + 1) * (degree + 2) * z_coordinate[degree + 2]
        for degree in range(5)
    ]
    y = series_add(
        series_sub(series_scale(z_coordinate[:5], 2), w[:5]),
        order6.stable_normalized_log(X, 4, NESTED_COMPOSITION_CAP),
    )
    D9 = [
        10 * b[degree] - z * (degree + 1) * (degree + 2) * y[degree + 2]
        for degree in range(3)
    ]
    v = series_add(
        series_sub(series_scale(y[:3], 2), z_coordinate[:3]),
        order6.stable_normalized_log(D9, 2, NESTED_COMPOSITION_CAP),
    )

    coordinates = {
        "J": J[0],
        "R": R[0],
        "S": S[0],
        "T": T[0],
        "U": U[0],
        "V": V[0],
        "W": W[0],
        "X": X[0],
        "D9": D9[0],
    }
    if not bool(all(value > 0 for value in coordinates.values())):
        raise RuntimeError("asymptotic stable-coordinate floor failed")
    floor = dimensionless_d9_floor()
    floor_arb = flint.arb(floor.numerator) / floor.denominator
    stable_y_second = 2 * y[2]
    stable_z_second = 2 * z_coordinate[2]
    positive_d9_second = flint.arb(max(flint.arb(0), flint.arb(D9[2].upper())))
    positive_term = 2 * positive_d9_second / floor_arb
    formula_upper = (
        flint.arb((2 * stable_y_second - stable_z_second).upper())
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
        "dimensionless_D9_floor": str(floor),
        "stable_y_second_scaled_ball": stable_y_second.str(50).replace("e", "E"),
        "stable_z_second_scaled_ball": stable_z_second.str(50).replace("e", "E"),
        "coordinate_D9_second_ball": D9[2].str(50).replace("e", "E"),
        "positive_phi_D9_second_upper": arb_upper_text(positive_term),
        "scaled_curvature_upper": arb_upper_text(formula_upper),
        "scaled_composed_ball": (2 * v[2]).str(50).replace("e", "E"),
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
            "t^2*v''<=2*t^2*y''-t^2*z''+phi(D9)*max(t^2*D9'',0)"
        ),
        "singular_term_reduction": (
            "phi(D9)<=1/D9 and t*D9>=F imply "
            "phi(D9)*max(t^2*D9'',0)<=2*max(D9_2,0)/F"
        ),
    }


def build_artifact() -> dict:
    ratios = exact_ratio_gates()
    interval = dimensionless_interval()
    rows = [
        RayRow(
            "co12ncarc_01_shifted_H_boxes",
            "exact_theorem_composition",
            "ready_to_apply",
            "Shifted normalized derivative boxes through H22 hold on every t+-8 collar for u>=20.",
            "|t^(r-1)H^(r)/(r-2)!|<C_r, r=2,...,22",
            "First Newman summand on the asymptotic saddle ray only.",
            ratios,
        ),
        RayRow(
            "co12ncarc_02_dimensionless_D9_floor",
            "exact_analytic_bound",
            "ready_to_apply",
            "The global order-eleven curvature theorem gives a strict dimensionless D9 floor.",
            "t*D9>=10/(2+3/t)-6000/t/(1-1/t^2)>499/100",
            "Uses the certified order-eleven first-summand theorem only.",
            {"dimensionless_D9_floor": interval["dimensionless_D9_floor"]},
        ),
        RayRow(
            "co12ncarc_03_dimensionless_box",
            "interval_analytic_theorem",
            "ready_to_apply",
            "One outward-rounded dimensionless box proves the asymptotic order-twelve ceiling.",
            "t^2*v_1''(t)<8000 for every mode u>=20",
            "Uniform interval theorem over z in (0,10^-30].",
            interval,
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order12_nested_curvature_asymptotic_ray_certificate",
        "date": "2026-07-22",
        "status": "rigorous order-twelve first-summand curvature theorem on u>=20",
        "proof_boundary": (
            "This artifact proves the displayed asymptotic first-summand ray. "
            "It does not prove the finite saddle ray, lower/compact ranges, "
            "global composition, full-kernel transfer, entry, PF-infinity, or RH."
        ),
        "source_contract": source_contract(),
        "ratio_gates": ratios,
        "dimensionless_interval": interval,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": 3,
            "ready_rows": 3,
            "asymptotic_ray_theorems": 1,
            "open_finite_ray_ranges": 0,
            "open_lower_compact_ranges": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order12_nested_curvature_asymptotic_ray_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order12_nested_curvature_asymptotic_ray_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    interval = artifact["dimensionless_interval"]
    lines = [
        "# Order-Twelve Nested Curvature Asymptotic-Ray Certificate",
        "",
        "Date: 2026-07-22",
        "",
        "Status: rigorous first-summand order-twelve theorem on `u>=20`.",
        "This is not a proof of order twelve, PF-infinity, RH, or `Lambda<=0`.",
        "",
        "```text",
        "t^2*v_1''(t)<8000 for every mode u>=20",
        "dimensionless D9 floor=" + interval["dimensionless_D9_floor"],
        "scaled upper=" + interval["scaled_curvature_upper"],
        "scaled margin=" + interval["scaled_margin_lower"],
        "```",
        "",
        "The finite ray `2001/1000<=u<=20` is certified separately. The lower",
        "and compact handoff remain open.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "wrote order-twelve asymptotic ray: scaled upper "
        f"{artifact['dimensionless_interval']['scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
