#!/usr/bin/env python3
"""Certify order-eleven first-summand curvature on the asymptotic saddle ray."""

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
import jensen_window_pf_compound_order9_nested_curvature_asymptotic_ray_certificate as order9  # noqa: E402
import jensen_window_pf_compound_order10_nested_curvature_asymptotic_ray_certificate as order10  # noqa: E402
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_scale,
    series_sub,
)
from jensen_window_pf_compound_order11_high_cumulant_coarse_corridor import (  # noqa: E402
    DEFAULT_OUT as SOURCE_HIGH_CUMULANTS,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_upper_text,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order11_nested_curvature_asymptotic_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order11_nested_curvature_asymptotic_ray_certificate.md"
)
SOURCE_ORDER10_GLOBAL = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_first_summand_curvature_certificate.json"
)
PRECISION_BITS = 384
INVERSE_T_CAP = order9.INVERSE_T_CAP
COLLAR_RADIUS = 8
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
EIGHTH_B_CAP = 2
FIRST_COMPOSITION_CAP = order9.FIRST_COMPOSITION_CAP
NESTED_COMPOSITION_CAP = order9.NESTED_COMPOSITION_CAP
SCALED_TARGET = 2000
ORDER_TEN_CURVATURE_CONSTANT = 4200


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
    global_order10 = load_json(SOURCE_ORDER10_GLOBAL)
    inherited = order10.source_contract()
    if high.get("exact", {}).get("ray_exact_corridor") != (
        "|kappa_r|*q^(r/2-1)/(r-2)!<1, r=19,20, u>=20"
    ):
        raise RuntimeError("order-eleven high-cumulant ray source changed")
    if (
        global_order10.get("status")
        != "rigorous global order-ten first-summand curvature theorem on t>=1251"
        or global_order10.get("theorem")
        != "z_1''(t)<=4200/t^2 for every real t>=1251"
    ):
        raise RuntimeError("global order-ten curvature source changed")
    return {
        "inherited_order10_asymptotic": inherited,
        "new_cumulants": high["exact"]["ray_exact_corridor"],
        "order10_global_curvature": (
            "z_1''(t)<=4200/t^2 for every real t>=1251"
        ),
        "high_cumulant_sha256": sha256(SOURCE_HIGH_CUMULANTS),
        "order10_global_sha256": sha256(SOURCE_ORDER10_GLOBAL),
    }


def raw_high_normalized_cap(order: int, cumulant_cap: int) -> Fraction:
    return order10.raw_high_normalized_cap(order, cumulant_cap)


def exact_ratio_gates() -> dict:
    inherited = order10.exact_ratio_gates()
    z = INVERSE_T_CAP
    eighth_ratio = Fraction(1) / (1 - COLLAR_RADIUS * z) ** 19
    if not eighth_ratio < NORMALIZED_CAP:
        raise RuntimeError("eighth-layer collar ratio gate failed")
    high_caps = {}
    for degree in (19, 20):
        raw = raw_high_normalized_cap(degree, 1)
        shifted = raw * eighth_ratio
        if not shifted < EIGHTH_B_CAP:
            raise RuntimeError(f"normalized derivative cap failed at {degree}")
        high_caps[str(degree)] = {
            "raw_cap": str(raw),
            "shifted_cap": str(shifted),
            "target_cap": EIGHTH_B_CAP,
            "cumulant_cap": 1,
        }
    return {
        "inherited": inherited,
        "eighth_order_upper_ratio": str(eighth_ratio),
        "upper_ratio_cap": str(NORMALIZED_CAP),
        "high_normalized_caps": high_caps,
        "all_strict": True,
    }


def dimensionless_x_floor() -> Fraction:
    z = INVERSE_T_CAP
    floor = (
        Fraction(9) / (2 + 3 * z)
        - ORDER_TEN_CURVATURE_CONSTANT * z / (1 - z**2)
    )
    if floor <= Fraction(449, 100):
        raise RuntimeError("asymptotic dimensionless X floor failed")
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
    ):
        b.append(order6.arb_interval(Fraction(-cap), Fraction(cap)))

    z_cap = flint.arb(INVERSE_T_CAP.numerator) / INVERSE_T_CAP.denominator
    z = z_cap / 2 + flint.arb(0, z_cap / 2)
    ell = order6.stable_normalized_log(b, 18, FIRST_COMPOSITION_CAP)
    J = [
        2 * b[degree]
        - z * (degree + 1) * (degree + 2) * ell[degree + 2]
        for degree in range(17)
    ]
    h = series_add(
        series_scale(ell[:17], 2),
        order6.stable_normalized_log(J, 16, NESTED_COMPOSITION_CAP),
    )
    R = [
        3 * b[degree]
        - z * (degree + 1) * (degree + 2) * h[degree + 2]
        for degree in range(15)
    ]
    q = series_add(
        series_sub(series_scale(h[:15], 2), ell[:15]),
        order6.stable_normalized_log(R, 14, NESTED_COMPOSITION_CAP),
    )
    S = [
        4 * b[degree]
        - z * (degree + 1) * (degree + 2) * q[degree + 2]
        for degree in range(13)
    ]
    p = series_add(
        series_sub(series_scale(q[:13], 2), h[:13]),
        order6.stable_normalized_log(S, 12, NESTED_COMPOSITION_CAP),
    )
    T = [
        5 * b[degree]
        - z * (degree + 1) * (degree + 2) * p[degree + 2]
        for degree in range(11)
    ]
    r = series_add(
        series_sub(series_scale(p[:11], 2), q[:11]),
        order6.stable_normalized_log(T, 10, NESTED_COMPOSITION_CAP),
    )
    U = [
        6 * b[degree]
        - z * (degree + 1) * (degree + 2) * r[degree + 2]
        for degree in range(9)
    ]
    s = series_add(
        series_sub(series_scale(r[:9], 2), p[:9]),
        order6.stable_normalized_log(U, 8, NESTED_COMPOSITION_CAP),
    )
    V = [
        7 * b[degree]
        - z * (degree + 1) * (degree + 2) * s[degree + 2]
        for degree in range(7)
    ]
    w = series_add(
        series_sub(series_scale(s[:7], 2), r[:7]),
        order6.stable_normalized_log(V, 6, NESTED_COMPOSITION_CAP),
    )
    W = [
        8 * b[degree]
        - z * (degree + 1) * (degree + 2) * w[degree + 2]
        for degree in range(5)
    ]
    z_coordinate = series_add(
        series_sub(series_scale(w[:5], 2), s[:5]),
        order6.stable_normalized_log(W, 4, NESTED_COMPOSITION_CAP),
    )
    X = [
        9 * b[degree]
        - z * (degree + 1) * (degree + 2) * z_coordinate[degree + 2]
        for degree in range(3)
    ]
    y = series_add(
        series_sub(series_scale(z_coordinate[:3], 2), w[:3]),
        order6.stable_normalized_log(X, 2, NESTED_COMPOSITION_CAP),
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
    }
    if not bool(all(value > 0 for value in coordinates.values())):
        raise RuntimeError("asymptotic stable-coordinate floor failed")
    floor = dimensionless_x_floor()
    floor_arb = flint.arb(floor.numerator) / floor.denominator
    stable_z_second = 2 * z_coordinate[2]
    stable_w_second = 2 * w[2]
    positive_x_second = flint.arb(max(flint.arb(0), flint.arb(X[2].upper())))
    positive_term = 2 * positive_x_second / floor_arb
    formula_upper = (
        flint.arb((2 * stable_z_second - stable_w_second).upper())
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
        "dimensionless_X_floor": str(floor),
        "stable_z_second_scaled_ball": stable_z_second.str(50).replace("e", "E"),
        "stable_w_second_scaled_ball": stable_w_second.str(50).replace("e", "E"),
        "coordinate_X_second_ball": X[2].str(50).replace("e", "E"),
        "positive_phi_X_second_upper": arb_upper_text(positive_term),
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
            "t^2*y''<=2*t^2*z''-t^2*w''+phi(X)*max(t^2*X'',0)"
        ),
        "singular_term_reduction": (
            "phi(X)<=1/X and t*X>=F imply "
            "phi(X)*max(t^2*X'',0)<=2*max(X_2,0)/F"
        ),
    }


def build_artifact() -> dict:
    ratios = exact_ratio_gates()
    interval = dimensionless_interval()
    rows = [
        RayRow(
            "co11ncarc_01_shifted_H_boxes",
            "exact_theorem_composition",
            "ready_to_apply",
            "Shifted normalized derivative boxes through H20 hold on every t+-8 collar for u>=20.",
            "|t^(r-1)H^(r)/(r-2)!|<C_r, r=2,...,20",
            "First Newman summand on the asymptotic saddle ray only.",
            ratios,
        ),
        RayRow(
            "co11ncarc_02_dimensionless_X_floor",
            "exact_analytic_bound",
            "ready_to_apply",
            "The global order-ten curvature theorem gives a strict dimensionless X floor.",
            "t*X>=9/(2+3/t)-4200/t/(1-1/t^2)>449/100",
            "Uses the certified order-ten first-summand theorem only.",
            {"dimensionless_X_floor": interval["dimensionless_X_floor"]},
        ),
        RayRow(
            "co11ncarc_03_dimensionless_box",
            "interval_analytic_theorem",
            "ready_to_apply",
            "One outward-rounded dimensionless box proves the asymptotic order-eleven ceiling.",
            "t^2*y_1''(t)<2000 for every mode u>=20",
            "Uniform interval theorem over z in (0,10^-30].",
            interval,
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order11_nested_curvature_asymptotic_ray_certificate",
        "date": "2026-07-18",
        "status": "rigorous order-eleven first-summand curvature theorem on u>=20",
        "proof_boundary": (
            "This artifact proves the displayed asymptotic first-summand ray. "
            "It does not prove the finite saddle ray, compact range, full-kernel "
            "transfer, entry, PF-infinity, or RH."
        ),
        "source_contract": source_contract(),
        "ratio_gates": ratios,
        "dimensionless_interval": interval,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": 3,
            "ready_rows": 3,
            "asymptotic_ray_theorems": 1,
            "open_finite_ray_ranges": 1,
            "open_compact_ranges": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order11_nested_curvature_asymptotic_ray_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order11_nested_curvature_asymptotic_ray_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    interval = artifact["dimensionless_interval"]
    lines = [
        "# Jensen-Window PF Compound Order-Eleven Nested Curvature Asymptotic-Ray Certificate",
        "",
        "Date: 2026-07-18",
        "",
        "Status: rigorous first-summand order-eleven theorem on `u>=20`.",
        "This is not a proof of the finite or compact ranges, full-kernel transfer, entry,",
        "PF-infinity, or RH.",
        "",
        "```text",
        "t^2*y_1''(t)<2000<6000 for every mode u>=20",
        f"dimensionless X floor={interval['dimensionless_X_floor']}",
        f"scaled upper={interval['scaled_curvature_upper']}",
        f"scaled margin to 2000={interval['scaled_margin_lower']}",
        "```",
        "",
        "The finite ray `2001/1000<=u<=20` and compact handoff remain open.",
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
        "wrote order-eleven nested asymptotic certificate: scaled upper "
        f"{artifact['dimensionless_interval']['scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
