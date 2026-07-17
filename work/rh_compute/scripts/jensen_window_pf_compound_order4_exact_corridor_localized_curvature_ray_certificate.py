#!/usr/bin/env python3
"""Prove the exact-corridor localized curvature bound on the ray u>=20."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
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
import sympy as sp  # noqa: E402

from jensen_window_pf_compound_order4_gaussian_cumulant_ray_target import (  # noqa: E402
    CORRIDOR_CAPS,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_rational,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.md"
)
SOURCE_CORRIDOR = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.json"
)
SOURCE_FINITE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate.json"
)
SOURCE_BRIDGE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_first_summand_curvature_bridge.json"
)

PRECISION_BITS = 256
COLLAR_MODE_START = 19
CENTRAL_MODE_START = 20
Q_FLOOR = 10**33
W_CAP = Fraction(1, 10**30)
X2_FLOOR = Fraction(97, 100)
X3_FLOOR = Fraction(24, 25)
ELL_SCALED_CAP = 30_000
ENVELOPE_SCALED_CAP = Fraction(1, 1_000)


@dataclass(frozen=True)
class RayRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def sympy_fraction(value: sp.Expr) -> Fraction:
    value = sp.cancel(value)
    numerator, denominator = sp.fraction(value)
    if not numerator.is_Integer or not denominator.is_Integer:
        raise RuntimeError(f"expression is not rational: {value}")
    return Fraction(int(numerator), int(denominator))


def shifted_positive_gate(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    substitutions: dict[sp.Symbol, sp.Expr],
) -> dict:
    numerator, denominator = sp.fraction(sp.factor(expression))
    shifted = sp.Poly(sp.expand(numerator.subs(substitutions)), *variables)
    coefficients = shifted.coeffs()
    if not coefficients or any(coefficient <= 0 for coefficient in coefficients):
        raise RuntimeError("shifted coefficient-positivity gate failed")
    return {
        "numerator_terms": len(shifted.terms()),
        "minimum_numerator_coefficient": str(min(coefficients)),
        "denominator": sp.sstr(sp.factor(denominator)),
    }


def potential_geometry() -> dict:
    """Certify the coarse ray geometry used by every normalized jet bound."""
    u, q, v, capital_q = sp.symbols("u q v Q", positive=True)
    t_numerator = (
        8 * q**2 * u
        + 400 * q * u**2
        - 30 * q * u
        - 2 * q
        - 600 * u**2
        + 15 * u
        + 3
    )
    t = t_numerator / (2 * (2 * q - 3))
    curvature_numerator = u * (
        64 * q**3 * u
        + 16 * q**3
        + 1408 * q**2 * u
        - 84 * q**2
        - 4560 * q * u
        + 120 * q
        + 3600 * u
        - 45
    )
    curvature = curvature_numerator / (4 * (2 * q - 3) ** 2)
    simple_t = (
        2 * u * q
        + 100 * u**2
        - sp.Rational(5, 2) * u
        - 4 * u * q / (2 * q - 3)
        - sp.Rational(1, 2)
    )
    if sp.cancel(t - simple_t) != 0:
        raise RuntimeError("potential slope decomposition failed")

    substitutions = {u: v + COLLAR_MODE_START, q: capital_q + Q_FLOOR}
    lower_t_gate = shifted_positive_gate(
        t - u * q, (v, capital_q), substitutions
    )
    lower_curvature_gate = shifted_positive_gate(
        curvature - sp.Rational(361, 100) * u**2 * q,
        (v, capital_q),
        substitutions,
    )

    flint.ctx.prec = PRECISION_BITS
    u19 = arb_rational(Fraction(COLLAR_MODE_START))
    u20 = arb_rational(Fraction(CENTRAL_MODE_START))
    q19 = flint.arb.pi() * (4 * u19).exp()
    t19 = potential_jet_arb(u19, 1)[1]
    t20 = potential_jet_arb(u20, 1)[1]
    if not bool(q19 > Q_FLOOR):
        raise RuntimeError("q(19) floor failed")
    if not bool(t20 - t19 > 2):
        raise RuntimeError("shifted t collar does not remain above u=19")
    if Q_FLOOR <= 10_000 * COLLAR_MODE_START:
        raise RuntimeError("q/u endpoint dominance failed")

    return {
        "mode_interval": "u>=19 for shifted jets; central u>=20",
        "q_endpoint_lower": arb_lower_text(q19),
        "q_floor": str(Q_FLOOR),
        "q_over_u_monotonicity": "(q/u)'=q*(4*u-1)/u^2>0",
        "q_dominance": "q>=10000*u",
        "t_decomposition": sp.sstr(simple_t),
        "t_bounds": "u*q<=t<=(201/100)*u*q",
        "curvature_bound": "V''(x_u)>=(361/100)*u^2*q",
        "slope_lower_gate": lower_t_gate,
        "curvature_lower_gate": lower_curvature_gate,
        "inverse_slope_cap": fraction_text(W_CAP),
        "collar_t20_minus_t19_lower": arb_lower_text(t20 - t19),
        "collar_consequence": (
            "If central u>=20 and |s|<=1, the mode belonging to t+s is >19."
        ),
    }


def normalized_h_boxes() -> dict:
    """Derive dimensionless H-jet boxes from the exact cumulant corridors."""
    rows: dict[str, dict] = {}
    for order in range(2, 9):
        corridor_cap = (
            Fraction(1001, 1000)
            if order == 2
            else Fraction(CORRIDOR_CAPS[order])
        )
        geometry_factor = (
            Fraction(201, 100) ** (order - 1)
            / Fraction(19, 10) ** order
        )
        endpoint_lower = (
            Fraction(1)
            - Fraction(order - 1, 2 * COLLAR_MODE_START * Q_FLOOR)
            - corridor_cap * geometry_factor / COLLAR_MODE_START
        )
        target = (
            X2_FLOOR
            if order == 2
            else X3_FLOOR if order == 3 else Fraction(0)
        )
        if endpoint_lower <= target:
            raise RuntimeError(f"normalized H floor failed at order {order}")
        rows[str(order)] = {
            "definition": (
                f"x_{order}=(-1)^{order}*t^{order - 1}*"
                f"H^({order})/{math.factorial(order - 2)}"
            ),
            "corridor_cap": fraction_text(corridor_cap),
            "geometry_factor_Dr": fraction_text(geometry_factor),
            "endpoint_lower": fraction_text(endpoint_lower),
            "required_floor": fraction_text(target),
            "floor_margin": fraction_text(endpoint_lower - target),
            "upper": "1",
        }
    return {
        "hurwitz_midpoint_bound": (
            "(t/(t+1/2))^(r-1)<=G_r=(r-1)t^(r-1)"
            "*zeta(r,t+1/2)<=1"
        ),
        "hurwitz_tangent_bound": "G_r>=1-(r-1)/(2t)",
        "cumulant_geometry": (
            "t^(r-1)q^(1-r/2)/V''^(r/2)<=D_r/u, "
            "D_r=(201/100)^(r-1)/(19/10)^r"
        ),
        "kappa2_cap": "kappa_2<=1001/1000",
        "boxes": rows,
        "consequence": (
            "0<x_r<=1 for 2<=r<=8, x_2>=97/100, x_3>=24/25"
        ),
    }


def logarithmic_defect_bounds() -> dict:
    """Bound ell=log(B)+R(B) and its derivatives through order six."""
    k = sp.symbols("k", integer=True, positive=True)
    r_bounds: dict[str, str] = {}
    for order in range(1, 7):
        start = (order + 1) // 2
        bound = (
            (sp.Rational(1, 2) if order == 1 else 0)
            + 2 ** (order + 1)
            * sp.summation(k ** (order - 1) / 36**k, (k, start, sp.oo))
        )
        bound_fraction = sympy_fraction(bound)
        if bound_fraction >= 1:
            raise RuntimeError(f"R derivative bound failed at order {order}")
        r_bounds[str(order)] = fraction_text(bound_fraction)

    log_caps: dict[str, dict] = {}
    for order in range(2, 7):
        log_cap = (
            Fraction(math.factorial(order) * 2 ** (order - 1))
            * Fraction(100, 97) ** order
        )
        composite_cap = math.factorial(order) * 2 ** (order - 1)
        ell_cap = log_cap + composite_cap * W_CAP
        if ell_cap >= ELL_SCALED_CAP:
            raise RuntimeError(f"scaled ell cap failed at order {order}")
        log_caps[str(order)] = {
            "scaled_log_B_cap": fraction_text(log_cap),
            "scaled_R_composite_cap": f"{composite_cap}/t",
            "scaled_ell_cap": str(ELL_SCALED_CAP),
        }

    sharp_log_second = Fraction(200, 97) - Fraction(576, 625)
    sharp_ell_second = sharp_log_second + 4 * W_CAP
    if sharp_ell_second >= Fraction(23, 20):
        raise RuntimeError("sharp ell second-derivative cap failed")

    collar_ratio = Fraction(1001, 1000)
    envelope_bound = (
        Fraction(ELL_SCALED_CAP, 12)
        * W_CAP**3
        * collar_ratio**6
    )
    if envelope_bound >= ENVELOPE_SCALED_CAP:
        raise RuntimeError("shifted derivative envelope cap failed")
    return {
        "defect": "R(B)=log((1-exp(-B))/B)",
        "series": (
            "R(B)=-B/2+sum_{k>=1}(-1)^(k+1)zeta(2k)"
            "*B^(2k)/(k*(2*pi)^(2k))"
        ),
        "series_majorant": (
            "zeta(2k)<2 and 2*pi>6 imply |R^(m)(B)|<1, "
            "1<=m<=6, 0<B<=1"
        ),
        "R_derivative_bounds": r_bounds,
        "bell_identity": (
            "B_{n,k}(1!,2!,...)=(n!/k!)*binomial(n-1,k-1)"
        ),
        "scaled_derivative_caps": log_caps,
        "sharp_log_second_upper": fraction_text(sharp_log_second),
        "sharp_ell_second_upper": "t^2*ell''<=23/20",
        "collar_ratio": "t/(t-1)<=1001/1000",
        "scaled_envelope_bound": fraction_text(envelope_bound),
        "scaled_envelope_cap": fraction_text(ENVELOPE_SCALED_CAP),
        "envelope_scaling": "t^(r+1)*E_r<1/1000, r=0,1,2",
    }


def localized_scaled_gates() -> dict:
    """Check the five rational gates entering the localized U expression."""
    a_minus_lower = (
        2 * X2_FLOOR
        - W_CAP * Fraction(23, 20)
        - ENVELOPE_SCALED_CAP
    )
    a_plus_upper = (
        Fraction(2)
        + W_CAP * ELL_SCALED_CAP
        + ENVELOPE_SCALED_CAP
    )
    c_plus_upper = (
        Fraction(4)
        + W_CAP * ELL_SCALED_CAP
        + ENVELOPE_SCALED_CAP
    )
    p_minus_lower = (
        2 * X3_FLOOR
        - W_CAP * ELL_SCALED_CAP
        - ENVELOPE_SCALED_CAP
    )
    targets = {
        "A_minus": (a_minus_lower, Fraction(193, 100), "lower"),
        "A_plus": (a_plus_upper, Fraction(201, 100), "upper"),
        "C_plus": (c_plus_upper, Fraction(401, 100), "upper"),
        "P_minus": (p_minus_lower, Fraction(191, 100), "lower"),
    }
    for name, (value, target, direction) in targets.items():
        if direction == "lower" and value <= target:
            raise RuntimeError(f"{name} lower gate failed")
        if direction == "upper" and value >= target:
            raise RuntimeError(f"{name} upper gate failed")
    z_upper = Fraction(201, 100) * W_CAP
    if z_upper >= Fraction(1, 1000):
        raise RuntimeError("localized z cap failed")
    return {
        "definitions": {
            "A_minus": "t*(j_0-E_0)",
            "A_plus": "t*(j_0+E_0)",
            "C_plus": "t^3*max(j_2+E_2,0)",
            "P_minus": "t^2*max(abs(j_1)-E_1,0)",
        },
        "gates": {
            "t2_ell_second_upper": "23/20",
            "A_minus_lower": "193/100",
            "A_plus_upper": "201/100",
            "C_plus_upper": "401/100",
            "P_minus_lower": "191/100",
            "j_plus_upper": fraction_text(z_upper),
            "j_plus_cap": "1/1000",
        },
        "exact_precomparison_values": {
            name: fraction_text(value)
            for name, (value, _target, _direction) in targets.items()
        },
        "positive_localized_gap": True,
    }


def scalar_curvature_comparison() -> dict:
    ell_term = 2 * Fraction(23, 20)
    phi_term = Fraction(401, 193)
    chi_term = Fraction(999, 1000) * Fraction(191, 201) ** 2
    upper = ell_term + phi_term - chi_term
    target = Fraction(7, 2)
    if upper >= target:
        raise RuntimeError("final localized curvature comparison failed")
    return {
        "phi_inequality": "z/(exp(z)-1)<=1",
        "chi_inequality": (
            "z^2*exp(z)/(exp(z)-1)^2>=exp(-z)>=1-z>=999/1000"
        ),
        "ell_term_upper": fraction_text(ell_term),
        "phi_term_upper": fraction_text(phi_term),
        "chi_term_lower": fraction_text(chi_term),
        "scaled_U_upper": fraction_text(upper),
        "scaled_U_upper_decimal": format(float(upper), ".16f"),
        "target": fraction_text(target),
        "margin_lower": fraction_text(target - upper),
        "margin_decimal": format(float(target - upper), ".16f"),
    }


def build_artifact() -> dict:
    corridor = json.loads(SOURCE_CORRIDOR.read_text(encoding="utf-8"))
    finite = json.loads(SOURCE_FINITE.read_text(encoding="utf-8"))
    if corridor.get("summary", {}).get("global_exact_corridors_closed") is not True:
        raise RuntimeError("exact corridor source is not closed")
    if (
        finite.get("summary", {}).get("finite_corridor_to_curvature_closed")
        is not True
    ):
        raise RuntimeError("finite corridor-to-curvature source is not closed")

    geometry = potential_geometry()
    normalized = normalized_h_boxes()
    defect = logarithmic_defect_bounds()
    localized = localized_scaled_gates()
    scalar = scalar_curvature_comparison()
    rows = [
        RayRow(
            id="co4eclcrc_01_exact_corridor_input",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="All seven exact alternating-factorial cumulant corridors hold at every central and shifted mode used on the ray.",
            formula="exact cumulant corridors hold globally for u>=2",
            proof_boundary="Exact cumulant input only.",
        ),
        RayRow(
            id="co4eclcrc_02_ray_geometry",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Coefficient-positive potential geometry gives a stable u>=19 collar around every central mode u>=20.",
            formula="q>=10^33; uq<=t<=(201/100)uq; V''>=(361/100)u^2q",
            proof_boundary="Potential geometry and shifted-mode collar only.",
            diagnostics=geometry,
        ),
        RayRow(
            id="co4eclcrc_03_normalized_H_boxes",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Midpoint Hurwitz-zeta bounds and the exact corridors place all seven normalized H derivatives in positive unit boxes.",
            formula="0<x_r<=1; x_2>=97/100; x_3>=24/25, 2<=r<=8",
            proof_boundary="Normalized H derivatives through order eight only.",
            diagnostics=normalized,
        ),
        RayRow(
            id="co4eclcrc_04_logarithmic_defect",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The convergent logarithmic-defect series and exact Bell identity control ell through order six, including every shifted Taylor envelope.",
            formula="t^2 ell''<=23/20; |t^r ell^(r)|<30000; t^(r+1)E_r<1/1000",
            proof_boundary="Localized logarithmic derivatives and envelopes only.",
            diagnostics=defect,
        ),
        RayRow(
            id="co4eclcrc_05_localized_scaled_gates",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Five dimensionless gates retain the positive curvature terms and the essential negative square term.",
            formula="A_->193/100; A_+<201/100; C_+<401/100; P_->191/100",
            proof_boundary="Inputs to the scalar localized-curvature comparison only.",
            diagnostics=localized,
        ),
        RayRow(
            id="co4eclcrc_06_scalar_comparison",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="Elementary exponential inequalities reduce the complete scaled localized curvature bound to one strict rational comparison.",
            formula="2*(23/20)+401/193-(999/1000)*(191/201)^2<7/2",
            proof_boundary="Localized U ceiling only.",
            diagnostics=scalar,
        ),
        RayRow(
            id="co4eclcrc_07_ray_theorem",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The exact cumulant corridors imply a positive localized gap and U(t)<7/(2t^2) for every central mode u>=20.",
            formula="j_0>E_0 and U(t)<7/(2t^2), u>=20",
            proof_boundary="Exact first-summand localized curvature ray only.",
        ),
        RayRow(
            id="co4eclcrc_08_global_corridor_to_U",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The finite cover and analytic ray close the exact-corridor-to-localized-curvature implication for every mode u>=2.",
            formula="exact corridors => j_0>E_0 and U(t)<7/(2t^2), u>=2",
            proof_boundary="Does not yet compose the compact interval, first-summand transfer, or full-kernel order-four entry.",
        ),
    ]
    source_paths = (SOURCE_CORRIDOR, SOURCE_FINITE, SOURCE_BRIDGE)
    return {
        "kind": "jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate",
        "date": "2026-07-13",
        "status": "rigorous exact-corridor localized-curvature ray and global corridor-to-U theorem",
        "proof_boundary": (
            "This artifact proves the localized U ceiling from the exact cumulant "
            "corridors on u>=20 and composes it with the finite u<=20 theorem. "
            "It does not by itself prove complete order-four entry, PF-infinity, "
            "the Riemann hypothesis, or Lambda<=0."
        ),
        "rows": [asdict(row) for row in rows],
        "geometry": geometry,
        "normalized_h_boxes": normalized,
        "logarithmic_defect": defect,
        "localized_scaled_gates": localized,
        "scalar_comparison": scalar,
        "summary": {
            "rows": len(rows),
            "exact_rows": len(rows),
            "open_analytic_rows": 0,
            "normalized_h_boxes": 7,
            "defect_derivative_bounds": 6,
            "localized_scaled_gates": 5,
            "positive_rational_comparisons": 1,
            "ray_corridor_to_curvature_closed": True,
            "global_corridor_to_curvature_closed": True,
        },
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for path in source_paths
        ],
        "source_hashes": {
            str(path.relative_to(REPO_ROOT)).replace("\\", "/"): sha256(path)
            for path in source_paths
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.py"
        ),
        "remaining_target": (
            "Compose compact curvature, global corridor-to-U, the first-summand "
            "tent transfer, and the full-kernel perturbation into contiguous "
            "order-four entry at lambda=-100."
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    geometry = artifact["geometry"]
    scalar = artifact["scalar_comparison"]
    lines = [
        "# Jensen-Window PF Compound Order-Four Exact-Corridor Localized-Curvature Ray Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous exact-corridor localized-curvature theorem on `u>=20`,",
        "composed with the finite theorem into a global corridor-to-`U` theorem on",
        "`u>=2`. This is not a proof of order-four entry, PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.py",
        "```",
        "",
        "## Ray Geometry",
        "",
        "Write `q=pi*exp(4u)`, `t=V'(x_u)`, and `a=V''(x_u)`. On the",
        "slightly enlarged mode ray `u>=19`, exact coefficient-positive gates and",
        "the monotonicity of `q/u` give",
        "",
        "```text",
        f"q >= {geometry['q_floor']}",
        "u*q <= t <= (201/100)*u*q",
        "a >= (361/100)*u^2*q",
        "1/t <= 10^-30.",
        "```",
        "",
        "The strict endpoint comparison",
        "",
        "```text",
        f"t(20)-t(19) > {geometry['collar_t20_minus_t19_lower']}",
        "```",
        "",
        "shows that every `t+s`, `|s|<=1`, belonging to a central mode `u>=20`",
        "has a shifted mode above `19`. Thus the same analytic boxes control all",
        "three Taylor envelopes without extrapolating a sampled cover.",
        "",
        "## Normalized H Jets",
        "",
        "For `2<=r<=8` put",
        "",
        "```text",
        "x_r=(-1)^r*t^(r-1)*H^(r)/(r-2)!.",
        "```",
        "",
        "Convex midpoint quadrature and the integral lower bound for the Hurwitz",
        "zeta function give",
        "",
        "```text",
        "(t/(t+1/2))^(r-1) <= (r-1)t^(r-1)zeta(r,t+1/2) <= 1.",
        "```",
        "",
        "Combining this with the seven exact cumulant corridors and the ray",
        "geometry proves",
        "",
        "```text",
        "0 < x_r <= 1       (2<=r<=8)",
        "x_2 >= 97/100",
        "x_3 >= 24/25.",
        "```",
        "",
        "## Logarithmic Defect",
        "",
        "With `B=H''` split",
        "",
        "```text",
        "ell=log(1-exp(-B))=log(B)+R(B)",
        "R(B)=log((1-exp(-B))/B).",
        "```",
        "",
        "The convergent product-series expansion, `zeta(2k)<2`, and `2*pi>6`",
        "give `|R^(m)(B)|<1` for `1<=m<=6`. The exact partial-Bell identity",
        "then yields",
        "",
        "```text",
        "t^2*ell'' <= 23/20",
        "|t^r*ell^(r)| < 30000       (2<=r<=6)",
        "t^(r+1)*E_r < 1/1000        (r=0,1,2).",
        "```",
        "",
        "## Curvature Comparison",
        "",
        "The scaled localized quantities satisfy",
        "",
        "```text",
        "A_-=t*(j_0-E_0)                 > 193/100",
        "A_+=t*(j_0+E_0)                 < 201/100",
        "C_+=t^3*max(j_2+E_2,0)          < 401/100",
        "P_-=t^2*max(abs(j_1)-E_1,0)     > 191/100.",
        "```",
        "",
        "For `0<z<=1/1000`,",
        "",
        "```text",
        "z/(exp(z)-1) <= 1",
        "z^2*exp(z)/(exp(z)-1)^2 >= exp(-z) >= 1-z >= 999/1000.",
        "```",
        "",
        "Crucially, the negative square term is retained. The final comparison is",
        "",
        "```text",
        "t^2*U(t)",
        " < 2*(23/20) + 401/193 - (999/1000)*(191/201)^2",
        f" = {scalar['scaled_U_upper']}",
        f" = {scalar['scaled_U_upper_decimal']}",
        " < 7/2,",
        f"margin > {scalar['margin_lower']}.",
        "```",
        "",
        "Therefore `j_0>E_0` and `U(t)<7/(2t^2)` for every central mode",
        "`u>=20`. The checked finite certificate supplies the same result on",
        "`2<=u<=20`, so the exact-corridor-to-localized-curvature implication is",
        "now complete for every `u>=2`.",
        "",
        "## Next Composition",
        "",
        "The next step is bookkeeping with mathematical content already isolated:",
        "compose the compact curvature theorem, this global theorem, the exact tent",
        "transfer, and the full-kernel perturbation into contiguous order-four entry",
        "at `lambda=-100`.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.md",
        "outputs/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate.md",
        "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md",
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
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "certified order-four exact-corridor curvature ray: "
        f"{summary['rows']} rows, {summary['exact_rows']} exact rows, "
        f"{summary['normalized_h_boxes']} normalized H boxes, "
        f"{summary['defect_derivative_bounds']} defect bounds, "
        f"{summary['localized_scaled_gates']} localized gates, "
        "global corridor-to-curvature closed"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
