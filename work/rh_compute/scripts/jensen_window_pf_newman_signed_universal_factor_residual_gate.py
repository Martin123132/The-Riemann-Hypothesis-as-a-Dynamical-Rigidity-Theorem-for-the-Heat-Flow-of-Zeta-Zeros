#!/usr/bin/env python3
"""Audit signed 9/5/1 Polya universal-factor residual decompositions."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp

import jensen_window_pf_newman_classical_three_block_residual_gate as prior


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

from flint import acb, arb, ctx  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_signed_universal_factor_residual_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_signed_universal_factor_residual_gate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fraction_ball(left: Fraction, right: Fraction) -> arb:
    midpoint = (left + right) / 2
    radius = (right - left) / 2
    midpoint_ball = arb(midpoint.numerator) / midpoint.denominator
    radius_ball = arb(radius.numerator) / radius.denominator
    return arb(midpoint_ball, radius_ball)


def derivative_discriminant_factor(beta: arb) -> arb:
    return ((32 * beta - 297) * beta + 1053) * beta - 1215


def quartic_discriminant_factor(beta: arb, gamma: arb) -> arb:
    """Evaluate disc(Q)/2^24 in nested form to limit dependency loss."""

    gamma2_coefficient = (-128 * beta - 48) * beta - 459
    gamma_coefficient = (
        (((16 * beta + 124) * beta - 384) * beta + 1080) * beta - 486
    )
    constant = (
        ((((-20 * beta + 64) * beta - 184) * beta + 432) * beta - 972)
        * beta
        + 729
    )
    return (
        (256 * gamma + gamma2_coefficient) * gamma + gamma_coefficient
    ) * gamma + constant


def build_exact() -> dict:
    beta, gamma, c, y = sp.symbols("beta gamma c y", real=True)
    q_poly = (
        256 * y**4
        - 576 * y**3
        + (432 + 16 * beta) * y**2
        - (120 + 20 * beta) * y
        + 9
        + 5 * beta
        + gamma
    )
    multiplier_polynomial = (
        sp.chebyshevt(9, c) + beta * sp.chebyshevt(5, c) + gamma * c
    )
    if sp.expand(multiplier_polynomial - c * q_poly.subs(y, c**2)) != 0:
        raise RuntimeError("9/5/1 Chebyshev reduction failed")

    derivative_factor = 32 * beta**3 - 297 * beta**2 + 1053 * beta - 1215
    derivative_discriminant = sp.factor(sp.discriminant(sp.diff(q_poly, y), y))
    if sp.expand(
        derivative_discriminant + sp.Integer(2) ** 22 * derivative_factor
    ) != 0:
        raise RuntimeError("critical-point discriminant failed")

    quartic_factor = (
        -20 * beta**5
        + 16 * beta**4 * gamma
        + 64 * beta**4
        + 124 * beta**3 * gamma
        - 184 * beta**3
        - 128 * beta**2 * gamma**2
        - 384 * beta**2 * gamma
        + 432 * beta**2
        - 48 * beta * gamma**2
        + 1080 * beta * gamma
        - 972 * beta
        + 256 * gamma**3
        - 459 * gamma**2
        - 486 * gamma
        + 729
    )
    quartic_discriminant = sp.factor(sp.discriminant(q_poly, y))
    if sp.expand(quartic_discriminant - sp.Integer(2) ** 24 * quartic_factor) != 0:
        raise RuntimeError("quartic discriminant failed")

    derivative_of_factor = sp.diff(derivative_factor, beta)
    derivative_of_factor_discriminant = sp.discriminant(derivative_of_factor, beta)
    if derivative_of_factor_discriminant != -51516:
        raise RuntimeError("monotonicity discriminant failed")
    if sp.factor(derivative_factor.subs(beta, sp.Rational(11, 5))) != sp.Rational(607, 125):
        raise RuntimeError("rational upper beta threshold failed")
    if sp.expand(sp.diff(q_poly, y).subs(y, 1) - (40 + 12 * beta)) != 0:
        raise RuntimeError("Q prime endpoint failed")
    if sp.expand(q_poly.subs(y, 1) - (1 + beta + gamma)) != 0:
        raise RuntimeError("Q endpoint failed")

    return {
        "multiplier_reduction": {
            "multiplier": "U_(beta,gamma)(z)=cosh(9z)+beta*cosh(5z)+gamma*cosh(z)",
            "reduction": (
                "U_(beta,gamma)(z)=cosh(z)*Q_(beta,gamma)(cosh(z)^2)"
            ),
            "quartic": f"Q_(beta,gamma)(y)={sp.sstr(q_poly)}",
            "imaginary_zero_equivalence": (
                "U_(beta,gamma) has only imaginary zeros iff every zero of "
                "Q_(beta,gamma) is real and lies in [0,1]."
            ),
            "reason": (
                "cosh(z)=0 is purely imaginary, and cosh(z)^2=r has only "
                "purely imaginary solutions exactly when 0<=r<=1."
            ),
        },
        "parameter_rectangle": {
            "lower_beta": (
                "If Q has roots in [0,1], Q'(1)=40+12*beta>=0, so "
                "beta>=-10/3>-17/5."
            ),
            "upper_beta": (
                "Rolle gives three real roots of Q', hence f(beta)<=0 for "
                "f=32beta^3-297beta^2+1053beta-1215. Since disc(f')=-51516, "
                "f is strictly increasing, and f(11/5)=607/125>0, beta<11/5."
            ),
            "lower_sum": (
                "Q(1)=1+beta+gamma=256*product_j(1-r_j)>=0, so "
                "s=beta+gamma>=-1."
            ),
            "upper_sum": (
                "Global residual nonnegativity at u=0 and the Arb origin "
                "bound give s<51/10."
            ),
            "rectangle": "-17/5<=beta<11/5, -1<=s=beta+gamma<51/10",
            "scope": (
                "A rational outer rectangle containing every signed 9/5/1 "
                "multiplier that passes the standard imaginary-zero test and "
                "leaves a globally nonnegative residual."
            ),
        },
        "discriminants": {
            "critical_factor": f"f(beta)={sp.sstr(derivative_factor)}",
            "critical_discriminant": (
                f"disc_y Q'={sp.sstr(derivative_discriminant)}"
            ),
            "quartic_factor": f"Delta(beta,gamma)={sp.sstr(quartic_factor)}",
            "quartic_discriminant": (
                f"disc_y Q=2^24*Delta={sp.sstr(quartic_discriminant)}"
            ),
            "exclusions": (
                "f(beta)>0 makes disc(Q')<0, while Delta(beta,gamma)<0 "
                "makes disc(Q)<0; either condition excludes four roots in [0,1]."
            ),
        },
        "laguerre_identity": {
            "identity": (
                "L[R_(beta,gamma)]=L[E]-beta*B[E,P5]-gamma*B[E,P1]+"
                "beta^2*L[P5]+beta*gamma*B[P5,P1]+gamma^2*L[P1]"
            ),
            "definitions": "E=H_0-P_(9/4), P5=P_(5/4), P1=P_(1/4)",
        },
    }


def build_interval_certificate() -> dict:
    ctx.dps = 170
    step = arb("1e-6")
    cauchy_radius = arb("0.04")
    box_radius = arb("0.1")
    spectral_points = (86, 122)

    def e_function(x: acb) -> acb:
        return prior.acb_h0(x) - prior.acb_gasper_block(x, arb("9/4"))

    def p5_function(x: acb) -> acb:
        return prior.acb_gasper_block(x, arb("5/4"))

    def p1_function(x: acb) -> acb:
        return prior.acb_gasper_block(x, arb("1/4"))

    rows: list[dict] = []
    coefficients_by_x: dict[int, tuple[arb, arb, arb, arb, arb, arb]] = {}
    for x in spectral_points:
        e_jet, e_maximum = prior.cauchy_finite_difference_jet(
            e_function, str(x), step, cauchy_radius, box_radius
        )
        p5_jet, p5_maximum = prior.cauchy_finite_difference_jet(
            p5_function, str(x), step, cauchy_radius, box_radius
        )
        p1_jet, p1_maximum = prior.cauchy_finite_difference_jet(
            p1_function, str(x), step, cauchy_radius, box_radius
        )
        raw_coefficients = prior.laguerre_coefficients(e_jet, p5_jet, p1_jet)
        coefficient_strings = [str(value) for value in raw_coefficients]
        coefficients = tuple(arb(value) for value in coefficient_strings)
        coefficients_by_x[x] = coefficients
        rows.append(
            {
                "x": str(x),
                "jets": {
                    "E": [str(value) for value in e_jet],
                    "P5": [str(value) for value in p5_jet],
                    "P1": [str(value) for value in p1_jet],
                },
                "cauchy_box_abs_upper": {
                    "E": str(e_maximum),
                    "P5": str(p5_maximum),
                    "P1": str(p1_maximum),
                },
                "laguerre_coefficients": {
                    name: value
                    for name, value in zip(
                        (
                            "c00",
                            "c_beta",
                            "c_gamma",
                            "c_beta2",
                            "c_beta_gamma",
                            "c_gamma2",
                        ),
                        coefficient_strings,
                    )
                },
            }
        )

    Box = tuple[Fraction, Fraction, Fraction, Fraction, int]

    def classify(box: Box) -> tuple[str, arb] | None:
        beta_left, beta_right, sum_left, sum_right, _ = box
        beta_box = fraction_ball(beta_left, beta_right)
        sum_box = fraction_ball(sum_left, sum_right)
        gamma_box = sum_box - beta_box
        value_86 = prior.evaluate_laguerre_polynomial(
            coefficients_by_x[86], beta_box, gamma_box
        )
        if value_86 < 0:
            return "laguerre_x86", value_86
        value_122 = prior.evaluate_laguerre_polynomial(
            coefficients_by_x[122], beta_box, gamma_box
        )
        if value_122 < 0:
            return "laguerre_x122", value_122
        critical = derivative_discriminant_factor(beta_box)
        if critical > 0:
            return "critical_discriminant", critical
        quartic = quartic_discriminant_factor(beta_box, gamma_box)
        if quartic < 0:
            return "quartic_discriminant", quartic
        return None

    stack: list[Box] = []
    for beta_index in reversed(range(56)):
        for sum_index in reversed(range(61)):
            stack.append(
                (
                    Fraction(-17, 5) + Fraction(beta_index, 10),
                    Fraction(-17, 5) + Fraction(beta_index + 1, 10),
                    Fraction(-1) + Fraction(sum_index, 10),
                    Fraction(-1) + Fraction(sum_index + 1, 10),
                    0,
                )
            )

    classifications = (
        "laguerre_x86",
        "laguerre_x122",
        "critical_discriminant",
        "quartic_discriminant",
    )
    counts = {name: 0 for name in classifications}
    depth_counts: dict[int, int] = {}
    closest_bounds: dict[str, object | None] = {name: None for name in classifications}
    digest = hashlib.sha256()
    leaves = 0
    maximum_depth = 0
    depth_limit = 14

    while stack:
        box = stack.pop()
        classification = classify(box)
        if classification is not None:
            kind, value = classification
            leaves += 1
            depth = box[4]
            maximum_depth = max(maximum_depth, depth)
            counts[kind] += 1
            depth_counts[depth] = depth_counts.get(depth, 0) + 1
            bound = value.lower() if kind == "critical_discriminant" else value.upper()
            previous = closest_bounds[kind]
            if previous is None:
                closest_bounds[kind] = bound
            elif kind == "critical_discriminant" and bound < previous:
                closest_bounds[kind] = bound
            elif kind != "critical_discriminant" and bound > previous:
                closest_bounds[kind] = bound
            beta_left, beta_right, sum_left, sum_right, _ = box
            digest.update(
                (
                    f"{fraction_text(beta_left)},{fraction_text(beta_right)},"
                    f"{fraction_text(sum_left)},{fraction_text(sum_right)},"
                    f"{depth}:{kind}\n"
                ).encode("ascii")
            )
            continue

        beta_left, beta_right, sum_left, sum_right, depth = box
        if depth >= depth_limit:
            raise RuntimeError(
                "unresolved adaptive box: "
                f"beta=[{beta_left},{beta_right}], sum=[{sum_left},{sum_right}]"
            )
        beta_midpoint = (beta_left + beta_right) / 2
        sum_midpoint = (sum_left + sum_right) / 2
        next_depth = depth + 1
        children: list[Box] = [
            (beta_left, beta_midpoint, sum_left, sum_midpoint, next_depth),
            (beta_left, beta_midpoint, sum_midpoint, sum_right, next_depth),
            (beta_midpoint, beta_right, sum_left, sum_midpoint, next_depth),
            (beta_midpoint, beta_right, sum_midpoint, sum_right, next_depth),
        ]
        stack.extend(reversed(children))

    if leaves != 4094 or maximum_depth != 6:
        raise RuntimeError("adaptive cover size drifted")
    expected_counts = {
        "laguerre_x86": 2329,
        "laguerre_x122": 1281,
        "critical_discriminant": 240,
        "quartic_discriminant": 244,
    }
    if counts != expected_counts:
        raise RuntimeError("adaptive classification counts drifted")

    for row in rows:
        x = int(row["x"])
        row["adaptive_leaf_count"] = counts[f"laguerre_x{x}"]

    return {
        "rigorous": True,
        "arithmetic": "python-flint Acb/Arb balls at 170 decimal digits",
        "derivative_method": (
            "Central differences use h=1e-6. Acb boxes of radius 0.1 bound "
            "each component on Cauchy circles of radius 0.04; the added f' "
            "and f'' errors are h^2*M/r^3 and 2*h^2*M/r^4."
        ),
        "step": "1e-6",
        "cauchy_radius": "0.04",
        "acb_box_radius": "0.1",
        "spectral_points": [str(value) for value in spectral_points],
        "rows": rows,
        "adaptive_cover": {
            "coordinates": "beta and s=beta+gamma",
            "rectangle": "[-17/5,11/5] x [-1,51/10]",
            "initial_cell_width": "1/10",
            "initial_beta_cells": 56,
            "initial_sum_cells": 61,
            "initial_boxes": 3416,
            "classification_order": list(classifications),
            "classification_logic": (
                "A leaf is accepted by L[R](86)<0, then L[R](122)<0, then "
                "f(beta)>0 (so disc Q'<0), then Delta(beta,gamma)<0 "
                "(so disc Q<0). Otherwise it is bisected in both coordinates."
            ),
            "leaf_count": leaves,
            "classification_counts": counts,
            "depth_counts": {str(key): value for key, value in sorted(depth_counts.items())},
            "maximum_depth": maximum_depth,
            "depth_limit": depth_limit,
            "closest_certified_bounds": {
                key: str(value) for key, value in closest_bounds.items()
            },
            "leaf_sha256": digest.hexdigest(),
            "unresolved_boxes": 0,
        },
        "conclusion": (
            "Every point in the necessary signed rectangle either has a "
            "strictly negative residual Laguerre value at x=86 or x=122, or "
            "is rigorously excluded from the Polya imaginary-zero multiplier "
            "cone by a negative critical or quartic discriminant."
        ),
        "scope": (
            "This exhausts signed 9/5/1 bases certified by the standard Polya "
            "universal-factor hypothesis when their residual kernel is globally "
            "nonnegative and is treated as an independent Laguerre-Polya block. "
            "It does not reject higher shifts or a coupled signed identity."
        ),
    }


def build_payload() -> dict:
    exact = build_exact()
    origin = prior.build_origin_certificate()
    interval = build_interval_certificate()
    rows = [
        GateRow(
            id="sufr_01_multiplier_quartic",
            role="exact_multiplier_reduction",
            readiness="available_exact",
            claim="The signed 9/5/1 multiplier condition is exactly a quartic root-location condition on [0,1].",
            formula=exact["multiplier_reduction"]["reduction"],
            proof_boundary="Exact standard universal-factor hypothesis only; no claim about every possible real-zero proof.",
            diagnostics=exact["multiplier_reduction"],
        ),
        GateRow(
            id="sufr_02_rational_rectangle",
            role="exact_interval_reduction",
            readiness="ready_to_apply",
            claim="Every signed universal-factor candidate with globally nonnegative residual lies in one rational outer rectangle.",
            formula=exact["parameter_rectangle"]["rectangle"],
            proof_boundary="Necessary outer bounds only; membership in the rectangle is not sufficient.",
            diagnostics={"parameter": exact["parameter_rectangle"], "origin": origin},
        ),
        GateRow(
            id="sufr_03_discriminant_exclusions",
            role="exact_route_guard",
            readiness="guard_validated",
            claim="Two exact discriminants exclude boxes that cannot contain an imaginary-zero multiplier.",
            formula=exact["discriminants"]["exclusions"],
            proof_boundary="Necessary real-rootedness exclusions for Q and Q' only.",
            diagnostics=exact["discriminants"],
        ),
        GateRow(
            id="sufr_04_bivariate_laguerre",
            role="exact_identity",
            readiness="available_exact",
            claim="The signed residual first Laguerre expression remains an exact bivariate quadratic.",
            formula=exact["laguerre_identity"]["identity"],
            proof_boundary="Exact algebra used by the adaptive interval cover.",
            diagnostics=exact["laguerre_identity"],
        ),
        GateRow(
            id="sufr_05_acb_spectral_jets",
            role="interval_method",
            readiness="interval_validated",
            claim="Acb Cauchy certificates enclose the two residual Laguerre quadratics used by the signed cover.",
            formula=interval["derivative_method"],
            proof_boundary="Finite spectral points x=86 and x=122 at t=0.",
            diagnostics={key: interval[key] for key in ("arithmetic", "spectral_points", "rows")},
        ),
        GateRow(
            id="sufr_06_adaptive_rectangle_cover",
            role="interval_theorem",
            readiness="interval_validated",
            claim="A 4,094-leaf adaptive rational cover classifies every point of the necessary signed rectangle.",
            formula=interval["conclusion"],
            proof_boundary="Complete for the stated rectangle and four classification tests; not higher-block structure.",
            diagnostics=interval["adaptive_cover"],
        ),
        GateRow(
            id="sufr_07_signed_universal_obstruction",
            role="exact_interval_composition",
            readiness="ready_to_apply",
            claim="No signed 9/5/1 base certified by the standard Polya factor theorem leaves a globally positive independently-LP residual.",
            formula=interval["scope"],
            proof_boundary="Rejects this decomposition architecture only; it does not prove or disprove RH.",
            diagnostics={"parameter": exact["parameter_rectangle"], "cover": interval["adaptive_cover"]},
        ),
        GateRow(
            id="sufr_08_coupled_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The surviving Gasper route needs higher shifts or one coupled signed square rather than independent LP summands.",
            formula="Control all mixed terms after modular cancellation in one matrix-valued identity.",
            proof_boundary="Open theorem target; not strict Laguerre positivity, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_signed_universal_factor_residual_gate",
        "date": "2026-07-11",
        "status": "exact and interval-certified signed universal-factor residual obstruction",
        "proof_boundary": (
            "This artifact reduces the signed 9/5/1 Polya multiplier hypothesis "
            "to a quartic root-location condition, derives a rational rectangle "
            "containing every globally positive residual candidate, and covers "
            "that rectangle by two Acb Laguerre tests and two exact discriminant "
            "exclusions. It closes an independent Laguerre-Polya residual for "
            "every signed base certified by this standard universal-factor route. "
            "Higher shifts and coupled signed identities remain open; this does "
            "not prove or disprove RH or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_classical_three_block_residual_gate.md",
            "outputs/jensen_window_pf_newman_gasper_residual_two_block_gate.md",
            "outputs/formal_core.md",
            "https://arxiv.org/abs/0801.2996",
            "https://arxiv.org/abs/1502.06844",
            "https://doi.org/10.1007/BF02565336",
        ],
        "exact": exact,
        "origin_certificate": origin,
        "interval_certificate": interval,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    interval = payload["interval_certificate"]
    cover = interval["adaptive_cover"]
    lines = [
        "# Jensen-Window PF Newman Signed Universal-Factor Residual Gate",
        "",
        "Date: 2026-07-11",
        "",
        "Status: exact and interval-certified signed universal-factor residual",
        "obstruction. This is not a proof or disproof of RH or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_newman_signed_universal_factor_residual_gate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_signed_universal_factor_residual_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_signed_universal_factor_residual_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_signed_universal_factor_residual_gate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF Newman signed universal-factor residual gate: 8 rows, 0 issues, 1 exact multiplier reduction, 1 rational parameter rectangle, 2 discriminant exclusions, 2 Acb spectral certificates, 4094 adaptive leaves, 3416 base boxes, maximum depth 6, 1 exhaustive signed universal-factor obstruction, 1 coupled handoff",
        "```",
        "",
        "## Multiplier Cone",
        "",
        "```text",
        exact["multiplier_reduction"]["multiplier"],
        exact["multiplier_reduction"]["reduction"],
        exact["multiplier_reduction"]["quartic"],
        exact["multiplier_reduction"]["imaginary_zero_equivalence"],
        "```",
        "",
        "Thus the standard Polya universal-factor hypothesis is a concrete",
        "quartic root-location condition, including signed coefficients.",
        "",
        "## Necessary Rectangle",
        "",
        "```text",
        exact["parameter_rectangle"]["lower_beta"],
        exact["parameter_rectangle"]["upper_beta"],
        exact["parameter_rectangle"]["lower_sum"],
        exact["parameter_rectangle"]["upper_sum"],
        exact["parameter_rectangle"]["rectangle"],
        "```",
        "",
        "The upper sum bound is the independently certified Xi-kernel origin",
        "inequality from the classical three-block gate.",
        "",
        "## Adaptive Certificate",
        "",
        "The residual Laguerre expression is",
        "",
        "```text",
        exact["laguerre_identity"]["identity"],
        "```",
        "",
        interval["derivative_method"],
        "The adaptive classification is:",
        "",
        "```text",
        f"initial boxes={cover['initial_boxes']}",
        f"adaptive leaves={cover['leaf_count']}",
        f"maximum depth={cover['maximum_depth']}",
        f"L(x=86) leaves={cover['classification_counts']['laguerre_x86']}",
        f"L(x=122) leaves={cover['classification_counts']['laguerre_x122']}",
        f"critical-discriminant leaves={cover['classification_counts']['critical_discriminant']}",
        f"quartic-discriminant leaves={cover['classification_counts']['quartic_discriminant']}",
        f"unresolved boxes={cover['unresolved_boxes']}",
        f"leaf sha256={cover['leaf_sha256']}",
        "```",
        "",
        exact["discriminants"]["exclusions"],
        "Every leaf is therefore either spectrally incompatible with an",
        "independently Laguerre-Polya residual or algebraically outside the",
        "standard imaginary-zero multiplier cone.",
        "",
        "## Boundary",
        "",
        interval["scope"],
        "The next route must add higher shifts or derive one genuinely coupled",
        "matrix square that retains and controls the mixed signs.",
        "",
        "References: https://arxiv.org/abs/0801.2996,",
        "https://arxiv.org/abs/1502.06844, and",
        "https://doi.org/10.1007/BF02565336.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(f"wrote Newman signed universal-factor residual gate: {len(payload['rows'])} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
