#!/usr/bin/env python3
"""Independently validate the signed universal-factor residual gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_signed_universal_factor_residual_gate as target


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_signed_universal_factor_residual_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_signed_universal_factor_residual_gate.md"

EXPECTED_IDS = [
    "sufr_01_multiplier_quartic",
    "sufr_02_rational_rectangle",
    "sufr_03_discriminant_exclusions",
    "sufr_04_bivariate_laguerre",
    "sufr_05_acb_spectral_jets",
    "sufr_06_adaptive_rectangle_cover",
    "sufr_07_signed_universal_obstruction",
    "sufr_08_coupled_handoff",
]


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = target.build_payload()
    if stored != rebuilt:
        issues.append("stored signed universal-factor payload differs from reconstruction")
    if stored.get("status") != (
        "exact and interval-certified signed universal-factor residual obstruction"
    ):
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 8:
        issues.append("expected 8 rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    expected_roles = {
        "exact_multiplier_reduction": 1,
        "exact_interval_reduction": 1,
        "exact_route_guard": 1,
        "exact_identity": 1,
        "interval_method": 1,
        "interval_theorem": 1,
        "exact_interval_composition": 1,
        "open_handoff": 1,
    }
    for role, expected in expected_roles.items():
        if sum(row.get("role") == role for row in rows) != expected:
            issues.append(f"expected {expected} row(s) with role {role}")

    exact = stored.get("exact", {})
    beta, gamma, c, y = sp.symbols("beta gamma c y", real=True)
    quartic = (
        256 * y**4
        - 576 * y**3
        + (432 + 16 * beta) * y**2
        - (120 + 20 * beta) * y
        + 9
        + 5 * beta
        + gamma
    )
    multiplier = sp.chebyshevt(9, c) + beta * sp.chebyshevt(5, c) + gamma * c
    if sp.expand(multiplier - c * quartic.subs(y, c**2)) != 0:
        issues.append("independent multiplier quartic reduction failed")
    reduction = exact.get("multiplier_reduction", {})
    for phrase in ("cosh(z)*Q", "iff every zero", "[0,1]"):
        joined = " ".join(str(value) for value in reduction.values())
        if phrase not in joined:
            issues.append(f"multiplier reduction missing: {phrase}")

    critical_factor = 32 * beta**3 - 297 * beta**2 + 1053 * beta - 1215
    critical_discriminant = sp.factor(sp.discriminant(sp.diff(quartic, y), y))
    if sp.expand(critical_discriminant + sp.Integer(2) ** 22 * critical_factor) != 0:
        issues.append("critical-point discriminant failed")
    quartic_factor = sp.factor(sp.discriminant(quartic, y) / sp.Integer(2) ** 24)
    if sp.Poly(quartic_factor, beta, gamma).total_degree() != 5:
        issues.append("quartic discriminant factor degree drifted")
    if sp.discriminant(sp.diff(critical_factor, beta), beta) != -51516:
        issues.append("critical-factor monotonicity failed")
    if sp.factor(critical_factor.subs(beta, sp.Rational(11, 5))) != sp.Rational(607, 125):
        issues.append("critical-factor rational threshold failed")
    discriminants = exact.get("discriminants", {})
    for phrase in ("disc(Q')<0", "disc(Q)<0", "four roots in [0,1]"):
        if phrase not in discriminants.get("exclusions", ""):
            issues.append(f"stored discriminant exclusion missing: {phrase}")

    parameter = exact.get("parameter_rectangle", {})
    if parameter.get("rectangle") != (
        "-17/5<=beta<11/5, -1<=s=beta+gamma<51/10"
    ):
        issues.append("signed parameter rectangle drifted")
    for phrase in ("Q'(1)=40+12*beta", "f(11/5)=607/125", "Q(1)=1+beta+gamma"):
        joined = " ".join(str(value) for value in parameter.values())
        if phrase not in joined:
            issues.append(f"parameter reduction missing: {phrase}")

    origin = stored.get("origin_certificate", {})
    try:
        phi_upper = target.arb(origin.get("phi_upper", "nan"))
        comparison_upper = target.arb(origin.get("comparison_upper", "nan"))
    except (TypeError, ValueError):
        phi_upper = comparison_upper = target.arb("nan")
        issues.append("origin bounds are unparseable")
    if not phi_upper < comparison_upper:
        issues.append("origin upper-sum certificate failed")

    interval = stored.get("interval_certificate", {})
    if interval.get("rigorous") is not True:
        issues.append("adaptive certificate is not marked rigorous")
    if interval.get("arithmetic") != (
        "python-flint Acb/Arb balls at 170 decimal digits"
    ):
        issues.append("adaptive arithmetic metadata drifted")
    if interval.get("spectral_points") != ["86", "122"]:
        issues.append("signed spectral points drifted")
    for phrase in ("h^2*M/r^3", "2*h^2*M/r^4", "Cauchy"):
        if phrase not in interval.get("derivative_method", ""):
            issues.append(f"Cauchy method missing: {phrase}")

    coefficient_names = (
        "c00",
        "c_beta",
        "c_gamma",
        "c_beta2",
        "c_beta_gamma",
        "c_gamma2",
    )
    coefficients_by_x: dict[int, tuple[target.arb, ...]] = {}
    certificate_rows = interval.get("rows", [])
    if [row.get("x") for row in certificate_rows] != ["86", "122"]:
        issues.append("signed coefficient rows drifted")
    for row in certificate_rows:
        try:
            x = int(row["x"])
            jets = {
                key: [target.arb(value) for value in row["jets"][key]]
                for key in ("E", "P5", "P1")
            }
            coefficients = tuple(
                target.arb(row["laguerre_coefficients"][name])
                for name in coefficient_names
            )
        except (KeyError, TypeError, ValueError):
            issues.append(f"unparseable signed spectral row at x={row.get('x')}")
            continue
        coefficients_by_x[x] = coefficients
        recomputed = target.prior.laguerre_coefficients(
            jets["E"], jets["P5"], jets["P1"]
        )
        if not all(left.overlaps(right) for left, right in zip(coefficients, recomputed)):
            issues.append(f"signed Laguerre coefficients do not match jets at x={x}")

    cover = interval.get("adaptive_cover", {})
    expected_counts = {
        "laguerre_x86": 2329,
        "laguerre_x122": 1281,
        "critical_discriminant": 240,
        "quartic_discriminant": 244,
    }
    if cover.get("initial_boxes") != 3416:
        issues.append("adaptive base-box count drifted")
    if cover.get("leaf_count") != 4094:
        issues.append("adaptive leaf count drifted")
    if cover.get("maximum_depth") != 6:
        issues.append("adaptive maximum depth drifted")
    if cover.get("classification_counts") != expected_counts:
        issues.append("adaptive classification counts drifted")
    if cover.get("unresolved_boxes") != 0:
        issues.append("adaptive cover has unresolved boxes")
    expected_depths = {
        "0": 3400,
        "1": 41,
        "2": 54,
        "3": 86,
        "4": 197,
        "5": 252,
        "6": 64,
    }
    if cover.get("depth_counts") != expected_depths:
        issues.append("adaptive depth histogram drifted")

    if set(coefficients_by_x) == {86, 122}:
        Box = tuple[Fraction, Fraction, Fraction, Fraction, int]

        def classify(box: Box) -> tuple[str, target.arb] | None:
            beta_left, beta_right, sum_left, sum_right, _ = box
            beta_box = target.fraction_ball(beta_left, beta_right)
            sum_box = target.fraction_ball(sum_left, sum_right)
            gamma_box = sum_box - beta_box
            value = target.prior.evaluate_laguerre_polynomial(
                coefficients_by_x[86], beta_box, gamma_box
            )
            if value < 0:
                return "laguerre_x86", value
            value = target.prior.evaluate_laguerre_polynomial(
                coefficients_by_x[122], beta_box, gamma_box
            )
            if value < 0:
                return "laguerre_x122", value
            value = target.derivative_discriminant_factor(beta_box)
            if value > 0:
                return "critical_discriminant", value
            value = target.quartic_discriminant_factor(beta_box, gamma_box)
            if value < 0:
                return "quartic_discriminant", value
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

        counts = {key: 0 for key in expected_counts}
        depths: dict[int, int] = {}
        digest = hashlib.sha256()
        leaves = 0
        maximum_depth = 0
        while stack:
            box = stack.pop()
            classified = classify(box)
            if classified is not None:
                kind, _ = classified
                counts[kind] += 1
                leaves += 1
                depth = box[4]
                maximum_depth = max(maximum_depth, depth)
                depths[depth] = depths.get(depth, 0) + 1
                beta_left, beta_right, sum_left, sum_right, _ = box
                digest.update(
                    (
                        f"{target.fraction_text(beta_left)},"
                        f"{target.fraction_text(beta_right)},"
                        f"{target.fraction_text(sum_left)},"
                        f"{target.fraction_text(sum_right)},{depth}:{kind}\n"
                    ).encode("ascii")
                )
                continue
            beta_left, beta_right, sum_left, sum_right, depth = box
            if depth >= 14:
                issues.append("independent replay reached unresolved depth limit")
                break
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

        if leaves != 4094 or counts != expected_counts or maximum_depth != 6:
            issues.append("independent adaptive replay counts failed")
        if {str(key): value for key, value in sorted(depths.items())} != expected_depths:
            issues.append("independent adaptive replay depths failed")
        if digest.hexdigest() != cover.get("leaf_sha256"):
            issues.append("adaptive leaf digest drifted")

    for name, value_text in cover.get("closest_certified_bounds", {}).items():
        try:
            value = target.arb(value_text)
        except (TypeError, ValueError):
            issues.append(f"unparseable closest adaptive bound: {name}")
            continue
        if name == "critical_discriminant":
            if not value > 0:
                issues.append("closest critical-discriminant bound is not positive")
        elif not value < 0:
            issues.append(f"closest negative adaptive bound failed: {name}")

    note = NOTE.read_text(encoding="utf-8")
    required_note = (
        "U_(beta,gamma)(z)=cosh(9z)",
        "-17/5<=beta<11/5",
        "initial boxes=3416",
        "adaptive leaves=4094",
        "maximum depth=6",
        "L(x=86) leaves=2329",
        "L(x=122) leaves=1281",
        "unresolved boxes=0",
        "This is not a proof or disproof of RH",
    )
    for phrase in required_note:
        if phrase not in note:
            issues.append(f"rendered note missing: {phrase}")
    for forbidden in ("therefore RH", "we have disproved RH", "Lambda <= 0 is proved"):
        if forbidden.lower() in note.lower():
            issues.append(f"rendered note contains forbidden promotion: {forbidden}")
    return issues


def main() -> int:
    issues = validate()
    for issue in issues:
        print(f"SIGNED-UNIVERSAL-FACTOR [{issue}]")
    print(
        "validated Jensen-window PF Newman signed universal-factor residual "
        f"gate: 8 rows, {len(issues)} issues, 1 exact multiplier reduction, "
        "1 rational parameter rectangle, 2 discriminant exclusions, 2 Acb "
        "spectral certificates, 4094 adaptive leaves, 3416 base boxes, "
        "maximum depth 6, 1 exhaustive signed universal-factor obstruction, "
        "1 coupled handoff"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
