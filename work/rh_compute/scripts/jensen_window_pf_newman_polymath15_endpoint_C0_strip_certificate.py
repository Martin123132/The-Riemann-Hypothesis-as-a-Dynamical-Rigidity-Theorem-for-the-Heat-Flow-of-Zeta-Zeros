#!/usr/bin/env python3
"""Certify coarse complex-strip bounds for the Riemann-Siegel C0 factor."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.md"
)
PRECISION_BITS = 256
R_OUTER = Fraction(1, 4)
R_INNER = Fraction(3, 20)
STRIP = Fraction(1, 100)
REAL_PAD = Fraction(1, 100)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits)


def bounds() -> dict:
    flint.ctx.prec = PRECISION_BITS
    pi = flint.arb.pi()
    sqrt_two = flint.arb(2).sqrt()

    near_exp = (3 * pi / 16).exp()
    near_cos = (pi / 8).cosh()
    near_numerator = near_exp + sqrt_two * near_cos
    near_denominator = sqrt_two
    near_c0 = near_numerator / near_denominator
    radius_gap = R_OUTER - R_INNER
    near_derivative = near_c0 / flint.arb(
        flint.fmpq(radius_gap.numerator, radius_gap.denominator)
    )
    if not near_numerator < 4:
        raise RuntimeError("near-removable numerator bound failed")
    if not near_c0 < 3:
        raise RuntimeError("near-removable C0 bound failed")
    if not near_derivative < 30:
        raise RuntimeError("near-removable C0 derivative bound failed")

    real_limit = flint.arb(1) + flint.arb(
        flint.fmpq(REAL_PAD.numerator, REAL_PAD.denominator)
    )
    strip = flint.arb(flint.fmpq(STRIP.numerator, STRIP.denominator))
    outside_exp = (pi * real_limit * strip).exp()
    outside_hyperbolic = (pi * strip / 2).cosh()
    outside_numerator = outside_exp + sqrt_two * outside_hyperbolic
    z_abs = (real_limit**2 + strip**2).sqrt()
    outside_numerator_prime = (
        pi * z_abs * outside_exp
        + pi / sqrt_two * outside_hyperbolic
    )
    outside_denominator = 2 * (pi / 10).sin()
    outside_denominator_prime = 2 * pi * (pi * strip).cosh()
    outside_c0 = outside_numerator / outside_denominator
    outside_derivative = (
        outside_numerator_prime / outside_denominator
        + outside_numerator
        * outside_denominator_prime
        / outside_denominator**2
    )
    if not outside_numerator < flint.arb(5) / 2:
        raise RuntimeError("outside numerator bound failed")
    if not outside_numerator_prime < flint.arb(28) / 5:
        raise RuntimeError("outside numerator derivative bound failed")
    if not outside_denominator_prime < flint.arb(63) / 10:
        raise RuntimeError("outside denominator derivative bound failed")
    if not outside_c0 < 5:
        raise RuntimeError("outside C0 bound failed")
    if not outside_derivative < 60:
        raise RuntimeError("outside C0 derivative bound failed")

    l0 = flint.arb(50)
    p_displacement = (-l0 / 2).exp() / (10 * l0)
    if not p_displacement < strip:
        raise RuntimeError("critical collar p-displacement failed")
    return {
        "precision_bits": PRECISION_BITS,
        "domain": "|Re(p)|<=101/100, |Im(p)|<=1/100",
        "removable_centers": ["-1/2", "1/2"],
        "near_outer_radius": "1/4",
        "near_inner_radius": "3/20",
        "near_numerator_ball": arb_text(near_numerator),
        "near_denominator_lower": "sqrt(2)",
        "near_C0_ball": arb_text(near_c0),
        "near_C0_lt": "3",
        "near_C0_prime_ball": arb_text(near_derivative),
        "near_C0_prime_lt": "30",
        "outside_numerator_ball": arb_text(outside_numerator),
        "outside_numerator_prime_ball": arb_text(outside_numerator_prime),
        "outside_denominator_lower_ball": arb_text(outside_denominator),
        "outside_denominator_prime_ball": arb_text(outside_denominator_prime),
        "outside_C0_ball": arb_text(outside_c0),
        "outside_C0_lt": "5",
        "outside_C0_prime_ball": arb_text(outside_derivative),
        "outside_C0_prime_lt": "60",
        "global_bound": "|C0(p)|<5 and |C0'(p)|<60",
        "L_min": "50",
        "critical_collar_p_displacement_ball": arb_text(p_displacement),
        "critical_collar_p_displacement_lt": "1/100",
    }


def build_exact() -> dict:
    return {
        "C0": (
            "C0(p)=(exp(pi*i*(p^2/2+3/8))-i*sqrt(2)*cos(pi*p/2))/"
            "(2*cos(pi*p))"
        ),
        "removable": (
            "The numerator vanishes at p=+-1/2, so both apparent poles are removable"
        ),
        "near_lower": (
            "For z=p0+w, |w|=1/4 and p0=+-1/2, "
            "|2*cos(pi*z)|=2*|sin(pi*w)|>=sqrt(2)"
        ),
        "near_proof": (
            "sin(pi*|Re w|)>=2*sqrt(2)*|Re w| and "
            "sinh(pi*|Im w|)>=2*sqrt(2)*|Im w|"
        ),
        "cauchy": (
            "A boundary bound |C0|<3 on radius 1/4 gives |C0'|<30 "
            "on radius 3/20"
        ),
        "outside_lower": (
            "Outside the two radius-3/20 disks, the strip width 1/100 forces "
            "dist(Re p, Z+1/2)>1/10 and |2*cos(pi*p)|>=2*sin(pi/10)"
        ),
        "global": (
            "On |Re p|<=101/100 and |Im p|<=1/100, "
            "|C0(p)|<5 and |C0'(p)|<60"
        ),
        "collar_map": (
            "For L>=50 and |z-x|<=1/L, the fixed-N map "
            "p_N(z)=1-2*(sqrt((z/2+pi*t/8)/(2*pi))-N) stays in this strip"
        ),
    }


def build_artifact() -> dict:
    exact = build_exact()
    interval = bounds()
    rows = [
        GateRow(
            id="np15c0sc_01_definition",
            role="exact_definition",
            readiness="ready_to_apply",
            claim="The endpoint factor is meromorphic only in appearance at the two half-integers.",
            formula=f"{exact['C0']}; {exact['removable']}",
            proof_boundary="Principal exponential and trigonometric functions.",
        ),
        GateRow(
            id="np15c0sc_02_near_denominator",
            role="analytic_bound",
            readiness="ready_to_apply",
            claim="A fixed outer circle around each removable point has a uniform denominator lower bound.",
            formula=f"{exact['near_lower']}; {exact['near_proof']}",
            proof_boundary="Exact elementary trigonometric and hyperbolic inequalities.",
        ),
        GateRow(
            id="np15c0sc_03_near_C0",
            role="interval_certificate",
            readiness="certified",
            claim="Arb bounds C0 on both removable-point outer circles.",
            formula="|C0|<3 on |p-(+-1/2)|=1/4",
            proof_boundary="Uses the analytic numerator bound and denominator lower bound.",
            diagnostics={
                "numerator": interval["near_numerator_ball"],
                "C0": interval["near_C0_ball"],
            },
        ),
        GateRow(
            id="np15c0sc_04_near_derivative",
            role="cauchy_certificate",
            readiness="certified",
            claim="Cauchy's estimate controls C0' through neighborhoods of both removable points.",
            formula=exact["cauchy"],
            proof_boundary="Inner radius 3/20, outer radius 1/4.",
            diagnostics={"C0_prime": interval["near_C0_prime_ball"]},
        ),
        GateRow(
            id="np15c0sc_05_outside_denominator",
            role="analytic_bound",
            readiness="ready_to_apply",
            claim="Away from the removable disks the denominator is uniformly separated from zero.",
            formula=exact["outside_lower"],
            proof_boundary="Covers the padded strip outside both inner disks.",
        ),
        GateRow(
            id="np15c0sc_06_outside_quotient",
            role="interval_certificate",
            readiness="certified",
            claim="Direct quotient differentiation bounds C0 and C0' on the outside region.",
            formula="|C0|<5 and |C0'|<60",
            proof_boundary="Arb evaluates conservative elementary numerator and denominator budgets.",
            diagnostics={
                "C0": interval["outside_C0_ball"],
                "C0_prime": interval["outside_C0_prime_ball"],
            },
        ),
        GateRow(
            id="np15c0sc_07_global_strip",
            role="proved_bound",
            readiness="ready_to_apply",
            claim="The near and outside covers give one rigorous complex-strip bound.",
            formula=exact["global"],
            proof_boundary="Closed padded strip only.",
        ),
        GateRow(
            id="np15c0sc_08_collar_map",
            role="analytic_transfer",
            readiness="ready_to_apply",
            claim="Every radius-1/L critical collar with L at least 50 maps into the certified p-strip.",
            formula=exact["collar_map"],
            proof_boundary="Fixed N and principal square-root branch.",
            diagnostics={
                "p_displacement": interval["critical_collar_p_displacement_ball"],
                "bound": interval["critical_collar_p_displacement_lt"],
            },
        ),
        GateRow(
            id="np15c0sc_09_lift_handoff",
            role="closure_record",
            readiness="closed",
            claim="The derivative input left open by the endpoint holomorphic-lift artifact is now explicit.",
            formula="complex-strip |C0|<5, |C0'|<60",
            proof_boundary="Closes the local C0 bound, not the corrected first-jet lower bound.",
        ),
        GateRow(
            id="np15c0sc_10_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="A correction-factor derivative bound does not imply Newman transversality.",
            formula="C0 strip control != T_L[J] lower bound",
            proof_boundary="Not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate",
        "date": "2026-07-17",
        "status": (
            "rigorous complex-strip bounds for the endpoint C0 factor and its "
            "first derivative; not corrected transversality, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This artifact proves |C0|<5 and |C0'|<60 on the stated complex "
            "strip and proves that L>=50 radius-1/L endpoint collars map into it. "
            "It does not by itself complete the full refined remainder constants, "
            "prove the corrected first-jet lower bound, exclude a positive Newman "
            "boundary, prove Lambda<=0, or prove RH."
        ),
        "exact": exact,
        "interval": interval,
        "rows": [asdict(row) for row in rows],
        "sources": [
            "outputs/jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.md",
            "https://arxiv.org/abs/1904.12438",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    interval = artifact["interval"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Endpoint C0 Strip Certificate",
            "",
            "Date: 2026-07-17",
            "",
            "Status: rigorous complex-strip bounds for the endpoint factor. This",
            "is not a proof of `Lambda <= 0` or RH; corrected transversality remains open.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.py",
            "```",
            "",
            "## Domain",
            "",
            "```text",
            exact["C0"],
            interval["domain"],
            "```",
            "",
            "The apparent poles at `p=+-1/2` are removable. They are handled by",
            "Cauchy's estimate rather than unstable quotient evaluation.",
            "",
            "## Removable Disks",
            "",
            "On each radius-`1/4` outer circle,",
            "",
            "```text",
            exact["near_lower"],
            exact["near_proof"],
            f"numerator bound = {interval['near_numerator_ball']} < 4",
            f"|C0| bound      = {interval['near_C0_ball']} < 3",
            "```",
            "",
            "The gap from outer radius `1/4` to inner radius `3/20` is `1/10`,",
            "so Cauchy gives",
            "",
            "```text",
            f"|C0'| <= {interval['near_C0_prime_ball']} < 30",
            "```",
            "",
            "## Outside Region",
            "",
            "Outside the two inner disks,",
            "",
            "```text",
            exact["outside_lower"],
            f"direct |C0| bound  = {interval['outside_C0_ball']} < 5",
            f"direct |C0'| bound = {interval['outside_C0_prime_ball']} < 60",
            "```",
            "",
            "Combining the covers proves",
            "",
            "```text",
            exact["global"],
            "```",
            "",
            "## Endpoint Collar",
            "",
            "For the critical scale,",
            "",
            "```text",
            exact["collar_map"],
            f"|Delta p| <= {interval['critical_collar_p_displacement_ball']} < 1/100",
            "```",
            "",
            "This closes the finite-dimensional derivative input in the endpoint",
            "holomorphic-lift reduction. The next step is to compose it with the",
            "published `e_A+e_B+e_C` bounds and produce explicit normalized `C1`",
            "remainder constants. The corrected first-jet lower bound remains open.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 endpoint C0 strip certificate: "
        f"{len(artifact['rows'])} rows, 4 Arb bounds, 2 removable disks, "
        "1 global strip theorem"
    )


if __name__ == "__main__":
    main()
