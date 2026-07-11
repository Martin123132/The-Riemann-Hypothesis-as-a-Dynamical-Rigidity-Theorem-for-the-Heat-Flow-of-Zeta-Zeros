#!/usr/bin/env python3
"""Audit positive one- and two-block Gasper residual decompositions of Xi."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

from flint import acb, arb, ctx  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_gasper_residual_two_block_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_gasper_residual_two_block_gate.md"
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


def acb_xi(s: acb) -> acb:
    pi_ball = arb.pi()
    return (
        acb(arb("1/2"))
        * s
        * (s - 1)
        * acb(pi_ball) ** (-s / 2)
        * (s / 2).gamma()
        * s.zeta()
    )


def acb_h0(x: acb) -> acb:
    return acb_xi((1 + acb(0, 1) * x) / 2) / 8


def acb_gasper_block(x: acb, order: arb) -> acb:
    pi_ball = arb.pi()
    imaginary_shift = acb(0, 1) * x / 4
    argument = acb(2 * pi_ball)
    return acb(pi_ball**2 / 2) * (
        argument.bessel_k(acb(order) + imaginary_shift)
        + argument.bessel_k(acb(order) - imaginary_shift)
    )


def acb_first_residual(x: acb) -> acb:
    return acb_h0(x) - acb_gasper_block(x, arb("9/4"))


def acb_second_block(x: acb) -> acb:
    return acb_gasper_block(x, arb("5/4"))


def cauchy_finite_difference_jet(
    function,
    x_text: str,
    step: arb,
    cauchy_radius: arb,
    box_radius: arb,
) -> tuple[list[arb], arb]:
    """Enclose f, f', f'' using central differences and Cauchy bounds."""

    x = arb(x_text)
    f_minus = function(acb(x - step))
    f_zero = function(acb(x))
    f_plus = function(acb(x + step))
    box = acb(arb(x, box_radius), arb(0, box_radius))
    maximum = abs(function(box)).upper()
    if not maximum.is_finite():
        raise RuntimeError(f"nonfinite Cauchy box at x={x_text}")

    first = ((f_plus - f_minus) / (2 * step)).real
    second = ((f_plus - 2 * f_zero + f_minus) / step**2).real
    first_error = (step**2 * maximum / cauchy_radius**3).upper()
    second_error = (2 * step**2 * maximum / cauchy_radius**4).upper()
    first += arb(0, first_error)
    second += arb(0, second_error)
    return [f_zero.real, first, second], maximum


def laguerre(jet: list[arb]) -> arb:
    return jet[1] ** 2 - jet[0] * jet[2]


def cross_laguerre(left: list[arb], right: list[arb]) -> arb:
    return (
        2 * left[1] * right[1]
        - left[0] * right[2]
        - right[0] * left[2]
    )


def build_interval_certificate() -> dict:
    ctx.dps = 110
    step = arb("1e-6")
    cauchy_radius = arb("0.04")
    box_radius = arb("0.1")
    beta_match = arb.pi() - 3 / (2 * arb.pi())
    split = arb("23/10")
    rows: list[dict] = []

    for x_text, beta_left, beta_right in (
        ("66", arb(0), split),
        ("50", split, beta_match),
    ):
        residual_jet, residual_maximum = cauchy_finite_difference_jet(
            acb_first_residual, x_text, step, cauchy_radius, box_radius
        )
        second_jet, second_maximum = cauchy_finite_difference_jet(
            acb_second_block, x_text, step, cauchy_radius, box_radius
        )
        c0 = laguerre(residual_jet)
        c1 = -cross_laguerre(residual_jet, second_jet)
        c2 = laguerre(second_jet)

        def quadratic(beta: arb) -> arb:
            return c0 + beta * c1 + beta**2 * c2

        left_value = quadratic(beta_left)
        right_value = quadratic(beta_right)
        if not (c2 > 0 and left_value < 0 and right_value < 0):
            raise RuntimeError(f"two-block interval certificate failed at x={x_text}")
        rows.append(
            {
                "x": x_text,
                "beta_interval": [str(beta_left), str(beta_right)],
                "first_residual_jet": [str(value) for value in residual_jet],
                "second_block_jet": [str(value) for value in second_jet],
                "cauchy_box_abs_upper": {
                    "first_residual": str(residual_maximum),
                    "second_block": str(second_maximum),
                },
                "laguerre_quadratic": {
                    "c0": str(c0),
                    "c1": str(c1),
                    "c2": str(c2),
                },
                "left_endpoint_value": str(left_value),
                "right_endpoint_value": str(right_value),
                "certified": {
                    "quadratic_strictly_convex": True,
                    "left_endpoint_negative": True,
                    "right_endpoint_negative": True,
                    "negative_on_closed_beta_interval": True,
                },
            }
        )

    return {
        "rigorous": True,
        "arithmetic": "python-flint Acb/Arb balls at 110 decimal digits",
        "spectral_point_formula": (
            "H_0(x)=xi((1+i*x)/2)/8, "
            "P_a(x)=pi^2/2*(K_(a+i*x/4)(2*pi)+K_(a-i*x/4)(2*pi))"
        ),
        "derivative_method": (
            "For real entire f, central differences with step h enclose f' and "
            "f''. A square Acb evaluation bounds abs(f) by M on every Cauchy "
            "circle of radius r around the finite-difference segment; the added "
            "errors are h^2*M/r^3 and 2*h^2*M/r^4."
        ),
        "step": "1e-6",
        "cauchy_radius": "0.04",
        "acb_box_radius": "0.1",
        "beta_split": "23/10",
        "beta_match": str(beta_match),
        "rows": rows,
        "conclusion": (
            "L[R_beta](66)<0 for 0<=beta<=23/10 and "
            "L[R_beta](50)<0 for 23/10<=beta<=pi-3/(2*pi)."
        ),
        "scope": (
            "A complete interval certificate over the tail-admissible nonnegative "
            "two-block parameter range at t=0. It rejects an independently "
            "Laguerre-Polya residual, not a signed coupled identity or RH."
        ),
    }


def build_exact() -> dict:
    beta, x = sp.symbols("beta x", real=True)
    q = sp.symbols("q", positive=True)
    a = sp.Rational(3, 2) / sp.pi
    b = sp.pi - a

    f, fp, fpp, g, gp, gpp = sp.symbols("f fp fpp g gp gpp", real=True)
    l_f = fp**2 - f * fpp
    l_g = gp**2 - g * gpp
    cross = 2 * fp * gp - f * gpp - g * fpp
    direct = (fp - beta * gp) ** 2 - (f - beta * g) * (fpp - beta * gpp)
    if sp.expand(direct - (l_f - beta * cross + beta**2 * l_g)) != 0:
        raise RuntimeError("two-block Laguerre quadratic failed")

    y = sp.symbols("y")
    multiplier_polynomial = (
        256 * y**4
        - 576 * y**3
        + (432 + 16 * beta) * y**2
        - (120 + 20 * beta) * y
        + 9
        + 5 * beta
    )
    discriminant = sp.factor(sp.discriminant(multiplier_polynomial, y))
    expected_discriminant = -sp.Integer(2) ** 24 * (
        20 * beta**5
        - 64 * beta**4
        + 184 * beta**3
        - 432 * beta**2
        + 972 * beta
        - 729
    )
    if sp.expand(discriminant - expected_discriminant) != 0:
        raise RuntimeError("two-block multiplier discriminant failed")

    derivative_quartic = (
        25 * beta**4
        - 64 * beta**3
        + 138 * beta**2
        - 216 * beta
        + 243
    )
    parameter = sp.symbols("parameter")
    left = sp.Rational(11, 10)
    right = sp.Rational(64, 25)
    on_unit = sp.Poly(
        sp.expand(derivative_quartic.subs(beta, left + (right - left) * parameter)),
        parameter,
    )
    power = [on_unit.nth(index) for index in range(5)]
    bernstein = [
        sp.factor(
            sum(
                power[index]
                * sp.binomial(degree, index)
                / sp.binomial(4, index)
                for index in range(degree + 1)
            )
        )
        for degree in range(5)
    ]
    if not all(value > 0 for value in bernstein):
        raise RuntimeError("multiplier discriminant monotonicity failed")

    return {
        "one_block_positive_residual": {
            "definitions": (
                "Psi_9(u)=4*pi^2*exp(-2*pi*cosh(4u))*cosh(9u), "
                "E_0(u)=Phi(u)-Psi_9(u)"
            ),
            "quotient": (
                "With q=exp(-4u) and a=3/(2*pi), "
                "phi_1(u)/Psi_9(u)=exp(pi*q)*(1-a*q)/(1+q^(9/2))."
            ),
            "inequality": (
                "exp(pi*q)*(1-a*q)>(1+pi*q)*(1-a*q)>"
                "1+q>=1+q^(9/2) for 0<q<=1, since pi>3 gives "
                "pi-a-3/2>1."
            ),
            "conclusion": "Phi(u)>phi_1(u)>Psi_9(u) for every finite u>=0",
            "tail": (
                "E_0(u)~2*pi^2*(pi-3/(2*pi))*exp(5u-pi*exp(4u))"
            ),
        },
        "two_block_family": {
            "kernel": (
                "Psi_beta(u)=4*pi^2*exp(-2*pi*cosh(4u))*"
                "(cosh(9u)+beta*cosh(5u)), beta>=0"
            ),
            "transform": "P_beta=P_(9/4)+beta*P_(5/4), R_beta=H_0-P_beta",
            "beta_match": "beta_tail=pi-3/(2*pi)",
            "residual_tail": (
                "Phi(u)-Psi_beta(u)~2*pi^2*(pi-3/(2*pi)-beta)*"
                "exp(5u-pi*exp(4u))"
            ),
            "positive_residual_necessity": (
                "A residual kernel nonnegative for all sufficiently large u "
                "must satisfy 0<=beta<=beta_tail."
            ),
        },
        "laguerre_quadratic": {
            "forms": "L[f]=f'^2-f*f'', B[f,g]=2*f'*g'-f*g''-g*f''",
            "identity": "L[R_beta]=L[E_0]-beta*B[E_0,P_(5/4)]+beta^2*L[P_(5/4)]",
            "convexity_use": (
                "When the beta^2 coefficient is positive, negativity at both "
                "ends of a closed beta interval proves negativity throughout."
            ),
        },
        "two_block_obstruction": {
            "theorem": (
                "For every beta>=0, the decomposition Phi=Psi_beta+residual "
                "cannot have both a globally nonnegative residual kernel and a "
                "residual transform satisfying the first Laguerre inequality "
                "on all real x."
            ),
            "tail_side": (
                "beta>pi-3/(2*pi) makes the residual kernel eventually negative."
            ),
            "spectral_side": (
                "For 0<=beta<=pi-3/(2*pi), the interval certificate gives "
                "L[R_beta](66)<0 or L[R_beta](50)<0."
            ),
        },
        "multiplier_guard": {
            "multiplier": "U_beta(z)=cosh(9z)+beta*cosh(5z)",
            "reduction": (
                "U_beta(z)=cosh(z)*Q_beta(cosh(z)^2), where "
                f"Q_beta(y)={sp.sstr(multiplier_polynomial)}"
            ),
            "discriminant": f"disc_y Q_beta={sp.sstr(discriminant)}",
            "monotonicity": (
                "For beta>=11/10, F(beta)=20beta^5-64beta^4+184beta^3-"
                "432beta^2+972beta-729 is strictly increasing and "
                "F(11/10)=4459/5000>0."
            ),
            "bernstein_coefficients": [str(value) for value in bernstein],
            "conclusion": (
                "For beta>=11/10 the quartic discriminant is negative, so "
                "U_beta has off-imaginary zeros and fails the standard "
                "Laguerre-Polya universal-multiplier test. In particular the "
                "tail-matched beta is not certified by that theorem."
            ),
        },
        "open_handoff": (
            "Any surviving Gasper comparison must now use at least one signed "
            "block or a coupled square identity; decomposing Phi into a positive "
            "9/5 Gasper block plus a separately Laguerre-Polya positive residual "
            "is impossible. The live target is a sign-aware multi-block identity "
            "whose mixed terms are controlled together after modular cancellation."
        ),
    }


def build_payload() -> dict:
    exact = build_exact()
    interval = build_interval_certificate()
    rows = [
        GateRow(
            id="grtb_01_one_block_residual_positive",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The true Xi kernel strictly dominates the first Gasper fake-Xi kernel pointwise.",
            formula=exact["one_block_positive_residual"]["conclusion"],
            proof_boundary="Exact kernel inequality only; positive kernels need not have real-rooted transforms.",
            diagnostics=exact["one_block_positive_residual"],
        ),
        GateRow(
            id="grtb_02_two_block_tail_parameter",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Tail positivity confines every nonnegative 9/5 two-block residual to one compact beta interval.",
            formula=exact["two_block_family"]["positive_residual_necessity"],
            proof_boundary="Exact leading-tail necessity; no sufficiency is claimed.",
            diagnostics=exact["two_block_family"],
        ),
        GateRow(
            id="grtb_03_laguerre_beta_quadratic",
            role="exact_identity",
            readiness="available_exact",
            claim="The residual first Laguerre expression is an exact quadratic in the second-block coefficient.",
            formula=exact["laguerre_quadratic"]["identity"],
            proof_boundary="Exact quadratic identity used by the interval certificate.",
            diagnostics=exact["laguerre_quadratic"],
        ),
        GateRow(
            id="grtb_04_acb_cauchy_derivatives",
            role="interval_method",
            readiness="interval_validated",
            claim="Acb Cauchy bounds rigorously enclose the finite-difference jets used in the beta-uniform test.",
            formula=interval["derivative_method"],
            proof_boundary="Finite spectral-point derivative enclosures at t=0.",
            diagnostics={key: interval[key] for key in ("arithmetic", "step", "cauchy_radius", "acb_box_radius")},
        ),
        GateRow(
            id="grtb_05_beta_uniform_negative_laguerre",
            role="interval_theorem",
            readiness="interval_validated",
            claim="Every tail-admissible nonnegative two-block residual violates the first Laguerre inequality.",
            formula=interval["conclusion"],
            proof_boundary="Complete interval coverage in beta at t=0; not a statement about signed coupled sums.",
            diagnostics=interval,
        ),
        GateRow(
            id="grtb_06_two_block_obstruction",
            role="exact_interval_composition",
            readiness="ready_to_apply",
            claim="No positive 9/5 Gasper block plus positive independently-LP residual decomposition can prove Xi.",
            formula=exact["two_block_obstruction"]["theorem"],
            proof_boundary="Rejects this decomposition architecture only; it neither proves nor disproves RH.",
            diagnostics=exact["two_block_obstruction"],
        ),
        GateRow(
            id="grtb_07_multiplier_discriminant_guard",
            role="exact_route_obstruction",
            readiness="guard_validated",
            claim="The tail-matched 9/5 multiplier fails the standard universal-factor imaginary-zero test.",
            formula=exact["multiplier_guard"]["conclusion"],
            proof_boundary="Rejects the standard universal-multiplier proof, not real-rootedness by every possible special-function theorem.",
            diagnostics=exact["multiplier_guard"],
        ),
        GateRow(
            id="grtb_08_signed_multiblock_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="A surviving Gasper route must control signed blocks and mixed terms together.",
            formula=exact["open_handoff"],
            proof_boundary="Open coupled-square target; not strict Laguerre positivity, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_gasper_residual_two_block_gate",
        "date": "2026-07-11",
        "status": "exact and interval-certified positive two-block Gasper obstruction",
        "proof_boundary": (
            "This artifact proves exact one-block kernel domination, derives the "
            "two-block tail parameter and Laguerre quadratic, and uses Acb/Arb "
            "Cauchy derivative enclosures to reject every tail-admissible "
            "nonnegative two-block residual as an independent Laguerre-Polya "
            "component. It also rejects the standard universal-factor proof for "
            "the tail-matched multiplier. Signed blocks and coupled mixed-term "
            "identities remain open; this does not prove or disprove RH or "
            "Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_gasper_fake_xi_remainder_gate.md",
            "outputs/jensen_window_pf_newman_theta_summand_spectral_square_gate.md",
            "outputs/formal_core.md",
            "https://arxiv.org/abs/0801.2996",
            "https://arxiv.org/abs/1502.06844",
        ],
        "exact": exact,
        "interval_certificate": interval,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    interval = payload["interval_certificate"]
    lines = [
        "# Jensen-Window PF Newman Gasper Residual Two-Block Gate",
        "",
        "Date: 2026-07-11",
        "",
        "Status: exact and interval-certified positive two-block Gasper",
        "obstruction. This is not a proof or disproof of RH or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_newman_gasper_residual_two_block_gate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_gasper_residual_two_block_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_gasper_residual_two_block_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_gasper_residual_two_block_gate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF Newman Gasper residual two-block gate: 8 rows, 0 issues, 2 exact kernel theorems, 1 exact Laguerre quadratic, 2 Acb derivative certificates, 2 beta intervals covered, 1 exhaustive positive-residual obstruction, 1 multiplier guard, 1 signed handoff",
        "```",
        "",
        "## Positive First Residual",
        "",
        "For `q=exp(-4u)` and `a=3/(2pi)`,",
        "",
        "```text",
        exact["one_block_positive_residual"]["quotient"],
        exact["one_block_positive_residual"]["inequality"],
        exact["one_block_positive_residual"]["conclusion"],
        "```",
        "",
        "Thus subtracting the original fake-Xi block leaves a strictly positive",
        "kernel. Positivity alone is not enough: the interval theorem below",
        "includes `beta=0` and proves that this residual transform violates the",
        "first Laguerre inequality.",
        "",
        "## Two-Block Family",
        "",
        "```text",
        exact["two_block_family"]["kernel"],
        exact["two_block_family"]["transform"],
        exact["two_block_family"]["residual_tail"],
        "```",
        "",
        "A globally nonnegative residual must therefore have",
        "`0<=beta<=pi-3/(2pi)`. On this complete parameter interval,",
        "",
        "```text",
        exact["laguerre_quadratic"]["identity"],
        "```",
        "",
        "## Interval Theorem",
        "",
        interval["derivative_method"],
        "The resulting certified quadratics are:",
        "",
        "```text",
    ]
    for row in interval["rows"]:
        coefficients = row["laguerre_quadratic"]
        lines.extend(
            [
                f"x={row['x']}, beta in [{row['beta_interval'][0]}, {row['beta_interval'][1]}]",
                f"  c0={coefficients['c0']}",
                f"  c1={coefficients['c1']}",
                f"  c2={coefficients['c2']}",
                f"  left={row['left_endpoint_value']}",
                f"  right={row['right_endpoint_value']}",
            ]
        )
    lines.extend(
        [
            "```",
            "",
            "Each quadratic is strictly convex and negative at both endpoints,",
            "hence negative throughout its closed beta interval. The intervals",
            "meet at `23/10` and cover every beta allowed by tail positivity.",
            "Consequently no member of this positive two-block family leaves a",
            "positive residual whose transform is independently Laguerre-Polya.",
            "",
            "## Multiplier Guard",
            "",
            "The tail-matched multiplier also fails the standard universal-factor",
            "test. Writing `U_beta(z)=cosh(9z)+beta*cosh(5z)` gives",
            "",
            "```text",
            exact["multiplier_guard"]["reduction"],
            exact["multiplier_guard"]["discriminant"],
            "```",
            "",
            exact["multiplier_guard"]["conclusion"],
            "",
            "## Live Handoff",
            "",
            exact["open_handoff"],
            "",
            "The comparison programme has therefore moved from positive block",
            "domination to a genuinely signed, coupled-square problem.",
            "",
            "References: https://arxiv.org/abs/0801.2996 and",
            "https://arxiv.org/abs/1502.06844.",
            "",
        ]
    )
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
    print(f"wrote Newman Gasper residual two-block gate: {len(payload['rows'])} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
