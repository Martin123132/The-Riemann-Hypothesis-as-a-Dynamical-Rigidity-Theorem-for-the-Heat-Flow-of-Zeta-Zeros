#!/usr/bin/env python3
"""Build an admissible strong-log-concavity countermodel at correlation order one."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
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
import sympy as sp  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.json"
)

CERTIFICATE_PRECISION_BITS = 256
CERTIFICATE_CUTOFF = 6
EXPLICIT_DELTA = sp.Rational(1, 10)
EXPLICIT_EPSILON = sp.Rational(1, 1000)
WITNESS_X = sp.Rational(21, 5)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.md"
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


def exact_algebra() -> dict:
    u, x = sp.symbols("u x", real=True)
    r = sp.symbols("r", positive=True)
    delta = sp.symbols("delta", nonnegative=True)
    epsilon = sp.symbols("epsilon", nonnegative=True)
    polynomial = 1 + u**4 / 10
    log_polynomial_second = sp.factor(sp.diff(sp.log(polynomial), u, 2))
    expected_log_polynomial_second = sp.factor(
        4 * sp.Rational(1, 10) * u**2 * (3 - u**4 / 10)
        / (1 + u**4 / 10) ** 2
    )
    if sp.simplify(log_polynomial_second - expected_log_polynomial_second) != 0:
        raise RuntimeError("polynomial log-curvature identity failed")

    curvature_constant = sp.simplify(2 - 6 / sp.sqrt(10))
    if not curvature_constant.is_positive:
        raise RuntimeError("strong-curvature constant did not prove positive")

    gaussian_transform = sp.sqrt(sp.pi) * sp.exp(-x**2 / 4)
    base_transform = sp.factor(
        gaussian_transform + sp.diff(gaussian_transform, x, 4) / 10
    )
    expected_transform = (
        sp.sqrt(sp.pi)
        * sp.exp(-x**2 / 4)
        * (x**4 - 12 * x**2 + 172)
        / 160
    )
    if sp.simplify(base_transform - expected_transform) != 0:
        raise RuntimeError("Gaussian-polynomial Fourier identity failed")

    laguerre = sp.factor(
        sp.diff(base_transform, x) ** 2
        - base_transform * sp.diff(base_transform, x, 2)
    )
    expected_laguerre = (
        sp.pi
        * sp.exp(-x**2 / 2)
        * (x**8 - 16 * x**6 + 440 * x**4 - 7680 * x**2 + 37840)
        / 51200
    )
    if sp.simplify(laguerre - expected_laguerre) != 0:
        raise RuntimeError("base Laguerre formula failed")
    witness = sp.simplify(laguerre.subs(x, 3))
    expected_witness = -sp.Rational(743, 51200) * sp.pi * sp.exp(
        -sp.Rational(9, 2)
    )
    if sp.simplify(witness - expected_witness) != 0 or not witness.is_negative:
        raise RuntimeError("negative Laguerre witness failed")

    root_log = -r - r**2 / 10 + sp.log(1 + r**2 / 10)
    root_second = sp.factor(sp.diff(root_log, r, 2))
    expected_root_second = -r**2 * (r**2 + 30) / (5 * (r**2 + 10) ** 2)
    if sp.simplify(root_second - expected_root_second) != 0:
        raise RuntimeError("root-variable concavity identity failed")

    y = sp.symbols("y", positive=True)
    theta_root_second = sp.factor(
        sp.diff(sp.cosh(4 * sp.sqrt(r)), r, 2)
    )
    expected_theta_root_second = (
        4 * sp.sqrt(r) * sp.cosh(4 * sp.sqrt(r))
        - sp.sinh(4 * sp.sqrt(r))
    ) / r ** sp.Rational(3, 2)
    if sp.simplify(theta_root_second - expected_theta_root_second) != 0:
        raise RuntimeError("theta-tail root-convexity identity failed")
    theta_numerator_derivative = sp.factor(
        sp.diff(4 * y * sp.cosh(4 * y) - sp.sinh(4 * y), y)
    )
    if theta_numerator_derivative != 16 * y * sp.sinh(4 * y):
        raise RuntimeError("theta-tail root-convexity sign identity failed")

    return {
        "family": {
            "definition": (
                "phi_(delta,epsilon)(u)=exp(-u^2-delta*u^4-"
                "epsilon*(cosh(4*u)-1))*(1+u^4/10), delta,epsilon>=0"
            ),
            "log_curvature_identity": (
                "(log phi_(delta,epsilon))''=-2-12*delta*u^2-"
                "16*epsilon*cosh(4*u)+4*(u^2/10)*(3-u^4/10)/"
                "(1+u^4/10)^2"
            ),
            "symbolic_log_polynomial_second": str(log_polynomial_second),
        },
        "strong_logconcavity": {
            "substitution": "z=u^2/sqrt(10)>=0",
            "elementary_bound": (
                "z*(3-z^2)/(1+z^2)^2<=3*z/(1+z^2)<=3/2"
            ),
            "curvature_bound": (
                "(log phi_(delta,epsilon))''<=-(2-6/sqrt(10))-"
                "12*delta*u^2-16*epsilon*cosh(4*u)"
            ),
            "uniform_constant": str(curvature_constant),
            "positivity_check": "3*sqrt(10)<10 because 90<100",
        },
        "admissibility": {
            "statement": (
                "For delta>0 or epsilon>0, phi_(delta,epsilon) is smooth, "
                "positive, even, strictly decreasing on (0,infinity), and "
                "admissible. At epsilon>0 it has theta-type double-exponential "
                "decay."
            ),
            "monotonicity": (
                "(log phi_delta)'(0)=0 and the uniform negative second-derivative "
                "bound imply (log phi_delta)'(u)<0 for u>0."
            ),
        },
        "explicit_root_concavity": {
            "delta": "1/10",
            "epsilon": "1/1000",
            "root_profile": (
                "log(phi_(1/10,1/1000)(sqrt(r)))=-r-r^2/10+"
                "log(1+r^2/10)-(cosh(4*sqrt(r))-1)/1000"
            ),
            "base_second_derivative": str(root_second),
            "theta_factor_second_derivative": str(theta_root_second),
            "theta_factor_positive_numerator_derivative": str(
                theta_numerator_derivative
            ),
            "strict_statement": (
                "d^2/dr^2 log(phi_(1/10,1/1000)(sqrt(r)))<0 for every r>0"
            ),
            "tail_bound": (
                "phi_(1/10,1/1000)(u)<=exp(1/1000)*(1+u^4/10)*"
                "exp(-u^2-u^4/10-exp(4*abs(u))/2000)"
            ),
        },
        "base_transform": {
            "definition": "F_0(x)=integral_R phi_0(u)*exp(i*x*u)du",
            "Gaussian_derivative_identity": (
                "F_0(x)=G(x)+(1/10)*G''''(x), G(x)=sqrt(pi)*exp(-x^2/4)"
            ),
            "closed_form": str(base_transform),
        },
        "negative_laguerre_witness": {
            "formula": str(laguerre),
            "point": "x=3",
            "value": str(witness),
            "decimal": format(float(witness.evalf(18)), ".17e"),
        },
        "weighted_correlation": {
            "definition": (
                "K_(1,delta,epsilon)(v)=integral_R "
                "phi_(delta,epsilon)(s+v)*phi_(delta,epsilon)(s-v)*s^2 ds"
            ),
            "source": "https://arxiv.org/abs/1309.0055",
            "identity": (
                "L_1[F_(delta,epsilon)](x)="
                "4*Fourier[K_(1,delta,epsilon)](2*x)"
            ),
            "witness": (
                "Fourier[K_(1,1/10,1/1000)](42/5)="
                "L_1[F_(1/10,1/1000)](21/5)/4<0"
            ),
            "conclusion": (
                "A smooth positive even uniformly strongly log-concave admissible "
                "kernel can have a non-positive-definite first weighted correlation."
            ),
        },
        "nonpromotion_decision": (
            "Even strict concavity of log(phi(sqrt(r))), uniform strong "
            "log-concavity, theta-type double-exponential decay, and admissibility "
            "of phi do not imply positive definiteness, let alone strict Fourier "
            "positivity, of K_(1). The Xi route must use arithmetic or modular "
            "structure beyond these generic shape and tail properties and the "
            "s^2-weighted correlation construction."
        ),
        "open_handoff": (
            "For Xi, seek a relative F_2/F_1 coercivity inequality or an "
            "s^2-weighted modular square that is absent from the countermodel."
        ),
        "symbolic_checks": {
            "curvature_constant": str(curvature_constant),
            "base_transform": str(base_transform),
            "laguerre": str(laguerre),
            "witness": str(witness),
        },
    }


def arb_text(value: flint.arb, digits: int = 55) -> str:
    return value.str(digits).replace("e", "E")


def explicit_countermodel_certificate(
    *, precision_bits: int = CERTIFICATE_PRECISION_BITS
) -> dict:
    """Certify a negative first Laguerre value at delta=1/10 and x=21/5."""
    old_precision = flint.ctx.prec
    flint.ctx.prec = precision_bits
    try:
        arb = flint.arb
        acb = flint.acb
        delta = acb(1) / 10
        epsilon = acb(1) / 1000
        witness_x = acb(21) / 5
        cutoff = acb(CERTIFICATE_CUTOFF)
        tolerance = arb("1e-55")

        def phi(u: flint.acb) -> flint.acb:
            return (
                -u * u
                - delta * u**4
                - epsilon * ((4 * u).cosh() - 1)
            ).exp() * (1 + u**4 / 10)

        values: list[flint.arb] = []
        tails: list[flint.arb] = []
        compact_values: list[flint.acb] = []
        for derivative_order in range(3):
            def integrand(u: flint.acb, analytic: bool) -> flint.acb:
                del analytic
                base = phi(u)
                if derivative_order == 0:
                    return 2 * base * (witness_x * u).cos()
                if derivative_order == 1:
                    return -2 * u * base * (witness_x * u).sin()
                return -2 * u * u * base * (witness_x * u).cos()

            compact = acb.integral(
                integrand,
                acb(0),
                cutoff,
                abs_tol=tolerance,
                rel_tol=tolerance,
                deg_limit=64,
                eval_limit=200000,
                depth_limit=64,
            )
            if not compact.imag.contains(0):
                raise RuntimeError("compact Fourier integral acquired an imaginary part")

            # For u>=6, exp(-u^2)*(1+u^4/10)<=2*u^4. Therefore the
            # full even-tail radius for derivative order j is bounded by
            # 4*integral_6^infinity u^(j+4)*exp(-u^4/10)du.
            power = 4 + derivative_order
            gamma_order = arb(power + 1) / 4
            gamma_argument = arb(CERTIFICATE_CUTOFF) ** 4 / 10
            tail = (
                arb(10) ** gamma_order
                * gamma_argument.gamma_upper(gamma_order)
            )
            value = compact.real + arb(0, tail.upper())
            compact_values.append(compact)
            tails.append(tail)
            values.append(value)

        function, first, second = values
        laguerre = first * first - function * second
        if not laguerre < 0:
            raise RuntimeError("explicit positive-delta Laguerre value was not negative")

        return {
            "arithmetic": "python-flint Acb integration and Arb tail balls",
            "precision_bits": precision_bits,
            "delta": "1/10",
            "epsilon": "1/1000",
            "x": "21/5",
            "correlation_frequency": "42/5",
            "cutoff": str(CERTIFICATE_CUTOFF),
            "absolute_tolerance": "1e-55",
            "F_ball": arb_text(function),
            "F_prime_ball": arb_text(first),
            "F_second_ball": arb_text(second),
            "laguerre_ball": arb_text(laguerre),
            "laguerre_upper": laguerre.upper().str(55, radius=False).replace("e", "E"),
            "compact_imaginary_parts": [
                arb_text(value.imag, 12) for value in compact_values
            ],
            "tail_radii": [arb_text(value, 18) for value in tails],
            "integration_tail_bound": (
                "2*integral_6^infinity u^j*phi(u)du<="
                "4*integral_6^infinity u^(j+4)*exp(-u^4/10)du"
            ),
            "theta_tail_bound": (
                "exp(-(cosh(4*u)-1)/1000)<="
                "exp(1/1000-exp(4*abs(u))/2000)"
            ),
            "passed": True,
        }
    finally:
        flint.ctx.prec = old_precision


def build_payload() -> dict:
    exact = exact_algebra()
    certificate = explicit_countermodel_certificate()
    rows = [
        GateRow(
            id="nwslc_01_countermodel_family",
            role="exact_family",
            readiness="available_exact",
            claim="A quartically damped Gaussian-polynomial family supplies admissible perturbations of an exactly transformable kernel.",
            formula=exact["family"]["definition"],
            proof_boundary="Generic family construction only; it is not the Xi kernel.",
            diagnostics=exact["family"],
        ),
        GateRow(
            id="nwslc_02_strong_curvature",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="Every member of the family has one uniform positive strong-log-concavity constant.",
            formula=exact["strong_logconcavity"]["curvature_bound"],
            proof_boundary="Exact elementary curvature estimate only.",
            diagnostics=exact["strong_logconcavity"],
        ),
        GateRow(
            id="nwslc_03_root_variable_concavity",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The explicit delta=1/10, epsilon=1/1000 kernel has the same strict root-variable log-concavity property known for Phi.",
            formula=exact["explicit_root_concavity"]["strict_statement"],
            proof_boundary="Exact generic shape theorem; it is not an arithmetic identity for Xi.",
            diagnostics=exact["explicit_root_concavity"],
        ),
        GateRow(
            id="nwslc_04_explicit_admissible_kernel",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The explicit delta=1/10, epsilon=1/1000 member is admissible with theta-type double-exponential decay.",
            formula=exact["admissibility"]["statement"],
            proof_boundary="Kernel-class membership only; no transform sign follows.",
            diagnostics=exact["admissibility"],
        ),
        GateRow(
            id="nwslc_05_base_fourier_transform",
            role="exact_identity",
            readiness="available_exact",
            claim="The delta-zero Fourier transform is an explicit Gaussian times quartic polynomial.",
            formula=exact["base_transform"]["closed_form"],
            proof_boundary="Exact Gaussian differentiation identity at delta=0.",
            diagnostics=exact["base_transform"],
        ),
        GateRow(
            id="nwslc_06_endpoint_laguerre_witness",
            role="exact_counterexample",
            readiness="guard_validated",
            claim="The first Laguerre expression is strictly negative at one exact frequency.",
            formula=exact["negative_laguerre_witness"]["value"],
            proof_boundary="Exact delta-zero benchmark; the explicit admissible witness is independently interval-certified.",
            diagnostics=exact["negative_laguerre_witness"],
        ),
        GateRow(
            id="nwslc_07_explicit_arb_witness",
            role="interval_certificate",
            readiness="interval_validated",
            claim="Acb quadrature with explicit Arb tails certifies a negative first Laguerre value for the root-concave theta-tail admissible kernel.",
            formula=certificate["laguerre_ball"],
            proof_boundary="One rigorous scalar witness at delta=1/10 and x=21/5.",
            diagnostics=certificate,
        ),
        GateRow(
            id="nwslc_08_weighted_correlation_failure",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Even root-variable concavity plus strong log-concavity does not make the first s^2-weighted correlation positive definite.",
            formula=exact["weighted_correlation"]["witness"],
            proof_boundary="Explicit generic admissible countermodel; it does not decide the Xi first correlation.",
            diagnostics={
                "correlation": exact["weighted_correlation"],
                "certificate": certificate,
            },
        ),
        GateRow(
            id="nwslc_09_xi_specific_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The surviving Xi route must use structure absent from generic weighted correlations.",
            formula=exact["open_handoff"],
            proof_boundary="Open Xi-specific coercivity or modular-square target, not RH or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate",
        "date": "2026-07-17",
        "status": (
            "explicit root-concave theta-tail admissible weighted-correlation "
            "countermodel to generic shape-and-tail promotion"
        ),
        "proof_boundary": (
            "This artifact gives an explicit smooth positive even admissible kernel "
            "whose log(phi(sqrt(r))) is strictly concave, whose ordinary log is "
            "uniformly strongly concave, whose tail is theta-type double "
            "exponential, and whose first s^2-weighted correlation has a "
            "rigorously negative Fourier value. The sign uses Acb quadrature with "
            "explicit Arb tails. It closes generic shape-tail-plus-weighting "
            "promotion but does not determine the Xi correlation, prove Wiener "
            "density, Lambda<=0, or RH."
        ),
        "exact": exact,
        "certificate": certificate,
        "rows": [asdict(row) for row in rows],
        "sources": [
            "https://arxiv.org/abs/1309.0055",
            "outputs/jensen_window_pf_newman_positive_time_strong_logconcavity_gate.md",
            "outputs/jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.md",
            "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md",
        ],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    certificate = payload["certificate"]
    return f"""# Jensen-Window PF Newman Weighted Strong-Log-Concavity Countermodel Gate

Date: {payload['date']}

Status: {payload['status']}. This is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.json
python work/rh_compute/scripts/jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_weighted_strong_logconcavity_countermodel_gate.py
```

Current result:

```text
validated Jensen-window PF Newman weighted strong-log-concavity countermodel gate: 9 rows, 0 issues, 1 exact strong-curvature bound, 1 exact root-variable concavity theorem, 1 explicit theta-tail admissible kernel, 1 Gaussian endpoint identity, 1 exact endpoint witness, 1 Arb theta-tail witness, 1 weighted-correlation countermodel, 1 Xi-specific handoff
```

## Countermodel Family

```text
{exact['family']['definition']}
```

For `z=u^2/sqrt(10)>=0`,

```text
{exact['strong_logconcavity']['elementary_bound']}
{exact['strong_logconcavity']['curvature_bound']}
2-6/sqrt(10)>0
```

Thus every `delta>=0` member is uniformly strongly log-concave. For every
`delta>0`, quartic damping makes the kernel admissible and
super-Gaussian; evenness and the curvature bound make it strictly decreasing
on the positive half-line.

At the explicit values `delta=1/10`, `epsilon=1/1000`, it also has the
published Xi kernel's stronger root-variable shape and tail properties:

```text
{exact['explicit_root_concavity']['root_profile']}
base root curvature={exact['explicit_root_concavity']['base_second_derivative']}<0
d^2/dr^2 cosh(4*sqrt(r))={exact['explicit_root_concavity']['theta_factor_second_derivative']}>0
{exact['explicit_root_concavity']['tail_bound']}
```

## Exact Endpoint Witness

At `delta=0`, Gaussian differentiation gives

```text
{exact['base_transform']['closed_form']}
L_1[F_0](x)={exact['negative_laguerre_witness']['formula']}
L_1[F_0](3)={exact['negative_laguerre_witness']['value']}<0
```

This endpoint identity is a cross-check. The admissible root-concave member is
certified directly at `delta=1/10`, `epsilon=1/1000`, `x=21/5`. Acb
integration on `[0,6]` and upper incomplete-gamma tail balls give

```text
F={certificate['F_ball']}
F'={certificate['F_prime_ball']}
F''={certificate['F_second_ball']}
L_1[F]={certificate['laguerre_ball']}<0
```

## Weighted-Correlation Failure

For the first Csordas correlation,

```text
{exact['weighted_correlation']['definition']}
{exact['weighted_correlation']['identity']}
{exact['weighted_correlation']['witness']}
```

Hence a smooth positive even admissible kernel with uniform strong
log-concavity, strict root-variable log-concavity, and the Xi theta tail class
can have a first `s^2`-weighted correlation that is not positive definite.
This is stronger than a zeroth-correlation shape guard: even the stronger
Xi-style shape and tail properties plus the generic weighted construction do
not supply the missing sign.

Primary source for the correlation identity:
https://arxiv.org/abs/1309.0055.

## Nonpromotion Gate

{exact['nonpromotion_decision']}

## Live Handoff

{exact['open_handoff']}
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Jensen-window PF Newman weighted strong-log-concavity "
        "countermodel gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
