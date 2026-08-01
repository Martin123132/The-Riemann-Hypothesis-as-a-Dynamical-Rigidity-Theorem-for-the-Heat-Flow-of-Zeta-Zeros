#!/usr/bin/env python3
"""Certify a real-zero band and its degree-361 Jensen consequence."""

from __future__ import annotations

import argparse
import ctypes
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402
from flint import acb, arb  # noqa: E402


PRECISION_BITS = 192
TIME_ORDER = 2
IMAGINARY_ORDER = 4
TIME_LOWER = Fraction(0)
TIME_UPPER = Fraction(1, 5)
TIME_STEP = Fraction(1, 50)
IMAGINARY_LOWER = Fraction(0)
IMAGINARY_UPPER = Fraction(1)
IMAGINARY_STEP = Fraction(1, 20)
REAL_BOUNDARY = Fraction(38)
INTEGRATION_CUTOFF = Fraction(2)
RETAINED_THETA_TERMS = 4
FAR_TAIL_RADIUS = "1e-800"
ABS_TOL = "1e-45"

DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_real_zero_band_"
    "degree361_sector_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_newman_real_zero_band_"
    "degree361_sector_certificate.md"
)
CONTACT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_"
    "compact_transversality_interval_certificate.json"
)
DEGREE71_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_zero_slab_degree71_sector_certificate.json"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def request_below_normal_priority() -> bool:
    if os.name != "nt":
        return False
    try:
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        kernel32.SetPriorityClass.argtypes = (
            ctypes.c_void_p,
            ctypes.c_uint32,
        )
        kernel32.SetPriorityClass.restype = ctypes.c_int
        return bool(
            kernel32.SetPriorityClass(
                kernel32.GetCurrentProcess(),
                0x00004000,
            )
        )
    except (AttributeError, OSError):
        return False


def fraction_decimal(value: Fraction) -> str:
    getcontext().prec = 100
    result = Decimal(value.numerator) / Decimal(value.denominator)
    return format(result, "f")


def interval_ball(center: Fraction, radius: Fraction) -> arb:
    return arb(
        f"[{fraction_decimal(center)} +/- {fraction_decimal(radius)}]"
    )


def symmetric_ball(radius: Fraction) -> arb:
    return arb(f"[0 +/- {fraction_decimal(radius)}]")


def arb_power(value: arb, exponent: int) -> arb:
    result = arb(1)
    for _ in range(exponent):
        result *= value
    return result


def fraction_range(
    lower: Fraction,
    upper: Fraction,
    step: Fraction,
) -> list[tuple[Fraction, Fraction]]:
    result: list[tuple[Fraction, Fraction]] = []
    left = lower
    while left < upper:
        right = min(left + step, upper)
        result.append((left, right))
        left = right
    return result


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interval_record(value: arb, digits: int = 70) -> dict[str, object]:
    return {
        "ball": value.str(n=digits, more=True),
        "lower": value.lower().str(n=digits, more=True),
        "upper": value.upper().str(n=digits, more=True),
        "strictly_negative": bool(value.upper() < 0),
    }


class BoundaryCertifier:
    def __init__(self) -> None:
        flint.ctx.prec = PRECISION_BITS
        self.pi = acb.pi()
        self.zero = acb(0)
        self.cutoff = acb(fraction_decimal(INTEGRATION_CUTOFF))
        self.real_boundary = acb(fraction_decimal(REAL_BOUNDARY))
        self.abs_tol = arb(ABS_TOL)
        self.far_tail = arb(f"[0 +/- {FAR_TAIL_RADIUS}]")
        self.moment_cache: dict[
            tuple[int, Fraction, Fraction],
            arb,
        ] = {}
        self.coefficient_cache: dict[
            tuple[int, int, Fraction, Fraction],
            arb,
        ] = {}

    def phi_n(self, u: acb, n: int) -> acb:
        n_ball = acb(n)
        return (
            2 * self.pi * self.pi * n_ball**4 * (9 * u).exp()
            - 3 * self.pi * n_ball**2 * (5 * u).exp()
        ) * (-self.pi * n_ball**2 * (4 * u).exp()).exp()

    def retained_phi(self, u: acb) -> acb:
        return sum(
            (
                self.phi_n(u, n)
                for n in range(1, RETAINED_THETA_TERMS + 1)
            ),
            acb(0),
        )

    def omitted_theta_envelope(self, u: acb) -> acb:
        first = RETAINED_THETA_TERMS + 1
        ratio = (
            acb(first + 1) / first
        ) ** 4 * (
            -self.pi * (2 * first + 1)
        ).exp()
        return (
            2
            * self.pi
            * self.pi
            * first**4
            * (9 * u).exp()
            * (-self.pi * first**2 * (4 * u).exp()).exp()
            / (1 - ratio)
        )

    def integrate_real(self, func) -> arb:
        value = acb.integral(
            func,
            self.zero,
            self.cutoff,
            abs_tol=self.abs_tol,
            rel_tol=self.abs_tol,
            deg_limit=80,
            eval_limit=300000,
            depth_limit=80,
        )
        if not value.imag.contains(0):
            raise RuntimeError("directed real integral gained imaginary part")
        return value.real

    def finite_moment(
        self,
        moment: int,
        time: Fraction,
        imaginary: Fraction,
    ) -> arb:
        time_ball = acb(fraction_decimal(time))
        imaginary_ball = acb(fraction_decimal(imaginary))

        def retained(u: acb, analytic: bool) -> acb:
            del analytic
            return (
                u**moment
                * (time_ball * u * u).exp()
                * self.retained_phi(u)
                * (imaginary_ball * u).cosh()
            )

        def omitted(u: acb, analytic: bool) -> acb:
            del analytic
            return (
                u**moment
                * (time_ball * u * u).exp()
                * self.omitted_theta_envelope(u)
                * (imaginary_ball * u).cosh()
            )

        return (
            self.integrate_real(retained)
            + self.integrate_real(omitted)
            + self.far_tail
        )

    def moment_upper(
        self,
        moment: int,
        time: Fraction,
        imaginary: Fraction,
    ) -> arb:
        key = (moment, time, imaginary)
        cached = self.moment_cache.get(key)
        if cached is not None:
            return cached
        value = self.finite_moment(moment, time, imaginary)
        self.moment_cache[key] = value
        return value

    def coefficient(
        self,
        time_order: int,
        imaginary_order: int,
        time: Fraction,
        imaginary: Fraction,
    ) -> arb:
        key = (
            time_order,
            imaginary_order,
            time,
            imaginary,
        )
        cached = self.coefficient_cache.get(key)
        if cached is not None:
            return cached
        moment = 2 * time_order + imaginary_order
        time_ball = acb(fraction_decimal(time))
        imaginary_ball = acb(fraction_decimal(imaginary))
        use_cosh = imaginary_order % 2 == 0

        def retained(u: acb, analytic: bool) -> acb:
            del analytic
            hyperbolic = (
                (imaginary_ball * u).cosh()
                if use_cosh
                else (imaginary_ball * u).sinh()
            )
            return (
                u**moment
                * (time_ball * u * u).exp()
                * self.retained_phi(u)
                * (self.real_boundary * u).cos()
                * hyperbolic
            )

        value = self.integrate_real(retained)
        omitted = self.moment_upper(
            moment,
            time,
            imaginary,
        ) - self.integrate_retained_moment(
            moment,
            time,
            imaginary,
        )
        radius = omitted.abs_upper()
        result = value + arb(0, str(radius))
        self.coefficient_cache[key] = result
        return result

    def integrate_retained_moment(
        self,
        moment: int,
        time: Fraction,
        imaginary: Fraction,
    ) -> arb:
        time_ball = acb(fraction_decimal(time))
        imaginary_ball = acb(fraction_decimal(imaginary))

        def retained(u: acb, analytic: bool) -> acb:
            del analytic
            return (
                u**moment
                * (time_ball * u * u).exp()
                * self.retained_phi(u)
                * (imaginary_ball * u).cosh()
            )

        return self.integrate_real(retained)

    def certify_box(
        self,
        time_low: Fraction,
        time_high: Fraction,
        imaginary_low: Fraction,
        imaginary_high: Fraction,
    ) -> dict[str, object]:
        time_center = (time_low + time_high) / 2
        time_radius = (time_high - time_low) / 2
        imaginary_center = (
            imaginary_low + imaginary_high
        ) / 2
        imaginary_radius = (
            imaginary_high - imaginary_low
        ) / 2
        tau = symmetric_ball(time_radius)
        eta = symmetric_ball(imaginary_radius)
        polynomial = arb(0)
        for time_order in range(TIME_ORDER + 1):
            for imaginary_order in range(IMAGINARY_ORDER + 1):
                polynomial += (
                    self.coefficient(
                        time_order,
                        imaginary_order,
                        time_center,
                        imaginary_center,
                    )
                    * arb_power(tau, time_order)
                    * arb_power(eta, imaginary_order)
                    / (
                        math.factorial(time_order)
                        * math.factorial(imaginary_order)
                    )
                )

        time_remainder = (
            arb_power(
                arb(fraction_decimal(time_radius)),
                TIME_ORDER + 1,
            )
            / math.factorial(TIME_ORDER + 1)
            * self.moment_upper(
                2 * (TIME_ORDER + 1),
                time_high,
                imaginary_high,
            ).upper()
        )
        imaginary_remainder = arb(0)
        for time_order in range(TIME_ORDER + 1):
            imaginary_remainder += (
                arb_power(
                    arb(fraction_decimal(time_radius)),
                    time_order,
                )
                / math.factorial(time_order)
                * arb_power(
                    arb(fraction_decimal(imaginary_radius)),
                    IMAGINARY_ORDER + 1,
                )
                / math.factorial(IMAGINARY_ORDER + 1)
                * self.moment_upper(
                    2 * time_order + IMAGINARY_ORDER + 1,
                    time_center,
                    imaginary_high,
                ).upper()
            )
        remainder = (
            time_remainder + imaginary_remainder
        ).upper()
        enclosure = polynomial + arb(0, str(remainder))
        if not enclosure.upper() < 0:
            raise RuntimeError(
                "unresolved boundary box "
                f"t=[{time_low},{time_high}], "
                f"y=[{imaginary_low},{imaginary_high}]: "
                f"{enclosure}"
            )
        margin = -enclosure.upper()
        return {
            "t_low": str(time_low),
            "t_high": str(time_high),
            "y_low": str(imaginary_low),
            "y_high": str(imaginary_high),
            "enclosure": enclosure.str(n=70, more=True),
            "upper": enclosure.upper().str(n=70, more=True),
            "negative_margin_lower": margin.str(n=70, more=True),
            "certified": True,
        }

    def build(self) -> dict[str, object]:
        records: list[dict[str, object]] = []
        for time_low, time_high in fraction_range(
            TIME_LOWER,
            TIME_UPPER,
            TIME_STEP,
        ):
            for imaginary_low, imaginary_high in fraction_range(
                IMAGINARY_LOWER,
                IMAGINARY_UPPER,
                IMAGINARY_STEP,
            ):
                records.append(
                    self.certify_box(
                        time_low,
                        time_high,
                        imaginary_low,
                        imaginary_high,
                    )
                )
        weakest = min(
            records,
            key=lambda row: arb(row["negative_margin_lower"]),
        )
        return {
            "precision_bits": PRECISION_BITS,
            "time_taylor_order": TIME_ORDER,
            "imaginary_taylor_order": IMAGINARY_ORDER,
            "time_step": str(TIME_STEP),
            "imaginary_step": str(IMAGINARY_STEP),
            "integration_cutoff": str(INTEGRATION_CUTOFF),
            "retained_theta_terms": RETAINED_THETA_TERMS,
            "far_tail_radius": FAR_TAIL_RADIUS,
            "boxes": len(records),
            "unresolved_boxes": 0,
            "weakest_negative_margin_lower": weakest[
                "negative_margin_lower"
            ],
            "weakest_box": weakest,
            "records": records,
        }


def rows() -> list[GateRow]:
    return [
        GateRow(
            "nrzb361_01_theta_enclosure",
            "interval_input",
            "ready_to_apply",
            "Four theta summands plus an analytic geometric envelope enclose the exact Xi kernel.",
            "Phi=sum_(n=1)^4 Phi_n+E_ge5, E_ge5>=sum_(n>=5)Phi_n",
            "The envelope is used only inside the directed boundary calculation.",
        ),
        GateRow(
            "nrzb361_02_taylor_enclosure",
            "interval_method",
            "ready_to_apply",
            "A bivariate Taylor model encloses the real part on each rational parameter box.",
            "orders (t,y)=(2,4), steps (1/50,1/20)",
            "Absolute positive moments enclose both Taylor remainders.",
        ),
        GateRow(
            "nrzb361_03_vertical_boundary",
            "new_theorem",
            "ready_to_apply",
            "The vertical sides of the central zero rectangle are zero-free.",
            "Re H_t(38+iy)<0 for 0<=t<=1/5 and |y|<=1",
            "Evenness and conjugation give the side Re z=-38.",
        ),
        GateRow(
            "nrzb361_04_horizontal_boundary",
            "published_composition",
            "ready_to_apply",
            "The horizontal sides are zero-free.",
            "|Im z|<=sqrt(1-2t)<1 for t>0; H_0(z)=xi((1+iz)/2)/8 at t=0",
            "The endpoint uses the classical zero-free lines Re s=0,1.",
        ),
        GateRow(
            "nrzb361_05_top_time_reality",
            "published_composition",
            "ready_to_apply",
            "Every H_(1/5) zero is real.",
            "Lambda<=1/5",
            "This imports Platt-Trudgian Corollary 2.",
        ),
        GateRow(
            "nrzb361_06_no_collision",
            "imported_certificate",
            "ready_to_apply",
            "No real multiple zero occurs in the central band during the homotopy.",
            "(H_t(x),H_t'(x))!=(0,0), 0<=t<=1/5, |x|<=38",
            "This imports the 1,900-box compact transversality certificate.",
        ),
        GateRow(
            "nrzb361_07_real_zero_band",
            "new_theorem",
            "ready_to_apply",
            "Every zero in the central real-width band is real.",
            "H_t(z)=0 and |Re z|<=38 imply Im z=0",
            "Boundary exclusion, top-time reality, and no real collision give the homotopy.",
        ),
        GateRow(
            "nrzb361_08_sector_transfer",
            "exact_transfer",
            "ready_to_apply",
            "Every squared-variable zero lies in a narrow negative-axis sector.",
            "delta=2 atan(1/38), sin(delta)=76/1445",
            "Real central zeros map to the negative axis; nonreal zeros have |Re z|>38 and |Im z|<=1.",
        ),
        GateRow(
            "nrzb361_09_chasse_transfer",
            "published_composition",
            "ready_to_apply",
            "Chasse's sector theorem applies to F_t and all derivatives.",
            "|sin(delta)|^-2=2088025/5776",
            "The order-one-half and positive-coefficient inputs are inherited from the degree-71 certificate.",
        ),
        GateRow(
            "nrzb361_10_degree361",
            "new_theorem",
            "ready_to_apply",
            "Every shifted Xi Jensen polynomial through degree 361 is hyperbolic with negative zeros.",
            "P_(d,n,t) hyperbolic for 0<=d<=361, n>=0, 0<=t<=1/5",
            "This is a finite-degree theorem and does not include degree 362.",
        ),
        GateRow(
            "nrzb361_11_degree71_sharpening",
            "scope_audit",
            "ready_to_apply",
            "The real-zero-band theorem strictly sharpens the prior zero-slab cutoff.",
            "71<361",
            "The earlier theorem remains an independent direct zero-slab proof.",
        ),
        GateRow(
            "nrzb361_12_nonpromotion",
            "proof_boundary",
            "open",
            "No all-degree or Newman-bound conclusion follows from a finite sector cutoff.",
            "degree 362, cofinal degrees, PF-infinity, Lambda<=0, and RH remain open",
            "The next sector upgrade needs a wider certified real-zero band or an unbounded sequence.",
        ),
    ]


def build_payload() -> dict[str, object]:
    lowered = request_below_normal_priority()
    contact = json.loads(CONTACT_SOURCE.read_text(encoding="utf-8"))
    degree71 = json.loads(DEGREE71_SOURCE.read_text(encoding="utf-8"))
    if contact.get("certificate", {}).get("unresolved_boxes") != 0:
        raise RuntimeError("contact source has unresolved boxes")
    if degree71.get("summary", {}).get("maximum_proved_degree") != 71:
        raise RuntimeError("degree-71 source changed")
    certificate = BoundaryCertifier().build()
    radius = REAL_BOUNDARY
    sector_sine = 2 * radius / (radius * radius + 1)
    inverse_sine_square = Fraction(1, 1) / sector_sine**2
    cutoff = (
        inverse_sine_square.numerator
        // inverse_sine_square.denominator
    )
    if sector_sine != Fraction(76, 1445) or cutoff != 361:
        raise RuntimeError("sector arithmetic drifted")
    return {
        "kind": (
            "jensen_window_pf_newman_real_zero_band_"
            "degree361_sector_certificate"
        ),
        "date": "2026-07-25",
        "status": (
            "rigorous real-zero-band homotopy and "
            "published-sector composition theorem"
        ),
        "priority_lowered": lowered,
        "sources": {
            "compact_no_collision": {
                "path": str(CONTACT_SOURCE.relative_to(REPO_ROOT)).replace(
                    "\\",
                    "/",
                ),
                "sha256": sha256(CONTACT_SOURCE),
                "theorem": contact["exact"]["proved_rectangle"],
            },
            "degree71_sector": {
                "path": str(DEGREE71_SOURCE.relative_to(REPO_ROOT)).replace(
                    "\\",
                    "/",
                ),
                "sha256": sha256(DEGREE71_SOURCE),
                "maximum_proved_degree": 71,
            },
            "debruijn_strip": {
                "url": "https://arxiv.org/abs/1904.12438",
                "theorem": "Theorem 3.2",
            },
            "lambda_upper": {
                "url": "https://doi.org/10.1090/mcom/3595",
                "result": "Corollary 2",
                "statement": "Lambda<=1/5",
            },
            "chasse_sector": {
                "url": "https://doi.org/10.1080/17476933.2011.584250",
                "theorem": "Theorem 3.6 and derivative closure",
            },
        },
        "exact": {
            "kernel": {
                "summand": (
                    "Phi_n(u)=pi*n^2*exp(5u)"
                    "*(2*pi*n^2*exp(4u)-3)"
                    "*exp(-pi*n^2*exp(4u))"
                ),
                "omitted_envelope": (
                    "sum_(n>=5)Phi_n(u)<="
                    "2*pi^2*5^4*exp(9u-25*pi*exp(4u))/"
                    "(1-(6/5)^4*exp(-11*pi))"
                ),
                "far_tail": (
                    "all required moments on u>=2 are enclosed "
                    "by symmetric radius 1e-800"
                ),
                "far_tail_audit": {
                    "moment_range": "0<=m<=9",
                    "exp8_partial_sum": "47259/35",
                    "exponent": (
                        "u^2/5+(29/2)u-pi*exp(4u)"
                    ),
                    "endpoint_upper": "-140734/35",
                    "derivative_rule": (
                        "the exponent derivative is below -1 "
                        "for u>=2"
                    ),
                    "prefactor_upper": "40",
                    "integer_comparison": (
                        "40*10^800<2^4020"
                    ),
                    "conclusion": (
                        "40*exp(-140734/35)<10^-800"
                    ),
                },
            },
            "boundary_certificate": certificate,
            "vertical_boundary_theorem": (
                "Re H_t(38+iy)<0 for "
                "0<=t<=1/5 and |y|<=1"
            ),
            "horizontal_boundary": (
                "For t>0, Theorem 3.2 gives "
                "|Im z|<=sqrt(1-2t)<1. At t=0, "
                "H_0(z)=xi((1+iz)/2)/8 and xi has no zeros "
                "on Re s=0 or Re s=1."
            ),
            "homotopy": {
                "rectangle": "|Re z|<=38, |Im z|<=1",
                "top_time": (
                    "Lambda<=1/5 implies every H_(1/5) zero is real"
                ),
                "no_collision": (
                    "(H_t(x),H_t'(x))!=(0,0) for "
                    "0<=t<=1/5 and |x|<=38"
                ),
                "argument": (
                    "Boundary zero-freeness fixes the finite zero "
                    "multiset in the rectangle. Simple real roots "
                    "continue as real roots; leaving the real axis "
                    "would require a real multiple zero."
                ),
                "theorem": (
                    "H_t(z)=0 and |Re z|<=38 imply Im z=0 "
                    "for every 0<=t<=1/5"
                ),
            },
            "sector_transfer": {
                "zero_map": "H_t(z)=0 => F_t(-z^2)=0",
                "sector_angle": "delta=2*atan(1/38)",
                "sector_sine": str(sector_sine),
                "inverse_sine_square": str(inverse_sine_square),
                "degree_cutoff": cutoff,
                "theorem": (
                    "P_(d,n,t)(x)=sum_(j=0)^d "
                    "binom(d,j)A_(n+j)(t)x^j is hyperbolic "
                    "with real negative zeros for every "
                    "0<=d<=361, n>=0, and 0<=t<=1/5"
                ),
            },
        },
        "summary": {
            "rows": 12,
            "precision_bits": PRECISION_BITS,
            "boundary_boxes": certificate["boxes"],
            "unresolved_boxes": certificate["unresolved_boxes"],
            "vertical_zero_free_boundaries": 2,
            "real_zero_bands": 1,
            "all_shifts": True,
            "maximum_proved_degree": cutoff,
            "previous_maximum_degree": 71,
            "all_degree_theorems": 0,
        },
        "rows": [asdict(row) for row in rows()],
    }


def success_line(payload: dict[str, object]) -> str:
    summary = payload["summary"]
    certificate = payload["exact"]["boundary_certificate"]
    return (
        "validated Newman real-zero-band degree-361 sector certificate: "
        f"{summary['rows']} rows, 0 issues, "
        f"{summary['boundary_boxes']} certified boundary boxes, "
        f"{summary['unresolved_boxes']} unresolved, "
        "1 real-zero band, all shifts through degree 361, "
        "0 all-degree theorems, weakest negative margin "
        f"{certificate['weakest_negative_margin_lower']}"
    )


def render_note(payload: dict[str, object]) -> str:
    certificate = payload["exact"]["boundary_certificate"]
    weak = certificate["weakest_box"]
    return f"""# Newman Real-Zero-Band Degree-361 Sector Certificate

Date: 2026-07-25

Status: rigorous Arb/Taylor boundary certificate, zero-homotopy theorem,
and published-sector composition. This is not a proof of PF-infinity,
`Lambda<=0`, RH, or a Clay-prize conclusion.

```text
work/rh_compute/results/jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.py
```

Current result:

```text
{success_line(payload)}
```

## Vertical Boundary

For

```text
H_t(z)=integral_0^infinity exp(t*u^2)*Phi(u)*cos(z*u)du,
```

the real part on `z=38+iy` is

```text
Re H_t(38+iy)=integral_0^infinity
 exp(t*u^2)*Phi(u)*cos(38*u)*cosh(y*u)du.
```

A 192-bit outward-rounded Taylor certificate proves

```text
Re H_t(38+iy)<0
for every 0<=t<=1/5 and |y|<=1.
```

The partition has `{certificate['boxes']}` rational boxes, Taylor orders
`(t,y)=({TIME_ORDER},{IMAGINARY_ORDER})`, steps
`({TIME_STEP},{IMAGINARY_STEP})`, and zero unresolved boxes. The weakest
certified negative margin is

```text
{certificate['weakest_negative_margin_lower']}
```

on `t=[{weak['t_low']},{weak['t_high']}]`,
`y=[{weak['y_low']},{weak['y_high']}]`. Four theta summands are integrated
directly. The omitted `n>=5` terms use

```text
2*pi^2*5^4*exp(9u-25*pi*exp(4u))
 /(1-(6/5)^4*exp(-11*pi)),
```

For every required `0<=m<=9`, `u^m<=exp(9u/2)` and
`cosh(yu)<=exp(u)`. The partial sum
`sum_(k=0)^7 8^k/k!=47259/35`, together with `pi>3` and `pi^2<10`,
gives the directed audit

```text
u^2/5+(29/2)u-pi*exp(4u)<=-140734/35-(u-2),
40*10^800<2^4020,
```

so the analytic `u>=2` moment tail is strictly below `1e-800`.
Evenness and conjugation give the same zero-free conclusion on `Re z=-38`.

## Real-Zero Homotopy

Polymath Theorem 3.2 gives

```text
H_t(z)=0 => |Im z|<=sqrt(1-2t)<1, 0<t<=1/5.
```

At `t=0`, `H_0(z)=xi((1+iz)/2)/8`; the horizontal sides `|Im z|=1`
map to the classical zero-free lines `Re s=0,1`. Thus the complete
rectangle boundary is zero-free. Platt--Trudgian Corollary 2 gives
`Lambda<=1/5`, so every zero at the top time is real. The independent
1,900-box compact transversality certificate proves

```text
(H_t(x),H_t'(x))!=(0,0)
for 0<=t<=1/5 and |x|<=38.
```

The argument principle fixes the finite zero multiset in the rectangle.
Simple real roots continue as real roots, and a pair could leave the real
axis only through a real multiple zero. Therefore

```text
H_t(z)=0 and |Re z|<=38 imply Im z=0
for every 0<=t<=1/5.
```

## Sector Transfer

All remaining nonreal zeros have `|Re z|>38`, while the strip theorem
gives `|Im z|<=1`. Under `s=-z^2`, every zero of
`F_t(s)=2H_t(i*sqrt(s))` therefore lies in the negative-axis sector

```text
delta=2*atan(1/38),
sin(delta)=76/1445,
|sin(delta)|^(-2)=2088025/5776=361.500... .
```

Chasse's sector theorem, the order-one-half property, coefficient
positivity, and derivative closure give

```text
P_(d,n,t)(x)=sum_(j=0)^d binom(d,j)A_(n+j)(t)x^j
is hyperbolic with real negative zeros
for every 0<=d<=361, n>=0, and 0<=t<=1/5.
```

Primary theorem inputs:

- D. H. J. Polymath, Theorem 3.2:
  <https://arxiv.org/abs/1904.12438>
- Platt--Trudgian, Corollary 2:
  <https://doi.org/10.1090/mcom/3595>
- Matthew Chasse, Theorem 3.6:
  <https://doi.org/10.1080/17476933.2011.584250>

## Proof Boundary

This strictly sharpens the independent degree-71 zero-slab theorem.
The cutoff `361` is finite. Degree `362`, an unbounded cofinal degree
sequence, all-degree Jensen hyperbolicity, PF-infinity, `Lambda<=0`, RH,
and the Clay prize remain open. The next sector upgrade needs a wider
real-zero band or a genuinely unbounded sector mechanism.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(success_line(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
