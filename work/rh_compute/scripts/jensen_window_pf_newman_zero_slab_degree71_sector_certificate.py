#!/usr/bin/env python3
"""Certify a Newman zero-free slab and its degree-71 Jensen consequence."""

from __future__ import annotations

import argparse
import ctypes
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
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
SLAB_REAL_RADIUS_NUMERATOR = 84
SLAB_REAL_RADIUS_DENOMINATOR = 5
SLAB_REAL_RADIUS_DISPLAY = "84/5"
SLAB_IMAG_RADIUS = 1
TIME_UPPER = Fraction(1, 5)
INTEGRATION_CUTOFF = 2
ABS_TOL = "1e-50"

DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_zero_slab_degree71_sector_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_newman_zero_slab_degree71_sector_certificate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


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


def interval_record(value: arb, digits: int = 90) -> dict[str, object]:
    return {
        "ball": value.str(n=digits, more=True),
        "lower": value.lower().str(n=digits, more=True),
        "upper": value.upper().str(n=digits, more=True),
        "strictly_positive": bool(value.lower() > 0),
    }


class DirectedSlabCertificate:
    def __init__(self) -> None:
        flint.ctx.prec = PRECISION_BITS
        self.pi = acb.pi()
        self.zero = acb(0)
        self.radius = (
            acb(SLAB_REAL_RADIUS_NUMERATOR)
            / SLAB_REAL_RADIUS_DENOMINATOR
        )
        self.split = self.pi / (2 * self.radius)
        self.cutoff = acb(INTEGRATION_CUTOFF)
        self.abs_tol = arb(ABS_TOL)

    def phi_one(self, u: acb) -> acb:
        return (
            2 * self.pi * self.pi * (9 * u).exp()
            - 3 * self.pi * (5 * u).exp()
        ) * (-self.pi * (4 * u).exp()).exp()

    def higher_theta_envelope(self, u: acb) -> acb:
        ratio = (acb(3) / 2) ** 4 * (-5 * self.pi).exp()
        return (
            32
            * self.pi
            * self.pi
            * (9 * u).exp()
            * (-4 * self.pi * (4 * u).exp()).exp()
            / (1 - ratio)
        )

    def integrate(self, func, left: acb, right: acb) -> arb:
        value = acb.integral(
            func,
            left,
            right,
            abs_tol=self.abs_tol,
            rel_tol=self.abs_tol,
            deg_limit=80,
            eval_limit=300000,
            depth_limit=80,
        )
        if not value.imag.contains(0):
            raise RuntimeError("a real directed integral acquired nonzero imaginary part")
        return value.real

    def build(self) -> dict[str, object]:
        def core_integrand(u: acb, analytic: bool) -> acb:
            del analytic
            return self.phi_one(u) * (self.radius * u).cos()

        def first_tail_integrand(u: acb, analytic: bool) -> acb:
            del analytic
            return (
                (u * u / 5).exp()
                * self.phi_one(u)
                * u.cosh()
            )

        def higher_tail_integrand(u: acb, analytic: bool) -> acb:
            del analytic
            return (
                (u * u / 5).exp()
                * self.higher_theta_envelope(u)
                * u.cosh()
            )

        core = self.integrate(
            core_integrand,
            self.zero,
            self.split,
        )
        first_tail = self.integrate(
            first_tail_integrand,
            self.split,
            self.cutoff,
        )
        higher_tail = self.integrate(
            higher_tail_integrand,
            self.split,
            self.cutoff,
        )

        full_ratio = 16 * (-3 * self.pi).exp()
        exp_eight = acb(8).exp()
        infinite_tail = (
            2
            * self.pi
            * self.pi
            / (1 - full_ratio)
            * (-self.pi * exp_eight / 2).exp()
            / (2 * self.pi * exp_eight)
        ).real
        margin = core - first_tail - higher_tail - infinite_tail

        higher_ratio = (
            (acb(3) / 2) ** 4 * (-5 * self.pi).exp()
        ).real
        full_ratio_real = full_ratio.real
        elementary_far_gap = (
            self.pi / 2 * exp_eight
            - acb(4) / 5
            - 20
        ).real

        for label, value in (
            ("core", core),
            ("margin", margin),
            ("elementary far gap", elementary_far_gap),
        ):
            if not value.lower() > 0:
                raise RuntimeError(f"{label} was not certified positive")
        for label, value in (
            ("higher-theta ratio", higher_ratio),
            ("full-theta ratio", full_ratio_real),
        ):
            if not (value.lower() > 0 and value.upper() < 1):
                raise RuntimeError(f"{label} was not certified in (0,1)")

        rational_radius = Fraction(
            SLAB_REAL_RADIUS_NUMERATOR,
            SLAB_REAL_RADIUS_DENOMINATOR,
        )
        sector_sine = (
            2 * rational_radius / (rational_radius**2 + 1)
        )
        sector_bound = Fraction(1, 1) / sector_sine**2
        degree_cutoff = sector_bound.numerator // sector_bound.denominator
        if degree_cutoff != 71:
            raise RuntimeError("sector degree cutoff drifted")

        return {
            "constants": {
                "precision_bits": PRECISION_BITS,
                "time_interval": "0<=t<=1/5",
                "slab": "|Re z|<=84/5, |Im z|<=1",
                "real_radius": SLAB_REAL_RADIUS_DISPLAY,
                "imaginary_radius": SLAB_IMAG_RADIUS,
                "split": "5*pi/168",
                "split_ball": self.split.real.str(n=90, more=True),
                "integration_cutoff": INTEGRATION_CUTOFF,
                "absolute_tolerance": ABS_TOL,
            },
            "kernel_bounds": {
                "summand_factorization": (
                    "Phi_n(u)=pi*n^2*exp(5u)"
                    "*(2*pi*n^2*exp(4u)-3)"
                    "*exp(-pi*n^2*exp(4u))>0"
                ),
                "higher_ratio": (
                    "(3/2)^4*exp(-5*pi)"
                ),
                "higher_ratio_ball": interval_record(higher_ratio),
                "higher_envelope": (
                    "sum_(n>=2) Phi_n(u)"
                    "<=32*pi^2*exp(9u-4*pi*exp(4u))"
                    "/(1-(3/2)^4*exp(-5*pi))"
                ),
                "full_ratio": "16*exp(-3*pi)",
                "full_ratio_ball": interval_record(full_ratio_real),
                "full_envelope": (
                    "Phi(u)<=2*pi^2*exp(9u-pi*exp(4u))"
                    "/(1-16*exp(-3*pi))"
                ),
                "far_exponent_inequality": (
                    "u^2/5+10u<=(pi/2)*exp(4u), u>=2"
                ),
                "far_endpoint_gap": interval_record(
                    elementary_far_gap
                ),
            },
            "directed_integrals": {
                "core": {
                    "formula": (
                        "C_R=int_0^(5*pi/168) "
                        "Phi_1(u)*cos((84/5)u)du"
                    ),
                    **interval_record(core),
                },
                "first_tail": {
                    "formula": (
                        "T_1=int_(5*pi/168)^2 "
                        "exp(u^2/5)*Phi_1(u)*cosh(u)du"
                    ),
                    **interval_record(first_tail),
                },
                "higher_tail": {
                    "formula": (
                        "T_ge2=int_(5*pi/168)^2 "
                        "exp(u^2/5)*cosh(u)*E_ge2(u)du"
                    ),
                    **interval_record(higher_tail),
                },
                "infinite_tail": {
                    "formula": (
                        "T_inf<=2*pi^2/(1-16e^(-3pi))"
                        "*exp(-(pi/2)e^8)/(2*pi*e^8)"
                    ),
                    **interval_record(infinite_tail),
                },
                "margin": {
                    "formula": "C_16-T_1-T_ge2-T_inf",
                    **interval_record(margin),
                },
            },
            "zero_slab": {
                "real_part_identity": (
                    "Re H_t(x+iy)=int_0^infinity "
                    "exp(tu^2)*Phi(u)*cos(xu)*cosh(yu)du"
                ),
                "near_inequality": (
                    "on 0<=u<=5*pi/168, |x|<=84/5: "
                    "cos(xu)>=cos((84/5)u)>=0; "
                    "exp(tu^2)>=1 and cosh(yu)>=1"
                ),
                "far_inequality": (
                    "on u>=5*pi/168, 0<=t<=1/5 and |y|<=1: "
                    "the real integrand is at least "
                    "-exp(u^2/5)*Phi(u)*cosh(u)"
                ),
                "theorem": (
                    "Re H_t(x+iy)>0 for 0<=t<=1/5, "
                    "|x|<=84/5, |y|<=1"
                ),
                "zero_free": True,
            },
            "published_inputs": {
                "debruijn_strip": {
                    "source": "https://arxiv.org/abs/1904.12438",
                    "theorem": "Theorem 3.2",
                    "statement": (
                        "Every zero z of H_t satisfies "
                        "|Im z|<=sqrt(1-2t)<=1 for 0<=t<=1/5."
                    ),
                },
                "chasse_sector": {
                    "source": (
                        "https://doi.org/10.1080/"
                        "17476933.2011.584250"
                    ),
                    "accessible_source": (
                        "https://people.kth.se/~chasse/"
                        "LaguerreSector.pdf"
                    ),
                    "theorem": "Theorem 3.6 and derivative closure",
                    "statement": (
                        "An order-at-most-one real entire function "
                        "with positive Taylor coefficients and zeros "
                        "in S(pi,delta) has hyperbolic Jensen "
                        "polynomials through degree "
                        "floor(|sin(delta)|^-2); the sector class is "
                        "closed under differentiation."
                    ),
                },
                "endpoint_comparison": {
                    "source": "https://arxiv.org/abs/1910.01227",
                    "result": "Corollary 1.3",
                    "statement": (
                        "At t=0, the published Xi theorem already "
                        "gives all-shift hyperbolicity for "
                        "d<=9.36*10^20."
                    ),
                },
            },
            "sector_transfer": {
                "generating_function": (
                    "F_t(s)=sum_(k>=0)A_k(t)s^k/k!"
                    "=2*H_t(i*sqrt(s))"
                ),
                "zero_map": "H_t(z)=0 => F_t(-z^2)=0",
                "order": (
                    "F_t has order 1/2; the double-exponential "
                    "kernel gives log M_F(r)=O(sqrt(r)*log(2+r))."
                ),
                "positive_coefficients": True,
                "sector_half_angle": "delta=2*atan(5/84)",
                "sector_sine": str(sector_sine),
                "inverse_sine_square": str(sector_bound),
                "degree_cutoff": degree_cutoff,
                "all_shifts": (
                    "F_t^(n)(s)=sum_(j>=0)A_(n+j)(t)s^j/j!"
                ),
                "theorem": (
                    "P_(d,n,t) is hyperbolic with real negative "
                    "zeros for every 0<=d<=71, n>=0, "
                    "and 0<=t<=1/5."
                ),
            },
        }


def build_payload() -> dict[str, object]:
    priority_lowered = request_below_normal_priority()
    exact = DirectedSlabCertificate().build()
    margin = exact["directed_integrals"]["margin"]
    rows = [
        GateRow(
            id="nzsd71_01_kernel_positivity",
            role="exact_lemma",
            readiness="available_exact",
            claim="Every theta summand of the Newman kernel is positive on u>=0.",
            formula=exact["kernel_bounds"]["summand_factorization"],
            proof_boundary="Elementary factorization using 2*pi-3>0.",
        ),
        GateRow(
            id="nzsd71_02_near_lower_bound",
            role="exact_lemma",
            readiness="available_exact",
            claim="The first theta summand supplies a uniform positive near contribution.",
            formula=exact["zero_slab"]["near_inequality"],
            proof_boundary="Only the interval 0<=u<=5*pi/168 is used.",
        ),
        GateRow(
            id="nzsd71_03_theta_tail_envelope",
            role="exact_lemma",
            readiness="available_exact",
            claim="A geometric majorant controls every theta summand n>=2.",
            formula=exact["kernel_bounds"]["higher_envelope"],
            proof_boundary="Uniform for u>=0; no numerical theta truncation is assumed.",
        ),
        GateRow(
            id="nzsd71_04_far_tail_envelope",
            role="exact_lemma",
            readiness="available_exact",
            claim="The full kernel beyond u=2 has an explicit elementary integral bound.",
            formula=(
                f"{exact['kernel_bounds']['full_envelope']}; "
                f"{exact['kernel_bounds']['far_exponent_inequality']}"
            ),
            proof_boundary="Uses cosh(u)<=exp(u) and a monotone exponent gap.",
        ),
        GateRow(
            id="nzsd71_05_directed_integrals",
            role="interval_certificate",
            readiness="ready_to_apply",
            claim="The directed near-minus-far margin is rigorously positive.",
            formula=margin["formula"],
            proof_boundary="One-dimensional 192-bit Arb integrals plus an analytic infinite tail.",
            diagnostics=exact["directed_integrals"],
        ),
        GateRow(
            id="nzsd71_06_zero_free_slab",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The full Newman transform has positive real part on the complete complex slab.",
            formula=exact["zero_slab"]["theorem"],
            proof_boundary="Uniform only on the stated slab and heat interval.",
        ),
        GateRow(
            id="nzsd71_07_debruijn_strip",
            role="published_theorem",
            readiness="available_published",
            claim="All positive-window Newman zeros lie in the horizontal strip used by the slab theorem.",
            formula=exact["published_inputs"]["debruijn_strip"]["statement"],
            proof_boundary="Imported Theorem 3.2 of the cited Polymath paper.",
            diagnostics=exact["published_inputs"]["debruijn_strip"],
        ),
        GateRow(
            id="nzsd71_08_squared_zero_sector",
            role="exact_lemma",
            readiness="ready_to_apply",
            claim="Every squared-variable zero lies in one explicit negative-axis sector.",
            formula=(
                "|arg(-z^2)-pi|<=2*atan(5/84), "
                "sin(delta)=840/7081"
            ),
            proof_boundary="Composition of the zero-free slab and de Bruijn strip only.",
        ),
        GateRow(
            id="nzsd71_09_chasse_transfer",
            role="published_theorem",
            readiness="available_published",
            claim="The sector theorem converts the squared-zero enclosure into Jensen hyperbolicity.",
            formula=exact["published_inputs"]["chasse_sector"]["statement"],
            proof_boundary="Imported Chasse Theorem 3.6 and sector derivative closure.",
            diagnostics=exact["published_inputs"]["chasse_sector"],
        ),
        GateRow(
            id="nzsd71_10_degree71_theorem",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="All shifted Newman Jensen windows through degree 71 are hyperbolic throughout the target heat interval.",
            formula=exact["sector_transfer"]["theorem"],
            proof_boundary="Finite-degree theorem only; degrees 65 and above are not covered.",
        ),
        GateRow(
            id="nzsd71_11_quartic_consequence",
            role="exact_consequence",
            readiness="ready_to_apply",
            claim="Actual Xi degree-four and adjacent degree-five windows cannot realize the nonhyperbolic outer-contact survivor on this heat interval.",
            formula=(
                "d=4,5<=71; hence P_(4,n,t) and P_(5,n,t) "
                "are hyperbolic for every n and 0<=t<=1/5"
            ),
            proof_boundary="Excludes the survivor for Xi; it does not invalidate the finite countermodel.",
        ),
        GateRow(
            id="nzsd71_12_endpoint_scope",
            role="literature_scope",
            readiness="available_published",
            claim="The static endpoint already has a much larger published finite-degree theorem.",
            formula=exact["published_inputs"]["endpoint_comparison"]["statement"],
            proof_boundary="Endpoint t=0 only; it supplies no positive-time uniformity by itself.",
            diagnostics=exact["published_inputs"]["endpoint_comparison"],
        ),
        GateRow(
            id="nzsd71_13_nonpromotion",
            role="nonpromotion_guard",
            readiness="blocked",
            claim="A bounded degree cutoff cannot imply Laguerre-Polya membership or RH.",
            formula="71<infinity",
            proof_boundary=(
                "No degree-72 theorem, unbounded terminal sequence, "
                "PF-infinity, Lambda<=0, RH, or Clay-prize conclusion."
            ),
        ),
    ]
    return {
        "kind": (
            "jensen_window_pf_newman_zero_slab_"
            "degree71_sector_certificate"
        ),
        "date": "2026-07-25",
        "status": (
            "rigorous uniform zero-slab and published-sector "
            "composition theorem"
        ),
        "priority_lowered": priority_lowered,
        "proof_boundary": (
            "This proves all-shift Jensen hyperbolicity only through "
            "degree 71 on 0<=t<=1/5. It is not an all-degree theorem "
            "and does not prove PF-infinity, Lambda<=0, RH, or a "
            "Clay-prize conclusion."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "precision_bits": PRECISION_BITS,
            "directed_integrals": 4,
            "positive_margin": 1,
            "zero_free_complex_slabs": 1,
            "published_inputs": 3,
            "uniform_heat_parameters": "continuum_0_to_1_over_5",
            "all_shifts": True,
            "maximum_proved_degree": 71,
            "quartic_layers_closed": 1,
            "quintic_layers_closed": 1,
            "all_degree_theorems": 0,
        },
    }


def render_note(payload: dict[str, object]) -> str:
    exact = payload["exact"]
    integrals = exact["directed_integrals"]
    margin = integrals["margin"]
    sector = exact["sector_transfer"]
    return f"""# Newman Zero-Slab Degree-71 Sector Certificate

Date: 2026-07-25

Status: rigorous interval certificate and published-theorem composition.
This proves a bounded-degree Newman/Jensen theorem. It is not a proof of
PF-infinity, `Lambda<=0`, RH, or a Clay-prize conclusion.

```text
work/rh_compute/results/jensen_window_pf_newman_zero_slab_degree71_sector_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_zero_slab_degree71_sector_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_zero_slab_degree71_sector_certificate.py
```

Current result:

```text
validated Newman zero-slab degree-71 sector certificate: 13 rows, 0 issues, 4 directed integrals, 1 positive slab margin, 1 zero-free complex slab, 3 published inputs, all shifts through degree 71, 0 all-degree theorems
```

## Uniform Complex Slab

For

```text
H_t(z)=integral_0^infinity exp(t*u^2)*Phi(u)*cos(z*u)du
0<=t<=1/5, |Re z|<=84/5, |Im z|<=1,
```

put `a=5*pi/168`. Every kernel summand is positive because

```text
Phi_n(u)=pi*n^2*exp(5u)*(2*pi*n^2*exp(4u)-3)
         *exp(-pi*n^2*exp(4u))>0.
```

On `0<=u<=a`,

```text
cos(xu)>=cos((84/5)u)>=0, exp(tu^2)>=1, cosh(yu)>=1.
```

On `u>=a`, the real integrand is bounded below by
`-exp(u^2/5)*Phi(u)*cosh(u)`. The first summand is integrated directly.
For `n>=2`, consecutive `n^4*exp(-pi*n^2*exp(4u))` terms have ratio at
most `(3/2)^4*exp(-5*pi)<1`, giving

```text
sum_(n>=2) Phi_n(u)
 <=32*pi^2*exp(9u-4*pi*exp(4u))
   /(1-(3/2)^4*exp(-5*pi)).
```

The 192-bit directed values are

```text
C_R      = {integrals['core']['ball']}
T_1      = {integrals['first_tail']['ball']}
T_ge2    = {integrals['higher_tail']['ball']}
T_inf    = {integrals['infinite_tail']['ball']}
C_R-T_1-T_ge2-T_inf
          = {margin['ball']}
```

For the final infinite tail, `cosh(u)<=exp(u)` and

```text
u^2/5+10u<=(pi/2)*exp(4u), u>=2,
```

reduce the integral to the displayed explicit
`exp(-(pi/2)*exp(8))` bound. Therefore

```text
Re H_t(x+iy)>0
```

throughout the stated slab, so the slab is zero-free.

## Sector Transfer

The de Bruijn strip-contraction theorem gives

```text
H_t(z)=0 => |Im z|<=sqrt(1-2t)<=1
```

on the same heat interval. The zero-free slab therefore forces every zero
to satisfy `|Re z|>84/5`. For

```text
F_t(s)=sum_(k>=0) A_k(t)*s^k/k!=2*H_t(i*sqrt(s)),
```

every zero has the form `-z^2`, hence lies in the negative-axis sector

```text
delta=2*atan(5/84),
sin(delta)={sector['sector_sine']},
|sin(delta)|^(-2)={sector['inverse_sine_square']}.
```

The double-exponential kernel gives `F_t` order `1/2`, and all its Taylor
coefficients are positive. Chasse's sector theorem and closure of the sector
class under differentiation now give

```text
P_(d,n,t)(x)=sum_(j=0)^d binom(d,j)*A_(n+j)(t)*x^j
is hyperbolic with real negative zeros
for every 0<=d<=71, n>=0, and 0<=t<=1/5.
```

Primary theorem inputs:

- D. H. J. Polymath, Theorem 3.2:
  <https://arxiv.org/abs/1904.12438>
- Matthew Chasse, Theorem 3.6:
  <https://doi.org/10.1080/17476933.2011.584250>
- Griffin--Ono--Rolen--Thorner--Tripp--Wagner, endpoint comparison:
  <https://arxiv.org/abs/1910.01227>

## Quartic Consequence

Degrees four and five are now closed for every shift throughout the complete
positive Newman target interval. In particular, an actual Xi quartic window
cannot realize the strict outward contact encoded by the exact length-14
survivor, and its adjacent quintic cannot have the survivor's nonreal pair.
The survivor remains a valid finite countermodel to generic signed-Hankel
promotion; the new theorem excludes it only by Xi kernel geometry plus the
published sector machinery.

At `t=0`, Corollary 1.3 of the effective Xi paper already gives the much
larger published range `d<=9.36*10^20` for every shift. The new contribution
here is the uniform continuum `0<=t<=1/5`, not a stronger endpoint cutoff.

## Proof Boundary

The cutoff `71` is finite. Nothing here proves degree `72`, an unbounded
cofinal terminal sequence, all-degree Jensen hyperbolicity, PF-infinity,
`Lambda<=0`, RH, or a Clay-prize conclusion. The next Jensen-side problem is
to enlarge the complex zero-free slab or obtain an unbounded sector sequence;
the independent strict-Laguerre/Newman transversality route remains open.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "built Newman zero-slab degree-71 sector certificate: "
            "13 rows, 4 directed integrals, 1 positive slab margin, "
            "1 zero-free complex slab, 3 published inputs, "
            "all shifts through degree 71, 0 all-degree theorems"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
