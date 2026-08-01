#!/usr/bin/env python3
"""Build the fixed-shift inner/reciprocal-boundary separation gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_fixed_shift_inner_"
    "reciprocal_boundary_separation_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_fixed_shift_inner_"
    "reciprocal_boundary_separation_gate.md"
)


def row(
    row_id: str,
    role: str,
    status: str,
    statement: str,
    proof_boundary: str,
) -> dict[str, str]:
    return {
        "id": row_id,
        "role": role,
        "status": status,
        "statement": statement,
        "proof_boundary": proof_boundary,
    }


def build_payload() -> dict:
    rows = [
        row(
            "fsirbsg_01_parameters",
            "exact_definition",
            "available_exact",
            "Fix omega>0 and T>0 in the centered variable "
            "z=s-1/2.",
            "The construction is finite-dimensional and makes no claim "
            "about zeros of the actual xi function.",
        ),
        row(
            "fsirbsg_02_toy_completed_function",
            "exact_definition",
            "available_exact",
            "Define F(z)=((z-omega)^2+T^2)"
            "((z+omega)^2+T^2).",
            "A real even quartic with the functional-equation symmetries "
            "of one off-axis zero quartet.",
        ),
        row(
            "fsirbsg_03_symmetry",
            "exact_identity",
            "available_exact",
            "F(-z)=F(z) and F(conj(z))=conj(F(z)).",
            "Immediate from the two paired quadratic factors.",
        ),
        row(
            "fsirbsg_04_zero_quartet",
            "exact_identity",
            "available_exact",
            "The zeros of F are z=+/-omega+/-iT, all simple.",
            "omega,T>0 keep the four roots distinct.",
        ),
        row(
            "fsirbsg_05_fixed_shift_quotient",
            "exact_definition",
            "available_exact",
            "Set Q(z)=F(z-omega)/F(z+omega).",
            "This is the centered analogue of a completed fixed-shift "
            "quotient.",
        ),
        row(
            "fsirbsg_06_boundary_zero_cancellation",
            "exact_identity",
            "available_exact",
            "Both numerator and denominator contain the factor z^2+T^2.",
            "The boundary zeros z=+/-iT cancel exactly in Q.",
        ),
        row(
            "fsirbsg_07_quotient_simplification",
            "exact_identity",
            "available_exact",
            "Q(z)=(((z-2omega)^2+T^2)/"
            "((z+2omega)^2+T^2)).",
            "Cancel the common boundary factor from row 6.",
        ),
        row(
            "fsirbsg_08_pole_location",
            "exact_identity",
            "available_exact",
            "The poles of Q are -2omega+/-iT, strictly in Re(z)<0.",
            "Hence Q is analytic on the open right half-plane.",
        ),
        row(
            "fsirbsg_09_boundary_unimodularity",
            "exact_identity",
            "available_exact",
            "|Q(it)|=1 for every real t.",
            "The numerator and denominator in row 7 are conjugates on "
            "the imaginary axis.",
        ),
        row(
            "fsirbsg_10_rational_inner",
            "exact_consequence",
            "available_exact",
            "Q is a rational inner function on Re(z)>0.",
            "It is the product of the two half-plane Blaschke factors "
            "(z-a)/(z+conj(a)) with a=2omega+/-iT.",
        ),
        row(
            "fsirbsg_11_causal_difference",
            "exact_identity",
            "available_exact",
            "Q(z)-1=-8omega*z/((z+2omega)^2+T^2).",
            "Direct subtraction in row 7.",
        ),
        row(
            "fsirbsg_12_exact_hardy_energy",
            "exact_identity",
            "available_exact",
            "The boundary integral int_R |Q(it)-1|^2 dt=16*pi*omega, "
            "so ||Q-1||_(H2)^2=8omega under the 1/(2pi) convention.",
            "A two-Cauchy-kernel partial fraction evaluates the integral "
            "exactly; in particular the causal boundary energy is finite.",
        ),
        row(
            "fsirbsg_13_reciprocal_boundary_energy",
            "exact_definition",
            "available_exact",
            "Define I_(omega,T)=int_R "
            "dt/[(1+t^2)|F(omega+it)|^2].",
            "Any positive smooth weight at t=+/-T gives the same "
            "divergence conclusion.",
        ),
        row(
            "fsirbsg_14_shifted_line_zero",
            "exact_identity",
            "available_exact",
            "F(omega+iT)=F(omega-iT)=0, and both zeros are simple "
            "along the shifted vertical line.",
            "The first quadratic factor vanishes simply and the second "
            "is nonzero.",
        ),
        row(
            "fsirbsg_15_reciprocal_divergence",
            "exact_consequence",
            "available_exact",
            "I_(omega,T)=infinity.",
            "Near t=T the integrand is a positive constant times "
            "|t-T|^(-2).",
        ),
        row(
            "fsirbsg_16_fixed_shift_nonpromotion",
            "proof_guard",
            "guard_validated",
            "Finite causal energy or fixed-shift innerness does not "
            "generically bound reciprocal square energy on the shifted "
            "line.",
            "The same symmetric boundary zero is canceled in Q but is a "
            "pole of 1/F(omega+it). No universal fixed-shift norm "
            "inequality can promote the former to the latter.",
        ),
        row(
            "fsirbsg_17_cofinal_distinction",
            "proof_guard",
            "guard_validated",
            "The countermodel blocks a same-shift promotion only; it does "
            "not block the audited cofinal implication from full "
            "determinant/Burnol control to RH.",
            "One off-axis zero quartet can align with an exceptional "
            "cancellation shift, not with every member of a sequence "
            "omega_j->0.",
        ),
        row(
            "fsirbsg_18_open_arithmetic_target",
            "open_target",
            "open_target",
            "Obtain sum_N P_(alpha/2,N)/N<infinity from a direct "
            "Mobius/Mertens estimate, or prove a cofinal full-Burnol "
            "theorem; do not infer it from limiting fixed-shift innerness "
            "alone.",
            "The logarithmic tail-energy estimate, RH, PF-infinity, and "
            "Lambda<=0 all remain open.",
        ),
    ]
    return {
        "kind": (
            "jensen_window_pf_fixed_shift_inner_"
            "reciprocal_boundary_separation_gate"
        ),
        "date": "2026-07-23",
        "status": (
            "exact rational-inner/reciprocal-boundary separation "
            "countermodel with one open arithmetic target"
        ),
        "parameters": {
            "omega": "omega>0",
            "height": "T>0",
            "half_plane": "Re(z)>0",
        },
        "rows": rows,
        "audit": {
            "row_count": 18,
            "exact_reduction_count": 15,
            "proof_guard_count": 2,
            "open_arithmetic_target_count": 1,
            "rational_inner_proved": True,
            "causal_energy_finite": True,
            "reciprocal_boundary_energy_finite": False,
            "fixed_shift_promotion_valid": False,
            "cofinal_full_burnol_implication_rejected": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Fixed-Shift Inner/Reciprocal-Boundary Separation Gate

Date: 2026-07-23

Status: exact rational-inner/reciprocal-boundary separation countermodel
with one open arithmetic target. This is not a proof of RH, PF-infinity,
or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.json
python work/rh_compute/scripts/jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.py
python work/rh_compute/scripts/check_jensen_window_pf_fixed_shift_inner_reciprocal_boundary_separation_gate.py
```

## Symmetric Four-Zero Model

Fix `omega>0`, `T>0`, and use the centered variable `z=s-1/2`. Define

```text
F(z)
 :=((z-omega)^2+T^2)((z+omega)^2+T^2).              (FSIRBSG.1)
```

This real even quartic obeys

```text
F(-z)=F(z),
F(conj(z))=conj(F(z)),
```

and has the simple zero quartet

```text
z=+/-omega+/-iT.                                     (FSIRBSG.2)
```

It is a finite toy completed function with the functional-equation
symmetries of one off-axis quartet. It makes no assertion that the actual
xi function has such a zero.

## The Fixed-Shift Quotient Is Inner

Form

```text
Q(z):=F(z-omega)/F(z+omega).                         (FSIRBSG.3)
```

Before cancellation,

```text
F(z-omega)
 =((z-2omega)^2+T^2)(z^2+T^2),

F(z+omega)
 =(z^2+T^2)((z+2omega)^2+T^2).
```

The boundary-zero factor cancels, leaving

```text
Q(z)
 =((z-2omega)^2+T^2)/((z+2omega)^2+T^2).             (FSIRBSG.4)
```

Its only poles are `-2omega+/-iT`, strictly in the left half-plane. On the
imaginary axis the numerator and denominator are conjugates, so

```text
|Q(it)|=1.                                            (FSIRBSG.5)
```

More explicitly, (FSIRBSG.4) is the product of the two right-half-plane
Blaschke factors

```text
(z-a)/(z+conj(a)),
a=2omega+iT and 2omega-iT.
```

Thus `Q` is rational inner on `Re(z)>0`.

The causal difference is

```text
Q(z)-1
 =-8omega*z/((z+2omega)^2+T^2).                      (FSIRBSG.6)
```

On the boundary,

```text
|Q(it)-1|^2
 =64omega^2*t^2
  /[(((t-T)^2+4omega^2)((t+T)^2+4omega^2))].
```

Using

```text
t^2/[D_-(t)D_+(t)]
 =t/(4T)[1/D_-(t)-1/D_+(t)],

D_-(t)=(t-T)^2+4omega^2,
D_+(t)=(t+T)^2+4omega^2,
```

and the elementary Cauchy integrals gives

```text
integral_R |Q(it)-1|^2 dt=16*pi*omega,

||Q-1||_(H2)^2=8omega                                  (FSIRBSG.7)
```

under the `1/(2*pi)` Hardy-norm convention. The causal boundary energy is
finite.

## Reciprocal Energy On The Shifted Line Diverges

Now inspect the reciprocal line energy

```text
I_(omega,T)
 :=integral_R dt/[(1+t^2)|F(omega+it)|^2].            (FSIRBSG.8)
```

The shifted line passes through the two zeros

```text
F(omega+iT)=F(omega-iT)=0.                            (FSIRBSG.9)
```

They are simple. Near `t=T`,

```text
1/[(1+t^2)|F(omega+it)|^2]
 ~C_(omega,T)/|t-T|^2,

C_(omega,T)
 =1/[64omega^2*T^2*(omega^2+T^2)*(1+T^2)]>0.          (FSIRBSG.10)
```

Therefore

```text
I_(omega,T)=infinity.                                  (FSIRBSG.11)
```

The same boundary zero is removable in the fixed-shift quotient but is a
pole of the reciprocal shifted-line function.

## Proof-Safety Consequence

This model rejects the generic promotion

```text
fixed-shift innerness or finite limiting causal energy
  =>
finite reciprocal square energy on the same shifted line.              (false)
```

In the RH corpus, the right-hand energy is the spectral form of the
weighted-prefix/logarithmic stable-prefix target. Limiting Jordan/Burnol
innerness cannot supply it at an exceptional fixed shift through a
universal norm multiplier.

The scope is deliberately narrow. The countermodel does not reject the
audited cofinal implication from full determinant or full Burnol control
to RH. One zero quartet can align with one exceptional cancellation shift;
it cannot align with every member of a sequence `omega_j->0`.

The live routes remain:

```text
1. prove sum_N P_(alpha/2,N)/N<infinity by a direct
   Mobius/Mertens mean-square estimate; or

2. prove the full Burnol norm on one cofinal shift sequence,
   which already implies RH through the existing three-gate criterion.
```

Neither route is closed here. The logarithmic tail-energy estimate, RH,
PF-infinity, and `Lambda <= 0` all remain open.
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    args.result.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    args.result.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(), encoding="utf-8")
    print(
        "wrote fixed-shift inner/reciprocal-boundary separation gate: "
        f"{len(payload['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
