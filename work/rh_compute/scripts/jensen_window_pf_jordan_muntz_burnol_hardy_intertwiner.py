#!/usr/bin/env python3
"""Build the Jordan-Muntz/Burnol Hardy-intertwiner audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.md"
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
            "jmbh_01_fractional_part_kernel",
            "exact_identity",
            "available_exact",
            "R_0(y)=floor(y)-y=-{y}.",
            "Endpoint notation only; the infinite omega=0 Mobius tail is not "
            "asserted to converge.",
        ),
        row(
            "jmbh_02_fractional_power_abel_identity",
            "exact_identity",
            "available_exact",
            "R_omega(y)=-y^omega*{y}+omega*integral_0^y "
            "{u}u^(omega-1)du.",
            "Exact Stieltjes/Abel summation for 0<omega<1.",
        ),
        row(
            "jmbh_03_hardy_operator_definition",
            "exact_definition",
            "available_exact",
            "(T_omega f)(x)=x^omega*f(x)-omega*integral_0^x "
            "f(u)u^(omega-1)du.",
            "Definition only.",
        ),
        row(
            "jmbh_04_remainder_intertwining",
            "exact_identity",
            "available_exact",
            "R_omega=T_omega R_0.",
            "Exact kernel identity; it supplies no Mobius cancellation.",
        ),
        row(
            "jmbh_05_dilation_covariance",
            "exact_identity",
            "available_exact",
            "T_omega D_d=d^omega D_d T_omega for "
            "(D_d f)(x)=f(x/d).",
            "Exact change of variables.",
        ),
        row(
            "jmbh_06_finite_jordan_factorization",
            "exact_identity",
            "available_exact",
            "E_(omega,N)=T_omega sum_(d<=N)mu(d)d^(-2omega)"
            "R_0(x/d).",
            "Finite identity; no infinite interchange is used.",
        ),
        row(
            "jmbh_07_burnol_partial_sum",
            "exact_definition",
            "available_exact",
            "f_(epsilon,N)(t)=sum_(d<=N)mu(d)d^(-epsilon)"
            "{1/(d*t)}.",
            "This is the finite natural fractional-part approximant.",
        ),
        row(
            "jmbh_08_reciprocal_unitary",
            "exact_identity",
            "available_exact",
            "f(x)->f(1/t) is unitary from L2(dx/x^2) to L2(dt).",
            "Exact substitution x=1/t.",
        ),
        row(
            "jmbh_09_tail_hardy_operator",
            "exact_definition",
            "available_exact",
            "(H_star g)(t)=integral_t^infinity g(v)dv/v.",
            "Definition on its natural L2 domain.",
        ),
        row(
            "jmbh_10_exact_burnol_intertwiner",
            "exact_identity",
            "available_exact",
            "E_(omega,N)(1/t)=-(I-omega H_star)"
            "[t^(-omega)f_(2omega,N)(t)].",
            "Exact finite identity with epsilon=2omega.",
        ),
        row(
            "jmbh_11_hardy_norm",
            "classical_theorem",
            "source_backed",
            "||H_star||_(L2(dt)->L2(dt))=2.",
            "Classical Hardy inequality; not a zeta estimate.",
        ),
        row(
            "jmbh_12_two_sided_norm_equivalence",
            "exact_corollary",
            "available_exact",
            "(1-2omega)||t^(-omega)f_(2omega,N)||_2 "
            "<=||E_(omega,N)||_H<="
            "(1+2omega)||t^(-omega)f_(2omega,N)||_2.",
            "Uses 0<omega<1/2; constants degenerate only at the excluded "
            "endpoint.",
        ),
        row(
            "jmbh_13_invertibility",
            "exact_corollary",
            "available_exact",
            "I-omega H_star is boundedly invertible by its Neumann series "
            "for 0<omega<1/2.",
            "Follows from omega*||H_star||<1.",
        ),
        row(
            "jmbh_14_mellin_multiplier",
            "exact_identity",
            "available_exact",
            "M[H_star g](s)=M[g](s)/s and the intertwiner multiplier is "
            "(s-omega)/s.",
            "Initially in a common honest convergence strip.",
        ),
        row(
            "jmbh_15_transform_consistency",
            "exact_identity",
            "available_exact",
            "M[t^(-omega)f_(2omega,N)](s)="
            "-zeta(s-omega)M_N(s+omega)/(s-omega), so the exact "
            "intertwiner recovers zeta(s-omega)M_N(s+omega)/s.",
            "Finite Mobius sum only; analytic continuation is not used to "
            "move an infinite Mobius series.",
        ),
        row(
            "jmbh_16_cofinal_burnol_corollary",
            "theorem_candidate",
            "source_backed_candidate",
            "RH iff there are omega_j->0 such that "
            "sup_N||E_(omega_j,N)||_H<infinity for every j.",
            "The reverse direction uses the two-sided norm, pointwise "
            "convergence, Fatou, and Burnol's published sequence theorem. "
            "The RH direction uses Balazard-Saias/Burnol L2 convergence. "
            "Normalization and source transfer still require expert review.",
        ),
        row(
            "jmbh_17_nonpromotion_guards",
            "proof_guard",
            "guard_validated",
            "Hardy boundedness transfers the norm problem but does not bound "
            "the Burnol partial sums; compact-support Nyman sufficiency and "
            "scalar l2 arguments remain unavailable.",
            "Rejects treating an invertible operator identity as the missing "
            "arithmetic estimate.",
        ),
        row(
            "jmbh_18_open_natural_mollifier_gate",
            "open_target",
            "open_target",
            "Prove sup_N||t^(-omega)f_(2omega,N)||_2<infinity along one "
            "explicit cofinal omega sequence without assuming RH.",
            "Equivalent by exact constants to the open Jordan all-height "
            "natural-mollifier bound; RH and Lambda <= 0 remain unproved.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner",
        "date": "2026-07-23",
        "status": (
            "exact finite Hardy intertwiner with one source-backed cofinal "
            "theorem candidate and one open natural-mollifier gate"
        ),
        "rows": rows,
        "primary_sources": [
            {
                "title": (
                    "On an analytic estimate in the theory of the Riemann "
                    "Zeta function and a Theorem of Baez-Duarte"
                ),
                "url": "https://arxiv.org/abs/math/0202166",
                "use": (
                    "weighted fractional-part functions, Balazard-Saias "
                    "convergence, square-integrability-to-zero-cancellation, "
                    "and cofinal RH theorem"
                ),
            },
            {
                "title": (
                    "A general strong Nyman-Beurling Criterion for the "
                    "Riemann Hypothesis"
                ),
                "url": "https://arxiv.org/abs/math/0505453",
                "use": "general Muntz framework and compact-support guard",
            },
            {
                "title": (
                    "On monotonicity of certain weighted summatory "
                    "functions associated with L-functions"
                ),
                "url": "https://arxiv.org/abs/1204.1823",
                "use": "Jordan coefficients and shifted quotient",
            },
        ],
        "audit": {
            "row_count": 18,
            "exact_identity_or_corollary_count": 13,
            "classical_source_backed_count": 1,
            "source_backed_cofinal_candidate_count": 1,
            "nonpromotion_guard_count": 2,
            "open_natural_mollifier_gate_count": 1,
            "uniform_partial_norm_proved": False,
            "rh_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Jordan-Muntz/Burnol Hardy Intertwiner

Date: 2026-07-23

Status: exact finite operator bridge with one source-backed cofinal theorem
candidate and one open natural-mollifier gate. This is not a proof of RH or
`Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.json
python work/rh_compute/scripts/jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.py
python work/rh_compute/scripts/check_jensen_window_pf_jordan_muntz_burnol_hardy_intertwiner.py
```

## Fractional-Part Intertwiner

Put

```text
R_0(y)=floor(y)-y=-{y}.
```

Stieltjes integration by parts gives, for `0<omega<1`,

```text
R_omega(y)
 =-y^omega*{y}
  +omega*integral_0^y {u}u^(omega-1)du.             (JMBH.1)
```

Define

```text
(T_omega f)(x)
 =x^omega*f(x)
  -omega*integral_0^x f(u)u^(omega-1)du.            (JMBH.2)
```

Then (JMBH.1) is exactly

```text
R_omega=T_omega R_0.                                (JMBH.3)
```

For `D_d f(x)=f(x/d)`, a change of variables gives

```text
T_omega D_d=d^omega D_d T_omega.                    (JMBH.4)
```

Consequently the finite Jordan-Muntz partial sum factors as

```text
E_(omega,N)(x)
 =T_omega[
   sum_(d<=N)mu(d)d^(-2omega)R_0(x/d)].             (JMBH.5)
```

No infinite sum or limiting argument occurs in (JMBH.5).

## Burnol Coordinates

Use Burnol's finite natural fractional-part approximant

```text
f_(epsilon,N)(t)
 =sum_(d<=N)mu(d)d^(-epsilon){1/(d*t)}.             (JMBH.6)
```

The reciprocal substitution `x=1/t` is unitary:

```text
integral_0^infinity |F(x)|^2 dx/x^2
 =integral_0^infinity |F(1/t)|^2 dt.                (JMBH.7)
```

Let

```text
(H_star g)(t)=integral_t^infinity g(v)dv/v.         (JMBH.8)
```

Applying (JMBH.5) at `x=1/t`, with `epsilon=2omega`, gives the exact finite
intertwiner

```text
E_(omega,N)(1/t)
 =-(I-omega*H_star)
   [t^(-omega)f_(2omega,N)(t)].                     (JMBH.9)
```

This identifies the Jordan partial sums with the weighted functions used in
Burnol's discussion of Baez-Duarte's natural approximants.

## Two-Sided Norm

The classical tail Hardy operator satisfies

```text
||H_star||_(L2(dt)->L2(dt))=2.
```

One proof sets `G=H_star g`, uses `G'=-g/t`, and integrates
`integral G^2` by parts before Cauchy-Schwarz. Therefore, for
`0<omega<1/2`,

```text
(1-2omega)||t^(-omega)f_(2omega,N)||_2
 <=||E_(omega,N)||_H
 <=(1+2omega)||t^(-omega)f_(2omega,N)||_2.          (JMBH.10)
```

Moreover `I-omega H_star` is boundedly invertible by its Neumann series.
Thus uniform boundedness, Cauchy convergence, and convergence of either
finite family are equivalent to the corresponding property of the other.
The constants improve, rather than deteriorate, along `omega->0`.

## Mellin Check

In a common convergence strip,

```text
M[H_star g](s)=M[g](s)/s,

M[t^(-omega)f_(2omega,N)](s)
 =-zeta(s-omega)M_N(s+omega)/(s-omega).             (JMBH.11)
```

The multiplier of `I-omega H_star` is `(s-omega)/s`. Combining this with
the minus sign in (JMBH.9) recovers

```text
M[E_(omega,N)(1/t)](s)
 =zeta(s-omega)M_N(s+omega)/s,                      (JMBH.12)
```

which is the independently obtained Jordan-Muntz transform.

## Cofinal Consequence

Suppose `omega_j->0` and

```text
sup_N ||E_(omega_j,N)||_H<infinity
```

for every `j`. Equation (JMBH.10) uniformly bounds
`t^(-omega_j)f_(2omega_j,N)`. The infinite fractional-part series converges
pointwise for each positive epsilon, so Fatou gives the weighted Burnol
function in `L2`. On `(0,1)`, this implies the unweighted function is in
`L2`; on `(1,infinity)` it is an explicit constant times `1/t`. Burnol's
published sequence theorem then implies RH.

Conversely, under RH the Balazard-Saias estimate used by Burnol gives `L2`
convergence of these weighted natural partial sums. The bounded operator in
(JMBH.9) transfers convergence to `E_(omega,N)`. This yields the
source-backed cofinal theorem candidate

```text
RH
iff there are omega_j->0 such that
    sup_N ||E_(omega_j,N)||_H<infinity
    for every j.                                    (JMBH.13)
```

The normalization transfer and published-theorem composition in (JMBH.13)
still warrant independent expert review, but the route no longer depends
only on the new archimedean factorization.

## What Remains Open

The operator identity does not provide the arithmetic estimate. By
(JMBH.10), the live target is equivalently

```text
sup_N ||t^(-omega)f_(2omega,N)||_2<infinity          (JMBH.14)
```

on one explicit cofinal sequence, without assuming RH. Burnol explicitly
locates these natural Mobius approximants inside the Nyman programme; the
Balazard-Saias convergence estimate used there is conditional on RH.

Baez-Duarte's general strong Nyman theorem does not bypass (JMBH.14): the
audited sufficiency direction requires compact support, and `R_omega` is not
compactly supported. Hardy boundedness also cannot be mistaken for a bound
on the arithmetic input.

Primary sources:

- Jean-Francois Burnol, `On an analytic estimate in the theory of the Riemann Zeta function and a Theorem of Baez-Duarte`: https://arxiv.org/abs/math/0202166
- Luis Baez-Duarte, `A general strong Nyman-Beurling Criterion for the Riemann Hypothesis`: https://arxiv.org/abs/math/0505453
- Masatoshi Suzuki, `On monotonicity of certain weighted summatory functions associated with L-functions`: https://arxiv.org/abs/1204.1823

## Proof Boundary

Equations (JMBH.1)-(JMBH.12) are finite identities, unitary substitutions,
or consequences of the classical Hardy inequality. Equation (JMBH.13) is a
source-backed theorem candidate whose normalization chain requires external
review. Equation (JMBH.14) is open. RH, PF-infinity, and `Lambda <= 0`
remain unproved.
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.result.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.result.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(), encoding="utf-8")
    print(
        "wrote Jordan-Muntz/Burnol Hardy intertwiner: "
        f"{len(payload['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
