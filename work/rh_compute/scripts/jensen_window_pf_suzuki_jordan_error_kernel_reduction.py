#!/usr/bin/env python3
"""Build the Suzuki Jordan-error kernel reduction and absolute-bound gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_jordan_error_kernel_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_suzuki_jordan_error_kernel_reduction.md"
)


def entry(
    row_id: str,
    role: str,
    status: str,
    formula: str,
    boundary: str,
) -> dict[str, str]:
    return {
        "id": row_id,
        "role": role,
        "status": status,
        "formula": formula,
        "proof_boundary": boundary,
    }


def build_payload() -> dict:
    rows = [
        entry(
            "jek_01_dirichlet_convolution",
            "exact_identity",
            "available_exact",
            "c_omega=id^omega*(mu*id^(-omega)); "
            "c_omega(n)=sum_(dm=n)mu(d)d^(-omega)m^omega.",
            "Exact coefficient identity only.",
        ),
        entry(
            "jek_02_cumulative_mass",
            "definition",
            "available_exact",
            "C_omega(x)=sum_(n<=x)c_omega(n), with C_omega(x)=0 for x<1.",
            "Definition only.",
        ),
        entry(
            "jek_03_residue_main_term",
            "exact_identity",
            "available_exact",
            "A_omega=1/((1+omega)zeta(1+2omega)); "
            "A_omega*x^(1+omega) is the residue main term of C_omega.",
            "A residue main term is not a pointwise approximation with a "
            "signed remainder.",
        ),
        entry(
            "jek_04_error_definition",
            "definition",
            "available_exact",
            "E_omega(x)=C_omega(x)-A_omega*x^(1+omega).",
            "Definition of the signed summatory discrepancy.",
        ),
        entry(
            "jek_05_weight_derivative",
            "exact_identity",
            "available_exact",
            "W_(omega,k)(u)=-u*g_(omega,k)'(u); for k>=1, "
            "W_(omega,k)=g_(omega,k-1)+(1/2)g_(omega,k).",
            "Uses Suzuki's weight recursion and the level-zero weight.",
        ),
        entry(
            "jek_06_abel_formula",
            "exact_identity",
            "available_exact",
            "H_(omega,k)(x)=x^(-1/2)*integral_0^x "
            "C_omega(t)W_(omega,k)(t/x)dt/t.",
            "Stieltjes/Abel summation only; no sign follows because W is "
            "signed.",
        ),
        entry(
            "jek_07_archimedean_zero",
            "exact_identity",
            "available_exact",
            "integral_0^1 g_(omega,k)(u)u^omega du=0.",
            "The zero is the factor s-omega-1 at s=1+omega in the "
            "archimedean Mellin transform.",
        ),
        entry(
            "jek_08_derivative_moment_zero",
            "exact_identity",
            "available_exact",
            "integral_0^1 W_(omega,k)(u)u^omega du=0.",
            "Integration by parts is valid because the endpoint products "
            "vanish.",
        ),
        entry(
            "jek_09_main_term_annihilation",
            "exact_lemma",
            "available_exact",
            "integral_0^1 A_omega*(x*u)^(1+omega)"
            "*W_(omega,k)(u)du/u=0.",
            "The positive residue main term contributes exactly zero.",
        ),
        entry(
            "jek_10_error_kernel",
            "exact_lemma",
            "available_exact",
            "H_(omega,k)(x)=x^(-1/2)*integral_0^1 "
            "E_omega(x*u)W_(omega,k)(u)du/u.",
            "Exact reduction to a signed error convolution; it proves no "
            "positivity.",
        ),
        entry(
            "jek_11_power_sum_decomposition",
            "exact_identity",
            "available_exact",
            "E_omega(x)=sum_(d<=x)mu(d)d^(-omega)R_omega(x/d)"
            "-x^(1+omega)/(1+omega)*sum_(d>x)"
            "mu(d)d^(-(1+2omega)), where R_omega(y)="
            "sum_(m<=y)m^omega-y^(1+omega)/(1+omega).",
            "Both terms require Mobius cancellation and have no fixed sign.",
        ),
        entry(
            "jek_12_elementary_error_bound",
            "exact_bound",
            "available_exact",
            "E_omega(x)=O_omega(x^(1-omega)) for 0<omega<1/2.",
            "Absolute-value estimate only; it discards all Mobius "
            "cancellation.",
        ),
        entry(
            "jek_13_absolute_kernel_bound",
            "exact_bound",
            "available_exact",
            "The elementary E bound and endpoint size "
            "W_(omega,k)(u)=O_(omega,k)(u^(omega-1)) give "
            "H_(omega,k)(x)=O_(omega,k)"
            "(x^(1/2-omega)*(1+log x)).",
            "Far too large to establish eventual sign or the RH-conditional "
            "polylogarithmic scale.",
        ),
        entry(
            "jek_14_smoothing_barrier",
            "route_obstruction",
            "guard_validated",
            "Every finite k has the same near-zero power u^(omega-1), so "
            "logarithmic smoothing alone does not improve the exponent in "
            "the absolute-value bound.",
            "Smoothing may improve cancellation, but not through a naive "
            "absolute estimate.",
        ),
        entry(
            "jek_15_positive_main_term_rejected",
            "route_obstruction",
            "guard_validated",
            "The Jordan-totient residue main term cannot be used as a "
            "positive dominant term because the Suzuki kernel annihilates "
            "it exactly.",
            "Any proof claiming main-term domination fails this gate.",
        ),
        entry(
            "jek_16_open_cancellation_target",
            "open_arithmetic_gate",
            "open_target",
            "For one explicit omega_j->0 and selected k_j, prove that the "
            "signed convolution of E_(omega_j) with W_(omega_j,k_j) has one "
            "eventual sign, gaining beyond the "
            "O(x^(1/2-omega_j)log x) absolute barrier.",
            "This cancellation estimate is not proved.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_suzuki_jordan_error_kernel_reduction",
        "date": "2026-07-23",
        "status": "exact arithmetic reduction with one open cancellation gate",
        "proof_boundary": (
            "This artifact is not a proof of an eventual sign, RH, or "
            "Lambda<=0. It exactly removes the Jordan-totient residue main "
            "term and identifies the remaining signed Mobius-error "
            "convolution and its inadequate absolute bound."
        ),
        "rows": rows,
        "exact": {
            "coefficient_convolution": (
                "c_omega(n)=sum_(dm=n)mu(d)d^(-omega)m^omega"
            ),
            "main_coefficient": (
                "A_omega=1/((1+omega)zeta(1+2omega))"
            ),
            "weight_derivative": (
                "W_(omega,k)=-u*g_(omega,k)'"
            ),
            "abel_formula": (
                "H=x^(-1/2)integral_0^x C(t)W(t/x)dt/t"
            ),
            "moment_zero": "integral_0^1 u^omega W(u)du=0",
            "error_kernel": (
                "H=x^(-1/2)integral_0^1 E(xu)W(u)du/u"
            ),
            "elementary_error_bound": "E=O_omega(x^(1-omega))",
            "absolute_barrier": (
                "H=O_(omega,k)(x^(1/2-omega)(1+log x))"
            ),
        },
        "sources": [
            {
                "title": (
                    "On monotonicity of certain weighted summatory "
                    "functions associated with L-functions"
                ),
                "url": "https://arxiv.org/abs/1204.1823",
                "use": "coefficients, weights, recursion, and Mellin transform",
            },
            {
                "title": "Cofinal monotonicity hierarchy",
                "url": (
                    "outputs/"
                    "jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.md"
                ),
                "use": "RH-strength eventual-sign handoff",
            },
        ],
        "audit": {
            "row_count": len(rows),
            "exact_reduction_count": 11,
            "route_obstruction_count": 2,
            "open_arithmetic_gate_count": 1,
            "main_term_survives": False,
            "eventual_sign_proved": False,
            "rh_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note(payload: dict) -> str:
    return """# Jensen-Window PF Suzuki Jordan-Error Kernel Reduction

Date: 2026-07-23

Status: exact arithmetic reduction with one open cancellation gate. This is
not a proof of an eventual sign, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_suzuki_jordan_error_kernel_reduction.json
python work/rh_compute/scripts/jensen_window_pf_suzuki_jordan_error_kernel_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_suzuki_jordan_error_kernel_reduction.py
```

## Cumulative Jordan Mass

For `0<omega<1/2`,

```text
c_omega(n)
 =sum_(d*m=n)mu(d)d^(-omega)m^omega,

C_omega(x)=sum_(n<=x)c_omega(n),

A_omega=1/((1+omega)*zeta(1+2omega)),

E_omega(x)=C_omega(x)-A_omega*x^(1+omega).           (JEK.1)
```

The constant `A_omega` is the residue main-term coefficient coming from
`zeta(s-omega)/zeta(s+omega)` at `s=1+omega`.

For Suzuki's level-`k` weight, set

```text
W_(omega,k)(u)=-u*g_(omega,k)'(u).
```

The logarithmic smoothing recursion gives

```text
W_(omega,k)
 =g_(omega,k-1)+(1/2)g_(omega,k).                    (JEK.2)
```

Stieltjes summation, using `g_(omega,k)(1)=0`, yields

```text
H_(omega,k)(x)
 =x^(-1/2)*integral_0^x
  C_omega(t)W_(omega,k)(t/x)dt/t.                   (JEK.3)
```

## Exact Main-Term Annihilation

The archimedean Mellin transform of `g_(omega,k)` contains the factor
`s-omega-1`. At `s=1+omega`,

```text
integral_0^1 g_(omega,k)(u)u^omega du=0.
```

Integration by parts, with vanishing endpoint products, then gives

```text
integral_0^1 W_(omega,k)(u)u^omega du=0.             (JEK.4)
```

Therefore the entire positive residue main term in (JEK.1) contributes
exactly zero to (JEK.3):

```text
H_(omega,k)(x)
 =x^(-1/2)*integral_0^1
  E_omega(x*u)W_(omega,k)(u)du/u.                   (JEK.5)
```

This is the sharp arithmetic handoff. Eventual positivity cannot come from
ordinary Jordan-totient main-term domination; it is wholly a signed
summatory-error cancellation problem.

## Exact Mobius Decomposition

Put

```text
P_omega(y)=sum_(m<=y)m^omega,
R_omega(y)=P_omega(y)-y^(1+omega)/(1+omega).
```

Using the coefficient convolution and the absolutely convergent identity
`sum mu(d)d^(-(1+2omega))=1/zeta(1+2omega)` gives

```text
E_omega(x)
 =sum_(d<=x)mu(d)d^(-omega)R_omega(x/d)
  -x^(1+omega)/(1+omega)
   *sum_(d>x)mu(d)d^(-(1+2omega)).                   (JEK.6)
```

Both terms require Mobius cancellation and neither has a fixed sign.

If `R_omega(y)=-y^(1+omega)/(1+omega)` is retained for `0<y<1`, the
tail is absorbed and (JEK.6) is equivalently the full generalized Muntz
identity

```text
E_omega(x)
 =sum_(d>=1)mu(d)d^(-omega)R_omega(x/d).
```

The finite natural-dilation energy of this full form is developed in
`outputs/jensen_window_pf_jordan_muntz_causal_energy_bridge.md`.

## Absolute-Bound Barrier

The elementary power-sum estimate

```text
R_omega(y)=O_omega(y^omega)
```

and absolute summation in (JEK.6) give

```text
E_omega(x)=O_omega(x^(1-omega)).                     (JEK.7)
```

Every finite smoothing level has

```text
W_(omega,k)(u)=O_(omega,k)(u^(omega-1))
```

near zero. Substitution into (JEK.5), split at `u=1/x`, yields only

```text
H_(omega,k)(x)
 =O_(omega,k)(x^(1/2-omega)*(1+log x)).              (JEK.8)
```

This is far above the RH-conditional scale `1` for `k=1` or
`(log x)^(k-1)` for higher `k`. Increasing `k` does not improve the power in
this naive absolute estimate because the near-zero exponent of every
smoothed weight is the same.

## Surviving Target

For one explicit sequence `omega_j->0` and convenient smoothing levels
`k_j`, prove directly that

```text
integral_0^1
 E_(omega_j)(x*u)W_(omega_j,k_j)(u)du/u
```

has one eventual sign after the `x^(-1/2)` normalization. A successful
argument must use cancellation or order structure beyond (JEK.7); positive
residue-main-term domination and absolute-value estimates are ruled out.

## Sources

- Masatoshi Suzuki, `On monotonicity of certain weighted summatory functions associated with L-functions`: https://arxiv.org/abs/1204.1823
- `outputs/jensen_window_pf_suzuki_cofinal_monotonicity_hierarchy.md`

## Proof Boundary

Equations (JEK.1)-(JEK.8) are exact identities or elementary bounds. They
identify the signed Mobius-error convolution that must be controlled and
reject two inadequate proof routes. They do not prove the required
cancellation, eventual sign, RH, or Lambda<=0.
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
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Suzuki Jordan-error kernel reduction: "
        f"{payload['audit']['row_count']} rows, "
        f"{payload['audit']['route_obstruction_count']} route guards, "
        f"{payload['audit']['open_arithmetic_gate_count']} open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
