#!/usr/bin/env python3
"""Build the weighted fractional-part autocorrelation Gram bridge."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.md"
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
            "wfagb_01_parameters",
            "exact_definition",
            "available_exact",
            "alpha=2omega, a=(1-alpha)/2, and "
            "beta=(1+alpha)/2=1-a for 0<omega<1/2.",
            "Thus 0<alpha<1, 0<a<1/2, and 1/2<beta<1.",
        ),
        row(
            "wfagb_02_weighted_autocorrelation",
            "exact_definition",
            "available_exact",
            "A_alpha(lambda)=integral_(0,infinity) "
            "{x}{lambda*x}x^(alpha-2)dx.",
            "Absolute convergence follows from O(x^alpha) at zero "
            "and O(x^(alpha-2)) at infinity.",
        ),
        row(
            "wfagb_03_scaling_law",
            "exact_identity",
            "available_exact",
            "A_alpha(lambda)=lambda^(1-alpha)A_alpha(1/lambda).",
            "The substitution x=lambda*y in A_alpha(1/lambda).",
        ),
        row(
            "wfagb_04_log_stationary_kernel",
            "exact_identity",
            "available_exact",
            "C_alpha(u)=exp(-a*u)A_alpha(exp(u)) is real, positive, "
            "continuous, and even.",
            "Evenness is exactly the scaling law; pointwise positivity "
            "comes from the defining integral.",
        ),
        row(
            "wfagb_05_fractional_part_mellin",
            "exact_identity",
            "available_exact",
            "integral_(0,infinity){x}x^(s-1)dx=zeta(-s)/s "
            "for -1<Re(s)<0.",
            "Classical Mellin identity, obtained by unit-interval "
            "decomposition and continuation inside the convergence strip.",
        ),
        row(
            "wfagb_06_autocorrelation_mellin",
            "exact_identity",
            "available_exact",
            "integral_(0,infinity)A_alpha(lambda)lambda^(s-1)d(lambda)="
            "-zeta(-s)zeta(1-alpha+s)/[s(s+1-alpha)].",
            "Tonelli/Fubini is valid first on alpha-1<Re(s)<0; insert "
            "the fractional-part Mellin identity twice.",
        ),
        row(
            "wfagb_07_fourier_density",
            "exact_identity",
            "available_exact",
            "At s=-a+it, Phi_alpha(t)=|zeta(a+it)|^2/(a^2+t^2) "
            "and C_alpha(u)=(1/(2*pi))integral_R "
            "Phi_alpha(t)exp(-itu)dt.",
            "Mellin-to-Fourier substitution lambda=exp(u); the standard "
            "convexity bound makes Phi_alpha integrable.",
        ),
        row(
            "wfagb_08_positive_definiteness",
            "exact_implication",
            "available_exact",
            "Phi_alpha(t)>=0 makes C_alpha positive definite and gives "
            "|C_alpha(u)|<=C_alpha(0).",
            "Bochner/Fourier inversion and the two-by-two principal minor.",
        ),
        row(
            "wfagb_09_burnol_reciprocal_norm",
            "exact_identity",
            "available_exact",
            "For f_(alpha,N)(t)=sum_(d<=N)mu(d)d^(-alpha){1/(dt)}, "
            "integral t^(-alpha)|f_(alpha,N)(t)|^2dt equals "
            "integral x^(alpha-2)|sum_(d<=N)mu(d)d^(-alpha){x/d}|^2dx.",
            "The reciprocal substitution x=1/t.",
        ),
        row(
            "wfagb_10_cross_kernel",
            "exact_identity",
            "available_exact",
            "I_alpha(d,e)=integral x^(alpha-2){x/d}{x/e}dx="
            "d^(alpha-1)A_alpha(d/e)="
            "(de)^((alpha-1)/2)C_alpha(log(d/e)).",
            "Scale x=d*y and then apply the definition of C_alpha.",
        ),
        row(
            "wfagb_11_stationary_gram_form",
            "exact_identity",
            "available_exact",
            "The Burnol norm is Gamma_(alpha,N)=sum_(d,e<=N)"
            "mu(d)mu(e)(de)^(-beta)C_alpha(log(d/e)).",
            "Finite expansion of the reciprocal norm and the cross-kernel "
            "identity.",
        ),
        row(
            "wfagb_12_spectral_gram_form",
            "exact_identity",
            "available_exact",
            "Gamma_(alpha,N)=(1/(2*pi))integral_R Phi_alpha(t)"
            "|sum_(d<=N)mu(d)d^(-beta-it)|^2dt.",
            "Insert the Fourier representation of C_alpha into the finite "
            "Gram sum.",
        ),
        row(
            "wfagb_13_mellin_plancherel_match",
            "exact_equivalence",
            "available_exact",
            "The stationary Gram formula and the Burnol-Hardy "
            "Mellin-Plancherel formula are the same identity, with "
            "beta=1/2+omega and a=1/2-omega.",
            "This is an exact normalization cross-check, not a new "
            "uniform estimate.",
        ),
        row(
            "wfagb_14_diagonal_split",
            "exact_identity",
            "available_exact",
            "Gamma_(alpha,N)=Diag_(alpha,N)+Off_(alpha,N), where "
            "Diag=C_alpha(0)sum_(d<=N)mu(d)^2d^(-(1+alpha)).",
            "Separate d=e from d!=e in the real symmetric Gram sum.",
        ),
        row(
            "wfagb_15_uniform_diagonal_bound",
            "exact_inequality",
            "available_exact",
            "0<=Diag_(alpha,N)<=C_alpha(0)zeta(1+alpha) uniformly in N.",
            "Use mu(d)^2<=1 and absolute convergence at 1+alpha>1.",
        ),
        row(
            "wfagb_16_signed_off_diagonal",
            "exact_identity",
            "available_exact",
            "Off_(alpha,N)=2sum_(d<e<=N)mu(d)mu(e)(de)^(-beta)"
            "C_alpha(log(d/e)).",
            "All remaining difficulty is a signed multiplicative "
            "correlation, despite C_alpha being positive definite.",
        ),
        row(
            "wfagb_17_absolute_value_barrier",
            "proof_guard",
            "available_exact",
            "Positive definiteness only gives "
            "Gamma_(alpha,N)>=0; the bound |C_alpha(u)|<=C_alpha(0) "
            "gives merely Gamma_(alpha,N)=O(N^(1-alpha)).",
            "Since sum_(d<=N)d^(-beta)=O(N^(1-beta)) and "
            "2(1-beta)=1-alpha; this is not uniform.",
        ),
        row(
            "wfagb_18_source_scope",
            "literature_guard",
            "source_backed_guard",
            "The alpha=0 autocorrelation and Mellin formula are in "
            "Báez-Duarte et al.; the 0<alpha<1 weighted extension and "
            "its finite Mobius Gram reduction are derived here.",
            "Do not attribute the weighted uniform bound to the "
            "unweighted source.",
        ),
        row(
            "wfagb_19_nonpromotion_guard",
            "proof_guard",
            "guarded",
            "Kernel positivity, pointwise positivity, bounded diagonal, "
            "or an unweighted autocorrelation theorem cannot by themselves "
            "bound the signed off-diagonal uniformly.",
            "A uniform upper bound needs arithmetic cancellation, not only "
            "positive-definite kernel structure.",
        ),
        row(
            "wfagb_20_open_off_diagonal_gate",
            "open_target",
            "open_target",
            "Prove sup_N Gamma_(alpha,N)<infinity along one explicit "
            "cofinal alpha sequence by a weighted multiplicative "
            "large-sieve, bilinear Mobius, or equivalent estimate.",
            "This RH-strength off-diagonal gate remains open; no RH, "
            "PF-infinity, or Lambda<=0 conclusion is recorded.",
        ),
    ]
    return {
        "kind": (
            "jensen_window_pf_weighted_fractional_"
            "autocorrelation_gram_bridge"
        ),
        "date": "2026-07-23",
        "status": (
            "exact weighted autocorrelation Gram/spectral reduction "
            "with one open signed off-diagonal gate"
        ),
        "parameters": {
            "omega_range": "0<omega<1/2",
            "alpha": "2omega",
            "a": "(1-alpha)/2",
            "beta": "(1+alpha)/2",
        },
        "source_anchors": [
            "https://arxiv.org/abs/math/0306251",
            "https://arxiv.org/abs/math/0202166",
            "https://arxiv.org/abs/1305.4395",
        ],
        "rows": rows,
        "audit": {
            "row_count": 20,
            "exact_identity_count": 12,
            "normalization_crosscheck_count": 1,
            "proof_guard_count": 3,
            "literature_guard_count": 1,
            "open_signed_off_diagonal_gate_count": 1,
            "uniform_gram_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Weighted Fractional Autocorrelation Gram Bridge

Date: 2026-07-23

Status: exact weighted autocorrelation Gram/spectral reduction with one open
signed off-diagonal gate. This is not a proof of RH, PF-infinity, or
`Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.json
python work/rh_compute/scripts/jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.py
python work/rh_compute/scripts/check_jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.py
```

## Weighted Autocorrelation

Fix

```text
0 < omega < 1/2,
alpha = 2*omega,
a = (1-alpha)/2 = 1/2-omega,
beta = (1+alpha)/2 = 1/2+omega.
```

For `lambda>0`, define

```text
A_alpha(lambda)
  = integral_(0,infinity) {x}{lambda*x} x^(alpha-2) dx.       (WFAGB.1)
```

The integrand is `O(x^alpha)` at zero and `O(x^(alpha-2))` at
infinity, so (WFAGB.1) converges absolutely. Substitution in
`A_alpha(1/lambda)` gives

```text
A_alpha(lambda)
  = lambda^(1-alpha) A_alpha(1/lambda).                       (WFAGB.2)
```

The logarithmically stationary normalization is

```text
C_alpha(u) = exp(-a*u) A_alpha(exp(u)).                       (WFAGB.3)
```

Equation (WFAGB.2) says exactly that `C_alpha(-u)=C_alpha(u)`.

## Mellin And Fourier Forms

In `-1<Re(s)<0`,

```text
integral_(0,infinity) {x} x^(s-1) dx = zeta(-s)/s.            (WFAGB.4)
```

Using (WFAGB.4) twice and first working in the common absolute
convergence strip `alpha-1<Re(s)<0` gives

```text
integral_(0,infinity) A_alpha(lambda) lambda^(s-1) d(lambda)
  = -zeta(-s) zeta(1-alpha+s) / [s(s+1-alpha)].               (WFAGB.5)
```

On the central line `s=-a+it`, the two zeta factors are conjugates:

```text
Phi_alpha(t) = |zeta(a+it)|^2 / (a^2+t^2),                   (WFAGB.6)

C_alpha(u)
  = (1/(2*pi)) integral_R Phi_alpha(t) exp(-itu) dt.          (WFAGB.7)
```

The standard convexity bound for zeta makes `Phi_alpha` integrable
because `a>0`. In particular, `Phi_alpha>=0`, so `C_alpha` is positive
definite and

```text
|C_alpha(u)| <= C_alpha(0).                                  (WFAGB.8)
```

## Exact Burnol Gram Form

Set

```text
f_(alpha,N)(t)
  = sum_(d<=N) mu(d)d^(-alpha) {1/(d*t)}.
```

The reciprocal substitution `x=1/t` gives

```text
integral_(0,infinity) t^(-alpha)|f_(alpha,N)(t)|^2 dt
  = integral_(0,infinity) x^(alpha-2)
      |sum_(d<=N) mu(d)d^(-alpha){x/d}|^2 dx.                 (WFAGB.9)
```

For the cross-kernel,

```text
I_alpha(d,e)
  = integral_(0,infinity) x^(alpha-2){x/d}{x/e} dx
  = d^(alpha-1) A_alpha(d/e)
  = (de)^((alpha-1)/2) C_alpha(log(d/e)).                    (WFAGB.10)
```

Consequently, the common value in (WFAGB.9) is

```text
Gamma_(alpha,N)
  = sum_(d,e<=N) mu(d)mu(e)(de)^(-beta)
      C_alpha(log(d/e))                                      (WFAGB.11)

  = (1/(2*pi)) integral_R Phi_alpha(t)
      |sum_(d<=N) mu(d)d^(-beta-it)|^2 dt.                   (WFAGB.12)
```

This is exactly the Burnol-Hardy Mellin-Plancherel norm already used in
the fixed-shift bridge: `beta=1/2+omega` and `a=1/2-omega`. Thus the
reciprocal-cell, Mellin, and stationary-Gram formulations have now been
normalization-checked against one another.

## The Real Remaining Gate

Split (WFAGB.11) into

```text
Diag_(alpha,N)
  = C_alpha(0) sum_(d<=N) mu(d)^2 d^(-(1+alpha)),             (WFAGB.13)

Off_(alpha,N)
  = 2 sum_(d<e<=N) mu(d)mu(e)(de)^(-beta)
      C_alpha(log(d/e)).                                     (WFAGB.14)
```

The diagonal is already harmless:

```text
0 <= Diag_(alpha,N) <= C_alpha(0) zeta(1+alpha).              (WFAGB.15)
```

All growth is therefore a signed off-diagonal arithmetic question.
Positive definiteness gives a lower bound for the full Gram form, not the
uniform upper bound that the proof programme needs. Indeed, applying
(WFAGB.8) absolutely gives only

```text
Gamma_(alpha,N)
  <= C_alpha(0)(sum_(d<=N)d^(-beta))^2
  = O(N^(1-alpha)),                                          (WFAGB.16)
```

which diverges. Kernel positivity, pointwise positivity, the bounded
diagonal, and the unweighted autocorrelation theorem therefore cannot be
promoted to the missing estimate.

The exact open target is:

> **Weighted signed off-diagonal gate.** Prove
> `sup_N Gamma_(alpha,N)<infinity` along one explicit cofinal sequence
> `alpha->0+`, using a weighted multiplicative large-sieve, bilinear
> Mobius-cancellation theorem, or an equivalent estimate strong enough to
> control (WFAGB.14).

Equation (WFAGB.16) is only a nonuniform barrier estimate.
The weighted signed off-diagonal gate is open.

## Source Boundary

- [Báez-Duarte et al., *A general strong Nyman-Beurling criterion for the
  Riemann hypothesis*](https://arxiv.org/abs/math/0306251) supplies the
  `alpha=0` fractional-part autocorrelation/Mellin prototype.
- [Burnol, *A note on Nyman's equivalent formulation of the Riemann
  hypothesis*](https://arxiv.org/abs/math/0202166) supplies the weighted
  approximants and Hardy-space setting.
- [Báez-Duarte et al., later autocorrelation analysis](https://arxiv.org/abs/1305.4395)
  studies the unweighted autocorrelation further.

The `0<alpha<1` stationary kernel, its finite weighted Mobius Gram form,
and the diagonal/off-diagonal reduction are derived here. None of the
sources is cited as proving the open weighted uniform bound.
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
        "wrote weighted fractional autocorrelation Gram bridge: "
        f"{len(payload['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
