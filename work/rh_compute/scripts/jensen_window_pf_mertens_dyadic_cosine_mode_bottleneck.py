#!/usr/bin/env python3
"""Build the dyadic reciprocal-tail cosine-mode bottleneck reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.md"
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
            "mdcos_01_tail_definitions",
            "exact_definition",
            "available_exact",
            "For 0<alpha<1 set q_n=mu(n)n^(-1-alpha), "
            "r_N=sum_(n>N)q_n, and "
            "R_alpha=sum_(N>=1)N^alpha|r_N|^2.",
            "This is the reciprocal-Mobius tail coordinate from the "
            "checked anchor/logarithmic tail-energy reduction.",
        ),
        row(
            "mdcos_02_dyadic_block",
            "exact_definition",
            "available_exact",
            "For dyadic K define y_m^(K)=r_(K+m-1), 1<=m<=K, and "
            "B_K=sum_(N=K)^(2K-1)N^alpha|r_N|^2.",
            "The blocks [K,2K-1] partition the positive integers.",
        ),
        row(
            "mdcos_03_weight_comparison",
            "exact_inequality",
            "available_exact",
            "K^alpha||y^(K)||_2^2<=B_K<="
            "2^alpha K^alpha||y^(K)||_2^2, and "
            "R_alpha=sum_(K dyadic)B_K.",
            "Only the elementary comparison K<=N<2K is used.",
        ),
        row(
            "mdcos_04_tail_difference",
            "exact_identity",
            "available_exact",
            "For 1<=m<K, y_m^(K)-y_(m+1)^(K)=q_(K+m).",
            "The endpoint tail and q_(2K) add a constant to the whole "
            "block and therefore disappear under the first difference.",
        ),
        row(
            "mdcos_05_neumann_cosine_basis",
            "exact_definition",
            "available_exact",
            "Let phi_0(m)=K^(-1/2) and, for 1<=r<K, "
            "phi_r(m)=sqrt(2/K)cos(pi*r*(m-1/2)/K). "
            "These diagonalize D_K^*D_K with "
            "nu_r=4sin^2(pi*r/(2K)).",
            "D_K is the forward first-difference map on a K-vector. "
            "The basis and eigenvalues follow by direct trigonometric "
            "substitution, including both boundary rows.",
        ),
        row(
            "mdcos_06_parseval",
            "exact_identity",
            "available_exact",
            "With c_(K,r)=<y^(K),phi_r>, "
            "||y^(K)||_2^2=sum_(r=0)^(K-1)|c_(K,r)|^2.",
            "The discrete cosine vectors in row 5 are an orthonormal basis.",
        ),
        row(
            "mdcos_07_gradient_parseval",
            "exact_identity",
            "available_exact",
            "sum_(m=1)^(K-1)|q_(K+m)|^2="
            "sum_(r=1)^(K-1)nu_r|c_(K,r)|^2.",
            "Combine row 4 with ||D_K y||_2^2="
            "<D_K^*D_K y,y> and row 5.",
        ),
        row(
            "mdcos_08_constant_mode",
            "exact_identity",
            "available_exact",
            "c_(K,0)=K^(-1/2)sum_(m=1)^K r_(K+m-1)="
            "sqrt(K)r_(2K)+K^(-1/2)"
            "sum_(ell=1)^K ell*q_(K+ell).",
            "Every q_(K+ell) occurs in exactly ell tails in the block. "
            "This is the sole mode affected by adding an endpoint tail "
            "constant.",
        ),
        row(
            "mdcos_09_nonconstant_modes",
            "exact_identity",
            "available_exact",
            "For 1<=r<K, c_(K,r)="
            "[sqrt(2/K)/(2sin(pi*r/(2K)))]"
            "sum_(m=1)^(K-1)q_(K+m)sin(pi*r*m/K).",
            "Use c=nu_r^(-1)<D_K y,D_K phi_r> and the exact cosine "
            "difference. The endpoint tail and q_(2K) cancel.",
        ),
        row(
            "mdcos_10_high_mode_bound",
            "exact_inequality",
            "available_exact",
            "For 1<=R<K, K^alpha sum_(r=R)^(K-1)|c_(K,r)|^2"
            "<=K^(1-alpha)/(4R^2).",
            "Since sin x>=2x/pi on [0,pi/2], nu_R>=4R^2/K^2. "
            "Then use row 7 and |q_n|<=K^(-1-alpha) on the block.",
        ),
        row(
            "mdcos_11_summable_threshold",
            "exact_consequence",
            "available_exact",
            "If R_K=ceil(K^beta) with beta>(1-alpha)/2, then the "
            "dyadic high-mode series is finite. The universal choice "
            "R_K=ceil(sqrt(K)) gives the block bound (1/4)K^(-alpha).",
            "Insert R_K into row 10 and sum the resulting geometric "
            "series over dyadic K.",
        ),
        row(
            "mdcos_12_low_mode_equivalence",
            "exact_equivalence",
            "available_exact",
            "For R_K=ceil(sqrt(K)), R_alpha<infinity iff "
            "sum_(K dyadic)K^alpha"
            "sum_(0<=r<R_K)|c_(K,r)|^2<infinity.",
            "Rows 3 and 6 compare R_alpha with the complete cosine "
            "energy, while row 11 proves the omitted high modes summable "
            "without any Mobius cancellation.",
        ),
        row(
            "mdcos_13_cofinal_rh_criterion",
            "exact_equivalence",
            "available_exact",
            "For any fixed cofinal alpha_j->0 in (0,1), RH iff the "
            "low-mode series in row 12 is finite for every alpha_j.",
            "Compose row 12 with the checked cofinal R_alpha criterion. "
            "This is an equivalent arithmetic target, not a convergence "
            "proof.",
        ),
        row(
            "mdcos_14_dimension_reduction",
            "exact_consequence",
            "available_exact",
            "At scale K the unresolved tail energy is confined to "
            "ceil(sqrt(K)) smooth cosine coordinates instead of all K "
            "cutoffs; every higher coordinate is automatically summable.",
            "This is a deterministic spectral reduction. It does not "
            "estimate any of the surviving low coordinates.",
        ),
        row(
            "mdcos_15_lowest_oscillatory_mode",
            "exact_consequence",
            "available_exact",
            "The r=1 coordinate is a positive one-hump sine-weighted "
            "Mobius sum on (K,2K), amplified by a factor asymptotic to "
            "sqrt(2K)/pi.",
            "Apply row 9 with r=1. This mode is a smooth anchored block "
            "average, not a high-frequency oscillation.",
        ),
        row(
            "mdcos_16_davenport_linear_phase",
            "classical_theorem_step",
            "source_backed",
            "For every A>0, Davenport's theorem gives uniformly in theta "
            "the estimate |sum_(n<=X)mu(n)e(n theta)|"
            "<<_A X(log X)^(-A).",
            "Primary source: Davenport (1937). Green and Tao, Example 3 "
            "and Section 5, record the normalized uniform form and its "
            "ineffective constant.",
        ),
        row(
            "mdcos_17_davenport_partial_summation",
            "classical_theorem_step",
            "source_backed",
            "Partial summation gives |c_(K,0)|"
            "<<_(alpha,A)K^(1/2-alpha)(log K)^(-A) and, for 1<=r<K, "
            "|c_(K,r)|<<_(alpha,A)"
            "K^(1/2-alpha)/(r(log K)^A).",
            "For r>0 absorb the sine into the two phases "
            "theta=+/-r/(2K) before summation by parts. For r=0 use "
            "theta=0 to bound the infinite reciprocal tail.",
        ),
        row(
            "mdcos_18_davenport_power_deficit_guard",
            "proof_guard",
            "guard_validated",
            "Davenport's bounds yield at best a low-mode block majorant "
            "of order K^(1-alpha)(log K)^(-2A), which does not prove "
            "dyadic summability for any fixed 0<alpha<1.",
            "Arbitrarily many logarithms cannot absorb the positive power "
            "K^(1-alpha). This is a proof-method insufficiency statement, "
            "not a lower bound on the actual Mobius modes.",
        ),
        row(
            "mdcos_19_general_exponent_threshold",
            "conditional_consequence",
            "conditional_exact",
            "A uniform additive-twist estimate with sigma<1+alpha, "
            "sum_(n<=X)mu(n)e(n theta)<<X^sigma would give low-mode "
            "block energy <<K^(2sigma-1-alpha), and closes the dyadic "
            "series when sigma<(1+alpha)/2.",
            "This is the exact exponent bookkeeping obtained from rows "
            "8-9 by partial summation. It assumes the displayed "
            "unproved power estimate.",
        ),
        row(
            "mdcos_20_square_root_sufficient",
            "conditional_consequence",
            "conditional_exact",
            "Uniform square-root-scale cancellation "
            "<<_epsilon X^(1/2+epsilon), with epsilon<alpha/2, makes "
            "the low-mode block energy O(K^(-alpha+2epsilon)) and hence "
            "proves R_alpha<infinity.",
            "This is a sufficient conditional closure of row 12. No such "
            "uniform additive-twist theorem is proved here.",
        ),
        row(
            "mdcos_21_square_root_scope_guard",
            "proof_guard",
            "guard_validated",
            "The square-root-scale hypothesis in row 20 contains the "
            "theta=0 Mertens bound and is not asserted to follow from RH "
            "alone; it is used only to calibrate the missing exponent.",
            "A sufficient hypothesis must not be reported as an "
            "unconditional theorem or as an established equivalent form.",
        ),
        row(
            "mdcos_22_translation_average_guard",
            "proof_guard",
            "guard_validated",
            "Translation-averaged Mobius cancellation and almost-all "
            "short-interval estimates do not automatically control the "
            "fixed dyadic mean and lowest one-hump modes.",
            "The surviving coordinates are anchored at (K,2K), including "
            "the exceptional origin-linked tail mode; an averaging theorem "
            "needs an explicit de-averaging step.",
        ),
        row(
            "mdcos_23_full_burnol_guard",
            "proof_guard",
            "guard_validated",
            "The cosine reduction concerns R_alpha and the logarithmic "
            "stable-prefix component P, not the full Burnol energy Q.",
            "A cofinal full-Q theorem would still imply RH through the "
            "existing route; this reduction neither proves nor obstructs "
            "such a theorem.",
        ),
        row(
            "mdcos_24_literature_context",
            "literature_guard",
            "source_backed",
            "The source audit used Davenport's original uniform linear-"
            "phase estimate and Green-Tao's primary exposition, including "
            "their summation-by-parts and completion tools.",
            "Those sources support the logarithmic linear-phase estimate, "
            "not the square-root hypothesis or the low-mode convergence.",
        ),
        row(
            "mdcos_25_focused_source_search_guard",
            "literature_guard",
            "source_backed",
            "A focused search found nearby exponential-sum refinements but "
            "no unconditional result removing the theta=0 power deficit "
            "required by row 19 on a cofinal alpha sequence.",
            "This is a search report, not a novelty, priority, or "
            "exhaustiveness claim.",
        ),
        row(
            "mdcos_26_open_low_mode_gate",
            "open_target",
            "open_target",
            "Prove the low-mode series in row 12 for every member of one "
            "fixed cofinal alpha sequence, using Mobius-specific signed "
            "cancellation across the dyadic mean and smooth modes.",
            "This target is RH-equivalent by row 13. No low-mode gain, "
            "uniform full-Burnol bound, RH, PF-infinity, or Lambda<=0 is "
            "proved.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck",
        "date": "2026-07-23",
        "status": (
            "exact dyadic cosine-mode reduction with automatically "
            "summable high modes and one open RH-equivalent low-mode gate"
        ),
        "parameters": {
            "alpha_range": "0<alpha<1",
            "dyadic_blocks": "K=2^j",
            "universal_low_mode_cutoff": "R_K=ceil(sqrt(K))",
            "general_summable_cutoff": "R_K=ceil(K^beta), beta>(1-alpha)/2",
        },
        "source_anchors": [
            "https://doi.org/10.1093/qmath/os-8.1.313",
            "https://doi.org/10.5802/aif.2401",
            "https://arxiv.org/abs/math/0202166",
        ],
        "rows": rows,
        "audit": {
            "row_count": 26,
            "exact_reduction_count": 15,
            "classical_theorem_step_count": 2,
            "conditional_consequence_count": 2,
            "literature_guard_count": 2,
            "proof_guard_count": 4,
            "open_low_mode_gate_count": 1,
            "cosine_diagonalization_proved": True,
            "high_modes_automatically_summable": True,
            "low_mode_cofinal_rh_equivalence_proved": True,
            "davenport_closes_low_modes": False,
            "square_root_twist_bound_proved": False,
            "low_mode_gain_proved": False,
            "full_burnol_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Mertens Dyadic Cosine-Mode Bottleneck

Date: 2026-07-23

Status: exact dyadic cosine-mode reduction with automatically summable
high modes and one open RH-equivalent low-mode gate. This is not a proof
of RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.json
python work/rh_compute/scripts/jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.py
```

## Dyadic Tail Blocks

Fix `0<alpha<1` and put

```text
q_n:=mu(n)n^(-1-alpha),
r_N:=sum_(n>N)q_n,
R_alpha:=sum_(N>=1)N^alpha|r_N|^2.                    (MDCOS.1)
```

For a dyadic integer `K`, define the length-`K` tail vector and its
weighted block energy by

```text
y_m^(K):=r_(K+m-1),             1<=m<=K,
B_K:=sum_(N=K)^(2K-1)N^alpha|r_N|^2.                  (MDCOS.2)
```

The dyadic intervals partition the positive integers, and

```text
R_alpha=sum_(K dyadic)B_K,

K^alpha||y^(K)||_2^2
 <=B_K
 <=2^alpha K^alpha||y^(K)||_2^2.                      (MDCOS.3)
```

The adjacent-tail identity becomes

```text
y_m^(K)-y_(m+1)^(K)=q_(K+m),       1<=m<K.            (MDCOS.4)
```

Notice what (MDCOS.4) does cleanly: the future tail `r_(2K)` and the
coefficient `q_(2K)` add the same constant to every entry of `y^(K)`.
They therefore affect only the mean mode below.

## Exact Neumann Cosine Diagonalization

Let `D_K` be the first-difference map in (MDCOS.4). The orthonormal
Neumann cosine basis is

```text
phi_0(m):=K^(-1/2),

phi_r(m):=sqrt(2/K)
          cos(pi*r*(m-1/2)/K),       1<=r<K,           (MDCOS.5)

D_K^*D_K phi_r=nu_r phi_r,
nu_r:=4sin^2(pi*r/(2K)).
```

For

```text
c_(K,r):=<y^(K),phi_r>,
```

Parseval and gradient Parseval give the exact identities

```text
||y^(K)||_2^2
 =sum_(r=0)^(K-1)|c_(K,r)|^2,                          (MDCOS.6)

sum_(m=1)^(K-1)|q_(K+m)|^2
 =sum_(r=1)^(K-1)nu_r|c_(K,r)|^2.                     (MDCOS.7)
```

The mean coordinate retains the endpoint information:

```text
c_(K,0)
 =K^(-1/2)sum_(m=1)^K r_(K+m-1)

 =sqrt(K)r_(2K)
  +K^(-1/2)sum_(ell=1)^K ell*q_(K+ell).               (MDCOS.8)
```

Every nonconstant coordinate is a smooth sine-weighted Mobius sum:

```text
c_(K,r)
 =[sqrt(2/K)/(2sin(pi*r/(2K)))]
   sum_(m=1)^(K-1)
     mu(K+m)(K+m)^(-1-alpha)sin(pi*r*m/K),

1<=r<K.                                                (MDCOS.9)
```

This follows from
`c_(K,r)=nu_r^(-1)<D_K y^(K),D_K phi_r>`. In particular,
the endpoint tail cancels before any estimate is made.

## High Modes Are Free

For `1<=R<K`, monotonicity of `nu_r` and
`sin x>=2x/pi` on `[0,pi/2]` give

```text
nu_R>=4R^2/K^2.
```

Since `|mu(n)|<=1`,

```text
sum_(m=1)^(K-1)|q_(K+m)|^2<=K^(-1-2alpha).
```

Consequently,

```text
K^alpha sum_(r=R)^(K-1)|c_(K,r)|^2
 <=K^(1-alpha)/(4R^2).                                 (MDCOS.10)
```

Thus any cutoff

```text
R_K:=ceil(K^beta),       beta>(1-alpha)/2,             (MDCOS.11)
```

makes the omitted high-mode blocks dyadically summable. The universal
choice `R_K=ceil(sqrt(K))` gives

```text
K^alpha sum_(r>=R_K)|c_(K,r)|^2
 <=(1/4)K^(-alpha).
```

Combining this with (MDCOS.3) and (MDCOS.6) yields the exact reduction

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^alpha
  sum_(0<=r<ceil(sqrt(K)))|c_(K,r)|^2
 <infinity.                                             (MDCOS.12)
```

For any fixed cofinal sequence `alpha_j->0`, the previously checked
Mertens/tail criterion therefore gives

```text
RH
 iff
(MDCOS.12) holds for every alpha_j.                     (MDCOS.13)
```

This removes all but `O(sqrt(K))` smooth coordinates at scale `K`
without using cancellation. It does not estimate the remaining
coordinates.

## What The Lowest Modes Ask For

The first oscillatory coordinate is

```text
c_(K,1)
 =[sqrt(2/K)/(2sin(pi/(2K)))]
   sum_(m=1)^(K-1)q_(K+m)sin(pi*m/K).                  (MDCOS.14)
```

The sine window is positive on the whole block and the prefactor is
asymptotic to `sqrt(2K)/pi`. This is a one-hump anchored Mobius average,
not a rapidly oscillating mode. The zero mode is the block mean of the
tails.

Davenport's theorem states, for every `A>0`, uniformly in real `theta`,

```text
|sum_(n<=X)mu(n)e(n theta)|
 <<_A X(log X)^(-A).                                   (MDCOS.15)
```

Partial summation, with the sine in (MDCOS.9) split into the phases
`theta=+/-r/(2K)`, gives

```text
|c_(K,0)|
 <<_(alpha,A)K^(1/2-alpha)(log K)^(-A),

|c_(K,r)|
 <<_(alpha,A)
 K^(1/2-alpha)/(r(log K)^A),       1<=r<K.             (MDCOS.16)
```

Hence Davenport gives only

```text
K^alpha
sum_(0<=r<sqrt(K))|c_(K,r)|^2
 <<_(alpha,A)K^(1-alpha)(log K)^(-2A).                 (MDCOS.17)
```

The positive power `K^(1-alpha)` defeats every fixed logarithmic saving.
This does not prove that the actual low modes are large. It proves that
classical arbitrary-log linear-phase orthogonality does not close the
required dyadic series.

More generally, a hypothetical uniform estimate with `sigma<1+alpha`,

```text
|sum_(n<=X)mu(n)e(n theta)|<<X^sigma
```

would produce the block scale

```text
K^(2sigma-1-alpha).                                    (MDCOS.18)
```

It closes the dyadic sum for fixed `alpha` precisely at the level

```text
sigma<(1+alpha)/2.
```

In particular, the stronger unproved hypothesis

```text
sup_theta |sum_(n<=X)mu(n)e(n theta)|
 <<_epsilon X^(1/2+epsilon),

epsilon<alpha/2,                                       (MDCOS.19)
```

would make each low-mode block
`O(K^(-alpha+2epsilon))` and prove `R_alpha<infinity`.
Equation (MDCOS.19) is only a calibration of the missing exponent. It
contains the `theta=0` RH-scale Mertens estimate and is not claimed here
as a theorem or as a consequence of RH alone.

## Guards And Open Gate

Translation-averaged cancellation and almost-all short-interval theorems
do not automatically control the fixed dyadic mean and one-hump modes.
An explicit de-averaging theorem would be required. Likewise, this
reduction concerns `R_alpha` and the logarithmically summed
stable-prefix component `P`; it neither proves nor obstructs a cofinal
theorem for the full Burnol energy `Q`.

The primary linear-phase sources checked were:

- H. Davenport, *On Some Infinite Series Involving Arithmetical
  Functions (II)*:
  https://doi.org/10.1093/qmath/os-8.1.313
- B. Green and T. Tao, *Quadratic Uniformity of the Mobius Function*,
  especially Example 3, Section 5, and the summation-by-parts appendix:
  https://doi.org/10.5802/aif.2401

A focused search found later exponential-sum refinements but no
unconditional result removing the `theta=0` power deficit needed by
(MDCOS.18) on a cofinal alpha sequence. This is a search report, not a
novelty, priority, or exhaustiveness claim.

The sharpened open target is

```text
for every alpha in one fixed cofinal sequence alpha_j->0,

sum_(K dyadic)K^alpha
  sum_(0<=r<ceil(sqrt(K)))|c_(K,r)|^2
 <infinity.                                             (MDCOS.20)
```

This is RH-equivalent by (MDCOS.13). The exact spectral reduction proves
the high modes harmless and identifies the low smooth modes, but supplies
no Mobius-specific gain. RH, PF-infinity, and `Lambda <= 0` remain open.
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
        json.dumps(payload, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(), encoding="utf-8")
    print(
        "built Mertens dyadic cosine-mode bottleneck: "
        f"{len(payload['rows'])} rows, "
        f"{payload['audit']['exact_reduction_count']} exact reductions, "
        f"{payload['audit']['proof_guard_count']} proof guards, "
        f"{payload['audit']['open_low_mode_gate_count']} open low-mode gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
