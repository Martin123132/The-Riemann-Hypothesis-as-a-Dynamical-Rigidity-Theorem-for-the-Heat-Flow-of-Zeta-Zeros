#!/usr/bin/env python3
"""Build the ordinary-Mobius feature-Gram and vector Vaughan handoff."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.md"
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
            "movh_01_prior_equivalence",
            "exact_equivalence",
            "available_exact",
            "For each fixed 0<alpha<1, Lemma 11.22Y gives "
            "R_alpha<infinity iff "
            "sum_(K dyadic)K^(-2-alpha)H_K<infinity.",
            "This artifact begins from the independently checked "
            "ordinary-Mobius anchored criterion.",
        ),
        row(
            "movh_02_future_anchor",
            "exact_definition",
            "available_exact",
            "With s=1+alpha, beta_K=mu(2K)+"
            "sum_(2K<n<4K)((2K)/n)^s(4K-n)/(2K)mu(n).",
            "The future anchor is retained exactly rather than bounded "
            "or discarded.",
        ),
        row(
            "movh_03_anchor_coefficients",
            "exact_bound",
            "available_exact",
            "Every coefficient b_K(n)=((2K)/n)^s(4K-n)/(2K) in "
            "beta_K lies in [0,1].",
            "Bounded coefficients do not imply cancellation of beta_K.",
        ),
        row(
            "movh_04_feature_support",
            "exact_definition",
            "available_exact",
            "Let I_K={n:K<n<4K}; the feature vectors have K real "
            "coordinates indexed by 1<=t<=K.",
            "All arithmetic coefficients now lie in one compact "
            "three-to-one interval.",
        ),
        row(
            "movh_05_feature_vector",
            "exact_definition",
            "available_exact",
            "Set a_(K,t)(K+i)=1_(t<=i) for 1<=i<K and "
            "a_(K,t)(n)=b_K(n) for 2K<=n<4K.",
            "The current block is a nested suffix feature and the future "
            "anchor is the constant vector direction.",
        ),
        row(
            "movh_06_suffix_vector",
            "exact_identity",
            "available_exact",
            "B_(K,t)=sum_(n in I_K)mu(n)a_(K,t)(n)="
            "beta_K+sum_(i=t)^(K-1)mu(K+i).",
            "The endpoint t=K contains beta_K and an empty current "
            "suffix.",
        ),
        row(
            "movh_07_energy_norm",
            "exact_identity",
            "available_exact",
            "H_K=sum_(t=1)^K|B_(K,t)|^2=||B_K||_2^2.",
            "The open scalar suffix family is one finite-dimensional "
            "Hilbert norm.",
        ),
        row(
            "movh_08_vector_energy_criterion",
            "exact_equivalence",
            "available_exact",
            "R_alpha<infinity iff "
            "sum_(K dyadic)K^(-2-alpha)||B_K||_2^2<infinity.",
            "This is an equivalent RH-strength criterion for every fixed "
            "positive alpha in a cofinal sequence.",
        ),
        row(
            "movh_09_gram_kernel",
            "exact_definition",
            "available_exact",
            "Define L_K(n,m)=sum_(t=1)^K"
            "a_(K,t)(n)a_(K,t)(m).",
            "The kernel is a Gram kernel and is therefore positive "
            "semidefinite.",
        ),
        row(
            "movh_10_current_current_kernel",
            "exact_identity",
            "available_exact",
            "For n=K+i and m=K+j with 1<=i,j<K, "
            "L_K(n,m)=min(i,j).",
            "Count the coordinates t lying below both suffix endpoints.",
        ),
        row(
            "movh_11_current_future_kernel",
            "exact_identity",
            "available_exact",
            "For n=K+i and 2K<=m<4K, "
            "L_K(n,m)=i*b_K(m), symmetrically.",
            "The future anchor occupies the constant vector direction.",
        ),
        row(
            "movh_12_future_future_kernel",
            "exact_identity",
            "available_exact",
            "For 2K<=n,m<4K, "
            "L_K(n,m)=K*b_K(n)b_K(m).",
            "The future-future block has rank one.",
        ),
        row(
            "movh_13_kernel_envelope",
            "exact_bound",
            "available_exact",
            "For all n,m in I_K, 0<=L_K(n,m)<=K.",
            "After multiplying by K^(-2-alpha), every pair weight is at "
            "most K^(-1-alpha).",
        ),
        row(
            "movh_14_gram_expansion",
            "exact_identity",
            "available_exact",
            "H_K=sum_(n,m in I_K)mu(n)mu(m)L_K(n,m).",
            "The positive norm still contains a genuinely signed "
            "off-diagonal Mobius form.",
        ),
        row(
            "movh_15_diagonal",
            "exact_identity",
            "available_exact",
            "Delta_K=sum_(i=1)^(K-1)i*mu(K+i)^2+"
            "K*sum_(2K<=n<4K)b_K(n)^2mu(n)^2.",
            "This is the complete Gram diagonal, including the future "
            "anchor.",
        ),
        row(
            "movh_16_future_diagonal_envelope",
            "exact_bound",
            "available_exact",
            "Since b_K(n)<=(4K-n)/(2K), "
            "K*sum_(2K<=n<4K)b_K(n)^2<="
            "(2K+1)(4K+1)/12.",
            "The bound is an exact finite sum-of-squares envelope.",
        ),
        row(
            "movh_17_diagonal_bound",
            "exact_bound",
            "available_exact",
            "0<=Delta_K<=(14K^2+1)/12<=(5/4)K^2.",
            "The last inequality holds for every integer K>=1.",
        ),
        row(
            "movh_18_diagonal_summability",
            "exact_bound",
            "available_exact",
            "sum_(K dyadic)K^(-2-alpha)Delta_K<="
            "(5/4)sum_(K dyadic)K^(-alpha)<infinity.",
            "The diagonal is closed for each fixed alpha>0, not at "
            "alpha=0.",
        ),
        row(
            "movh_19_signed_offdiagonal",
            "exact_definition",
            "available_exact",
            "Let C_(vec,K)=sum_(n<m; n,m in I_K)"
            "mu(n)mu(m)L_K(n,m); then H_K=Delta_K+2C_(vec,K).",
            "No absolute values are introduced in the correlation sum.",
        ),
        row(
            "movh_20_signed_gram_criterion",
            "exact_equivalence",
            "available_exact",
            "R_alpha<infinity iff sup_J sum_(K dyadic<=2^J)"
            "K^(-2-alpha)C_(vec,K)<infinity.",
            "The summable diagonal and nonnegative partial energy make "
            "the upper boundedness criterion exact.",
        ),
        row(
            "movh_21_absolute_barrier",
            "proof_guard",
            "guard_validated",
            "The feature envelope gives "
            "||sum_n mu(n)a_K(n)||_2=O(K^(3/2)) and H_K=O(K^3), "
            "hence weighted scale O(K^(1-alpha)).",
            "Coefficientwise or rowwise absolute values lose one full "
            "power relative to the summable K^(-alpha) calibration.",
        ),
        row(
            "movh_22_vaughan_interval_split",
            "exact_definition",
            "available_exact",
            "Split I_K into X<n<=2X for X=K and X=2K; for each X set "
            "U_X=V_X=floor(sqrt(X)), with a zero feature at n=4K.",
            "Both pieces now have the standard finite Vaughan interval "
            "shape.",
        ),
        row(
            "movh_23_type_i_coefficients",
            "exact_definition",
            "available_exact",
            "Define A_X(d)=sum_(bc=d,b<=U_X,c<=V_X)mu(b)mu(c).",
            "These are the finite Type I coefficients in Green-Tao's "
            "Vaughan identity.",
        ),
        row(
            "movh_24_type_ii_coefficients",
            "exact_definition",
            "available_exact",
            "Define D_X(d)=sum_(c|d,c>V_X)mu(c).",
            "These are the finite Type II divisor coefficients.",
        ),
        row(
            "movh_25_type_i_vector",
            "exact_definition",
            "available_exact",
            "V_(I,K,X)=sum_(d<=U_X*V_X)A_X(d)"
            "sum_(X/d<w<=2X/d)a_K(dw).",
            "The inner sum is a deterministic vector and contains no "
            "collapsed Mertens target.",
        ),
        row(
            "movh_26_type_ii_vector",
            "exact_definition",
            "available_exact",
            "V_(II,K,X)=sum_(V_X<d<=2X/U_X)D_X(d)"
            "sum_(max(U_X,X/d)<w<=2X/d)mu(w)a_K(dw).",
            "The bilinear Mobius factor is retained inside the vector "
            "sum.",
        ),
        row(
            "movh_27_componentwise_vaughan",
            "source_backed_identity",
            "source_backed",
            "Green-Tao's finite Vaughan identity applied componentwise "
            "gives sum_(X<n<=2X)mu(n)a_K(n)="
            "-V_(I,K,X)+V_(II,K,X).",
            "Finite-dimensional componentwise use requires no "
            "vector-valued extension theorem.",
        ),
        row(
            "movh_28_full_vector_vaughan",
            "exact_identity",
            "available_exact",
            "With V_(r,K)=V_(r,K,K)+V_(r,K,2K), "
            "B_K=-V_(I,K)+V_(II,K).",
            "The signed Type I/II cancellation is preserved before "
            "taking a norm.",
        ),
        row(
            "movh_29_vector_vaughan_criterion",
            "exact_equivalence",
            "available_exact",
            "R_alpha<infinity iff sum_(K dyadic)K^(-2-alpha)"
            "||-V_(I,K)+V_(II,K)||_2^2<infinity.",
            "This is a noncircular vector-valued arithmetic handoff with "
            "bounded deterministic features.",
        ),
        row(
            "movh_30_separate_norm_guard",
            "proof_guard",
            "guard_validated",
            "Bounds for ||V_(I,K)||_2 and ||V_(II,K)||_2 separately are "
            "sufficient only if they retain the required power.",
            "They are not equivalent to the signed-difference criterion "
            "and may erase its decisive cancellation.",
        ),
        row(
            "movh_31_coordinatewise_guard",
            "proof_guard",
            "guard_validated",
            "Even a uniform coordinate estimate B_(K,t)=o(K) yields "
            "only H_K=o(K^3).",
            "The target needs average-square cancellation at essentially "
            "K^2 scale, not merely sublinear cancellation per suffix.",
        ),
        row(
            "movh_32_finite_validation",
            "finite_reproduction",
            "finite_check_passed",
            "Mobius and arbitrary finite-vector checks reproduce the "
            "features, explicit Gram blocks, diagonal bound, signed "
            "expansion, and both interval Vaughan vectors.",
            "Finite validation proves the identities only and supplies "
            "no asymptotic Mobius estimate.",
        ),
        row(
            "movh_33_all_interval_log_guard",
            "literature_guard",
            "source_backed",
            "Matomaki-Shao-Tao-Teravainen give H*log^(-A)X "
            "cancellation for Mobius in every interval when "
            "H>=X^(5/8+epsilon).",
            "Squaring this scalar bound over O(K) long suffix endpoints "
            "still gives the K^3*log^(-2A) scale; the theorem does not "
            "state the required K^(2+epsilon) vector square estimate.",
        ),
        row(
            "movh_34_open_vector_gate",
            "open_theorem_target",
            "open",
            "Prove the vector Vaughan criterion in row 29 for every "
            "member of one fixed cofinal sequence alpha_j->0 without "
            "assuming RH.",
            "No vector Type I/II power gain, full Burnol bound, RH, "
            "PF-infinity, or Lambda<=0 conclusion is supplied here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_mertens_ordinary_vector_vaughan_handoff",
        "date": "2026-07-23",
        "status": (
            "exact ordinary-Mobius feature Gram and pre-square vector "
            "Vaughan reduction with one open signed gate"
        ),
        "parameters": {
            "alpha_range": "0<alpha<1",
            "dyadic_blocks": "K=2^j",
            "feature_support": "K<n<4K",
            "feature_dimension": "K",
            "vaughan_intervals": "X=K and X=2K",
        },
        "source_anchors": [
            "outputs/jensen_window_pf_mertens_affine_tent_bridge_handoff.md",
            "https://doi.org/10.5802/aif.2401",
            "https://doi.org/10.1017/fmp.2023.28",
        ],
        "rows": rows,
        "audit": {
            "row_count": 34,
            "exact_reduction_count": 28,
            "literature_guard_count": 2,
            "proof_guard_count": 3,
            "finite_validation_count": 1,
            "open_signed_gate_count": 1,
            "ordinary_feature_vector_proved": True,
            "explicit_gram_kernel_proved": True,
            "diagonal_summable": True,
            "componentwise_vaughan_proved": True,
            "pre_square_signed_difference_preserved": True,
            "vector_type_i_ii_gain_proved": False,
            "full_burnol_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Mertens Ordinary-Vector Vaughan Handoff

Date: 2026-07-23

Status: exact ordinary-Mobius feature Gram and pre-square vector Vaughan
reduction with one open signed gate. This is not a proof of RH,
PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.py
```

## Input From The Affine-Tent Lemma

Fix `0<alpha<1`, put `s=1+alpha`, and let `K` run over dyadic
integers. Formal Lemma 11.22Y proves

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^(-2-alpha)H_K<infinity,               (MOVH.1)
```

where

```text
beta_K
 :=mu(2K)
   +sum_(2K<n<4K)
     ((2K)/n)^s*(4K-n)/(2K)*mu(n),

H_K
 :=sum_(t=1)^K
   |beta_K+sum_(i=t)^(K-1)mu(K+i)|^2.                (MOVH.2)
```

Write

```text
b_K(n):=((2K)/n)^s*(4K-n)/(2K),     2K<=n<4K.
```

Then `0<=b_K(n)<=1`.

## Bounded Feature Vector

For `K<n<4K` define the vector `a_K(n)` in `R^K` by

```text
a_(K,t)(K+i):=1_(t<=i),        1<=i<K,
a_(K,t)(n):=b_K(n),            2K<=n<4K,
                                      1<=t<=K.       (MOVH.3)
```

The anchored suffix vector is exactly

```text
B_K:=sum_(K<n<4K)mu(n)a_K(n),

B_(K,t)
 =beta_K+sum_(i=t)^(K-1)mu(K+i),

H_K=||B_K||_2^2.                                    (MOVH.4)
```

Thus (MOVH.1) becomes

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^(-2-alpha)||B_K||_2^2<infinity.      (MOVH.5)
```

All arithmetic coefficients in this formulation are ordinary Mobius
values multiplied by deterministic numbers in `[0,1]`.

## Explicit Gram Kernel

Let

```text
L_K(n,m):=<a_K(n),a_K(m)>.
```

For `1<=i,j<K` and `2K<=n,m<4K`,

```text
L_K(K+i,K+j)=min(i,j),
L_K(K+i,n)=i*b_K(n),
L_K(n,m)=K*b_K(n)b_K(m).                             (MOVH.6)
```

Consequently `L_K` is positive semidefinite and
`0<=L_K(n,m)<=K`. The exact expansion is

```text
H_K
 =sum_(n,m in (K,4K))mu(n)mu(m)L_K(n,m)
 =Delta_K+2C_(vec,K),                                (MOVH.7)
```

where

```text
Delta_K
 :=sum_(i=1)^(K-1)i*mu(K+i)^2
   +K*sum_(2K<=n<4K)b_K(n)^2mu(n)^2,

C_(vec,K)
 :=sum_(K<n<m<4K)mu(n)mu(m)L_K(n,m).                 (MOVH.8)
```

Since `b_K(n)<=(4K-n)/(2K)`,

```text
K*sum_(2K<=n<4K)b_K(n)^2
 <=K*sum_(r=1)^(2K)(r/(2K))^2
 =(2K+1)(4K+1)/12.
```

Therefore

```text
0<=Delta_K
 <=K(K-1)/2+(2K+1)(4K+1)/12
 =(14K^2+1)/12
 <=(5/4)K^2.                                        (MOVH.9)
```

The normalized diagonal is dyadically summable for every fixed
`alpha>0`. Hence

```text
R_alpha<infinity
 iff
sup_J sum_(K dyadic<=2^J)
 K^(-2-alpha)C_(vec,K)<infinity.                    (MOVH.10)
```

This is still a signed criterion. The envelope
`||a_K(n)||_2<=sqrt(K)` and `O(K)` supported coefficients give only
`||B_K||_2=O(K^(3/2))`, hence `H_K=O(K^3)` and normalized scale
`O(K^(1-alpha))`. Absolute values therefore lose one full power.

## Vaughan Before Squaring

The vector identity admits a cleaner handoff than the scalar
off-diagonal expansion. Split the support into the standard intervals
`X<n<=2X` for `X=K,2K`. Put

```text
U_X=V_X=floor(sqrt(X)),

A_X(d)
 :=sum_(bc=d, b<=U_X, c<=V_X)mu(b)mu(c),

D_X(d)
 :=sum_(c|d, c>V_X)mu(c).                            (MOVH.11)
```

Define the finite vectors

```text
V_(I,K,X)
 :=sum_(d<=U_X*V_X)A_X(d)
   sum_(X/d<w<=2X/d)a_K(dw),

V_(II,K,X)
 :=sum_(V_X<d<=2X/U_X)D_X(d)
   sum_(max(U_X,X/d)<w<=2X/d)mu(w)a_K(dw).           (MOVH.12)
```

The inequalities in the inner sums mean integer ranges. Green and Tao's
finite Vaughan identity applies separately to every one of the `K`
coordinates, so no infinite-dimensional extension is being invoked:

```text
sum_(X<n<=2X)mu(n)a_K(n)
 =-V_(I,K,X)+V_(II,K,X).                             (MOVH.13)
```

Set

```text
V_(r,K):=V_(r,K,K)+V_(r,K,2K),       r in {I,II}.
```

The feature at `n=4K` is zero, and the two intervals cover all nonzero
features. Therefore

```text
B_K=-V_(I,K)+V_(II,K),                               (MOVH.14)
```

and the live theorem is exactly

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^(-2-alpha)
 ||-V_(I,K)+V_(II,K)||_2^2
 <infinity.                                          (MOVH.15)
```

This is Vaughan before squaring: the signed Type I/II cancellation
survives inside the norm. Separate bounds for the two vector norms may
be useful if they save the required power, but they are not equivalent
to (MOVH.15). Even a uniform coordinate estimate `B_(K,t)=o(K)` gives
only `H_K=o(K^3)`, still above the required essentially quadratic scale.

## All-Interval Logarithmic Guard

Matomaki, Shao, Tao, and Teravainen prove, in every interval of length
`H>=X^(5/8+epsilon)`, a Mobius bound of size
`H*log^(-A)X` for arbitrary fixed `A`. This is much stronger than an
unquantified `o(H)` statement, but it is still scalar. Applying it
coordinatewise to the `O(K)` long suffixes and squaring gives at best

```text
K^3*log^(-2A)K,                                     (MOVH.16)
```

up to harmless shorter-suffix and smooth-anchor terms. No logarithmic
choice of `A` turns this into `K^(2+epsilon)`. The cited theorem does
not give the vector correlation across nested suffix endpoints required
by (MOVH.15).

## Open Gate

For every member of one fixed cofinal sequence `alpha_j->0`, prove

```text
sum_(K dyadic)K^(-2-alpha_j)
 ||-V_(I,K)+V_(II,K)||_2^2
 <infinity                                            (MOVH.17)
```

without assuming RH. A viable estimate must save essentially one power
over coefficientwise control while preserving the signed Type I/II
difference and the vector correlation across suffix endpoints. No such
estimate, full Burnol bound, RH, PF-infinity, or `Lambda<=0` conclusion
is supplied here.

Primary source:

- Green and Tao, *Quadratic Uniformity of the Mobius Function*,
  Lemma 4.1:
  https://doi.org/10.5802/aif.2401
- Matomaki, Shao, Tao, and Teravainen, *Higher Uniformity of
  Arithmetic Functions in Short Intervals I. All Intervals*:
  https://doi.org/10.1017/fmp.2023.28
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
        "built Mertens ordinary-vector Vaughan handoff: "
        f"{len(payload['rows'])} rows, "
        f"{payload['audit']['exact_reduction_count']} exact reductions, "
        f"{payload['audit']['proof_guard_count']} proof guards, "
        f"{payload['audit']['open_signed_gate_count']} open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
