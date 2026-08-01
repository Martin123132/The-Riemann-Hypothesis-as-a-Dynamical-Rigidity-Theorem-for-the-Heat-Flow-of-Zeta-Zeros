#!/usr/bin/env python3
"""Build the Mertens planar Abel and mixed-variation handoff."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_planar_abel_handoff.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_mertens_planar_abel_handoff.md"
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
            "pah_01_inherited_signed_criterion",
            "exact_equivalence",
            "available_exact",
            "Corollary 11.22Z.4 reduces the reciprocal-tail criterion "
            "to the dyadic signed off-diagonal O_(alpha,K).",
            "The high modes and uniformly summable diagonal remain closed.",
        ),
        row(
            "pah_02_rectangle",
            "exact_definition",
            "available_exact",
            "Use a_K=K+1<=n<=b_K=4K-2 and "
            "1<=h<=H_K=3K-2.",
            "This rectangle contains the triangular support n+h<=4K-1.",
        ),
        row(
            "pah_03_zero_extended_array",
            "exact_definition",
            "available_exact",
            "a_(K)(n,h)=mu(n)mu(n+h) when n+h<=4K-1 and is zero "
            "otherwise.",
            "The zero extension introduces no arithmetic terms outside "
            "the original off-diagonal.",
        ),
        row(
            "pah_04_zero_extended_kernel",
            "exact_definition",
            "available_exact",
            "Wtilde_(K)(n,h)=G_(alpha,K,R)(n,n+h) on n+h<=4K-1 "
            "and is zero otherwise.",
            "Also set Wtilde(b_K+1,h)=Wtilde(n,H_K+1)=0 for Abel "
            "boundary terms.",
        ),
        row(
            "pah_05_planar_prefix",
            "exact_definition",
            "available_exact",
            "S_K(x,H)=sum_(n=K+1)^x sum_(h=1)^H a_K(n,h).",
            "This is a joint base/shift prefix, not a collection of "
            "separately absolutized shift rows.",
        ),
        row(
            "pah_06_mixed_difference",
            "exact_definition",
            "available_exact",
            "Delta W(n,h)=W(n,h)-W(n+1,h)-W(n,h+1)+W(n+1,h+1).",
            "The n increment moves both ordinary kernel arguments because "
            "m=n+h.",
        ),
        row(
            "pah_07_double_abel",
            "exact_identity",
            "available_exact",
            "O_(alpha,K)=sum_(n=K+1)^(4K-2) "
            "sum_(h=1)^(3K-2) S_K(n,h) Delta Wtilde_K(n,h).",
            "This is finite two-dimensional Abel summation with forward "
            "zero boundary conditions.",
        ),
        row(
            "pah_08_mixed_variation",
            "exact_definition",
            "available_exact",
            "V_(alpha,K)=sum_(n,h)|Delta Wtilde_K(n,h)|.",
            "The absolute value is applied only to kernel curvature after "
            "the exact planar identity.",
        ),
        row(
            "pah_09_dirichlet_kernel",
            "exact_identity",
            "available_exact",
            "D_R(u)=sum_(r<=R)cos((2r-1)u)=sin(2Ru)/(2sin u).",
            "The quotient is interpreted continuously at its endpoints.",
        ),
        row(
            "pah_10_dirichlet_l1",
            "exact_bound",
            "available_exact",
            "L_R:=integral_0^pi |D_R(u)|du is at most "
            "(pi/2)(1+log(2R)).",
            "Use symmetry and min(R,pi/(4u)) on 0<u<=pi/2.",
        ),
        row(
            "pah_11_cc_coordinate",
            "exact_identity",
            "available_exact",
            "For n=K+i and m=n+h=K+j in the current block, "
            "W=(2/N)[C_R(h*pi/N)-C_R((i+j)*pi/N)].",
            "The first cosine term is annihilated by the mixed difference.",
        ),
        row(
            "pah_12_cc_curvature",
            "exact_identity",
            "available_exact",
            "If i+h<=K-2 and s=2i+h, Delta W equals "
            "(2/N) integral_0^(pi/N) integral_0^(2pi/N) "
            "D_R(s*pi/N+u+v)dvdu.",
            "This follows from C_R''=-D_R and two finite differences.",
        ),
        row(
            "pah_13_cc_variation",
            "exact_bound",
            "available_exact",
            "The total current/current interior mixed variation is at "
            "most 4*pi*K*L_R/N^2.",
            "For each s there are at most K pairs; the pi/N integration "
            "strips tile with bounded overlap.",
        ),
        row(
            "pah_14_future_convexity",
            "exact_derivative",
            "available_exact",
            "b_K(x)=(2K)^alpha(4K-x)x^(-1-alpha) is decreasing and "
            "strictly convex on [2K,4K].",
            "Its first two derivatives have the required signs for "
            "0<alpha<1.",
        ),
        row(
            "pah_15_future_decrements",
            "exact_monotonicity",
            "available_exact",
            "With b_x=b_K(x), b_x=0 for x>=4K, and d_x=b_x-b_(x+1), "
            "one has d_x>=d_(x+1)>=0 and sum_(x>=2K)d_x=1.",
            "Convexity gives decreasing forward decrements and the sum "
            "telescopes from b_(2K)=1.",
        ),
        row(
            "pah_16_initial_decrement",
            "exact_bound",
            "available_exact",
            "d_(2K)<=(2+alpha)/(2K)<3/(2K).",
            "Integrate -b_K' over [2K,2K+1] and use convexity.",
        ),
        row(
            "pah_17_endpoint_column",
            "exact_monotonicity",
            "available_exact",
            "p_i=P_(K,R)(i,K) satisfies 0=p_0<p_1<...<p_K<=rho_K.",
            "This is the inherited positive nested-interval formula.",
        ),
        row(
            "pah_18_endpoint_strip",
            "exact_bound",
            "available_exact",
            "sum_(i=1)^(K-2)|P(i,K)-P(i,K-1)|<pi^2/(2N).",
            "The lower and upper sine-integral boundary strips are "
            "disjoint and S_R is positive.",
        ),
        row(
            "pah_19_second_argument_transition",
            "exact_bound",
            "available_exact",
            "Across m=2K-1, Delta W=P(i,K-1)-p_i-p_(i+1)d_(2K), "
            "and the total absolute contribution is below "
            "pi^2/(2N)+(3/2)rho_K.",
            "This keeps the current/current and current/future pieces "
            "joined; zeroing them separately would create a false O(1) "
            "boundary loss.",
        ),
        row(
            "pah_20_current_future_curvature",
            "exact_bound",
            "available_exact",
            "For current n=K+i and future m, "
            "Delta W=p_i d_m-p_(i+1)d_(m+1); its total absolute "
            "contribution is at most (5/2)rho_K.",
            "Use monotonicity of p_i and d_m before summing the two "
            "telescoping factors.",
        ),
        row(
            "pah_21_first_argument_transition",
            "exact_bound",
            "available_exact",
            "For n=2K-1, Delta W=p_(K-1)d_m-P(K,K)d_(m+1), and the "
            "whole transition row costs at most 2rho_K.",
            "Both decrement sums are at most one.",
        ),
        row(
            "pah_22_future_future_curvature",
            "exact_bound",
            "available_exact",
            "For 2K<=n<m, Delta W=P(K,K)[b_n(d_m-d_(m+1))"
            "+d_n d_(m+1)] is nonnegative and its total mass is at "
            "most 2rho_K.",
            "Both future sums telescope and b_n<=1.",
        ),
        row(
            "pah_23_total_mixed_variation",
            "exact_bound",
            "available_exact",
            "V_(alpha,K)<3*pi^2*(1+log(2R))/K.",
            "Add rows 13 and 19-22, then use row 10, "
            "rho_K<pi^2/(4K), and N>2K.",
        ),
        row(
            "pah_24_planar_maximum_bound",
            "exact_bound",
            "available_exact",
            "With M_K=max_(x,H)|S_K(x,H)|, "
            "|O_(alpha,K)|<=V_(alpha,K)M_K.",
            "This is the total-variation consequence of the exact double "
            "Abel identity.",
        ),
        row(
            "pah_25_curvature_energy",
            "exact_definition",
            "available_exact",
            "E_(alpha,K)=sum_(n,h)|Delta Wtilde_K(n,h)|"
            "|S_K(n,h)|^2.",
            "This weighted energy averages planar prefixes only where "
            "the kernel has mixed curvature.",
        ),
        row(
            "pah_26_curvature_cauchy",
            "exact_bound",
            "available_exact",
            "|O_(alpha,K)|^2<=V_(alpha,K)E_(alpha,K).",
            "Apply weighted Cauchy-Schwarz to row 7.",
        ),
        row(
            "pah_27_energy_sufficient_condition",
            "conditional_calibration",
            "conditional_only",
            "If E_(alpha,K)=O_epsilon(K^(1+epsilon)) on dyadic K "
            "for every epsilon>0, then the weighted off-diagonal series "
            "is absolutely summable for that alpha.",
            "Choose epsilon<2alpha in rows 23 and 26; this estimate is "
            "not proved.",
        ),
        row(
            "pah_28_maximum_sufficient_condition",
            "conditional_calibration",
            "conditional_only",
            "The stronger estimate M_K=O_epsilon(K^(1+epsilon)) for "
            "every epsilon>0 closes every fixed positive alpha.",
            "Choose epsilon<alpha in rows 23 and 24; this maximal estimate "
            "is not proved.",
        ),
        row(
            "pah_29_full_triangle_identity",
            "exact_identity",
            "available_exact",
            "At x=4K-2 and H=3K-2, "
            "S_K=(A_K^2-Q_K)/2 with "
            "A_K=sum_(K<n<4K)mu(n) and Q_K=sum_(K<n<4K)mu(n)^2.",
            "The terminal planar prefix contains every unordered pair "
            "in the block exactly once.",
        ),
        row(
            "pah_30_maximum_strength_guard",
            "proof_guard",
            "guard_active",
            "The maximal condition in row 28 already forces "
            "A_K=O_epsilon(K^(1/2+epsilon)) through row 29.",
            "It contains RH-scale block-Mertens input and must not be "
            "advertised as an unconditional shortcut.",
        ),
        row(
            "pah_31_axiswise_countermodel",
            "countermodel_guard",
            "guard_active",
            "For K=L^2, an L-block diagonal K by K zero-one array with "
            "L by L all-one blocks has every row and column prefix at "
            "most sqrt(K), but total rectangle mass K^(3/2).",
            "Separate square-root prefix bounds do not imply planar "
            "square-root cancellation, even abstractly.",
        ),
        row(
            "pah_32_literature_guard",
            "literature_guard",
            "guard_active",
            "Averaged Chowla gives an o(HX) terminal average and the "
            "all-interval nilsequence theorem gives scalar logarithmic "
            "decay; neither supplies the K^(1+epsilon) curvature energy "
            "or planar maximum above.",
            "The endpoint/maximal quantifiers and a full power of "
            "cancellation remain missing.",
        ),
        row(
            "pah_33_planar_vaughan_identity",
            "exact_identity",
            "available_exact",
            "For the triangular band operator T_(x,H), "
            "S_K=<u_I,T u_I>+<u_II,T u_II>"
            "-<u_I,T u_II>-<u_II,T u_I> when mu=-u_I+u_II.",
            "T need not be self-adjoint; both cross terms must remain in "
            "the signed expression.",
        ),
        row(
            "pah_34_open_planar_gate",
            "open_gate",
            "open",
            "Prove the curvature-energy estimate in row 27, or directly "
            "bound the signed planar pairing in row 7, for every member "
            "of one cofinal positive-alpha sequence.",
            "No such joint Mobius/Vaughan estimate is currently available.",
        ),
        row(
            "pah_35_finite_validation",
            "finite_validation",
            "validated_finite",
            "Independent finite checks reproduce the double Abel identity, "
            "all five mixed-variation pieces, the global bound, the "
            "curvature Cauchy inequality, and the countermodel.",
            "Finite checks audit the algebra and constants; they do not "
            "prove the all-scale arithmetic gate.",
        ),
        row(
            "pah_36_proof_boundary",
            "proof_guard",
            "guard_active",
            "The logarithmic mixed-variation theorem is exact, but it "
            "does not estimate its Mobius planar prefixes.",
            "No curvature-energy gain, full Burnol bound, RH, PF-infinity, "
            "or Lambda<=0 conclusion is supplied.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_mertens_planar_abel_handoff",
        "date": "2026-07-24",
        "status": (
            "exact planar Abel identity and logarithmic mixed-kernel "
            "variation with one open arithmetic energy gate"
        ),
        "parameters": {
            "alpha_range": "0<alpha<1",
            "dyadic_blocks": "K=2^j",
            "cutoff": "R_(alpha,K)=ceil(K^(1-alpha/2))",
            "ambient_size": "N=2K+1",
            "rectangle": "K+1<=n<=4K-2, 1<=h<=3K-2",
            "mixed_variation_bound": (
                "V_(alpha,K)<3*pi^2*(1+log(2R))/K"
            ),
            "conditional_energy_target": (
                "E_(alpha,K)=O_epsilon(K^(1+epsilon))"
            ),
        },
        "source_anchors": [
            "outputs/jensen_window_pf_mertens_shift_kernel_variation_handoff.md",
            "outputs/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.md",
            "https://doi.org/10.2140/ant.2015.9.2167",
            "https://doi.org/10.1017/fmp.2023.28",
            "https://doi.org/10.5802/aif.2401",
        ],
        "rows": rows,
        "audit": {
            "row_count": 36,
            "exact_reduction_count": 28,
            "conditional_calibration_count": 2,
            "literature_guard_count": 1,
            "proof_guard_count": 3,
            "finite_validation_count": 1,
            "open_planar_gate_count": 1,
            "double_abel_proved": True,
            "mixed_variation_bound_proved": True,
            "curvature_cauchy_proved": True,
            "axiswise_countermodel_proved": True,
            "curvature_energy_proved": False,
            "joint_mobius_gain_proved": False,
            "full_burnol_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return """# Jensen-Window PF Mertens Planar Abel Handoff

Date: 2026-07-24

Status: exact planar Abel identity and logarithmic mixed-kernel
variation with one open arithmetic energy gate.
This is not a proof of RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_planar_abel_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_planar_abel_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_abel_handoff.py
```

## Exact Planar Abel Identity

Retain the notation of Corollary 11.22Z.4. Thus

```text
N=2K+1,
R=ceil(K^(1-alpha/2)),
rho_K=pi^2*K/N^2,
O_(alpha,K)=sum_(K<n<m<4K)mu(n)mu(m)G(n,m).
```

Use the rectangle

```text
a_K:=K+1<=n<=b_K:=4K-2,
1<=h<=H_K:=3K-2.
```

Extend both the arithmetic array and the kernel by zero:

```text
a_K(n,h):=
  mu(n)mu(n+h),  n+h<=4K-1,
  0,              n+h>4K-1,

Wtilde_K(n,h):=
  G(n,n+h),       n+h<=4K-1,
  0,              n+h>4K-1.                  (PAH.1)
```

Also set `Wtilde_K(b_K+1,h)=Wtilde_K(n,H_K+1)=0`. Define

```text
S_K(x,H)
 :=sum_(n=K+1)^x sum_(h=1)^H a_K(n,h),        (PAH.2)

Delta Wtilde_K(n,h)
 :=Wtilde_K(n,h)-Wtilde_K(n+1,h)
   -Wtilde_K(n,h+1)+Wtilde_K(n+1,h+1).        (PAH.3)
```

Two finite Abel summations give the exact identity

```text
O_(alpha,K)
 =sum_(n=K+1)^(4K-2)sum_(h=1)^(3K-2)
   S_K(n,h)Delta Wtilde_K(n,h).               (PAH.4)
```

No absolute value over the base point or shift has been taken in
(PAH.4).

## Current-Current Curvature

Put

```text
C_R(u)=sum_(r=1)^R cos((2r-1)u)/(2r-1)^2,

D_R(u):=sum_(r=1)^R cos((2r-1)u)
       =sin(2Ru)/(2sin u),

L_R:=integral_0^pi |D_R(u)|du.                (PAH.5)
```

On `(0,pi/2]`,

```text
|D_R(u)|<=min(R,pi/(4u)).
```

Symmetry and a split at `u=pi/(4R)` therefore give

```text
L_R<=(pi/2)(1+log(2R)).                       (PAH.6)
```

For `n=K+i`, `m=n+h=K+j`, and `s=i+j=2i+h`, the
current-current formula is

```text
Wtilde_K(n,h)
 =(2/N)[C_R(h*pi/N)-C_R(s*pi/N)].
```

The first term disappears under the mixed difference. Whenever
`i+h<=K-2`,

```text
Delta Wtilde_K(n,h)
 =(2/N)integral_0^(pi/N)integral_0^(2pi/N)
   D_R(s*pi/N+u+v)dvdu.                       (PAH.7)
```

For each `s` there are at most `K` pairs `(i,h)`. For fixed `v`, the
`pi/N` strips indexed by `s` are disjoint. Hence the total interior
current-current contribution to mixed variation is at most

```text
4*pi*K*L_R/N^2.                               (PAH.8)
```

The logarithm in (PAH.6) is the ordinary `L1` cost of the finite
Dirichlet kernel. No arithmetic estimate has entered.

## Transition And Future Curvature

For real `2K<=x<=4K`, write

```text
b_K(x)=(2K)^alpha(4K-x)x^(-1-alpha).
```

It is decreasing and strictly convex because

```text
b_K'(x)
 =(2K)^alpha[
   -4K(1+alpha)x^(-2-alpha)+alpha*x^(-1-alpha)
  ]<0,

b_K''(x)
 =(2K)^alpha(1+alpha)x^(-3-alpha)
  [4K(2+alpha)-alpha*x]>0.                    (PAH.9)
```

Set `b_x=b_K(x)` on `[2K,4K]`, extend it by zero for `x>=4K`, and
put `d_x=b_x-b_(x+1)`. Then

```text
d_x>=d_(x+1)>=0,
sum_(x>=2K)d_x=1,
d_(2K)<=(2+alpha)/(2K)<3/(2K).               (PAH.10)
```

Also retain `p_i=P_(K,R)(i,K)` and `c=P_(K,R)(K,K)`. Then

```text
0=p_0<p_1<...<p_K=c<=rho_K.                  (PAH.11)
```

The sine-integral boundary strips are disjoint, so

```text
sum_(i=1)^(K-2)|P(i,K)-P(i,K-1)|
 <pi^2/(2N).                                  (PAH.12)
```

All remaining mixed cells fall into four explicit families.

At the second-argument transition `m=2K-1`,

```text
Delta W
 =P(i,K-1)-p_i-p_(i+1)d_(2K),

sum_i|Delta W|
 <pi^2/(2N)+(3/2)rho_K.                      (PAH.13)
```

For current `n=K+i` and future `m>=2K`,

```text
Delta W=p_i d_m-p_(i+1)d_(m+1),

sum_(i,m)|Delta W|<=(5/2)rho_K.              (PAH.14)
```

At the first-argument transition `n=2K-1`,

```text
Delta W=p_(K-1)d_m-c*d_(m+1),

sum_m|Delta W|<=2rho_K.                       (PAH.15)
```

Finally, for `2K<=n<m`,

```text
Delta W
 =c[b_n(d_m-d_(m+1))+d_n d_(m+1)]>=0,

sum_(n,m)Delta W<=2rho_K.                     (PAH.16)
```

The joined treatment at `m=2K` is essential. Extending the
current-current and current-future blocks separately by zero would
manufacture a false `O(1)` interface variation that is absent from
the actual kernel.

## Logarithmic Mixed Variation

Define

```text
V_(alpha,K)
 :=sum_(n=K+1)^(4K-2)sum_(h=1)^(3K-2)
   |Delta Wtilde_K(n,h)|.
```

Adding (PAH.8) and (PAH.13)-(PAH.16) gives

```text
V_(alpha,K)
 <=4*pi*K*L_R/N^2+pi^2/(2N)+8rho_K.
```

Using (PAH.6), `N>2K`, and `rho_K<pi^2/(4K)` yields the convenient
uniform bound

```text
V_(alpha,K)
 <3*pi^2*(1+log(2R))/K.                       (PAH.17)
```

Thus, with

```text
M_K:=max_(K+1<=x<=4K-2,1<=H<=3K-2)|S_K(x,H)|,
```

(PAH.4) implies

```text
|O_(alpha,K)|
 <3*pi^2*(1+log(2R))*M_K/K.                  (PAH.18)
```

This is the first joint base/shift Abel interface: separate
square-root estimates are no longer summed over the other axis.

## Curvature-Weighted Energy Gate

The maximal condition in (PAH.18) is stronger than necessary. Define

```text
E_(alpha,K)
 :=sum_(n,h)|Delta Wtilde_K(n,h)||S_K(n,h)|^2. (PAH.19)
```

Weighted Cauchy-Schwarz applied only after (PAH.4) gives

```text
|O_(alpha,K)|^2
 <=V_(alpha,K)E_(alpha,K).                    (PAH.20)
```

Consequently, for one fixed `alpha>0`, the conditional estimate

```text
E_(alpha,K)=O_epsilon(K^(1+epsilon))          (PAH.21)
```

on dyadic `K`, for every `epsilon>0`, implies absolute summability of
`sum_K K^(-alpha)O_(alpha,K)`: choose `epsilon<2alpha` in
(PAH.17) and (PAH.20).

The stronger planar maximum

```text
M_K=O_epsilon(K^(1+epsilon))                  (PAH.22)
```

also suffices, now by choosing `epsilon<alpha` in (PAH.18). Both
(PAH.21) and (PAH.22) have the square-root scale for a planar region
with `O(K^2)` terms. Neither estimate is proved.

The energy gate (PAH.21) is the narrower target: it samples prefixes
only against the mixed-curvature measure of the actual kernel.

## Strength And Countermodel Guards

The terminal prefix is

```text
S_K(4K-2,3K-2)
 =1/2[A_K^2-Q_K],

A_K:=sum_(K<n<4K)mu(n),
Q_K:=sum_(K<n<4K)mu(n)^2.                    (PAH.23)
```

Therefore the maximal condition (PAH.22) already implies
`A_K=O_epsilon(K^(1/2+epsilon))`. It contains RH-scale block-Mertens
input and is not an unconditional shortcut. The curvature energy
(PAH.21) does not isolate that single corner in the same way.

Separate row and column square-root prefixes still do not imply a
planar square-root estimate. Let `K=L^2` and form a `K` by `K`
zero-one matrix from `L` diagonal all-one blocks of size `L` by `L`.
Every row and column prefix is at most `L=sqrt(K)`, while the full
rectangle sum is

```text
L^3=K^(3/2).                                  (PAH.24)
```

This abstract countermodel blocks promotion from the two
one-dimensional interfaces of Corollary 11.22Z.4.

Matomaki-Radziwill-Tao averaged Chowla has terminal scale `o(HX)`.
It neither provides the varying two-parameter prefixes in (PAH.2)
nor a fixed power saving from `K^2` to `K^(1+epsilon)`. The
all-interval Mobius/nilsequence theorem supplies logarithmic decay
for scalar linear tests, not the quadratic curvature energy
(PAH.19). Both are therefore power- or quantifier-short here.

## Vaughan Form And Open Gate

Let `T_(x,H)` be the finite triangular band operator whose quadratic
form is `S_K(x,H)`. Inserting the finite Vaughan identity
`mu=-u_I+u_II` before any square gives

```text
S_K(x,H)
 =<u_I,T u_I>+<u_II,T u_II>
  -<u_I,T u_II>-<u_II,T u_I>.                (PAH.25)
```

The operator need not be self-adjoint, so its two cross terms must
both remain. This is a rational-frequency, large-sieve, or
Hilbert-valued Vaughan interface; it is not permission to estimate
the four terms separately.

The exact surviving criterion is

```text
sup_J sum_(K dyadic<=2^J)K^(-alpha)
 sum_(n,h)S_K(n,h)Delta Wtilde_K(n,h)
 <infinity                                    (PAH.26)
```

for every member of one fixed cofinal positive-`alpha` sequence.
The concrete sufficient theorem target is (PAH.21), or a direct
signed estimate of (PAH.26) that is weaker.

The planar Abel identity, logarithmic mixed-variation estimate, and
Vaughan expansion are exact. The curvature-energy estimate, weighted
joint Mobius gain, full Burnol bound, RH, PF-infinity, and
`Lambda<=0` remain open.

The separate finite float64 diagnostic

```text
outputs/jensen_window_pf_mertens_planar_curvature_energy_scout.md
work/rh_compute/results/jensen_window_pf_mertens_planar_curvature_energy_scout.json
python work/rh_compute/scripts/check_jensen_window_pf_mertens_planar_curvature_energy_scout.py
```

tests the actual-Mobius quantities through `K=1024`. It is a
falsification scout only and is not used in the proof of (PAH.17).

Primary sources:

- Matomaki, Radziwill, and Tao, *An Averaged Form of Chowla's
  Conjecture*:
  https://doi.org/10.2140/ant.2015.9.2167
- Matomaki, Shao, Tao, and Teravainen, *Higher Uniformity of
  Arithmetic Functions in Short Intervals I. All Intervals*:
  https://doi.org/10.1017/fmp.2023.28
- Green and Tao, *Quadratic Uniformity of the Mobius Function*,
  Lemma 4.1:
  https://doi.org/10.5802/aif.2401
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
        "built Mertens planar Abel handoff: "
        "36 rows, 28 exact reductions, 3 proof guards, 1 open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
