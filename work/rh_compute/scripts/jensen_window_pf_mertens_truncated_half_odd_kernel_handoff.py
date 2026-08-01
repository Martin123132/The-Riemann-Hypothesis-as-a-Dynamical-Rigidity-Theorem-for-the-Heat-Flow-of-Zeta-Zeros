#!/usr/bin/env python3
"""Build the truncated half-odd kernel and signed Vaughan handoff."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.md"
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
            "thok_01_reduced_input",
            "exact_equivalence",
            "available_exact",
            "Corollary 11.22Z.2 reduces R_alpha<infinity to the "
            "low half-odd energy with "
            "R_(alpha,K)=ceil(K^(1-alpha/2)).",
            "The complementary high modes are already summable and "
            "are not reopened.",
        ),
        row(
            "thok_02_notation",
            "exact_definition",
            "available_exact",
            "Put N=2K+1, q_r=2r-1, theta_r=q_r*pi/N, and "
            "phi_r(i)=2*N^(-1/2)sin(i*theta_r).",
            "For 1<=r<=K these are the mixed-boundary orthonormal "
            "eigenvectors from Corollary 11.22Z.1.",
        ),
        row(
            "thok_03_kernel",
            "exact_definition",
            "available_exact",
            "P_(K,R)(i,j):=sum_(r<=R)"
            "phi_r(i)phi_r(j)/q_r^2.",
            "The cutoff satisfies 1<=R<=K.",
        ),
        row(
            "thok_04_quadratic_identity",
            "exact_identity",
            "available_exact",
            "sum_(r<=R)|<x,phi_r>|^2/q_r^2="
            "x^T P_(K,R)x.",
            "No absolute value or Type I/II separation is introduced.",
        ),
        row(
            "thok_05_cosine_sum",
            "exact_definition",
            "available_exact",
            "C_R(u):=sum_(r<=R)cos(q_r*u)/q_r^2.",
            "This finite trigonometric polynomial is even.",
        ),
        row(
            "thok_06_cosine_kernel",
            "exact_identity",
            "available_exact",
            "P_(K,R)(i,j)=2/N*[C_R((i-j)pi/N)-"
            "C_R((i+j)pi/N)].",
            "This is product-to-sum applied before any estimate.",
        ),
        row(
            "thok_07_sine_primitive",
            "exact_definition",
            "available_exact",
            "S_R(u):=sum_(r<=R)sin(q_r*u)/q_r=-C_R'(u).",
            "The identity is termwise finite differentiation.",
        ),
        row(
            "thok_08_odd_dirichlet",
            "exact_identity",
            "available_exact",
            "S_R'(u)=sum_(r<=R)cos(q_r*u)="
            "sin(2Ru)/(2sin u), with continuous values at multiples "
            "of pi.",
            "This is the finite odd Dirichlet identity.",
        ),
        row(
            "thok_09_integral_kernel",
            "exact_identity",
            "available_exact",
            "P_(K,R)(i,j)=2/N*integral_"
            "(|i-j|pi/N)^((i+j)pi/N)S_R(u)du.",
            "Both endpoints lie in [0,pi) for 1<=i,j<=K.",
        ),
        row(
            "thok_10_sine_symmetry",
            "exact_identity",
            "available_exact",
            "S_R(pi-u)=S_R(u) for 0<=u<=pi.",
            "Every q_r is odd.",
        ),
        row(
            "thok_11_sine_positivity",
            "exact_positivity",
            "available_exact",
            "S_R(u)>0 for 0<u<pi.",
            "On (0,pi/2], pair successive half-waves in "
            "integral_0^u sin(2Rv)/(2sin v)dv using the decreasing "
            "weight 1/sin v; use symmetry on the other half.",
        ),
        row(
            "thok_12_entrywise_positivity",
            "exact_positivity",
            "available_exact",
            "P_(K,R)(i,j)>0 for all 1<=i,j<=K.",
            "The integration interval in row 9 has positive length.",
        ),
        row(
            "thok_13_psd",
            "exact_positivity",
            "available_exact",
            "P_(K,R) is positive semidefinite.",
            "It is a sum of positive rank-one matrices.",
        ),
        row(
            "thok_14_rank",
            "exact_identity",
            "available_exact",
            "rank P_(K,R)=R.",
            "The first R mixed-boundary eigenvectors are orthonormal.",
        ),
        row(
            "thok_15_nesting",
            "exact_comparison",
            "available_exact",
            "P_(K,R+1)-P_(K,R)="
            "phi_(R+1)phi_(R+1)^T/q_(R+1)^2 is rank-one PSD.",
            "Nesting is Loewner monotonicity, not necessarily "
            "entrywise monotonicity.",
        ),
        row(
            "thok_16_infinite_odd_cosine",
            "exact_identity",
            "available_exact",
            "sum_(r>=1)cos((2r-1)u)/(2r-1)^2="
            "pi^2/8-pi*u/4 for 0<=u<=pi.",
            "Subtract the even modes from the standard absolutely "
            "convergent cosine series.",
        ),
        row(
            "thok_17_min_kernel_completion",
            "exact_identity",
            "available_exact",
            "P_(K,infinity)(i,j)=pi^2*N^(-2)min(i,j).",
            "Insert row 16 into the cosine-kernel formula.",
        ),
        row(
            "thok_18_tail_psd",
            "exact_comparison",
            "available_exact",
            "0<=P_(K,R)<=pi^2*N^(-2)C_K in Loewner order, "
            "where C_K(i,j)=min(i,j).",
            "The omitted infinite tail is an absolutely convergent "
            "sum of rank-one PSD matrices.",
        ),
        row(
            "thok_19_full_lower_comparison",
            "exact_comparison",
            "available_exact",
            "4*N^(-2)C_K<=P_(K,K).",
            "Compare the common eigenvalues using "
            "sin(theta/2)>=theta/pi.",
        ),
        row(
            "thok_20_full_upper_comparison",
            "exact_comparison",
            "available_exact",
            "P_(K,K)<=pi^2*N^(-2)C_K.",
            "Compare the common eigenvalues using "
            "sin(theta/2)<=theta/2; this also follows from row 18.",
        ),
        row(
            "thok_21_trace",
            "exact_identity",
            "available_exact",
            "tr P_(K,R)=sum_(r<=R)q_r^(-2)<pi^2/8.",
            "The trace is uniformly bounded because the modes are "
            "orthonormal.",
        ),
        row(
            "thok_22_diagonal_envelope",
            "exact_bound",
            "available_exact",
            "0<P_(K,R)(i,i)<=pi^2*i/N^2.",
            "Take the ith diagonal entry in row 18.",
        ),
        row(
            "thok_23_pullback_maps",
            "exact_definition",
            "available_exact",
            "For K<n<4K, set tau_K(n)=n-K and w_K(n)=1 on "
            "K<n<2K, while tau_K(n)=K and w_K(n)=b_K(n) on "
            "2K<=n<4K.",
            "The future affine interval is collapsed only to the "
            "already-defined terminal anchor.",
        ),
        row(
            "thok_24_feature_pullback",
            "exact_identity",
            "available_exact",
            "g_(K,r)(n)=w_K(n)phi_r(tau_K(n)).",
            "This reproduces both feature cases in Corollary 11.22Z.1.",
        ),
        row(
            "thok_25_ordinary_kernel",
            "exact_definition",
            "available_exact",
            "G_(alpha,K,R)(n,m):=w_K(n)w_K(m)"
            "P_(K,R)(tau_K(n),tau_K(m)).",
            "G is a positive semidefinite and entrywise positive "
            "ordinary-support kernel.",
        ),
        row(
            "thok_26_ordinary_energy",
            "exact_identity",
            "available_exact",
            "L_(alpha,K):=sum_(r<=R)|X_(K,r)|^2/q_r^2="
            "sum_(n,m)mu(n)mu(m)G_(alpha,K,R)(n,m).",
            "The sums are finite on K<n,m<4K.",
        ),
        row(
            "thok_27_diagonal_offdiagonal",
            "exact_identity",
            "available_exact",
            "L_(alpha,K)=D_(alpha,K)+2O_(alpha,K), with "
            "D=sum_n mu(n)^2G(n,n) and "
            "O=sum_(n<m)mu(n)mu(m)G(n,m).",
            "All unresolved growth is isolated in the signed "
            "off-diagonal.",
        ),
        row(
            "thok_28_current_diagonal",
            "exact_bound",
            "available_exact",
            "The current-support diagonal is less than pi^2/8.",
            "Use mu^2<=1 and the trace bound.",
        ),
        row(
            "thok_29_future_tent_square",
            "exact_bound",
            "available_exact",
            "sum_(2K<=n<4K)b_K(n)^2<="
            "(2K+1)(4K+1)/(12K).",
            "Drop the factor (2K/n)^(1+alpha)<=1 and sum the "
            "squared affine tent exactly.",
        ),
        row(
            "thok_30_uniform_diagonal",
            "exact_bound",
            "available_exact",
            "0<=D_(alpha,K)<7pi^2/24 uniformly in alpha, K, and R.",
            "The future diagonal is at most "
            "pi^2(4K+1)/[12(2K+1)]<pi^2/6.",
        ),
        row(
            "thok_31_offdiagonal_equivalence",
            "exact_equivalence",
            "available_exact",
            "For fixed alpha, R_alpha<infinity iff "
            "sup_J sum_(K dyadic<=2^J)K^(-alpha)"
            "O_(alpha,K)<infinity.",
            "The high modes and the uniformly bounded per-scale "
            "diagonal are both dyadically summable.",
        ),
        row(
            "thok_32_endpoint_formula",
            "exact_identity",
            "available_exact",
            "p_i:=P_(K,R)(i,K)=4/N sum_(r<=R)"
            "(-1)^(r-1)sin(i*theta_r)cos(theta_r/2)/q_r^2.",
            "This retains the exact endpoint phase.",
        ),
        row(
            "thok_33_endpoint_integral",
            "exact_identity",
            "available_exact",
            "p_i=2/N integral_((K-i)pi/N)^((K+i)pi/N)"
            "S_R(u)du.",
            "This is row 9 with j=K.",
        ),
        row(
            "thok_34_endpoint_monotonicity",
            "exact_positivity",
            "available_exact",
            "0=p_0<p_1<...<p_K.",
            "Each successive endpoint interval strictly contains the "
            "previous one and S_R is positive.",
        ),
        row(
            "thok_35_endpoint_abel",
            "exact_identity",
            "available_exact",
            "If A_t=sum_(i=t)^(K-1)mu(K+i), then "
            "sum_(i=1)^(K-1)mu(K+i)p_i="
            "sum_(t=1)^(K-1)A_t(p_t-p_(t-1)).",
            "This is exact finite suffix Abel summation.",
        ),
        row(
            "thok_36_endpoint_abel_bound",
            "exact_bound",
            "available_exact",
            "The endpoint average has absolute value at most "
            "pi^2*K*N^(-2)max_t|A_t|.",
            "The positive Abel increments total p_(K-1)<=p_K and "
            "row 22 bounds p_K.",
        ),
        row(
            "thok_37_vaughan_coefficients",
            "exact_definition",
            "available_exact",
            "Collect the two interval Vaughan sums as "
            "T_(q,K,r)=sum_n u_(q,K)(n)g_(K,r)(n), q in {I,II}.",
            "The collected coefficients preserve every divisor and "
            "inner-variable incidence.",
        ),
        row(
            "thok_38_vaughan_identity",
            "source_backed_identity",
            "available_exact",
            "mu(n)=-u_(I,K)(n)+u_(II,K)(n) on the feature support.",
            "This is Green-Tao's finite Vaughan identity on "
            "(K,2K] and (2K,4K], with the zero feature at 4K harmless.",
        ),
        row(
            "thok_39_signed_gram_expansion",
            "exact_identity",
            "available_exact",
            "L=<u_I,G u_I>+<u_II,G u_II>-"
            "2<u_I,G u_II>.",
            "The Type I/II cross term remains inside the exact "
            "positive-semidefinite kernel form.",
        ),
        row(
            "thok_40_cross_term_guard",
            "proof_guard",
            "guard_validated",
            "Replacing row 39 by separate absolute Type I and Type II "
            "bounds is not an equivalent step.",
            "The missing arithmetic cancellation may occur in the "
            "cross term.",
        ),
        row(
            "thok_41_rank_guard",
            "proof_guard",
            "guard_validated",
            "For R<K, P_(K,R) has rank R and admits no positive "
            "Loewner lower bound by the full min kernel.",
            "Kernel positivity, trace control, and generic large-sieve "
            "bounds do not supply Mobius cancellation.",
        ),
        row(
            "thok_42_davenport_guard",
            "literature_guard",
            "source_backed",
            "Davenport-scale arbitrary logarithmic cancellation in "
            "the endpoint anchor or local suffixes does not close the "
            "K^(1-alpha) power deficit.",
            "The endpoint Abel identity changes the geometry, not the "
            "required power scale.",
        ),
        row(
            "thok_43_conditional_endpoint",
            "conditional_calibration",
            "conditional",
            "If both |beta_K| and max_t|A_t| are "
            "O_epsilon(K^(1/2+epsilon)), then the anchor-current cross "
            "term is O_epsilon(K^(2epsilon)).",
            "Both square-root estimates are unproved, and the "
            "current-current and future-future off-diagonals remain.",
        ),
        row(
            "thok_44_finite_reproduction",
            "finite_reproduction",
            "finite_check_passed",
            "Finite arbitrary-vector and Mobius checks reproduce the "
            "trigonometric kernel, positivity, nesting, min-kernel "
            "comparison, diagonal bound, Abel identity, and signed "
            "Vaughan Gram expansion.",
            "Finite checks validate identities only and supply no "
            "asymptotic Mobius estimate.",
        ),
        row(
            "thok_45_open_gate",
            "open_theorem_target",
            "open",
            "Prove the weighted low-kernel off-diagonal bound in row "
            "31, equivalently the joint signed Vaughan Gram estimate "
            "in row 39, for one cofinal sequence alpha_j->0.",
            "No such Mobius-specific bilinear estimate, full Burnol "
            "bound, RH, PF-infinity, or Lambda<=0 conclusion is "
            "supplied.",
        ),
    ]
    return {
        "kind": (
            "jensen_window_pf_mertens_"
            "truncated_half_odd_kernel_handoff"
        ),
        "date": "2026-07-24",
        "status": (
            "exact truncated half-odd kernel, uniformly summable "
            "diagonal, and signed off-diagonal/Vaughan reduction with "
            "one open arithmetic gate"
        ),
        "parameters": {
            "alpha_range": "0<alpha<1",
            "dyadic_blocks": "K=2^j",
            "cutoff": "R_(alpha,K)=ceil(K^(1-alpha/2))",
            "ambient_size": "N=2K+1",
            "mode_numbers": "q_r=2r-1",
            "uniform_diagonal_bound": "7pi^2/24",
        },
        "source_anchors": [
            "outputs/jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.md",
            "outputs/jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.md",
            "https://doi.org/10.5802/aif.2401",
            "https://doi.org/10.1093/qmath/os-8.1.313",
        ],
        "rows": rows,
        "audit": {
            "row_count": 45,
            "exact_reduction_count": 39,
            "literature_guard_count": 1,
            "proof_guard_count": 2,
            "conditional_calibration_count": 1,
            "finite_validation_count": 1,
            "open_signed_gate_count": 1,
            "kernel_formula_proved": True,
            "entrywise_positivity_proved": True,
            "min_kernel_completion_proved": True,
            "uniform_diagonal_bound_proved": True,
            "endpoint_abel_identity_proved": True,
            "signed_vaughan_gram_identity_proved": True,
            "offdiagonal_mobius_gain_proved": False,
            "full_burnol_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Mertens Truncated Half-Odd Kernel Handoff

Date: 2026-07-24

Status: exact truncated half-odd kernel, uniformly summable diagonal,
and signed off-diagonal/Vaughan reduction with one open arithmetic
gate. This is not a proof of RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.py
```

## Exact Truncated Kernel

Fix `0<alpha<1`, let `K` be dyadic, and put

```text
R=R_(alpha,K):=ceil(K^(1-alpha/2)),
N:=2K+1,
q_r:=2r-1,
theta_r:=q_r*pi/N,

phi_r(i):=2/sqrt(N)*sin(i*theta_r),       1<=i<=K.
                                                        (THOK.1)
```

Corollary 11.22Z.2 proves that the reciprocal-tail criterion is
equivalent to summability of

```text
L_(alpha,K)
 :=sum_(r=1)^R |X_(K,r)|^2/q_r^2.          (THOK.2)
```

Define

```text
P_(K,R)(i,j)
 :=sum_(r=1)^R phi_r(i)phi_r(j)/q_r^2.      (THOK.3)
```

For every `x in R^K`,

```text
sum_(r=1)^R |<x,phi_r>|^2/q_r^2
 =x^T P_(K,R)x.                             (THOK.4)
```

Put

```text
C_R(u):=sum_(r=1)^R cos(q_r*u)/q_r^2,
S_R(u):=sum_(r=1)^R sin(q_r*u)/q_r.
```

Product-to-sum gives the exact finite formula

```text
P_(K,R)(i,j)
 =2/N[
    C_R((i-j)pi/N)-C_R((i+j)pi/N)
   ].                                       (THOK.5)
```

Moreover,

```text
C_R'(u)=-S_R(u),

S_R'(u)
 =sum_(r=1)^R cos(q_r*u)
 =sin(2R*u)/(2sin u),                       (THOK.6)
```

where the last quotient is interpreted continuously. Hence

```text
P_(K,R)(i,j)
 =2/N integral_(|i-j|pi/N)^((i+j)pi/N)
      S_R(u)du.                              (THOK.7)
```

The sine polynomial is strictly positive:

```text
S_R(u)>0,                         0<u<pi.    (THOK.8)
```

Indeed, `S_R(pi-u)=S_R(u)`. On `(0,pi/2]`, write

```text
S_R(u)=integral_0^u sin(2R*v)/(2sin v)dv.
```

Pair each negative half-wave with the preceding positive half-wave.
The amplitude `1/sin v` is decreasing there, so every completed pair
is nonnegative; a terminal partial half-wave has the same property.
The first positive half-wave makes the inequality strict. Symmetry
handles `(pi/2,pi)`.

It follows from (THOK.7) that

```text
P_(K,R)(i,j)>0                    (THOK.9)
```

entrywise. Independently, (THOK.3) proves that `P_(K,R)` is positive
semidefinite, has rank `R`, and is nested in Loewner order:

```text
P_(K,R+1)-P_(K,R)
 =phi_(R+1)phi_(R+1)^T/q_(R+1)^2 >=0.       (THOK.10)
```

## Completion To The Min Kernel

For `0<=u<=pi`, the full odd cosine series is

```text
sum_(r>=1)cos((2r-1)u)/(2r-1)^2
 =pi^2/8-pi*u/4.                             (THOK.11)
```

Inserting (THOK.11) into (THOK.5) gives the exact infinite completion

```text
P_(K,infinity)(i,j)
 =pi^2/N^2*min(i,j).                         (THOK.12)
```

If `C_K(i,j)=min(i,j)`, the omitted tail is an absolutely convergent
sum of rank-one positive semidefinite matrices. Therefore

```text
0<=P_(K,R)<=pi^2/N^2*C_K                     (THOK.13)
```

in Loewner order. At the full finite rank `R=K`, comparison in the
common sine eigenbasis gives

```text
4/N^2*C_K
 <=P_(K,K)
 <=pi^2/N^2*C_K.                             (THOK.14)
```

The lower bound uses `sin(theta/2)>=theta/pi`; the upper bound uses
`sin(theta/2)<=theta/2`. No lower comparison with `C_K` is possible
when `R<K`, because `P_(K,R)` then has rank `R`.

Orthonormality also gives

```text
tr P_(K,R)
 =sum_(r=1)^R 1/q_r^2
 <pi^2/8,

0<P_(K,R)(i,i)<=pi^2*i/N^2.                  (THOK.15)
```

## Ordinary-Mobius Pullback

Retain

```text
b_K(n)
 :=((2K)/n)^(1+alpha)*(4K-n)/(2K),
                                      2K<=n<4K.
```

For `K<n<4K`, define

```text
tau_K(n):=
  n-K,                              K<n<2K,
  K,                               2K<=n<4K,

w_K(n):=
  1,                                K<n<2K,
  b_K(n),                           2K<=n<4K.
```

Then the ordinary-Mobius feature from Corollary 11.22Z.1 is exactly

```text
g_(K,r)(n)=w_K(n)phi_r(tau_K(n)).             (THOK.16)
```

Pull the truncated kernel back to the ordinary support:

```text
G_(alpha,K,R)(n,m)
 :=w_K(n)w_K(m)
   P_(K,R)(tau_K(n),tau_K(m)).                 (THOK.17)
```

This kernel is positive semidefinite and entrywise positive. The low
energy has the exact finite expansion

```text
L_(alpha,K)
 =sum_(K<n,m<4K)
   mu(n)mu(m)G_(alpha,K,R)(n,m)

 =D_(alpha,K)+2O_(alpha,K),                    (THOK.18)

D_(alpha,K)
 :=sum_(K<n<4K)mu(n)^2G_(alpha,K,R)(n,n),

O_(alpha,K)
 :=sum_(K<n<m<4K)mu(n)mu(m)G_(alpha,K,R)(n,m).
```

The current diagonal is less than `pi^2/8` by (THOK.15). On the future
interval,

```text
sum_(2K<=n<4K)b_K(n)^2
 <=sum_(ell=1)^(2K)(ell/(2K))^2
 =(2K+1)(4K+1)/(12K).                         (THOK.19)
```

Using `P_(K,R)(K,K)<=pi^2*K/N^2`, the future diagonal is less than
`pi^2/6`. Thus

```text
0<=D_(alpha,K)<7pi^2/24                       (THOK.20)
```

uniformly in `alpha`, `K`, and `R`. Its normalized dyadic series is
automatically finite. Combining (THOK.18)-(THOK.20) with Corollary
11.22Z.2 proves

```text
R_alpha<infinity
 iff
sup_J sum_(K dyadic<=2^J)
 K^(-alpha)O_(alpha,K)<infinity.              (THOK.21)
```

This is an exact low-mode, ordinary-Mobius, signed off-diagonal
criterion. Kernel positivity does not estimate that signed sum.

## Endpoint Abel Structure

Let

```text
p_i:=P_(K,R)(i,K),                    0<=i<=K,
p_0:=0.
```

The exact endpoint phase gives

```text
p_i
 =4/N sum_(r=1)^R
   (-1)^(r-1)sin(i*theta_r)cos(theta_r/2)/q_r^2.
                                                        (THOK.22)
```

The integral form is

```text
p_i
 =2/N integral_((K-i)pi/N)^((K+i)pi/N)
      S_R(u)du.                              (THOK.23)
```

Since `S_R>0`, these intervals are strictly nested:

```text
0=p_0<p_1<...<p_K.                            (THOK.24)
```

Put

```text
A_(K,t):=sum_(i=t)^(K-1)mu(K+i).
```

Finite suffix Abel summation now gives

```text
sum_(i=1)^(K-1)mu(K+i)p_i
 =sum_(t=1)^(K-1)
   A_(K,t)(p_t-p_(t-1)),                     (THOK.25)
```

and therefore

```text
|sum_(i=1)^(K-1)mu(K+i)p_i|
 <=pi^2*K/N^2*max_t|A_(K,t)|.                (THOK.26)
```

Thus the current-anchor cross term is an exact positive Abel average,
not an uncontrolled collection of modes. This geometry does not by
itself supply the missing power: Davenport-scale logarithmic bounds
remain one power short. If both `|beta_K|` and
`max_t|A_(K,t)|` were conditionally
`O_epsilon(K^(1/2+epsilon))`, the cross term would be
`O_epsilon(K^(2epsilon))`; both inputs are unproved and the other
off-diagonal blocks would still remain.

## Signed Vaughan Kernel Form

For `X in {K,2K}`, retain `U_X=V_X=floor(sqrt(X))` and the coefficients
`A_X(d),D_X(d)` from Lemma 11.22Z. Collect the finite Vaughan
incidences as

```text
u_(I,X)(n)
 :=1_(X<n<=2X)
   sum_(d|n,d<=U_X*V_X)A_X(d),

u_(II,X)(n)
 :=1_(X<n<=2X)
   sum_(d|n,V_X<d<=2X/U_X,n/d>U_X)
   D_X(d)mu(n/d).                              (THOK.27)
```

Sum the two disjoint intervals and write the resulting coefficients as
`u_(I,K),u_(II,K)`. The feature at `4K` is zero. Green and Tao's finite
identity gives

```text
mu(n)=-u_(I,K)(n)+u_(II,K)(n)                 (THOK.28)
```

on the feature support, and

```text
T_(q,K,r)
 =sum_(K<n<4K)u_(q,K)(n)g_(K,r)(n),
                                      q in {I,II}.
```

Consequently the entire truncated square is

```text
L_(alpha,K)
 =<u_I,G u_I>
  +<u_II,G u_II>
  -2<u_I,G u_II>.                             (THOK.29)
```

The cross term is part of the exact identity. Replacing (THOK.29) by
separate absolute Type I and Type II estimates is not equivalent and
may discard the only available cancellation.

## Open Gate

For every member of one fixed cofinal sequence `alpha_j->0`, prove the
weighted off-diagonal bound (THOK.21), equivalently the joint signed
Vaughan Gram estimate (THOK.29), without assuming RH. The diagonal,
high modes, kernel algebra, and endpoint Abel geometry are closed.
The Mobius-specific off-diagonal gain, full Burnol bound, RH,
PF-infinity, and `Lambda<=0` remain open.

Primary sources:

- Davenport, *On Some Infinite Series Involving Arithmetical
  Functions (II)*:
  https://doi.org/10.1093/qmath/os-8.1.313
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
        "built Mertens truncated half-odd kernel handoff: "
        f"{len(payload['rows'])} rows, "
        f"{payload['audit']['exact_reduction_count']} exact reductions, "
        f"{payload['audit']['proof_guard_count']} proof guards, "
        f"{payload['audit']['open_signed_gate_count']} open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
