#!/usr/bin/env python3
"""Build the Mertens weighted-prefix and affine-defect reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.md"
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
            "mwpad_01_weighted_prefix_definitions",
            "exact_definition",
            "available_exact",
            "For 0<alpha<1 define M(x)=sum_(n<=x)mu(n), "
            "A_alpha(x)=sum_(n<=x)mu(n)n^(-alpha), "
            "I_alpha=int M(x)^2*x^(-2-alpha)dx, and "
            "J_alpha=int A_alpha(x)^2*x^(alpha-2)dx.",
            "All integrals start at 1; beta=(1+alpha)/2 and "
            "lambda=(1-alpha)/2.",
        ),
        row(
            "mwpad_02_forward_abel",
            "exact_identity",
            "available_exact",
            "A_alpha(x)=x^(-alpha)M(x)+alpha*int_1^x "
            "M(u)u^(-alpha-1)du.",
            "Stieltjes integration by parts; no zero-free region is used.",
        ),
        row(
            "mwpad_03_inverse_abel",
            "exact_identity",
            "available_exact",
            "M(x)=x^alpha*A_alpha(x)-alpha*int_1^x "
            "A_alpha(u)u^(alpha-1)du.",
            "The inverse Stieltjes identity for mu(n)=n^alpha "
            "[mu(n)n^(-alpha)].",
        ),
        row(
            "mwpad_04_logarithmic_coordinates",
            "exact_definition",
            "available_exact",
            "For x=e^t put m(t)=e^(-beta*t)M(e^t) and "
            "a(t)=e^(-lambda*t)A_alpha(e^t); then "
            "||m||_2^2=I_alpha and ||a||_2^2=J_alpha.",
            "A unitary logarithmic change of variables on t>=0.",
        ),
        row(
            "mwpad_05_forward_volterra",
            "exact_identity",
            "available_exact",
            "a=m+alpha*V_lambda*m, where "
            "V_c f(t)=int_0^t exp(-c(t-u))f(u)du.",
            "Equation 2 in logarithmic coordinates.",
        ),
        row(
            "mwpad_06_inverse_volterra",
            "exact_identity",
            "available_exact",
            "m=a-alpha*V_beta*a.",
            "Equation 3 in logarithmic coordinates.",
        ),
        row(
            "mwpad_07_sharp_norm_sandwich",
            "exact_inequality",
            "available_exact",
            "I_alpha<=J_alpha<=((1+alpha)/(1-alpha))^2*I_alpha.",
            "The forward Fourier multiplier is "
            "(beta+it)/(lambda+it), whose modulus lies in "
            "[1,beta/lambda].",
        ),
        row(
            "mwpad_08_weighted_prefix_mellin",
            "exact_identity",
            "available_exact",
            "int_1^infinity A_alpha(x)x^(-s-1)dx="
            "1/[s*zeta(s+alpha)] initially for Re(s)>1-alpha.",
            "Interchange the finite-prefix step function and the integral.",
        ),
        row(
            "mwpad_09_mellin_plancherel",
            "classical_theorem_step",
            "source_backed",
            "If the energies are finite, Mellin-Plancherel writes J_alpha "
            "with denominator lambda^2+t^2 and I_alpha with denominator "
            "beta^2+t^2 against the same |zeta(beta+it)|^(-2).",
            "Hardy-space continuation from the initial half-plane; "
            "finiteness excludes a reciprocal-zeta pole on the boundary.",
        ),
        row(
            "mwpad_10_shared_reciprocal_zeta_spectrum",
            "exact_identity",
            "available_exact",
            "J_alpha=(1/(2pi))*int_R "
            "dt/[(lambda^2+t^2)|zeta(beta+it)|^2], while "
            "I_alpha has beta^2+t^2 in the first factor.",
            "The same reciprocal-zeta boundary data appear in both "
            "coordinates; only a scalar Hardy multiplier changes.",
        ),
        row(
            "mwpad_11_weighted_prefix_cell_identity",
            "exact_identity",
            "available_exact",
            "J_alpha=(1/(1-alpha))*sum_(k>=1)A_alpha(k)^2*"
            "[k^(alpha-1)-(k+1)^(alpha-1)].",
            "A_alpha(x)=A_alpha(k) on every half-open cell [k,k+1).",
        ),
        row(
            "mwpad_12_weighted_prefix_cell_comparison",
            "exact_inequality",
            "available_exact",
            "For E_alpha=sum k^(alpha-2)A_alpha(k)^2, "
            "2^(alpha-2)E_alpha<=J_alpha<=E_alpha.",
            "Compare x^(alpha-2) on each unit cell.",
        ),
        row(
            "mwpad_13_mertens_cell_comparison",
            "exact_inequality",
            "available_exact",
            "For M_alpha=sum M(k)^2/k^(2+alpha), "
            "2^(-2-alpha)M_alpha<=I_alpha<=M_alpha.",
            "The companion unit-cell comparison in the Mertens coordinate.",
        ),
        row(
            "mwpad_14_discrete_energy_equivalence",
            "exact_equivalence",
            "available_exact",
            "2^(-2-alpha)M_alpha<=E_alpha<="
            "2^(2-alpha)*((1+alpha)/(1-alpha))^2*M_alpha.",
            "Rows 7, 12, and 13; the constants are explicit and nonsharp "
            "only because of unit-cell comparison.",
        ),
        row(
            "mwpad_15_weighted_prefix_rh_criterion",
            "exact_equivalence",
            "available_exact",
            "RH iff E_(alpha_j)<infinity for every member of any fixed "
            "cofinal sequence alpha_j->0 in (0,1).",
            "Compose row 14 with the already checked cofinal "
            "OU/Mertens criterion; this is an equivalent target, not a proof.",
        ),
        row(
            "mwpad_16_burnol_coefficient_match",
            "exact_identity",
            "available_exact",
            "With alpha=2omega, A_alpha(k) is exactly the weighted Mobius "
            "prefix M_alpha(k) used in the Burnol tail-discrepancy reduction.",
            "No normalization or deleted weight is hidden in the comparison.",
        ),
        row(
            "mwpad_17_first_block_floor_collapse",
            "exact_identity",
            "available_exact",
            "For N<k<=2N, H_(omega,N)(k)=A_alpha(k)-A_alpha(N).",
            "Every omitted divisor satisfies N<d<=k and hence floor(k/d)=1.",
        ),
        row(
            "mwpad_18_first_block_affine_defect",
            "exact_identity",
            "available_exact",
            "For N<k<=2N, D_(omega,N)(k)="
            "A_alpha(k)-A_alpha(N)-k*r_(omega,N).",
            "Insert row 17 into the exact finite-tail discrepancy.",
        ),
        row(
            "mwpad_19_weighted_tail_abel",
            "exact_identity",
            "available_exact",
            "r_(omega,N)=-A_alpha(N)/(N+1)+"
            "sum_(n>N)A_alpha(n)/(n(n+1)).",
            "Discrete Abel summation; A_alpha(n)/n->0 follows "
            "unconditionally from alpha>0.",
        ),
        row(
            "mwpad_20_residual_weight_comparison",
            "exact_inequality",
            "available_exact",
            "If V_N=sum_(N<k<=2N)|D_N(k)|^2 and "
            "E_D=sum k^(alpha-2)|D_N(k)|^2, then "
            "2^(alpha-2)N^(alpha-2)V_N<=E_D<=N^(alpha-2)V_N.",
            "Dyadic comparison for the decreasing weight k^(alpha-2).",
        ),
        row(
            "mwpad_21_affine_gram",
            "exact_identity",
            "available_exact",
            "For W_q=sum_(N<k<=2N)k^(alpha-2+q), the weighted energy of "
            "A_alpha(N)+k*r_N is A_alpha(N)^2*W_0+"
            "2A_alpha(N)r_N*W_1+r_N^2*W_2.",
            "The Gram determinant W_0*W_2-W_1^2 is strictly positive.",
        ),
        row(
            "mwpad_22_three_term_block_bound",
            "exact_inequality",
            "available_exact",
            "F_alpha(N)^(1/2)<=N^((alpha-2)/2)V_N^(1/2)+"
            "N^((alpha-1)/2)|A_alpha(N)|+"
            "2^(alpha/2)N^((1+alpha)/2)|r_N|.",
            "Weighted triangle inequality plus W_0<=N^(alpha-1) and "
            "W_2<=2^alpha*N^(1+alpha).",
        ),
        row(
            "mwpad_23_anchor_gauge_invariance",
            "exact_identity",
            "available_exact",
            "Adding a constant C to A_alpha(N) and every A_alpha(k) on "
            "N<k<=2N leaves D_N(k) unchanged when r_N is fixed.",
            "The first q-block is a high-pass observation with one "
            "unobserved constant anchor mode.",
        ),
        row(
            "mwpad_24_constant_tail_countermodel",
            "proof_guard",
            "guard_validated",
            "For any finitely supported coefficient sequence whose prefix "
            "is the nonzero constant C after N, r_N=0 and D_N(k)=0 for "
            "all k>N, but F_alpha(N)=C^2*W_0>0.",
            "This is an operator countermodel, not a surrogate for mu; it "
            "proves that post-prefix residual and tail data alone cannot "
            "bound the absolute weighted-prefix block.",
        ),
        row(
            "mwpad_25_tail_rate_does_not_fix_anchor",
            "proof_guard",
            "guard_validated",
            "Even the forced rate r_N=O(N^(-(1+alpha)/2)) does not control "
            "N^(alpha-1)A_alpha(N)^2 without an additional anchor theorem.",
            "The constant-tail countermodel has r_N=0 and arbitrary anchor.",
        ),
        row(
            "mwpad_26_correlation_one_dimensional_collapse",
            "exact_identity",
            "available_exact",
            "sum_(h<X)sum_(Y<=X-h)C_h(Y)="
            "sum_(n=2)^X (X-n+1)mu(n)M(n-1).",
            "Reorder the triangle by the larger index n=m+h.",
        ),
        row(
            "mwpad_27_discrete_energy_difference",
            "exact_identity",
            "available_exact",
            "M(n)^2-M(n-1)^2=2mu(n)M(n-1)+mu(n)^2, and triangular "
            "summation of this identity recovers sum_(k<=X)M(k)^2.",
            "The one-dimensional collapse is exact but still contains the "
            "unknown Mertens prefix endogenously.",
        ),
        row(
            "mwpad_28_vaughan_source",
            "literature_guard",
            "source_backed",
            "Green and Tao, Lemma 4.1, give Vaughan's exact Type I/II "
            "decomposition for sums of mu(n)f(n) on a dyadic interval.",
            "Only the published algebraic identity is imported; its "
            "bounded-test estimates are not promoted to the endogenous test.",
        ),
        row(
            "mwpad_29_vaughan_type_i_ii_identity",
            "exact_identity",
            "available_exact",
            "For 1<=U,V<=N, sum_(N<n<=2N)mu(n)f(n)=-TI(f)+TII(f), "
            "with a_d=sum_(bc=d,b<=U,c<=V)mu(b)mu(c) and "
            "b_d=sum_(c|d,c>V)mu(c).",
            "Finite divisor reindexing; independently checked for arbitrary "
            "finite test vectors.",
        ),
        row(
            "mwpad_30_endogenous_test_guard",
            "proof_guard",
            "guard_validated",
            "After row 26 the Vaughan test is "
            "f_X(n)=(X-n+1)M(n-1), not an exogenous bounded function.",
            "The inverse theorem for a bounded external f cannot be quoted "
            "without retaining the Mertens dependence.",
        ),
        row(
            "mwpad_31_cauchy_schwarz_recycling_guard",
            "proof_guard",
            "guard_validated",
            "sum_(n<=X)|f_X(n)|^2<=X^2*sum_(k<=X)M(k)^2, so generic "
            "Cauchy-Schwarz in the Type I/II forms recycles the target "
            "energy on the right-hand side.",
            "No small absorbable coefficient is supplied by the identity "
            "alone; logarithmic Mobius orthogonality is one power short.",
        ),
        row(
            "mwpad_32_precollapse_bounded_shift_route",
            "exact_identity",
            "available_exact",
            "Before collapsing h and Y, Vaughan applies with the bounded "
            "test f_(h,Y)(m)=mu(m+h)1_(m<=Y) on each dyadic m-block, "
            "giving T(X)=B_1(X)-TI_X+TII_X exactly.",
            "The explicit aggregate avoids the endogenous-test defect, but "
            "no O_epsilon(X^(2+epsilon)) estimate for -TI_X+TII_X is "
            "asserted.",
        ),
        row(
            "mwpad_33_open_anchor_gate",
            "open_target",
            "open_target",
            "Control the constant anchor mode "
            "N^((alpha-1)/2)A_alpha(N) in a summable dyadic theorem, using "
            "the limiting Burnol energy or another arithmetic input.",
            "Post-prefix V_N and the reciprocal-tail rate do not supply it.",
        ),
        row(
            "mwpad_34_open_type_i_ii_gate",
            "open_target",
            "open_target",
            "Prove an O_epsilon(X^(2+epsilon)) signed aggregate estimate "
            "for the pre-collapse Vaughan Type I/II forms over h and Y, "
            "without rowwise absolute values or RH.",
            "This is the live noncircular bilinear handoff; no RH, "
            "PF-infinity, or Lambda<=0 conclusion is recorded.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction",
        "date": "2026-07-23",
        "status": (
            "exact weighted-prefix, first-block affine-defect, and Vaughan "
            "handoff reduction with two open arithmetic gates"
        ),
        "parameters": {
            "alpha_range": "0<alpha<1",
            "beta": "(1+alpha)/2",
            "lambda": "(1-alpha)/2",
            "burnol_match": "alpha=2omega",
        },
        "source_anchors": [
            "https://doi.org/10.5802/aif.2401",
            "https://arxiv.org/abs/math/0202166",
            "https://arxiv.org/abs/2508.00388",
        ],
        "rows": rows,
        "audit": {
            "row_count": 34,
            "exact_reduction_count": 24,
            "classical_theorem_step_count": 1,
            "literature_guard_count": 1,
            "proof_guard_count": 4,
            "open_handoff_gate_count": 2,
            "weighted_prefix_rh_equivalence_proved": True,
            "first_block_anchor_defect_proved": True,
            "vaughan_identity_closes_gate": False,
            "anchor_gate_proved": False,
            "signed_type_i_ii_gate_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Mertens Weighted-Prefix/Affine-Defect Reduction

Date: 2026-07-23

Status: exact weighted-prefix, first-block affine-defect, and Vaughan
handoff reduction with two open arithmetic gates. This is not a proof of
RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.json
python work/rh_compute/scripts/jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_weighted_prefix_affine_defect_reduction.py
```

## Weighted Prefix Is The Same Hardy Energy

Fix `0<alpha<1` and set

```text
M(x):=sum_(n<=x)mu(n),
A_alpha(x):=sum_(n<=x)mu(n)n^(-alpha),

I_alpha:=integral_1^infinity M(x)^2*x^(-2-alpha)dx,
J_alpha:=integral_1^infinity A_alpha(x)^2*x^(alpha-2)dx,

beta:=(1+alpha)/2,
lambda:=(1-alpha)/2.                                      (MWPADR.1)
```

Stieltjes integration by parts in both directions gives

```text
A_alpha(x)
 =x^(-alpha)M(x)
  +alpha*integral_1^x M(u)u^(-alpha-1)du,                  (MWPADR.2)

M(x)
 =x^alpha*A_alpha(x)
  -alpha*integral_1^x A_alpha(u)u^(alpha-1)du.             (MWPADR.3)
```

With `x=e^t`, define

```text
m(t):=e^(-beta*t)M(e^t),
a(t):=e^(-lambda*t)A_alpha(e^t).
```

Then `||m||_2^2=I_alpha`, `||a||_2^2=J_alpha`, and

```text
a=m+alpha*V_lambda*m,
m=a-alpha*V_beta*a,

(V_c f)(t):=integral_0^t exp(-c(t-u))f(u)du.                (MWPADR.4)
```

The forward Fourier multiplier is

```text
(beta+it)/(lambda+it).
```

Its modulus lies between `1` and `beta/lambda`. Therefore

```text
I_alpha
 <=J_alpha
 <=((1+alpha)/(1-alpha))^2*I_alpha.                         (MWPADR.5)
```

This is an exact norm comparison, not a generic estimate for either
arithmetic function.

The Mellin transform makes the shared spectrum explicit:

```text
integral_1^infinity A_alpha(x)x^(-s-1)dx
 =1/[s*zeta(s+alpha)]                                      (MWPADR.6)
```

initially for `Re(s)>1-alpha`. If the energies are finite,
Mellin-Plancherel gives

```text
J_alpha
 =(1/(2*pi))*integral_R
   dt/[(lambda^2+t^2)|zeta(beta+it)|^2],

I_alpha
 =(1/(2*pi))*integral_R
   dt/[(beta^2+t^2)|zeta(beta+it)|^2].                     (MWPADR.7)
```

Thus both coordinates contain the same reciprocal-zeta boundary data.

## Discrete Energy

Since both summatory functions are constant on unit cells,

```text
J_alpha
 =(1/(1-alpha))*sum_(k>=1)A_alpha(k)^2
   [k^(alpha-1)-(k+1)^(alpha-1)].                          (MWPADR.8)
```

Put

```text
E_alpha:=sum_(k>=1)k^(alpha-2)A_alpha(k)^2,
M_alpha:=sum_(k>=1)M(k)^2/k^(2+alpha).
```

Unit-cell comparison and (MWPADR.5) give

```text
2^(-2-alpha)M_alpha
 <=E_alpha
 <=2^(2-alpha)*((1+alpha)/(1-alpha))^2*M_alpha.             (MWPADR.9)
```

Consequently, for any fixed cofinal sequence `alpha_j->0`,

```text
RH
 iff
E_(alpha_j)<infinity for every j.                          (MWPADR.10)
```

This is a new exact coordinate for the existing RH-equivalent energy, not
a proof that it is finite.

## First q-Block And Its Missing Mode

Set `alpha=2omega`. Then `A_alpha` is exactly the weighted Mobius prefix in
the Burnol tail-discrepancy reduction. For `N<k<=2N`,
every omitted divisor has `floor(k/d)=1`, so

```text
H_(omega,N)(k)=A_alpha(k)-A_alpha(N),

D_(omega,N)(k)
 =A_alpha(k)-A_alpha(N)-k*r_(omega,N).                     (MWPADR.11)
```

Discrete Abel summation also gives

```text
r_(omega,N)
 =-A_alpha(N)/(N+1)
  +sum_(n>N)A_alpha(n)/(n(n+1)).                           (MWPADR.12)
```

Define

```text
V_N:=sum_(N<k<=2N)|D_(omega,N)(k)|^2,
E_D(N):=sum_(N<k<=2N)k^(alpha-2)|D_(omega,N)(k)|^2,
F_alpha(N):=sum_(N<k<=2N)k^(alpha-2)|A_alpha(k)|^2.
```

Then

```text
2^(alpha-2)N^(alpha-2)V_N
 <=E_D(N)
 <=N^(alpha-2)V_N.                                        (MWPADR.13)
```

For

```text
W_q(N):=sum_(N<k<=2N)k^(alpha-2+q),
```

the missing affine line has the exact Gram

```text
||A_alpha(N)+k*r_N||_w^2
 =A_alpha(N)^2*W_0
  +2A_alpha(N)r_N*W_1
  +r_N^2*W_2,                                              (MWPADR.14)
```

where `W_0*W_2-W_1^2>0`. In particular,

```text
F_alpha(N)^(1/2)
 <=N^((alpha-2)/2)V_N^(1/2)
   +N^((alpha-1)/2)|A_alpha(N)|
   +2^(alpha/2)N^((1+alpha)/2)|r_N|.                       (MWPADR.15)
```

The first q-block is invariant under adding the same constant to
`A_alpha(N)` and every later prefix in the block. This is a genuine
unobserved anchor mode, not a notational artifact.

For a generic finitely supported coefficient sequence, take its prefix to
be a nonzero constant `C` from `N` onward. Then

```text
r_N=0,
D_N(k)=0 for every k>N,
F_alpha(N)=C^2*W_0(N)>0.                                   (MWPADR.16)
```

The exact tail identity (MWPADR.12) still holds. Hence no operator
inequality can bound the absolute prefix energy from `V_N` and the forced
reciprocal-tail rate alone. A separate arithmetic theorem must control the
constant anchor `N^((alpha-1)/2)A_alpha(N)`.

## Anchored Correlation Collapses, But Does Not Decouple

The signed two-dimensional correlation from the preceding reduction has
the exact one-dimensional form

```text
sum_(h=1)^(X-1)sum_(Y=1)^(X-h)C_h(Y)
 =sum_(n=2)^X (X-n+1)mu(n)M(n-1).                          (MWPADR.17)
```

This follows by using the larger index `n=m+h`. It is the discrete energy
increment identity

```text
M(n)^2-M(n-1)^2
 =2mu(n)M(n-1)+mu(n)^2.                                    (MWPADR.18)
```

Thus the collapse removes one summation variable but leaves the unknown
Mertens prefix inside the test function.

## Vaughan Audit

Green and Tao's Mobius Vaughan identity gives, for `1<=U,V<=N`,

```text
sum_(N<n<=2N)mu(n)f(n)
 =-sum_(d<=U*V)a_d
    sum_(N/d<w<=2N/d)f(d*w)

  +sum_(V<d<=2N/U)b_d
    sum_(max(U,N/d)<w<=2N/d)mu(w)f(d*w),                   (MWPADR.19)

a_d:=sum_(b*c=d, b<=U, c<=V)mu(b)mu(c),
b_d:=sum_(c|d, c>V)mu(c).
```

The identity is finite and valid for an arbitrary test vector `f`.
Applied after (MWPADR.17), however, the test is

```text
f_X(n):=(X-n+1)M(n-1).
```

It is endogenous, and

```text
sum_(n<=X)|f_X(n)|^2
 <=X^2*sum_(k<=X)M(k)^2.                                   (MWPADR.20)
```

Therefore the standard Cauchy-Schwarz elimination of the Type I/II
coefficients feeds the target energy back into the estimate. Vaughan's
identity alone supplies no absorbable small factor and does not prove RH.

There is one noncircular handoff worth retaining. Before collapsing the
shift sums, apply (MWPADR.19) to each dyadic `m`-block with

```text
f_(h,Y)(m)=mu(m+h)1_(m<=Y),
```

which is externally bounded by one. Let `N` run over powers of two below
`X`, choose `1<=U_N,V_N<=N`, and define

```text
a_(N,d):=sum_(b*c=d,b<=U_N,c<=V_N)mu(b)mu(c),
b_(N,d):=sum_(c|d,c>V_N)mu(c),

TI_X
 :=sum_N sum_(h<X) sum_(Y<=X-h)
   sum_(d<=U_N*V_N) a_(N,d)
   sum_(N/d<w<=min(2N,Y)/d) mu(d*w+h),

TII_X
 :=sum_N sum_(h<X) sum_(Y<=X-h)
   sum_(V_N<d<=2N/U_N) b_(N,d)
   sum_(max(U_N,N/d)<w<=min(2N,Y)/d)
      mu(w)mu(d*w+h).                                      (MWPADR.21)
```

Empty sums vanish. The dyadic blocks `(N,2N]` cover `2<=m<=X`; the
omitted `m=1` boundary is

```text
B_1(X):=sum_(n=2)^X mu(n)(X-n+1)=O(X^2).
```

Termwise use of (MWPADR.19) gives the exact aggregate identity

```text
sum_(h<X)sum_(Y<=X-h)C_h(Y)
 =B_1(X)-TI_X+TII_X.                                      (MWPADR.22)
```

The remaining obligation is therefore the concrete signed estimate

```text
-TI_X+TII_X
 =O_epsilon(X^(2+epsilon)).                                (MWPADR.23)
```

The summation must retain cancellation jointly in the shift and terminal
variables. Rowwise absolute values restore the previously identified cubic
loss.

## Open Handoff

Two gates remain:

```text
1. Control N^((alpha-1)/2)A_alpha(N) in a summable dyadic
   theorem, using the limiting Burnol energy or another arithmetic input.

2. Prove (MWPADR.23) for the pre-collapse Vaughan Type I/II
   aggregate without assuming RH.
```

Neither gate is proved here. The reduction proves no full Burnol bound,
RH, PF-infinity, or `Lambda <= 0`.

## Source Boundary

- [Green and Tao, Quadratic uniformity of the Mobius
  function](https://doi.org/10.5802/aif.2401), Lemma 4.1, supplies the
  finite Vaughan Type I/II identity. Their bounded-test estimates are not
  attributed to the endogenous Mertens test.
- [Burnol](https://arxiv.org/abs/math/0202166) supplies the published
  weighted natural-approximant setting. The first-block and anchor-defect
  reductions above are derived here.
- [Das and Manna](https://arxiv.org/abs/2508.00388) provide modern
  Hardy/Copson context. The sharp scalar multiplier comparison
  (MWPADR.5) is derived directly here.
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
        "wrote Mertens weighted-prefix/affine-defect reduction: "
        f"{len(payload['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
