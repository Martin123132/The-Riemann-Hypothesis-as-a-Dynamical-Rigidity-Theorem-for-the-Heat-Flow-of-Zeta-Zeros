#!/usr/bin/env python3
"""Build the centered Mertens-bridge and Vaughan handoff reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_centered_bridge_vaughan_handoff.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_centered_bridge_vaughan_handoff.md"
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
            "mcbvh_01_tail_setup",
            "exact_definition",
            "available_exact",
            "For 0<alpha<1 set q_n=mu(n)n^(-1-alpha), "
            "r_N=sum_(n>N)q_n, and "
            "R_alpha=sum_(N>=1)N^alpha|r_N|^2.",
            "These are the checked reciprocal-Mobius tail coordinates.",
        ),
        row(
            "mcbvh_02_dyadic_dct",
            "exact_definition",
            "available_exact",
            "For dyadic K let y_j=r_(K+j-1), "
            "c_(K,r)=<y,phi_r>, and "
            "E_alpha=sum_(K dyadic)K^alpha"
            "sum_(0<=r<K)|c_(K,r)|^2.",
            "The Neumann DCT basis is that of Lemma 11.22V.",
        ),
        row(
            "mcbvh_03_common_kernels",
            "exact_definition",
            "available_exact",
            "Define h_(K,0)(n)=K^(-1/2)min((n-K)_+,K), and for "
            "1<=r<K define h_(K,r)(K+m)="
            "sqrt(2/K)sin(pi*r*m/K)/[2sin(pi*r/(2K))] "
            "for 1<=m<K, with zero values elsewhere.",
            "The mean kernel includes every future dyadic block; each "
            "nonconstant kernel is supported on (K,2K).",
        ),
        row(
            "mcbvh_04_kernel_identity",
            "exact_identity",
            "available_exact",
            "For every 0<=r<K, c_(K,r)="
            "sum_(n>=1)mu(n)n^(-1-alpha)h_(K,r)(n).",
            "Count the number of tails containing q_n for r=0 and use "
            "the checked sine formula for r>0.",
        ),
        row(
            "mcbvh_05_analysis_operator",
            "exact_definition",
            "available_exact",
            "Let A_alpha x have coordinates "
            "K^(alpha/2)sum_n x_n n^(-1-alpha)h_(K,r)(n) "
            "over dyadic K and 0<=r<K.",
            "For x_n=mu(n), ||A_alpha x||_2^2=E_alpha.",
        ),
        row(
            "mcbvh_06_cosine_energy_equivalence",
            "exact_equivalence",
            "available_exact",
            "R_alpha<infinity iff E_alpha<infinity.",
            "On K<=N<2K, the weights N^alpha and K^alpha differ by at "
            "most the fixed factor 2^alpha.",
        ),
        row(
            "mcbvh_07_mean_row_bound",
            "exact_inequality",
            "available_exact",
            "With F_(K,r)(n)=n^(-1-alpha)h_(K,r)(n), "
            "||F_(K,0)||_2^2<=[1+2^(-1-2alpha)/(1+2alpha)]K^(-2alpha).",
            "Split K<n<=2K from n>2K and compare the latter decreasing "
            "series with its integral.",
        ),
        row(
            "mcbvh_08_nonconstant_row_bound",
            "exact_inequality",
            "available_exact",
            "For 1<=r<K, ||F_(K,r)||_2^2<="
            "K^(-2alpha)/(4r^2).",
            "Use sum_(m=1)^(K-1)sin^2(pi*r*m/K)=K/2 and "
            "sin(pi*r/(2K))>=r/K.",
        ),
        row(
            "mcbvh_09_finite_trace",
            "exact_inequality",
            "available_exact",
            "sum_(K dyadic)K^alpha sum_(0<=r<K)"
            "||F_(K,r)||_2^2<infinity.",
            "Rows 7-8 give at most "
            "[1+2^(-1-2alpha)/(1+2alpha)+pi^2/24]"
            "sum_K K^(-alpha).",
        ),
        row(
            "mcbvh_10_signed_gram_expansion",
            "exact_identity",
            "available_exact",
            "For a finite dyadic cutoff J, E_alpha(J)=D_alpha(J)+"
            "2O_alpha(J), where D is the m=n part and "
            "O=sum_(m<n)mu(m)mu(n)G_(alpha,J)(m,n) for the Gram kernel "
            "G=sum_(K<=2^J,r)K^alpha F_(K,r)(m)F_(K,r)(n).",
            "Each finite-J row is l1, so the displayed expansion is "
            "absolutely justified before taking J to infinity.",
        ),
        row(
            "mcbvh_11_offdiagonal_criterion",
            "exact_equivalence",
            "available_exact",
            "E_alpha<infinity iff sup_J O_alpha(J)<infinity.",
            "The diagonal D_alpha(J) converges by row 9 and E_alpha(J) "
            "is a monotone nonnegative partial energy.",
        ),
        row(
            "mcbvh_12_trace_class_interpretation",
            "exact_consequence",
            "available_exact",
            "A_alpha is Hilbert-Schmidt from l2(N) to the dyadic DCT "
            "coefficient space, while the RH input mu is not in l2.",
            "The squared Hilbert-Schmidt norm is exactly the trace in "
            "row 9; this operator fact alone does not evaluate A_alpha mu.",
        ),
        row(
            "mcbvh_13_linf_operator_guard",
            "proof_guard",
            "guard_validated",
            "A_alpha is not bounded from l_infinity to l2.",
            "The zero-mean first-mode sequence of Lemma 11.22W gives "
            "|x_n|=|q*_n n^(1+alpha)|<=1 but divergent output energy. "
            "It is not the Mobius sequence.",
        ),
        row(
            "mcbvh_14_local_mertens_path",
            "exact_definition",
            "available_exact",
            "For 0<=m<K define S_(K,m)=M(K+m)-M(K), "
            "T_(K,m)=sum_(ell=1)^m mu(K+ell)(K+ell)^(-1-alpha).",
            "Both paths start at zero and T is the weighted path from "
            "Lemma 11.22W.",
        ),
        row(
            "mcbvh_15_quotient_norm",
            "exact_definition",
            "available_exact",
            "On K-vectors define ||[x]||_K=min_a||x-a1||_2; "
            "equivalently ||[x]||_K^2="
            "sum_m|x_m-K^(-1)sum_j x_j|^2.",
            "This is the Euclidean quotient norm modulo constants.",
        ),
        row(
            "mcbvh_16_extended_abel_operator",
            "exact_definition",
            "available_exact",
            "With a_m=(K+m)^(-1-alpha), define "
            "(A_Kx)_0=a_1x_0 and "
            "(A_Kx)_m=a_1x_0+sum_(j=1)^m a_j(x_j-x_(j-1)).",
            "For S_(K,0)=0, A_K S_K=T_K exactly.",
        ),
        row(
            "mcbvh_17_constant_preservation",
            "exact_identity",
            "available_exact",
            "A_K(c1)=a_1 c1, so A_K induces an invertible operator on "
            "the quotient by constants.",
            "All a_m are nonzero; the increment relation "
            "Delta(A_Kx)_m=a_m Delta x_m gives invertibility.",
        ),
        row(
            "mcbvh_18_abel_operator_upper",
            "exact_inequality",
            "available_exact",
            "||A_Kx||_2<=(2+alpha)K^(-1-alpha)||x||_2.",
            "Use the diagonal-plus-prefix Abel matrix, "
            "a_m<=K^(-1-alpha), "
            "a_m-a_(m+1)<=(1+alpha)K^(-2-alpha), and prefix norm <=K.",
        ),
        row(
            "mcbvh_19_abel_inverse",
            "exact_identity",
            "available_exact",
            "For b_m=1/a_m, (A_K^(-1)y)_0=b_1y_0 and "
            "(A_K^(-1)y)_m=b_m y_m+"
            "sum_(j=1)^(m-1)(b_j-b_(j+1))y_j.",
            "This is finite summation of "
            "Delta x_m=b_m Delta y_m.",
        ),
        row(
            "mcbvh_20_abel_inverse_upper",
            "exact_inequality",
            "available_exact",
            "||A_K^(-1)y||_2<="
            "2^alpha(3+alpha)K^(1+alpha)||y||_2.",
            "Use b_m<=(2K)^(1+alpha), "
            "b_(m+1)-b_m<=(1+alpha)(2K)^alpha, and prefix norm <=K.",
        ),
        row(
            "mcbvh_21_centered_abel_equivalence",
            "exact_inequality",
            "available_exact",
            "[2^alpha(3+alpha)]^(-1)K^(-1-alpha)||[S_K]||_K"
            "<=||[T_K]||_K<="
            "(2+alpha)K^(-1-alpha)||[S_K]||_K.",
            "Rows 17-20 descend to quotient norms; no DCT-commutator "
            "estimate is needed.",
        ),
        row(
            "mcbvh_22_centered_dct_identity",
            "exact_identity",
            "available_exact",
            "||[T_K]||_K^2=sum_(r=1)^(K-1)|tau_(K,r)|^2.",
            "This is the centered-path Parseval identity from "
            "Lemma 11.22W.",
        ),
        row(
            "mcbvh_23_high_modes_summable",
            "exact_inequality",
            "available_exact",
            "For R_K=ceil(sqrt(K)), "
            "sum_K K^alpha sum_(r>=R_K)|tau_(K,r)|^2<infinity.",
            "For r>0, tau=-c; apply the checked high-mode bound from "
            "Lemma 11.22V.",
        ),
        row(
            "mcbvh_24_centered_bridge_equivalence",
            "exact_equivalence",
            "available_exact",
            "sum_K K^alpha sum_(1<=r<R_K)|tau_(K,r)|^2<infinity "
            "iff L_alpha:=sum_K K^(-2-alpha)||[S_K]||_K^2<infinity.",
            "Add the summable high modes, then apply the two-sided "
            "quotient estimate in row 21.",
        ),
        row(
            "mcbvh_25_two_component_criterion",
            "exact_equivalence",
            "available_exact",
            "R_alpha<infinity iff M_alpha+L_alpha<infinity, where "
            "M_alpha=sum_K K^alpha|sqrt(K)r_K-tau_(K,0)|^2.",
            "Compose row 24 with the exact two-component criterion of "
            "Lemma 11.22W.",
        ),
        row(
            "mcbvh_26_cofinal_rh_criterion",
            "exact_equivalence",
            "available_exact",
            "For any fixed cofinal alpha_j->0, RH iff "
            "M_(alpha_j)<infinity and L_(alpha_j)<infinity for every j.",
            "This is a criterion only; neither component is estimated.",
        ),
        row(
            "mcbvh_27_pairwise_variance",
            "exact_identity",
            "available_exact",
            "||[S_K]||_K^2=K^(-1)"
            "sum_(0<=i<j<K)|S_(K,j)-S_(K,i)|^2.",
            "Use the standard finite variance/pairwise-distance identity.",
        ),
        row(
            "mcbvh_28_brownian_bridge_kernel",
            "exact_identity",
            "available_exact",
            "Writing x_i=mu(K+i), 1<=i<K, "
            "||[S_K]||_K^2=sum_(i,j)x_i x_j B_K(i,j), with "
            "B_K(i,j)=min(i,j)-ij/K="
            "min(i,j)(K-max(i,j))/K.",
            "The kernel is the discrete Brownian-bridge Green matrix "
            "L^*(I-K^(-1)11^T)L.",
        ),
        row(
            "mcbvh_29_bridge_diagonal_offdiagonal",
            "exact_identity",
            "available_exact",
            "||[S_K]||_K^2=Delta_K+2C_K, where "
            "Delta_K=sum_(i=1)^(K-1)mu(K+i)^2 i(K-i)/K and "
            "C_K=sum_(h=1)^(K-2)sum_(i=1)^(K-1-h)"
            "mu(K+i)mu(K+i+h)i(K-i-h)/K.",
            "Separate i=j and write j=i+h.",
        ),
        row(
            "mcbvh_30_bridge_diagonal_bound",
            "exact_inequality",
            "available_exact",
            "0<=Delta_K<=(K^2-1)/6 and "
            "sum_K K^(-2-alpha)Delta_K<infinity.",
            "Use |mu|<=1 and "
            "sum_(i=1)^(K-1)i(K-i)/K=(K^2-1)/6.",
        ),
        row(
            "mcbvh_31_bridge_offdiagonal_criterion",
            "exact_equivalence",
            "available_exact",
            "L_alpha<infinity iff "
            "sup_J sum_(K dyadic<=2^J)K^(-2-alpha)C_K<infinity.",
            "The diagonal series converges by row 30 and the partial "
            "bridge energies are nonnegative and monotone.",
        ),
        row(
            "mcbvh_32_vaughan_source",
            "literature_guard",
            "source_backed",
            "Green and Tao, Lemma 4.1, give the finite Vaughan Type I/II "
            "identity for sum_(K<n<=2K)mu(n)f(n).",
            "Only the exact divisor identity is imported; no bound for "
            "the present shifted-Mobius tests is attributed to the source.",
        ),
        row(
            "mcbvh_33_shifted_bridge_test",
            "exact_definition",
            "available_exact",
            "For 1<=h<=K-2 and i=n-K define "
            "g_(K,h)(n)=mu(n+h)i(K-i-h)/K when "
            "1<=i<=K-1-h, and zero otherwise; |g_(K,h)|<=K/4.",
            "Then C_K=sum_h sum_(K<n<=2K)mu(n)g_(K,h)(n).",
        ),
        row(
            "mcbvh_34_vaughan_coefficients",
            "exact_definition",
            "available_exact",
            "Choose 1<=U_K,V_K<=K and set "
            "a_(K,d)=sum_(bc=d,b<=U_K,c<=V_K)mu(b)mu(c), "
            "b_(K,d)=sum_(c|d,c>V_K)mu(c).",
            "These are the coefficients in the finite Vaughan identity.",
        ),
        row(
            "mcbvh_35_exact_vaughan_handoff",
            "exact_identity",
            "available_exact",
            "With TI_K=sum_h sum_(d<=U_KV_K)a_(K,d)"
            "sum_(K/d<w<=2K/d)g_(K,h)(dw), and the corresponding "
            "TII_K=sum_h sum_(V_K<d<=2K/U_K)b_(K,d)"
            "sum_(max(U_K,K/d)<w<=2K/d)mu(w)g_(K,h)(dw), "
            "one has C_K=-TI_K+TII_K.",
            "Apply Vaughan termwise in h; invalid endpoint values vanish "
            "by the definition of g.",
        ),
        row(
            "mcbvh_36_signed_type_i_ii_criterion",
            "exact_equivalence",
            "available_exact",
            "L_alpha<infinity iff "
            "sup_J sum_(K<=2^J)K^(-2-alpha)(-TI_K+TII_K)<infinity.",
            "Substitute row 35 into the exact off-diagonal criterion.",
        ),
        row(
            "mcbvh_37_separate_bound_guard",
            "proof_guard",
            "guard_validated",
            "The target in row 36 retains the signed Type I/II "
            "difference; replacing it by |TI_K|+|TII_K| is not an "
            "equivalent reduction.",
            "The algebraic cancellation in Vaughan is part of the "
            "original Mobius sum and no absorbable separate bound is proved.",
        ),
        row(
            "mcbvh_38_averaged_chowla_guard",
            "literature_guard",
            "source_backed",
            "Averaged Chowla gives the natural unweighted average "
            "o(K^2) after collapsing the redundant shift average, but "
            "the Brownian-bridge weight is O(K), leaving o(K^3), one "
            "full power above the K^(2+epsilon) scale required here.",
            "This is exponent bookkeeping for the published averaged "
            "theorem, not a lower bound on C_K or a novelty claim.",
        ),
        row(
            "mcbvh_39_mean_and_scope_guard",
            "proof_guard",
            "guard_validated",
            "The centered bridge removes constants and therefore does "
            "not control M_alpha; both components in row 25 remain "
            "necessary. The reduction concerns R_alpha, not the full "
            "Burnol energy Q.",
            "No centered-path estimate may be promoted to RH without the "
            "affine mean series or an independently sufficient full-Q theorem.",
        ),
        row(
            "mcbvh_40_open_two_component_gate",
            "open_target",
            "open_target",
            "For every alpha in one fixed cofinal sequence, prove both "
            "the affine mean summability M_alpha<infinity and the signed "
            "Brownian-bridge/Vaughan criterion in row 36.",
            "No Type I/II power gain, mean estimate, full Burnol bound, "
            "RH, PF-infinity, or Lambda<=0 is proved.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_mertens_centered_bridge_vaughan_handoff",
        "date": "2026-07-23",
        "status": (
            "exact trace/off-diagonal and centered Brownian-bridge "
            "Vaughan reduction with one open two-component gate"
        ),
        "parameters": {
            "alpha_range": "0<alpha<1",
            "dyadic_blocks": "K=2^j",
            "low_mode_cutoff": "R_K=ceil(sqrt(K))",
            "vaughan_cutoffs": "1<=U_K,V_K<=K",
        },
        "source_anchors": [
            "https://doi.org/10.5802/aif.2401",
            "https://arxiv.org/abs/1503.05121",
            "https://doi.org/10.2140/ant.2015.9.2167",
        ],
        "rows": rows,
        "audit": {
            "row_count": 40,
            "exact_reduction_count": 34,
            "literature_guard_count": 2,
            "proof_guard_count": 3,
            "open_two_component_gate_count": 1,
            "global_gram_diagonal_summable": True,
            "analysis_operator_hilbert_schmidt_on_l2": True,
            "analysis_operator_bounded_linf_to_l2": False,
            "centered_abel_quotient_equivalence_proved": True,
            "brownian_bridge_kernel_proved": True,
            "bridge_vaughan_handoff_proved": True,
            "averaged_chowla_closes_bridge": False,
            "signed_type_i_ii_gain_proved": False,
            "affine_mean_series_proved": False,
            "full_burnol_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Mertens Centered-Bridge Vaughan Handoff

Date: 2026-07-23

Status: exact trace/off-diagonal and centered Brownian-bridge Vaughan
reduction with one open two-component gate. This is not a proof of RH,
PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_centered_bridge_vaughan_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_centered_bridge_vaughan_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_centered_bridge_vaughan_handoff.py
```

## One Global Coefficient Kernel

Fix `0<alpha<1` and put

```text
q_n:=mu(n)n^(-1-alpha),
r_N:=sum_(n>N)q_n,
R_alpha:=sum_(N>=1)N^alpha|r_N|^2.                    (MCBVH.1)
```

For dyadic `K`, let

```text
y_j:=r_(K+j-1),                       1<=j<=K,
c_(K,r):=<y,phi_r>,                   0<=r<K.
```

Define coefficient kernels by

```text
h_(K,0)(n)
 :=K^(-1/2)min((n-K)_+,K),                            (MCBVH.2)

h_(K,r)(K+m)
 :=sqrt(2/K)
   sin(pi*r*m/K)/[2sin(pi*r/(2K))],
                                      1<=m<K, 1<=r<K, (MCBVH.3)
```

and set the nonconstant kernels to zero elsewhere. Counting how many
tails contain each `q_n`, together with Lemma 11.22V, gives

```text
c_(K,r)
 =sum_(n>=1)mu(n)n^(-1-alpha)h_(K,r)(n),
                                      0<=r<K.          (MCBVH.4)
```

Thus the dyadic cosine energy

```text
E_alpha
 :=sum_(K dyadic)K^alpha
   sum_(r=0)^(K-1)|c_(K,r)|^2                         (MCBVH.5)
```

is the squared norm of one analysis operator `A_alpha` applied to the
sequence `mu(n)`. The dyadic weight comparison from Lemma 11.22V gives

```text
R_alpha<infinity iff E_alpha<infinity.                (MCBVH.6)
```

Let

```text
F_(K,r)(n):=n^(-1-alpha)h_(K,r)(n).
```

The mean row splits at `2K`:

```text
||F_(K,0)||_2^2
 <=[1+2^(-1-2alpha)/(1+2alpha)]K^(-2alpha).           (MCBVH.7)
```

For `r>0`, finite sine orthogonality gives

```text
sum_(m=1)^(K-1)|h_(K,r)(K+m)|^2
 =1/[4sin^2(pi*r/(2K))].
```

Since `sin(pi*r/(2K))>=r/K`,

```text
||F_(K,r)||_2^2
 <=K^(-2alpha)/(4r^2).                                (MCBVH.8)
```

Consequently

```text
sum_(K dyadic)K^alpha
 sum_(r=0)^(K-1)||F_(K,r)||_2^2

 <=[1+2^(-1-2alpha)/(1+2alpha)+pi^2/24]
   sum_(K dyadic)K^(-alpha)
 <infinity.                                           (MCBVH.9)
```

This proves that `A_alpha` is Hilbert-Schmidt from `l2(N)` to the
dyadic coefficient space. It does not close the problem because
`mu` is not in `l2`.

For a finite dyadic cutoff `J`, define

```text
G_(alpha,J)(m,n)
 :=sum_(K dyadic<=2^J)K^alpha
   sum_(r=0)^(K-1)F_(K,r)(m)F_(K,r)(n).
```

Every row entering a fixed `J` is in `l1`, so expansion before passage
to the limit is legitimate:

```text
E_alpha(J)
 =D_alpha(J)+2O_alpha(J),                              (MCBVH.10)

D_alpha(J)
 :=sum_n mu(n)^2 G_(alpha,J)(n,n),

O_alpha(J)
 :=sum_(m<n)mu(m)mu(n)G_(alpha,J)(m,n).
```

The diagonal converges by (MCBVH.9), while `E_alpha(J)` is monotone and
nonnegative. Hence

```text
E_alpha<infinity
 iff
sup_J O_alpha(J)<infinity.                             (MCBVH.11)
```

The RH-strength content is therefore purely signed and off-diagonal.
There is no generic `l_infinity -> l2` operator theorem available:
the scalar construction in Lemma 11.22W has

```text
x*_n:=q*_n n^(1+alpha),      |x*_n|<=1,
```

but its first-mode output energy diverges. That sequence is not
Mobius; it is a proof-method guard.

## Centered Abel Quotient

For `0<=m<K`, define the ordinary and weighted local paths

```text
S_(K,m):=M(K+m)-M(K),

T_(K,m):=sum_(ell=1)^m
          mu(K+ell)(K+ell)^(-1-alpha).                 (MCBVH.12)
```

Both start at zero. On K-vectors use the quotient norm modulo constants

```text
||[x]||_K
 :=min_a||x-a1||_2,

||[x]||_K^2
 =sum_(m=0)^(K-1)
   |x_m-K^(-1)sum_j x_j|^2.                            (MCBVH.13)
```

With `a_m=(K+m)^(-1-alpha)`, extend the Abel map to every K-vector:

```text
(A_Kx)_0:=a_1x_0,

(A_Kx)_m
 :=a_1x_0+sum_(j=1)^m a_j(x_j-x_(j-1)),
                                      1<=m<K.          (MCBVH.14)
```

For the anchored Mertens path, `A_K S_K=T_K`. More importantly,

```text
A_K(c1)=a_1c1.                                        (MCBVH.15)
```

So `A_K` descends to an invertible map on the quotient by constants.
The diagonal-plus-prefix formula gives

```text
||A_Kx||_2
 <=(2+alpha)K^(-1-alpha)||x||_2.                       (MCBVH.16)
```

Writing `b_m=1/a_m`, the inverse is

```text
(A_K^(-1)y)_0=b_1y_0,

(A_K^(-1)y)_m
 =b_my_m+sum_(j=1)^(m-1)(b_j-b_(j+1))y_j,             (MCBVH.17)
```

and

```text
||A_K^(-1)y||_2
 <=2^alpha(3+alpha)K^(1+alpha)||y||_2.                 (MCBVH.18)
```

Because constants map to constants, the same bounds descend to the
quotient:

```text
[2^alpha(3+alpha)]^(-1)
 K^(-1-alpha)||[S_K]||_K

 <=||[T_K]||_K

 <=(2+alpha)K^(-1-alpha)||[S_K]||_K.                  (MCBVH.19)
```

This resolves the commutator concern from Lemma 11.22W.
No commutator estimate is needed: the Abel map itself is uniformly invertible after
passing to the correct quotient.

The centered DCT identity is

```text
||[T_K]||_K^2
 =sum_(r=1)^(K-1)|tau_(K,r)|^2.                        (MCBVH.20)
```

All modes `r>=ceil(sqrt(K))` are already dyadically summable. Therefore
(MCBVH.19)-(MCBVH.20) give

```text
sum_(K dyadic)K^alpha
 sum_(1<=r<ceil(sqrt(K)))|tau_(K,r)|^2
 <infinity

iff

L_alpha
 :=sum_(K dyadic)K^(-2-alpha)||[S_K]||_K^2
 <infinity.                                            (MCBVH.21)
```

Retain the affine mean series

```text
M_alpha
 :=sum_(K dyadic)K^alpha
   |sqrt(K)r_K-tau_(K,0)|^2.                           (MCBVH.22)
```

The exact two-component criterion is now

```text
R_alpha<infinity
 iff
M_alpha+L_alpha<infinity.                              (MCBVH.23)
```

For every member of any fixed cofinal sequence `alpha_j->0`,
(MCBVH.23) is equivalent to RH. Neither term has been proved finite.

## Discrete Brownian-Bridge Kernel

The finite pairwise-variance identity gives

```text
||[S_K]||_K^2
 =K^(-1)sum_(0<=i<j<K)
   |S_(K,j)-S_(K,i)|^2.                                (MCBVH.24)
```

Writing `x_i:=mu(K+i)`, `1<=i<K`, this is exactly

```text
||[S_K]||_K^2
 =sum_(i,j=1)^(K-1)x_i x_j B_K(i,j),                   (MCBVH.25)

B_K(i,j)
 :=min(i,j)-ij/K
  =min(i,j)(K-max(i,j))/K.
```

`B_K` is the discrete Brownian-bridge Green matrix. Separating its
diagonal and writing `j=i+h` gives

```text
||[S_K]||_K^2=Delta_K+2C_K,                            (MCBVH.26)

Delta_K
 :=sum_(i=1)^(K-1)
   mu(K+i)^2 i(K-i)/K,

C_K
 :=sum_(h=1)^(K-2)
   sum_(i=1)^(K-1-h)
   mu(K+i)mu(K+i+h)i(K-i-h)/K.                         (MCBVH.27)
```

The diagonal is harmless:

```text
0<=Delta_K<=(K^2-1)/6,

sum_(K dyadic)K^(-2-alpha)Delta_K<infinity.            (MCBVH.28)
```

As the bridge energies are nonnegative, the complete nonconstant gate
is therefore

```text
L_alpha<infinity
 iff
sup_J sum_(K dyadic<=2^J)K^(-2-alpha)C_K
 <infinity.                                            (MCBVH.29)
```

This is an anchored but base-averaged, triangularly weighted two-point
Mobius correlation.

## Pre-Collapse Vaughan Handoff

For `1<=h<=K-2`, put `i=n-K` and define

```text
g_(K,h)(n)
 :=mu(n+h)i(K-i-h)/K
   1_(1<=i<=K-1-h).                                    (MCBVH.30)
```

Then

```text
|g_(K,h)(n)|<=K/4,

C_K
 =sum_(h=1)^(K-2)
   sum_(K<n<=2K)mu(n)g_(K,h)(n).                       (MCBVH.31)
```

Choose `1<=U_K,V_K<=K` and define the Vaughan coefficients

```text
a_(K,d)
 :=sum_(bc=d,b<=U_K,c<=V_K)mu(b)mu(c),

b_(K,d)
 :=sum_(c|d,c>V_K)mu(c).                               (MCBVH.32)
```

Green and Tao's finite Vaughan identity gives

```text
TI_K
 :=sum_h sum_(d<=U_KV_K)a_(K,d)
   sum_(K/d<w<=2K/d)g_(K,h)(d*w),

TII_K
 :=sum_h sum_(V_K<d<=2K/U_K)b_(K,d)
   sum_(max(U_K,K/d)<w<=2K/d)
   mu(w)g_(K,h)(d*w),

C_K=-TI_K+TII_K.                                      (MCBVH.33)
```

The exact live Type I/II target is

```text
sup_J sum_(K dyadic<=2^J)
 K^(-2-alpha)(-TI_K+TII_K)
 <infinity.                                            (MCBVH.34)
```

The signed difference must be retained. Replacing it with
`|TI_K|+|TII_K|` is not an equivalent handoff.

The elementary absolute scale is still one full power too large:

```text
|C_K|
 <=sum_(h,i)i(K-i-h)/K
 <=K^3/8.                                             (MCBVH.35)
```

After multiplication by `K^(-2-alpha)`, this gives the divergent scale
`K^(1-alpha)`. The averaged Chowla theorem of Matomaki, Radziwill, and
Tao gives the natural unweighted average `o(K^2)` after collapsing the
redundant common-shift average. Multiplication by the Brownian-bridge
weight `O(K)` still gives only `o(K^3)`, not the
`O_epsilon(K^(2+epsilon))` scale required here. This is a power audit,
not a lower bound on the actual correlation.

## Open Gate

The quotient reduction removes only local constants. It does not control
the affine mean series `M_alpha`. Thus the surviving cofinal gate has
two noninterchangeable pieces:

```text
1. Prove M_alpha<infinity.

2. Prove the signed Brownian-bridge/Vaughan criterion (MCBVH.34).
```

A bound

```text
||[S_K]||_K^2=O_epsilon(K^(2+epsilon))
```

with `epsilon<alpha` is a sufficient calibration for the second piece,
but is not proved here. The exact reduction concerns `R_alpha`, not the
full Burnol energy `Q`. No Type I/II power gain, affine mean estimate,
full-Q bound, RH, PF-infinity, or `Lambda <= 0` is proved.

Primary sources:

- Green and Tao, *Quadratic Uniformity of the Mobius Function*,
  Lemma 4.1 for the finite Vaughan identity:
  https://doi.org/10.5802/aif.2401
- Matomaki, Radziwill, and Tao, *An Averaged Form of Chowla's
  Conjecture*:
  https://doi.org/10.2140/ant.2015.9.2167
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
        "built Mertens centered-bridge Vaughan handoff: "
        f"{len(payload['rows'])} rows, "
        f"{payload['audit']['exact_reduction_count']} exact reductions, "
        f"{payload['audit']['proof_guard_count']} proof guards, "
        f"{payload['audit']['open_two_component_gate_count']} open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
