#!/usr/bin/env python3
"""Build the local Mertens-path cosine transference reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_local_path_cosine_transference.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_local_path_cosine_transference.md"
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
            "mlpath_01_tail_definitions",
            "exact_definition",
            "available_exact",
            "For 0<alpha<1 set q_n=mu(n)n^(-1-alpha), "
            "r_N=sum_(n>N)q_n, and "
            "R_alpha=sum_(N>=1)N^alpha|r_N|^2.",
            "These are the checked reciprocal-Mobius tail coordinates.",
        ),
        row(
            "mlpath_02_dyadic_tail_vector",
            "exact_definition",
            "available_exact",
            "For dyadic K let y_j=r_(K+j-1), 1<=j<=K, and let "
            "c_(K,r)=<y,phi_r> in the Neumann DCT basis.",
            "The basis phi_r and its high-mode estimate are those of "
            "Lemma 11.22V.",
        ),
        row(
            "mlpath_03_weighted_local_path",
            "exact_definition",
            "available_exact",
            "Define T_(K,m)=sum_(ell=1)^m q_(K+ell)="
            "r_K-r_(K+m), 0<=m<=K, with T_(K,0)=0, and set "
            "z_j=T_(K,j-1), 1<=j<=K.",
            "T is the origin-anchored reciprocal-weighted Mobius path "
            "inside the dyadic annulus.",
        ),
        row(
            "mlpath_04_affine_tail_path",
            "exact_identity",
            "available_exact",
            "The tail vector is exactly y_j=r_K-z_j.",
            "This keeps the incoming tail level and the local path in "
            "one affine identity; neither term may be discarded.",
        ),
        row(
            "mlpath_05_path_dct",
            "exact_definition",
            "available_exact",
            "Let tau_(K,r)=<z,phi_r>, 0<=r<K.",
            "This is an orthonormal DCT of the K-vector "
            "(T_(K,0),...,T_(K,K-1)).",
        ),
        row(
            "mlpath_06_affine_dct_transference",
            "exact_identity",
            "available_exact",
            "c_(K,0)=sqrt(K)r_K-tau_(K,0), while "
            "c_(K,r)=-tau_(K,r) for 1<=r<K.",
            "The nonconstant DCT vectors have zero sum; the constant "
            "vector has inner product sqrt(K) with the all-ones vector.",
        ),
        row(
            "mlpath_07_summation_by_parts",
            "exact_identity",
            "available_exact",
            "For 1<=r<K, c_(K,r)=-sqrt(2/K)"
            "sum_(m=1)^(K-1)T_(K,m)"
            "cos(pi*r*(m+1/2)/K).",
            "Discrete summation by parts in the sine formula of "
            "Lemma 11.22V cancels its 1/sin(pi*r/(2K)) factor exactly.",
        ),
        row(
            "mlpath_08_block_parseval",
            "exact_identity",
            "available_exact",
            "||y||_2^2=|sqrt(K)r_K-tau_(K,0)|^2+"
            "sum_(r=1)^(K-1)|tau_(K,r)|^2.",
            "Combine row 6 with DCT Parseval.",
        ),
        row(
            "mlpath_09_centered_path_variance",
            "exact_identity",
            "available_exact",
            "sum_(r=1)^(K-1)|tau_(K,r)|^2="
            "sum_(m=0)^(K-1)|T_(K,m)|^2-"
            "K^(-1)|sum_(m=0)^(K-1)T_(K,m)|^2.",
            "The nonconstant modes are exactly the centered local-path "
            "energy, not merely bounded by it.",
        ),
        row(
            "mlpath_10_deleted_dct_gram",
            "exact_identity",
            "available_exact",
            "After deleting the constant DCT row and the fixed column "
            "T_(K,0)=0, the (K-1)-square transform has Gram matrix "
            "I_(K-1)-K^(-1)11^T.",
            "Restrict the full identity "
            "sum_(r>=1)phi_r phi_r^T=I-K^(-1)11^T.",
        ),
        row(
            "mlpath_11_deleted_dct_spectrum",
            "exact_consequence",
            "available_exact",
            "The deleted transform has singular values 1 with "
            "multiplicity K-2 and K^(-1/2) with multiplicity 1.",
            "The all-ones direction has Gram eigenvalue 1/K; its "
            "orthogonal complement has eigenvalue 1.",
        ),
        row(
            "mlpath_12_low_mode_path_criterion",
            "exact_equivalence",
            "available_exact",
            "With R_K=ceil(sqrt(K)), R_alpha<infinity iff "
            "sum_(K dyadic)K^alpha("
            "|sqrt(K)r_K-tau_(K,0)|^2+"
            "sum_(1<=r<R_K)|tau_(K,r)|^2)<infinity.",
            "Rows 6 and 8 rewrite the exact low-mode criterion of "
            "Lemma 11.22V; its high modes remain automatically summable.",
        ),
        row(
            "mlpath_13_cofinal_rh_criterion",
            "exact_equivalence",
            "available_exact",
            "For any fixed cofinal alpha_j->0 in (0,1), RH is "
            "equivalent to the path criterion in row 12 for every "
            "alpha_j.",
            "This composes row 12 with the checked cofinal R_alpha "
            "criterion and does not prove the displayed series finite.",
        ),
        row(
            "mlpath_14_two_component_split",
            "exact_consequence",
            "available_exact",
            "The surviving target has two components: the affine mean "
            "defect sqrt(K)r_K-tau_(K,0), and the first O(sqrt(K)) "
            "nonconstant DCT coefficients of the local path.",
            "Control of only one component does not algebraically "
            "control the other.",
        ),
        row(
            "mlpath_15_local_mertens_path",
            "exact_definition",
            "available_exact",
            "Let S_(K,m)=M(K+m)-M(K)="
            "sum_(ell=1)^m mu(K+ell), and "
            "a_m=(K+m)^(-1-alpha).",
            "S is the ordinary local Mertens-increment path on the same "
            "anchored annulus.",
        ),
        row(
            "mlpath_16_abel_forward",
            "exact_identity",
            "available_exact",
            "For 1<=m<=K, T_(K,m)=a_m S_(K,m)+"
            "sum_(ell=1)^(m-1)S_(K,ell)(a_ell-a_(ell+1)).",
            "This is finite Abel summation applied to "
            "sum a_ell(S_ell-S_(ell-1)).",
        ),
        row(
            "mlpath_17_abel_inverse",
            "exact_identity",
            "available_exact",
            "With b_m=1/a_m, S_(K,m)=b_m T_(K,m)+"
            "sum_(ell=1)^(m-1)T_(K,ell)(b_ell-b_(ell+1)).",
            "Sum (T_ell-T_(ell-1))/a_ell. The correction is negative "
            "because b_m is increasing.",
        ),
        row(
            "mlpath_18_path_norm_upper",
            "exact_inequality",
            "available_exact",
            "For the vectors 1<=m<K, "
            "||T_K||_2<=(2+alpha)K^(-1-alpha)||S_K||_2.",
            "Use a_m<=K^(-1-alpha), "
            "a_ell-a_(ell+1)<=(1+alpha)K^(-2-alpha), "
            "and the strict-prefix operator norm <=K.",
        ),
        row(
            "mlpath_19_path_norm_lower",
            "exact_inequality",
            "available_exact",
            "For the same vectors, "
            "||T_K||_2>=[2^alpha(3+alpha)]^(-1)"
            "K^(-1-alpha)||S_K||_2.",
            "Apply row 17 with b_m<=(2K)^(1+alpha), "
            "b_(ell+1)-b_ell<=(1+alpha)(2K)^alpha, and prefix norm <=K.",
        ),
        row(
            "mlpath_20_projection_scope_guard",
            "proof_guard",
            "guard_validated",
            "Rows 18-19 compare full local paths only; the triangular "
            "Abel transform does not commute with projection onto the "
            "first O(sqrt(K)) DCT modes.",
            "A full-path norm equivalence cannot be promoted to a "
            "low-mode estimate without a separate commutator or signed "
            "arithmetic argument.",
        ),
        row(
            "mlpath_21_generic_countermodel",
            "exact_definition",
            "available_exact",
            "For fixed alpha and every dyadic K>=4 set "
            "A_K=2^(-2-alpha)pi^(-1)K^(-alpha), "
            "T*_(K,m)=A_K sin(2pi*m/K), and "
            "q*_(K+m)=T*_(K,m)-T*_(K,m-1).",
            "Set q*_1=0. The disjoint annuli (K,2K] define one global "
            "scalar sequence; it is not the Mobius sequence.",
        ),
        row(
            "mlpath_22_countermodel_envelope",
            "exact_inequality",
            "available_exact",
            "|q*_n|<=n^(-1-alpha), and the sequence is absolutely "
            "summable because each dyadic block has l1 norm <=2pi*A_K.",
            "Use |sin u-sin v|<=|u-v| and n<=2K; "
            "sum_(K dyadic)K^(-alpha)<infinity.",
        ),
        row(
            "mlpath_23_countermodel_zero_means",
            "exact_identity",
            "available_exact",
            "Every block has sum_(K<n<=2K)q*_n=0; hence r*_K=0 at "
            "every dyadic endpoint, tau*_(K,0)=0, and c*_(K,0)=0.",
            "T*_(K,K)=T*_(K,0)=0 and the complete discrete sine sum "
            "over m=0,...,K-1 vanishes.",
        ),
        row(
            "mlpath_24_countermodel_first_mode",
            "exact_identity",
            "available_exact",
            "|c*_(K,1)|=A_K sqrt(2/K) cos(pi/K)/2 times "
            "[csc(pi/(2K))+csc(3pi/(2K))].",
            "Use product-to-sum and the finite sine-sum formula on "
            "sin(2pi*m/K)cos(pi*(m+1/2)/K).",
        ),
        row(
            "mlpath_25_countermodel_divergence",
            "exact_consequence",
            "available_exact",
            "For K>=4, |c*_(K,1)|>=A_K sqrt(K)/pi and therefore "
            "K^alpha|c*_(K,1)|^2>="
            "2^(-4-2alpha)pi^(-4)K^(1-alpha); the dyadic series "
            "diverges although row 23 holds.",
            "Use cos(pi/K)>=sqrt(2)/2 and sin x<=x. Since alpha<1, "
            "the lower-bound terms do not even tend to zero.",
        ),
        row(
            "mlpath_26_generic_method_barrier",
            "proof_guard",
            "guard_validated",
            "Coefficient size, absolute convergence, zero dyadic "
            "endpoint tails, and zero block mean cannot control the "
            "first nonconstant mode by a generic Hilbert-space or "
            "large-sieve argument.",
            "Rows 21-25 are a generic scalar countermodel only. They "
            "show that Mobius-specific signed cancellation is necessary "
            "but do not describe the actual Mobius coefficients.",
        ),
        row(
            "mlpath_27_linear_phase_sources",
            "literature_guard",
            "source_backed",
            "Davenport gives the uniform arbitrary-log estimate for "
            "sum_(n<=X)mu(n)e(n theta); Zhang records this and quotes "
            "power refinements under a zero-free half-plane hypothesis "
            "for every Dirichlet L-function.",
            "These sources support the logarithmic unconditional bound "
            "and the stated conditional calibrations, not the open path "
            "series.",
        ),
        row(
            "mlpath_28_three_quarter_guard",
            "proof_guard",
            "guard_validated",
            "At the quoted all-Dirichlet zero-free endpoint a=1/2, the "
            "Baker-Harman/Zhang uniform exponent is 3/4+epsilon; the "
            "threshold from Lemma 11.22V then closes only "
            "alpha>1/2+2epsilon, not a cofinal alpha->0 sequence.",
            "This is exponent bookkeeping for that theorem, not a claim "
            "that stronger conditional or future methods are impossible.",
        ),
        row(
            "mlpath_29_rational_zero_density_guard",
            "literature_guard",
            "source_backed",
            "The unconditional Maier-Sankaranarayanan rational estimate "
            "quoted by Zhang retains Dirichlet-L zero-density terms; the "
            "focused source audit found no supported simplification to "
            "the uniform cofinal-closing power needed here.",
            "This is a focused search report, not a novelty, priority, "
            "or exhaustive-literature claim.",
        ),
        row(
            "mlpath_30_vaughan_handoff_guard",
            "proof_guard",
            "guard_validated",
            "A viable Vaughan or bilinear attack must act on mu before "
            "rowwise absolute values and preserve signed coupling across "
            "the low modes and scales; bounding the full local S path "
            "after row 19 merely recycles the Mertens-energy target.",
            "The exact transference identifies where a pre-collapse "
            "arithmetic gain must enter but supplies no such gain.",
        ),
        row(
            "mlpath_31_scope_guard",
            "proof_guard",
            "guard_validated",
            "Translation averages, almost-all short-interval theorems, "
            "and estimates for only the nonzero rational phases do not "
            "automatically control both components in row 14; this "
            "reduction also concerns R_alpha rather than the full Burnol "
            "energy Q.",
            "Each promotion requires an explicit de-averaging, mean-mode, "
            "or full-Q implication already absent from the proof.",
        ),
        row(
            "mlpath_32_open_two_component_gate",
            "open_target",
            "open_target",
            "For every alpha in one fixed cofinal sequence, prove "
            "dyadic summability of the affine mean defects and the "
            "first ceil(sqrt(K))-1 nonconstant local-path DCT modes in "
            "row 12, using Mobius-specific signed cancellation.",
            "This remains RH-equivalent. No low-mode power gain, RH, "
            "PF-infinity, or Lambda<=0 is proved.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_mertens_local_path_cosine_transference",
        "date": "2026-07-23",
        "status": (
            "exact local Mertens-path DCT transference with an isolated "
            "generic first-mode barrier and one open two-component gate"
        ),
        "parameters": {
            "alpha_range": "0<alpha<1",
            "dyadic_blocks": "K=2^j",
            "low_mode_cutoff": "R_K=ceil(sqrt(K))",
            "path_indices": "0<=m<K",
        },
        "source_anchors": [
            "https://doi.org/10.1093/qmath/os-8.1.313",
            "https://doi.org/10.5802/aif.2401",
            "https://doi.org/10.1112/jlms/s2-43.2.193",
            "https://doi.org/10.46298/hrj.2005.151",
            "https://arxiv.org/abs/2204.04613",
        ],
        "rows": rows,
        "audit": {
            "row_count": 32,
            "exact_reduction_count": 24,
            "literature_guard_count": 2,
            "proof_guard_count": 5,
            "open_two_component_gate_count": 1,
            "affine_dct_transference_proved": True,
            "centered_path_variance_proved": True,
            "local_mertens_abel_pair_proved": True,
            "full_path_norm_equivalence_proved": True,
            "low_projection_norm_equivalence_proved": False,
            "generic_zero_mean_first_mode_barrier_proved": True,
            "mobius_low_mode_gain_proved": False,
            "full_burnol_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Mertens Local-Path Cosine Transference

Date: 2026-07-23

Status: exact local Mertens-path DCT transference with an isolated
generic first-mode barrier and one open two-component gate. This is not a proof
of RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_local_path_cosine_transference.json
python work/rh_compute/scripts/jensen_window_pf_mertens_local_path_cosine_transference.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_local_path_cosine_transference.py
```

## Exact Local-Path Transform

Fix `0<alpha<1` and retain

```text
q_n:=mu(n)n^(-1-alpha),
r_N:=sum_(n>N)q_n,
R_alpha:=sum_(N>=1)N^alpha|r_N|^2.                    (MLPATH.1)
```

For dyadic `K`, set

```text
y_j:=r_(K+j-1),                                      1<=j<=K,

T_(K,m):=sum_(ell=1)^m q_(K+ell)
        =r_K-r_(K+m),                                0<=m<=K,

z_j:=T_(K,j-1).                                      (MLPATH.2)
```

Thus `T_(K,0)=0` and, exactly,

```text
y_j=r_K-z_j.                                         (MLPATH.3)
```

Let `phi_r` be the orthonormal Neumann DCT basis from Lemma 11.22V:

```text
phi_0(j)=K^(-1/2),

phi_r(j)=sqrt(2/K)cos(pi*r*(j-1/2)/K),      1<=r<K.
```

Write

```text
c_(K,r):=<y,phi_r>,
tau_(K,r):=<z,phi_r>.
```

Because every nonconstant DCT vector has zero sum, (MLPATH.3) gives the
affine transference

```text
c_(K,0)=sqrt(K)r_K-tau_(K,0),
c_(K,r)=-tau_(K,r),                         1<=r<K.  (MLPATH.4)
```

Discrete summation by parts gives the equivalent nonconstant formula

```text
c_(K,r)
 =-sqrt(2/K) sum_(m=1)^(K-1)
    T_(K,m)cos(pi*r*(m+1/2)/K),             1<=r<K.  (MLPATH.5)
```

The factor `1/sin(pi*r/(2K))` in the sine-weight formula from
Lemma 11.22V has canceled exactly. The low modes are therefore the
ordinary cosine coefficients of the reciprocal-weighted local Mertens
path, not generic high-denominator oscillatory sums.

Parseval now reads

```text
||y||_2^2
 =|sqrt(K)r_K-tau_(K,0)|^2
  +sum_(r=1)^(K-1)|tau_(K,r)|^2,                     (MLPATH.6)

sum_(r=1)^(K-1)|tau_(K,r)|^2
 =sum_(m=0)^(K-1)|T_(K,m)|^2
  -K^(-1)|sum_(m=0)^(K-1)T_(K,m)|^2.                 (MLPATH.7)
```

So the nonconstant spectrum is exactly the centered path variance. If
the constant DCT row and the fixed column `T_(K,0)=0` are deleted, the
remaining square transform `U_K` satisfies

```text
U_K^*U_K=I_(K-1)-K^(-1)11^T.                         (MLPATH.8)
```

Its singular values are `1` with multiplicity `K-2` and `K^(-1/2)`
with multiplicity one. The weak direction is the near-constant local
path. It cannot be separated from the affine mean defect in (MLPATH.4)
without losing the exact cancellation.

## RH-Equivalent Two-Component Gate

Let `R_K=ceil(sqrt(K))`. Lemma 11.22V already proves that every mode
`r>=R_K` is dyadically summable from coefficient size alone. Equations
(MLPATH.4)-(MLPATH.6) sharpen its criterion to

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^alpha
 (
   |sqrt(K)r_K-tau_(K,0)|^2
   +sum_(1<=r<R_K)|tau_(K,r)|^2
 )
<infinity.                                            (MLPATH.9)
```

For every member of any one fixed cofinal sequence `alpha_j->0`,
(MLPATH.9) is equivalent to RH by the checked reciprocal-tail
criterion. It has two genuinely different pieces:

```text
affine mean defect:
  sqrt(K)r_K-tau_(K,0);

nonconstant local-path spectrum:
  tau_(K,r),                         1<=r<ceil(sqrt(K)).
```

The first contains the incoming global tail and local path mean. The
second contains the smooth shape of the anchored path. Controlling one
does not algebraically control the other.

## Transfer To Ordinary Local Mertens Increments

Define

```text
S_(K,m):=M(K+m)-M(K)
        =sum_(ell=1)^m mu(K+ell),

a_m:=(K+m)^(-1-alpha),
b_m:=1/a_m.                                           (MLPATH.10)
```

Finite Abel summation and its inverse give

```text
T_(K,m)
 =a_m S_(K,m)
  +sum_(ell=1)^(m-1)
    S_(K,ell)(a_ell-a_(ell+1)),                       (MLPATH.11)

S_(K,m)
 =b_m T_(K,m)
  +sum_(ell=1)^(m-1)
    T_(K,ell)(b_ell-b_(ell+1)).                       (MLPATH.12)
```

For the vectors indexed by `1<=m<K`, the strict-prefix matrix has
operator norm at most `K`. Since

```text
a_m<=K^(-1-alpha),
a_ell-a_(ell+1)<=(1+alpha)K^(-2-alpha),

b_m<=(2K)^(1+alpha),
b_(ell+1)-b_ell<=(1+alpha)(2K)^alpha,
```

(MLPATH.11)-(MLPATH.12) imply the explicit full-path comparison

```text
[2^alpha(3+alpha)]^(-1)
 K^(-1-alpha)||S_K||_2
 <=||T_K||_2
 <=(2+alpha)K^(-1-alpha)||S_K||_2.                   (MLPATH.13)
```

This does not compare their first `sqrt(K)` DCT projections. The Abel
matrix is triangular and does not commute with the low-mode projector.
Using (MLPATH.13) after a generic Cauchy-Schwarz step therefore recycles
the full Mertens-energy target rather than proving the missing
low-mode gain.

## Zero-Mean First-Mode Countermodel

The following exact scalar construction proves that coefficient size,
absolute convergence, endpoint-tail cancellation, and mean cancellation
cannot close the nonconstant modes generically.

For this fixed `alpha` and every dyadic `K>=4`, put

```text
A_K:=2^(-2-alpha)pi^(-1)K^(-alpha),

T*_(K,m):=A_K sin(2pi*m/K),                 0<=m<=K,

q*_(K+m):=T*_(K,m)-T*_(K,m-1),             1<=m<=K.
                                                            (MLPATH.14)
```

Set `q*_1=0`. The annuli `(K,2K]` partition all remaining indices. For
`n=K+m`,

```text
|q*_n|
 <=2pi*A_K/K
 =2^(-1-alpha)K^(-1-alpha)
 <=n^(-1-alpha).                                    (MLPATH.15)
```

The block `l1` norm is at most `2pi*A_K`, so the global sequence is
absolutely summable. Moreover,

```text
sum_(K<n<=2K)q*_n=T*_(K,K)-T*_(K,0)=0.
```

All later block totals also vanish. Hence, at every dyadic endpoint,

```text
r*_K=0,
tau*_(K,0)=0,
c*_(K,0)=0.                                          (MLPATH.16)
```

Nevertheless the first nonconstant coefficient is

```text
|c*_(K,1)|
 =A_K sqrt(2/K) cos(pi/K)/2
  [csc(pi/(2K))+csc(3pi/(2K))].                      (MLPATH.17)
```

Indeed, product-to-sum gives

```text
sum_(m=0)^(K-1)
 sin(2pi*m/K)cos(pi*(m+1/2)/K)

 =cos(pi/K)/2
  [csc(pi/(2K))+csc(3pi/(2K))].
```

For `K>=4`, `cos(pi/K)>=sqrt(2)/2` and `sin x<=x`, so

```text
|c*_(K,1)|>=A_K sqrt(K)/pi,

K^alpha|c*_(K,1)|^2
 >=2^(-4-2alpha)pi^(-4)K^(1-alpha).                  (MLPATH.18)
```

The dyadic series diverges although every dyadic tail endpoint and
every mean mode vanish. This is deliberately not a Mobius model. It
proves a method barrier: the actual proof must use arithmetic signed
cancellation in the nonconstant local path, not only its coefficient
envelope or the scalar tail level.

## Source And Method Audit

Davenport's theorem gives, uniformly in `theta`,

```text
sum_(n<=X)mu(n)e(n theta)<<_A X(log X)^(-A).
```

That arbitrary-log estimate is power-short for (MLPATH.9). The focused
primary-source audit also checked:

- R. C. Baker and G. Harman, *Exponential Sums Formed with the Mobius
  Function*:
  https://doi.org/10.1112/jlms/s2-43.2.193
- H. Maier and A. Sankaranarayanan, *On an Exponential Sum Involving
  the Mobius Function*:
  https://doi.org/10.46298/hrj.2005.151
- W. Zhang, *On an Exponential Sum Related to the Mobius Function*:
  https://arxiv.org/abs/2204.04613

Zhang records the unconditional Davenport estimate and quotes
Baker-Harman power estimates under the hypothesis that every Dirichlet
`L(s,chi)` is zero-free in `Re(s)>a`. At the endpoint `a=1/2`, the
quoted uniform exponent is `3/4+epsilon`. Lemma 11.22V requires

```text
sigma<(1+alpha)/2,
```

so `sigma=3/4+epsilon` closes only
`alpha>1/2+2epsilon`, never a cofinal `alpha->0` sequence. The
unconditional rational estimate of Maier-Sankaranarayanan, as quoted
in Zhang's Lemma 2.1, retains Dirichlet-L zero-density terms. The audit
found no source-backed simplification yielding the required uniform
cofinal-closing power. This is a focused search report, not a novelty,
priority, or exhaustive-literature claim.

A viable Vaughan or bilinear attack must therefore enter before
rowwise absolute values, preserve signed coupling over the local modes
and dyadic scales, and treat the affine mean defect at the same time.
Translation-averaged and almost-all interval theorems need an explicit
de-averaging step. Estimates only for nonzero phases do not treat the
mean defect. Finally, (MLPATH.9) concerns `R_alpha`, not the complete
Burnol energy `Q`.

The open gate is:

```text
for every alpha in one fixed cofinal sequence alpha_j->0,

sum_(K dyadic)K^alpha
 (
   |sqrt(K)r_K-tau_(K,0)|^2
   +sum_(1<=r<ceil(sqrt(K)))|tau_(K,r)|^2
 )
<infinity.                                            (MLPATH.19)
```

This is RH-equivalent. The exact transference and countermodel sharpen
where new arithmetic input must act, but supply no such gain. RH,
PF-infinity, and `Lambda <= 0` remain open.
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
        "built Mertens local-path cosine transference: "
        f"{len(payload['rows'])} rows, "
        f"{payload['audit']['exact_reduction_count']} exact reductions, "
        f"{payload['audit']['proof_guard_count']} proof guards, "
        f"{payload['audit']['open_two_component_gate_count']} open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
