#!/usr/bin/env python3
"""Build the affine-tent and Brownian-bridge unified Vaughan handoff."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_affine_tent_bridge_handoff.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_affine_tent_bridge_handoff.md"
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
            "matbh_01_tail_setup",
            "exact_definition",
            "available_exact",
            "For 0<alpha<1 set q_n=mu(n)n^(-1-alpha), "
            "r_N=sum_(n>N)q_n, and "
            "R_alpha=sum_(N>=1)N^alpha|r_N|^2.",
            "These are the reciprocal-Mobius tail coordinates from "
            "Lemmas 11.22T-X.",
        ),
        row(
            "matbh_02_affine_block_sum",
            "exact_identity",
            "available_exact",
            "For dyadic K, u_K:=sum_(N=K)^(2K-1)r_N="
            "sqrt(K)c_(K,0).",
            "The affine mean coefficient is exactly one dyadic block sum "
            "of the tail sequence.",
        ),
        row(
            "matbh_03_infinite_mean_kernel",
            "exact_identity",
            "available_exact",
            "u_K=sum_n q_n w_K(n), where "
            "w_K(n)=min((n-K)_+,K).",
            "Count the block tails containing q_n; w_K is constant for "
            "n>=2K.",
        ),
        row(
            "matbh_04_mean_energy",
            "exact_identity",
            "available_exact",
            "M_alpha=sum_(K dyadic)K^(alpha-1)|u_K|^2.",
            "This is Lemma 11.22X.15 after u_K=sqrt(K)c_(K,0).",
        ),
        row(
            "matbh_05_scale_filter",
            "exact_definition",
            "available_exact",
            "Define v_K:=u_K-(1/2)u_(2K).",
            "The coefficient 1/2 is forced by cancellation of the two "
            "constant kernel plateaux beyond 4K.",
        ),
        row(
            "matbh_06_tent_kernel",
            "exact_definition",
            "available_exact",
            "Define W_K(n)=n-K for K<n<=2K, "
            "W_K(n)=(4K-n)/2 for 2K<n<4K, and W_K(n)=0 otherwise.",
            "W_K is a nonnegative compact tent with maximum K.",
        ),
        row(
            "matbh_07_compact_localization",
            "exact_identity",
            "available_exact",
            "v_K=sum_(K<n<4K)mu(n)n^(-1-alpha)W_K(n).",
            "This follows from W_K=w_K-(1/2)w_(2K); the infinite tail "
            "cancels exactly.",
        ),
        row(
            "matbh_08_normalized_sequences",
            "exact_definition",
            "available_exact",
            "Put p=(alpha-1)/2, d_K=K^p u_K, and e_K=K^p v_K.",
            "Then M_alpha=sum_K|d_K|^2 and the filtered affine energy is "
            "E_aff=sum_K|e_K|^2.",
        ),
        row(
            "matbh_09_filter_multiplier",
            "exact_definition",
            "available_exact",
            "Let rho_alpha=2^(-(1+alpha)/2), so 0<rho_alpha<1.",
            "The strict contraction is uniform over the dyadic scale "
            "index for fixed alpha.",
        ),
        row(
            "matbh_10_normalized_filter_identity",
            "exact_identity",
            "available_exact",
            "e_K=d_K-rho_alpha*d_(2K).",
            "Use rho_alpha*2^p=1/2.",
        ),
        row(
            "matbh_11_boundary_decay",
            "exact_inequality",
            "available_exact",
            "The coefficient envelope gives |d_K|=O_alpha("
            "K^((1-alpha)/2)) and hence "
            "rho_alpha^m d_(2^m K)=O_alpha(2^(-alpha*m)).",
            "This supplies the boundary condition needed to invert the "
            "forward scale filter without assuming the target energy.",
        ),
        row(
            "matbh_12_filter_inverse",
            "exact_identity",
            "available_exact",
            "d_K=sum_(m>=0)rho_alpha^m e_(2^m K).",
            "Iterate row 10 and use the boundary decay in row 11.",
        ),
        row(
            "matbh_13_filter_norm_bounds",
            "exact_inequality",
            "available_exact",
            "(1-rho_alpha)||d||_2<=||e||_2<="
            "(1+rho_alpha)||d||_2.",
            "The upper bound is the shift triangle inequality; the lower "
            "bound follows from the geometric-series inverse.",
        ),
        row(
            "matbh_14_affine_energy_equivalence",
            "exact_equivalence",
            "available_exact",
            "M_alpha<infinity iff "
            "E_aff:=sum_(K dyadic)|e_K|^2<infinity.",
            "The infinite affine mean kernel has therefore been replaced "
            "by compact tent tests without estimating it.",
        ),
        row(
            "matbh_15_affine_signed_expansion",
            "exact_identity",
            "available_exact",
            "|e_K|^2=D_aff,K+2C_aff,K, where "
            "D_aff,K=K^(alpha-1)sum_n mu(n)^2 n^(-2-2alpha)W_K(n)^2 "
            "and C_aff,K is the corresponding sum over K<n<m<4K.",
            "This is the diagonal/off-diagonal expansion of the compact "
            "tent square.",
        ),
        row(
            "matbh_16_affine_diagonal_bound",
            "exact_inequality",
            "available_exact",
            "0<=D_aff,K<=3K^(-alpha), so "
            "sum_(K dyadic)D_aff,K<infinity.",
            "There are fewer than 3K terms, W_K<=K, and n>K.",
        ),
        row(
            "matbh_17_affine_offdiagonal_criterion",
            "exact_equivalence",
            "available_exact",
            "M_alpha<infinity iff "
            "sup_J sum_(K<=2^J)C_aff,K<infinity.",
            "Use rows 14-16 and monotonicity of the partial filtered "
            "affine energies.",
        ),
        row(
            "matbh_18_bridge_energy_recall",
            "exact_definition",
            "available_exact",
            "Let E_br,K=K^(-2-alpha)||[S_K]||_K^2, where "
            "S_(K,m)=M(K+m)-M(K) and [.] is the quotient modulo constants.",
            "The sum of E_br,K is L_alpha from Lemma 11.22X.",
        ),
        row(
            "matbh_19_bridge_signed_expansion",
            "exact_identity",
            "available_exact",
            "E_br,K=D_br,K+2C_br,K, with "
            "D_br,K=K^(-2-alpha)Delta_K and "
            "C_br,K=K^(-2-alpha)C_K.",
            "Delta_K and C_K are the Brownian-bridge diagonal and "
            "off-diagonal terms from Lemma 11.22X.",
        ),
        row(
            "matbh_20_bridge_diagonal_bound",
            "exact_inequality",
            "available_exact",
            "0<=D_br,K<=K^(-alpha)/6.",
            "Use Delta_K<=(K^2-1)/6.",
        ),
        row(
            "matbh_21_unified_local_energy",
            "exact_definition",
            "available_exact",
            "Define E_loc=sum_(K dyadic)(|e_K|^2+E_br,K).",
            "Rows 14 and 18 give R_alpha<infinity iff "
            "E_loc<infinity.",
        ),
        row(
            "matbh_22_combined_shift_kernel",
            "exact_definition",
            "available_exact",
            "For h>=1 define G_(K,h)(n) as "
            "K^(alpha-1)[n(n+h)]^(-1-alpha)W_K(n)W_K(n+h) "
            "on K<n<n+h<4K, plus "
            "K^(-3-alpha)(n-K)(2K-n-h) "
            "on K<n<n+h<2K.",
            "The first summand is the affine tent kernel and the second "
            "is the normalized Brownian-bridge kernel.",
        ),
        row(
            "matbh_23_combined_support",
            "exact_identity",
            "available_exact",
            "G_(K,h)(n)=0 unless K<n<n+h<4K.",
            "Both formerly separate open components are now local on at "
            "most two adjacent dyadic annuli.",
        ),
        row(
            "matbh_24_kernel_amplitude",
            "exact_inequality",
            "available_exact",
            "0<=G_(K,h)(n)<=(5/4)K^(-1-alpha).",
            "The tent product contributes at most K^(-1-alpha); the "
            "bridge product contributes at most K^(-1-alpha)/4.",
        ),
        row(
            "matbh_25_unified_signed_expansion",
            "exact_identity",
            "available_exact",
            "E_loc(J)=D_loc(J)+2O_loc(J), where "
            "O_loc(J)=sum_(K<=2^J)sum_h sum_n "
            "mu(n)mu(n+h)G_(K,h)(n).",
            "Combine rows 15 and 19 before taking a cutoff limit.",
        ),
        row(
            "matbh_26_unified_diagonal_bound",
            "exact_inequality",
            "available_exact",
            "The per-scale unified diagonal is at most "
            "(19/6)K^(-alpha), hence D_loc(infinity)<infinity.",
            "Add the bounds 3K^(-alpha) and K^(-alpha)/6.",
        ),
        row(
            "matbh_27_unified_signed_criterion",
            "exact_equivalence",
            "available_exact",
            "R_alpha<infinity iff sup_J O_loc(J)<infinity.",
            "The partial local energies are nonnegative and monotone, and "
            "their diagonal converges.",
        ),
        row(
            "matbh_28_two_interval_split",
            "exact_identity",
            "available_exact",
            "Split each n-sum in O_loc into K<n<=2K and "
            "2K<n<4K.",
            "This places every compact test into the standard "
            "(X,2X] form with X=K or X=2K.",
        ),
        row(
            "matbh_29_vaughan_source",
            "literature_guard",
            "source_backed",
            "Green and Tao, Lemma 4.1, supply the finite Vaughan "
            "Type I/II identity on each interval (X,2X].",
            "Only the exact divisor identity is imported; no estimate "
            "for the endogenous shifted-Mobius tests is attributed to it.",
        ),
        row(
            "matbh_30_vaughan_tests",
            "exact_definition",
            "available_exact",
            "For interval a in {0,1}, define "
            "g_(K,h,a)(n)=mu(n+h)G_(K,h)(n) restricted to "
            "(2^a K,2^(a+1)K].",
            "Summing mu(n)g_(K,h,a)(n) over h and a recovers the "
            "unified off-diagonal O_loc scale by scale.",
        ),
        row(
            "matbh_31_exact_vaughan_handoff",
            "exact_identity",
            "available_exact",
            "With the standard finite Vaughan coefficients on each "
            "interval, O_loc,K=-TI_loc,K+TII_loc,K.",
            "Apply the identity termwise in h and in the two interval "
            "pieces, then retain the signed difference.",
        ),
        row(
            "matbh_32_final_type_i_ii_criterion",
            "exact_equivalence",
            "available_exact",
            "R_alpha<infinity iff "
            "sup_J sum_(K<=2^J)(-TI_loc,K+TII_loc,K)<infinity.",
            "This is one compact signed Type I/II target for both the "
            "affine mean and centered bridge.",
        ),
        row(
            "matbh_33_absolute_value_guard",
            "proof_guard",
            "guard_validated",
            "Replacing -TI_loc,K+TII_loc,K by the sum of separate "
            "absolute values is not an equivalent reduction.",
            "There are O(K^2) shift/base pairs of amplitude "
            "O(K^(-1-alpha)), giving the divergent scale K^(1-alpha).",
        ),
        row(
            "matbh_34_averaged_chowla_guard",
            "literature_guard",
            "source_backed",
            "The natural averaged-Chowla o(K^2) pair scale becomes only "
            "o(K^(1-alpha)) after the normalized kernel, one full power "
            "above the diagonal K^(-alpha) calibration.",
            "This is a power audit of the published averaged theorem, "
            "not a lower bound on the actual unified correlation.",
        ),
        row(
            "matbh_35_localization_scope_guard",
            "proof_guard",
            "guard_validated",
            "The invertible scale filter localizes M_alpha but supplies "
            "no Mobius cancellation estimate.",
            "Compact support must not be reported as a Type I/II gain, "
            "a full Burnol bound, or a proof of RH.",
        ),
        row(
            "matbh_36_commutator_retirement_guard",
            "proof_guard",
            "guard_validated",
            "The Abel/low-DCT commutator is not a remaining obligation: "
            "the centered quotient handles nonconstant modes and the "
            "scale filter handles the affine mean.",
            "Future work should attack the unified signed kernel rather "
            "than reintroduce the retired commutator.",
        ),
        row(
            "matbh_37_finite_validation",
            "finite_reproduction",
            "finite_check_passed",
            "Arbitrary finite coefficient tests verify the tail count, "
            "tent cancellation, normalized filter, inverse truncations, "
            "and norm bounds; Mobius tests verify both signed expansions "
            "and the two-interval Vaughan residual.",
            "Finite validation checks the algebra only and is not "
            "asymptotic evidence for the open signed estimate.",
        ),
        row(
            "matbh_38_weighted_bridge_coordinate",
            "exact_definition",
            "available_exact",
            "Let T_(K,m)=sum_(i=1)^m q_(K+i), 0<=m<K, and retain "
            "the weighted centered energy K^alpha||[T_K]||_K^2.",
            "Together with the affine mean this is the exact DCT "
            "decomposition before transfer to ordinary Mertens increments.",
        ),
        row(
            "matbh_39_future_anchor",
            "exact_definition",
            "available_exact",
            "Define Z_K=Kq_(2K)+sum_(2K<n<4K)(4K-n)q_n/2; then "
            "v_K=sum_(i=1)^(K-1)i*q_(K+i)+Z_K.",
            "The q_(2K) endpoint belongs to the affine tent and not to "
            "the length-K centered path.",
        ),
        row(
            "matbh_40_constant_direction_completion",
            "exact_identity",
            "available_exact",
            "For 1<=i,j<K, K^(alpha-1)ij+"
            "K^alpha[min(i,j)-ij/K]=K^alpha min(i,j).",
            "The affine rank-one square exactly restores the constant "
            "direction removed by the Brownian bridge.",
        ),
        row(
            "matbh_41_tail_lattice_square",
            "exact_identity",
            "available_exact",
            "|e_K|^2+K^alpha||[T_K]||_K^2="
            "K^alpha sum_(t=1)^K|Z_K/K+"
            "sum_(i=t)^(K-1)q_(K+i)|^2.",
            "This is finite variance completion, including the terminal "
            "row t=K whose suffix sum is zero.",
        ),
        row(
            "matbh_42_positive_compact_criterion",
            "exact_equivalence",
            "available_exact",
            "R_alpha<infinity iff sum_(K dyadic)K^alpha "
            "sum_(t=1)^K|Z_K/K+sum_(i=t)^(K-1)q_(K+i)|^2<infinity.",
            "Use dyadic tail-energy comparison, DCT Parseval, and the "
            "invertible affine scale filter.",
        ),
        row(
            "matbh_43_positive_square_scope",
            "exact_identity",
            "available_exact",
            "Every row in the positive criterion uses only q_n with "
            "K<n<4K.",
            "The Cholesky form is compact and nonnegative but does not "
            "bound its own Mobius suffix sums.",
        ),
        row(
            "matbh_44_ordinary_future_anchor",
            "exact_definition",
            "available_exact",
            "Put s=1+alpha, z_K=Z_K/K, and beta_K=(2K)^s*z_K; "
            "then beta_K=mu(2K)+sum_(2K<n<4K)"
            "((2K)/n)^s(4K-n)/(2K)*mu(n).",
            "This is one explicit compact future anchor, not an "
            "unproved short-interval cancellation estimate.",
        ),
        row(
            "matbh_45_bounded_anchor_weights",
            "exact_bound",
            "available_exact",
            "Every coefficient of mu(n) in beta_K lies in [0,1] and "
            "the support is 2K<=n<4K.",
            "Coefficient boundedness gives only the trivial "
            "|beta_K|=O(K); it does not itself close the criterion.",
        ),
        row(
            "matbh_46_full_suffix_abel_pair",
            "exact_identity",
            "available_exact",
            "With a_i=(K+i)^(-s) for i<K, a_K=(2K)^(-s), "
            "A_(K,t)=sum_(i=t)^K a_i*x_i and "
            "B_(K,t)=sum_(i=t)^K x_i, where x_K=beta_K, finite Abel "
            "summation gives mutually inverse triangular maps.",
            "The terminal choice a_K*beta_K=z_K is essential; omitting "
            "the future anchor would break the exact inverse.",
        ),
        row(
            "matbh_47_full_suffix_norm_transfer",
            "exact_bound",
            "available_exact",
            "(4sqrt(2))^(-1)K^(-s)||B_K||_2<="
            "||A_K||_2<=sqrt(6)K^(-s)||B_K||_2.",
            "The constants follow from explicit row/column sums for "
            "the Abel matrix and its inverse and are uniform for "
            "0<alpha<1.",
        ),
        row(
            "matbh_48_ordinary_anchored_criterion",
            "exact_equivalence",
            "available_exact",
            "R_alpha<infinity iff sum_(K dyadic)K^(-2-alpha)"
            "sum_(t=1)^K|beta_K+sum_(i=t)^(K-1)mu(K+i)|^2<infinity.",
            "This removes all power weights from the current-block "
            "suffixes while preserving one bounded-coefficient future "
            "anchor.",
        ),
        row(
            "matbh_49_ordinary_calibration",
            "sufficient_condition",
            "available_exact",
            "If H_K:=sum_(t=1)^K|beta_K+"
            "sum_(i=t)^(K-1)mu(K+i)|^2=O_epsilon(K^(2+epsilon)) "
            "for some epsilon<alpha, then the alpha criterion holds.",
            "The coefficientwise estimate gives only H_K=O(K^3), "
            "which is one full power short of this calibration.",
        ),
        row(
            "matbh_50_all_interval_guard",
            "literature_guard",
            "source_backed",
            "Matomaki and Teravainen prove o(H) cancellation for every "
            "interval of length H=x^theta when theta>0.55.",
            "An o(H) pointwise bound has no fixed square-root-scale power "
            "and does not make the K-term suffix-square family summable.",
        ),
        row(
            "matbh_51_almost_all_guard",
            "literature_guard",
            "source_backed",
            "Matomaki-Radziwill and the higher-uniformity short-interval "
            "theorems give strong cancellation outside small exceptional "
            "sets of base points.",
            "The present endpoints are one fixed dyadic sequence; no "
            "source-backed de-averaging theorem excludes every one of "
            "those points from the exceptional sets.",
        ),
        row(
            "matbh_52_positive_alpha_guard",
            "proof_guard",
            "guard_validated",
            "All criteria are imposed separately for fixed alpha>0; "
            "the diagonal majorant is sum_(K dyadic)K^(-alpha).",
            "At alpha=0 this majorant is not summable. Cofinality "
            "alpha_j->0 does not permit setting alpha=0 or demanding a "
            "uniform limit in alpha.",
        ),
        row(
            "matbh_53_open_signed_gate",
            "open_theorem_target",
            "open",
            "Prove the unified signed criterion in row 32 for every "
            "member of one fixed cofinal sequence alpha_j->0 without "
            "assuming RH.",
            "No signed Type I/II power gain, full Burnol bound, RH, "
            "PF-infinity, or Lambda<=0 conclusion is supplied here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_mertens_affine_tent_bridge_handoff",
        "date": "2026-07-23",
        "status": (
            "exact affine scale-localization and unified local Vaughan "
            "reduction with one open signed gate"
        ),
        "parameters": {
            "alpha_range": "0<alpha<1",
            "dyadic_blocks": "K=2^j",
            "scale_filter": "u_K-(1/2)u_(2K)",
            "filter_contraction": "rho_alpha=2^(-(1+alpha)/2)",
            "combined_support": "K<n<n+h<4K",
        },
        "source_anchors": [
            "https://doi.org/10.5802/aif.2401",
            "https://arxiv.org/abs/1503.05121",
            "https://doi.org/10.2140/ant.2015.9.2167",
            "https://arxiv.org/abs/1911.09076",
            "https://arxiv.org/abs/1501.04585",
            "https://doi.org/10.1007/s00222-026-01408-6",
        ],
        "rows": rows,
        "audit": {
            "row_count": 53,
            "exact_reduction_count": 44,
            "literature_guard_count": 4,
            "proof_guard_count": 4,
            "open_signed_gate_count": 1,
            "affine_infinite_tail_localized": True,
            "scale_filter_invertible": True,
            "affine_bridge_kernel_unified": True,
            "combined_diagonal_summable": True,
            "two_interval_vaughan_handoff_proved": True,
            "positive_tail_lattice_square_proved": True,
            "ordinary_mobius_anchored_criterion_proved": True,
            "alpha_zero_passage_allowed": False,
            "averaged_chowla_closes_gate": False,
            "signed_type_i_ii_gain_proved": False,
            "full_burnol_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Mertens Affine-Tent/Bridge Handoff

Date: 2026-07-23

Status: exact affine scale-localization and unified local Vaughan
reduction with one open signed gate. This is not a proof of RH,
PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_affine_tent_bridge_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_affine_tent_bridge_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_affine_tent_bridge_handoff.py
```

## Affine Mean As A Scale Difference

Fix `0<alpha<1` and set

```text
q_n:=mu(n)n^(-1-alpha),
r_N:=sum_(n>N)q_n,
R_alpha:=sum_(N>=1)N^alpha|r_N|^2.
```

For dyadic `K`, write

```text
u_K
 :=sum_(N=K)^(2K-1)r_N
  =sqrt(K)c_(K,0).                                    (MATBH.1)
```

The affine mean series from Lemma 11.22X is

```text
M_alpha
 =sum_(K dyadic)K^(alpha-1)|u_K|^2.                   (MATBH.2)
```

Counting the tails containing each coefficient gives

```text
u_K=sum_n q_n w_K(n),

w_K(n):=min((n-K)_+,K).                               (MATBH.3)
```

The plateau of `w_K` is the nonlocal part. Apply the forced scale
difference

```text
v_K:=u_K-(1/2)u_(2K).                                 (MATBH.4)
```

Define the compact tent

```text
W_K(n)
 :=
  n-K,          K<n<=2K,
  (4K-n)/2,     2K<n<4K,
  0,            otherwise.                            (MATBH.5)
```

Since `W_K=w_K-(1/2)w_(2K)`, the infinite tail cancels exactly:

```text
v_K
 =sum_(K<n<4K)mu(n)n^(-1-alpha)W_K(n).                (MATBH.6)
```

Put

```text
p:=(alpha-1)/2,
d_K:=K^p u_K,
e_K:=K^p v_K,
rho_alpha:=2^(-(1+alpha)/2).                          (MATBH.7)
```

Then

```text
e_K=d_K-rho_alpha*d_(2K).                             (MATBH.8)
```

The coefficient envelope gives

```text
|d_K|=O_alpha(K^((1-alpha)/2)),

rho_alpha^m d_(2^m K)=O_alpha(2^(-alpha*m))->0.       (MATBH.9)
```

Thus iteration of (MATBH.8) is legitimate without assuming the target:

```text
d_K=sum_(m>=0)rho_alpha^m e_(2^m K).                 (MATBH.10)
```

The unilateral shift and its geometric inverse give

```text
(1-rho_alpha)||d||_2
 <=||e||_2
 <=(1+rho_alpha)||d||_2.                              (MATBH.11)
```

Consequently

```text
M_alpha<infinity
 iff
E_aff:=sum_(K dyadic)|e_K|^2<infinity.                (MATBH.12)
```

This is an exact localization, not an estimate for `M_alpha`.

## Signed Affine Tent

Expanding the compact square,

```text
|e_K|^2=D_(aff,K)+2C_(aff,K),                         (MATBH.13)

D_(aff,K)
 :=K^(alpha-1)
   sum_(K<n<4K)
   mu(n)^2 n^(-2-2alpha)W_K(n)^2,

C_(aff,K)
 :=K^(alpha-1)
   sum_(K<n<m<4K)
   mu(n)mu(m)(nm)^(-1-alpha)W_K(n)W_K(m).             (MATBH.14)
```

Since there are fewer than `3K` terms and `W_K<=K`,

```text
0<=D_(aff,K)<=3K^(-alpha),

sum_(K dyadic)D_(aff,K)<infinity.                     (MATBH.15)
```

It follows that

```text
M_alpha<infinity
 iff
sup_J sum_(K dyadic<=2^J)C_(aff,K)<infinity.          (MATBH.16)
```

The affine mean has therefore become a compact signed off-diagonal
problem on `(K,4K)`.

## Unified Local Kernel

Retain the Brownian-bridge energy from Lemma 11.22X:

```text
E_(br,K)
 :=K^(-2-alpha)||[S_K]||_K^2
  =D_(br,K)+2C_(br,K),                                (MATBH.17)

D_(br,K):=K^(-2-alpha)Delta_K,
C_(br,K):=K^(-2-alpha)C_K.
```

Its diagonal satisfies

```text
0<=D_(br,K)<=K^(-alpha)/6.                            (MATBH.18)
```

Define

```text
E_loc
 :=sum_(K dyadic)(|e_K|^2+E_(br,K)).                  (MATBH.19)
```

Lemmas 11.22X and (MATBH.12) give the exact equivalence

```text
R_alpha<infinity iff E_loc<infinity.                  (MATBH.20)
```

For `h>=1`, define

```text
G_(K,h)(n)
 :=
 K^(alpha-1)[n(n+h)]^(-1-alpha)
 W_K(n)W_K(n+h)
 1_(K<n<n+h<4K)

 +K^(-3-alpha)(n-K)(2K-n-h)
 1_(K<n<n+h<2K).                                     (MATBH.21)
```

The first line is the affine-tent kernel and the second is the
normalized Brownian-bridge kernel. They have the same scale:

```text
0<=G_(K,h)(n)<=(5/4)K^(-1-alpha),                     (MATBH.22)
```

and both vanish unless `K<n<n+h<4K`. For a finite dyadic cutoff,

```text
E_loc(J)=D_loc(J)+2O_loc(J),                          (MATBH.23)

O_loc(J)
 :=sum_(K<=2^J)
   sum_h sum_n
   mu(n)mu(n+h)G_(K,h)(n).
```

The diagonal obeys

```text
D_(loc,K)
 <=(3+1/6)K^(-alpha)
 =(19/6)K^(-alpha),                                  (MATBH.24)
```

so it is dyadically summable. Hence

```text
R_alpha<infinity
 iff
sup_J O_loc(J)<infinity.                              (MATBH.25)
```

This replaces the two open pieces of Lemma 11.22X by one compact signed
correlation criterion.

## Two-Interval Vaughan Handoff

Split the base variable into

```text
K<n<=2K,
2K<n<4K.                                              (MATBH.26)
```

For `a in {0,1}`, let

```text
g_(K,h,a)(n)
 :=mu(n+h)G_(K,h)(n)
   1_(2^a K<n<=2^(a+1)K).                             (MATBH.27)
```

Each piece now has the standard `(X,2X]` form. Applying the finite
Vaughan identity of Green and Tao separately to the two intervals and
then summing over `h` gives

```text
O_(loc,K)=-TI_(loc,K)+TII_(loc,K).                    (MATBH.28)
```

Therefore the single live criterion is

```text
sup_J sum_(K dyadic<=2^J)
 (-TI_(loc,K)+TII_(loc,K))
 <infinity.                                           (MATBH.29)
```

The signed difference is essential. Separate absolute values see
`O(K^2)` shift/base pairs of size `O(K^(-1-alpha))` and restore the
divergent scale `K^(1-alpha)`. The natural `o(K^2)` averaged-Chowla
pair scale likewise gives only `o(K^(1-alpha))`, one full power above
the summable diagonal calibration `K^(-alpha)`.

## Positive Compact Tail Square

There is a simpler Cholesky form before transferring the centered path
to ordinary Mertens increments. Retain

```text
T_(K,m)
 :=sum_(i=1)^m q_(K+i),                    0<=m<K,

Z_K
 :=K*q_(2K)
   +sum_(2K<n<4K)(4K-n)q_n/2.                         (MATBH.30)
```

The endpoint `q_(2K)` belongs to the affine tent, while the centered
length-`K` path uses only `q_(K+1),...,q_(2K-1)`. Thus

```text
v_K
 =sum_(i=1)^(K-1)i*q_(K+i)+Z_K.                       (MATBH.31)
```

The affine rank-one kernel exactly fills the constant direction removed
by the Brownian bridge:

```text
K^(alpha-1)i*j
 +K^alpha[min(i,j)-i*j/K]
 =K^alpha min(i,j),                 1<=i,j<K.          (MATBH.32)
```

Finite variance completion therefore gives the exact positive identity

```text
|e_K|^2+K^alpha||[T_K]||_K^2

 =K^alpha sum_(t=1)^K
   |Z_K/K+sum_(i=t)^(K-1)q_(K+i)|^2.                  (MATBH.33)
```

For `t=K` the suffix sum is empty. DCT Parseval, dyadic weight
comparison, and the invertible scale filter now yield the alternative
one-gate criterion

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^alpha
 sum_(t=1)^K
 |Z_K/K+sum_(i=t)^(K-1)q_(K+i)|^2
 <infinity.                                           (MATBH.34)
```

Every coefficient in (MATBH.34) lies in `(K,4K)`. This positive local
tail-lattice square is equivalent to the signed-kernel target, but it
does not estimate the Mobius suffixes.

## Ordinary-Mobius Anchored Form

The power weights on the current-block suffixes can be removed exactly.
Put `s=1+alpha`, `z_K=Z_K/K`, and

```text
beta_K
 :=(2K)^s*z_K
  =mu(2K)
   +sum_(2K<n<4K)
     ((2K)/n)^s*(4K-n)/(2K)*mu(n).                   (MATBH.35)
```

Every coefficient in `beta_K` lies in `[0,1]`, and its support is
`2K<=n<4K`. For `1<=i<K` set

```text
a_i:=(K+i)^(-s),       x_i:=mu(K+i),
a_K:=(2K)^(-s),        x_K:=beta_K,

A_(K,t):=sum_(i=t)^K a_i*x_i
        =z_K+sum_(i=t)^(K-1)(K+i)^(-s)mu(K+i),

B_(K,t):=sum_(i=t)^K x_i
        =beta_K+sum_(i=t)^(K-1)mu(K+i).              (MATBH.36)
```

Writing `A_K` and `B_K` for these length-`K` vectors, finite Abel
summation gives the mutually inverse triangular maps

```text
A_(K,t)
 =a_t*B_(K,t)
  +sum_(i=t+1)^K(a_i-a_(i-1))*B_(K,i),

B_(K,t)
 =a_t^(-1)*A_(K,t)
  +sum_(i=t+1)^K(a_i^(-1)-a_(i-1)^(-1))*A_(K,i).
                                                            (MATBH.37)
```

The first matrix has absolute row sums at most `2K^(-s)` and column
sums at most `3K^(-s)`. Its inverse has absolute row sums at most
`4K^s` and column sums at most `8K^s`. Schur's test therefore proves
the uniform two-sided estimate

```text
(1/(4sqrt(2)))*K^(-s)||B_K||_2
 <=||A_K||_2
 <=sqrt(6)*K^(-s)||B_K||_2.                          (MATBH.38)
```

Combining (MATBH.34) with (MATBH.38) yields the ordinary-Mobius
anchored criterion

```text
R_alpha<infinity
 iff
sum_(K dyadic)K^(-2-alpha)
 sum_(t=1)^K
 |beta_K+sum_(i=t)^(K-1)mu(K+i)|^2
 <infinity.                                           (MATBH.39)
```

Thus, if

```text
H_K
 :=sum_(t=1)^K
   |beta_K+sum_(i=t)^(K-1)mu(K+i)|^2
 =O_epsilon(K^(2+epsilon)),       epsilon<alpha,      (MATBH.40)
```

then the criterion holds at that `alpha`. The coefficientwise bound
gives only `H_K=O(K^3)`, so cancellation sufficient to save almost one
full power is still required.

## Quantifier And Critical Guards

The strongest directly relevant short-interval theorems do not have the
required calibration. Matomaki and Teravainen prove `o(H)` cancellation
in every interval of length `H=x^theta` for `theta>0.55`; an unquantified
`o(H)` bound is still far above the square-root scale needed after
summing the suffix squares. Matomaki-Radziwill and the newer
higher-uniformity results give strong logarithmic cancellation outside
small exceptional sets of base points. The present family uses one fixed
dyadic base point at each scale, and no cited theorem de-averages those
exceptional sets along that sequence.

The parameter boundary is also strict:

```text
sum_(K dyadic)K^(-alpha)<infinity
```

holds for every fixed `alpha>0` and fails at `alpha=0`. The cofinal
criterion asks separately for each positive `alpha_j->0`; it does not
permit setting `alpha=0` or passing to a bound uniform in `alpha`.

## Open Gate

The scale filter retires the nonlocal affine tail and the quotient
retires the Abel/low-DCT commutator. Neither operation supplies the
missing arithmetic cancellation. The remaining task is:

```text
Prove (MATBH.29), equivalently (MATBH.34) and (MATBH.39), for every
member of one fixed cofinal sequence alpha_j->0
without assuming RH.                                  (MATBH.41)
```

A viable argument must preserve cancellation between the Type I and
Type II pieces and across dyadic scales. Rowwise absolute values,
post-collapse Vaughan, generic operator bounds, and direct promotion of
averaged Chowla remain insufficient. A signed Type I/II power gain, the
full Burnol bound, RH, PF-infinity, and `Lambda <= 0` all remain open.

Primary sources:

- Green and Tao, *Quadratic Uniformity of the Mobius Function*,
  Lemma 4.1:
  https://doi.org/10.5802/aif.2401
- Matomaki, Radziwill, and Tao, *An Averaged Form of Chowla's
  Conjecture*:
  https://doi.org/10.2140/ant.2015.9.2167
- Matomaki and Teravainen, *On the Mobius Function in All Short
  Intervals*:
  https://arxiv.org/abs/1911.09076
- Matomaki and Radziwill, *Multiplicative Functions in Short
  Intervals*:
  https://arxiv.org/abs/1501.04585
- Matomaki, Radziwill, Shao, Tao, and Teravainen, *Higher Uniformity
  of Arithmetic Functions in Short Intervals II*:
  https://doi.org/10.1007/s00222-026-01408-6
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
        "built Mertens affine-tent/bridge handoff: "
        f"{len(payload['rows'])} rows, "
        f"{payload['audit']['exact_reduction_count']} exact reductions, "
        f"{payload['audit']['proof_guard_count']} proof guards, "
        f"{payload['audit']['open_signed_gate_count']} open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
