#!/usr/bin/env python3
"""Build the Mertens shift-kernel variation and joint-cancellation handoff."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_shift_kernel_variation_handoff.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_shift_kernel_variation_handoff.md"
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
            "skv_01_inherited_reduction",
            "exact_equivalence",
            "available_exact",
            "Corollary 11.22Z.3 reduces the low half-odd criterion "
            "to the signed off-diagonal O_(alpha,K).",
            "The summable high modes and diagonal are not reopened.",
        ),
        row(
            "skv_02_shift_support",
            "exact_definition",
            "available_exact",
            "For 1<=h<=3K-2 put I_(K,h)=[K+1,4K-1-h] "
            "intersect Z.",
            "Every ordered pair K<n<m<4K occurs once as m=n+h.",
        ),
        row(
            "skv_03_shift_kernel",
            "exact_definition",
            "available_exact",
            "W_(K,h)(n):=G_(alpha,K)(n,n+h) on I_(K,h).",
            "The inherited ordinary-support kernel G is positive "
            "semidefinite and entrywise positive.",
        ),
        row(
            "skv_04_shift_decomposition",
            "exact_identity",
            "available_exact",
            "O_(alpha,K)=sum_(h=1)^(3K-2) Q_(K,h), where "
            "Q_(K,h)=sum_(n in I_(K,h))"
            "mu(n)mu(n+h)W_(K,h)(n).",
            "No absolute value is taken over h.",
        ),
        row(
            "skv_05_cc_support",
            "exact_support",
            "available_exact",
            "I_h^(cc)=[K+1,2K-1-h], nonempty exactly when "
            "1<=h<=K-2.",
            "Both entries lie in the current interval.",
        ),
        row(
            "skv_06_cf_support",
            "exact_support",
            "available_exact",
            "I_h^(cf)=[max(K+1,2K-h),"
            "min(2K-1,4K-1-h)].",
            "The interval is interpreted as empty when its lower "
            "endpoint exceeds its upper endpoint.",
        ),
        row(
            "skv_07_ff_support",
            "exact_support",
            "available_exact",
            "I_h^(ff)=[2K,4K-1-h], nonempty exactly when "
            "1<=h<=2K-1.",
            "Both entries lie in the future interval.",
        ),
        row(
            "skv_08_block_partition",
            "exact_partition",
            "available_exact",
            "I_(K,h) is the disjoint consecutive union of its "
            "nonempty cc, cf, and ff blocks.",
            "There are at most two adjacent block joins.",
        ),
        row(
            "skv_09_cc_formula",
            "exact_identity",
            "available_exact",
            "On I_h^(cc), with i=n-K, "
            "W_(K,h)(n)=P_(K,R)(i,i+h).",
            "Here 1<=i<i+h<=K-1.",
        ),
        row(
            "skv_10_cf_formula",
            "exact_identity",
            "available_exact",
            "On I_h^(cf), with i=n-K and m=n+h, "
            "W_(K,h)(n)=p_i*b_K(m), where p_i=P_(K,R)(i,K).",
            "The boundary value b_K(2K)=1 is included.",
        ),
        row(
            "skv_11_ff_formula",
            "exact_identity",
            "available_exact",
            "On I_h^(ff), W_(K,h)(n)="
            "P_(K,R)(K,K)b_K(n)b_K(n+h).",
            "Both future weights lie in (0,1].",
        ),
        row(
            "skv_12_future_monotonicity",
            "exact_monotonicity",
            "available_exact",
            "b_K(x)=(2K/x)^(1+alpha)(4K-x)/(2K) is strictly "
            "decreasing on [2K,4K).",
            "Its logarithmic derivative is "
            "-(1+alpha)/x-1/(4K-x)<0.",
        ),
        row(
            "skv_13_uniform_height",
            "exact_bound",
            "available_exact",
            "With rho_K=pi^2*K/(2K+1)^2, "
            "0<W_(K,h)(n)<=rho_K<pi^2/(4K).",
            "PSD Cauchy-Schwarz and "
            "P(i,i)<=pi^2*i/(2K+1)^2 give the bound.",
        ),
        row(
            "skv_14_cc_monotonicity",
            "exact_monotonicity",
            "available_exact",
            "For fixed h the cc weight P(i,i+h) is strictly "
            "increasing in i.",
            "Its sine-integral lower endpoint is fixed and its upper "
            "endpoint strictly increases inside (0,pi).",
        ),
        row(
            "skv_15_cc_variation",
            "exact_bound",
            "available_exact",
            "The discrete variation of W on I_h^(cc) is less than "
            "rho_K.",
            "A positive increasing sequence has variation equal to "
            "its endpoint difference.",
        ),
        row(
            "skv_16_cf_variation",
            "exact_bound",
            "available_exact",
            "The discrete variation of p_i*b_K(K+i+h) on "
            "I_h^(cf) is less than 2rho_K.",
            "The first factor increases, the second decreases, and "
            "the product variation is bounded by the sum of their "
            "two monotone variations.",
        ),
        row(
            "skv_17_ff_variation",
            "exact_bound",
            "available_exact",
            "The discrete variation of W on I_h^(ff) is less than "
            "rho_K.",
            "The product b_K(n)b_K(n+h) is strictly decreasing.",
        ),
        row(
            "skv_18_join_variation",
            "exact_bound",
            "available_exact",
            "Each nonempty cc/cf or cf/ff join contributes at most "
            "rho_K to variation.",
            "Both one-sided values lie in [0,rho_K].",
        ),
        row(
            "skv_19_fixed_shift_variation",
            "exact_bound",
            "available_exact",
            "Var_(n in I_(K,h)) W_(K,h)(n)<6rho_K.",
            "Add rows 15-18 over at most three blocks and two joins.",
        ),
        row(
            "skv_20_shift_prefix",
            "exact_definition",
            "available_exact",
            "C_(K,h)(x):=sum_(n=K+1)^x mu(n)mu(n+h), and "
            "A_(K,h):=max_(K+1<=x<=4K-1-h)|C_(K,h)(x)|.",
            "This is a fixed-origin maximal two-point correlation.",
        ),
        row(
            "skv_21_fixed_shift_abel",
            "exact_identity",
            "available_exact",
            "Q_(K,h)=C_(K,h)(b_h)W_h(b_h)+"
            "sum_(n=K+1)^(b_h-1)C_(K,h)(n)"
            "[W_h(n)-W_h(n+1)], where b_h=4K-1-h.",
            "This is finite Abel summation on one shift row.",
        ),
        row(
            "skv_22_fixed_shift_bound",
            "exact_bound",
            "available_exact",
            "|Q_(K,h)|<7rho_K*A_(K,h)<"
            "7pi^2*A_(K,h)/(4K).",
            "The terminal height costs rho_K and row 19 costs "
            "6rho_K.",
        ),
        row(
            "skv_23_fixed_base_path",
            "exact_identity",
            "available_exact",
            "For n=K+i<2K, h first traces P(i,i+h) through "
            "the current block up to P(i,K), then traces "
            "p_i*b_K(n+h) down the future block.",
            "The value at n+h=2K agrees because b_K(2K)=1.",
        ),
        row(
            "skv_24_kernel_row_variation",
            "exact_bound",
            "available_exact",
            "For fixed i, the variation of P_(K,R)(i,j) over "
            "i<j<=K is less than pi^2/(2K+1).",
            "Successive sine-integral boundary strips are disjoint; "
            "their total is at most "
            "4(2K+1)^(-1) integral_0^pi S_R<pi^2/(2K+1).",
        ),
        row(
            "skv_25_current_future_tail",
            "exact_bound",
            "available_exact",
            "For fixed current n, the variation after n+h reaches "
            "2K is at most p_i<=rho_K.",
            "The future factor b_K(n+h) decreases from one.",
        ),
        row(
            "skv_26_future_base_monotonicity",
            "exact_monotonicity",
            "available_exact",
            "For fixed future n, W_(K,h)(n) is decreasing in h.",
            "Only the factor b_K(n+h) varies.",
        ),
        row(
            "skv_27_fixed_base_variation",
            "exact_bound",
            "available_exact",
            "For every K<n<4K-1, "
            "Var_h W_(K,h)(n)<pi^2/(2K+1)+rho_K.",
            "Rows 23-26 cover current and future base points.",
        ),
        row(
            "skv_28_base_prefix",
            "exact_definition",
            "available_exact",
            "B_(K,n)(H):=mu(n)sum_(h=1)^H mu(n+h), "
            "1<=H<=4K-1-n.",
            "This is one local Mertens interval multiplied by mu(n).",
        ),
        row(
            "skv_29_fixed_base_abel",
            "exact_identity",
            "available_exact",
            "The fixed-n row sum has the exact Abel expansion in the "
            "prefixes B_(K,n)(H) and successive h-differences of W.",
            "No absolute value is inserted before the identity.",
        ),
        row(
            "skv_30_fixed_base_bound",
            "exact_bound",
            "available_exact",
            "The absolute fixed-n row is less than "
            "pi^2/K times max_H|sum_(h<=H)mu(n+h)|.",
            "The terminal height plus row 27 is "
            "pi^2/(2K+1)+2rho_K<pi^2/K.",
        ),
        row(
            "skv_31_future_l1",
            "exact_bound",
            "available_exact",
            "B_1:=sum_(2K<=n<4K)b_K(n)<=(2K+1)/2.",
            "Use b_K(4K-ell)<=ell/(2K) and sum 1<=ell<=2K.",
        ),
        row(
            "skv_32_future_l2",
            "exact_bound",
            "available_exact",
            "B_2:=sum_(2K<=n<4K)b_K(n)^2<="
            "(2K+1)(4K+1)/(12K).",
            "Use the same pointwise envelope and sum the squares.",
        ),
        row(
            "skv_33_total_cc_mass",
            "exact_bound",
            "available_exact",
            "The total cc kernel mass is at most "
            "rho_K*binom(K-1,2)<pi^2*K/8.",
            "There are binom(K-1,2) current-current pairs.",
        ),
        row(
            "skv_34_total_cf_mass",
            "exact_bound",
            "available_exact",
            "The total cf kernel mass is at most "
            "rho_K*(K-1)B_1<pi^2*K/4.",
            "Use p_i<=rho_K and row 31.",
        ),
        row(
            "skv_35_total_ff_mass",
            "exact_bound",
            "available_exact",
            "The total ff kernel mass is at most "
            "rho_K*B_1^2/2<=pi^2*K/8.",
            "Use sum_(n<m)b_K(n)b_K(m)<=B_1^2/2.",
        ),
        row(
            "skv_36_total_mass",
            "exact_bound",
            "available_exact",
            "sum_(K<n<m<4K)G_(alpha,K)(n,m)<pi^2*K/2.",
            "Add rows 33-35.",
        ),
        row(
            "skv_37_per_shift_mass",
            "exact_bound",
            "available_exact",
            "For every h, sum_(n in I_(K,h))W_(K,h)(n)"
            "<13pi^2/24.",
            "The cc, cf, and ff blocks cost respectively less than "
            "pi^2/8, pi^2/4, and pi^2/6.",
        ),
        row(
            "skv_38_trivial_barrier",
            "exact_bound",
            "available_exact",
            "|O_(alpha,K)|<pi^2*K/2, so the normalized scale is "
            "K^(1-alpha).",
            "This recovers the one-power absolute barrier and is not "
            "summable for 0<alpha<1.",
        ),
        row(
            "skv_39_axiswise_square_root_guard",
            "proof_guard",
            "guard_active",
            "Even uniform square-root-plus-epsilon bounds for every "
            "fixed-shift prefix or every fixed-base Mertens interval "
            "give only O(K^(1/2+epsilon)) after summing the other "
            "axis absolutely.",
            "That closes only alpha>1/2+epsilon and cannot reach a "
            "cofinal alpha sequence tending to zero.",
        ),
        row(
            "skv_40_averaged_chowla_guard",
            "literature_guard",
            "guard_active",
            "The natural averaged-Chowla terminal-correlation scale "
            "o(K^2), even if promoted to the maximal prefixes in row "
            "20, gives only o(K) through row 22.",
            "The cited theorem does not supply this maximal promotion, "
            "and o(K) still has no fixed power saving sufficient for "
            "the cofinal criterion.",
        ),
        row(
            "skv_41_absolute_shift_sufficient",
            "conditional_calibration",
            "conditional_only",
            "The strong variational estimate "
            "sum_(h<=3K-2)A_(K,h)=O_epsilon(K^(1+epsilon)) "
            "would imply O_(alpha,K)=O_epsilon(K^epsilon).",
            "This absolute-shift condition is far stronger than "
            "random square-root size and is not claimed.",
        ),
        row(
            "skv_42_joint_arithmetic_gate",
            "open_gate",
            "open",
            "Prove the weighted cumulative signed sum "
            "sup_J sum_(K<=2^J)K^(-alpha)sum_h Q_(K,h)<infinity "
            "for each member of one cofinal positive-alpha sequence.",
            "Cancellation must remain joint in n, h, Vaughan type, "
            "and possibly K.",
        ),
        row(
            "skv_43_joint_vaughan_identity",
            "exact_identity",
            "available_exact",
            "The same sum is exactly "
            "<u_I,G u_I>+<u_II,G u_II>-2<u_I,G u_II>.",
            "Rows 21 and 29 are theorem interfaces, not permission to "
            "separate the three Vaughan terms.",
        ),
        row(
            "skv_44_finite_validation",
            "finite_validation",
            "validated_finite",
            "Independent finite checks verify the shift partition, "
            "both variation bounds, both Abel identities, mass bounds, "
            "and direct/shift quadratic agreement.",
            "Finite checks audit the algebra and do not prove the open "
            "all-scale Mobius estimate.",
        ),
        row(
            "skv_45_proof_boundary",
            "proof_guard",
            "guard_active",
            "The shift representation and O(1/K) bounded variation "
            "do not prove cancellation across shifts.",
            "No weighted joint Mobius gain, full Burnol bound, RH, "
            "PF-infinity, or Lambda<=0 conclusion is supplied.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_mertens_shift_kernel_variation_handoff",
        "date": "2026-07-24",
        "status": (
            "exact shift decomposition, two-axis bounded variation, "
            "and joint-cancellation handoff with one open arithmetic gate"
        ),
        "parameters": {
            "alpha_range": "0<alpha<1",
            "dyadic_blocks": "K=2^j",
            "cutoff": "R_(alpha,K)=ceil(K^(1-alpha/2))",
            "ambient_size": "N=2K+1",
            "kernel_height": "rho_K=pi^2*K/N^2<pi^2/(4K)",
            "fixed_shift_variation": "less than 6rho_K",
            "fixed_base_variation": "less than pi^2/N+rho_K",
            "total_offdiagonal_mass": "less than pi^2*K/2",
        },
        "source_anchors": [
            "outputs/jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.md",
            "outputs/jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.md",
            "https://doi.org/10.2140/ant.2015.9.2167",
            "https://doi.org/10.1017/fmp.2023.28",
            "https://doi.org/10.5802/aif.2401",
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
            "shift_decomposition_proved": True,
            "fixed_shift_variation_proved": True,
            "fixed_base_variation_proved": True,
            "mass_bounds_proved": True,
            "axiswise_routes_rejected": True,
            "joint_mobius_gain_proved": False,
            "full_burnol_bound_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
    }


def render_note() -> str:
    return r"""# Jensen-Window PF Mertens Shift-Kernel Variation Handoff

Date: 2026-07-24

Status: exact shift decomposition, two-axis bounded variation, and
joint-cancellation handoff with one open arithmetic gate.
This is not a proof of RH, PF-infinity, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_mertens_shift_kernel_variation_handoff.json
python work/rh_compute/scripts/jensen_window_pf_mertens_shift_kernel_variation_handoff.py
python work/rh_compute/scripts/check_jensen_window_pf_mertens_shift_kernel_variation_handoff.py
```

## Exact Shift Decomposition

Retain the notation of Corollary 11.22Z.3:

```text
N=2K+1,
R=ceil(K^(1-alpha/2)),
P=P_(K,R),
G=G_(alpha,K,R),
rho_K:=pi^2*K/N^2<pi^2/(4K).               (SKV.1)
```

For `1<=h<=3K-2`, put

```text
I_(K,h):=[K+1,4K-1-h] intersect Z,
W_(K,h)(n):=G(n,n+h),

Q_(K,h):=sum_(n in I_(K,h))
 mu(n)mu(n+h)W_(K,h)(n).                    (SKV.2)
```

Every pair `K<n<m<4K` occurs exactly once with `h=m-n`, so

```text
O_(alpha,K)=sum_(h=1)^(3K-2)Q_(K,h).        (SKV.3)
```

This identity preserves the sign across shifts.

## Three Consecutive Blocks

The shift support is the disjoint consecutive union

```text
I_h^(cc)=[K+1,2K-1-h],

I_h^(cf)=[
 max(K+1,2K-h),
 min(2K-1,4K-1-h)
],

I_h^(ff)=[2K,4K-1-h],                       (SKV.4)
```

with an interval omitted when empty. With `i=n-K`, `m=n+h`, and
`p_i=P(i,K)`,

```text
W_h(n)=
  P(i,i+h),                    n in I_h^(cc),
  p_i*b_K(m),                  n in I_h^(cf),
  P(K,K)b_K(n)b_K(n+h),        n in I_h^(ff). (SKV.5)
```

The future weight

```text
b_K(x)=(2K/x)^(1+alpha)*(4K-x)/(2K)
```

is strictly decreasing on `[2K,4K)`. PSD Cauchy-Schwarz and the
diagonal estimate from Corollary 11.22Z.3 give

```text
0<W_h(n)<=rho_K<pi^2/(4K).                   (SKV.6)
```

## Fixed-Shift Variation

On the current-current block, the sine-integral formula becomes

```text
P(i,i+h)
 =2/N integral_(h*pi/N)^((2i+h)*pi/N)S_R(u)du. (SKV.7)
```

Since `S_R>0`, this is strictly increasing in `i`, and its variation
is less than `rho_K`. On the current-future block, `p_i` increases
while `b_K(K+i+h)` decreases. For positive monotone factors,

```text
Var(p_i*b_K(K+i+h))<2rho_K.                  (SKV.8)
```

The future-future block is decreasing and has variation less than
`rho_K`. Each of the at most two block joins costs at most `rho_K`.
Consequently,

```text
Var_(n in I_(K,h))W_(K,h)(n)<6rho_K.         (SKV.9)
```

Define the fixed-origin maximal correlation

```text
C_(K,h)(x):=sum_(n=K+1)^x mu(n)mu(n+h),

A_(K,h):=max_(K+1<=x<=4K-1-h)|C_(K,h)(x)|.  (SKV.10)
```

Finite Abel summation gives, with `b_h=4K-1-h`,

```text
Q_(K,h)
 =C_(K,h)(b_h)W_h(b_h)
  +sum_(n=K+1)^(b_h-1)
    C_(K,h)(n)[W_h(n)-W_h(n+1)],             (SKV.11)

|Q_(K,h)|
 <7rho_K*A_(K,h)
 <7pi^2/(4K)*A_(K,h).                        (SKV.12)
```

Thus each shift weight is a clean bounded-variation test of a local
two-point Mobius correlation.

## Fixed-Base Variation

Fix a current base point `n=K+i`. As `h` grows, the second position
first traces

```text
P(i,i+1),...,P(i,K),
```

then traces `p_i*b_K(n+h)` down the future interval. Successive
sine-integral boundary strips are disjoint, and

```text
integral_0^pi S_R(u)du
 =2sum_(r<=R)1/(2r-1)^2<pi^2/4.
```

Therefore

```text
Var_(i<j<=K)P(i,j)<pi^2/N.                  (SKV.13)
```

The future tail costs at most `p_i<=rho_K`. A future base point has a
decreasing `h`-row. Uniformly,

```text
Var_h W_(K,h)(n)<pi^2/N+rho_K.               (SKV.14)
```

Put

```text
B_(K,n)(H):=mu(n)sum_(h=1)^H mu(n+h).
```

A second Abel summation gives the exact fixed-base expansion and the
bound

```text
|sum_(h=1)^(4K-1-n)
  mu(n)mu(n+h)W_(K,h)(n)|

 <pi^2/K*
  max_(H<=4K-1-n)|sum_(h<=H)mu(n+h)|.        (SKV.15)
```

Thus the other coordinate is an `O(1/K)` bounded-variation test of
local Mertens intervals.

## Mass Calibration

The future weights satisfy

```text
B_1:=sum_(2K<=n<4K)b_K(n)<=(2K+1)/2,

B_2:=sum_(2K<=n<4K)b_K(n)^2
 <=(2K+1)(4K+1)/(12K).                      (SKV.16)
```

The total current-current, current-future, and future-future masses
are respectively less than

```text
pi^2*K/8,  pi^2*K/4,  pi^2*K/8.
```

Hence

```text
sum_(K<n<m<4K)G(n,m)<pi^2*K/2,

|O_(alpha,K)|<pi^2*K/2.                     (SKV.17)
```

For a single shift, PSD Cauchy-Schwarz, (SKV.16), and the trace bound
give

```text
sum_(n in I_(K,h))W_(K,h)(n)<13pi^2/24.     (SKV.18)
```

The absolute normalized scale remains `K^(1-alpha)`, one full power
above the summable target.

## Why Axiswise Cancellation Is Insufficient

Suppose, conditionally, every maximal fixed-shift correlation in
(SKV.10), or every local Mertens maximum in (SKV.15), had the
square-root-plus-epsilon size `O(K^(1/2+epsilon))`. Summing the other
coordinate absolutely still gives

```text
O_(alpha,K)=O(K^(1/2+epsilon)).              (SKV.19)
```

After multiplication by `K^(-alpha)`, this closes only
`alpha>1/2+epsilon`, not a cofinal sequence tending to zero.

Matomaki, Radziwill, and Tao prove an averaged Chowla theorem whose
natural terminal two-point scale is `o(K^2)` when both averaging
ranges have size `K`. Their theorem does not provide the varying
maximal endpoints in (SKV.10). Even granting that promotion, (SKV.12)
would give only `o(K)`, with no fixed power saving sufficient for the
cofinal criterion.

Similarly, scalar all-interval logarithmic cancellation applied one
base point at a time yields at best `K*log^(-A)K` after the `1/K`
normalization. Logarithms do not replace the missing power.

The stronger condition

```text
sum_(h<=3K-2)A_(K,h)
 =O_epsilon(K^(1+epsilon))                   (SKV.20)
```

would suffice, but it is stronger than random square-root behavior
after absolute summation and is not proposed as the likely route.

## Joint Gate

The exact remaining criterion is

```text
sup_J sum_(K dyadic<=2^J)
 K^(-alpha)sum_(h=1)^(3K-2)Q_(K,h)
 <infinity                                   (SKV.21)
```

for every member of one fixed cofinal sequence `alpha_j->0`.
Equivalently, before any axiswise absolute value, the scale is

```text
<u_I,G u_I>+<u_II,G u_II>-2<u_I,G u_II>.    (SKV.22)
```

The new bounded-variation formulas are theorem interfaces for a
genuinely joint estimate in base point, shift, Vaughan type, and
possibly dyadic scale. They do not themselves provide that estimate.

The weighted joint Mobius gain, full Burnol bound, RH, PF-infinity,
and `Lambda<=0` remain open.

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
        "built Mertens shift-kernel variation handoff: "
        "45 rows, 39 exact reductions, 2 proof guards, 1 open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
