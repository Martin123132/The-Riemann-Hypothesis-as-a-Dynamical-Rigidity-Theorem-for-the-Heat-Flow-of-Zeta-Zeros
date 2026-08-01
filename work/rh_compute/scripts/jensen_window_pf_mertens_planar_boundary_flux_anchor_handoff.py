#!/usr/bin/env python3
"""Build the joined planar boundary-flux and anchor handoff."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.md"
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


def build_rows() -> list[dict[str, str]]:
    return [
        row(
            "pbfah_01_parameters",
            "exact_definition",
            "available_exact",
            "N=2K+1, q_r=2r-1, theta_r=q_r*pi/N, and "
            "phi_r(i)=2*N^(-1/2)sin(i*theta_r).",
            "The retained modes are 1<=r<=R; no asymptotic estimate is "
            "asserted.",
        ),
        row(
            "pbfah_02_ordinary_feature",
            "exact_definition",
            "available_exact",
            "With e_r=phi_r(K), f_r(K+i)=phi_r(i) for 1<=i<K, "
            "f_r(n)=b_n e_r for 2K<=n<4K, and f_r=0 outside the "
            "support, G(n,m)=sum_r f_r(n)f_r(m)/q_r^2.",
            "This is the truncated ordinary-Mobius pullback from the "
            "parent handoff.",
        ),
        row(
            "pbfah_03_global_flux",
            "exact_identity",
            "available_exact",
            "For m=n+h and delta f_r(t)=f_r(t)-f_r(t+1), "
            "Delta W(n,h)=sum_r[f_r(n)delta f_r(m)-"
            "f_r(n+1)delta f_r(m+1)]/q_r^2.",
            "The identity is finite and modewise; it supplies no Mobius "
            "cancellation.",
        ),
        row(
            "pbfah_04_current_interior",
            "exact_identity",
            "available_exact",
            "For n=K+i and m=K+j<=2K-2, the global flux specializes "
            "to the signed current-current curvature used in "
            "Corollary 11.22Z.7.",
            "This recovers the anti-diagonal component but does not "
            "separate it safely from its terminal boundary.",
        ),
        row(
            "pbfah_05_second_transition",
            "exact_identity",
            "available_exact",
            "At m=2K-1 the mode numerator is "
            "phi_r(i)[phi_r(K-1)-e_r]-phi_r(i+1)e_r d_(2K).",
            "Summing with q_r^(-2) reproduces "
            "P(i,K-1)-p_i-p_(i+1)d_(2K).",
        ),
        row(
            "pbfah_06_current_future",
            "exact_identity",
            "available_exact",
            "For current n=K+i and future m, the mode numerator is "
            "e_r[(phi_r(i)-phi_r(i+1))d_m+"
            "phi_r(i+1)(d_m-d_(m+1))].",
            "Both terms must remain signed before summing cells or modes.",
        ),
        row(
            "pbfah_07_first_transition",
            "exact_identity",
            "available_exact",
            "At n=2K-1 and future m, the mode numerator is "
            "e_r[phi_r(K-1)-e_r]d_m+"
            "e_r^2(d_m-d_(m+1)).",
            "Summing with q_r^(-2) reproduces "
            "p_(K-1)d_m-c d_(m+1).",
        ),
        row(
            "pbfah_08_future_future",
            "exact_identity",
            "available_exact",
            "For 2K<=n<m, the mode numerator is "
            "e_r^2[b_n(d_m-d_(m+1))+d_n d_(m+1)].",
            "Its geometric sign is nonnegative, but its Mobius prefix "
            "pairing is not thereby estimated.",
        ),
        row(
            "pbfah_09_endpoint_symmetry",
            "exact_identity",
            "available_exact",
            "For odd q_r, phi_r(K+1)=phi_r(K)=e_r.",
            "This is a finite trigonometric endpoint identity.",
        ),
        row(
            "pbfah_10_endpoint_recurrence",
            "exact_identity",
            "available_exact",
            "The sine recurrence gives "
            "phi_r(K-1)=(2cos(theta_r)-1)e_r.",
            "No small-angle approximation is used.",
        ),
        row(
            "pbfah_11_smoothed_endpoint_jump",
            "exact_identity",
            "available_exact",
            "phi_r(K-1)-e_r=-4sin^2(theta_r/2)e_r="
            "-(pi^2 q_r^2/N^2)sinc^2(theta_r/2)e_r.",
            "The q_r^2 factor cancels the kernel denominator only in "
            "this endpoint jump.",
        ),
        row(
            "pbfah_12_second_transition_smoothing",
            "exact_identity",
            "available_exact",
            "The first part of the m=2K-1 transition equals "
            "-(pi^2/N^2)sum_r phi_r(i)e_r "
            "sinc^2(theta_r/2); the remaining term carries d_(2K).",
            "This smoothing does not control the signed prefix paired "
            "with the complete transition.",
        ),
        row(
            "pbfah_13_first_transition_smoothing",
            "exact_identity",
            "available_exact",
            "The first part of the n=2K-1 transition equals "
            "-(pi^2/N^2)sum_r e_r^2 d_m "
            "sinc^2(theta_r/2); the other part carries "
            "d_m-d_(m+1).",
            "The formula is exact but is not an arithmetic estimate.",
        ),
        row(
            "pbfah_14_lower_interface_join",
            "exact_identity",
            "available_exact",
            "At m=2K-1 the zero-extended CC boundary "
            "phi_r(i)phi_r(K-1) joins the CF boundary "
            "-phi_r(i)e_r-phi_r(i+1)e_r d_(2K) exactly.",
            "Estimating the two artificial block boundaries separately "
            "can create a loss absent from the actual kernel.",
        ),
        row(
            "pbfah_15_left_interface_join",
            "exact_identity",
            "available_exact",
            "At n=2K-1 the zero-extended CF boundary "
            "phi_r(K-1)e_r d_m joins the FF boundary "
            "-e_r^2 d_(m+1) exactly.",
            "The join removes an artificial interface split, not the "
            "genuine current boundary contribution.",
        ),
        row(
            "pbfah_16_block_split",
            "exact_definition",
            "available_exact",
            "Write u_i=mu(K+i), v_m=mu(m), and split O into "
            "O_CC+O_CF+O_FF on the current and future supports.",
            "The split is algebraic; estimates must retain any stated "
            "cross-block cancellation.",
        ),
        row(
            "pbfah_17_current_mode_data",
            "exact_definition",
            "available_exact",
            "Set x_r=sum_(i<K)u_i phi_r(i), "
            "D_(C,r)=sum_(i<K)u_i^2 phi_r(i)^2, "
            "gamma=sum_i u_i p_i, and c=sum_r e_r^2/q_r^2.",
            "The scalar gamma is the weighted endpoint projection of "
            "the current mode vector.",
        ),
        row(
            "pbfah_18_future_scalars",
            "exact_definition",
            "available_exact",
            "Set beta=sum_(2K<=m<4K)v_m b_m and "
            "Q=sum_(2K<=m<4K)v_m^2 b_m^2.",
            "Square-root cancellation for beta is not assumed.",
        ),
        row(
            "pbfah_19_block_compression",
            "exact_identity",
            "available_exact",
            "O_CC=(1/2)sum_r(x_r^2-D_(C,r))/q_r^2, "
            "O_CF=beta*gamma, and O_FF=(c/2)(beta^2-Q).",
            "The future interval is rank one only after its weighted "
            "Mobius sum beta is retained.",
        ),
        row(
            "pbfah_20_boundary_prefix",
            "exact_definition",
            "available_exact",
            "Put S_i^bd=S_K(K+i,K-1-i), "
            "S_i^sh=S_K(K+i,K-2-i), "
            "U_r=sum_i S_i^bd phi_r(i), and "
            "V_r=sum_i S_i^sh phi_r(i+1).",
            "The shoulder is forced because the mixed stencil reaches "
            "m+2=2K on the top current-interior anti-diagonal.",
        ),
        row(
            "pbfah_21_boundary_spectrum",
            "exact_identity",
            "available_exact",
            "B_CC=sum_i S_i^bd P(i,K-1)="
            "sum_r phi_r(K-1)U_r/q_r^2, while the embedded CF shoulder "
            "A_CF=sum_i S_i^sh p_(i+1)=sum_r e_r V_r/q_r^2.",
            "B_CC and A_CF must be joined before absolute values.",
        ),
        row(
            "pbfah_22_current_reconstruction",
            "exact_identity",
            "available_exact",
            "O_CC=O_(CC,int)+B_CC-A_CF.",
            "Omitting A_CF gives the zero-extended CC interior, not the "
            "actual joined curvature used in Corollary 11.22Z.7.",
        ),
        row(
            "pbfah_23_joined_direct_decomposition",
            "exact_identity",
            "available_exact",
            "O=O_(CC,int)+B_CC-A_CF+beta*gamma+"
            "(c/2)(beta^2-Q).",
            "This is the lossless signed join of every planar cell.",
        ),
        row(
            "pbfah_24_remainder_modes",
            "exact_definition",
            "available_exact",
            "Define Z_r=phi_r(K-1)U_r-e_r V_r+"
            "beta e_r x_r+(e_r^2/2)(beta^2-Q).",
            "The shoulder, boundary, and anchor terms must be summed in "
            "Z_r before any modewise absolute value.",
        ),
        row(
            "pbfah_25_remainder_projection",
            "exact_identity",
            "available_exact",
            "B_CC-A_CF+beta*gamma+(c/2)(beta^2-Q)="
            "sum_r Z_r/q_r^2.",
            "This closes the formerly unnamed transition/future remainder "
            "as an exact projection, not as an estimate.",
        ),
        row(
            "pbfah_26_remainder_cauchy",
            "exact_bound",
            "available_exact",
            "|sum_r Z_r/q_r^2|^2<=(sum_r q_r^(-2))"
            "(sum_r |Z_r|^2/q_r^2)<(pi^2/8)"
            "sum_r |Z_r|^2/q_r^2.",
            "The joined Z_r retains the cancellation between B_CC and "
            "the embedded shoulder A_CF.",
        ),
        row(
            "pbfah_27_remainder_square_diagnostic",
            "conditional_calibration",
            "guard_active",
            "sum_r|Z_r|^2/q_r^2=O_epsilon(K^epsilon) is sufficient "
            "for a subpower remainder.",
            "It is a secondary sufficient target; the parent joined "
            "low-energy criterion remains the lossless theorem gate.",
        ),
        row(
            "pbfah_28_current_energy",
            "exact_definition",
            "available_exact",
            "Set E_C=sum_r x_r^2/q_r^2 and "
            "E_perp=E_C-gamma^2/c.",
            "Here c>0; E_perp is the squared distance from the endpoint "
            "mode direction.",
        ),
        row(
            "pbfah_29_orthogonal_anchor_completion",
            "exact_identity",
            "available_exact",
            "The complete low-mode energy is "
            "L=sum_r(x_r+beta e_r)^2/q_r^2="
            "E_perp+c(beta+gamma/c)^2.",
            "Both terms are nonnegative, so cancellation between them "
            "cannot prove the positive energy criterion.",
        ),
        row(
            "pbfah_30_anchor_direction",
            "exact_identity",
            "available_exact",
            "The future rank-one anchor can shift only the endpoint "
            "direction gamma/c; it leaves E_perp unchanged.",
            "It does not algebraically cancel the joined boundary flux "
            "B_CC-A_CF.",
        ),
        row(
            "pbfah_31_joined_energy_gate",
            "theorem_target",
            "open",
            "For every member of one fixed cofinal positive-alpha "
            "sequence, prove dyadic summability of "
            "K^(-alpha)[E_perp+c(beta+gamma/c)^2].",
            "This is the earlier full low-mode energy gate in joined "
            "coordinates; no estimate is supplied here.",
        ),
        row(
            "pbfah_32_boundary_cancellation_witness",
            "countermodel_gate",
            "guard_validated",
            "For K=9, P_infinity(i,j)=pi^2 min(i,j)/19^2, and "
            "u=(-1,-1,-1,-1,-1,-1,1,1), one has "
            "O_CC=O_(CC,int)=0 and "
            "B_CC=A_CF=198pi^2/19^2.",
            "This bounded synthetic sequence proves that omitting the "
            "shoulder manufactures a nonzero boundary remainder.",
        ),
        row(
            "pbfah_33_shoulder_join_guard",
            "proof_guard",
            "guard_active",
            "The Corollary 11.22Z.7 Y_(K,r) square-function remains a "
            "valid and potentially useful condition for the actual "
            "joined interior; the remaining theorem must use B_CC-A_CF.",
            "No estimate may replace the shoulder by zero or bound "
            "B_CC and A_CF separately.",
        ),
        row(
            "pbfah_34_mobius_finite_audit",
            "finite_validation",
            "validated_finite",
            "For alpha=1/2 and K=16,...,1024, direct Mobius audits "
            "reproduce O_CC=O_(CC,int)+B_CC-A_CF and the full block "
            "compression; at K=1024 O_(CC,int)=2.086e-5 while "
            "B_CC=-184.3328081 and A_CF=-184.1213147.",
            "The table is finite diagnostic evidence only and proves no "
            "asymptotic cancellation.",
        ),
        row(
            "pbfah_35_proof_boundary",
            "proof_guard",
            "guard_active",
            "The flux, joins, block compression, boundary spectrum, "
            "orthogonal completion, and countermodel are exact.",
            "No Mobius energy estimate, reciprocal-tail theorem, RH, "
            "PF-infinity, Lambda<=0, or Clay-prize result is proved.",
        ),
    ]


MOBIUS_AUDIT_ROWS = [
    {
        "K": 16,
        "R": 8,
        "O_CC": 0.0112021851505,
        "B_CC": -1.64707693737,
        "A_CF": -1.65999311032,
        "joined_boundary": 0.0129161729502,
        "O_CC_int": -0.00171398779971,
        "anchor_correction": -0.606047689794,
        "O_total": -0.594845504643,
    },
    {
        "K": 32,
        "R": 14,
        "O_CC": -0.258233497456,
        "B_CC": -3.3956430504,
        "A_CF": -3.13823518152,
        "joined_boundary": -0.257407868874,
        "O_CC_int": -0.00082562858183,
        "anchor_correction": -0.338634157266,
        "O_total": -0.596867654722,
    },
    {
        "K": 64,
        "R": 23,
        "O_CC": -0.285714586905,
        "B_CC": -9.87073687071,
        "A_CF": -9.5849702148,
        "joined_boundary": -0.285766655914,
        "O_CC_int": 5.20690090173e-05,
        "anchor_correction": -0.289114323711,
        "O_total": -0.574828910616,
    },
    {
        "K": 128,
        "R": 39,
        "O_CC": -0.308845847746,
        "B_CC": -20.950184064,
        "A_CF": -20.6412849618,
        "joined_boundary": -0.30889910222,
        "O_CC_int": 5.32544736735e-05,
        "anchor_correction": -0.277062237592,
        "O_total": -0.585908085338,
    },
    {
        "K": 256,
        "R": 64,
        "O_CC": -0.337635715345,
        "B_CC": -47.8276172227,
        "A_CF": -47.4901484097,
        "joined_boundary": -0.337468812975,
        "O_CC_int": -0.000166902369834,
        "anchor_correction": -0.268082212656,
        "O_total": -0.605717928,
    },
    {
        "K": 512,
        "R": 108,
        "O_CC": -0.311239977615,
        "B_CC": -88.8508867578,
        "A_CF": -88.5396892323,
        "joined_boundary": -0.311197525512,
        "O_CC_int": -4.24521027469e-05,
        "anchor_correction": -0.295508209027,
        "O_total": -0.606748186642,
    },
    {
        "K": 1024,
        "R": 182,
        "O_CC": -0.211472500994,
        "B_CC": -184.332808071,
        "A_CF": -184.12131471,
        "joined_boundary": -0.21149336131,
        "O_CC_int": 2.08603166243e-05,
        "anchor_correction": -0.402560286509,
        "O_total": -0.614032787502,
    },
]


NOTE_BODY = r"""
# Jensen-Window PF Mertens Planar Boundary-Flux/Anchor Handoff

Date: 2026-07-24

Status: exact joined boundary-flux and rank-one anchor reduction with
one open Mobius energy gate. This is not a proof of RH, PF-infinity,
`Lambda<=0`, or a Clay-prize result.

```text
work/rh_compute/results/
jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.json
python work/rh_compute/scripts/
jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.py
python work/rh_compute/scripts/
check_jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.py
```

## One Modewise Flux Law

Retain

```text
N=2K+1, q_r=2r-1, theta_r=q_r*pi/N,
phi_r(i)=2/sqrt(N) sin(i theta_r), e_r=phi_r(K).
```

Extend the ordinary feature by zero outside `K<n<4K`:

```text
f_r(K+i)=phi_r(i),                    1<=i<K,
f_r(n)=b_n e_r,                       2K<=n<4K,
delta f_r(t)=f_r(t)-f_r(t+1).
```

Then

```text
G(n,m)=sum_(r<=R) f_r(n)f_r(m)/q_r^2.
```

For `m=n+h`, direct expansion of the four terms in the mixed
difference gives the global finite identity

```text
Delta W(n,h)
 =sum_(r<=R) 1/q_r^2 [
    f_r(n)delta f_r(m)
    -f_r(n+1)delta f_r(m+1)
  ].                                                   (PBFA.1)
```

This single law reproduces all five regions from Corollary 11.22Z.5.
For current `n=K+i`, its four noninterior specializations have mode
numerators

```text
m=2K-1:
 phi_r(i)[phi_r(K-1)-e_r]
 -phi_r(i+1)e_r d_(2K),                               (PBFA.2)

m>=2K:
 e_r[(phi_r(i)-phi_r(i+1))d_m
     +phi_r(i+1)(d_m-d_(m+1))],                       (PBFA.3)

n=2K-1, m>=2K:
 e_r[phi_r(K-1)-e_r]d_m
 +e_r^2(d_m-d_(m+1)),                                 (PBFA.4)

2K<=n<m:
 e_r^2[b_n(d_m-d_(m+1))+d_n d_(m+1)].                (PBFA.5)
```

Summing (PBFA.2)-(PBFA.5) with `q_r^(-2)` gives exactly the two
transition, current-future, and future-future formulas already proved
in Corollary 11.22Z.5.

## Endpoint Smoothing And Exact Joins

Oddness of `q_r` and the sine recurrence give

```text
phi_r(K+1)=phi_r(K)=e_r,
phi_r(K-1)=(2cos(theta_r)-1)e_r,

phi_r(K-1)-e_r
 =-4sin^2(theta_r/2)e_r
 =-(pi^2 q_r^2/N^2)sinc^2(theta_r/2)e_r.              (PBFA.6)
```

Thus the endpoint jump in each transition cancels the `q_r^2`
kernel denominator exactly. The other transition terms carry either
`d_(2K)`, `d_m-d_(m+1)`, or both. This is genuine geometric smoothing,
but it does not estimate the signed Mobius prefixes paired with those
cells.

The zero-extended block boundaries explain why the join must precede
absolute values. At `m=2K-1`, the CC and CF mode numerators are

```text
phi_r(i)phi_r(K-1),
-phi_r(i)e_r-phi_r(i+1)e_r d_(2K).
```

Their sum is (PBFA.2). At `n=2K-1`, the CF and FF numerators are

```text
phi_r(K-1)e_r d_m,
-e_r^2 d_(m+1),
```

whose sum is (PBFA.4). The joins remove artificial interface jumps.
They do not erase the genuine terminal current boundary.

## Lossless Block Compression

Write

```text
u_i=mu(K+i),                           1<=i<K,
v_m=mu(m),                             2K<=m<4K,

x_r=sum_(i=1)^(K-1)u_i phi_r(i),
D_(C,r)=sum_(i=1)^(K-1)u_i^2 phi_r(i)^2,

beta=sum_(m=2K)^(4K-1)v_m b_m,
Q=sum_(m=2K)^(4K-1)v_m^2 b_m^2,

gamma=sum_(i=1)^(K-1)u_i p_i
     =sum_(r<=R)e_r x_r/q_r^2,
c=P(K,K)=sum_(r<=R)e_r^2/q_r^2.
```

The three direct blocks are exactly

```text
O_CC=(1/2)sum_(r<=R)(x_r^2-D_(C,r))/q_r^2,
O_CF=beta gamma,
O_FF=(c/2)(beta^2-Q).                              (PBFA.7)
```

For the terminal current prefixes and the top interior shoulder set

```text
S_i^bd:=S_K(K+i,K-1-i),                  1<=i<=K-2,
S_i^sh:=S_K(K+i,K-2-i),                  1<=i<=K-3,
U_r:=sum_(i=1)^(K-2)S_i^bd phi_r(i),
V_r:=sum_(i=1)^(K-3)S_i^sh phi_r(i+1),

B_CC:=sum_(i=1)^(K-2)S_i^bd P(i,K-1)
    =sum_(r<=R)phi_r(K-1)U_r/q_r^2,

A_CF:=sum_(i=1)^(K-3)S_i^sh p_(i+1)
    =sum_(r<=R)e_r V_r/q_r^2.                       (PBFA.8)
```

The shoulder is not optional: on the top cell `m=2K-2`, the fourth
point of the mixed stencil is `m+2=2K`, where the joined kernel equals
the current-future anchor because `b_(2K)=1`. Finite double Abel
summation on the current triangle therefore gives

```text
O_CC=O_(CC,int)+B_CC-A_CF.                         (PBFA.9)
```

Consequently every transition and future cell is closed into the
lossless direct formula

```text
O
 =O_(CC,int)+B_CC-A_CF+beta gamma+(c/2)(beta^2-Q).
                                                               (PBFA.10)
```

Equivalently, with

```text
Z_r
 :=phi_r(K-1)U_r-e_r V_r+beta e_r x_r
   +(e_r^2/2)(beta^2-Q),
```

one has

```text
O=O_(CC,int)+sum_(r<=R)Z_r/q_r^2.                 (PBFA.11)
```

This is the compatible signed endpoint-spectral coordinate requested
by Corollary 11.22Z.7. It retains cancellation among all boundary and
future terms. Cauchy gives the valid sufficient estimate

```text
|sum_r Z_r/q_r^2|^2
 <=(sum_r q_r^(-2))(sum_r |Z_r|^2/q_r^2)
 <(pi^2/8)sum_r |Z_r|^2/q_r^2.                    (PBFA.12)
```

The subtraction `phi_r(K-1)U_r-e_r V_r` is essential. Bounding the two
fluxes separately restores a large artificial interface loss.

## What The Future Anchor Can Cancel

Put

```text
E_C:=sum_(r<=R)x_r^2/q_r^2,
E_perp:=E_C-gamma^2/c.
```

Orthogonal projection onto the endpoint mode vector gives the exact
completion

```text
L_(alpha,K)
 :=sum_(r<=R)(x_r+beta e_r)^2/q_r^2

 =E_perp+c(beta+gamma/c)^2.                       (PBFA.13)
```

Both terms on the right are nonnegative. The future rank-one anchor
can shift only the endpoint direction `gamma/c`; it leaves `E_perp`
unchanged. In particular, it does not algebraically cancel the actual
joined current boundary flux `B_CC-A_CF`.
The correct joined theorem target remains dyadic summability of

```text
K^(-alpha)[E_perp+c(beta+gamma/c)^2]               (PBFA.14)
```

on one fixed cofinal positive-`alpha` sequence. This is exactly the
earlier low-mode energy gate in orthogonal coordinates, not a proof of
that gate.

## Exact Shoulder-Omission Witness

The shoulder is forced even in the simplest completed kernel. Use

```text
P_(K,infinity)(i,j)=pi^2 min(i,j)/N^2
```

at `K=9`, `N=19`, and take

```text
(u_1,...,u_8)=(-1,-1,-1,-1,-1,-1,+1,+1).
```

After suppressing the common factor `pi^2/19^2`,

```text
sum_(i<j)i u_i u_j=0,
(S_1^bd,...,S_7^bd)=(3,6,9,10,9,6,5),
(S_1^sh,...,S_6^sh)=(4,8,10,10,8,4),

sum_(i=1)^7 i S_i^bd
 =sum_(i=1)^6 (i+1)S_i^sh
 =198.
```

Therefore

```text
O_CC=0,
O_(CC,int)=0,
B_CC=198pi^2/19^2,
A_CF=198pi^2/19^2.                                (PBFA.15)
```

This exact bounded-sign witness proves that dropping `A_CF` manufactures
a nonzero remainder from a zero current block. The square-function
target for `Y_(K,r)` in Corollary 11.22Z.7 remains a valid condition for
the actual joined interior. Its complementary target must use
`B_CC-A_CF`, never the two terms separately.

## Finite Mobius Diagnostic

At `alpha=1/2`, direct ordinary-Mobius audits show the same cancellation
pattern at the tested scales:
"""


NOTE_END = r"""

The table is deterministic finite evidence only. It does not prove an
asymptotic bound, but it sharply confirms the shoulder join: at
`K=1024`, `O_(CC,int)` is approximately `2.086e-5`, while `B_CC` and
`A_CF` are approximately `-184.3328` and `-184.1213`. Their joined
difference, not either term, reconstructs the modest current block.

## Proof Boundary

The global flux identity, all four remaining cell formulas, endpoint
smoothing, both interface joins, block compression, boundary/shoulder
spectrum, orthogonal anchor completion, and the `K=9` shoulder witness
are exact and independently checked. The finite Mobius table is
diagnostic only.

No estimate for (PBFA.14), no full reciprocal-tail theorem, no RH,
no PF-infinity, no `Lambda<=0`, and no Clay-prize result is proved.
"""


def build_payload() -> dict:
    rows = build_rows()
    return {
        "kind": (
            "jensen_window_pf_mertens_planar_"
            "boundary_flux_anchor_handoff"
        ),
        "date": "2026-07-24",
        "status": (
            "exact joined boundary-flux and rank-one anchor reduction "
            "with one open Mobius energy gate"
        ),
        "source": {
            "planar_parent": (
                "outputs/jensen_window_pf_mertens_planar_abel_handoff.md"
            ),
            "antidiagonal_parent": (
                "outputs/"
                "jensen_window_pf_mertens_planar_"
                "incidence_antidiagonal_handoff.md"
            ),
            "kernel_parent": (
                "outputs/"
                "jensen_window_pf_mertens_"
                "truncated_half_odd_kernel_handoff.md"
            ),
            "formal_core": "outputs/formal_core.md",
        },
        "parameters": {
            "normalization": "N=2K+1",
            "mode_cutoff": "R=ceil(K^(1-alpha/2))",
            "current_support": "K<n<2K",
            "future_support": "2K<=n<4K",
            "future_weight": (
                "b_n=(2K)^alpha(4K-n)n^(-1-alpha)"
            ),
        },
        "rows": rows,
        "mobius_audit": {
            "alpha": 0.5,
            "rows": MOBIUS_AUDIT_ROWS,
            "status": "finite_diagnostic_only",
        },
        "summary": {
            "row_count": 35,
            "exact_reduction_count": 29,
            "conditional_calibration_count": 1,
            "proof_guard_count": 2,
            "countermodel_count": 1,
            "finite_validation_count": 1,
            "open_joint_energy_gate_count": 1,
            "global_flux_proved": True,
            "all_remaining_cells_projected": True,
            "interface_joins_proved": True,
            "shoulder_join_proved": True,
            "block_compression_proved": True,
            "boundary_cancellation_witness_proved": True,
            "orthogonal_anchor_completion_proved": True,
            "interior_projection_estimate_proved": False,
            "joint_mobius_energy_proved": False,
            "full_reciprocal_tail_proved": False,
            "rh_proved": False,
            "pf_infinity_proved": False,
            "lambda_le_zero_proved": False,
        },
        "proof_boundary": (
            "Exact joined algebra, one countermodel, and finite "
            "diagnostics only. The joint Mobius energy estimate remains "
            "open."
        ),
    }


def render_note(payload: dict) -> str:
    lines = [
        NOTE_BODY.strip(),
        "",
        "| K | R | O_CC | O_CC,int | B_CC | A_CF | B_CC-A_CF | "
        "Anchor correction | O total |",
    ]
    lines.append(
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
    )
    for item in payload["mobius_audit"]["rows"]:
        lines.append(
            "| {K} | {R} | {O_CC:.9g} | {O_CC_int:.9g} | "
            "{B_CC:.9g} | {A_CF:.9g} | {joined_boundary:.9g} | "
            "{anchor_correction:.9g} | "
            "{O_total:.9g} |".format(**item)
        )
    lines.extend([NOTE_END.strip(), "", "## Claim Ledger", ""])
    lines.append("| ID | Role | Status | Statement | Boundary |")
    lines.append("|---|---|---|---|---|")
    for item in payload["rows"]:
        cells = [
            f"`{item['id']}`",
            f"`{item['role']}`",
            f"`{item['status']}`",
            item["statement"],
            item["proof_boundary"],
        ]
        cells = [cell.replace("|", "\\|") for cell in cells]
        lines.append("| " + " | ".join(cells) + " |")
    summary = payload["summary"]
    lines.extend(
        [
            "",
            "Summary:",
            "",
            f"- rows: {summary['row_count']}",
            f"- exact reductions: {summary['exact_reduction_count']}",
            f"- conditional calibrations: "
            f"{summary['conditional_calibration_count']}",
            f"- proof guards: {summary['proof_guard_count']}",
            f"- countermodels: {summary['countermodel_count']}",
            f"- finite validations: {summary['finite_validation_count']}",
            f"- open joint energy gates: "
            f"{summary['open_joint_energy_gate_count']}",
            "",
        ]
    )
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote planar boundary-flux/anchor handoff: "
        f"{len(payload['rows'])} rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
