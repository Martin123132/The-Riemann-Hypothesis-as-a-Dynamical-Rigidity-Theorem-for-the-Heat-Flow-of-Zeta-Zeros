# Selector increment in paired fold coordinates

Date: 2026-08-13

Status: exact selector-cocycle and paired-sector embedding identity proved;
the increment remainder is not bounded

Fix the upper odd endpoint `B=5122423` and let the lower selector advance from
`C=159577` to `C+2=159579`.  If `I_m^(C)` is the finite-Poisson coefficient
before the advance, the integer character is unchanged by the unit shift, so

```text
I_m^(C)(x)-I_m^(C+2)(x)=S_m(x),
K_m^(C)(t)-K_m^(C+2)(t)=Z_m(t),
Z_m(t)=integral_0^(1/2) W_t(x)S_m(x)dx.              (SI1)
```

The endpoint half-current difference is `Z_H`.  The completed symmetric
Poisson theorem therefore gives the exact selector cocycle

```text
P_C=Q_half^(C)-Q_half^(C+2)
   =Z_H+lim_sym sum_m Z_m,                            (SI2)
```

where `P_C` is the single removed source term after the same half-domain
weight and integration.  This identity keeps the endpoint, zero, negative,
and outer-positive sectors attached.

For the fold-owned target block define

```text
F_C=sum_(r=39695)^39894
 [K_r^(C)+K_-r^(C)-tauhat_r].                         (SI3)
```

The target carrier is independent of the selector and cancels exactly in the
state difference.  Hence

```text
D_fold=F_C-F_(C+2)
      =sum_(r=39695)^39894(Z_r+Z_-r).         (SI4)
```

This is the natural location of the one-cell strip: it is an increment of
the paired residual, not a summand of one fixed selector state.

Let `Z_R=sum_(m=39696)^40094Z_m` be the raw positive
399-mode roster.  Exact coefficient cancellation leaves

```text
E_raw=D_fold-Z_R
 =Z_39695+sum_(r=39695)^39894Z_-r
  -sum_(m=39895)^40094Z_m. (SI5)
```

The common positive intersection `39696..39894`
contains exactly 199 modes.  Thus the true raw mismatch is one lower edge,
the 200 negative target modes, and the 200 outer-positive modes.  The earlier
399-to-200 event allocation remains useful internal branch bookkeeping, but
it is not the Fourier-sector embedding (SI5).

There is also a cancellation-safe completed form of the raw mismatch:

```text
C_out=P_C-D_fold,       R_raw=P_C-Z_R,
D_fold-Z_R=R_raw-C_out=E_raw.                         (SI6)
```

Equation (SI6) is cancellation-safe: `R_raw` and `C_out` are completed
objects before their difference is estimated.  It does **not** insert the
corrected selected quantity from Section 11.403.  That quantity is

```text
C_corr^K=e^(i beta^3)2sqrt(2)W_corr,
W_corr=sum_m a_m(G_ex,m-G_opp,m),                     (SI7)
```

The factor `e^(i beta^3)` is mandatory: it is the common source carrier
suppressed in the common-profile coefficients.  Since
`exp(i*pi*(C^2-1)/8)=1`, its paired physical projection is exactly the
quantity in Section 11.403.  The term still carries the nonconstant
leading-defect weights `a_m`, whereas `Z_R` is the unweighted raw strip
roster.  Equal physical units do not make these objects equal.  The next
exact task is to decompose `Z_R` into (SI7) plus its compulsory weighted
companion before attempting a fixed-state sign theorem.

The physical increment is

```text
Delta Q_fold=2 Re[exp(-i*pi/8)D_fold].                (SI8)
```

The `pi/8` is inherited from odd-square Kummer reflection.  The selector
cocycle and coefficient cancellation introduce no new `pi` and no fitted
constant.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_increment_paired_sector_embedding_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_increment_paired_sector_embedding_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_increment_paired_sector_embedding_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_increment_paired_sector_embedding_gate.py
```

No weighted-defect-to-fixed-state embedding, no bound for `E_raw`, no signed
fold-increment theorem, no fixed-state fold residual bound, no all-corridor telescope, complete `Q_K-T` or
`T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
