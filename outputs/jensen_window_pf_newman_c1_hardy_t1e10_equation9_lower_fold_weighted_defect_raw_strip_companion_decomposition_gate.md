# Weighted defect and raw-strip companion decomposition

Date: 2026-08-13

Status: exact carrier-restored insertion into the selector increment proved;
the compulsory joint companion is not bounded

Let

```text
kappa_C=e^(i beta^3)2sqrt(2),       beta^3=pi*C^2/8,
Z_m=kappa_C G_ex,m.                                  (WD1)
```

The factor `kappa_C` restores both the exact equation-(9) half-mode
normalization and the common source carrier suppressed in the common-profile
Fourier coefficients.  Extend the certified leading-defect weights by the
existing event-zero convention `a_39894=0`, and put

```text
c_m=kappa_C a_m(G_ex,m-G_opp,m),
b_m=kappa_C[(1-a_m)G_ex,m+a_m G_opp,m].              (WD2)
```

Direct coefficient algebra, with no estimate, gives

```text
Z_m=c_m+b_m.                                         (WD3)
```

Thus on the augmented roster `R=39696..40094`,

```text
Z_R=sum_R Z_m=C_corr+B_comp,
C_corr=sum_R c_m=kappa_C W_corr,
B_comp=sum_R b_m.                                    (WD4)
```

The event-zero mode has `c_39894=0` and
`b_39894=Z_39894`; it has not disappeared.  Combining (WD4) with the exact
paired-sector identity of Section 11.405 proves

```text
D_fold=C_corr+B_comp+E_raw,                           (WD5)
P_C=C_corr+B_comp+E_raw+C_out.                       (WD6)
```

Equation (WD6) is the completed one-cell source identity.  It retains the
coefficient companion, lower edge, negative target modes, outer-positive
modes, and completed outside sector before any absolute value is taken.

The paired physical projection splits exactly as

```text
Delta Q_fold=R_corr+R_joint,
R_corr =2Re[e^(-i*pi/8)C_corr]=4sqrt(2)Re W_corr,
R_joint=2Re[e^(-i*pi/8)(B_comp+E_raw)].              (WD7)
```

The second equality uses
`exp(i*pi*(C^2-1)/8)=1`, with `(C^2-1)/8` even.  The existing interval
certificate therefore gives

```text
R_corr<[-0.00012633020694136396822650036607776931238404392877833648447949440001439331656143490 +/- 2.31e-84]<-1.263e-4,
|R_corr|<[0.00023129817306407708082063016174763530876611806139822690924052796930526831295442464 +/- 2.19e-84]<2.314e-4. (WD8)
```

There is a sharper cancellation.  Split `C_corr=C_inner+C_outer` at event
zero and define

```text
J_pair=Z_39695+Z_-39695
       +sum_(m=39696)^39894(b_m+Z_-m).               (WD9)
```

Using `Z_m=c_m+b_m` on the 200 outer-positive modes gives

```text
B_comp+E_raw=J_pair-C_outer,
D_fold=C_inner+J_pair,
P_C=C_inner+J_pair+C_out.                            (WD10)
```

Thus the outer weighted correction cancels exactly from the operative fold
increment.  The signed theorem (WD8) concerns `C_inner+C_outer`; it does not
by itself determine the sign of `C_inner`.  The next quantitative step is to
certify the inner corrected projection and then estimate `J_pair` as one
paired object.  Bounding `B_comp` and `E_raw` separately is not licensed.

Pi provenance: `beta^3=pi*C^2/8` and `pi/8` are inherited from the exact
Kummer phase and odd-square half-domain reflection.  No fitted scale or
geometric insertion is used.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_defect_raw_strip_companion_decomposition_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_defect_raw_strip_companion_decomposition_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_defect_raw_strip_companion_decomposition_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_weighted_defect_raw_strip_companion_decomposition_gate.py
```

No bound for `R_joint`, signed fold increment, fixed-state fold residual,
all-corridor telescope, complete `Q_K-T` or `T_upper`, height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
