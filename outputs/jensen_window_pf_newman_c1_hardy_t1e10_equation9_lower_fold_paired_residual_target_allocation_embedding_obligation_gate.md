# Fold-owned target allocation and embedding obligation

Date: 2026-08-13

Status: exact 399-to-200 allocation and four-term assembly identity proved;
the analytic selector-to-paired-residual embedding remains open

The exact paired residual contains the half-domain fold block

```text
F_fold=sum_(r=39695)^39894(K_r+K_-r-tauhat_r).       (FA1)
```

Its target roster has 200 modes.  The corrected selector calculation instead
uses the 399-mode augmented roster `39696..40094`.  Exact
event ordering gives a disjoint allocation of that roster to the target:

```text
A_39894={39894},
A_(39894-j)={39894-j,39894+j},       1<=j<=198,
A_39695={40093,40094}.                              (FA2)
```

Thus (FA2) consists of one event-zero group, 198 reflected pairs, and one
selector-exchange edge.  Its 200 target labels are exactly
`39695..39894`, and its 399 augmented modes are exactly
`39696..40094`, with no overlap or omission.  This is an allocation theorem,
not an analytic equality between a group and its target mode.

Let `Z_m` denote the augmented selector contribution after an exact
embedding into the paired-residual coordinates has been supplied, and put

```text
E_emb=sum_(r=39695)^39894(K_r+K_-r)
      -sum_(m=39696)^40094 Z_m.                       (FA3)
```

For `m!=39894`, restore the suppressed common source carrier and define

```text
c_m=e^(i beta^3)2sqrt(2)a_m(G_ex,m-G_opp,m),
C_corr=sum_(m!=39894)c_m=e^(i beta^3)2sqrt(2)W_corr.
```

Its projected value is exactly the quantity certified in Section 11.403.
Pure finite algebra then gives

```text
C_corr=sum_(m!=39894)c_m,
E_0=Z_39894-tauhat_39894,
A_tar=sum_(r=39695)^39893
 [sum_(m in A_r)(Z_m-c_m)-tauhat_r],

F_fold=E_emb+C_corr+E_0+A_tar.                        (FA4)
```

An independent exact-rational test verifies (FA4) on all 399 augmented and
200 target slots.  Applying the common half-domain projection gives

```text
Q_fold=2 Re[e^(-i*pi/8)F_fold].                       (FA5)
```

Section 11.403 fixes the scale and negative sign of the projection of
`C_corr`; it does not show `E_emb=0`.  In fact the existing completed
selector theorem retains the zero, negative, and outer-positive complement
and endpoint half-current as one nonlocal object.  The paired-target theorem
also states explicitly that `tauhat_r` is an algebraic allocation, not a
modewise half-Gamma identity.  Therefore neither completion may be localized
to (FA1) by index counting alone.

The event-zero gate controls only the exact-minus-beta-minus-four transport
at mode 39894; it does not bound `E_0`.  No current gate bounds `A_tar` or
gives an exact source formula for `E_emb`.  The first decisive obligation is
to derive (FA3) from the original paired finite-Poisson current and the
completed one-cell selector identity, preserving all complement sectors.
Only then should `E_emb+E_0+A_tar` be bounded jointly with the signed
`C_corr` projection.

Pi provenance: `pi/8` in (FA5) is the exact odd-square half-Kummer
reflection phase.  The allocation (FA2)--(FA4) introduces no new `pi` and no
fitted constant.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_residual_target_allocation_embedding_obligation_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_residual_target_allocation_embedding_obligation_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_residual_target_allocation_embedding_obligation_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_paired_residual_target_allocation_embedding_obligation_gate.py
```

No selector-to-paired-residual embedding, event-zero target residual,
target-allocation bound, complete fold-owned residual, complete `Q_K-T` or
`T_upper`, all-corridor or height-uniform theorem, `Lambda<=0`, PF-infinity,
RH, or prize-level conclusion is proved.
