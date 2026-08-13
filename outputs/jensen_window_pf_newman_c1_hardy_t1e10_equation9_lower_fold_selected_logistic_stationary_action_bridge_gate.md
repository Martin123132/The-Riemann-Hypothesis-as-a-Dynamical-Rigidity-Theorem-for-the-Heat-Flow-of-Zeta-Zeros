# Selected fold/logistic stationary-action bridge

Date: 2026-08-13

Status: exact all-mode top-corridor stationary carrier map certified; the
uniform finite-integral amplitude match, complete ordinary splice, `T_upper`,
and RH remain open

Put

```text
s=4m/C,  mu=(t*-t)/t*=lambda/beta^2,  beta^3=t*.       (SA1)
```

Solving the exact logistic and selected canonical joint stationary equations
gives

```text
q_e=(1/2)log((1-mu)/s^2),
y_e/beta^2=((s-1)^2-mu)/s,
z_c/beta=1-sqrt(2s+mu-1),
y_c/beta^2=2[s-sqrt(2s+mu-1)].                         (SA2)
```

Substitution into the two stationary actions simplifies exactly to

```text
Theta_e-Theta_c=beta^3 F(s,mu),                        (SA3)

6F=3s^2-8sR+12s-4mu R+6mu log s-3mu log(1-mu)
   +9mu+4R-6log s+3log(1-mu)-11,
R=sqrt(2s+mu-1).                                       (SA4)
```

This is the finite-height carrier bridge: the exact logistic carrier is the
selected canonical carrier multiplied by `exp(i*beta^3*F)`.  It is not an
asymptotic phase declaration.

At the selector center, symbolic differentiation proves

```text
F(s,0)=(s-1)^5/20-(s-1)^6/8+13(s-1)^7/56+...;          (SA5)
```

the derivatives through order four vanish and the fifth derivative is `6`.
The first height derivative has the lower-order but still rigid tangency

```text
F_mu(s,0)=-(s-1)^3/6+3(s-1)^4/8-27(s-1)^5/40+... .    (SA6)
```

All `399` modes `m=39696,...,40094` were enclosed on
the full top corridor `t*-pi/16<=t<=t*`.  Both stationary normal coordinates
remain inside the one-cell strip, with

```text
0 < y_e <= 115.7846195742020540819794405251741409301757812500000000000000000000000000000000000000000000000000000,
0 < y_c <= [115.7853374149729823009277131705244766524074069105512803421828717194283298592154790119897362826789336 +/- 5.15e-100],
Y=[116.3473994433307293304405864127926592709684040724544140144822587448142837461736192118000150309760753 +/- 5.02e-98].                            (SA7)
```

The uniform action and multiplier bounds are

```text
max_m sup_t |Theta_e-Theta_c|
 < [0.001553939890085833702727967179932554524723632822141271740851090546738296584448362628543893575103730246 +/- 3.86e-12] < 0.001555,
|exp(i[Theta_e-Theta_c])-1| <= |Theta_e-Theta_c|.       (SA8)
```

The maximum occurs at mode `40094`.  Its center defect is
`[0.001553935828165629707118148436169587307247550749305690297814652620619088839989719469114436462330482971 +/- 1.86e-90]` and is rigorously nonzero, so
the multiplier cannot be replaced by an identity.  It must remain inside the
199-pair projection before summation or absolute values.

The exact scalar amplitude and the two Hessian determinants simplify at the
same saddles.  Their leading stationary-coefficient ratio is

```text
R_A(s,mu)=[(1-mu)(2s+mu-1)]^(1/4)/sqrt(s),             (SA9)
L_exact/L_canonical=R_A exp(i*beta^3 F).               (SA10)
```

At `mu=0`, `R_A=1-(s-1)^2/4+(s-1)^3/2+...`.  Across the
full top corridor, the exact inequality

```text
|R_A exp(i*Delta)-1| <= |R_A-1|+R_A|Delta|
```

gives the carrier-scaled bounds

```text
max_m sup_t |R_A exp(i*Delta)-1|/sqrt(m)
 < [7.791527784240833546356591078618780699116541654802858829498291015625000000000000000000000000000000000e-6 +/- 4.37e-14] < 7.8e-6,
physical < [3.507419795971723210618594366244496537629308893611396327102485483498293716936664475003239965201315984e-6 +/- 1.98e-14].
                                                               (SA11)
```

This consumes less than the existing local normalized headroom
`[1.684762986976495868506729268152917255910224328644042473174483323829871500000000000000000000000000000e-5 +/- 3.40e-21]` and leaves at least
`[9.056102041924121624183921674257025768134916762471474817766256578142465000000000000000000000000000000e-6 +/- 1.68e-106]` for the uniform
endpoint/Fresnel and amplitude-variation remainder.  These are local
one-carrier diagnostics, not a 399-term triangle estimate and not an aggregate
`T_upper` budget.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_logistic_stationary_action_bridge_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_logistic_stationary_action_bridge_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_logistic_stationary_action_bridge_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_logistic_stationary_action_bridge_gate.py
```

No uniform finite-integral selected-branch amplitude theorem, grouped 199-pair error bound, all-400-
corridor continuation, complete `Q_K-T` or `T_upper`, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
