# Fourth Theta-Summand Tangency Homotopy Scout

Date: 2026-08-03

Status: high-precision branch diagnostic. This is not a proof of RH or contact exclusion.

```text
work/rh_compute/results/jensen_window_pf_newman_theta_fourth_summand_tangency_homotopy_scout.json
python work/rh_compute/scripts/jensen_window_pf_newman_theta_fourth_summand_tangency_homotopy_scout.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_fourth_summand_tangency_homotopy_scout.py
```

## Exact Homotopy

```text
F_lambda(t,x)=H_(1,t)(x)+H_(2,t)(x)+H_(3,t)(x)+lambda H_(4,t)(x),
F_lambda=0,  partial_x F_lambda=0.
```

At a contact, heat flow gives the exact Jacobian and response

```text
J=[[-F_xx,0],[-F_xxx,F_xx]],  det J=-F_xx^2,
(t',x')=-J^(-1)(H_4,H_4').
```

## Pi Provenance

pi is the canonical Jacobi-theta constant in sum_n exp(-pi*n^2*y), here y=exp(4u). It is fixed by the theta modular normalization inherited by Xi; it is not inserted from an arbitrary circle or polygon. Replacing it changes the theta kernel and no longer represents the same zeta problem.

## Continuation

| lambda | t | x | dt/dlambda | dx/dlambda | F_xx |
|---:|---:|---:|---:|---:|---:|
| 0.0 | 0.432004524261 | 135.56163204268 | -1.91549903 | 1.342158327 | -1.3214387e-21 |
| 0.02 | 0.393303528925 | 135.58873204827 | -1.954873095 | 1.368008013 | -1.2946571e-21 |
| 0.04 | 0.353798344418 | 135.6163591904 | -1.995935873 | 1.394881936 | -1.2678563e-21 |
| 0.06 | 0.313454089225 | 135.644534567 | -2.03879971 | 1.422842224 | -1.2410355e-21 |
| 0.08 | 0.272233534017 | 135.67328056875 | -2.083587231 | 1.451956109 | -1.214194e-21 |
| 0.1 | 0.230096884161 | 135.70262098636 | -2.130432556 | 1.482296453 | -1.1873308e-21 |
| 0.12 | 0.187001536264 | 135.73258112904 | -2.17948269 | 1.51394235 | -1.160445e-21 |
| 0.14 | 0.142901804885 | 135.76318795563 | -2.230899119 | 1.54697979 | -1.1335357e-21 |
| 0.16 | 0.0977486149058 | 135.79447022004 | -2.284859654 | 1.581502423 | -1.1066017e-21 |
| 0.18 | 0.0514891541925 | 135.82645863284 | -2.341560568 | 1.617612415 | -1.079642e-21 |
| 0.2 | 0.00406648019228 | 135.85918604105 | -2.401219076 | 1.65542142 | -1.0526552e-21 |
| 0.22 | -0.0445809271325 | 135.89268762879 | -2.464076233 | 1.69505169 | -1.0256401e-21 |

## Boundary Crossing

The tracked branch reaches `t=0` at `lambda=0.201691679922127687` and `x=135.86198925655846`.
The transverse derivative is `dt/dlambda=-2.4064083881484`.
Thus this branch exits nonnegative heat time after only about one
fifth of the true fourth coefficient has been switched on.

## Interpretation

This supplies a concrete drifting-root mechanism. The fixed-root
tangency-chain theorem does not cover it, because the contact location
moves and the fourth-summand response is order one despite its small mass.
The next serious theorem is an interval branch-and-count argument, not a
claim that this one tracked branch exhausts the four-summand problem.

## Proof Boundary

This is a high-precision, analytically tailed point continuation, not interval continuation. It tracks one branch only. It does not exclude another four-summand branch, prove no re-entry for lambda>the crossing, control the complete theta sum, exclude an Xi contact, prove Lambda<=0, RH, or a prize-level conclusion.
