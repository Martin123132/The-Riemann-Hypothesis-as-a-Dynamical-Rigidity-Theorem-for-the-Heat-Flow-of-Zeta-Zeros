# Exact common-profile initial slope

Date: 2026-08-13

Status: rigorous exact-minus-beta-minus-four initial slope on the top selector
corridor; this is not yet the full profile or a proof of `T_upper` or RH

After removing the common quadratic phase, write the exact transformed fold as

```text
H_ex(lambda,y)=int_R (1+y/(2beta^2)) cosh(z/beta)^(-3/2)
 exp(i{beta^3[z/beta-tanh(z/beta)]-lambda*z-beta*y*tanh(z/beta)}) dz. (IS1)
```

At `y=0`, direct differentiation gives

```text
d_y H_ex=int_R cosh(z/beta)^(-3/2) exp(i Phi_0)
 [1/(2beta^2)-i beta tanh(z/beta)] dz.                 (IS2)
```

The exact and beta-minus-four derivatives share the carrier
`exp(i[z^3/3-lambda*z])`.  Subtracting their ratio multipliers before taking
absolute values, then using a degree-17 tanh polynomial on `0<=r<=12`, proves

```text
sup_(0<=lambda<=pi/(16beta)) |d_y(H_ex-H4)(lambda,0)|
 <= abs([-5.5287800957085793122985763322674881918496126811619178281999605932129114076136882931034364178393e-22 +/- 2.18e-24])
 < 5.56e-22.                                             (IS3)
```

The appended compact-replacement and contour-tail radius is
`[2.3144403817744743910823600119361744808045874052217964208801416849950048517028849592092911675913e-28 +/- 2.32e-116]`.  The independently derived
beta-minus-six coefficient gives

```text
d_y H6(lambda,0)=[-5.5287563502829390685218300270665205844532925088468823949627789677080155552346599900569668846348e-22 +/- 2.39e-27],
d_y(H_ex-H4-H6)(lambda,0)
 =[-2.3745425640243776746305200967607396320172315035433237181625504895852379028303046469533204543168e-27 +/- 2.18e-24].           (IS4)
```

The left Airy ray is minus the conjugate of the right ray, so (IS3) is a
two-ray real integral, not a one-sided estimate.  Pi enters only through the
Airy-ray angle `pi/6`, `beta^3=pi*C^2/8`, and the selector width
`lambda_max=pi/(16beta)`.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_profile_initial_slope_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_profile_initial_slope_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_profile_initial_slope_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_profile_initial_slope_gate.py
```

No full exact-minus-beta-minus-four profile bound, weighted finite-height
remainder, complete `T_upper`, height-uniform theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
