# Jensen-Window PF Compound Order-Four Eventual-Tail Forward-Invariance Reduction

Date: 2026-07-13

Status: exact finite-confinement forward-invariance reduction with one
open uniform tail. This is not a proof of unconditional order-four
invariance, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_eventual_tail_forward_invariance_reduction.py
```

## Available Inputs

The completed entry and order-three propagation theorems give

```text
H_(4,n)(-100)>0 for every n>=0
H_(3,n)(lambda)<0 on -100<=lambda<=0
Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n>0
```

## Alternative Tail Hypothesis

Instead of bounding every effective diagonal coefficient, assume only

```text
exists N such that Q_n(lambda)>0 for every n>=N and -100<=lambda<=0.
```

The endpoint theorem currently proves only

```text
exists N_H4: H_(4,n)(0)>0 for n>=N_H4.
```

That endpoint statement does not yet control the interior heat interval.

## Finite Confinement

For each fixed `n`, put

```text
E_n(lambda)=exp(integral_(-100)^lambda b_n(s)ds)>0
```

Variation of constants gives

```text
Q_n(lambda)=E_n(lambda)*(Q_n(-100)+integral_(-100)^lambda E_n(s)^(-1)*a_n(s)*Q_(n+1)(s)ds).
```

Under the uniform tail hypothesis, `Q_N` is positive throughout the
interval. Since `Q_n(-100)>0`, `E_n>0`, and `a_n>0`, backward induction
from `N-1` to `0` proves

```text
uniform eventual order-four positivity on [-100,0] implies H_(4,n)(lambda)>0 for every n>=0 and -100<=lambda<=0.
```

This route needs no supremum bound on `beta_n+alpha_n*(n+2)/(n+1)`.
Its sole open input is a compact-heat, uniform version of the available
lambda-zero eventual asymptotic theorem.

```text
outputs/jensen_window_pf_compound_order4_m100_entry_certificate.md
outputs/jensen_window_pf_compound_order4_forward_flow_reduction.md
outputs/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.md
outputs/signed_hankel_jensen_dependency_graph.md
```
