# Jensen-Window PF Compound Order-Eight Uniform-Heat Forward Invariance Certificate

Date: 2026-07-13

Status: signed contiguous and arbitrary-column order-eight theorem on
`-100<=lambda<=0`. This is not a proof of the all-degree Jensen
bridge, RH, or `Lambda <= 0`; it proves one fixed compound layer,
not PF-infinity.

```text
work/rh_compute/results/jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate.py
```

## Cooperative Propagation

The endpoint theorem supplies

```text
Q_(8,n)(-100)=H_(8,n)(-100)>0 for every integer n>=0
```

The compact-uniform eventual theorem and exact signed flow are

```text
there exists N_8 such that Q_(8,n)(lambda)>0 for every n>=N_8 and -100<=lambda<=0
Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n=(4*n+58)*Q_(7,n)/Q_(7,n+1)>0, b_n=((4*n+58)/(4*n+54))*(log Q_(7,n+1))'
```

The uniform tail confines a hypothetical first loss to finitely many
shifts. Variation of constants on that finite set gives

```text
Q_n(lambda)=E_n(lambda)*(Q_n(-100)+integral_(-100)^lambda E_n(s)^(-1)*a_n(s)*Q_(n+1)(s)ds)
```

whose integrand is nonnegative under backward induction from the
positive tail. Therefore

```text
Q_(8,n)(lambda)=H_(8,n)(lambda)>0 for every n>=0 and -100<=lambda<=0
```

## Arbitrary Columns

The signed contiguous layers through order eight satisfy the
initial-minor hypotheses pointwise in `lambda`. The fixed-order
Gasca-Pena transfer therefore gives

```text
epsilon_8*R_(8,n)(j_1,j_2,j_3,j_4,j_5,j_6,j_7,j_8;lambda)>0 for every n>=0, 0<=j_1<j_2<j_3<j_4<j_5<j_6<j_7<j_8, and -100<=lambda<=0
```

Thus contiguous and arbitrary-column order eight are both closed on
the full heat interval. The remaining programme is genuinely all-order:

```text
extend the zeta-specific endpoint and cooperative-flow programme to order nine while deriving a uniform-in-order mechanism
```

```text
outputs/jensen_window_pf_compound_order8_m100_entry_certificate.md
outputs/jensen_window_pf_compound_order8_uniform_tail_flow_reduction.md
outputs/jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate.md
outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
