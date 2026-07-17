# Jensen-Window PF Compound Order-Five Uniform-Heat Forward Invariance Certificate

Date: 2026-07-13

Status: contiguous and arbitrary-column order-five theorem on `-100<=lambda<=0`.
This is not a proof of the all-degree Jensen bridge, RH, or `Lambda <= 0`; it proves one fixed compound layer, not PF-infinity.

```text
work/rh_compute/results/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_uniform_heat_forward_invariance_certificate.py
```

## Cooperative Propagation

The endpoint theorem now supplies

```text
H_(5,n)(-100)>0 for every integer n>=0
```

The compact-uniform eventual theorem and exact normalized flow are

```text
there exists N_5 such that H_(5,n)(lambda)>0 for every n>=N_5 and -100<=lambda<=0
Q_n'=a_n*Q_(n+1)+b_n*Q_n, a_n=(4*n+34)*H_(4,n)/H_(4,n+1)>0, b_n=((4*n+34)/(4*n+30))*(log H_(4,n+1))'
```

The uniform tail confines a hypothetical first loss to finitely many shifts. Variation of constants on that finite set gives

```text
Q_n(lambda)=E_n(lambda)*(Q_n(-100)+integral_(-100)^lambda E_n(s)^(-1)*a_n(s)*Q_(n+1)(s)ds)
```

whose integrand is nonnegative under backward induction from the positive tail. Therefore

```text
H_(5,n)(lambda)>0 for every n>=0 and -100<=lambda<=0
```

## Arbitrary Columns

The signed contiguous layers through order five satisfy the initial-minor hypotheses pointwise in `lambda`. The fixed-order Gasca-Pena transfer therefore gives

```text
R_(5,n)(j_1,j_2,j_3,j_4,j_5;lambda)>0 for every n>=0, 0<=j_1<j_2<j_3<j_4<j_5, and -100<=lambda<=0
```

Thus contiguous and arbitrary-column order five are both closed on the full heat interval. The next fixed-order task is

```text
derive a cancellation-preserving contiguous order-six entry coordinate and flow reduction
```

```text
outputs/jensen_window_pf_compound_order5_m100_entry_certificate.md
outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md
outputs/jensen_window_pf_order4_noncontiguous_total_positivity_transfer.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
