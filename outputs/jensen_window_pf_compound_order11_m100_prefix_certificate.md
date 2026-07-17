# Order-Eleven Lambda=-100 Prefix Certificate

Date: 2026-07-16

Status: rigorous finite endpoint theorem. This is not a proof of the
analytic order-eleven tail, heat propagation, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order11_m100_prefix_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order11_m100_prefix_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order11_m100_prefix_certificate.py
```

## Exact Reduction

```text
epsilon_11=(-1)^55=-1, Q_(11,n)=-H_(11,n)
Q_(11,n)*Q_(9,n+2)=Q_(10,n+1)^2-Q_(10,n)*Q_(10,n+2)
M_n=Q_(10,n+1)^2/(Q_(10,n)*Q_(10,n+2))-1
```

The four rows `n=0,1,2,3` lie outside the positive `Q10` cone, so
their condensation numerators are enclosed directly through the stable
chain and audited with independent `11x11` determinants. For
`4<=n<=1242`, the relative `Q10` margin is rigorously positive.

## Finite Theorem

```text
Q_(10,n+1)^2-Q_(10,n)*Q_(10,n+2)>0 for 0<=n<=3
M_n>0 for every 4<=n<=1242
Q_(11,n)(-100)>0 for every 0<=n<=1242
minimum relative margin at n=1242: [0.0050653409726380927142089621927090981269206001522366632107599374870969111066134147 +/- 2.76E-83]
```

There are `1243` positive rows, zero negative rows, and zero
inconclusive rows. The remaining endpoint handoff begins exactly at
`n=1243` and requires the continuous order-eleven curvature theorem.
