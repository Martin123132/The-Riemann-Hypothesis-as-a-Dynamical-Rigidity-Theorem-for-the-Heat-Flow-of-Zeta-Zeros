# Jensen-Window PF Compound Order-Nine Lambda=-100 Tail Curvature Reduction

Date: 2026-07-13

Status: exact endpoint-tail reduction with one open sixth-nested
stable curvature ceiling. This is not a proof of order-nine entry,
PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order9_m100_tail_curvature_reduction.py
```

## Canonical Normalization

```text
U(t)=6*B(t)-r(t-1)+2*r(t)-r(t+1); s(t)=2*r(t)-p(t)+log(1-exp(-U(t)))
V(t)=7*B(t)-s(t-1)+2*s(t)-s(t+1)
w(t)=2*s(t)-r(t)+log(1-exp(-V(t)))
Q_(8,n)=A_(n+7)^8*exp(w(n+7))
```

## Exact Sign Reduction

```text
E_n=log(Q_(8,n)*Q_(8,n+2)/Q_(8,n+1)^2)=8*log(x_k)+Y_k, Y_k=w(k-1)-2*w(k)+w(k+1), k=n+8
M_n=exp(-E_n)-1
Q_(9,n)>0 iff M_n>0 iff E_n<0
k=n+8, so n>=1241 iff k>=1249
```

## Tail Arithmetic

```text
-8*log(x_k)>=8*d_k>=1004/(125*(2*k+1)), k>=320
Y_k<=4900/k^2 for every real/integer k>=1249
4900/k^2<1004/(125*(2*k+1)), k>=1249
1004*m**2 + 1282992*m + 35603504>0 for m>=0,
coefficients=['1004', '1282992', '35603504']
```

The scalar ceiling would prove the complete `n>=1241` endpoint tail
and splice it to the finite prefix.

## Open Input

```text
Prove Y_k<=4900/k^2 for every k>=1249.
```

A continuous proof requires a common `t+-7` cover and `H` derivatives
through order sixteen before tent integration.

```text
outputs/jensen_window_pf_compound_order9_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order9_uniform_tail_flow_reduction.md
outputs/formal_core.md
```
