# Jensen-Window PF Compound Order-Eight Lambda=-100 Tail Curvature Reduction

Date: 2026-07-13

Status: exact endpoint-tail reduction with one open fifth-nested
stable curvature ceiling. This is not a proof of order-eight entry,
PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order8_m100_tail_curvature_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order8_m100_tail_curvature_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_m100_tail_curvature_reduction.py
```

## Canonical Normalization

The completed order-seven hierarchy supplies

```text
T(t)=5*B(t)-p(t-1)+2*p(t)-p(t+1); r(t)=2*p(t)-q(t)+log(1-exp(-T(t)))
```

Signed condensation adds one gap and one stable logarithm:

```text
U(t)=6*B(t)-r(t-1)+2*r(t)-r(t+1)
s(t)=2*r(t)-p(t)+log(1-exp(-U(t)))
Q_(7,n)=A_(n+6)^7*exp(s(n+6))
```

The centered Q_6(n+1)^2 scale is A_(n+6)^12. Division by
H_5(n+2)=A_(n+6)^5 exp(p(n+6)) leaves the seventh power.
No raw near-cancelling determinant is used.

## Exact Sign Reduction

For k=n+7,

```text
E_n=log(Q_(7,n)*Q_(7,n+2)/Q_(7,n+1)^2)=7*log(x_k)+W_k, W_k=s(k-1)-2*s(k)+s(k+1), k=n+7
M_n=exp(-E_n)-1
Q_(8,n)>0 iff M_n>0 iff E_n<0
```

The rigorous prefix proves Q_(8,n)(-100)>0 for 0<=n<=1242, so
the missing tail starts at n=1243, equivalently k=1250.

## Tail Arithmetic

The completed coefficient-defect theorem gives

```text
-7*log(x_k)>=7*d_k>=1757/(250*(2*k+1)), k>=320
```

A buffered sufficient curvature theorem is

```text
W_k<=4300/k^2 for every real/integer k>=1250
```

and the comparison is exact:

```text
4300/k^2<1757/(250*(2*k+1)), k>=1250
1757*m**2 + 2242500*m + 56737500>0 for m>=0,
coefficients=['1757', '2242500', '56737500']
```

Consequently the scalar ceiling would prove the complete n>=1243
endpoint tail and splice it to the finite prefix.

## Open Input

```text
Prove W_k<=4300/k^2 for every k>=1250.
```

A continuous proof will require one common t+-6 cover and H
derivatives through order fourteen before tent integration.

```text
outputs/jensen_window_pf_compound_order8_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order8_uniform_tail_flow_reduction.md
outputs/formal_core.md
```
