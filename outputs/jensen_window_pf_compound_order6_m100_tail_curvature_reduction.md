# Jensen-Window PF Compound Order-Six Lambda=-100 Tail Curvature Reduction

Date: 2026-07-13

Status: exact endpoint-tail reduction with one open third-nested stable
curvature ceiling. This is not a proof of order-six entry, PF-infinity,
RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order6_m100_tail_curvature_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order6_m100_tail_curvature_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_m100_tail_curvature_reduction.py
```

## Canonical Normalization

The apparently complicated positive monomial in the stable H5
factorization collapses exactly:

```text
A_n^5*rho_n^20*x_(n+1)^15*x_(n+2)^10*x_(n+3)^5=A_(n+4)^5
S(t)=4*B(t)-q(t-1)+2*q(t)-q(t+1)
p(t)=2*q(t)-h(t)+log(1-exp(-S(t)))
H_(5,n)=A_(n+4)^5*exp(p(n+4))
```

Thus order six adds exactly one stable logarithm to the completed
order-five hierarchy; no raw near-cancelling determinant is used.

## Exact Sign Reduction

For `k=n+5`,

```text
D_n=log(H_(5,n)*H_(5,n+2)/H_(5,n+1)^2)=5*log(x_k)+P_k, P_k=p(k-1)-2*p(k)+p(k+1), k=n+5
K_n=exp(-D_n)-1
Q_(6,n)>0 iff K_n>0 iff D_n<0
```

The rigorous prefix proves `Q_(6,n)(-100)>0` for `0<=n<=316`,
so the missing tail starts at `n=317`, equivalently `k=322`.

## Tail Arithmetic

The completed defect theorem gives

```text
-5*log(x_k)>=5*d_k>=251/(50*(2*k+1)), k>=320
```

A deliberately loose sufficient curvature theorem is

```text
P_k<=320/k^2 for every real/integer k>=322
```

and the comparison is exact:

```text
320/k^2<251/(50*(2*k+1)), k>=322
251*m**2 + 129644*m + 15704684>0 for m>=0,
coefficients=['251', '129644', '15704684']
```

Consequently the displayed scalar ceiling would prove the entire
`n>=317` endpoint tail and splice it to the finite prefix.

## Open Input

```text
Prove P_k<=320/k^2 for every k>=322.
```

The companion first/full bridge separates that task into a continuous
first-summand curvature theorem and an explicit kernel-transfer error.

```text
outputs/jensen_window_pf_compound_order6_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.md
outputs/formal_core.md
```
