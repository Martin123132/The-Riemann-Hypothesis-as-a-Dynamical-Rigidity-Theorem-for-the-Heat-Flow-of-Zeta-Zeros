# Jensen-Window PF Compound Order-Seven Lambda=-100 Tail Curvature Reduction

Date: 2026-07-13

Status: exact endpoint-tail reduction with one open fourth-nested
stable curvature ceiling. This is not a proof of order-seven entry,
PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order7_m100_tail_curvature_reduction.py
```

## Canonical Normalization

The completed hierarchy through order six is

```text
x=exp(-B), d=1-x, g=d^2-x^2*d(t-1)*d(t+1), h=log(g)
f=g^2-x^3*g(t-1)*g(t+1), q=log(f/d)
S(t)=4*B(t)-q(t-1)+2*q(t)-q(t+1), p(t)=2*q(t)-h(t)+log(1-exp(-S(t)))
```

Signed condensation adds exactly one stable gap and logarithm:

```text
T(t)=5*B(t)-p(t-1)+2*p(t)-p(t+1)
r(t)=2*p(t)-q(t)+log(1-exp(-T(t)))
Q_(6,n)=A_(n+5)^6*exp(r(n+5))
```

Indeed the centered `H_5(n+1)^2` scale is `A_(n+5)^10`; division
by `H_4(n+2)=A_(n+5)^4 exp(q(n+5))` leaves the displayed sixth
power. No raw near-cancelling determinant is used.

## Exact Sign Reduction

For `k=n+6`,

```text
E_n=log(Q_(6,n)*Q_(6,n+2)/Q_(6,n+1)^2)=6*log(x_k)+R_k, R_k=r(k-1)-2*r(k)+r(k+1), k=n+6
L_n=exp(-E_n)-1
Q_(7,n)>0 iff L_n>0 iff E_n<0
```

The rigorous prefix proves `Q_(7,n)(-100)>0` for `0<=n<=314`,
so the missing tail starts at `n=315`, equivalently `k=321`.

## Tail Arithmetic

The completed defect theorem gives

```text
-6*log(x_k)>=6*d_k>=753/(125*(2*k+1)), k>=320
```

A deliberately buffered sufficient curvature theorem is

```text
R_k<=900/k^2 for every real/integer k>=321
```

and the comparison is exact:

```text
900/k^2<753/(125*(2*k+1)), k>=321
753*m**2 + 258426*m + 5252373>0 for m>=0,
coefficients=['753', '258426', '5252373']
```

Consequently the displayed scalar ceiling would prove the entire
`n>=315` endpoint tail and splice it to the finite prefix.

## Open Input

```text
Prove R_k<=900/k^2 for every k>=321.
```

The complete-to-first-summand transfer is now proved in
`outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md`.
It reduces the remaining task to the continuous theorem
`r_1''(t)<=600/t^2` on `t>=320`. The shifted-jet compact certificate now
proves this on `320<=t<=1000`, and the nested-curvature compact certificate
extends it through `t=V'(2)`. The finite and asymptotic saddle rays `u>=2`
remain open. The extra stable layer requires a common `t+-5` cover and
potential derivatives through order twelve.

```text
outputs/jensen_window_pf_compound_order7_m100_prefix_certificate.md
outputs/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.md
outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md
outputs/jensen_window_pf_compound_order7_shifted_jet_t320_t1000_certificate.md
outputs/jensen_window_pf_compound_order7_nested_curvature_compact_certificate.md
outputs/formal_core.md
```
