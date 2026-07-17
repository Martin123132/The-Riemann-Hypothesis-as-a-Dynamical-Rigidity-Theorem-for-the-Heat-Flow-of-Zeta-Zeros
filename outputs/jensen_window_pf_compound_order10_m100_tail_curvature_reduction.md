# Jensen-Window PF Compound Order-Ten Lambda=-100 Tail Curvature Reduction

Date: 2026-07-16

Status: exact positive-tail reduction with one open seventh-nested
curvature ceiling. This does not remove the four rigorous negative
endpoint shifts and is not a proof of RH or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order10_m100_tail_curvature_reduction.json
python work/rh_compute/scripts/jensen_window_pf_compound_order10_m100_tail_curvature_reduction.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order10_m100_tail_curvature_reduction.py
```

## Canonical Coordinate

```text
W(t)=8*B(t)-w(t-1)+2*w(t)-w(t+1)
z(t)=2*w(t)-s(t)+log(1-exp(-W(t)))
Q_(9,n)=A_(n+8)^9*exp(z(n+8))
```

## Order-Ten Sign

```text
E_n=log(Q_(9,n)*Q_(9,n+2)/Q_(9,n+1)^2)=9*log(x_k)+Z_k, Z_k=z(k-1)-2*z(k)+z(k+1), k=n+9
L_n=exp(-E_n)-1
Q_(10,n)>0 iff L_n>0 iff E_n<0
n>=1243 iff k=n+9>=1252
```

## Tail Budget

```text
-9*log(x_k)>=9*d_k>=2259/(250*(2*k+1)), k>=320
Z_k<=5500/k^2 for every integer k>=1252
5500/k^2<2259/(250*(2*k+1)), k>=1252
2259*m**2 + 2906536*m + 96616536>0 for m>=0,
coefficients=['2259', '2906536', '96616536']
```

Thus the full-kernel ceiling would prove `Q_(10,n)(-100)>0` for
every `n>=1243`. It would meet the rigorous finite block
`4<=n<=1242` exactly and give positivity for every `n>=4`.

## Boundary And Handoff

The endpoint signs at `n=0,1,2,3` remain negative. Conditional on
the positive `n>=4` splice, the existing cooperative system and
uniform eventual tail propagate `Q_(10,n)(lambda)>0` for every
`n>=4` and `-100<=lambda<=0`. The four delayed-entry crossings
remain a separate finite-shift theorem.

```text
outputs/jensen_window_pf_compound_order10_m100_finite_splice_certificate.md
outputs/jensen_window_pf_compound_order9_m100_entry_certificate.md
outputs/jensen_window_pf_compound_order9_first_summand_curvature_certificate.md
outputs/jensen_window_pf_all_order_endpoint_heat_reduction.md
outputs/formal_core.md
```
