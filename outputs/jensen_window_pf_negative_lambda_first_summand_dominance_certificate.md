# Jensen-Window PF Negative-Lambda First-Summand Dominance Certificate

Date: 2026-07-10

Status: all-k analytic first-summand dominance certificate. This is not
a proof of cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_first_summand_dominance_certificate`.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_dominance_certificate.py
```

Current result:

```text
validated Jensen-window PF negative-lambda first-summand dominance certificate: 10 rows, 0 issues, 4 exact rows, 5 interval rows, 15 positive analytic gates, 1 open dominant-wall row, 0 ready-to-apply rows
```

## Exact Kernel Ratio

For `q=pi*exp(4u)`, divide the nth positive Newman summand by the
first:

```text
r_n(q)=n^2*(2*n^2*q-3)/(2*q-3)*exp(-(n^2-1)*q).
d_q log(r_n)=-(n^2-1)*(1+6/((2*n^2*q-3)*(2*q-3)))<0.
```

Thus every tail ratio and their sum `epsilon(u)` decrease with `u`.
A finite Arb sum through `n=20` plus a geometric tail gives

```text
epsilon(0) = [0.002176061157485405504580960975113297499158520899631047145786721045566703 +/- 2.56E-73]
epsilon(0) < 0.0022.
```

## Adaptive Split

For every integer `k>=300`, set

```text
a(k)=log(k)/8
b(k)=a(k)+1/10
c(k)=b(k)+1/100.
```

If `S_k` is the logarithm of the `n=1`, `lambda=-100` moment
integrand, then

```text
S_k'(u)=2*k/u-200*u+5+8*q/(2*q-3)-4*q
S_k''(u)=-2*k/u^2-200-16*q-96*q/(2*q-3)^2<0.
```

Endpoint bounds and monotone comparison functions prove
`S_k'(a(k))>100` and `S_k'(c(k))>0` for the complete half-line
`k>=300`. Concavity therefore gives

```text
P_k(u<a(k)) <= exp(S_k(a)-S_k(b))/(0.01*S_k'(a)).
epsilon(0)*P_k(u<a(k)) < k^(-6).
```

On the complementary region, the decreasing kernel ratio and
`q(a(k))=pi*sqrt(k)` give

```text
epsilon(a(k)) <= 17*exp(-3*pi*sqrt(k)) < k^(-6).
```

All 15 endpoint and monotonic-propagation gates are Arb-positive.
The half-line propagation is exact: with `L=log(k)`, the comparison
ratios have logarithmic derivatives

```text
1/2+1/(L-alpha)-2/(L+beta)>0
1  +1/(L-alpha)-2/(L+beta)>0
```

for the three pairs `(alpha,beta)=(1,0),(3/25,22/25),`
`(1/5,4/5)` and every `L>=log(300)`. Hence each split endpoint
comparison only strengthens with `k`. Also

```text
d_k[log(17)-3*pi*sqrt(k)+6*log(k)]
  =(6-(3*pi/2)*sqrt(k))/k<0.
```

## All-k Consequence

Let `M_k^(1)` be the first-summand moment and `M_k` the full kernel
moment. The two regions compose to prove

```text
M_k=M_k^(1)*(1+delta_k), 0<=delta_k<=2/k^6, k>=300.
```

For the adjacent log wall

```text
L_k=log(x_(k+1)/x_k),
|L_k-L_k^(1)|<=16/(k-1)^6, k>=301.
```

The shifted `v>=3/2` far-tail source from the earlier shift lemma is
therefore discharged. It remains to prove a quantitative first-summand
wall, for example

```text
L_k^(1)>=1/(4*k^2), k>=301.
```

That bound would dominate the certified perturbation by many orders of
magnitude and splice to the repaired finite collar.

```text
outputs/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.md
outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

Exact monotonicity of every original-variable tail ratio, strict concavity of the first-summand tilted integrand, and 15 positive Arb endpoint gates prove that at lambda=-100 the complete n>=2 moment tail satisfies 0<=delta_k<=2/k^6 for every integer k>=300. Consequently the full adjacent log wall differs from the n=1 wall by at most 16/(k-1)^6 for k>=301. The shifted far-tail source is discharged; the only remaining tail inequality is a quantitative dominant n=1 saddle wall.
