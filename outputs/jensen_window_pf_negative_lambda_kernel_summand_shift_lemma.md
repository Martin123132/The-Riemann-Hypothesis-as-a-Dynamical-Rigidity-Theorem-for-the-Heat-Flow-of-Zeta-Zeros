# Jensen-Window PF Negative-Lambda Kernel Summand-Shift Lemma

Date: 2026-07-10

Status: exact kernel-shift lemma with compact interval certificate. This
is not a proof of cone entry, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_negative_lambda_kernel_summand_shift_lemma`.

Machine-readable result:

```text
work/rh_compute/results/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.json
```

Generator and checker:

```text
python work/rh_compute/scripts/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.py
python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.py
```

Current result:

```text
validated Jensen-window PF negative-lambda kernel summand-shift lemma: 8 rows, 0 issues, 6 exact rows, 1 compact interval row, 1 open far-tail row, 0 ready-to-apply rows
```

## Exact Shift

Write the positive Newman-kernel summands as

```text
phi_n(u)=(2*pi^2*n^4*exp(9u)-3*pi*n^2*exp(5u))*exp(-pi*n^2*exp(4u)).
```

For `a_n=(log n)/2`, direct substitution gives

```text
phi_n(u)=n^(-1/2)*phi_1(u+a_n).
```

Therefore the nth contribution to the raw moment at `lambda=-T` is

```text
M_(k,n)(-T)=2*n^(-1/2)*integral_(a_n)^infinity
  (v-a_n)^(2k)*exp(-T*(v-a_n)^2)*phi_1(v)dv.
```

Relative to the first-summand integrand at the same `v`,

```text
rho_(n,k,T)(v)=n^(-1/2)*(1-a_n/v)^(2k)
                   *exp(T*(2*a_n*v-a_n^2)).
d_v log rho=2*k*a_n/(v*(v-a_n))+2*T*a_n>0.
```

Thus `rho` increases in `v` and decreases in `k`.

## Compact Certificate

At `T=100`, `k>=300`, and `v<=3/2`, only `n=2..20` can occur before
the split because `a_21>3/2`. Arb certifies

```text
sum_(n=2)^20 rho_(n,300,100)(3/2) = [2.121833358028468684385733322599482971419953276137963993084126740175644E-29 +/- 2.63E-100]
upper bound < 2.122e-29
```

Hence the complete shifted `n=2..20` compact contribution is below
`2.122e-29` times the first-summand moment for every `k>=300`.

## Remaining Tail

The revised `lambda=-100` tail theorem now has two explicit pieces:

```text
1. prove the adjacent-k wall for the dominant n=1 moment when k>=300;
2. bound all shifted contributions on v>=3/2 below the resulting margin.
```

The k300 Arb cache supplies the finite collar through `k=300`; this lemma
does not yet supply either infinite-tail piece.

## 2026-07-10 Later Handoff

The later certificate
`outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md`
uses a stronger adaptive split in the original variable and proves the
complete all-k bound `0<=delta_k<=2/k^6`. Thus the shifted `v>=3/2` far-tail
requirement named above is discharged. The only remaining tail theorem is the
quantitative first-summand saddle wall recorded in
`outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md`.

```text
outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md
outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md
outputs/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.md
outputs/signed_hankel_jensen_dependency_graph.md
```

Summary:

The exact identity phi_n(u)=n^(-1/2)phi_1(u+(log n)/2) converts the full Newman kernel into shifted copies of one dominant summand. For lambda=-100, k>=300, and shifted variable v<=3/2, Arb bounds the complete n=2..20 relative contribution below 2.122e-29; n>=21 begins entirely beyond the split. The remaining tail theorem is now the dominant n=1 adjacent-wall inequality plus a v>=3/2 far-tail perturbation bound.
