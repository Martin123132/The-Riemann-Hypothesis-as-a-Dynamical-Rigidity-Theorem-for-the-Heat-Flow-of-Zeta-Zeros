# Jensen-Window PF Lambda-Zero First-Summand Dominance Transfer

Date: 2026-07-13

Status: exact heat-parameter dominance transfer with an open lambda-zero
first-summand curvature tail. This is not a proof of all-shift order-four
positivity, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_lambda0_first_summand_dominance_transfer.json
python work/rh_compute/scripts/jensen_window_pf_lambda0_first_summand_dominance_transfer.py
python work/rh_compute/scripts/check_jensen_window_pf_lambda0_first_summand_dominance_transfer.py
```

## Monotone Heat Transfer

Write `T=-lambda`, so `0<=T<=100`, and normalize the first theta
summand to the probability measure

```text
dmu_(k,T)(u)=Z_(k,T)^(-1)*u^(2k)*exp(-T*u^2)*Phi_1(u)du
```

For the decreasing pointwise higher-summand ratio `epsilon`,

```text
delta_k(T)=E_(mu_(k,T))[epsilon(U)]
delta_k'(T)=-Cov_(mu_(k,T))(epsilon(U),U^2)
2*Cov(f(U),g(U))=E[(f(U)-f(V))*(g(U)-g(V))]
```

Since `epsilon(u)` decreases while `u^2` increases, the two-copy
integrand is nonpositive. Therefore

```text
delta_k'(T)>=0 on 0<=T<=100
0<=delta_k(T)<=delta_k(100)<=2/k^6, k>=300, 0<=T<=100
```

The last inequality uses the certified all-index theorem at `T=100`.
It follows in particular that

```text
0<=delta_k(0)<=2/k^6, every integer k>=300
```

## Order-Four Perturbation

The direct Arb prefix ends at `n=500`; the first unproved tail index is
`n=501`, hence `k=n+3=504`. The same exact Lipschitz calculation now gives

```text
a_j=2*((j-1)^(-6)+2*j^(-6)+(j+1)^(-6)); |B_j-B_j^(1)|<=a_j
|P_n-P_n^(1)|<=10112*(k+1)^2/(k-3)^6<=2/(5*k^2), k=n+3>=504
```

After writing `k=504+n`, the final rational comparison has the
coefficient-positive numerator

```text
n**6 + 3006*n**5 + 3739735*n**4 + 2464014980*n**3 + 906416911135*n**2 + 176398167034206*n + 14175790723241001.
```

The remaining analytic target is now isolated as

```text
K_1(t)<=(7/2)/t^2 for every real t>=503 at lambda=0.
```

This transfer does not supply that curvature theorem or promote the
lambda-zero finite prefix to an all-shift order-four theorem.

```text
outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md
outputs/arb_xi_lambda0_order4_prefix_certificate.md
outputs/jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate.md
outputs/signed_hankel_jensen_dependency_graph.md
```
