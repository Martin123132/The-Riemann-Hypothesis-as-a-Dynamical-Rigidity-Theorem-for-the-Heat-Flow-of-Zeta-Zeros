# Jensen-Window PF Cubic Lambda=-100 Tail-Entry Certificate

Date: 2026-07-10

Status: all-k lambda=-100 cubic Jensen entry theorem with open forward-uniform
tail. This is not a proof of higher-degree Jensen hyperbolicity, PF-infinity,
RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_cubic_m100_tail_entry_certificate.json
python work/rh_compute/scripts/jensen_window_pf_cubic_m100_tail_entry_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_cubic_m100_tail_entry_certificate.py
```

Current result:

```text
validated Jensen-window PF cubic lambda=-100 tail-entry certificate: 10 rows, 0 issues, 4074 negative-skewness blocks, 1 analytic negative-skewness ray, 318 prefix margins, 1 all-k cubic tail, 1 full cubic entry theorem, 1 open forward-uniform tail
```

## Negative Skewness

The compact paired engine certifies `kappa_3,t(2*log(U))<0` on all `4074`
mode blocks through `u=5`. Its closest outward-rounded upper bound is

```text
-9.13150910388589240652708356293943489101299308784061615122009E-6
```

On `u>=5`, put `q=pi*exp(4u)`. Exact component bounds give

```text
V''<=5*u^2*q,
V'''>=7*u^3*q,
alpha=V'''/(V'')^(3/2)>=7/(5^(3/2)*sqrt(q)).
```

The analytic ray theorem already proves `|kappa_3+alpha|<=120/q` and
`q>=10^9`; squaring reduces dominance to `49*q>1800000`. Hence the
third cumulant is negative for every real `t>=318`.

## Two-Sided Adjacent Wall

The exact Gamma/cumulant identity and first-summand dominance now give

```text
1/(4*k^2)-16/(k-1)^6 <= L_k
L_k <= -log(1-4/(2*k+1)^2)+16/(k-1)^6,
k>=319.
```

Since `x_k->1`, summing the lower wall proves

```text
-log(x_m)=sum_(j=m)^infinity L_j,
d_m=1-x_m>=1/(5*m+1), m>=320.
```

The latter estimate is exact: after clearing denominators and setting
`n=m-320`, every coefficient is positive.

## Cubic Tail

Using `d_k>=d_(k+1)` and the upper adjacent wall gives

```text
q_(k+1)-q_k<=(100000/99999)^2*(5*k+6)^(3/2)/k^2<1
```

At `k=319` the explicit upper ball is `[0.629526936581100872034573873588795950729998253934931248714367264075224560815668364090713063263859264831287335392570724809083958566675428179651740356280886 +/- 6.77e-154]`.
After squaring and setting `n=k-319`, the final quartic has positive
coefficients

```text
[625120211416829925149906901121, 9139382719646770221315447036, 49046657772633715774210566, 115094896076559489601276, 99996000059999600001]
```

so the bound remains below one for every `k>=319` and tends to zero.
Together with the 318 repaired prefix margins, this proves every shifted
degree-3 Jensen polynomial is hyperbolic at `lambda=-100`.

## Remaining Handoff

Forward cubic propagation now has one explicit analytic requirement:

```text
sup_(lambda in [-100,L]) |q_(k+1)(lambda)-q_k(lambda)| -> 0
for every finite L.
```

The entry theorem and inward boundary algebra are closed. Uniformity of the
spatial tail along the forward heat interval, and every higher-degree minor
cone, remain open.
