# Jensen-Window PF Mellin Multiplier Power-Sum Obstruction

Date: 2026-07-10

Status: interval-certified continuous-interpolation obstruction. This is not a
proof of PF-infinity, Jensen hyperbolicity, RH, or `Lambda <= 0`, and it does
not rule out a multiplier product asserted only at integer coefficient indices.

```text
work/rh_compute/results/jensen_window_pf_mellin_multiplier_power_sum_obstruction.json
python work/rh_compute/scripts/jensen_window_pf_mellin_multiplier_power_sum_obstruction.py
python work/rh_compute/scripts/check_jensen_window_pf_mellin_multiplier_power_sum_obstruction.py
```

Current result:

```text
validated Jensen-window PF Mellin multiplier power-sum obstruction: 11 log moments, 9 power sums, 6 Hankel determinants, 3 negative, 0 inconclusive, 1 continuous route ruled out, 0 discrete routes ruled out, 0 issues
```

## Exact Reduction

For `Re(s)>-1/2`, define

```text
M(s)=integral_0^infinity u^(2s)*Phi(u)du,
A(s)=2*sqrt(pi)*M(s)/(4^s*Gamma(s+1/2)),
B(s)=A(s)*A(0)^(s-1)/A(1)^s.
```

The duplication formula gives `A(k)=A_k` at every nonnegative integer. For
`m>=2`, geometric normalization is linear in `s`, so

```text
d^m log B(0)=kappa_m(2*log U)-psi_(m-1)(1/2).
```

If this interpolation were a continuous elementary multiplier product, then

```text
p_m=(-1)^(m-1)*d^m log B(0)/(m-1)!=sum_j alpha_j^(-m).
```

Every shifted Hankel matrix `[p_(r+i+j)]` would therefore be positive semidefinite.

## Arb Certificate

Arb integrates the finite `n<=10` kernel after `u=exp(-x)` on nine
panels through `x=200`, and after `u=exp(t)` through `u=8`. A common
`1e-45` radius covers the geometric `n`-tail, omitted small-`u` tail, and
double-exponential large-`u` tail for log moments through order 10.

Selected candidate power sums:

```text
p_2 = [0.08866959360452186625898917050448368042996 +/- 3.61e-42]
p_4 = [0.001895020035266623043978694514287596850 +/- 3.80e-40]
p_8 = [3.33434016766525080576860469494e-6 +/- 3.96e-36]
```

Separated Hankel failures:

```text
shift=2, size=4: [-2.588644974358276189292582822e-19 +/- 8.65e-47]
shift=4, size=3: [-2.5739894194580514748560442e-17 +/- 6.00e-43]
shift=7, size=2: [-1.1686954903466295579809210e-14 +/- 2.08e-40]
```

The shift-2 size-4 determinant is strictly negative, so the natural Mellin
interpolation is not such a continuous multiplier product.

## Proof Boundary

The original target asks for equality of coefficient sequences at integer
indices. Equality there does not by itself force equality with this Mellin
interpolation between the integers. A Carlson-type or other uniqueness theorem
with verified growth hypotheses would be needed for that promotion, and none is
claimed here. The discrete counting-measure route therefore remains open, while
its most direct continuous-index strengthening is rigorously closed.
