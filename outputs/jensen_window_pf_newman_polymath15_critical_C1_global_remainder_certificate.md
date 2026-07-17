# Jensen-Window PF Newman Polymath-15 Critical C1 Global Remainder Certificate

Date: 2026-07-17

Status: global corrected `C1` remainder on the `L>=50` critical ray,
including every cutoff transition. This is not a proof of
`Lambda <= 0` or RH; corrected transversality remains open.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.py
```

## Adjacent Ratio

The exact C0 recurrence reduces a cutoff jump to one main block minus
one endpoint block. On a doubled transition collar put

```text
T=z/2+pi*t/8, a=sqrt(T/(2*pi)), h=1/2-i*pi*t/8, s=i*T+h, delta=log(n/a), r=(n-a)/a
r0=alpha(iT)-(log(a)+i*pi/4), r1=alpha(s)-(log(a)+i*pi/4), R_M=log(M0(s))-log(M0(iT))-h*alpha(iT)
```

The recurrence phase has the exact residual

```text
Psi=a^2*(2*log(1+r)-2*r+r^2)=O(a^2*r^3)
```

and direct cancellation gives

```text
log(main_n/endpoint_n)=h*r0+R_M-h*delta+t/4*(i*pi*(r1-delta)/2+(r1-delta)^2)-i*pi*Psi
```

The elementary bounds saved in the artifact yield

```text
|h| <= [0.5371714271922468021254818134070188143898179360655684061240855679759843 +/- 3.97e-71] < 3/5
log-ratio constant = [1.784000000000000000000000000000000000000000000000000000000000000000000 +/- 3e-74] < 3
|log(main_n/endpoint_n)|<3/|T| and |main_n/endpoint_n-1|<6/|T|
```

A separate direct complex-point audit compares the block quotient
with the exponential of the displayed exact residual. Across
`33` points its maximum relative discrepancy is
`1.48706884775618692828577350123e-87`, while the largest
sampled `|log ratio|*|T|` is `0.819771075241446295806733521638`.
This finite audit only guards signs and branches; the uniform bound
above comes from the analytic estimates.

## Transition Bound

Combining the ratio estimate with the explicit endpoint size gives

```text
|endpoint_n|/A_t(x)<50*exp(-L/4) on the doubled transition collar
direct adjacent constant = [190.9859317102744029226605160470172344413515748885477384972008128706762 +/- 4.29e-68] < 200
|Delta corrected lift|/A_t(x)<1000*exp(-5L/4)
```

The saved constant `1000` is deliberately loose. At `L=50` its
effective contribution to the `exp(-3L/4)` collar constant is only
`[1.388794386496402059466176374608685691039976038020505558354779996088137e-8 +/- 1.87e-78]`. Therefore

```text
2300*exp(-3L/4)+1000*exp(-5L/4)<2500*exp(-3L/4) for L>=50
For every critical radius-1/L disk, including cutoff crossings, sup |R_hat|/A_t(x)<2500*exp(-3L/4)
```

## Global C1 Budget

Cauchy's estimate now applies with one fixed analytic cutoff on every
center disk, whether or not the prescribed cutoff changes inside it:

```text
|r|<2500*exp(-3L/4), |r'|<5000*L*exp(-3L/4), r^2+(r'/L)^2<32000000*exp(-3L/2)
```

The approximation and cutoff sides of the `L>=50` critical problem
are therefore closed. The remaining asymptotic target is exactly

```text
Prove T_L[J]>32000000*exp(-3L/2) for the corrected finite main
```

That arithmetic lower bound is not proved here; neither are
positive-boundary exclusion, `Lambda <= 0`, or RH.
