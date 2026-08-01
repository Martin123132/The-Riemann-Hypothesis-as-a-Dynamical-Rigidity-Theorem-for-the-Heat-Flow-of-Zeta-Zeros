# Theta Square-Root to Corrected Riemann-Siegel C1 Transfer

Date: 2026-07-25

Status: exact transfer and quantitative architecture guard.
This is not retained separation and not a proof of `Lambda<=0`
or RH.

## Common Normalization

On

```text
L=log(x/(4*pi))>=50, 0<t<=1/5, tL<=25, h>=0,
Z_t=H_t/A_t,
O_(h,t)=J_(N_(F,h),t)^F/(16*x^4*A_t),
```

the exact real-axis normalizer satisfies

```text
A_0/exp(-pi*x/8)
 =pi^(1/4)*(1+x^2)^(7/8)
  *exp((x*atan(1/x)-1)/4)/32.
```

Elementary bounds on `alpha` give

```text
A_t>=A_0>=exp(-pi*x/8)*x^(7/4)/32,
|(log A_t)'|<L/2.
```

The certified prefactor margin at the coarse threshold is
`[0.3313335154949063760847634430497333607831618194069044056338910099840727 +/- 4.01e-71]`.

## Normalized Ordinary Tail

Converting `J_t=16*x^4*H_t` and differentiating exactly gives

```text
|Z_t-O_(h,t)|<E_F0,
|(Z_t-O_(h,t))'|/L<E_F1,
E_F0=exp(-3h)/(25*x^(23/4)),
E_F1=E_F0*(1/2+(1+4/x)/L).
```

The final factor is below `53/100`. Therefore the direct
ordinary sufficient target is

```text
T_L[O_h]=O_h^2+(O_h'/L)^2
  >2*exp(-6h)/(625*x^(23/2)).
```

Any point satisfying this inequality cannot be a double zero
of `H_t`.

## Corrected-Main Transfer

The two finite mains are not termwise identical. They satisfy

```text
Z_t=O_(h,t)+e_F=J_hat_(N,t)+r_RS,
O_(h,t)-J_hat_(N,t)=r_RS-e_F.
```

Combining the two independently certified remainders yields

```text
|O_(h,t)-J_hat_(N,t)|<2501*exp(-3L/4),
|O_(h,t)'-J_hat_(N,t)'|/L<5001*exp(-3L/4),
scaled first-jet distance<5600*exp(-3L/4).
```

## Architecture Consequence

The transfer is exact, but the certified envelopes have very
different scales:

```text
ordinary amplitude: O(exp(-3h-23L/4)),
corrected-RS amplitude: O(exp(-3L/4)).
```

Thus the current corrected-RS envelope costs about `exp(5L+3h)`
relative to the ordinary one. This is a comparison of proved
upper bounds, not a claim that the actual RS error is that large.
Moreover, the current bulk estimate already has

```text
e_A+e_B<1000*exp(-3L/4).
```

Therefore adding higher classical endpoint `C_k` corrections alone
cannot close the certified five-exponent gap. Such a route would
have to upgrade the bulk heat-flow saddle approximation and the
endpoint expansion together.
A proof must therefore establish either the tiny direct ordinary
first-jet target or the existing corrected-main target
`T_L[J_hat]>32000000*exp(-3L/2)`. The phrase 'equivalent after
transfer' is valid algebraically but not as equality of current
quantitative proof burdens.
