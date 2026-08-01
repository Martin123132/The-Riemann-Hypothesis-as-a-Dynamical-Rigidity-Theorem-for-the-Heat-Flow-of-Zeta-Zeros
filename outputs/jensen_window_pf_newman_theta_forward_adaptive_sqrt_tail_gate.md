# Newman Theta Forward Adaptive Square-Root Tail Gate

Date: 2026-07-25

Status: exact cofinal omitted-tail theorem. This is not a
retained-separation theorem and not a proof of `Lambda<=0` or RH.

## Adaptive Ordinary Count

For `x>=245`, set

```text
g(x)=x/8+2*log(1+x),
K_h(x)=ceil(sqrt(g(x)+h)), h>=0,
N_(F,h)(x)=K_h(x)-1.
```

The ordinary positive-half-line theta series is retained through
`n=N_F(x)`. Its summands are not individually even, so this proof
uses the direct `L1` and first-moment tail; it does not integrate
individual summands by parts.

## Gaussian Tail

For `p=0,1` and `K>=7`,

```text
E_p(K)<=5*p!*exp(1/80)*pi*K^2*exp(-pi*K^2)
         /[4*theta_K*(1-rho_K)],
theta_K=1-121/(80*pi*K^2),
rho_K=exp(2/K-pi*(2K+1)).
```

The defining ceiling gives

```text
exp(-pi*K_h^2+pi*x/8)
  <=exp(-pi*h)*(1+x)^(-2*pi).
```

The ceiling bound `K_h^2<=2(g+h)+2`, together with
`g+1>42` and `pi-1/42>3`, gives

```text
[1+h/(g+1)]*exp(-pi*h)<=exp(-3h).
```

## Cofinal Envelopes

With

```text
C=20*exp(1/80)*pi/[theta_7*(1-rho_7)],
V(x)=C*x^4*(2g+2)/(1+x)^(2*pi),
D(x)=C*x^3*(x+4)*(2g+2)/(1+x)^(2*pi),
```

both envelopes decrease on `x>=245`. At the left endpoint,

| envelope | upper ratio to `exp(-pi*x/8)` |
|---|---:|
| value | `[0.01873595260145674757868387975376218425555027166839555374884467219544046 +/- 3.63e-72]` |
| derivative | `[0.01904184570515400060037667779055830154951843936910405258556050357822316 +/- 2.30e-72]` |

Both are strictly below `1/50`.

## Theorem

For every `0<=t<=1/5`, `x>=245`, and `h>=0`,

```text
|J_t-J_(N_(F,h),t)^F|<exp(-3h-pi*x/8)/50,
|J_t'-(J_(N_(F,h),t)^F)'|<exp(-3h-pi*x/8)/50.
```

The count is on the saddle scale:

```text
For fixed h, N_(F,h)(x)/sqrt(x/(4*pi))->sqrt(pi/2).
```

This replaces the `x^(3/4)` modular count on the omitted-error
side with only about `sqrt(pi/2)` times the Riemann-Siegel count.
For any `eta>0`, taking
`h>=max(0,log(1/(50*eta))/3)` makes both errors smaller than
`eta*exp(-pi*x/8)`. Thus the approximation can follow a margin
that vanishes as `t->0`; retained first-jet separation remains open.
