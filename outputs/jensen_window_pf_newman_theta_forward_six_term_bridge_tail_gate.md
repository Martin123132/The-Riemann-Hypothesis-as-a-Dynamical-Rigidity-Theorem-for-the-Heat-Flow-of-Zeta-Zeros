# Newman Theta Forward Six-Term Bridge Tail Gate

Date: 2026-07-25

Status: exact finite-bridge omitted-tail theorem. This is not a
retained-separation theorem and not a proof of `Lambda<=0` or RH.

## Why Ordinary Terms Work Here

For `u>=0`,

```text
Phi(u)=sum_(n>=1)phi_n(u),
phi_n(u)>0.
```

The ordinary summands are not even. Accordingly, this gate does not
reuse the ninefold modular integration-by-parts estimate. It bounds
the raw positive-half-line `L1` tail and its first moment directly.

Define

```text
S_(6,t)^F(x)
 =sum_(n=1)^6 integral_0^infinity
   exp(tu^2)phi_n(u)cos(xu)du.
```

## Explicit Tail

For `K>=6`, `p in {0,1}`, and

```text
E_p(K)=sum_(n>=K) integral_0^infinity
       u^p exp(u^2/5)phi_n(u)du,
```

the directed majorant is

```text
E_p(K)<=5*p!*exp(1/80)*pi*K^2*exp(-pi*K^2)/[4*theta_K*(1-rho_K)]; theta_K=1-121/(80*pi*K^2), rho_K=exp(2/K-pi*(2K+1)), p in {0,1}
```

This uses the coefficient norm five of `P_0(X)=2X-3`. No
cancellation or zero information enters the bound.

## Direct C1 Contract

With `K=7`,

```text
|J-J_6^F|<=16*x^4*E_0(7),
|J'-(J_6^F)'|<=64*x^3*E_0(7)+16*x^4*E_1(7).
```

After division by `exp(-pi*x/8)`, every positive multiplier is
increasing because

```text
d/dx log(exp(pi*x/8)*x^s)=pi/8+s/x>0.
```

Thus `x=245` is the worst endpoint throughout `38<=x<=245`.

## Directed Witnesses

| retained | first omitted | J error / gamma at 245 | J' error / gamma at 245 |
|---:|---:|---:|---:|
| 5 | 6 | `[388060.7167916033030525796948738226863094375369078980252670927776314124 +/- 3.85e-65]` | `[394396.4019637111120820095674431912199634691701635371767000249046131498 +/- 1.57e-65]` |
| 6 | 7 | `[9.645876016252800903561127821835024777240708998165622402398980787313357e-13 +/- 1.98e-83]` | `[9.803359706314071122394778888313963957277292002217306033458555983840922e-13 +/- 2.13e-83]` |
| 7 | 8 | `[4.302009913278438859342572536685387837218011277150707285729765403704792e-33 +/- 2.83e-103]` | `[4.372246809821760310107349231161883965172591053104188221006986063357115e-33 +/- 4.22e-103]` |

Five retained terms are an explicit failure guard. Six terms
put both ratios below `10^-11`; seven terms give a much smaller
cross-check.

## Theorem

For every `0<=t<=1/5` and `38<=x<=245`,

```text
|J_t-J_(6,t)^F|<10^-11*exp(-pi*x/8),
|J_t'-(J_(6,t)^F)'|<10^-11*exp(-pi*x/8).
```

The next finite task is therefore a retained six-term interval
cover on

```text
[1/1035,1/155] x [38,69]
union [1/1035,1/5] x [69,245].
```

That retained cover and the cofinal retained theorem remain open.
