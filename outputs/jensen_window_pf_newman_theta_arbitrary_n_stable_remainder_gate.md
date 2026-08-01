# Newman Theta Arbitrary-N Stable Remainder Gate

Date: 2026-07-24

Status: exact arbitrary-`N` forward-remainder theorem.
This is not the cofinal retained-separation theorem and not a proof
of `Lambda<=0`, RH, or a Clay-prize result.

## Exact Split

For every integer `N>=1`,

```text
r_N=sum_(n>N)b_n
   =sum_(n>N)phi_n+sum_(n<=N)delta_n
delta_n=(erfc(3*sinh(4u))/2)*(phi_n(u)-phi_n(-u)).
```

For every `M>=N`, differentiation through order nine gives

```text
delta_n(u)=phi_n(u)-b_n(u)=(erfc(3*sinh(4u))/2)*(phi_n(u)-phi_n(-u)); r_N^(q)(u)=sum_(n=1)^N D^q delta_n(u)+sum_(n=N+1)^M D^q phi_n(u)+T_(M+1,q)(u), T_(K,q)=sum_(n>=K)D^q phi_n, N>=1, M>=N, K=M+1
```

Taking `M=N` removes the former fixed cap:

```text
r_N^(q)(u)=sum_(n=1)^N D^q delta_n(u)+T_(N+1,q)(u)
```

## Uniform Tail

With `X=pi*n^2*exp(4u)`,

```text
D^q phi_n=pi*n^2*exp(5u)*P_q(X)*exp(-X)
P_0=2X-3
P_(q+1)=(5-4X)P_q+4X P_q'.
```

If `A_q` is the coefficient 1-norm of `P_q`, then for every
`K>=2`, `n>=K`, `u>=0`, and `0<=q<=9`,

```text
|D^q phi_n(u)|
 <=A_q*pi^(q+2)*n^(2q+4)*exp(-pi*n^2).
rho_(q,K)=exp((2q+4)/K-pi*(2K+1))<1.
B_(q,K)=A_q*pi^(q+2)*K^(2q+4)*exp(-pi*K^2)/(1-exp((2q+4)/K-pi*(2K+1)))
sup_(u>=0)|T_(K,q)(u)|<=B_(q,K).
```

The two uniform endpoint checks reduce exactly to
`16*pi-45>48-45>0` and `11-5*pi<11-15=-4`.

## Witness Scale

| K | q | log10(B_(q,K)) |
|---:|---:|---:|
| 2 | 9 | `18.3907534327205064106618448609` |
| 3 | 9 | `15.4389430648130260743924124925` |
| 4 | 9 | `8.63696060630065286403137083898` |
| 8 | 9 | `-50.2304444735559946975241769770` |
| 13 | 9 | `-188.851187590021686145095154645` |
| 33 | 9 | `-1435.17688019795261023756419028` |
| 65 | 9 | `-5707.38433867713777784213579508` |

The table is only a scale check. The theorem is the symbolic
all-`K` inequality above.

## Cofinal Boundary

For an adaptive retained count, the forward arithmetic tail is
now explicit as `B_(q,N+1)`. The unresolved terms are the full
finite switch-defect derivative sum and, more importantly, a
uniform retained value-or-derivative lower margin. Those are
separate gates and are not inferred from Q31 or this theorem.
