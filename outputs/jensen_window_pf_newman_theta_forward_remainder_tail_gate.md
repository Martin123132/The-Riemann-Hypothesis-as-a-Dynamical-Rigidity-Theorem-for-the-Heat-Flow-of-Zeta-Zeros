# Newman Theta Forward-Remainder Tail Gate

Date: 2026-07-24

Status: exact stable-remainder and compact arithmetic-tail theorem.
This is not a complete derivative budget and not a proof of
`Lambda<=0`, RH, or a Clay-prize result.

## Stable Remainder

The exact modular partition gives

```text
r_N=Phi-sum_(n<=N)b_n
   =sum_(n>=1)phi_n-sum_(n<=N)b_n.
delta_n=phi_n-b_n
       =(erfc(3*sinh(4u))/2)*(phi_n(u)-phi_n(-u)).
delta_n(u)=(erfc(3*sinh(4u))/2)*(phi_n(u)-phi_n(-u)); r_N^(q)(u)=D^q[sum_(n=N+1)^12 phi_n(u)+sum_(n=1)^N delta_n(u)]+tau_q(u), tau_q(u)=sum_(n>=13)D^q phi_n(u), 1<=N<=12
```

This representation preserves modular cancellation while leaving
only the rapidly convergent forward theta tail to bound.

## Derivative Recurrence

With `X=pi*n^2*exp(4u)`,

```text
D^q phi_n=pi*n^2*exp(5u)*P_q(X)*exp(-X)
P_0=2X-3
P_(q+1)=(5-4X)P_q+4X P_q'
```

For `A_q` equal to the coefficient 1-norm of `P_q`, `n>=13`,
and `u>=0`,

```text
|D^q phi_n(u)|
 <=A_q*pi^(q+2)*n^(2q+4)*exp(-pi*n^2).
```

## Explicit Tail Constants

| q | A_q | log10(B_q) |
|---:|---:|---:|
| 0 | 5 | `-224.430560640319554134782614384` |
| 1 | 53 | `-220.680218197746976501171604594` |
| 2 | 661 | `-216.859256030554317917393746384` |
| 3 | 9445 | `-212.979218950474298371287584553` |
| 4 | 151269 | `-209.049634399525175370365652954` |
| 5 | 2672213 | `-205.077476685661153682202861769` |
| 6 | 51450613 | `-201.067920627348149193750603841` |
| 7 | 1069789829 | `-197.024976138966201909388726545` |
| 8 | 23844232389 | `-192.951854680324361201803837325` |
| 9 | 566257358709 | `-188.851187590021686145095154645` |

Here `sup_(u>=0)|tau_q(u)|<=B_q` for the omitted forward
series `n>=13`.

## Compact Matrix Correction

For every matrix weight on `I=[0,11/5]`,

```text
M_(p,b)=integral_I W_(p,b)
Q_f=integral_I W_(p,b)f^2
sqrt(M_(p,b)Q_f)+B_q M_(p,b)
```

The remaining obligation is an Arb matrix for the stable finite
remainder plus an explicit `u>11/5` tail. No retained first-jet
lower separation or Newman conclusion is claimed.
