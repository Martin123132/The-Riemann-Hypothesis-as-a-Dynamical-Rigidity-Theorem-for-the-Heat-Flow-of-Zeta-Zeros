# Newman Theta Stable-Remainder Arb Quadratic Matrix

Date: 2026-07-24

Status: complete rigorous full-arithmetic compact matrix.
This is not an outer-tail certificate and not a proof of
`Lambda<=0`, RH, or a Clay-prize result.

## Stable Identity

```text
s(u)=erfc(3*sinh(4u))/2=1-omega(u)
delta_n(u)=s(u)*(phi_n(u)-phi_n(-u))=phi_n(u)-b_n(u)
f_N(u)=sum_(n=N+1)^12 phi_n(u)+sum_(n=1)^N delta_n(u)
tau(u)=sum_(n>=13)phi_n(u)
r_N(u)=f_N(u)+tau(u)
```

The finite integrand evaluates the small switch defect directly,
rather than subtracting two independently accumulated full sums.

## Matrix

```text
M_(p,b)=integral_0^(11/5) u^p exp(u^2/5)P_b(1/5,u)du
Q_(N,p,b,q)=integral_0^(11/5) u^p exp(u^2/5)P_b(1/5,u)[D^q f_N(u)]^2du
integral W|D^q r_N|<=sqrt(M_(p,b)Q_(N,p,b,q))+B_q*M_(p,b)
sum_(b=0)^9 binom(9,b) entry_(N,0,b,9-b)
sum_(b=0)^9 binom(9,b) entry_(N,1,b,9-b)
9*sum_(b=0)^8 binom(8,b) entry_(N,0,b,8-b)
```

The cache is append-only, hash-chained, and fsynced after every
Arb integral.

## Progress

Cache rows: `223/223`.
Complete retained counts: `[4, 5, 6, 7, 8, 9, 10]`.
Missing tasks: `0`.

## Certified Compact Budgets

| N | full d0 upper | full d1 upper |
|---:|---:|---:|
| 4 | `[2.311446228803606412062389673565557072551 +/- 2.14e-40]` | `[1.481199334082046836491114484070924167398 +/- 2.02e-40]` |
| 5 | `[3.578912303690706823227733723563769919738e-5 +/- 3.43e-46]` | `[2.377618136103753562674753318243528584736e-5 +/- 4.5e-48]` |
| 6 | `[2.156724837370598223767145199565110175041e-10 +/- 3.93e-51]` | `[1.474981610581005128085444638524413039291e-10 +/- 4.46e-50]` |
| 7 | `[5.775616164915660570751792075824312574621e-16 +/- 1.55e-56]` | `[4.046480799196019817481453212762728080373e-16 +/- 4.33e-56]` |
| 8 | `[7.558991977401040858882533204137291574703e-22 +/- 3.36e-62]` | `[5.406154074721168308587992953598875573732e-22 +/- 3.95e-62]` |
| 9 | `[5.194670455410591619566561472474149887870e-28 +/- 5.84e-69]` | `[3.782308254512628363063829384710824400067e-28 +/- 1.88e-68]` |
| 10 | `[1.982963848396864385776534656519998107370e-34 +/- 4.57e-74]` | `[1.466800142804524406183633608817826322741e-34 +/- 2.47e-74]` |

## Non-Promotion Guard

Completed rows cover the entire arithmetic remainder on
`0<=u<=11/5`. The outer interval `u>11/5`, retained first-jet
lower separation, and transition-cell theorem remain open.
