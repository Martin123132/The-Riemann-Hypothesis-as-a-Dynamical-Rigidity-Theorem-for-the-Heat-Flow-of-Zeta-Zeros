# Newman Theta Modular-Tail Arb Quadratic Refinement

Date: 2026-07-24

Status: partial resumable rigorous finite compact matrix.
This is not an infinite-tail certificate and not a proof of
`Lambda<=0`, RH, or a Clay-prize result.

## Matrix

For `m=9`, weighted Cauchy-Schwarz reduces the compact pieces of
`d0` and `d1` to 29 analytic quadratic integrals per retained `N`.
The cache is append-only, hash-chained, and fsynced after every
integral.

```text
M_(p,b)=integral_0^(11/5) u^p exp(u^2/5)P_b(1/5,u)du
Q_(N,p,b,q)=integral_0^(11/5) u^p exp(u^2/5)P_b(1/5,u)[D^q sum_(n=N+1)^24 b_n(u)]^2du
sum_(b=0)^9 binom(9,b) sqrt(M_(0,b)Q_(N,0,b,9-b))
sum_(b=0)^9 binom(9,b) sqrt(M_(1,b)Q_(N,1,b,9-b))
9*sum_(b=0)^8 binom(8,b) sqrt(M_(0,b)Q_(N,0,b,8-b))
```

## Progress

Cache rows: `79/136`.
Complete retained counts: `[7, 8]`.
Missing tasks: `57`.

## Certified Compact Budgets

| N | d0 upper | d1 upper |
|---:|---:|---:|
| 7 | `[5.775616164915660570751792075510261002765e-16 +/- 3.92e-56]` | `[4.046480799196019817481453212589053038097e-16 +/- 7.36e-57]` |
| 8 | `[7.558991977401040858882490304203121677054e-22 +/- 4.84e-62]` | `[5.406154074721168308587968382225149559683e-22 +/- 3.12e-62]` |

## Non-Promotion Guard

Every row covers only `n=N+1..24` and `0<=u<=11/5`.
The arithmetic tail `n>=25`, outer interval `u>11/5`, retained
first-jet lower balls, and transition-cell theorem remain open.
