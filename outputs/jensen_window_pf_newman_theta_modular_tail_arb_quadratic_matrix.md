# Newman Theta Modular-Tail Arb Quadratic Matrix

Date: 2026-07-24

Status: complete rigorous finite compact matrix.
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

Cache rows: `223/223`.
Complete retained counts: `[4, 5, 6, 7, 8, 9, 10]`.
Missing tasks: `0`.

## Certified Compact Budgets

| N | d0 upper | d1 upper |
|---:|---:|---:|
| 4 | `[2.311446228803606412062423872378569481940 +/- 3.54e-40]` | `[1.481199334082046836491260234937470906161 +/- 3.44e-40]` |
| 5 | `[3.578912303690782413909208426422165234438e-5 +/- 3.50e-45]` | `[2.377618136103799725793656934898333934474e-5 +/- 2.47e-46]` |
| 6 | `[2.156724869008086962213757515243563440621e-10 +/- 3.25e-50]` | `[1.474982029958602760829431976673771326220e-10 +/- 4.70e-50]` |
| 7 | `[5.879536661916348072733220123296501393868e-16 +/- 3.54e-56]` | `[5.003701572871956160161867949879959353078e-16 +/- 1.35e-56]` |
| 8 | `[7.448670578364535386802462777645216158452e-18 +/- 1.97e-58]` | `[2.287819193327009786997766725938843647199e-17 +/- 1.33e-57]` |
| 9 | `[2.116147365805095197074245093326076197656e-17 +/- 1.04e-57]` | `[7.598243383666986441557266853224076565498e-17 +/- 3.13e-57]` |
| 10 | `[7.876421143304176606333708440141244118827e-19 +/- 2.25e-59]` | `[3.571073894852620912304873161001730472769e-18 +/- 1.24e-58]` |

## Non-Promotion Guard

Every row covers only `n=N+1..24` and `0<=u<=11/5`.
The arithmetic tail `n>=25`, outer interval `u>11/5`, retained
first-jet lower balls, and transition-cell theorem remain open.
