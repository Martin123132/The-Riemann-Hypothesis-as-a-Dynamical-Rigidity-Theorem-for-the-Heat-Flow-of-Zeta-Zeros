# Small-t characteristic Euler--Maclaurin remainder gate

Date: 2026-08-25

Status: interval certificate; finite parameter transport remains open.

## Scope

This gate treats the two length-`496284` affine currents between the
`x=2/5` normalization wall and the first positive third-floor walls.  It does
not yet build the physical parameter-cell cover.

## Exact characteristic phases

At `s=1/3`, after removing an integer multiple of `k`, both currents have

```text
W(p,b,t) = sum_(k=0)^n (p+k) e(t k^2 + (b+2pt)k),
e(u)=exp(2 pi i u),   0 <= t <= 1/(2n).
```

The right pair is `(p,b)=(95747/6,-1/3)`; the left pair is
`(p,b)=(47873/3,-2/3)`.  Their phase derivatives stay in the exact
strips

```text
right: [-1/3, 2080879/2977698],
left:  [-2/3, 544156/1488849].
```

Each strip crosses zero once and has absolute ceiling below `0.7`.  Thus the
apparent small-`t` Mordell singularity has disappeared from the finite current;
only one ordinary quadratic stationary point migrates into each interval.

## Explicit remainder

Writing `q'(x)=b+2t(p+x)`, direct differentiation gives

```text
d^r e(q(x))/dx^r = e(q(x)) E_r(q'(x),t),
H_r(D,T) = r! sum_(j=0)^(floor(r/2))
           (2 pi D)^(r-2j) (2 pi T)^j / (j! (r-2j)!),
|E_r| <= H_r(D,T).
```

Since the affine weight is `p+x`,

```text
integral |d^m((p+x)e(q(x)))/dx^m| dx
 <= (pn+n^2/2) H_m + mn H_(m-1).
```

The Fourier series of the periodic Bernoulli polynomial gives the exact
Euler--Maclaurin remainder coefficient

```text
2 zeta(2M)/(2 pi)^(2M).
```

For `M=48` this proves the uniform bounds

| side | whole-interval EM remainder |
|---|---:|
| right | `< 0.000301311815099289837405899250555307844479102641` |
| left | `< 0.00000327515689685535581025799256593700192752294242` |

Both are below `0.001`.  The production run also evaluates the exact Fresnel
primitive and all endpoint corrections through `B_96`, then overlaps the
result with direct phase-reset sums at `t/t_* = 0, 1/16, 1/4, 1/2, 1` on
both sides.

## Pi provenance

Every `pi` here has one of two explicit origins.  The factors `2 pi` in the
derivatives and Fresnel primitive come from the declared character
`e(u)=exp(2 pi i u)`.  The denominator `(2 pi)^(2M)` comes from the Fourier
series of the periodic Bernoulli polynomial.  No geometric circle or polygon
constant is inserted into the theta-current model.

## Source audit

Kuznetsov's truncated-theta algorithm identifies Euler--Maclaurin summation as
the small-parameter branch when `|tau|<1/n`.  Our whole interval satisfies
`0<=t<=1/(2n)<1/n`.  The paper supplies the algorithmic regime cue only; all
characteristic identities, derivative envelopes, constants, and interval
checks in this gate are derived here.

## Boundary

An exact characteristic reduction at s=1/3, one stationary-crossing inventory, and a uniform order-48 Euler--Maclaurin remainder below 0.001 for each of the two finite child currents on 0<=t<=1/(2n), plus five direct-sum overlaps per side only. The endpoint/Fresnel partial sum is pointwise; no finite t- or s-cell transport theorem, endpoint-complete physical source cover, physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.
