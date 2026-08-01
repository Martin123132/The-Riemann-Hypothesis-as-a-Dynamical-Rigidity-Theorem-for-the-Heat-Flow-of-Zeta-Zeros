# Newman Theta Q31 Margin Geometry Audit

Date: 2026-07-24

Status: rigorous finite audit and an unproved cofinal candidate.
This is not an extrapolation of Q31 and not a proof of `Lambda<=0`,
RH, or a Clay-prize result.

## What The Weakest Box Means

The complete cache has `1240` records and
`1566` certified leaves. The weakest reported
tail-relative ratio occurs on

```text
t in [423/3100,227/1550]
x in [121/2,61]
J_7 lower endpoint = [-0.9275539323943026077312401498631990174681850823180591442 +/- 2.65e-56]
J_7 upper endpoint = [0.9514626588416926849323340688868009825318149176819408558 +/- 2.66e-56]
J_7' lower endpoint = [0.001271684214987597582727291505044551671012136278112135435 +/- 4.22e-58]
J_7' upper endpoint = [0.6759502762263201094845827602550445516710121362781121354 +/- 3.56e-56]
ratio lower = [145477302928566486117.5587962758389478448188984993262392 +/- 3.93e-35]
```

The retained value interval crosses zero, while the retained
derivative is strictly positive. The ratio is enormous because the
full-tail derivative bar is around `10^-23`; it is not evidence for
an intrinsic margin of size `10^20`.

## Saddle Coordinate

Set

```text
s=sqrt(x/(4*pi)),  x=4*pi*s^2,
partial_s=4*sqrt(pi*x)*partial_x.
```

This compares a value with the change produced by moving one unit
in saddle count. On the fixed Q31 cover, endpoint enclosures give

```text
min max(|J_7|_lower,4*sqrt(pi*x_low)*|J_7'|_lower)
  = [0.007144288599064160132663990479637952513446064310084250408 +/- 3.10e-58]
weakest box: t=[1/155,
51/3100], x=[135/2,
271/4].
```

That is a finite diagnostic only.

## Why Q31 Cannot Set The Cofinal Scaling

Q31 keeps `N=7`. Even the `kappa=1` absolute-tail count
`ceil((1+x)^(3/4))` runs from 16 to 25 across this slab. Thus Q31
does not sample the adaptive-count regime; it can only anchor it.

## Candidate Gate

```text
s=sqrt(x/(4*pi)); G_N(t,x)=max(|J_N(t,x)|,|partial_s J_N(t,x)|)/A_N(t,x)=max(|J_N|,4*sqrt(pi*x)*|J_N'|)/A_N; E_N=max(E_(0,N),4*sqrt(pi*x)*E_(1,N))/A_N; prove G_N>E_N on cells split at every exact n_*(a,t,x)=k surface for a in {5,9} and at every adaptive-count jump; x=4*pi*k^2 is only the leading marker
```

The amplitude `A_N` must come from an explicit Xi-specific
saddle or corrected Riemann-Siegel theorem. The partition must split
at every exact `n_*(a,t,x)=k` surface for `a=5,9` and every
adaptive-count jump. The simpler `x=4*pi*k^2` locations are only
leading markers. The partition must terminate
after a proved asymptotic region. Until those pieces exist, this is
a precise proof-search target rather than a theorem.
