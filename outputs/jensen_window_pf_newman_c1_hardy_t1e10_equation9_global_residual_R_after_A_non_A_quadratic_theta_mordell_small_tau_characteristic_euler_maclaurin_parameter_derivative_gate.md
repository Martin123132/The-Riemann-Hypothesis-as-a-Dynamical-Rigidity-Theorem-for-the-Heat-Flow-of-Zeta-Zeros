# Small-t characteristic Euler--Maclaurin parameter derivatives

Date: 2026-08-25

Status: interval derivative certificate; complete-source transport remains open.

## Exact derivative identities

For

```text
W(p,b,t)=sum_(k=0)^n (p+k)e(t k^2+(b+2pt)k),
```

the exact characteristic derivatives at fixed `s` and fixed `t` are

```text
partial_t W=2 pi i sum (p+k)(k^2+2pk)e(q)
           =2 pi i sum (y^3-p^2 y)e(q),

partial_s W=epsilon[sum e(q)
           +4 pi i(1+t)sum(y^2-py)e(q)],
y=p+k,
epsilon=+1 right, -1 left.                         (PD1)
```

The second identity uses `p_s=epsilon` and `b_s=2epsilon`.  It retains the
amplitude derivative and phase derivative together.

## Uniform Euler--Maclaurin remainders

Applying the Section 11.473 Hermite envelope to each polynomial derivative
through order `2M=160` proves

| side | `partial_t W` remainder | `partial_s W` remainder |
|---|---:|---:|
| right | `< 0.0273239048896268960087141408621391747146844864` | `< 0.000000135854701287693782861949699444781280988081562` |
| left | `< 0.0000145786238874247882894462022274062462656729622` | `< 7.24848685667169314994612865950393207348234625e-11` |

The exact Fresnel moments through degree three and the Bernoulli endpoint
corrections through `B_160` overlap direct weighted 496284-term sums at
`t/t_* = 0,1/16,1/4,1/2,1` on both sides.  The independent checker uses
512 bits, reverse Hermite/Bernoulli/block order, and reset blocks of length
13 instead of 16.

## Transport-scale diagnostic

Using only the pointwise derivative magnitudes, a current variation budget
of `0.001` would suggest half-widths no larger than

```text
minimum t diagnostic: 4.42342986709169656168011248122789584819055316e-24,
minimum s diagnostic: 7.66918032743106123704364481335285968963908344e-19.              (PD2)
```

These are diagnostics, not certified cell widths: they do not include
derivative variation across a cell.  Their purpose is architectural.  A
small value proves that the characteristic current should not be transported
in isolation.  Its complete grouped source derivative must be audited next.
If that derivative remains microscopic, the source and pole-free kernel must
be integrated together before norms instead of enclosed by near-constant
value boxes.

## Pi provenance

Every `pi` in (PD1) and the endpoint formulas comes from differentiating the
declared character `e(v)=exp(2 pi i v)`.  The Euler--Maclaurin remainder
denominator comes from the Fourier series of the periodic Bernoulli
polynomial.  No geometric or fitted occurrence is introduced.

## Boundary

Exact partial_t W and partial_s W identities at s=1/3, uniform order-80 Euler--Maclaurin remainder certificates for both derivatives on 0<=t<=1/(2n), and five direct weighted-sum overlaps per side only. Reported transport half-widths are pointwise scale diagnostics, not certified cells. No differentiated joined Mordell endpoint or outer source assembly, finite t- or s-cover, physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.
