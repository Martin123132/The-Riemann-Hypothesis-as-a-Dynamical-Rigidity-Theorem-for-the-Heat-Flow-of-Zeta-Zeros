# First nonzero fixed-selector lower-fold height cell

Date: 2026-08-11

Status: rigorous local height-cell enclosure validated; not a proof of a
selector-transition join or a complete source-height theorem

Put

```text
S(t)=2pi sum_m G_64(d_m;lambda(t)),
lambda(t)=[pi*C^2/8-t]/beta,
C=159577.
```

On the selector-stable cell

```text
|t-10^10| <= 0.01,
[9999999999.990000000000000000000000000000000000000000000000000000000000 +/- 3.06e-61] <= t <= [10000000000.01000000000000000000000000000000000000000000000000000000000 +/- 1.62e-61],
```

the Airy arguments lie in

```text
[-69.10994768107896592736337592278572333035237947644572864317321132594660 +/- 1.84e-64]
    <= -lambda(t)-y <=
[-5.109938397904705329248462412183606171580254564285152942184929613014258 +/- 1.84e-64],
0 <= y <= 64.
```

An Arb interval subdivision into 2048 covering panels
proves `|Ai(x)|<0.54` throughout that interval.  Since

```text
G_tt=(-lambda*G-i*G')/beta^2,
```

the integral representation gives the uniform grouped bound

```text
|S''(t)| < 0.146.
```

Taylor's theorem therefore proves, for every real height in the cell,

```text
|S(t)-S(t0)-(t-t0)S'(t0)| < 0.0000073,
|S(t)-S(t0)|                 < 0.001313.
```

This is the first certified nonzero-radius height theorem for the canonical
lower-fold lattice.  The radius is intentionally small: it establishes the
local mechanism and a reusable checker before any large subdivision campaign.

Proof boundary: this certificate covers only `t=10^10 +/- 0.01` with the odd
selector fixed.  It does not cross a selector fold, join adjacent rosters,
or cover the source interval.  It does not prove `T_upper`, `Lambda<=0`, RH,
or a prize-level result.
