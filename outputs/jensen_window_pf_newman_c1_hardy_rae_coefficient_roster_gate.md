# Hardy coefficient-to-physical-RAE roster gate

Date: 2026-08-06

Status: exact inversion and rigorous 32-pivot fixture recovery; not a proof and continuous q-cell coverage remains open

## What was proved

For either source branch, write the saved cubic phase as

```text
F(x) = a1*x + a2*x^2 + a3*x^3.
```

Before binary128 rounding the accepted source has

```text
branch +: a2 = xx/2 + 2*c,  a3 =  4*c/3,
branch -: a2 = xx/2 - 2*c,  a3 = -4*c/3.
```

Consequently both branches obey the same exact cancellation

```text
xx = 2*a2 - 3*a3.
```

The source also defines `xx=1/(pc-1)`.  If `r=rae/a>1`, then
`pc=2(r^2+r*sqrt(r^2-1))-1=(r+sqrt(r^2-1))^2`.  Since
`(r+s)(r-s)=1`, this gives the exact inverse

```text
pc = 1 + 1/xx,
r  = (sqrt(pc) + 1/sqrt(pc))/2.
```

Eight rational fixtures check these identities without floating-point arithmetic.

## Where pi comes from

The source does not insert an unexplained decimal.  Its ordinary real constant is
`p=4*atan(1)`.  For the more accurate central scale, PARI supplies `2*pi`, the
code divides that value by `2`, and then computes

```text
a = sqrt(8*t/pi).
```

Thus this is the usual circle constant characterized by `pi=4*atan(1)`, with an
independent high-precision PARI evaluation used by the physical scale.  It is not
an arbitrary circle or polygon chosen to fit the data.

## Independent physical roster

The checkpoint immediately before block `20` records
`RN1=656853`.  The completed block records `MT=209`
and `212` sums.  The accepted loop therefore gives

```text
rae_j = 656853 + 209 + 2
        + 2*(209+1)*(j-1)
      = 657064 + 420*(j-1).
```

The telemetry stores the first `32` values.  From
each of their two coefficient branches independently, the inverse above recovers
the same unique integer pivot.  The recovered roster starts at
`657064` and ends at `670084`.

At 320-bit Arb precision:

* maximum saved-coefficient versus ideal-formula gap:
  `2.23209802317547015246360663786087355339012722652143037e-28`;
* maximum recovered-`rae` error:
  `1.50092581422481707279726872703229150907968275712941097e-26`;
* minimum distance of the source linear modulo choice from a half-integer:
  `0.00526593489108541845375616382250543835190112579921387570`;
* maximum exact difference between the two branch reconstructions of `xx`:
  `7.1794121037062591972291734514340896542281320233263195045E-37`.

The calculation was repeated at 192 and 320 bits; all 64 modulo lifts and all
32 recovered integer pivots were identical.

## Boundary

This gate proves a discrete physical-coordinate adapter for the 32 saved block-20
pivots.  It does **not** yet prove that an interval of physical `rae` values maps
inside any saved q-selector cell, does not cover the remaining 180 pivots of the
block, and does not certify the numerical values returned by the source `PSI` or
`ERF` routines.  Those are separate obligations, so this is not a Riemann
Hypothesis or prize-level conclusion.

The next falsifiable test is to interval-evaluate the physical coefficient curve
between adjacent roster points and compare its four coordinate widths with the
already certified q-cell radii.  A failed containment is evidence that subdivision
or a different coordinate enclosure is required; it must not be relabelled as
coverage.
