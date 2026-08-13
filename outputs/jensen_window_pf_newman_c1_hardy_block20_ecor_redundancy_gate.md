# Hardy block-20 `ecor` redundancy gate

Date: 2026-08-06

Status: exact source algebra; the commented leading cubic phase correction is already present in `gc`; not a proof of the contour estimate or RH

## Source calculation

Write

```text
a = i-phi1,  x = 2 phi2,  b = phi3.
```

For a cubic parent the accepted source forms the saddle approximation

```text
c_src = a/x - 3 b a^2/x^3
```

in both endpoint loops at lines [2427, 2470] and
[2435, 2478].  It then evaluates the full cubic Legendre
phase `g(c)=a c-x c^2/2-b c^3` at that `c` (source locations
[2441, 2486]--[2446, 2491]).  Exact expansion gives

```text
g(c_src) = a^2/(2x) - b a^3/x^3 + 9 b^2 a^4/(2x^5)
           - 27 b^3 a^5/x^7 + 27 b^4 a^6/x^9.          (1)
```

Since `x=2 phi2` and the source phase is `2*pi*g`, the coefficient of
`b^2 a^4` already present in (1) is

```text
9*pi/(32*phi2^5).                                      (2)
```

This is exactly the coefficient in the commented `ecor` assignment at line
2465.  Adding a separate factor
`exp(i*ecor*a^4)` to the current `gc` calculation would therefore double the
coefficient in (2), not restore an absent term.

The source assigns `ecor=0` at line 2466, but
`ecor` has no executable use after that assignment.  Its only executable
occurrences are the declaration and zero assignment.  The nearby statement
about `exp(I*ecor)` is a stale comment relative to this source body.

## Boundary

This exact identity removes the displayed `ecor` coefficient from the list of
plausible missing leading corrections for the accepted evaluator.  It does
not prove that the saddle approximation is exact, control terms beyond the
source `c` model, reconstruct the unavailable Maple implementation, bound the
contour or finite-`ip` remainders, or establish a height-uniform recurrence,
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
