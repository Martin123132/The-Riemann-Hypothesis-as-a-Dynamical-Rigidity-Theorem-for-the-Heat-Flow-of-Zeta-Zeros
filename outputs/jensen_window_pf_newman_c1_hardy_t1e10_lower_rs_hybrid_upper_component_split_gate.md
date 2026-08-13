# Lower-RS versus hybrid-upper residual split

Date: 2026-08-09

Status: finite interval diagnostic validated; not a proof of RH

For each of the fifteen saved heights, Arb reconstructs the exact mathematical
version of the source's visible lower component,

```text
L(t)=2 sum_(n=1)^621 cos(theta(t)-t log n)/sqrt(n)+C_0(t),
```

where `C_0(t)` is the leading Riemann--Siegel endpoint formula used by the
source.  The exact complementary upper target is defined by `U(t)=Z(t)-L(t)`.
The source values are read separately from the last stage-1 `ZP` checkpoint,
the stage-2 `rszsum`, and the final rounded sum.

The residual therefore splits identically as

```text
S-Z = (R_source-L) + (P_source-U) + (S-P_source-R_source).
```

This is an exact identity; the interval calculation only encloses its
components.

The 70- and 110-digit component enclosures overlap on all
15 outputs, and every reconstructed identity interval
contains zero.  The independent checker repeats the calculation at higher
precision.

Lower-component dominance count: 0.
Hybrid-upper dominance count: 15.
Final-addition dominance count: 0.
Outputs whose hybrid-upper discrepancy alone rigorously exceeds `0.005`:
15.

The largest lower, upper, and final-addition discrepancy balls are recorded in
the machine artifact.  This finite split localizes the saved error; it does not
turn `U(t)` into a proved hybrid formula and does not identify which analytic
term inside block zero, blocks 1--19, `H(t)`, or the outer truncations is
responsible.  The next step must refine the dominant component from an exact
identity with explicit remainders rather than fit the observed residual.
