# Classical upper-main versus hybrid-ZP split

Date: 2026-08-09

Status: finite interval diagnostic validated; not a proof of RH

The exact complementary upper target from Section 11.297 is split once more:

```text
U(t)=M_upper(t)+rho_RS(t),
M_upper(t)=2 sum_(n=622)^39894 cos(theta(t)-t log n)/sqrt(n).
```

Here `rho_RS=U-M_upper` is defined exactly.  It contains all classical
correction and remainder content left after the source-style leading endpoint
term already included in the lower component.  The source hybrid residual
therefore satisfies the exact identity

```text
P_source-U=(P_source-M_upper)-rho_RS.
```

The 70- and 110-digit Arb enclosures overlap for all fifteen upper main sums
and remainder balls; the independent checker repeats them at higher precision.
All 15 source
hybrid-minus-classical-main discrepancies exceed `0.005`.

Maximum exact remaining correction/remainder:

```text
[1.3072692374137841682761465547411376966304627619316165121701537784213069382284574418877819314609368924252976213e-9 +/- 9.10e-96]
```

Minimum source-hybrid versus classical-main discrepancy:

```text
[0.0066363171935595117189810598088943535341522286888172335760663213272673584147518827962222923550608457701116390500 +/- 3.59e-97]
```

The minimum discrepancy-to-remainder ratio is

```text
[5076473.1576552409116810744202439145693397769399638233410684402275279173536456900466783288364854960783938349687 +/- 3.53e-80]
```

Thus the tiny exact Riemann--Siegel correction/remainder left outside the
classical upper main sum cannot explain the saved `ZP` failure.  The error is
inside the source's approximation to the upper main sum itself.  This finite
localization does not identify which hybrid block, transition, normalization,
or asymptotic replacement creates it and has no RH implication.
