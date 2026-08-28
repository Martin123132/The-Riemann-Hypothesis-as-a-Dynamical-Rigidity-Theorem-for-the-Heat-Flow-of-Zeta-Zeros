# Exact Kummer-to-Hardy truncation bridge target

Date: 2026-08-27

Status: exact target isolated; not a proof of the bridge value

Use the corrected half-domain normalization from the repair gate and define

```text
U=Z-L,
rho_RS=U-T,
Delta_KU=Q_K-U.                                      (KB1)
```

Here `U` is the independently enclosed complementary Hardy upper component,
`T` is the exact classical upper main, and `Q_K` is the corrected finite
equation-(9) Kummer roster.  These definitions give the exact identities

```text
Q_K-T=Delta_KU+rho_RS,
R_KGamma=Delta_KU+rho_RS-(G-T),
J_Z=Delta_KU+rho_RS-(G-T)-A_transition.               (KB2)
```

At `t=10^10`, Arb gives

```text
rho_RS=[-1.3072691049744271751888591526308660287074286085960356843173854653718142590663242700912616358030642389000000000e-9 +/- 2.41e-97].                        (KB3)
```

Transporting the already certified sufficient `Q_K-T` corridor through
(KB2) yields the rigorous displayed target

```text
-0.0663407986927308950<Delta_KU<-0.0072743686927308951    (KB4)
```

The full outward-rounded endpoints retained in the artifact are

```text
lower: [-0.06634079869273089502557282481114084736913397129257139140396431568261453462818574093367572990873812380 +/- 3.07e-102]
upper: [-0.007274368692730895025572824811140847369133971292571391403964315682614534628185740933675729908738604597 +/- 6.36e-104]
```

Thus the stronger route does not ask us to prove `Q_K` is almost equal to
`U`.  It requires a definite negative finite-roster bias: even `Delta_KU=0`
misses the upper wall by more than `[0.0072743686927308950255728248111408473691339712925713914039643156826145346281857409 +/- 3.37e-83]`.

The saved hybrid value is deliberately not substituted for `Q_K`.  As a
route diagnostic only, its hybrid-minus-exact-upper ball lies above the
required upper wall by at least
`[0.00029135199153005446802741961708662959307117850104262984456587947166456944600688036 +/- 3.09e-84]`.  This reinforces the need
for an exact truncation/continuation theorem rather than telemetry transfer.

The direct fixed-precision Arb call to `1F1` is also rejected as a production
evaluator.  All three corrected-label scouts contain zero with astronomical
enclosures; the lower endpoint already has

```text
[+/- 6.02e+6821866891].
```

This is a diagnostic of the naive representation, not an impossibility
theorem for scaled Kummer functions, steepest descent, or a common contour.
Standard contiguous relations shift Kummer parameters, whereas the odd-label
roster keeps `(a,b)` fixed and changes
`z_alpha=i*pi*alpha^2/4`; moreover
`z_(alpha+2)-z_alpha=i*pi(alpha+1)`.  The exponential phase repeats on the
odd lattice, but `1F1(a;b;z)` is not periodic in `z`.

The selected analytic object is therefore one exact continuation bridge for
`Delta_KU`: lower omitted labels, the analytically continued upper tail,
the exact `H(t)` normalization, and any lower transition term must be
assembled with signs intact before taking absolute values.

Proof boundary: exact saved-height identities and sufficient interval target,
plus fixed-precision conditioning diagnostics only.  No enclosure of
`Delta_KU`, `Q_K`, or `J_Z`, no non-A bound, all-height theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
