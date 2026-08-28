# Non-A physical-transform ownership ledger

Date: 2026-08-26

Status: exact physical ownership certified; one joined Kummer transform remains open

Use the inherited physical functional

```text
P_t[F]=2(pi/(32t))^(1/4)
 Re[e^(-i*pi/8) integral_0^1 W_t(x)F(x)dx].          (PT1)
```

For `a=3/4-it/2`, `b=3/4+it/2`, Euler's Kummer integral gives, label by
label,

```text
K_t(alpha)=2(pi/(32t))^(1/4)
 Re[e^(-i*pi/8) alpha B(a,b)
                    1F1(a;3/2;i*pi*alpha^2/4)],

Q_K=sum_(alpha=A,A+2,...,B)K_t(alpha)
   =P_t[sum_(n=0)^L f_x(n)].                          (PT2)
```

Thus the complete source has an exact physical transform: it is the finite
equation-(9) Kummer roster `Q_K`.  It does not have a certified numerical
interval here.  The implemented or published hybrid telemetry cannot be
substituted for it.

The finite carrier ownership now closes without overlap.  The target
`P` block transforms to `G`; modes `39895..39936` transform to `G_extra`;
and the paired `A` block transforms to `A_endpoint`.  The previously
certified identity `A_transition=A_endpoint+G_extra` therefore gives

```text
J_Z=P_t[Z]=Q_K-G-A_transition=R_KGamma-A_transition. (PT3)
```

The remaining native physical subtractions retain their original signs:

```text
R_nonA=J_Z-E_Btr,win-E_outer-I_(A,42).                (PT4)
```

Reassembling the signed local and nonlocal first B currents with every
remaining absolute allowance gives the new two-sided enclosure

```text
E_Btr,win=[-0.0001334644064788904440605731217335137745361029174306820826121883351835273065327100000000000000000000000 +/- 1.49e-6],
-0.000134948<E_Btr,win<-0.00013198.                  (PT5)
```

Together with `|E_outer|<8e-10` and `|I_(A,42)|<0.00364`, exact interval
arithmetic gives the sufficient transformed target

```text
-0.02966668<J_Z<0.02939975,

or, more simply,

|J_Z|<0.02939975.              (PT6)
```

Either statement implies `|R_nonA|<0.0331747039947`.  In terms of the
unextracted Gamma residual, the same sufficient corridor is

```text
-0.0663408<R_KGamma<-0.00727437. (PT7)
```

The sub-`5e-20` Gamma-to-classical correction leaves the same displayed
rounded corridor for `Q_K-T`.  These are sufficient conditions for the
stronger absolute non-A route, not necessary conditions for the original
one-sided upper theorem.

The selected next object is therefore the single joined transform `J_Z`,
not six separate pointwise integrals and not the source-hybrid telemetry.
An admissible next method must evaluate or enclose `J_Z` in its physical
normalization while preserving the source/Gamma/A cancellation.  Direct
event-aware `x` panels remain only a fallback.

Pi provenance: every `pi` in (PT1)--(PT7) is inherited from the equation-(9)
physical normalization, Kummer quadratic phase, Gaussian Abel regulator, or
the already certified Gamma/A/B identities.  No geometric or fitted
occurrence is introduced.

Proof boundary: exact saved-height physical-transform ownership, a derived
two-sided B-trace interval, and sufficient scalar target arithmetic only.
No numerical enclosure of `J_Z`, non-A bound, joined `R_after_A`, `R_Dir`,
`Q_K-T`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
