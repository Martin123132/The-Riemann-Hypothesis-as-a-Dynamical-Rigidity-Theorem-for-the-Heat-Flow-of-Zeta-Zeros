# Fresnel-to-boundary B-face truncation dictionary

Date: 2026-08-13

Status: exact algebra plus saved-height interval correction; not a complete
B-face estimate

The three direct sign-adapted Fresnel terms and the first three repeated
`z`-boundary currents are two nearby, but not identical, truncations.  For
one paired positive index `m`, exact rational algebra gives

```text
C_0,m+C_1,m+C_2,m = F_0:2,m + D_m,                  (TD1)

D_m=48*x^2*(B^4*x^4+40*B^2*m^2*x^2+80*m^4)
    /[pi^4*(B*x-2m)^5*(B*x+2m)^5].                  (TD2)
```

Therefore a direct Fresnel remainder cannot be attached to the boundary
truncation as though `D_m=0`; instead

```text
|exact-(C_0+C_1+C_2)| <= |R_F,m|+|D_m|.             (TD3)
```

Modes 621 and 622 are removed before summation.  Endpoint-monotone
denominators, a finite Arb sum through `100000`, and an analytic
`m^-6` tail prove uniformly on `|xi|<=70`

```text
sum_(m notin {621,622})|D_m|
 <= [1.473382220896900982381436728578127689315566567946299748218399913058673969628811570087359703502695308e-12 +/- 1.18e-106].             (TD4)
```

After the exact equation-(9) projection factor, window width, and maximum
Kummer weight are restored, this costs only

```text
[2.008623530541314263038626341400090849584077132853930785636956278176743558106204164082269591186308018e-20 +/- 2.06e-114] < 2.1e-20.               (TD5)
```

TD5 is numerically negligible but logically required.  It repairs the
dictionary between the direct Fresnel remainder and the boundary-current
budget without changing the existing `1.491e-6` higher-current ceiling.

Pi provenance: every power of `pi` in TD1--TD5 comes from the exact Fresnel
phase, its integration-by-parts recurrence, or the equation-(9) physical
normalization.  No fitted constant is introduced.

Proof boundary: the nonlocal truncation-dictionary correction on the one
saved-height B window only.  No local 621/622 replacement, outside-window
tail, complete B estimate, A-fold splice, complete `T_upper`, height-uniform
theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
