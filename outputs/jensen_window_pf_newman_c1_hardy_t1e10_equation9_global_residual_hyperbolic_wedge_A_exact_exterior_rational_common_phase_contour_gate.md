# Rational common-phase contour certificate for the exact A exterior

Date: 2026-08-14

Status: rigorous four-term rational common-phase value and full exact-exterior
enclosure after the certified endpoint replacement; not a complete A theorem

Put `q=log((1-x)/x)`.  The four-term replacement from Section 11.438 gives

```text
I_4=sum_m integral_(q0_m)^infinity G_m(q)exp(i(Phi(q)-Phi_*))dq,
G_m(q)=-exp(3q/4)(b_0+b_1+b_2+b_3)/(2pi i).          (RC1)
```

The compact interval ends at the exact rational splice `q=0.0083`,
where `Phi-Phi_*=[200.2564258576398451275253787789496345252602556158712800 +/- 2.51e-48]`.  Splitting at all 84 exact
cutoffs and summing the active modes before integration, 376
Taylor panels of order 36 give

```text
I_compact=[0.7773861590620306329482711937195171625036705660071043078 +/- 5.31e-11]
          +i [-0.7708849002837852545593855942947741759676302665815269438 +/- 5.31e-11].             (RC2)
```

Every panel integrates its Arb power series exactly.  Its omitted analytic
tail is bounded by Cauchy's estimate on a disk of twice the panel radius;
the total compact Cauchy error is
`[5.3019026921526616371242469296064622885824000070460e-11 +/- 4.58e-56]`.

For the far interval, shift to `q=s+i*0.001`.  No endpoint pole
or logistic singularity lies in the contour strip: the minimum real
cutoff-to-pole margin is `[0.0036920709732345963151411464148382534186363037135412 +/- 3.65e-53]`.
The finite vertical leg gives

```text
I_vertical=[0.039682859566687547015343321218921269429676773844610 +/- 9.91e-14]
           +i [-0.041707058899746012931644950525587451269775557609097 +/- 9.91e-14].                     (RC3)
```

On the horizontal ray, monotonicity of
`R/(1+2R cos(eta)+R^2)` and the exact rational amplitude imply

```text
Im(Phi(s+i eta)-Phi_*) >= [80.190446337452177256162256151936206043443883320530 +/- 1.49e-49],
integral_(q_split)^infinity |G(s+i eta)exp(i(Phi-Phi_*))|ds
 <= [1.3705217217812944276404284412719285254858153320418e-31 +/- 1.30e-81].                          (RC4)
```

Restoring the common phase, paper rotation, and equation-(9) normalization
gives

```text
four-term rational exterior = [-0.002401534018481420957842956095803618141597237555197487995 +/- 3.38e-13],
endpoint replacement error <= [2.282352550043819570691992735073827166694472774698517842e-15 +/- 2.33e-70],
full exact A exterior        = [-0.002401534018481420957842956095803618141597237555197487995 +/- 3.40e-13].      (RC5)
```

An independent checker increases precision and Taylor order and decreases
both panel widths.  Its physical and canonical enclosures overlap the saved
balls and have smaller radii.

Pi provenance: every `pi` in (RC1)--(RC5) comes from the exact endpoint
hierarchy, common outer phase, paper rotation, and equation-(9) normalization.
No fitted constant is introduced.

Proof boundary: the full exact exterior for the finite A-transition roster
39853..39936 at t=10^10 only.  No compact exact-minus-affine transformed-
amplitude value, complete A endpoint block, `R_Dir` estimate, complete
`Q_K-T` or `T_upper`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
