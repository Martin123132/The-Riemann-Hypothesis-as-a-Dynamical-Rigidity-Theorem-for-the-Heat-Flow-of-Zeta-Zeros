# Signed saved-height one-sided residual closure gate

Date: 2026-08-28

Status: interval certificate for the original saved-height one-sided residual;
not an all-height proof, not a proof of `Lambda<=0`, and not an RH proof.

Two already certified exact identities give

```text
R_after_A=I_(A,42)+R_nonA,
R_nonA=J_Z-E_Btr,win-E_outer-I_(A,42).
```

Collecting coefficients before any interval norm cancels the translated
A-face channel exactly:

```text
R_after_A=J_Z-E_Btr,win-E_outer.                    (SC1)
```

The same identity follows independently from

```text
R_KGamma=J_Z+A_transition,
R_KGamma=E_Btr,win+E_outer+A_transition+R_after_A. (SC2)
```

Thus both `I_(A,42)` and `A_transition` have coefficient zero in the final
saved-height `R_after_A` formula.  No sign or quadrature estimate for the
translated A-face is needed for this closure.

The signed B trace is rebuilt from the raw local principal, signed nonlocal
first current, and every remaining certified absolute allowance:

```text
E_Btr,win=[-0.000133464406478890444060573121733513774536102917430682082612188335183527307 +/- 1.49e-6],
-0.000134948<E_Btr,win<-0.00013198.                (SC3)
```

Together with

```text
J_Z=[0.0334441839950159192085266113281250000000000000000000000000000000000000000 +/- 9.69e-7],
|E_outer|< [6.56548602237922141651026664754486292307352108262691963788397295997613088e-10 +/- 4.50e-82],
```

outward-rounded Arb arithmetic gives

```text
R_after_A=[0.0335776484014948096525871844498585137745361029174306820826121883351835273 +/- 2.46e-6].            (SC4)
```

Its upper endpoint is below both `0.0368147039947` and
`0.0368061039947`.  The corresponding rigorous margins are
`[0.00323460341064151104184801422505125671420291075444431791738781166481647269 +/- 6.94e-18]` and
`[0.00322600341064151354331371546848396107349589476323338041738781166481647269 +/- 6.94e-18]`.

Restoring the independently certified projector-completed A transition gives

```text
R_Dir=[-0.00309647639627883706537890664906199823786932306198112566050784988570576608 +/- 2.46e-6].                    (SC5)
```

This is below zero, hence below both inherited `R_Dir` sufficient targets.
The alternate construction `R_KGamma-E_Btr,win-E_outer` overlaps (SC5).

Finally, the sub-`5e-20` Gamma-to-classical correction gives

```text
Q_K-T=[-0.00322994080275772750943947977079551201240542597941180774312003822088929339 +/- 9.69e-7].              (SC6)
```

Its upper endpoint is negative with margin
`[0.00322897280275224181756717080874841194526675410441180774312003822088929339 +/- 3e-79]` and therefore is also
below the working `8.6e-6` threshold.  Rejoining the B trace, outer ball,
`R_Dir`, and Gamma correction with dependency loss still remains negative;
its margin is `[0.00322600443459119253354079483701388817207339472941180774312003822088929339 +/- 3.71e-75]`.

Pi provenance: this gate introduces no new occurrence of `pi`.  All inherited
occurrences belong to the already certified Fresnel, Fourier, Kummer, Gamma,
or Riemann--Siegel normalizations.

Proof boundary: this certifies the original signed one-sided residual and
`Q_K-T<0` only at the single saved height `t=10^10`, for the exact finite
roster and common-regulator identities pinned here.  It does not prove the
stronger absolute non-A bound, an all-height finite split-contour transport
theorem, equation-(4) as a global infinite-series identity, `Lambda<=0`,
PF-infinity, RH, or a prize-level conclusion.
