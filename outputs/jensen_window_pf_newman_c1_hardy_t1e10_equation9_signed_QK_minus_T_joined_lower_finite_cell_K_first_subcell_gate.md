# Joined lower finite-cell K on the first height subcell

Date: 2026-08-28

Status: rigorous lower finite-cell component certificate; the ordinary
complements and complete `K_T` are not yet enclosed

On

```text
I_1=[10^10-10^-4,10^10+10^-4],
```

write the lower finite cell as

```text
V_L(t)=integral_0^L F(y,t)dy,       L=621.5,
K_L=(d_t+i theta'-H'/H)V_L.                         (LK1)
```

Since `d_t F=-i log(y)F`, the complete Hardy operator is inserted before the
reversed 2481423-label endpoint collapse:

```text
K_L=integral_0^L [i(theta'-log(y))-H'/H]F(y,t)dy.   (LK2)
```

For `d_r=p/y+y-q_++r` and `T(f)=(f/d_r)'`, exact linearity gives

```text
T^n[(C-i log(y))y^(-1/2)]
 =C T^n[y^(-1/2)]-i T^n[log(y)y^(-1/2)],
C=i theta'-H'/H.                                   (LK3)
```

The second recurrence contains only monomials
`p^j y^(e/2)(log y)^ell d_r^(-k)` with `ell=0,1`.  Near zero its discarded
integral uses the exact elementary moments

```text
integral_0^S y^(b-1)dy=S^b/b,
integral_0^S y^(b-1)|log y|dy
 =S^b(-log(S)/b+1/b^2),             0<S<=1,
 =S^b(log(S)/b-1/b^2)+2/b^2,        S>=1.          (LK4)
```

On the upper slabs the complete signed complex coefficient of each
denominator power is interval-evaluated before its modulus is taken.  At
order 18, split `300`, and 4096 slabs the direct remainder is

```text
[9.557776483835090311765406221501451061119691683058037573292e-16 +/- 3.80e-74].
```

The direct operator enclosure is

```text
K_L: real [-5.550521035144591024700789451841786798604412850259680048077e-5 +/- 1.03e-8]
     imag [1.575591933754295400329504103333343296818185409077174744126e-5 +/- 3.61e-8].                   (LK5)
```

It overlaps the deliberately wider assembly

```text
V_L'+(i theta'-H'/H)V_L:
 real [-5.550521035144591024700789451841786798604412850259680048077e-5 +/- 4.16e-8]
 imag [1.575591933754295400329504103333343296818185409077174744126e-5 +/- 1.48e-7].
```

After the common Hardy projection and stable positive `H` transport, the
lower finite-cell component contributes

```text
Hardy_t[K_L]/H=[-6.682826256110274698585271835327148437500000000000000000000e-5 +/- 1.81e-7] < 0.             (LK6)
```

The independent checker raises precision from 384 to 448 bits, changes the
endpoint order `18 -> 20`, split `300 -> 280`, and slab count `4096 -> 6144`.
It repeats a changed three-label direct logarithmic-chart quadrature and
requires direct and assembled operator forms to overlap.

The next cancellation-preserving coordinate is

```text
K_(L+O)=(d_t+i theta'-H'/H)
        [V_L-T_lower-T_upper+Gamma_defect].         (LK7)
```

At `L`, the lower complement, finite source roster, and upper complement
partition the complete half-integer lattice.  Their endpoint currents must
therefore be joined before a final norm; (LK6) is not substituted for that
future joined calculation.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
half-integer Fourier spacing, Riemann--Siegel phase, Gamma normalization, or
`p=t/(2*pi)`.  No fitted geometric constant supplies `pi`.

Proof boundary: (LK1)--(LK6) certify only the lower finite-cell contribution
to `K_T` on `I_1`.  They do not enclose either ordinary complementary tail,
the complete lower-plus-ordinary join (LK7), the positive-real tail
derivative, tiny target correction, complete `K_T`, a wider `Q_K-T` sign
interval, the full event cell, wall handoffs, an all-height theorem,
`Lambda<=0`, RH, or a prize-level conclusion.
