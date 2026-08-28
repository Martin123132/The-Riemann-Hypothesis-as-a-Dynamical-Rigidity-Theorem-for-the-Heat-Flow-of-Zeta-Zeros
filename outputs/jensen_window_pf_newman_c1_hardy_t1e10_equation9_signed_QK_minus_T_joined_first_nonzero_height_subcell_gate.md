# First nonzero joined Q_K-T height subcell

Date: 2026-08-28

Status: rigorous first nonzero-radius height subcell with `Q_K-T<0`
certified; extension to the maximal event cell remains open.

The exact A-free identity from Section 11.508 is

```text
Q_K(t)-T(t)=Hardy_t[Y_T(t)]/H(t),
Y_T=S_W-H D_T.                                      (HS1)
```

On the fixed-roster interval

```text
I_1=[10^10-10^(-4),10^10+10^(-4)],                 (HS2)
```

the same packet is rebuilt before projection as

```text
Y_T=U_unowned+O_join+I_transition
    -C_G D_(39853..39894)+(C_G-H)D_(622..39894).    (HS3)
```

Every scalar height occurrence in the lower cell, upper quarter arc,
ordinary complementary packet, transition packet, finite Dirichlet blocks,
phase, and prefactor is evaluated with Arb on the full box.  The Section
11.507 event atlas proves that every roster and endpoint used by (HS3) is
unchanged on `I_1`.

Two large-expression dependency losses are removed without approximation.
First,

```text
log H(t)-log H(t_0)=integral_(t_0)^t H'(u)/H(u) du. (HS4)
```

The direct and duplication formulas for `H'/H` overlap on the box, giving the
stable prefactor enclosure

```text
H(I_1)=[1.000000000000000000000312499999999999841801075230765990 +/- 3.51e-28].                        (HS5)
```

Second, the upper-arc action has one maximizer `delta_*(t)`.  The exact
envelope identity

```text
d A_max(t)/dt=delta_*(t)                            (HS6)
```

transports its point value across the box and gives

```text
A_max(I_1)=[0.01979967368031121930651718904437006522489619007936360954 +/- 1.44e-8].            (HS7)
```

All original polynomial, omitted-numerator, compact-tail, and positive-tail
error guards remain active.  Representative component radii are

```text
U_unowned: real [8.940865063777891919016838e-5 +/- 7.38e-31],
           imag [0.0005090024442324647679924965 +/- 9.53e-30],
O_join:    real [4.075355967358973430236802e-6 +/- 4.95e-31],
           imag [4.018557639540176751324907e-6 +/- 5.48e-32],
transition join:
           real [0.0006036017166479723528027534 +/- 4.85e-29],
           imag [0.0006331650974971125833690166 +/- 4.74e-29].
```

After the single final Hardy projection and division by the positive `H` box,

```text
Q_K-T=[-0.002620053190184989944100379943847656250000000000000000000 +/- 2.58e-3],
lower=-0.005192062995774904265999794006347656250000,
upper=-4.804338459507562220096588134765625000000e-5<0.                               (HS8)
```

The strict upper-endpoint negativity margin is
`4.00531862396746873855590820312500000000000000e-5`.  The independent checker covers the
same interval by two adjacent half-boxes at 448 bits and changes the upper-arc
configuration, grouped panel counts, transition jet order, ordinary-tail
recurrence order, slab counts, and clustering power.

Pi provenance: every `pi` in (HS1)--(HS8) comes from the inherited
Riemann-Siegel phase, Gaussian/Fresnel kernel, Fourier characters, gamma
normalization, or the already-certified event equations.  No fitted geometric
constant is introduced.

Proof boundary: a rigorous fixed-roster sign theorem only on the compact
interval (HS2).  This does not extend the sign to the rest of the maximal
event cell, cross either event wall, establish an all-height theorem, prove
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
