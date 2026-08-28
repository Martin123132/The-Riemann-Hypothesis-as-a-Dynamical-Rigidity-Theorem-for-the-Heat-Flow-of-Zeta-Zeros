# Rationally closed 240-label upper-complement collar

Date: 2026-08-27

Status: the least rationally closed endpoint collar containing all 230 upper
saddles is rigorously enclosed; both remaining upper tails are nonstationary;
not a proof of the complete joined packet or RH.

At `L=621.5`, every even consecutive half-integer roster has zero grouped
amplitude.  At `B=L+1/16=621.5625`, the Fourier step has exact order 16.
Therefore the least multiple of 16 containing the 230 stationary labels is

```text
240=15*16.
```

The extended collar has exact zeros at both endpoints:

```text
sum_(j=0)^239 (-1)^j=0,
sum_(j=0)^239 exp(2*pi*i*j*621.5625)=0.            (RC1)
```

It contains the 230 stationary labels and ten additional nonstationary labels.
The grouped collar is

```text
P_240=sum_(q=q_+)^(q_++239)
      integral_L^B y^(-s)exp(-i*pi*y^2+2*pi*i*q*y)dy.           (RC2)
```

Thirty-two Arb panels at 384 bits give

```text
P_240 = [0.00030782796494112583871038043632907191040644603185258356848647300656075359934231000 +/- 9.76e-57]
      + i [0.00028016069584300095386042085295876417097966110951237517727265227761890702061954502 +/- 9.76e-57],

|P_240| = [0.00041623079114240153776910839305849149486077044551163656236875082522620543223326049 +/- 1.38e-56].
```

The exact remainder split is

```text
T_upper=P_240+R_B+R_L,                            (RC3)

R_B=sum_(q=q_+)^(q_++239) integral_B^a (...),
R_L=sum_(q>=q_++240) integral_L^a (...).
```

Every phase in both remainders is nonstationary.  The minimum `R_B` gap at
`B` is `[27.756031598419002372836641297342834233954904672518891224228313374799248240762753 +/- 3.23e-79]`.  The minimum
`R_L` launch gap at `L` is
`[10.320323486150138636148471208977723899505382840728114765180066020955201503860175 +/- 1.61e-79]`, replacing the previous
near-endpoint gap `0.320323...` by a gap greater than ten.

Pi provenance: every `pi` comes from the inherited Fresnel phase, exact
half-integer Fourier spacing, or `p=t/(2*pi)`.  The numbers `1/16`, `16`, and
`240` come exactly from the rational endpoint and its root-of-unity order; no
fitted geometric constant is used.

Proof boundary: exact dual-endpoint cancellation, least closed collar length,
production gap census, and rigorous grouped complex enclosure of `P_240` only.
No quantitative enclosure of `R_B`, `R_L`, the lower complementary tail,
complete ordinary or joined packet, `J_Z`, or `D_K` is proved, nor is any
non-A, all-height, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
