# Offset-20 complete `K_T` monotonicity bridge

Date: 2026-08-28

Status: rigorous two-configuration local derivative certificate and leftward
sign extension; not a proof of the all-height theorem or RH

Let

```text
D_20=[10^10+20-0.006,10^10+20+0.006].              (MB1)
```

The complete ownership identity is

```text
K_T=K_(L+O)+(K_tr+K_U)+K_tail+K_corr,
(Q_K-T)'=Hardy_t[K_T]/H.                            (MB2)
```

All five components are evaluated on the same Arb height box before the final
projection.  The production component sum gives

```text
(Q_K-T)' in [0.002545273240684764460052793108960924216727489649125110160450140000780610 +/- 2.50e-3],
lower endpoint = [5.238088851911015840942152692967421672748964912511016045e-5 +/- 1.41e-61] > 0.                 (MB3)
```

An independent configuration reverses the ordinary expansion path and changes
the duplication, transition, endpoint, and upper-arc settings.  It gives

```text
(Q_K-T)' in [0.002545273240684764460052793310909315953306511834527622873000490000218029 +/- 2.50e-3],
lower endpoint = [4.598662781063467296535282018665970330651183452762287300e-5 +/- 4.91e-61] > 0.                 (MB4)
```

The direct and component-sum projections overlap within both configurations;
the production and independent direct projections, component projections, and
complete complex `K_T` boxes also overlap.  Hence `Q_K-T` is strictly increasing
throughout `D_20`.

## Sign extension

The earlier value gate proves `Q_K-T<0` on
`I_20=[10^10+20-0.00012,10^10+20+0.00012]`.  Put
`a=10^10+20+0.00012`.  For every `t` between the left endpoint of `D_20` and
`a`, the fundamental theorem of calculus gives

```text
Q_K(t)-T(t)
 =Q_K(a)-T(a)-integral_t^a (Q_K-T)'(u) du
 <=Q_K(a)-T(a)<0.                                   (MB5)
```

Therefore

```text
Q_K(t)-T(t)<0 for every
t in [10000000019.994,
      10000000020.00012].            (MB6)
```

This extends the offset-20 sign island leftward by `0.00588` height units.  The
remaining open gap to the right endpoint `10000000000.0001` of the first
subcell is `19.9939` height units.

## Radius diagnostic

Production remains positive at radius
`0.0061`, but at radius
`0.01` its direct and
component enclosures both straddle zero.  This is an interval-enclosure width
obstruction only; it proves no pointwise derivative zero or sign change.

## Where pi comes from

Here `pi` is the ordinary circle constant, circumference divided by diameter.
It enters through the exact Fourier kernel `exp(2*pi*i*q*y)`, the quadratic
kernel `exp(-i*pi*y^2)`, Gamma normalization, and the height coordinate
`p=t/(2*pi)`.  It is not fitted from a plot, polygon, or numerical pattern.

## Proof boundary

This gate proves `(Q_K-T)'>0` only on `D_20` and `Q_K-T<0` only on the
left-anchored interval (MB6).  It does not bridge the remaining gap, certify
the whole event cell, cross an event wall, prove an all-height theorem,
`Lambda<=0`, PF-infinity, RH, or a Clay-prize conclusion.
