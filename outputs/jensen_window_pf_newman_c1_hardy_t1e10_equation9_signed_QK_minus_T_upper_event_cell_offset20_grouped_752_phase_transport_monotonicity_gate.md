# Offset-20 grouped-752 phase-transport monotonicity gate

Date: 2026-08-28

Status: rigorous two-configuration local derivative and sign-extension
certificate; not a proof of the all-height theorem or RH

Let `t_c=10^10+20`, `h=t-t_c`, and

```text
D_20^phase=[t_c-0.07,t_c+0.07].                     (PT1)
```

The old direct grouped-collar integrator evaluated 192 height-dependent panel
phases separately.  This forgot their common rotation.  On panel `j`, choose
`s_j` as the midpoint of its log-height interval and put
`s_0=(log L+log C)/2`.  For the exact grouped integrand,

```text
F_t(y)=F_(t_c)(y) exp(-i h log y)
      =exp(-i h s_0) exp(-i h(s_j-s_0))F_(t_c)(y)+E_j. (PT2)
```

If `epsilon_j=max|log y-s_j|` and
`B_j=752*(y_right-y_left)/sqrt(y_left)`, then the finite 752-term geometric
sum gives

```text
|E_j|<=|h| epsilon_j B_j,                           (PT3)
```

because `|exp(-i x)-exp(-i y)|<=|x-y|` for real `x,y`.  The log-weighted
integral has the companion bound `log(C)|h|epsilon_j B_j`.  The common Hardy
phase is enclosed without losing `h` by

```text
theta(t_c)+h(theta'(D_20^phase)-s_0),               (PT4)
```

which follows from the real mean-value theorem.  Thus (PT2)--(PT4) are a
rigorous representation change, not sampled interpolation.

The 192-panel, 384-bit production configuration gives

```text
(Q_K-T)' in [0.002437108602634907583837059793662644647767919620084953180594311561448800 +/- 2.18e-3],
lower endpoint = [0.0002613876520172198071522987462993633977679196200849531806 +/- 5.70e-60] > 0.                  (PT5)
```

A 256-panel, 448-bit independent configuration integrates panels in reverse
order and changes the underlying lower-ordinary, transition, endpoint, and arc
settings.  It gives

```text
(Q_K-T)' in [0.002437108602634907583837059793662644647767919620084953180594311561464550 +/- 2.18e-3],
lower endpoint = [0.0002575619025719788625283516027446758977679196200849531806 +/- 5.70e-60] > 0.                  (PT6)
```

The derivative, lower-ordinary, and complete complex `K_T` enclosures overlap.
Each phase-transport result also overlaps its deliberately broader direct
interval evaluation.  Hence `(Q_K-T)'>0` throughout `D_20^phase`.

The prior value gate gives `Q_K-T<0` at
`a=10^10+20+0.00012`.  Therefore

```text
Q_K(t)-T(t)<0 for every
t in [10000000019.93,
      10000000020.00012].            (PT7)
```

The exact remaining gap to the first subcell is
`19.9299` height units.  At radius
`0.1`, the production
phase-transport enclosure

```text
[0.002322161277561463066395034147191846236874110051551234334965556627101350 +/- 3.18e-3]         (PT8)
```

straddles zero.  This is an enclosure-width diagnostic only and proves no
pointwise derivative zero or sign change.

## Where pi comes from

`pi` is the ordinary circumference-to-diameter constant.  It enters through
the exact Fourier kernel `exp(2*pi*i*q*y)`, quadratic kernel
`exp(-i*pi*y^2)`, Gamma normalization, and `p=t/(2*pi)`.  It is not fitted
from a plot, polygon, or numerical pattern.

## Proof boundary

This gate proves the complete derivative only on `D_20^phase` and negativity
only on (PT7).  It does not bridge the remaining gap, cover the event cell,
cross an event wall, prove an all-height theorem, `Lambda<=0`, PF-infinity,
RH, or a Clay-prize conclusion.
