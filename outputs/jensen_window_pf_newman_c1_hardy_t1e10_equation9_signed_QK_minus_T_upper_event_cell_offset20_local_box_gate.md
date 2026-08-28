# Offset-20 local `Q_K-T` sign box

Date: 2026-08-28

Status: rigorous two-configuration local interval certificate; disconnected
from the first height subcell and not a maximal-event-cell or RH theorem

Let

```text
I_20=[10000000019.99987999999984822352416813373565673828125000,
      10000000020.00012000000015177647583186626434326171875000]
    =[10^10+20-0.00012, 10^10+20+0.00012].           (LB1)
```

The event-atlas wall enclosures place all of `I_20` strictly inside the same
fixed-roster open cell as `t=10^10`.  The complete joined contour packet is
therefore evaluated without changing labels or crossing an event wall.

The production packet gives

```text
Q_K-T in [-0.002606182806630386039614677429199218750000000000000000000000000000000000 +/- 2.58e-3],
upper endpoint = -2.725495869526639580726623535156250000000000000000000000e-5 < 0.                    (LB2)
```

An independent configuration changes transition panels, endpoint-remainder
slabs, ordinary expansion order, and transition-jet order.  It gives

```text
Q_K-T in [-0.002606182784802513197064399719238281250000000000000000000000000000000000 +/- 2.59e-3],
upper endpoint = -2.547956682974472641944885253906250000000000000000000000e-5 < 0.                    (LB3)
```

The two `Q_K-T`, complete-packet, and joined-transition boxes overlap.  Thus

```text
Q_K(t)-T(t)<0 for every t in I_20.                    (LB4)
```

## Radius barrier

At the same center, every analytic guard also passes at radius
`0.00013`, but the direct Arb enclosure is

```text
[-0.002606206373457098379731178283691406250000000000000000000000000000000000 +/- 2.78e-3].                   (LB5)
```

It straddles zero.  This is an interval-dependency/width obstruction only; it
does not imply a pointwise sign change.

## Where pi comes from

Here `pi` is the ordinary circle constant, circumference divided by diameter.
It is not fitted or inserted from a plot.  It enters the exact Fourier kernel
`exp(2*pi*i*q*y)` and quadratic contour kernel `exp(-i*pi*y^2)`.  Differentiating
their phases produces the natural height coordinate `p=t/(2*pi)` used in the
stationary-point equations.

## Proof boundary

`I_20` is a second certified local island.  There remains an uncovered gap of
almost 20 height units between the first box around `10^10` and this box.  No
derivative sign on that gap, connected interval from `I_1`, maximal-event-cell
coverage, event-wall handoff, all-height theorem, RH, or Clay-prize conclusion
is asserted.
