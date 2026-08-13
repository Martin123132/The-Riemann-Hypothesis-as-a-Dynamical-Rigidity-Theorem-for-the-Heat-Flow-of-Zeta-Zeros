# Extreme-event corrected normal form

Date: 2026-08-12

Status: rigorous beta^-2 corrected diagnostic at both extreme event
detunings; not a proof of uniform 399-event propagation

At an event `t=tau_m=beta^3-beta*d_m^2`, expansion of the exact common
hyperbolic chart gives

```text
Theta_exact-Theta_0=beta^-2 r_1+O(beta^-4),
A_exact=1+beta^-2 a_1+O(beta^-4),

a_1=y/2-3z^2/4,
r_1=-2z^5/15+z^3y/3-zy^2/4.                         (ECN1)
```

The corrected comparison integrand is therefore

```text
exp(i Theta_0)[1+beta^-2(a_1+i r_1)].                (ECN2)
```

The inner quadratic exponential is integrated with rigorous moments through
`y^2`; real, connector, and lifted-contour Taylor tails are enclosed on the
common `R=14`, `Im z=1`, `|Re z|<=21` contour.  A separate polynomial-
Gaussian estimate controls the corrected far tail.  At the two extreme
event parameters:

```text
event 398 (mode 39695): normalized=[2.5066972175164844316717959252227565514203888596966862678527832031250000000000000e-5 +/- 3.06e-14], physical=[1.1284102022978546396912054848062079148718339638435281813144683837890625000000000e-5 +/- 1.38e-14]
event 397 (mode 40093): normalized=[2.5455545543935088995428014935207272628758801147341728210449218750000000000000000e-5 +/- 4.70e-13], physical=[1.1459021494942538177045201000514573763666703598573803901672363281250000000000000e-5 +/- 2.12e-13]
```

Both recover the prototype targets: `False`.
Both rigorously exceed those targets: `True`.
This tests the first corrected comparison only at the two extreme event
points.  It does not enclose intermediate detunings or certify continuous
height on those cells.  It does not propagate all 399 events.  It does not
establish complete `T_upper` or prove `Lambda<=0` or RH.  It does not
establish a prize-level conclusion.
