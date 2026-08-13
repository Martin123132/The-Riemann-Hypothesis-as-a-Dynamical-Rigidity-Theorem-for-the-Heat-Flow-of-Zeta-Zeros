# Extreme-event second-order normal form

Date: 2026-08-12

Status: complete beta^-4 comparison certified at both extreme event
detunings; not a proof of uniform 399-event propagation

Write the exact amplitude and phase relative to the event Airy chart as

```text
A=1+beta^-2 a_1+beta^-4 a_2+O(beta^-6),
Theta=Theta_0+beta^-2 r_1+beta^-4 r_2+O(beta^-6),

a_1=y/2-3z^2/4,
r_1=-2z^5/15+z^3y/3-zy^2/4,
a_2=13z^4/32-3z^2y/8,
r_2=17z^7/315-2z^5y/15+z^3y^2/12.               (ESN1)
```

The complete second-order comparison multiplier is

```text
1+beta^-2(a_1+i r_1)
 +beta^-4(a_2+i r_2+i a_1 r_1-r_1^2/2).             (ESN2)
```

Its degree-four normal polynomial is integrated signed with moments
`H_0,...,H_4`.  The common event contour and a degree-ten polynomial-
Gaussian far-tail bound give

```text
event 398 (mode 39695): normalized=[3.5143444321317218737021903507411479949951171875000000000000000000000000000000000e-8 +/- 3.06e-14], physical=[1.5820108167918789376926724798977375030517578125000000000000000000000000000000000e-8 +/- 1.38e-14]
event 397 (mode 40093): normalized=[3.5707269002704933313907531555742025375366210937500000000000000000000000000000000e-8 +/- 4.74e-13], physical=[1.6073918460146607856131595326587557792663574218750000000000000000000000000000000e-8 +/- 2.14e-13]
```

Both extreme points are below the prototype targets.  This is a rigorous
two-point endpoint result only.  It does not enclose the continuous detuning
interval or intermediate discrete events.  It does not extend exact height
transport to every corrected cell.  It does not establish complete `T_upper`
or prove `Lambda<=0` or RH.  It does not establish a prize-level conclusion.
