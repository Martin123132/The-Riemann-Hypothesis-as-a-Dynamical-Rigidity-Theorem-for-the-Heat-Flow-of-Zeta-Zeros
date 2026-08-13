# Endpoint-driven transport for symmetric Fourier pairs

Date: 2026-08-13

Status: exact-lemma certificate; not a proof of the Volterra operator bound

Let

```text
f_x(u)=(A+2u)exp[i*pi*x(A+2u)^2/4],
P_m(x)=I_m(x)+I_-m(x).                                (VT1)
```

The canonical interpolation obeys the exact differential identity

```text
partial_u^2 f_x
 =6 i*pi*x f_x+4 i*pi*x^2 partial_x f_x.              (VT2)
```

Integrating (VT2) against the symmetric Fourier character and using the
integer endpoints gives a first-order transport equation:

```text
4 i*pi*x^2 P_m'(x)
 +(4*pi^2*m^2+6 i*pi*x)P_m(x)=2 Delta f_x',           (VT3)

Delta f_x'=e^(i*pi*x*B^2/4)(2+i*pi*x*B^2)
            -e^(i*pi*x*A^2/4)(2+i*pi*x*A^2).
```

At `x=0`, the interpolation is linear, and direct integration gives
`P_m(0)=0` for every nonzero integer `m`.  The integrating factor for (VT3)
is

```text
mu_m(x)=x^(3/2)exp(i*pi*m^2/x).                        (VT4)
```

Therefore the complete pair coefficient, including all endpoint/Fresnel
cancellation, is exactly

```text
P_m(x)=x^(-3/2)exp(-i*pi*m^2/x)/(2 i*pi)
 integral_0^x z^(-1/2)exp(i*pi*m^2/z)Delta f_z' dz.   (VT5)
```

This is not an asymptotic expansion.  The driver vanishes at `z=0`, so its
integrand is locally `O(z^(1/2))` and the lower endpoint is ordinary.  After
inserting (VT5) into the half-Kummer integral, absolute local integrability
permits the triangular order to be exchanged:

```text
K_m+K_-m=1/(2 i*pi) integral_0^(1/2)
 z^(-1/2)exp(i*pi*m^2/z)Delta f_z' G_m(z)dz,

G_m(z)=integral_z^(1/2) W_t(x)x^(-3/2)
                  exp(-i*pi*m^2/x)dx.                 (VT6)
```

Equations (VT5)--(VT6) are a cancellation-preserving replacement for the
unsafe sum of bare `A` and `B` endpoint exponentials.  The finite-roster
interior has not been discarded: it is encoded exactly by the Volterra
propagator from the endpoint driver.

The geometry is also exact.  The endpoint-driver phase has stationary point
`z_D=2m/D`; it aligns with the outer Morse saddle precisely when

```text
D=2m+t/(pi*m).                                        (VT7)
```

At the saved height, (VT7) reproduces

```text
B lower root: 621.5560036102371463889412039565834444235567629462757775115581974443559,
A roots:      39852.3913862312389618431367868261478442698900177001428220670814997234,
              39936.1086137687610381568632131738521557301099822998571779329185002766.       (VT8)
```

The half-boundary corner occurs at

```text
A/4=39894.25,
sqrt(8t/pi)=159576.9121605730711759784239737527473903434524659738630663703318682632,
A-sqrt(8t/pi)=0.08783942692882402157602624725260965654753402613693362966813173682964028.    (VT9)
```

Thus the interesting midpoint symmetry is real but has a disciplined role:
it gives exact boundary data for (VT3)--(VT6).  It does not replace the
continuous quadratic chirp by a linear interpolation.

The next quantitative problem is now an oscillatory Volterra-operator bound
for (VT6), with separate local charts only where the two phases align in
(VT7).  This formulation automatically keeps the `+m/-m` cancellation that
the tangent B profile lacked.

Pi provenance: every `pi` in (VT1)--(VT9) is differentiated or completed
from the equation-(9) Kummer phase and integer Fourier characters.  No
geometric or fitted constant is introduced.

Proof boundary: exact PDE, pair ODE, Volterra solution, locally justified
triangular exchange, and characteristic geometry only.  No quantitative
Volterra norm, nonlinear B crossing, A-fold splice, complete `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
