# Characteristic saddle migration and contour geometry

Date: 2026-08-11

Status: exact contour geometry validated; not a proof of the contour-tail
integral bound or the outer saddle join

For fixed `y`, the canonical `z` phase is

```text
phi_a(z)=z^3/3-a*z,       a=lambda+y,
z_+-=+/-sqrt(a).                                           (CG1)
```

At the previous radius `R=4`, a saddle crosses the boundary when

```text
y=16-lambda
 =[10.890056960508164371694080832515335249033682979634559207320929530479738006951859443149906077384149044060698049 +/- 1.69e-104].              (CG2)
```

Because (CG2) lies inside `0<=y<=64`, the measured `|z|>4` contribution is
not purely a nonstationary tail.  It contains the stationary connector for
larger `y`.

Choose `R=9`.  The largest saddle on this normal range is

```text
sqrt(lambda+64)=[8.3132390221556745293520374740162070113233622063143660330145200430796663162629018073581594379976614297567297654 +/- 1.02e-105]<9,      (CG3)
```

with clearance `[0.68676097784432547064796252598379298867663779368563396698547995692033368373709819264184056200233857024327023463 +/- 1.02e-105]>0.68`.

The outgoing rays are

```text
z=9+exp(i*pi/6)s,
z=-9+exp(5i*pi/6)s,       s>=0.                          (CG4)
```

and exact expansion on either ray gives

```text
Im phi_a(z)
 =[81-a]s/2+(9sqrt(3))s^2/2+s^3/3.                      (CG5)
```

Uniformly for `0<=y<=64`,

```text
81-a >= [11.890056960508164371694080832515335249033682979634559207320929530479738006951859443149906077384149044060698049 +/- 1.69e-104],
(81-a)/2 >= [5.9450284802540821858470404162576676245168414898172796036604647652398690034759297215749530386920745220303490243 +/- 8.41e-105]>5.9. (CG6)
```

Thus both rays decay immediately and cubically at infinity.  The integrand is
entire in `z`, so the real tails outside `9` may be deformed to (CG4); the
connector `4<|z|<=9` must instead remain part of the stationary core.

The exact finite-`t` comparison can still be continued on the expanded
rectangle.  The conservative phase and relative-amplitude errors obey

```text
|R_phase| <= [0.059615404731365909993317952480708409615772794066338317235454203957417308721391325503965566365720139396565688900 +/- 4.10e-111] <0.06,
|A/(sqrt(2)/pi)-1|
 <= [1.3088057420087649447106292535067770040433775321372522976933399295668057876675816755547305636517899422616548977e-5 +/- 3.92e-111] <1.4e-5. (CG7)
```

These are admissibility bounds, not a replacement for direct
cancellation-preserving quadrature.

Pi provenance: the ray angles are the cubic Airy steepest-descent angles,
while `lambda` and the finite-`t` amplitude retain the Kummer and integer
Fourier--Poisson normalization already traced.  No geometric fit is used.

Proof boundary: exact contour algebra and source-height interval geometry
only.  No connector quadrature, contour-tail integral constant, `y>64`
partition, complete `T_upper`, `Lambda<=0`, RH, or prize-level conclusion is
proved.
