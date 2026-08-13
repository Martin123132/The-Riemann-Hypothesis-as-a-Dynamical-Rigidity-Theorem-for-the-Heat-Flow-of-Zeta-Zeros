# R=9 outgoing Airy-ray quadrature

Date: 2026-08-11

Status: direct saved-height outgoing-ray enclosure validated; not a proof of
a height-uniform contour theorem or the `y>64` saddle join

For `a=lambda+y`, deform the positive and negative real tails to

```text
z_+(s)=9+exp(i*pi/6)s,
z_-(s)=-9+exp(5i*pi/6)s,       s>=0.                    (R9R1)
```

The negative real tail has the opposite outward orientation, so the complete
ray integral is

```text
I_ray= exp(i*pi/6) integral_0^infinity F(z_+(s))ds
      -exp(5i*pi/6) integral_0^infinity F(z_-(s))ds.    (R9R2)
```

Here `F` retains the complete 84-mode kernel inside the `y` integral.  The
phase is combined before exponentiation; no exponentially large closed-form
Fresnel terms are subtracted numerically.

On either ray,

```text
|exp(i phi_(lambda+y))|
 =exp(-[(81-lambda-y)s/2+(9sqrt(3))s^2/2+s^3/3]).       (R9R3)
```

After taking absolute values only for `s>2`, convexity gives the two-ray
envelope

```text
2*84*64*exp(-q(2))/q'(2)
 =[3.59375889562187124835772160686e-18 +/- 7.26e-42] <1e-14,   (R9R4)

q(s)=(81-lambda-64)s/2+(9sqrt(3))s^2/2+s^3/3.          (R9R5)
```

Thirty-two resumable interval panels on `0<=s<=2` give

```text
I_ray,trunc=[-0.144582281914288298389553059137 +/- 7.09e-12]
            +i*[-0.00391902734872929159093178210949 +/- 7.09e-12].        (R9R6)
```

Adding (R9R4) componentwise gives the complete enclosure

```text
I_ray=[-0.144582281914288298389553059137 +/- 7.09e-12]
      +i*[-0.00391902734872929159093178210949 +/- 7.09e-12].     (R9R7)
```

Both components of (R9R7) overlap the independently obtained full-line minus
`R=9` real-axis ball.  This closes the direct contour check without using that
independent value as numerical input to the ray integration.

Pi provenance: the angles in (R9R1) are forced by the cubic Airy phase.  The
remaining pi factors retain the Kummer/Fourier--Poisson normalization; no
circle, polygon, or fitted geometric constant is introduced.

Proof boundary: rigorous direct contour quadrature and analytic `s>2` tail at
one saved height only.  No height-uniform ray theorem, `y>64` endpoint/interior
partition, complete `T_upper`, `Lambda<=0`, RH, or prize-level conclusion is
proved.
