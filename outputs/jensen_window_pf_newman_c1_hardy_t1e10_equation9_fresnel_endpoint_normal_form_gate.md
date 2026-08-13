# Exact Fresnel endpoint normal form for the portcullis modes

Date: 2026-08-10

Status: exact endpoint-uniform alpha normal form validated; not a proof of the quantitative x-stationary remainder

For each positive Poisson mode, complete the square in the alpha integral and
introduce

```text
z_C(m,x)=sqrt(x/2)(C-2m/x), C in {A,B},
F(z)=exp(i*pi/4)/sqrt(2)
     erf(exp(-i*pi/4)*sqrt(pi/2)z).
```

Then `F'(z)=exp(i*pi*z^2/2)` and the mode is exactly

```text
F_m(A,B;x)=exp(-i*pi*m^2/x)*sqrt(2/x)
            [F(z_B)-F(z_A)],                            (FN1)

J_m=(-1)^m/x{[E_m(B)-E_m(A)]/(i*pi)+m F_m(A,B;x)}.    (FN2)
```

No asymptotic expansion has been used in (FN1)--(FN2).  At the reduced
Kummer saddle

```text
x_m=2*pi*m^2/(t+2*pi*m^2),
z_C*=m sqrt(pi/(t+2*pi*m^2))[C-2m-t/(pi*m)].             (FN3)
```

The remaining x phase is

```text
psi_m(x)=-pi*m^2/x+(t/2)log((1-x)/x),
psi_m''(x_m)=-(2*pi*m**2 + t)**4/(8*pi**2*m**4*t)<0.               (FN4)
```

Thus the alpha integral is a smooth Fresnel endpoint transition, while the
reduced x saddle remains nondegenerate.  These two facts do not by themselves
exclude a characteristic boundary tangent in the coupled x/Fresnel integral;
that composite question must be tested separately.

At `t=10^10`, all 84 modes `39853..39936` lie in one certified compact chart:

```text
[0.00022750894296694990809559692404036293914315357882116410965099627379429790500158717305687866190481159006741483041 +/- 1.11e-106] <= z_A* <=
[0.043918284430471208236568666879227136537014420731701157675047926192404798109785259085750350552102576016421131307 +/- 1.11e-106],

z_B* >= [2480138.8143857504635233249542825428297102819312099114926173913418685905387959935446594103017278016880437874277 +/- 5.05e-104].
```

Using `F(+infinity)=(1+i)/2` and
`|F(+infinity)-F(z)|<=2/(pi*z)` for `z>0`, every one of these
84 factors satisfies

```text
|[F(z_B*)-F(z_A*)]-(1+i)/2|
 <= [0.043918541117624634868984033259106670070103648831783030548299706563710949369340990910746913440142704777456777788 +/- 1.11e-106] < 0.044. (FN5)
```

This explains the 42 classical lower endpoint terms: they are nearly
half-saddle Fresnel contributions, not failed or divergent saddles.  The exact
positive-mode partition is

```text
1..621             lower nonstationary B-tail,
622..39852         lower interior,
39853..39894       lower A-endpoint transition,
39895..39936       reflected A-endpoint transition,
39937..2560588     upper interior,
2560589..infinity  upper nonstationary B-tail.
```

The zero and negative modes have no positive-alpha stationary point.  The
next theorem must estimate the x integral in these six positive ranges plus
the nonpositive modes, keeping the boundary and Fresnel pieces in (FN2)
grouped.

Proof boundary: exact endpoint-normal-form algebra and a certified finite
transition atlas at `t=10^10`.  No explicit x-stationary remainder,
height-uniform source error, `Lambda<=0`, RH, or prize-level conclusion is
proved.
