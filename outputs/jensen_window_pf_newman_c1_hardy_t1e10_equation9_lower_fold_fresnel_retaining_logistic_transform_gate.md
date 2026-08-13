# Exact Fresnel-retaining logistic-Morse transform

Date: 2026-08-12

Status: exact positive-mode transform certified; this is not a proof of an
Airy/logistic remainder comparison or 399-event propagation theorem

For a positive finite-Poisson mode, define

```text
q_D(m,x)=sqrt(x/2)(D-2m/x),       D in {C,B},
F'(q)=exp(i*pi*q^2/2),            F(0)=0,
```

and keep the endpoint current and Fresnel term paired:

```text
P_m(x)= [exp(i*pi*q_B^2/2)-exp(i*pi*q_C^2/2)]/(i*pi)
       +m sqrt(2/x)[F(q_B)-F(q_C)].                    (FL1)
```

Completing the square in the exact mode integral gives

```text
J_m(x)=(-1)^m exp(-i*pi*m^2/x) P_m(x)/x.               (FL2)
```

No approximation has entered (FL1)--(FL2).

Put

```text
r=t/(2*pi*m^2),       v=(1-x)/(r*x),
x=1/(1+r*v),
s=sgn(v-1)sqrt(2[v-1-log(v)]).                         (FL3)
```

For

```text
psi_m(x)=-pi*m^2/x+(t/2)log((1-x)/x),
```

the phase is globally exact:

```text
psi_m(x)-psi_m(x_m)=-t*s^2/4.                          (FL4)
```

If `K_m` denotes the Kummer-weighted `x` integral of `J_m`, reversing the
monotone `x(s)` orientation gives

```text
K_m=(-1)^m exp(i*psi_m(x_m))
    integral_R exp(-i*t*s^2/4) A_m(s) ds,              (FL5)

A_m(s)=r*x(s)*(dv/ds)*[x(s)(1-x(s))]^(-1/4) P_m(x(s)),
dv/ds=v*s/(v-1),             (dv/ds)|_(s=0)=1.         (FL6)
```

Equations (FL5)--(FL6) are the required Fresnel-retaining logistic chart.
The universal Morse phase is Gaussian, while the endpoint/Fresnel current
remains exact inside the amplitude and is summed before absolute values.

## Event regularity

At `tau_m=pi*m(C-2m)`,

```text
x_m=2m/C,       q_C(m,x_m)=0.                          (FL7)
```

For the nearest event, `m=39894`,

```text
x_m=[0.4999968667163814334145898218414934483039535772699073174705628004035669300713762008309468156438584508 +/- 3.56e-101],
q_C=[9.363323371298775135582115452008733351814599525865910576241412384112335337790759144117852033355661050e-97 +/- 3.72e-96],
q_B=[2481415.224985782684224518129617620840289245227135651492989529376158664479066332836005433747371783561 +/- 2.39e-94],
r*x_m*[x_m(1-x_m)]^(-1/4)
 =[0.7071112123256777031645955082146022816916880477231222365949880387944855147256459786707552174331662772 +/- 1.51e-100].
```

Thus at the crossing, the lower endpoint contributes the finite values
`exp(i*pi*q_C^2/2)=1` and `F(q_C)=0`.  The apparent characteristic defect is

```text
Delta=[-6.266567237133170820356317013103392092845460185365058874399192866139857247598338106368712283098438782e-6 +/- 2.06e-101],
```

but neither (FL1) nor (FL5)--(FL6) divides by it.

## Next comparison

The valid prototype target is now precise: retain mode `39894` exactly on
its `pi/16` buffer, and on the overlap compare the grouped Airy-Fresnel
integral with (FL5) using the same `P_m` current before taking absolute
values.  The full-interior replacement `F(q_B)-F(q_C) -> 1+i` is excluded by
the preceding Fresnel-retention gate.

## Boundary

This is an exact coordinate and current identity.  It does not bound the
variation of `A_m`, prove an Airy-to-logistic remainder, propagate through
399 events, establish complete `T_upper`, prove `Lambda<=0`, RH, or a
prize-level conclusion.
