# Full exact A-exterior common-phase reduction

Date: 2026-08-14

Status: exact one-dimensional reduction and certified cutoff geometry, with
nonrigorous scale telemetry; not an exterior value or complete A bound

For

```text
F_m(z)=pi m^2/z+pi A^2 z/4,
h_A(z)=z^(-1/2)(2+i pi A^2 z),
H_m(x)=integral_0^x h_A(z)exp(iF_m(z))dz,             (XR1)
```

put `alpha=exp(-i*pi/4)sqrt(pi)A/2`,
`beta=exp(-i*pi/4)sqrt(pi)m`, `X=sqrt(x)`, and
`w_+/-=beta/X +/- alpha X`.  Differentiating the exact inverse-Gaussian
erfc primitive gives

```text
H_m(x)=4sqrt(x)exp(iF_m(x))
 -2sqrt(pi)beta(-1)^(mA)[erfc(w_-)+erfc(w_+)].        (XR2)
```

The phase-stripped endpoint tail `B_m=exp(-iF_m)H_m` obeys

```text
B_m'+iF_m'B_m=h_A.                                   (XR3)
```

Since every exterior cutoff lies before the endpoint saddle, define

```text
b_0=h_A/(iF_m'),
b_(j+1)=-b_j'/(iF_m').                               (XR4)
```

This gives an exact finite hierarchy plus one explicit oscillatory
remainder.  At the 84 cutoff points, the maxima of `|b_0|,...,|b_3|` are

```text
[764.2294305113607894710054242156947998475897810123903792780329832548107 +/- 1.16e-68],
[0.01122818975505151250696597487708410083543408437918844421666611210235773 +/- 1.85e-72],
[4.935333096449782921147980581214820564433178567820281495326604650545145e-7 +/- 1.16e-77],
[3.615523237262133425186793862670614095851727838502636699436997331055288e-11 +/- 3.47e-81].       (XR5)
```

The endpoint phase cancels the mode-dependent outer term exactly.  Hence the
full convergent exterior is the single common-phase family

```text
C_A,ext^exact=-1/(2pi i) sum_(m=39853)^39936
 integral_0^(x0_m) x^(-7/4)(1-x)^(-1/4)B_m(x)
 exp(i Phi_A(x))dx,                                  (XR6)

Phi_A(x)=pi A^2 x/4+(t/2)log((1-x)/x),
x0_m={1+[t/(2pi m^2)](1+0.0037)}^(-1).             (XR7)
```

Its unique stationary point is
`x_*=[0.499475380364729741276539059975136114155171359502937676758769515653551522212 +/- 5.80e-77]`.  The signed Morse endpoint
`w0=sgn(x0-x_*)sqrt(2[Phi_A(x0)-Phi_A(x_*)])` runs from
`[-10.55088896860266921802726608598352738766947446953811197267676806252211 +/- 3.26e-69]` to `[1.095334032285792364311797927904263153061057949912153416127995196554505 +/- 4.98e-70]` and is positive exactly for
modes `39927..39936`.

Diagnostic only: a 180-node-per-mode Gauss pilot on `w>=-10`, a three-term
endpoint hierarchy, and 5 outer boundary terms
give physical partial estimates

```text
-0.0023881597576624815847689516466433381
-0.0024016189542409348569034552089469742
-0.0024015622018266558833286391168445659
-0.0024015335665506138029390894152305718
-0.0024015338135894309700554353233689293       (XR8)
```

The final displayed scale, `-0.0024015338135894309700554353233689293`, is not an
interval enclosure.  It only shows that the exact exterior is likely a
signed carrier component at a much larger scale than the final `R_Dir`
budget.  Promotion requires interval bounds for the endpoint-hierarchy
remainder, stationary-region quadrature, and the final outer-IBP remainder.

Pi provenance: every `pi` in (XR1)--(XR8) comes from the original
equation-(9) triangle, exact endpoint/outer phases, and paper normalization.
No fitted constant is introduced.

Proof boundary: exact erfc primitive, phase-stripped ODE, one-dimensional
common-phase reduction, cutoff/Morse roster, and pointwise cutoff hierarchy
only.  The diagnostic pilot is not a certificate.  No exact exterior value,
uniform hierarchy remainder, outer-tail bound, compact nonlinear-amplitude
bound, complete A theorem, `R_Dir` estimate, RH, or prize-level conclusion is
proved.
