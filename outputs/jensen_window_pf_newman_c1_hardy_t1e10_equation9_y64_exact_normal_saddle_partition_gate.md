# Exact y=64 normal-saddle partition

Date: 2026-08-11

Status: exact saved-height endpoint/Fresnel/Morse partition validated; not a
proof of uniform estimates for the three resulting integrals

For a transition mode `m`, the complete finite-`t` dependence on the normal
coordinate `y` is

```text
Theta_m(z,y)=constant-(d_m+beta*tanh(z/beta))*y
             +x(z)y^2/(2beta),
x(z)=[1-tanh(z/beta)]/2.                                (Y64P1)
```

It has positive curvature `x(z)/beta` and the unique normal saddle

```text
y_m(z)=beta[d_m+beta*tanh(z/beta)]/x(z).                (Y64P2)
```

At `Y0=64`, define

```text
g_m(z)=d_m+beta*tanh(z/beta)-64x(z)/beta.               (Y64P3)
```

Then `y_m(z)>=64` exactly when `g_m(z)>=0`, and

```text
g_m'(z)=sech(z/beta)^2[1+32/beta^2]>0.                 (Y64P4)
```

Thus every mode has one crossing.  More generally `g_m(z)=c` has the closed
solution

```text
z_m(c)=beta*atanh([c-d_m+32/beta]/[beta+32/beta]).      (Y64P5)
```

No fitted root finder is used.

The normalized Fresnel coordinate at the boundary is

```text
q_m(64,z)=-g_m(z)sqrt(beta/x(z)).                       (Y64P6)
```

Choose the common buffer

```text
delta=[0.06106374974066852196960112275697700327066298956976831494680221728038747586393433087668391952505892023 +/- 9.39e-102].               (Y64P7)
```

It guarantees `|q_m(64,z)|>=4` whenever `|g_m(z)|>=delta` throughout
`|z|<=9`.  The exact disjoint partition is therefore

```text
g_m(z)<=-delta:             nonstationary outer chart,
-delta<g_m(z)<delta:        endpoint/Fresnel chart,
g_m(z)>=delta:              ordinary Morse chart.       (Y64P8)
```

Assign equality to the two outer charts as displayed.  The 84 modewise bands
all lie inside

```text
[-2.300848454368048408275273263077647154811009046913125013081567602452624926431631092799171896572534703 +/- 3.28e-98] <= z <=
[2.303552766850944527056958976058364256501322959784794517154637035423609327524839776767735300413389493 +/- 4.05e-98]                              (Y64P9)
```

and hence well inside `|z|<9`.  The largest exact normal saddle on the R=9
rectangle is

```text
[48698.05236996770386052579173796680146981224638421403087861099339243375037747434592455992317508809410 +/- 3.70e-95] <50000,         (Y64P10)
```

whereas the upper endpoint is at

```text
y_B=(B-C)/sigma=[288706996.6214686250366005408176047861534849591702779769326445954590745496019976345317832168981947148 +/- 8.95e-92]>1e8.            (Y64P11)
```

Its minimum normalized clearance exceeds one million, so no transition-mode
saddle in `|z|<=9` can be assigned accidentally to the upper endpoint chart.
The exact crossing differs from the canonical threshold
`64/(2beta)-d_m` by less than `2e-5` for every mode.

Pi provenance: `beta`, `sigma`, and `d_m` retain the Kummer and integer
Fourier--Poisson normalization of Section 11.314.  The four-width buffer is a
dimensionless Fresnel-coordinate choice, not a geometric fit.

Proof boundary: exact finite-`t` normal square completion and a disjoint
saved-height chart assignment only.  No uniform integral estimate on any
chart, join to all ordinary interior modes, complete `T_upper`, `Lambda<=0`,
RH, or prize-level conclusion is proved.
