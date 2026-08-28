# Certified affine A-corner wedge values

Date: 2026-08-13

Status: rigorous canonical scalar and affine wedge evaluation at modes 39894
and 39895; not a curved-face, amplitude-remainder, or global residual bound

For a quadratic phase `A R^2+B R`, complete the phase before deformation.
The lower-indented full-line integral is

```text
L_full(A,B)=i pi[1+erf(e^(i sgn(A)pi/4)B/(2sqrt(|A|)))]. (CE1)
```

Subtract its right ray
`R=r0+e^(i sgn(A)pi/4)y`.  Along that ray the quadratic modulus is exactly
`exp(-|A|y^2+Ly)`, so the omitted tail has an explicit monotone Gaussian
majorant.  Combining the rays according to (NU3)--(NU4) retains the Abel
half-jump and avoids epsilon extrapolation.

The certified scalar wedges are

```text
C_39894=[0.55049140553459663112859950915520916939514170453882768431637604447162677336343298 +/- 2.98e-32]
       +i [0.51629103497222081174590329766744384157723201692609454873505773010637774928159094 +/- 2.98e-32],

C_39895=[0.22614941108099700601176040873911234153214197291793806380663987155676996657476770 +/- 4.78e-33]
       +i [0.15430715853068163107687694878365332420574747590804558243996880432469501004838343 +/- 4.78e-33]. (CE2)
```

The exact affine moments from Section 11.428 give

```text
C_aff=C+i mu_S partial_h C-i(lambda_P+mu_S rho)partial_a C. (CE3)
```

Hence

```text
C_aff,39894=[0.55115455318147109192380669333650930558785340212997628799345999364156127281076714 +/- 2.98e-32]
           +i [0.51595349243531182980598371367242479463406528694263694460257419305072874245833444 +/- 2.98e-32],

C_aff,39895=[0.22614030418561043001638213324418544153500580727428556237440107288193002408853150 +/- 4.78e-33]
           +i [0.15430902651539819474224065397253725318119507144361064529800653354874726698181585 +/- 4.78e-33]. (CE4)
```

At mode 39894 the mandatory affine correction is

```text
[0.00066314764687446079520718418130013619271169759114860367708394916993449944733415945 +/- 5.96e-32]
+i [-0.00033754253690898193991958399501904694316672998345760413248353705564900682325649128 +/- 5.96e-32]. (CE5)
```

Its real component exceeds `6e-4`, so a scalar-only corner carrier would not
be an admissible approximation at the target precision.  These are values of
the canonical tangent wedge only.  The exact paired-triangle phase carrier,
curved-face strip, transformed-amplitude remainder, all other modes, and the
global projector orientation have not been summed here.

Pi provenance: (CE1) is the exact Gaussian rotation of the equation-(9)
hyperbolic phase; (CE3) uses the Fourier Gaussian moments.  No fitted
constant or geometric surrogate is introduced.

Proof boundary: rigorous canonical scalar/affine wedge values and explicit
steepest-ray tail bounds at two saved-height modes only.  No exact curved-face
or transformed-amplitude remainder, `R_Dir` estimate, complete `Q_K-T` or
`T_upper`, all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
