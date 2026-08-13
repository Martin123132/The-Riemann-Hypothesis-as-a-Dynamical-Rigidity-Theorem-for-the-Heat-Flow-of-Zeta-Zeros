# Hardy block-20 exceptional-mode recombination pilot

Date: 2026-08-07

Status: rigorous four-branch finite-input identity and quadrature certificate; not a height-uniform bound and not a proof of RH

## Why the reciprocal poles disappear

The factor `2*pi` comes from the Fourier phase `exp(2*pi*i*n*u)` in the
Poisson endpoint rays.  The paper's exceptional quadratic models retain its
global coefficient `Phi2`; they do not use the exact local Taylor coefficient
`Phi2+3*Phi3*x`.  Let `M` denote a single endpoint-ray mode and let `Z_quad`
denote the positive-`Phi1` quadratic zero-mode half-ray.  Direct evaluation of
the paper formulas gives

```text
W2 = M_b,lin(N)-M_b,quad(N),
W3 = conjugate(M_c,lin(0)-M_c,quad(0)),
W4 = Z_quad(0)+i*exp(-2*pi*i*F(N))/(2*pi*F'(N)).       (1)
```

Therefore the upper endpoint pair is exactly

```text
R_b-W2
 = G_b-[M_b,exact(N)-M_b,quad(N)],                    (2)
```

where `G_b` contains only full-minus-linear generic modes.  In the W3 branch,

```text
R_c-W3
 = G_c-conjugate[M_c,exact(0)-M_c,quad(0)].           (3)
```

Equations (2)--(3) contain neither `1/delta`, for
`delta=ceil(F'(N))-F'(N)`, nor `1/eta`, for `eta=1+Phi1`.  In the W4 branch the
lower pair instead splits into the finite differences
`Z_exact-Z_quad` and `R_c-i*exp(-2*pi*i*F(N))/(2*pi*F'(N))`.

## Arb certificate

The exceptional exact modes are integrated directly on the already certified
cubic rays, with explicit polynomial-times-exponential tail balls.  The 70-
and 110-digit enclosures overlap, and (1)--(3) plus the complete correction
identity enclose zero on every pilot branch:

```text
chain   1 (W4): model gap <= 1.35540799423525169891116790343507111787918457651082E-108, recombination gap <= 3.53657251307850457530679300077736115781590342521667E-14
chain   2 (W3): model gap <= 1.08832907926713622818471796546948422713406208308813E-108, recombination gap <= 3.29901692435632227194675092363240764825604856014252E-14
chain   3 (W3): model gap <= 9.10238723219679221266584609046797927747050380946744E-109, recombination gap <= 3.29681105468530302130686360229105957841966301202774E-14
chain  36 (W4): model gap <= 1.85675202920745081348041965437585405169225241384569E-106, recombination gap <= 4.99380188630967080699873950067058103741146624088287E-14
```

Aggregate bounds:

```text
maximum model identity gap       <= 1.85675202920745081348041965437585405169225241384569E-106
maximum recombination gap        <= 4.99380188630967080699873950067058103741146624088287E-14
maximum exceptional-mode tail    <= 7.48577985778680634169871981009110712649984275922310E-60
maximum four-piece triangle/correction ratio <= 7.98464023111844078340892734108799427643761676441522E+2
```

This removes the apparent small-gap singularities at the identity level.  It
does not yet supply a height-uniform majorant for the finite differences, an
all-374 replay, recursive accumulation, outer Hardy control, `Lambda<=0`,
PF-infinity, RH, or a prize-level conclusion.
