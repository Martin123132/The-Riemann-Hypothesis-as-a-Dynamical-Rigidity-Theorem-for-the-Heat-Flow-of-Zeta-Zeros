# Shared-contour exceptional-mode majorants

Date: 2026-08-07

Status: rigorous symbolic W2/W3 inequalities with an all-374 finite-roster application; not a complete height-uniform correction theorem and not a proof of RH

## Common-contour lemma

After the reciprocal terms have canceled, integration by parts writes the W2
exceptional difference as a difference of phase integrals.  Interpolate
between the paper quadratic phase and the exact cubic phase.  For `Phi3>=0`,
deform both integrals to `arg(u)=7*pi/12`; for `Phi3<0`, use
`arg(u)=3*pi/4`.  Along the entire homotopy the quadratic decay is at least
`B_* sigma_b r^2`, where

```text
B_* = min(Phi2, Phi2+3*Phi3*N),
sigma_b = 1/2  if Phi3>=0,
sigma_b = 1    if Phi3<0.
```

The phase derivative with respect to the homotopy parameter has magnitude at
most `2*pi*|Phi3|*(3*N*r^2+r^3)`.  Dropping the nonnegative linear and cubic
decay terms and integrating the two Gaussian moments gives

```text
|M_b,exact-M_b,quad|
 <= 2*pi*|Phi3| [3*N*sqrt(pi)/(4*a_b^(3/2)) + 1/(2*a_b^2)],
a_b = 2*pi*B_*sigma_b.                                 (1)
```

For W3 use `arg(u)=pi/6` when `Phi3>=0` and `arg(u)=5*pi/12`
when `Phi3<0`.  Since `sin(2*arg(u))>=1/2`, the same homotopy gives

```text
|M_c,exact-M_c,quad| <= pi*|Phi3|/(pi*Phi2)^2.          (2)
```

Here `pi` is forced by the Fourier normalization `exp(2*pi*i*n*u)`.  It is
not an inserted geometric fitting constant.  Neither (1) nor (2) contains
`delta=ceil(F'(N))-F'(N)` or `eta=1+Phi1`; both remain finite as either gap
tends to zero.

## Certified application

Arb verifies (1) on all 374 W2 modes and (2) on all
165 W3 modes:

```text
minimum W2 slack  >= 8.19247881923105036292636473276466283865476232141554E-5  (chain 36)
minimum W3 slack  >= 8.76288601452327450134485487549190873236719965790528E-5  (chain 2)
maximum W2 bound  <= 3.85998126991787118824804275008923069770211445186728E-3
maximum W3 bound  <= 9.93827660502072829705591177351013867627221056762071E-5
worst W2 bound/actual ratio <= 2.03010770351914576976214770600847211044706438540839E+4
worst W3 bound/actual ratio <= 3.38760464934580613018770388065972673310359447484023E+5
```

The large worst ratios show that these are conservative theorem bounds, not
fitted approximations.  They solve the small-gap singularity problem for W2
and W3, but do not yet bound the W4 zero-mode pair or the coupled generic
endpoint tails.  They imply no recurrence budget, outer Hardy estimate,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
