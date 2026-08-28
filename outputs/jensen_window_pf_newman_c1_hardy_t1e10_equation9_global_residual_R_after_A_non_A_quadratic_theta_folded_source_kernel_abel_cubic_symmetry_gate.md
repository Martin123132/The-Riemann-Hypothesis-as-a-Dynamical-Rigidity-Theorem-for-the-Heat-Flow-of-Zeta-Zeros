# Folded source-kernel Abel limit and cubic symmetry

Date: 2026-08-25

Status: exact Abel/symmetry reduction and two rigorous wall derivative rows;
the physical non-A quadrature remains open.

## Abel-limit decomposition

For the extended notch `U={622,...,39936}`, put

```text
C_U,epsilon(s)=sum_(m in Z)[1-chi_U(m)]
               exp(-pi epsilon m^2)exp(-2 pi i m s).             (FSK1)
```

Away from the integer lattice, Poisson summation makes the full theta term
vanish as `epsilon` decreases to zero.  On the complete cell, however, that
theta term is a periodic approximate identity.  Therefore

```text
lim integral_0^1 F_x(s)C_U,epsilon(s)ds
 =[F_x(0)+F_x(1)]/2-sum_(m=622)^39936 hat F_x(m),                 (FSK2)

hat F_x(m)=integral_0^L f_x(u)exp(-2 pi i m u)du.                (FSK3)
```

The endpoint half-sum in (FSK2) is mandatory.  Integrating the pointwise
interior limit through `s=0` would silently discard the theta boundary layer.

With `phi_m(y)=x y^2/4-m y`, `A=159577`, and `B=5122421`, the finite coefficient is

```text
hat F_x(m)=(-1)^m{[exp(i pi phi_m(B))-exp(i pi phi_m(A))]/(i pi x)
                  +(m/x)integral_A^B exp(i pi phi_m(y))dy}.      (FSK4)
```

Equation (FSK4) is written for `x>0`.  At `x=0` its two displayed pieces
must not be bounded separately; use the continuous direct-integral limit.

Thus no microscopic `s`-cell cover is required merely to remove the Abel
regulator: the exact handoff is a joined endpoint term plus 39,315 oscillatory
Fresnel coefficients.  Those coefficients must still be summed and integrated
in `x` without splitting their cancellation.

## Cubic symmetry

The notch length is

```text
Q=39936-622+1=39315=3*13105.
```

For `omega=exp(-2 pi i/3)`, the pointwise interior kernel has the exact factor

```text
C_U,0(1/3+u)
 =-omega exp(-2 pi i 622u)
   [1-exp(-2 pi i 39315u)]/[1-omega exp(-2 pi i u)].             (FSK5)
```

Consequently

```text
C_U,0(1/3)=0,
C_U,0'(1/3)=13105*pi*(-sqrt(3)+3i),
|C_U,0'(1/3)|=2*pi*13105*sqrt(3).                                (FSK6)
```

The factor in (FSK5), rather than a near-constant value box, is the natural
interior quadrature coordinate.

## Complete-source wall audit

Exact-rational phase resets at 352 bits give:

| side | `x` | `|F|` lower | `|partial_t F|` lower | `|partial_s F|` lower | product `s` diagnostic width |
|---|---:|---:|---:|---:|---:|
| right | `992567/2481417` | `2.60588197510086679777881211452933315603886538e9` | `4.48103023211304199751582640749685852977845508e20` | `5.41729548776711291439378709941051862418475131e15` | `2.69071523116409066758880619695133350168029342e-18` |
| left | `992567/2481418` | `5.21730125656717748664863016528256567589944535e9` | `3.59363807035235088498978900192681874835506292e21` | `2.17040719113376564723398780249809864091414253e16` | `1.34392973995780689013292614939131460805635000e-18` |

Both direct 2,481,422-term source balls overlap the independently certified
grouped first-wall source boxes.  At `s=1/3`, the interior-limit product and
its `t` derivative vanish exactly because of (FSK6), but

```text
partial_s(F C)=F C_U,0'(1/3)                                    (FSK7)
```

remains large.  The cubic zero removes pointwise `t` stiffness on the symmetry
line; it does not create a macroscopic two-dimensional value box.

## Decision

Use (FSK2)--(FSK5) to split the next proof into a rigorously retained theta
boundary layer and an interior oscillatory Dirichlet/Fresnel assembly.  Do not
enumerate near-constant source or product boxes.  Preserve the finite 39,315
mode sum and its joined endpoint term before norms, then build an altered-
partition `x`-quadrature certificate.

## Pi provenance

Every `pi` above comes from the declared Gaussian Abel/Fourier kernel or the
original quadratic Kummer character.  The cubic root is selected by the exact
integer identity `39315=3*13105`; no geometric or fitted occurrence is added.

## Boundary

The exact complete-cell Abel-limit identity for the folded source factor, its finite Fresnel-mode handoff, the interior cubic-root factorization at s=1/3, and rigorous complete-source value and first-derivative balls at the two first walls only. The reported variation widths are diagnostics, not certified cells. No uniform complete-source derivative theorem, boundary-layer numerical estimate, finite 39315-mode signed bound, physical x quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.
