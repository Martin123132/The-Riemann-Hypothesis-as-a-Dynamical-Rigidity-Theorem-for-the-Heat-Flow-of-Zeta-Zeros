# Newman Zero-Slab Degree-71 Sector Certificate

Date: 2026-07-25

Status: rigorous interval certificate and published-theorem composition.
This proves a bounded-degree Newman/Jensen theorem. It is not a proof of
PF-infinity, `Lambda<=0`, RH, or a Clay-prize conclusion.

```text
work/rh_compute/results/jensen_window_pf_newman_zero_slab_degree71_sector_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_zero_slab_degree71_sector_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_zero_slab_degree71_sector_certificate.py
```

Current result:

```text
validated Newman zero-slab degree-71 sector certificate: 13 rows, 0 issues, 4 directed integrals, 1 positive slab margin, 1 zero-free complex slab, 3 published inputs, all shifts through degree 71, 0 all-degree theorems
```

## Uniform Complex Slab

For

```text
H_t(z)=integral_0^infinity exp(t*u^2)*Phi(u)*cos(z*u)du
0<=t<=1/5, |Re z|<=84/5, |Im z|<=1,
```

put `a=5*pi/168`. Every kernel summand is positive because

```text
Phi_n(u)=pi*n^2*exp(5u)*(2*pi*n^2*exp(4u)-3)
         *exp(-pi*n^2*exp(4u))>0.
```

On `0<=u<=a`,

```text
cos(xu)>=cos((84/5)u)>=0, exp(tu^2)>=1, cosh(yu)>=1.
```

On `u>=a`, the real integrand is bounded below by
`-exp(u^2/5)*Phi(u)*cosh(u)`. The first summand is integrated directly.
For `n>=2`, consecutive `n^4*exp(-pi*n^2*exp(4u))` terms have ratio at
most `(3/2)^4*exp(-5*pi)<1`, giving

```text
sum_(n>=2) Phi_n(u)
 <=32*pi^2*exp(9u-4*pi*exp(4u))
   /(1-(3/2)^4*exp(-5*pi)).
```

The 192-bit directed values are

```text
C_R      = [0.0250137178322222592781720466155290738354880163433099023243796069742164755707699343706888691 +/- 5.98e-55]
T_1      = [0.0250013637245988698624665825257790936376379193863078792549355807003372801796452099498413696 +/- 6.62e-58]
T_ge2    = [1.26239731337832931418386021631717220860018844644733550153248964299810354830020759218822325e-7 +/- 4.24e-59]
T_inf    = [2.81229244156541662833207459697003328528302017357678380448767119539977295514227473578212149e-2037 +/- 7.91e-2091]
C_R-T_1-T_ge2-T_inf
          = [1.22278678920515827740457037283484806292369381573783358938692646700682598895904006843418812e-5 +/- 5.99e-55]
```

For the final infinite tail, `cosh(u)<=exp(u)` and

```text
u^2/5+10u<=(pi/2)*exp(4u), u>=2,
```

reduce the integral to the displayed explicit
`exp(-(pi/2)*exp(8))` bound. Therefore

```text
Re H_t(x+iy)>0
```

throughout the stated slab, so the slab is zero-free.

## Sector Transfer

The de Bruijn strip-contraction theorem gives

```text
H_t(z)=0 => |Im z|<=sqrt(1-2t)<=1
```

on the same heat interval. The zero-free slab therefore forces every zero
to satisfy `|Re z|>84/5`. For

```text
F_t(s)=sum_(k>=0) A_k(t)*s^k/k!=2*H_t(i*sqrt(s)),
```

every zero has the form `-z^2`, hence lies in the negative-axis sector

```text
delta=2*atan(5/84),
sin(delta)=840/7081,
|sin(delta)|^(-2)=50140561/705600.
```

The double-exponential kernel gives `F_t` order `1/2`, and all its Taylor
coefficients are positive. Chasse's sector theorem and closure of the sector
class under differentiation now give

```text
P_(d,n,t)(x)=sum_(j=0)^d binom(d,j)*A_(n+j)(t)*x^j
is hyperbolic with real negative zeros
for every 0<=d<=71, n>=0, and 0<=t<=1/5.
```

Primary theorem inputs:

- D. H. J. Polymath, Theorem 3.2:
  <https://arxiv.org/abs/1904.12438>
- Matthew Chasse, Theorem 3.6:
  <https://doi.org/10.1080/17476933.2011.584250>
- Griffin--Ono--Rolen--Thorner--Tripp--Wagner, endpoint comparison:
  <https://arxiv.org/abs/1910.01227>

## Quartic Consequence

Degrees four and five are now closed for every shift throughout the complete
positive Newman target interval. In particular, an actual Xi quartic window
cannot realize the strict outward contact encoded by the exact length-14
survivor, and its adjacent quintic cannot have the survivor's nonreal pair.
The survivor remains a valid finite countermodel to generic signed-Hankel
promotion; the new theorem excludes it only by Xi kernel geometry plus the
published sector machinery.

At `t=0`, Corollary 1.3 of the effective Xi paper already gives the much
larger published range `d<=9.36*10^20` for every shift. The new contribution
here is the uniform continuum `0<=t<=1/5`, not a stronger endpoint cutoff.

## Proof Boundary

The cutoff `71` is finite. Nothing here proves degree `72`, an unbounded
cofinal terminal sequence, all-degree Jensen hyperbolicity, PF-infinity,
`Lambda<=0`, RH, or a Clay-prize conclusion. The next Jensen-side problem is
to enlarge the complex zero-free slab or obtain an unbounded sector sequence;
the independent strict-Laguerre/Newman transversality route remains open.
