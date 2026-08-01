# Newman Real-Zero-Band Degree-361 Sector Certificate

Date: 2026-07-25

Status: rigorous Arb/Taylor boundary certificate, zero-homotopy theorem,
and published-sector composition. This is not a proof of PF-infinity,
`Lambda<=0`, RH, or a Clay-prize conclusion.

```text
work/rh_compute/results/jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.py
```

Current result:

```text
validated Newman real-zero-band degree-361 sector certificate: 12 rows, 0 issues, 200 certified boundary boxes, 0 unresolved, 1 real-zero band, all shifts through degree 361, 0 all-degree theorems, weakest negative margin [1.413787535739456082478734678236082408014047530929237461735354861030476e-5 +/- 8.48e-76]
```

## Vertical Boundary

For

```text
H_t(z)=integral_0^infinity exp(t*u^2)*Phi(u)*cos(z*u)du,
```

the real part on `z=38+iy` is

```text
Re H_t(38+iy)=integral_0^infinity
 exp(t*u^2)*Phi(u)*cos(38*u)*cosh(y*u)du.
```

A 192-bit outward-rounded Taylor certificate proves

```text
Re H_t(38+iy)<0
for every 0<=t<=1/5 and |y|<=1.
```

The partition has `200` rational boxes, Taylor orders
`(t,y)=(2,4)`, steps
`(1/50,1/20)`, and zero unresolved boxes. The weakest
certified negative margin is

```text
[1.413787535739456082478734678236082408014047530929237461735354861030476e-5 +/- 8.48e-76]
```

on `t=[9/50,1/5]`,
`y=[19/20,1]`. Four theta summands are integrated
directly. The omitted `n>=5` terms use

```text
2*pi^2*5^4*exp(9u-25*pi*exp(4u))
 /(1-(6/5)^4*exp(-11*pi)),
```

For every required `0<=m<=9`, `u^m<=exp(9u/2)` and
`cosh(yu)<=exp(u)`. The partial sum
`sum_(k=0)^7 8^k/k!=47259/35`, together with `pi>3` and `pi^2<10`,
gives the directed audit

```text
u^2/5+(29/2)u-pi*exp(4u)<=-140734/35-(u-2),
40*10^800<2^4020,
```

so the analytic `u>=2` moment tail is strictly below `1e-800`.
Evenness and conjugation give the same zero-free conclusion on `Re z=-38`.

## Real-Zero Homotopy

Polymath Theorem 3.2 gives

```text
H_t(z)=0 => |Im z|<=sqrt(1-2t)<1, 0<t<=1/5.
```

At `t=0`, `H_0(z)=xi((1+iz)/2)/8`; the horizontal sides `|Im z|=1`
map to the classical zero-free lines `Re s=0,1`. Thus the complete
rectangle boundary is zero-free. Platt--Trudgian Corollary 2 gives
`Lambda<=1/5`, so every zero at the top time is real. The independent
1,900-box compact transversality certificate proves

```text
(H_t(x),H_t'(x))!=(0,0)
for 0<=t<=1/5 and |x|<=38.
```

The argument principle fixes the finite zero multiset in the rectangle.
Simple real roots continue as real roots, and a pair could leave the real
axis only through a real multiple zero. Therefore

```text
H_t(z)=0 and |Re z|<=38 imply Im z=0
for every 0<=t<=1/5.
```

## Sector Transfer

All remaining nonreal zeros have `|Re z|>38`, while the strip theorem
gives `|Im z|<=1`. Under `s=-z^2`, every zero of
`F_t(s)=2H_t(i*sqrt(s))` therefore lies in the negative-axis sector

```text
delta=2*atan(1/38),
sin(delta)=76/1445,
|sin(delta)|^(-2)=2088025/5776=361.500... .
```

Chasse's sector theorem, the order-one-half property, coefficient
positivity, and derivative closure give

```text
P_(d,n,t)(x)=sum_(j=0)^d binom(d,j)A_(n+j)(t)x^j
is hyperbolic with real negative zeros
for every 0<=d<=361, n>=0, and 0<=t<=1/5.
```

Primary theorem inputs:

- D. H. J. Polymath, Theorem 3.2:
  <https://arxiv.org/abs/1904.12438>
- Platt--Trudgian, Corollary 2:
  <https://doi.org/10.1090/mcom/3595>
- Matthew Chasse, Theorem 3.6:
  <https://doi.org/10.1080/17476933.2011.584250>

## Proof Boundary

This strictly sharpens the independent degree-71 zero-slab theorem.
The cutoff `361` is finite. Degree `362`, an unbounded cofinal degree
sequence, all-degree Jensen hyperbolicity, PF-infinity, `Lambda<=0`, RH,
and the Clay prize remain open. The next sector upgrade needs a wider
real-zero band or a genuinely unbounded sector mechanism.
