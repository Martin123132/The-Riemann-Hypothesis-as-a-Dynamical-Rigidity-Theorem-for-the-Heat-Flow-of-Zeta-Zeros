# Jensen-Window PF Newman Polymath-15 Dominant-Saddle Cell-Interior Certificate

Date: 2026-07-16

Status: strict first Laguerre positivity for the exact Newman heat-flow
function on dominant-saddle fixed-cutoff cell interiors. This is not a
proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate.py
```

The scalar approximation is imported from D. H. J. Polymath,
[Effective approximation of heat flow evolution of the Riemann xi function](https://arxiv.org/abs/1904.12438),
Theorem 1.3. All derivative and curvature transfers below are explicit.

## Region

```text
0<t<=1/2, L=log(x/(4*pi)), t*L>=25, N=floor(sqrt(x/(4*pi)+t/16)), and x is at distance at least 1/2 from both boundaries of its fixed-N cell
4*pi*(N^2-t/16)<=x<4*pi*((N+1)^2-t/16), dist(x,cell boundary)>=1/2
```

Since `t<=1/2`, the first condition forces `L>=50` and
`x>[65152931198542130352545.16179758607386787326595894428541600421898494750 +/- 2.12e-48]`, in particular `x>10^22`.
This deliberately large threshold buys a clean theorem; it does not
indicate failure below it.

## Raw Collar

With `N` fixed on the radius-`1/4` disk, define

```text
S_(N,t)(z)=B_t(z)*f_(N,t)(z), R_(N,t)(z)=H_t(z)-S_(N,t)(z)
S_(N,t)(x)=A_t(x)*P_(N,t)(x), A_t=|B_t|, Z_t-P_(N,t)=R_(N,t)/A_t
```

The main sum is holomorphic and real on the real axis. Schwarz
reflection therefore transfers the upper-half Polymath-15 bound to
the lower semicircle. The explicit published errors give

```text
positive eC exponent correction < 1/100
raw eA+eB / A(center) < 1/200000
raw eC0 / A(center) < 3/25000
sup_(|z-x|=1/4)|R_(N,t)(z)|/A_t(x)<1/8000
```

The saved outward-rounded combined bound is
`[0.0001144926418905000756380979082439183439073632620634334355034308779338375 +/- 4.16e-74]`.

## Derivatives

Cauchy's estimate first gives

```text
|R|/A<=delta, |R'|/A<=4*delta, |R''|/A<=32*delta
```

For `ell=(log A)'` and `V=-(log A)''`, exact differentiation of
`e=R/A=Z-P` gives

```text
e'=R'/A-ell*e
e''=R''/A-2*ell*R'/A+(ell^2+V)*e
eps0<=delta, eps1<=0.35*L*delta, eps2<=0.13*L^2*delta, delta=1/8000
```

The arithmetic ray theorem supplies the companion main-jet caps

```text
|P|<=2.256, |P'|<=0.59L, |P''|<=0.154L^2
```

## Exact Heat Flow

Substitution into the exact curvature perturbation identity gives
the rational budget

```text
finite-main margin / L^2 = 19/250
remainder error / L^2  < 1/1000
exact-H margin / L^2   = 3/40
```

Therefore

```text
C_t[Z_t](x)>3/40*L^2 and L_t(x)=A_t(x)^2*C_t[Z_t](x)>0
```

## Remaining Gap

The exact high-frequency theorem now fails only in two explicitly
identified places: width-one neighborhoods around the discrete cutoff
transitions, and the near-zero-time layer `0<tL<25`. The first is a
local adjacent-`N` problem; the second retains the full RH-level
arithmetic coupling.
