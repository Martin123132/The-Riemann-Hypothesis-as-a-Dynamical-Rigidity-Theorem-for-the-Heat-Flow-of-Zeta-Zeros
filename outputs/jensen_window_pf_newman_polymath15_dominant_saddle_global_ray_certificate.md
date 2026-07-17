# Jensen-Window PF Newman Polymath-15 Dominant-Saddle Global Ray Certificate

Date: 2026-07-17

Status: strict first Laguerre positivity for the exact Newman heat-flow
function on the complete dominant-saddle ray, including every cutoff
transition. This is not a proof of `Lambda <= 0` or RH.

```text
work/rh_compute/results/jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate.json
python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate.py
```

The effective approximation is imported from D. H. J. Polymath,
[Effective approximation of heat flow evolution of the Riemann xi function](https://arxiv.org/abs/1904.12438),
Theorem 1.3.

## Region

```text
0<t<=1/2, L=log(x/(4*pi)), t*L>=25, N_x=floor(sqrt(x/(4*pi)+t/16))
```

The prior certificate proved the exact theorem whenever a radius-`1/4`
disk stayed inside one cutoff cell. This certificate removes that
restriction.

## One-Block Repair

Consecutive boundaries have spacing

```text
x_(N+1)-x_N=4*pi*(2N+1); for L>=50 a radius-1/4 disk crosses at most one cutoff boundary
```

At `L=50` the relevant adjacent index is already enclosed below by
`[72004899336.38587252416135146612615791522353381339527873622138644723206 +/- 6.92e-61]`. Thus a radius-`1/4` disk
crosses at most one boundary. Keep the center cutoff fixed throughout
the disk. Where the theorem prescribes its neighbor, the two main sums
differ by one block and

```text
|f_(N_z,t)(z)-f_(N_x,t)(z)|<=3*n^(-7/2), n>=exp(L/2)-1
```

After restoring the raw normalizer, Arb gives

```text
adjacent raw jump / A(center) = [8.751811692390317453095950416291923816104918478522244230782868344244592e-37 +/- 4.19e-107]
|S_(N_z,t)(z)-S_(N_x,t)(z)|/A_t(x)<10^-30
```

This is negligible beside the existing fixed-cell collar. Their
outward-rounded sum is

```text
[0.0001144926418905000756380979082439192190885325010951787450984725071262191 +/- 5.21e-74] < 1/8000
```

so all previous Cauchy, normalized-jet, and curvature budgets apply
unchanged at the transition itself.

## Global Ray

The exact rational margins remain

```text
finite-main margin / L^2 = 19/250
remainder error / L^2  < 1/1000
exact-H margin / L^2   = 3/40
```

Therefore, with no cutoff exception,

```text
For every real x on the ray, C_t[H_t/A_t](x)>3/40*L^2 and L_t(x)>0
```

## Remaining Gap

The high-frequency dominant-saddle branch is now closed. The remaining
global target is exactly the nonuniform layer

```text
0<t*log(x/(4*pi))<25
```

where the finite Dirichlet terms no longer form a small perturbation
of the first saddle and the RH-level arithmetic cancellation must be
retained.
