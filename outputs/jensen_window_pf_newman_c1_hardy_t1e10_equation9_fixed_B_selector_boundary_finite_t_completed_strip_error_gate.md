# Fixed-B selector-boundary finite-t completed-strip error

Date: 2026-08-12

Status: rigorous one-height completed-strip comparison; this is not a proof of
a nonzero height cell, the ordinary-Morse join, `T_upper`, `Lambda<=0`, or RH

At the exact selector boundary

```text
A=159577, B=5122423,
t*=pi*A^2/8=[10000011009.0425883324999665107390608929017717861783444351676067861157838771817261874891645 +/- 3.27e-80],
beta=t*^(1/3)=[2154.43548064036880514538712430952510213643988046170017933872008509078888609877630125920487 +/- 1.17e-86].                         (FT1)
```

The exact identity `hY=2pi` makes every complete symmetric mode sum a Fourier
midpoint.  Adding the endpoint half-current evaluates the completed strip at
`y=0`.  This applies to the exact finite-`t` strip and to the canonical
Fourier object before either is estimated.  After suppressing their common
carrier, the comparison is therefore only

```text
J_beta=int_R cosh(z/beta)^(-3/2)
       exp(i*beta^3[z/beta-tanh(z/beta)])dz,
J_0=2pi*Ai(0).                                         (FT2)
```

Thus the 399-mode core and the zero/negative/outer-positive complement are
not estimated separately.  Their exact completion from Section 11.333 is
used first.

For the positive half-line deform to `z=exp(i*pi/6)r`.  With
`rho=r/beta`, the exact decay is

```text
D(rho)=rho/2-sin(rho)/(cosh(sqrt(3)rho)+cos(rho)).     (FT3)
```

Interval subdivision plus an analytic tail proves

```text
D(rho)>=rho^3/4,  0<=rho<=1,
D(rho)>=rho/5,    rho>=1.                              (FT4)
```

The compact proof uses 256 panels on `1/2<=rho<=1`
and 128 panels on `1<=rho<=6/5`; the small range is
handled by the degree-17 tanh polynomial and its Cauchy tail.  Also
`|cosh(rho exp(i*pi/6))|>=1`, so no amplitude growth is discarded.

Twenty-four interval panels on `0<=r<=12`, followed by the explicit
polynomial and contour-tail bounds, give

```text
J_beta-J_0
 =[1.21310104929970121481949162693097801772702871840046930443976781381094935040935440624574124e-15 +/- 2.93e-29].       (FT5)
```

The total error radius appended to the direct quadrature is

```text
[1.92760602634080897875175554388423514132159494709658484176136177638403287125007705133930109e-29 +/- 1.77e-117].                (FT6)
```

The apparent `beta^(-2)` correction cancels exactly after the two Airy rays
are recombined.  The first formal nonzero integrated term is

```text
-(9/560)*2pi*Ai'(0)/beta^4
 =[1.21310170301935625353011202762802769279599339278848656321623218664090887065250360264068247e-15 +/- 2.38e-104],              (FT7)
```

and the rigorous result divided by (FT7) is
`[0.999999461117189587953363933325780926741899179432948585599660873413085937500000000000000000 +/- 2.42e-14]`.

Restoring the transformed Kummer scale gives

```text
canonical completed strip
 =[116.832681847458611358570975865305759328684820915749960391047032569169702056774263809084607 +/- 9.54e-88],
exact finite-t completed strip
 =[116.832681847458674894412145923617298122502411496750206730658088167017053241429538252379402 +/- 1.54e-27],
finite-t minus canonical
 =[6.35358411700583115387938175905810002463396110555978473511846552744432947942965977354078551e-14 +/- 1.54e-27]. (FT8)
```

After the paper's equation-(9) normalization, the certified difference is

```text
[1.12465323837736911257466151907972231972557091782126906396733440983103588284330454391441735e-16 +/- 2.72e-30] < 1.2e-16. (FT9)
```

There is no hidden phase in (FT9): `A=159577` is `1 mod 8`, so
`exp(i[t*-pi/8])=1` exactly.

This closes the exact finite-`t` versus canonical comparison at the single
selector-boundary height.  It does not yet provide a uniform neighborhood in
height or join the completed strip to the ordinary-Morse modes.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate.py
```

No nonzero selector-height cell, complete `T_upper`, `Lambda<=0`, RH, or
prize-level conclusion is proved.
