# Jensen-Window PF Compound Order-Four Exact-Corridor Localized-Curvature Ray Certificate

Date: 2026-07-13

Status: rigorous exact-corridor localized-curvature theorem on `u>=20`,
composed with the finite theorem into a global corridor-to-`U` theorem on
`u>=2`. This is not a proof of order-four entry, PF-infinity,
RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_exact_corridor_localized_curvature_ray_certificate.py
```

## Ray Geometry

Write `q=pi*exp(4u)`, `t=V'(x_u)`, and `a=V''(x_u)`. On the
slightly enlarged mode ray `u>=19`, exact coefficient-positive gates and
the monotonicity of `q/u` give

```text
q >= 1000000000000000000000000000000000
u*q <= t <= (201/100)*u*q
a >= (361/100)*u^2*q
1/t <= 10^-30.
```

The strict endpoint comparison

```text
t(20)-t(19) > 6.84140403277229664153888231548284194636517090552223519946531E+36
```

shows that every `t+s`, `|s|<=1`, belonging to a central mode `u>=20`
has a shifted mode above `19`. Thus the same analytic boxes control all
three Taylor envelopes without extrapolating a sampled cover.

## Normalized H Jets

For `2<=r<=8` put

```text
x_r=(-1)^r*t^(r-1)*H^(r)/(r-2)!.
```

Convex midpoint quadrature and the integral lower bound for the Hurwitz
zeta function give

```text
(t/(t+1/2))^(r-1) <= (r-1)t^(r-1)zeta(r,t+1/2) <= 1.
```

Combining this with the seven exact cumulant corridors and the ray
geometry proves

```text
0 < x_r <= 1       (2<=r<=8)
x_2 >= 97/100
x_3 >= 24/25.
```

## Logarithmic Defect

With `B=H''` split

```text
ell=log(1-exp(-B))=log(B)+R(B)
R(B)=log((1-exp(-B))/B).
```

The convergent product-series expansion, `zeta(2k)<2`, and `2*pi>6`
give `|R^(m)(B)|<1` for `1<=m<=6`. The exact partial-Bell identity
then yields

```text
t^2*ell'' <= 23/20
|t^r*ell^(r)| < 30000       (2<=r<=6)
t^(r+1)*E_r < 1/1000        (r=0,1,2).
```

## Curvature Comparison

The scaled localized quantities satisfy

```text
A_-=t*(j_0-E_0)                 > 193/100
A_+=t*(j_0+E_0)                 < 201/100
C_+=t^3*max(j_2+E_2,0)          < 401/100
P_-=t^2*max(abs(j_1)-E_1,0)     > 191/100.
```

For `0<z<=1/1000`,

```text
z/(exp(z)-1) <= 1
z^2*exp(z)/(exp(z)-1)^2 >= exp(-z) >= 1-z >= 999/1000.
```

Crucially, the negative square term is retained. The final comparison is

```text
t^2*U(t)
 < 2*(23/20) + 401/193 - (999/1000)*(191/201)^2
 = 3011223637/866377000
 = 3.4756504812570048
 < 7/2,
margin > 21095863/866377000.
```

Therefore `j_0>E_0` and `U(t)<7/(2t^2)` for every central mode
`u>=20`. The checked finite certificate supplies the same result on
`2<=u<=20`, so the exact-corridor-to-localized-curvature implication is
now complete for every `u>=2`.

## Next Composition

The next step is bookkeeping with mathematical content already isolated:
compose the compact curvature theorem, this global theorem, the exact tent
transfer, and the full-kernel perturbation into contiguous order-four entry
at `lambda=-100`.

```text
outputs/jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem.md
outputs/jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate.md
outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md
outputs/formal_core.md
```
