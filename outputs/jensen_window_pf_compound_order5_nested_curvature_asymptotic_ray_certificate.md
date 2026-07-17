# Jensen-Window PF Compound Order-Five Nested Curvature Asymptotic-Ray Certificate

Date: 2026-07-13

Status: rigorous nested first-summand order-five curvature theorem on
the complete `u>=20` ray. This is not a proof of full order-five
entry by itself, PF-infinity, RH, or `Lambda <= 0`.
This is not by itself a proof of the full compact-to-kernel composition.

```text
work/rh_compute/results/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_nested_curvature_asymptotic_ray_certificate.py
```

## Shifted Normalized H Boxes

The completed order-four ray theorem gives, throughout the enlarged
`u>=19` collar,

```text
x_r=(-1)^r*t^(r-1)*H^(r)/(r-2)!
0<x_r<=1 for 2<=r<=8, x_2>=97/100, x_3>=24/25,
1/t<=10^-30.
```

Since `t(20)-t(19)>6.84e36`, every point within `t+-3` of a
central `u>=20` point remains in that collar. Exact rational arithmetic
also gives

```text
(1-3*10^-30)^(-7)=1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000/999999999999999999999999999979000000000000000000000000000188999999999999999999999999999055000000000000000000000000002834999999999999999999999999994897000000000000000000000000005102999999999999999999999999997813<1001/1000,
(1+3*10^-30)^(-2)=1000000000000000000000000000000000000000000000000000000000000/1000000000000000000000000000006000000000000000000000000000009>999/1000.
```

For `b_r=t^(r+1)B^(r)/r!`, unit-mass tent averaging therefore yields

```text
b_0 in [969/1000,1001/1000],
b_1 in [-1001/1000,-959/1000],
(-1)^r*b_r in [0,1001/1000], 2<=r<=6.
```

## Analytic Stable Logs

Put `z=1/t` and scale the local variable by `y=s/t`. Then the Taylor
jet of `B` is `z*b`. For

```text
R_0(w)=log((1-exp(-w))/w),
```

the existing convergent product series proves `|R_0^(k)(w)|<1`
through order six whenever `0<=w<=1/1000`. If the normalized Taylor
coefficients of `v` are bounded by `C`, the exact partial-Bell identity
gives the coefficient correction

```text
coefficient r error <=2^(r-1)*C*10^-30 for r>=1 and <=C*10^-30 for r=0
```

The proof uses `C=2` at the defect layer and `C=100` at each nested
stable layer. Arb verifies both coefficient caps and keeps `z*C<1/1000`.
The constants `log(z)` disappear from every derivative used below.

## Dimensionless Interval

The normalized recursion is

```text
J_y,r=z*j_r, j_r=2*b_r-z*(r+1)*(r+2)*ell_(r+2)
R_y,r=z*r_r, r_r=3*b_r-z*(r+1)*(r+2)*h_(r+2)
q=2*h-ell+log(z)+log(r)+R(z*r)
```

A single outward-rounded interval evaluation over the complete box
`0<=z<=10^-30` gives

```text
j_0 lower=1.93799999964423477649688720703123724792404705658555030822753E+0
r_0 lower=2.90699999952455982565879821777339924377214116975665092468261E+0
t^2*q_1'' ball=[+/- 9.16]
t^2*q_1'' upper=9.15835270172636955976486206054574636687947832827158772903837E+0<10,
margin below 10=8.41647283372469246387481689454253633120521671728412270961638E-1.
```

This is a uniform interval theorem, not evaluation at a representative
large value of `t`. Consequently

```text
t^2*q_1''(t)<10<60 for every mode u>=20.
```

Together with the compact and finite-ray certificates, this covers all
real `t>=320`. The full composition is recorded separately.

```text
outputs/jensen_window_pf_compound_order5_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order5_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order5_first_summand_curvature_bridge.md
outputs/signed_hankel_jensen_bridge_target.md
outputs/formal_core.md
```
