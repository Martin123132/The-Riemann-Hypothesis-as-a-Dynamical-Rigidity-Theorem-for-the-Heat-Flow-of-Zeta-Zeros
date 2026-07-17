# Jensen-Window PF Compound Order-Six Nested Curvature Asymptotic-Ray Certificate

Date: 2026-07-13

Status: rigorous nested first-summand order-six curvature theorem on
the complete `u>=20` ray. This certificate is not a proof of full
order-six entry, PF-infinity, RH, or `Lambda <= 0`.

```text
work/rh_compute/results/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.py
python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_nested_curvature_asymptotic_ray_certificate.py
```

## Shifted H Boxes

The completed order-four ray theorem gives signed unit boxes through
`H^(8)`. The exact order-six coarse corridor adds

```text
|kappa_r|*q^(r/2-1)/(r-2)!<50000, r=9,10.
```

The proved geometry

```text
t^(r-1)q^(1-r/2)/V''^(r/2)<=D_r/u,
D_r=(201/100)^(r-1)/(19/10)^r, u>=19,
```

and the midpoint Hurwitz ceiling give, throughout the `t+-4` collar,

```text
raw normalized H9 cap=2665436245700681801/1226213251560200,
raw normalized H10 cap=535739197040069879801/232980517796438000,
|b_7|<2500, |b_8|<2500.
```

The low derivatives retain

```text
b_0 in [969/1000,1001/1000],
b_1 in [-1001/1000,-959/1000],
(-1)^r*b_r in [0,1001/1000], 2<=r<=6.
```

## Dimensionless Stable Recursion

Put `z=1/t` and use the scaled local variable `y=s/t`. Each stable
logarithm is evaluated as

```text
log(1-exp(-z*v))=log(z)+log(v)+R(z*v),
R(w)=log((1-exp(-w))/w).
```

The convergent defect series and exact partial-Bell identity give

```text
coefficient r error <=2^(r-1)*C*10^-30 for r>=1 and <=C*10^-30 for r=0
```

for every layer. A single outward-rounded interval evaluation over
`0<=z<=10^-30` then gives

```text
j_0 lower=1.93799999964423477649688720703123724792404705658555030822753E+0
r_0 lower=2.90699999952455982565879821777339924377214116975665092468261E+0
s_0 lower=3.87599999928846955299377441406242348754462460055947303771972E+0
t^2*p_1'' ball=[+/- 22.8]
t^2*p_1'' upper=2.27683610696345567703247070312471784171958323101989922910682E+1<100<200,
margin below 100=7.72316389005631208419799804687528215828041676898010077089318E+1.
```

This is a uniform interval theorem, not evaluation at a representative
large value of `t`. Together with the compact and finite-ray
certificates it covers every real `t>=321`.

```text
outputs/jensen_window_pf_compound_order6_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order6_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order6_first_summand_curvature_bridge.md
outputs/formal_core.md
```
