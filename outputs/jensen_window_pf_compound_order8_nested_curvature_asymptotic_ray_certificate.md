# Jensen-Window PF Compound Order-Eight Nested Curvature Asymptotic-Ray Certificate

Date: 2026-07-13

Status: rigorous nested first-summand order-eight curvature theorem on
the complete `u>=20` ray. This certificate is not a proof of full
order-eight entry, PF-infinity, RH, or `Lambda <= 0`.

## Shifted H Boxes

The completed order-four ray theorem gives signed unit boxes through
`H^(8)`. Exact coarse corridors through order fourteen give

```text
|kappa_r|*q^(r/2-1)/(r-2)!<50000, r=9,10,
|kappa_r|*q^(r/2-1)/(r-2)!<700001, r=11,12,
|kappa_r|*q^(r/2-1)/(r-2)!<1, r=13,14.
```

The exact ray geometry gives, throughout the `t+-6` collar,

```text
raw normalized H13 cap=84249300895685402333762642401/79900668578288412100000000000,
raw normalized H14 cap=16055202125671593335986291122601/15181127029874798299000000000000,
|b_7|,|b_8|<2500, |b_9|,|b_10|<40000, |b_11|,|b_12|<2.
```

## Dimensionless Stable Recursion

Every stable logarithm uses the convergent defect series and the
exact partial-Bell identity in the form

```text
log(1-exp(-z*v))=log(z)+log(v)+R(z*v),
R(w)=log((1-exp(-w))/w).
coefficient r error <=2^(r-1)*C*10^-30 for r>=1 and <=C*10^-30 for r=0
```

One outward-rounded interval evaluation over `0<=z<=10^-30` gives

```text
j_0 lower=1.93799999964423477649688720703123724792404705658555030822753E+0
r_0 lower=2.90699999952455982565879821777339924377214116975665092468261E+0
s_0 lower=3.87599999928846955299377441406242348754462460055947303771972E+0
t_0 lower=4.84499999911058694124221801757799747924115508794784545898437E+0
u_0 lower=5.81399999904911965131759643554668371886155661195516586303710E+0
t^2*s_1'' ball=[+/- 1.35E+2]
t^2*s_1'' upper=1.34489839907184243202209472656240124460170043090398255726938E+2<200<4000,
margin below 200=6.55101598543971776962280273437598755398299569096017442730626E+1.
```

This is a uniform interval theorem, not evaluation at a representative
large value of `t`. Together with the compact and finite-ray
certificates it covers every real `t>=999` for the first summand.

```text
outputs/jensen_window_pf_compound_order8_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order8_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order8_first_summand_curvature_bridge.md
outputs/formal_core.md
```
