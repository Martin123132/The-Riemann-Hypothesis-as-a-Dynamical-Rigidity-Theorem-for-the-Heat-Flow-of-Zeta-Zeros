# Jensen-Window PF Compound Order-Seven Nested Curvature Asymptotic-Ray Certificate

Date: 2026-07-14

Status: rigorous nested first-summand order-seven curvature theorem on
the complete `u>=20` ray. This certificate is not a proof of full
order-seven entry, PF-infinity, RH, or `Lambda <= 0`.

## Shifted H Boxes

The completed order-four ray theorem gives signed unit boxes through
`H^(8)`. Exact coarse corridors add

```text
|kappa_r|*q^(r/2-1)/(r-2)!<50000, r=9,10,
|kappa_r|*q^(r/2-1)/(r-2)!<700001, r=11,12.
```

The exact ray geometry and Hurwitz midpoint ceiling give, throughout
the `t+-5` collar,

```text
raw normalized H11 cap=75348045616352460399362502001/2213314919066161000000000,
raw normalized H12 cap=15144932822422734812500862902201/420529834622570590000000000,
|b_7|,|b_8|<2500, |b_9|,|b_10|<40000.
```

## Dimensionless Stable Recursion

Put `z=1/t`. Every stable logarithm is evaluated as

```text
log(1-exp(-z*v))=log(z)+log(v)+R(z*v),
R(w)=log((1-exp(-w))/w).
```

The convergent defect series and exact partial-Bell identity give

```text
coefficient r error <=2^(r-1)*C*10^-30 for r>=1 and <=C*10^-30 for r=0
```

for every layer. A single outward-rounded interval evaluation over
`0<=z<=10^-30` gives

```text
j_0 lower=1.93799999964423477649688720703123724792404705658555030822753E+0
r_0 lower=2.90699999952455982565879821777339924377214116975665092468261E+0
s_0 lower=3.87599999928846955299377441406242348754462460055947303771972E+0
t_0 lower=4.84499999911058694124221801757799747924115508794784545898437E+0
t^2*r_1'' ball=[+/- 55.6]
t^2*r_1'' upper=5.55401953430641442537307739257756068343861012881347792016427E+1<100<600,
margin below 100=4.44598045973312109708786010742243931656138987118652207983573E+1.
```

This is a uniform interval theorem, not evaluation at a representative
large value of `t`. Together with the preceding certificates it covers
every real `t>=320` for the first Newman summand.

```text
outputs/jensen_window_pf_compound_order7_nested_curvature_compact_certificate.md
outputs/jensen_window_pf_compound_order7_nested_curvature_finite_ray_certificate.md
outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md
outputs/formal_core.md
```
