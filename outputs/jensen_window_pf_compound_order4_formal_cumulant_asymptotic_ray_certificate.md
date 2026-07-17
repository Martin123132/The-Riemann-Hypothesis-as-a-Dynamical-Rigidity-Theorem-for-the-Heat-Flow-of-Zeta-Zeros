# Jensen-Window PF Order-Four Formal Cumulant Asymptotic-Ray Certificate

Date: 2026-07-13

Status: exact analytic theorem for the epsilon-six formal cumulant ray.
This is not a proof of the exact cumulant ray, order-four entry,
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.py
```

## Leading Model

For the q-leading potential jets,

```text
P_0=1; P_(r+1)=(u/2)*P_r'+2u*P_r
corridor floor+1/(10u) <= F_r <= corridor ceiling-1/(10u)
```

After `u=20+v`, clearing positive denominators gives coefficient-positive
numerators for both sides of all seven buffered corridor inequalities.
Odd cumulant ratios are squared only after their positive numerators are
identified. The same method proves `0<L_r^(infinity)<3/2` through order eight.

## Exact Potential Transfer

With `q=pi*exp(4u)`, exact symbolic jets and the substitution
`u=20+v, q=10^35+Q` prove

```text
|V^(r)-q*P_r(u)|<=1000*u^8, r=2,...,8
eta(u)=250*u^6/q<=10^-24
|L_r-L_r^(infinity)|<=14*eta(u)
```

Arb gives `q(20)>100000000000000000000000000000000000` with lower enclosure
`1.74063785791258146612293033454152703544959046008510652604992E+35`. The logarithmic derivative
`d log(q/u^7)/du=4-7/u` is positive on the ray.

Exact coefficient norms give

```text
|R_r^[6]-F_r|<=22000000*u^6/q<1/(20u), u>=20
```

This consumes less than half the `1/(10u)` leading buffer. Therefore the
exact epsilon-six formal cumulant polynomial satisfies every candidate
corridor for all `u>=20`.

## Remaining Boundary

Together with the 1.8-million-block compact formal certificate, the formal
model is now closed for every `u>=2`. What remains is not formal algebra:
prove that the exact standardized-density cumulants differ from the formal
polynomial by less than the available scaled margins, using a central
expansion remainder and two adaptive tails.

```text
outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md
outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md
outputs/formal_core.md
```
