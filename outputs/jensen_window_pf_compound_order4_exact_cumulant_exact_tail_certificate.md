# Jensen-Window PF Order-Four Exact Cumulant Exact-Tail Certificate

Date: 2026-07-13

Status: exact-density two-tail theorem with one open central residual.
This is not a proof of the exact cumulant corridors, order-four entry,
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate.py
```

## Curvature Geometry

The exact curvature formula and two positive-coefficient substitutions prove

```text
V''(x_u)<=5*u^2*q                  (u>=2, q>=9000),
V''(x_v)>=(39/10)*v^2*q_v          (v>=39/20, q_v>=7200).
```

For `Y(q)=1+sqrt(32*log(q))`, the adaptive shift obeys

```text
u_-/u>39/40,
q_-/q>4/5,
W_u''(y)>59319/100000>1/2          (|y|<=Y).
```

Hence `W_u(+-Y)>=Y^2/4`, with outward slopes at least `Y/2`.
The inherited global potential geometry propagates those slopes through
both complete tails.

## Exact Tails

Uniformly for `|z|<=1`, either exact-density tail satisfies

```text
exact tail <exp(Y-Y^2/4)=exp(3/4+sqrt(32*log(q))/2)*q^-8
each exact tail <1/(500000*q^3)
each exact tail <1/(300000*u*q^3)
```

Together with the formal-tail theorem, all four tails in the complex-disk
partition decomposition are closed.

## Remaining Boundary

Only the central exact-minus-formal residual remains:

```text
finite: <1/(500000*q^3),             2<=u<=20,
ray:    <1/(300000*u*q^3),           u>=20.
```

This central residual is an open target. No exact cumulant corridor is
promoted by the two-tail theorem.

```text
outputs/jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md
outputs/formal_core.md
```
