# Jensen-Window PF Order-Four Exact Cumulant Complex-Disk Contract

Date: 2026-07-13

Status: exact complex-disk reduction with open partition residual.
This is not a proof of the exact cumulant corridors, order-four entry,
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.py
```

## Exact Factorization

At the mode, put

```text
W_u(y)=V(x_u+y/sqrt(a_u))-V(x_u)-t_u*y/sqrt(a_u)
A_u(z)=exp(-z^2/2)/sqrt(2*pi)*integral_R exp(z*y-W_u(y))dy
K_u(z)=z^2/2+log A_u(z)-log A_u(0)
```

Exact graded Gaussian algebra constructs `P_u^[10](z)` and
`S_u^[10](z)` and audits all 70 cumulant coefficients through epsilon ten.
The global normalized-jet caps give

```text
|P_u^[10](z)-1|<1/100, |z|<=1, u>=2.
```

Thus the formal partition has no unit-disk zero. Exact positive-majorant
recurrences prove

```text
|log P^[10]-S^[10]|<1/(7500*q^3)
|log P^[10]-S^[10]|<1/(100000*u*q^3)
```

## Cauchy Reduction

The relative partition logarithm lemma and Cauchy's estimate reduce the
simultaneous cumulant problem to the following sufficient targets:

```text
sup_|z|<=1 |A_u(z)-P_u^[10](z)|<1/(100000*q^3), 2<=u<=20,
sup_|z|<=1 |A_u(z)-P_u^[10](z)|<1/(20000*u*q^3), u>=20.
```

The resulting finite scaled cumulant cap is `532/61875`
against `9/1000`; the remaining margin is `199/495000`.
The ray cap is `7693/1237500/u` against `1/(100u)`,
leaving `2341/618750/u`.

## Remaining Boundary

The targets above are not yet proved. Their advantage is structural: the
epsilon-ten subtraction occurs at partition level before differentiation,
so the order-seven/eight cancellations are preserved automatically. The
next proof must split the exact central residual, exact left tail, exact
right tail, and formal Gaussian tails without returning to independent raw
moment boxes.

```text
outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_asymptotic_ray_certificate.md
outputs/formal_core.md
```
