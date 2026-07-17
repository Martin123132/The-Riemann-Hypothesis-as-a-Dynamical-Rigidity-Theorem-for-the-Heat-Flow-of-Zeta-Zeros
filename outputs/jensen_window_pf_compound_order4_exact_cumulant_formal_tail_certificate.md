# Jensen-Window PF Order-Four Exact Cumulant Formal-Tail Certificate

Date: 2026-07-13

Status: exact formal two-tail theorem with open exact-density components.
This is not a proof of the exact cumulant corridors, order-four entry,
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_exact_cumulant_formal_tail_certificate.py
```

## Formal Weight

The epsilon-ten formal density factor has degree thirty in `y`. The global
normalized-jet caps give the exact coefficient majorant

```text
sum_(n,m)|[y^m]E_n|=32408700190646835653/24385536000000000000<4/3.
```

Use the adaptive cutoff

```text
Y(q)=1+sqrt(32*log(q))
exp(-Y^2/2+Y)=exp(1/2)*q^-16
Y(q)<2*q^(1/4).
```

Exact exponential-tail integration through polynomial degree thirty has
maximum ratio below `3/2`. Including the complex factors on `|z|<=1` gives

```text
formal tail <=2^32*q^(-35/4) for |z|<=1
each formal tail <1/(500000*q^3)
each formal tail <1/(100000*u*q^3)
```

Both formal tails are therefore closed inside the unit-disk partition
contract.

## Remaining Boundary

Three exact-density components remain: the central exact-minus-formal
residual, the exact left tail, and the exact right tail. Their sufficient
allocations are

```text
finite: each <1/(500000*q^3), 2<=u<=20,
ray:    each <1/(300000*u*q^3), u>=20.
```

These are open targets. No exact-density tail or cumulant corridor is
promoted by this formal-tail theorem.

```text
outputs/jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract.md
outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md
outputs/formal_core.md
```
