# Jensen-Window PF Order-Four Formal Cumulant Next-Parity Certificate

Date: 2026-07-13

Status: exact formal epsilon-eight and next-parity coefficient theorem.
This is not a proof of the exact cumulant ray, order-four entry,
PF-infinity, RH, or `Lambda <= 0`.

Artifact kind: `jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate`.

```text
work/rh_compute/results/jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.json
python work/rh_compute/scripts/jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate.py
```

## Exact Audit

The tilted-Gaussian recurrence was extended from `epsilon^6` to
`epsilon^8`, with normalized potential jets `L_3,...,L_10`. All 42
coefficients at epsilon powers one through six agree exactly with the
stored formal-cumulant target.

Reflection parity gives

```text
[epsilon^j] kappa_r=0 when j is not congruent to r modulo 2.
```

## First Omitted Layer

After applying the corridor normalization, the epsilon-eight extension is
organized as

```text
orders 2,3,4: scaled(kappa_r^[8]-kappa_r^[6])=q^-3*C_r
orders 5,6:   scaled(kappa_r^[8]-kappa_r^[6])=q^-2*C_r
orders 7,8:   scaled(kappa_r^[8]-kappa_r^[6])=q^-1*C_r
```

The exact rational polynomials `C_r(L_3,...,L_10)` are stored in the JSON
artifact. In particular, the previously observed order-seven and
order-eight `q^-1` discrepancy is the predicted first omitted parity term,
not a failure of the epsilon-six algebra.

## Remaining Boundary

No exact-density remainder has been promoted. The next gates are to certify
the seven explicit coefficient functions on `2<=u<=20` and `u>=20`, then
bound the central remainder beyond epsilon eight and both adaptive tails.

```text
outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md
outputs/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.md
outputs/formal_core.md
```
